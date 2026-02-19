"""Star detection and plate solving pipeline."""

import math
import multiprocessing
import pathlib
import re
import time

import numpy as np
from astropy.wcs import WCS
from PIL import Image


# Quad angular size midpoints (arcmin) for each series and scale
_QUAD_ARCMIN: dict[tuple[str, int], float] = {}

# 4100 series: quad_arcmin = 30 * 2^((scale-10)/2)
for _s in range(7, 20):
    _QUAD_ARCMIN[("4100", _s)] = 30.0 * (2 ** ((_s - 10) / 2.0))

# 5200 series: known angular ranges per scale
_5200_RANGES = {
    0: (2.0, 2.8), 1: (2.8, 4.0), 2: (4.0, 5.6), 3: (5.6, 8.0),
    4: (8.0, 11.0), 5: (11.0, 16.0), 6: (16.0, 22.0),
}
for _s, (_lo, _hi) in _5200_RANGES.items():
    _QUAD_ARCMIN[("5200", _s)] = (_lo + _hi) / 2.0


def _parse_index_file(path: pathlib.Path) -> tuple[str, int] | None:
    """Return (series, scale) for an index file, or None if unrecognized."""
    name = path.name
    m = re.match(r"index-41(\d{2})\.fits$", name)
    if m:
        return ("4100", int(m.group(1)))
    m = re.match(r"index-52(\d{2})-\d{2}\.fits$", name)
    if m:
        return ("5200", int(m.group(1)))
    return None


def _filter_index_files(
    index_files: list[pathlib.Path],
    width: int,
    height: int,
    scale_low: float | None,
    scale_high: float | None,
) -> list[pathlib.Path]:
    """Pre-filter index files to only those whose quad scale overlaps the expected FOV."""
    if scale_low is None and scale_high is None:
        return index_files

    lo = scale_low or 0.1
    hi = scale_high or 1000.0

    diag_px = math.hypot(width, height)
    short_px = min(width, height)

    # Quads should be between ~5% and ~100% of the image shortest dimension
    min_quad_arcmin = short_px * lo / 60.0 * 0.05
    max_quad_arcmin = diag_px * hi / 60.0

    filtered = []
    for f in index_files:
        parsed = _parse_index_file(f)
        if parsed is None:
            filtered.append(f)  # Unknown format — keep to be safe
            continue
        quad_arcmin = _QUAD_ARCMIN.get(parsed)
        if quad_arcmin is None:
            filtered.append(f)
            continue
        if min_quad_arcmin <= quad_arcmin <= max_quad_arcmin:
            filtered.append(f)

    if not filtered:
        return index_files
    return filtered


def _parse_5200_tile(path: pathlib.Path) -> int | None:
    """Extract healpix tile number from a 5200-series filename, or None."""
    m = re.match(r"index-52\d{2}-(\d{2})\.fits$", path.name)
    return int(m.group(1)) if m else None


# Tile priority: loaded once from bundled JSON
_TILE_PRIORITY: dict[int, int] | None = None


def _get_tile_priority() -> dict[int, int]:
    """Load tile priority mapping {tile_number: rank} from bundled JSON."""
    global _TILE_PRIORITY
    if _TILE_PRIORITY is None:
        import json
        priority_path = pathlib.Path(__file__).parent / "data" / "tile_priority.json"
        if priority_path.exists():
            with open(priority_path) as f:
                data = json.load(f)
            _TILE_PRIORITY = {tile: rank for rank, tile in enumerate(data["tile_order"])}
        else:
            _TILE_PRIORITY = {}
    return _TILE_PRIORITY


def _sort_index_files_by_priority(index_files: list[pathlib.Path]) -> list[pathlib.Path]:
    """Sort index files for fastest solving of typical astrophotography.

    Two priority axes:
    1. Scale: middle-out from the most common focal lengths (85-600mm FF).
       4100 center=12 (60' quads), 5200 center=4 (8-11' quads).
    2. Sky region: object-dense healpix tiles first (from tile_priority.json).

    4100 files (all-sky) come first, sorted by scale priority.
    5200 files follow, sorted by tile priority then scale priority.
    """
    tile_priority = _get_tile_priority()

    # Middle-out scale ordering: distance from the "sweet spot"
    # 4100 scales 7-19: center at 12, order: 12,11,13,10,14,9,15,8,16,7,17,18,19
    # 5200 scales 0-6:  center at 4,  order: 4,3,5,2,6,1,0
    _4100_CENTER = 12
    _5200_CENTER = 4

    files_4100 = []
    files_5200 = []
    files_other = []

    for f in index_files:
        parsed = _parse_index_file(f)
        if parsed is None:
            files_other.append(f)
        elif parsed[0] == "4100":
            files_4100.append(f)
        else:
            files_5200.append(f)

    # Sort 4100 by distance from center scale (middle-out)
    files_4100.sort(key=lambda f: (
        abs((_parse_index_file(f)[1] if _parse_index_file(f) else 0) - _4100_CENTER),
        _parse_index_file(f)[1] if _parse_index_file(f) else 0,
    ))

    # Sort 5200: primary = tile priority, secondary = scale distance from center
    npix = max(tile_priority.values(), default=0) + 1 if tile_priority else 0
    files_5200.sort(key=lambda f: (
        tile_priority.get(_parse_5200_tile(f) or -1, npix),
        abs((_parse_index_file(f)[1] if _parse_index_file(f) else 0) - _5200_CENTER),
    ))

    return files_4100 + files_5200 + files_other


def estimate_plate_scale_from_exif(
    image_path: pathlib.Path,
) -> tuple[float, float] | None:
    """Estimate plate scale range from EXIF data.

    Returns (scale_low, scale_high) in arcsec/pixel, or None if insufficient data.
    """
    img = Image.open(image_path)
    exif = img.getexif()
    if not exif:
        return None

    width, height = img.size

    # Method 1: FocalLength + FocalPlaneXResolution (most precise)
    focal_length = exif.get(37386)  # FocalLength in mm
    fp_x_res = exif.get(41486)  # FocalPlaneXResolution
    fp_res_unit = exif.get(41488)  # FocalPlaneResolutionUnit

    if focal_length and fp_x_res and fp_res_unit:
        fl_mm = float(focal_length)
        unit_to_mm = {2: 25.4, 3: 10.0, 4: 1.0, 5: 0.001}
        if fl_mm > 0 and fp_res_unit in unit_to_mm:
            pixel_size_mm = unit_to_mm[fp_res_unit] / float(fp_x_res)
            plate_scale = 206.265 * pixel_size_mm / fl_mm  # arcsec/pixel
            return (plate_scale * 0.7, plate_scale * 1.3)

    # Method 2: FocalLengthIn35mmFilm (less precise but widely available)
    fl_35mm = exif.get(41989)  # FocalLengthIn35mmFilm in mm
    if fl_35mm and int(fl_35mm) > 0:
        fov_h_rad = 2 * math.atan(18.0 / int(fl_35mm))
        plate_scale = math.degrees(fov_h_rad) * 3600.0 / width  # arcsec/pixel
        return (plate_scale * 0.6, plate_scale * 1.5)

    return None


def load_image(image_path: str | pathlib.Path) -> np.ndarray:
    """Load a JPG/PNG image as a grayscale float32 numpy array."""
    img = Image.open(image_path)
    if img.mode != "L":
        img = img.convert("L")
    return np.array(img, dtype=np.float32)


def detect_stars(
    image: np.ndarray, max_stars: int = 200, detection_threshold: float = 5.0
) -> np.ndarray:
    """Detect stars using sep (Source Extractor Python).

    Returns an (N, 2) array of (x, y) pixel coordinates, sorted by brightness
    (brightest first).
    """
    import sep

    # sep requires C-contiguous array
    data = np.ascontiguousarray(image)

    # Background subtraction
    bkg = sep.Background(data)
    data_sub = data - bkg

    # Extract sources
    objects = sep.extract(data_sub, detection_threshold, err=bkg.globalrms)

    if len(objects) == 0:
        return np.empty((0, 2), dtype=np.float64)

    # Sort by flux (brightest first)
    flux_order = np.argsort(-objects["flux"])
    objects = objects[flux_order]

    # Take top N stars
    objects = objects[:max_stars]

    # Return (x, y) positions
    return np.column_stack([objects["x"], objects["y"]])


def _solve_worker(
    index_files: list[pathlib.Path],
    stars: np.ndarray,
    size_hint_args: tuple[float, float] | None,
    position_hint_args: tuple[float, float, float] | None,
    parity: str,
    worker_id: int,
    result_queue: multiprocessing.Queue,
    cancel_event: multiprocessing.Event,
) -> None:
    """Solver worker for parallel solving. Runs in a subprocess.

    Each worker gets a round-robin subset of index files and processes
    all star slices on its subset.
    """
    import astrometry

    parity_map = {
        "both": astrometry.Parity.BOTH,
        "normal": astrometry.Parity.NORMAL,
        "flip": astrometry.Parity.FLIP,
    }

    size_hint = None
    if size_hint_args is not None:
        size_hint = astrometry.SizeHint(
            lower_arcsec_per_pixel=size_hint_args[0],
            upper_arcsec_per_pixel=size_hint_args[1],
        )

    position_hint = None
    if position_hint_args is not None:
        position_hint = astrometry.PositionHint(
            ra_deg=position_hint_args[0],
            dec_deg=position_hint_args[1],
            radius_deg=position_hint_args[2],
        )

    def logodds_callback(logodds_list):
        if cancel_event.is_set():
            return astrometry.Action.STOP
        if logodds_list[0] > 20.0:
            return astrometry.Action.STOP
        return astrometry.Action.CONTINUE

    solution_parameters = astrometry.SolutionParameters(
        parity=parity_map[parity],
        logodds_callback=logodds_callback,
    )

    try:
        with astrometry.Solver(index_files) as solver:
            solution = solver.solve(
                stars=stars,
                size_hint=size_hint,
                position_hint=position_hint,
                solution_parameters=solution_parameters,
            )

        if solution.has_match():
            match = solution.best_match()
            result_queue.put({
                "worker_id": worker_id,
                "wcs_header": dict(match.astropy_wcs().to_header(relax=True)),
                "center_ra_deg": match.center_ra_deg,
                "center_dec_deg": match.center_dec_deg,
                "scale_arcsec_per_pixel": match.scale_arcsec_per_pixel,
                "logodds": match.logodds,
            })
    except Exception:
        pass  # Worker failed silently; main process handles timeout


def _solve_parallel(
    index_files: list[pathlib.Path],
    stars: np.ndarray,
    size_hint_args: tuple[float, float] | None,
    position_hint_args: tuple[float, float, float] | None,
    parity: str,
    workers: int,
) -> dict | None:
    """Run plate solving across multiple processes with round-robin file distribution.

    Index files are distributed round-robin so each worker gets a mix of
    scales and sky regions. Each worker processes all star slices on its
    subset of files.

    Returns the result dict from the first worker to find a match, or None.
    """
    # Round-robin distribution: worker i gets files i, i+N, i+2N, ...
    file_chunks: list[list[pathlib.Path]] = [[] for _ in range(workers)]
    for i, f in enumerate(index_files):
        file_chunks[i % workers].append(f)

    result_queue = multiprocessing.Queue()
    cancel_event = multiprocessing.Event()

    processes = []
    for i in range(workers):
        if not file_chunks[i]:
            continue
        p = multiprocessing.Process(
            target=_solve_worker,
            args=(
                file_chunks[i], stars, size_hint_args, position_hint_args,
                parity, i, result_queue, cancel_event,
            ),
        )
        processes.append(p)

    for p in processes:
        p.start()

    # Wait for first result or all processes to finish
    result = None
    alive = set(range(len(processes)))
    while alive:
        try:
            result = result_queue.get(timeout=0.5)
            break
        except Exception:
            pass
        for i in list(alive):
            if not processes[i].is_alive():
                alive.discard(i)

    # Signal remaining workers to stop and terminate them
    cancel_event.set()
    for p in processes:
        p.join(timeout=2.0)
        if p.is_alive():
            p.terminate()
            p.join(timeout=2.0)
            if p.is_alive():
                p.kill()
                p.join(timeout=1.0)

    return result


def solve(
    image_path: str | pathlib.Path,
    index_dir: pathlib.Path,
    scale_low: float | None = None,
    scale_high: float | None = None,
    ra_hint: float | None = None,
    dec_hint: float | None = None,
    radius_hint: float | None = None,
    max_stars: int = 1000,
    verbose: bool = False,
    parity: str = "both",
    use_exif: bool = True,
    workers: int = 1,
) -> tuple[WCS, dict]:
    """Plate-solve an image.

    Args:
        image_path: Path to the image file.
        index_dir: Directory containing astrometry.net index files.
        scale_low: Lower bound of plate scale in arcsec/pixel.
        scale_high: Upper bound of plate scale in arcsec/pixel.
        ra_hint: RA hint in degrees.
        dec_hint: Dec hint in degrees.
        radius_hint: Search radius in degrees (requires ra_hint and dec_hint).
        max_stars: Maximum number of stars to detect.
        verbose: Print progress info.
        parity: Image parity - "both", "normal", or "flip".
        use_exif: Try to auto-detect plate scale from EXIF data.
        workers: Number of parallel solver processes (1 = single-process).

    Returns:
        Tuple of (WCS object, metadata dict with center_ra, center_dec, scale, logodds).

    Raises:
        FileNotFoundError: If image or index files not found.
        RuntimeError: If solving fails.
    """
    import astrometry

    t_total = time.monotonic()
    image_path = pathlib.Path(image_path)
    if not image_path.exists():
        raise FileNotFoundError(f"Image not found: {image_path}")

    # Collect index files
    index_files = sorted(index_dir.rglob("*.fits"))
    if not index_files:
        raise FileNotFoundError(
            f"No index files found in {index_dir}. Run 'my-astrometry setup' first."
        )

    print(f"  Loading image ({image_path.name})...")
    image = load_image(image_path)
    height, width = image.shape
    print(f"  Image size: {width} x {height}")

    print(f"  Detecting stars...")
    t0 = time.monotonic()
    stars = detect_stars(image, max_stars=max_stars)
    if len(stars) == 0:
        raise RuntimeError("No stars detected in the image.")
    print(f"  Found {len(stars)} stars ({time.monotonic() - t0:.1f}s)")

    # Auto-detect plate scale from EXIF if no explicit scale hints given
    if use_exif and scale_low is None and scale_high is None:
        exif_scale = estimate_plate_scale_from_exif(image_path)
        if exif_scale is not None:
            scale_low, scale_high = exif_scale
            print(f"  EXIF auto-detected scale: {scale_low:.2f}-{scale_high:.2f} arcsec/px")

    # Build hints (as plain tuples for pickling to worker processes)
    size_hint_args = None
    if scale_low is not None or scale_high is not None:
        size_hint_args = (scale_low or 0.1, scale_high or 1000.0)

    position_hint_args = None
    if ra_hint is not None and dec_hint is not None:
        position_hint_args = (ra_hint, dec_hint, radius_hint or 10.0)

    if size_hint_args or position_hint_args:
        hints = []
        if size_hint_args:
            hints.append(f"scale {size_hint_args[0]:.2f}-{size_hint_args[1]:.2f}\"/px")
        if position_hint_args:
            hints.append(f"pos ({position_hint_args[0]:.1f}, {position_hint_args[1]:.1f}) r={position_hint_args[2]:.1f}")
        print(f"  Hints: {', '.join(hints)}")

    # Pre-filter index files to matching scales, then sort by sky region priority
    all_count = len(index_files)
    index_files = _filter_index_files(index_files, width, height,
                                      size_hint_args[0] if size_hint_args else None,
                                      size_hint_args[1] if size_hint_args else None)
    index_files = _sort_index_files_by_priority(index_files)
    if len(index_files) < all_count:
        parsed = [_parse_index_file(f) for f in index_files]
        series_counts = {}
        for p in parsed:
            if p:
                series_counts[p[0]] = series_counts.get(p[0], 0) + 1
        parts = [f"{count} {series}" for series, count in sorted(series_counts.items())]
        print(f"  Using {len(index_files)}/{all_count} index files ({', '.join(parts)})")

    if verbose:
        import logging
        logging.getLogger().setLevel(logging.INFO)

    effective_workers = min(workers, len(index_files))

    if effective_workers > 1:
        print(f"  Plate solving ({len(index_files)} index files, {effective_workers} workers)...")
        t0 = time.monotonic()

        result = _solve_parallel(
            index_files, stars, size_hint_args, position_hint_args,
            parity, effective_workers,
        )
        solve_time = time.monotonic() - t0

        if result is None:
            print(f"  Failed after {solve_time:.1f}s")
            raise RuntimeError(
                "Plate solving failed — no match found. "
                "Try providing scale hints (--scale-low, --scale-high) "
                "or position hints (--ra, --dec, --radius)."
            )

        from astropy.io.fits import Header
        header = Header()
        for k, v in result["wcs_header"].items():
            header[k] = v
        wcs = WCS(header, relax=True)
        total_time = time.monotonic() - t_total
        metadata = {
            "center_ra_deg": result["center_ra_deg"],
            "center_dec_deg": result["center_dec_deg"],
            "scale_arcsec_per_pixel": result["scale_arcsec_per_pixel"],
            "logodds": result["logodds"],
            "image_width": width,
            "image_height": height,
            "num_stars_detected": len(stars),
            "solve_time_s": solve_time,
            "total_time_s": total_time,
        }
        print(f"  Solved in {solve_time:.1f}s (total {total_time:.1f}s, worker {result['worker_id']})")
        return wcs, metadata

    # Single-process path (workers=1)
    parity_map = {
        "both": astrometry.Parity.BOTH,
        "normal": astrometry.Parity.NORMAL,
        "flip": astrometry.Parity.FLIP,
    }
    solution_parameters = astrometry.SolutionParameters(
        parity=parity_map[parity],
        logodds_callback=lambda logodds_list: (
            astrometry.Action.STOP
            if logodds_list[0] > 20.0
            else astrometry.Action.CONTINUE
        ),
    )

    print(f"  Plate solving ({len(index_files)} index files)...")
    t0 = time.monotonic()

    with astrometry.Solver(index_files) as solver:
        solution = solver.solve(
            stars=stars,
            size_hint=size_hint_args and astrometry.SizeHint(
                lower_arcsec_per_pixel=size_hint_args[0],
                upper_arcsec_per_pixel=size_hint_args[1],
            ),
            position_hint=position_hint_args and astrometry.PositionHint(
                ra_deg=position_hint_args[0],
                dec_deg=position_hint_args[1],
                radius_deg=position_hint_args[2],
            ),
            solution_parameters=solution_parameters,
        )

    solve_time = time.monotonic() - t0

    if not solution.has_match():
        print(f"  Failed after {solve_time:.1f}s")
        raise RuntimeError(
            "Plate solving failed — no match found. "
            "Try providing scale hints (--scale-low, --scale-high) "
            "or position hints (--ra, --dec, --radius)."
        )

    match = solution.best_match()
    wcs = match.astropy_wcs()

    total_time = time.monotonic() - t_total

    metadata = {
        "center_ra_deg": match.center_ra_deg,
        "center_dec_deg": match.center_dec_deg,
        "scale_arcsec_per_pixel": match.scale_arcsec_per_pixel,
        "logodds": match.logodds,
        "image_width": width,
        "image_height": height,
        "num_stars_detected": len(stars),
        "solve_time_s": solve_time,
        "total_time_s": total_time,
    }

    print(f"  Solved in {solve_time:.1f}s (total {total_time:.1f}s)")

    return wcs, metadata
