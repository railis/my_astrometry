"""Star detection and plate solving pipeline."""

import math
import pathlib
import re
import time

import numpy as np
from astropy.wcs import WCS
from PIL import Image


def _index_scale_number(path: pathlib.Path) -> int | None:
    """Extract the scale number from an index filename like index-4112.fits."""
    m = re.search(r"index-\d\d(\d\d)\.fits", path.name)
    if m:
        return int(m.group(1))
    return None


def _filter_index_files(
    index_files: list[pathlib.Path],
    width: int,
    height: int,
    scale_low: float | None,
    scale_high: float | None,
) -> list[pathlib.Path]:
    """Pre-filter index files to only those relevant for the given scale range.

    Each 4100-series index at scale N covers quads of angular size ~30*sqrt(2)^(N-10) arcmin.
    We keep only index files whose quad scale overlaps with the expected image FOV.
    """
    if scale_low is None and scale_high is None:
        return index_files

    lo = scale_low or 0.1
    hi = scale_high or 1000.0

    # Image FOV range in arcmin for the given plate scale bounds
    diag_px = math.hypot(width, height)
    short_px = min(width, height)

    # Quads should be between ~5% and ~100% of the image shortest dimension
    min_quad_arcmin = short_px * lo / 60.0 * 0.05
    max_quad_arcmin = diag_px * hi / 60.0

    filtered = []
    for f in index_files:
        scale_n = _index_scale_number(f)
        if scale_n is None:
            # Can't determine scale — keep it to be safe
            filtered.append(f)
            continue
        # Approximate quad angular size for this scale
        quad_arcmin = 30.0 * (2 ** ((scale_n - 10) / 2.0))
        if min_quad_arcmin <= quad_arcmin <= max_quad_arcmin:
            filtered.append(f)

    # If filtering removed everything, fall back to all files
    if not filtered:
        return index_files
    return filtered


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


def solve(
    image_path: str | pathlib.Path,
    index_dir: pathlib.Path,
    scale_low: float | None = None,
    scale_high: float | None = None,
    ra_hint: float | None = None,
    dec_hint: float | None = None,
    radius_hint: float | None = None,
    max_stars: int = 200,
    verbose: bool = False,
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

    # Build hints
    size_hint = None
    if scale_low is not None or scale_high is not None:
        size_hint = astrometry.SizeHint(
            lower_arcsec_per_pixel=scale_low or 0.1,
            upper_arcsec_per_pixel=scale_high or 1000.0,
        )

    position_hint = None
    if ra_hint is not None and dec_hint is not None:
        position_hint = astrometry.PositionHint(
            ra_deg=ra_hint,
            dec_deg=dec_hint,
            radius_deg=radius_hint or 10.0,
        )

    if size_hint or position_hint:
        hints = []
        if size_hint:
            hints.append(f"scale {size_hint.lower_arcsec_per_pixel}-{size_hint.upper_arcsec_per_pixel}\"/px")
        if position_hint:
            hints.append(f"pos ({position_hint.ra_deg:.1f}, {position_hint.dec_deg:.1f}) r={position_hint.radius_deg:.1f}")
        print(f"  Hints: {', '.join(hints)}")

    # Pre-filter index files to matching scales
    all_count = len(index_files)
    index_files = _filter_index_files(index_files, width, height, scale_low, scale_high)
    if len(index_files) < all_count:
        scales = sorted(s for f in index_files if (s := _index_scale_number(f)) is not None)
        print(f"  Using {len(index_files)}/{all_count} index files (scales {min(scales)}-{max(scales)})")

    if verbose:
        import logging

        logging.getLogger().setLevel(logging.INFO)

    solution_parameters = astrometry.SolutionParameters(
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
            size_hint=size_hint,
            position_hint=position_hint,
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
