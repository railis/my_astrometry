"""Load and filter catalog data (OpenNGC, Messier, constellations)."""

import csv
import json
import pathlib
import warnings

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

DATA_DIR = pathlib.Path(__file__).parent / "data"

# Common name aliases → Messier number (lowercase)
_COMMON_ALIASES = {
    "crab": "M1", "crab nebula": "M1",
    "andromeda": "M31", "andromeda galaxy": "M31",
    "triangulum": "M33", "triangulum galaxy": "M33",
    "orion": "M42", "orion nebula": "M42",
    "pleiades": "M45",
    "whirlpool": "M51", "whirlpool galaxy": "M51",
    "ring": "M57", "ring nebula": "M57",
    "dumbbell": "M27", "dumbbell nebula": "M27",
    "sombrero": "M104", "sombrero galaxy": "M104",
    "eagle": "M16", "eagle nebula": "M16",
    "lagoon": "M8", "lagoon nebula": "M8",
    "trifid": "M20", "trifid nebula": "M20",
    "omega": "M17", "omega nebula": "M17",
    "beehive": "M44", "beehive cluster": "M44",
    "hercules": "M13", "hercules cluster": "M13", "great cluster": "M13",
    "pinwheel": "M101", "pinwheel galaxy": "M101",
    "sunflower": "M63", "sunflower galaxy": "M63",
    "black eye": "M64", "black eye galaxy": "M64",
    "bode": "M81", "bode's galaxy": "M81",
    "cigar": "M82", "cigar galaxy": "M82",
    "owl": "M97", "owl nebula": "M97",
    "wild duck": "M11", "wild duck cluster": "M11",
}


def lookup_object(
    name: str,
    openngc_path: pathlib.Path | None = None,
) -> tuple[float, float, str] | None:
    """Look up an object by name and return (ra_deg, dec_deg, canonical_name).

    Accepts formats like: M42, m42, NGC1976, ngc 1976, IC434, ic 434,
    or common names like "orion nebula", "crab".

    Returns None if not found.
    """
    raw = name.strip()
    normalized = raw.lower().replace(" ", "")

    # Check common name aliases first
    alias_match = _COMMON_ALIASES.get(raw.lower().strip())
    if alias_match:
        normalized = alias_match.lower()

    # Try Messier catalog (bundled, always available)
    if normalized.startswith("m") and normalized[1:].isdigit():
        canon = f"M{int(normalized[1:])}"
        for obj in load_messier_catalog():
            if obj["name"] == canon:
                display = canon
                if obj.get("common_name"):
                    display += f" ({obj['common_name']})"
                return obj["ra_deg"], obj["dec_deg"], display

    # Try NGC/IC from OpenNGC
    if openngc_path and openngc_path.exists():
        target = None
        if normalized.startswith("ngc") and normalized[3:].isdigit():
            target = f"NGC{int(normalized[3:]):04d}"
        elif normalized.startswith("ic") and normalized[2:].isdigit():
            target = f"IC{int(normalized[2:]):04d}"

        if target:
            for obj in load_openngc(openngc_path):
                if obj["name"] == target:
                    display = target
                    if obj.get("messier"):
                        display = f"{obj['messier']} ({target})"
                    return obj["ra_deg"], obj["dec_deg"], display

    return None


def _parse_hms_to_deg(hms_str: str) -> float | None:
    """Parse 'HH:MM:SS.s' to degrees."""
    if not hms_str or hms_str.strip() == "":
        return None
    parts = hms_str.strip().split(":")
    if len(parts) != 3:
        return None
    try:
        h, m, s = float(parts[0]), float(parts[1]), float(parts[2])
        return (h + m / 60 + s / 3600) * 15.0
    except ValueError:
        return None


def _parse_dms_to_deg(dms_str: str) -> float | None:
    """Parse '+DD:MM:SS.s' to degrees."""
    if not dms_str or dms_str.strip() == "":
        return None
    s = dms_str.strip()
    sign = -1 if s.startswith("-") else 1
    s = s.lstrip("+-")
    parts = s.split(":")
    if len(parts) != 3:
        return None
    try:
        d, m, sec = float(parts[0]), float(parts[1]), float(parts[2])
        return sign * (d + m / 60 + sec / 3600)
    except ValueError:
        return None


def load_messier_catalog() -> list[dict]:
    """Load the bundled Messier catalog."""
    path = DATA_DIR / "messier.csv"
    objects = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mag = None
            if row.get("magnitude"):
                try:
                    mag = float(row["magnitude"])
                except ValueError:
                    pass
            size = None
            if row.get("major_axis_arcmin"):
                try:
                    size = float(row["major_axis_arcmin"])
                except ValueError:
                    pass
            objects.append({
                "name": row["name"],
                "ra_deg": float(row["ra_deg"]),
                "dec_deg": float(row["dec_deg"]),
                "type": row["type"],
                "common_name": row.get("common_name", ""),
                "magnitude": mag,
                "major_axis_arcmin": size,
            })
    return objects


def load_openngc(openngc_path: pathlib.Path) -> list[dict]:
    """Load the OpenNGC catalog CSV.

    Returns a list of dicts with name, ra_deg, dec_deg, type, messier.
    """
    objects = []
    with open(openngc_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            name = row.get("Name", "").strip()
            if not name:
                continue

            ra_str = row.get("RA", "")
            dec_str = row.get("Dec", "")
            ra = _parse_hms_to_deg(ra_str)
            dec = _parse_dms_to_deg(dec_str)
            if ra is None or dec is None:
                continue

            obj_type = row.get("Type", "").strip()
            messier = row.get("M", "").strip()
            mag_str = row.get("V-Mag", "") or row.get("B-Mag", "")
            mag = None
            if mag_str:
                try:
                    mag = float(mag_str)
                except ValueError:
                    pass

            size = None
            size_str = row.get("MajAx", "").strip()
            if size_str:
                try:
                    size = float(size_str)
                except ValueError:
                    pass

            objects.append({
                "name": name,
                "ra_deg": ra,
                "dec_deg": dec,
                "type": obj_type,
                "messier": f"M{messier}" if messier else "",
                "magnitude": mag,
                "major_axis_arcmin": size,
            })
    return objects


def filter_objects_in_fov(
    objects: list[dict],
    wcs: WCS,
    width: int,
    height: int,
    margin: float = 0.05,
) -> list[dict]:
    """Filter catalog objects to those visible in the image FOV.

    Args:
        objects: List of objects with ra_deg and dec_deg.
        wcs: WCS solution.
        width: Image width in pixels.
        height: Image height in pixels.
        margin: Fraction of image size to add as margin (for labels).

    Returns:
        List of objects with added 'x' and 'y' pixel coordinates.
    """
    if not objects:
        return []

    ras = np.array([o["ra_deg"] for o in objects])
    decs = np.array([o["dec_deg"] for o in objects])
    coords = SkyCoord(ra=ras * u.deg, dec=decs * u.deg)

    # Convert to pixel coordinates (suppress convergence warnings for
    # objects far from the image center — results will be NaN, filtered below)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        pixels = wcs.world_to_pixel(coords)
    xs, ys = pixels[0], pixels[1]

    # Filter to those within image bounds (with margin)
    mx = width * margin
    my = height * margin
    visible = []
    for i, obj in enumerate(objects):
        x, y = float(xs[i]), float(ys[i])
        if np.isnan(x) or np.isnan(y):
            continue
        if -mx <= x <= width + mx and -my <= y <= height + my:
            result = dict(obj)
            result["x"] = x
            result["y"] = y
            visible.append(result)

    return visible


def get_catalog_objects(
    wcs: WCS,
    width: int,
    height: int,
    openngc_path: pathlib.Path | None = None,
    mag_limit: float | None = None,
    include_ngc: bool = True,
) -> dict[str, list[dict]]:
    """Get all catalog objects visible in the FOV, categorized.

    Returns:
        Dict with keys 'messier', 'ic', 'ngc' each containing a list of objects.
    """
    result = {"messier": [], "ic": [], "ngc": []}

    # Always load bundled Messier catalog (more reliable names)
    messier_objects = load_messier_catalog()
    result["messier"] = filter_objects_in_fov(messier_objects, wcs, width, height)

    # Load OpenNGC for IC/NGC objects if available
    if openngc_path and openngc_path.exists():
        all_ngc = load_openngc(openngc_path)

        # Build set of Messier names to avoid duplicates
        messier_names = set()
        for obj in all_ngc:
            if obj["messier"]:
                messier_names.add(obj["name"])

        # Separate IC and NGC (non-Messier)
        ic_objects = []
        ngc_objects = []
        for obj in all_ngc:
            if obj["name"] in messier_names:
                continue  # Already in Messier list
            if mag_limit is not None and obj["magnitude"] is not None:
                if obj["magnitude"] > mag_limit:
                    continue
            if obj["name"].startswith("IC"):
                ic_objects.append(obj)
            elif include_ngc:
                ngc_objects.append(obj)

        result["ic"] = filter_objects_in_fov(ic_objects, wcs, width, height)
        result["ngc"] = filter_objects_in_fov(ngc_objects, wcs, width, height)

    return result


def load_constellation_lines() -> dict[str, list[list[list[float]]]]:
    """Load constellation stick figure data.

    Returns:
        Dict mapping IAU abbreviation to list of line segments.
        Each segment is [[ra1_deg, dec1_deg], [ra2_deg, dec2_deg]].
    """
    path = DATA_DIR / "constellations.json"
    with open(path) as f:
        return json.load(f)


def load_constellation_names() -> dict[str, str]:
    """Load constellation IAU abbreviation to full name mapping."""
    path = DATA_DIR / "constellation_names.json"
    with open(path) as f:
        return json.load(f)


def get_constellation_lines_in_fov(
    wcs: WCS, width: int, height: int
) -> list[dict]:
    """Get constellation line segments near the FOV in sky coordinates.

    Returns RA/Dec coordinates (not pixel coords).  The annotator uses WCSAxes
    to project and clip these onto the image automatically.

    Returns:
        List of dicts with 'abbr', 'name', and 'segments' (RA/Dec in degrees).
    """
    all_lines = load_constellation_lines()
    names = load_constellation_names()

    # Compute field center and angular radius for pre-filtering
    center_sky = wcs.pixel_to_world(width / 2, height / 2)
    corner_sky = wcs.pixel_to_world(0, 0)
    fov_radius_deg = center_sky.separation(corner_sky).deg
    max_angular_dist = fov_radius_deg * 2.0

    visible_constellations = []

    for abbr, segments in all_lines.items():
        visible_segments = []
        for segment in segments:
            ra1, dec1 = segment[0]
            ra2, dec2 = segment[1]

            coord1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg)
            coord2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg)

            # Only include segments where both endpoints are near the FOV
            if (center_sky.separation(coord1).deg > max_angular_dist
                    or center_sky.separation(coord2).deg > max_angular_dist):
                continue

            visible_segments.append({
                "ra1": ra1, "dec1": dec1,
                "ra2": ra2, "dec2": dec2,
            })

        if visible_segments:
            visible_constellations.append({
                "abbr": abbr,
                "name": names.get(abbr, abbr),
                "segments": visible_segments,
            })

    return visible_constellations
