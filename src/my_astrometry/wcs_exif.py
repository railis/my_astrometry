"""Serialize/deserialize WCS solutions to/from image EXIF UserComment."""

import json
import pathlib

from astropy.io.fits import Header
from astropy.wcs import WCS
from PIL import Image

# EXIF sub-IFD tag number (ExifOffset / ExifIFDPointer)
_EXIF_IFD = 0x8769
# UserComment tag in the EXIF sub-IFD
_USER_COMMENT_TAG = 0x9286
# Prefix to identify our data in UserComment
_MAGIC_PREFIX = "my_astrometry_wcs:"
# Current serialization format version
_FORMAT_VERSION = "0.1.0"


def save_solution(
    image_path: pathlib.Path,
    wcs: WCS,
    metadata: dict,
) -> None:
    """Save a WCS solution and solve metadata into the image's EXIF UserComment.

    For JPEG files, re-saves with quality='keep' and subsampling='keep' to
    preserve quality.  Preserves all existing EXIF data.

    Args:
        image_path: Path to the image file (JPEG, PNG, or TIFF).
        wcs: Astropy WCS object from plate solving.
        metadata: Solve metadata dict (center coords, scale, logodds, etc.).
    """
    image_path = pathlib.Path(image_path)
    img = Image.open(image_path)

    # Build the JSON payload
    wcs_header = {}
    for card in wcs.to_header().cards:
        if card.keyword:
            wcs_header[card.keyword] = card.value

    payload = {
        "version": _FORMAT_VERSION,
        "wcs_header": wcs_header,
        "metadata": metadata,
    }
    json_str = _MAGIC_PREFIX + json.dumps(payload, separators=(",", ":"))

    # Read existing EXIF and inject UserComment
    exif = img.getexif()
    exif_ifd = exif.get_ifd(_EXIF_IFD)
    exif_ifd[_USER_COMMENT_TAG] = json_str

    # Re-save with EXIF, preserving quality for JPEG
    save_kwargs = {"exif": exif.tobytes()}
    suffix = image_path.suffix.lower()
    if suffix in (".jpg", ".jpeg"):
        save_kwargs["quality"] = "keep"
        save_kwargs["subsampling"] = "keep"

    img.save(image_path, **save_kwargs)


def load_solution(
    image_path: pathlib.Path,
) -> tuple[WCS, dict] | None:
    """Load a WCS solution and metadata from the image's EXIF UserComment.

    Returns:
        Tuple of (WCS, metadata_dict), or None if no solution is stored.
    """
    image_path = pathlib.Path(image_path)
    try:
        img = Image.open(image_path)
        exif = img.getexif()
        exif_ifd = exif.get_ifd(_EXIF_IFD)
    except Exception:
        return None

    user_comment = exif_ifd.get(_USER_COMMENT_TAG)
    if not user_comment or not isinstance(user_comment, str):
        return None

    if not user_comment.startswith(_MAGIC_PREFIX):
        return None

    json_str = user_comment[len(_MAGIC_PREFIX):]
    try:
        payload = json.loads(json_str)
    except json.JSONDecodeError:
        return None

    wcs_header_dict = payload.get("wcs_header")
    metadata = payload.get("metadata", {})
    if wcs_header_dict is None:
        return None

    header = Header(wcs_header_dict)
    wcs = WCS(header)
    return wcs, metadata


def has_solution(image_path: pathlib.Path) -> bool:
    """Check if an image has a WCS solution stored in its EXIF.

    Lightweight check that avoids full WCS reconstruction.
    """
    image_path = pathlib.Path(image_path)
    try:
        img = Image.open(image_path)
        exif = img.getexif()
        exif_ifd = exif.get_ifd(_EXIF_IFD)
        user_comment = exif_ifd.get(_USER_COMMENT_TAG)
        return (
            isinstance(user_comment, str)
            and user_comment.startswith(_MAGIC_PREFIX)
        )
    except Exception:
        return False
