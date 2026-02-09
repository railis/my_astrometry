"""Extract RA/Dec coordinates, FOV, and plate scale from a WCS solution."""

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


def _format_ra_hms(coord: SkyCoord) -> str:
    """Format RA as HH:MM:SS.ss."""
    return coord.ra.to_string(unit="hour", sep=":", precision=2)


def _format_dec_dms(coord: SkyCoord) -> str:
    """Format Dec as +DD:MM:SS.s."""
    return coord.dec.to_string(unit="degree", sep=":", precision=1, alwayssign=True)


def get_field_info(wcs: WCS, width: int, height: int) -> dict:
    """Extract full field information from a WCS solution.

    Args:
        wcs: Astropy WCS object from plate solving.
        width: Image width in pixels.
        height: Image height in pixels.

    Returns:
        Dict with center, corners, FOV, plate scale info.
    """
    # Center of the image
    center = wcs.pixel_to_world(width / 2, height / 2)

    # Four corners (0-indexed pixel coordinates)
    top_left = wcs.pixel_to_world(0, 0)
    top_right = wcs.pixel_to_world(width - 1, 0)
    bottom_left = wcs.pixel_to_world(0, height - 1)
    bottom_right = wcs.pixel_to_world(width - 1, height - 1)

    # FOV: angular separation across the diagonal and along each axis
    fov_horizontal = top_left.separation(top_right)
    fov_vertical = top_left.separation(bottom_left)
    fov_diagonal = top_left.separation(bottom_right)

    # Plate scale from pixel scale
    # Use two adjacent pixels near center
    p1 = wcs.pixel_to_world(width / 2, height / 2)
    p2 = wcs.pixel_to_world(width / 2 + 1, height / 2)
    plate_scale = p1.separation(p2)

    return {
        "center": {
            "ra_deg": center.ra.deg,
            "dec_deg": center.dec.deg,
            "ra_hms": _format_ra_hms(center),
            "dec_dms": _format_dec_dms(center),
        },
        "corners": {
            "top_left": {"ra_deg": top_left.ra.deg, "dec_deg": top_left.dec.deg},
            "top_right": {"ra_deg": top_right.ra.deg, "dec_deg": top_right.dec.deg},
            "bottom_left": {"ra_deg": bottom_left.ra.deg, "dec_deg": bottom_left.dec.deg},
            "bottom_right": {"ra_deg": bottom_right.ra.deg, "dec_deg": bottom_right.dec.deg},
        },
        "fov": {
            "horizontal_deg": fov_horizontal.deg,
            "vertical_deg": fov_vertical.deg,
            "diagonal_deg": fov_diagonal.deg,
            "horizontal_arcmin": fov_horizontal.arcmin,
            "vertical_arcmin": fov_vertical.arcmin,
        },
        "plate_scale_arcsec_per_pixel": plate_scale.arcsec,
    }


def print_field_info(info: dict) -> None:
    """Print field information in a readable format."""
    c = info["center"]
    fov = info["fov"]
    ps = info["plate_scale_arcsec_per_pixel"]

    print(f"  Center:      RA {c['ra_hms']}  Dec {c['dec_dms']}")
    print(f"               ({c['ra_deg']:.5f}, {c['dec_deg']:.5f}) deg")
    print(f"  FOV:         {fov['horizontal_arcmin']:.1f}' x {fov['vertical_arcmin']:.1f}'")
    print(f"               ({fov['horizontal_deg']:.3f} x {fov['vertical_deg']:.3f} deg)")
    print(f"  Plate scale: {ps:.3f} arcsec/pixel")
