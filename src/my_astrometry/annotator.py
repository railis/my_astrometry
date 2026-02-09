"""Overlay rendering: constellation lines, RA/Dec grid, Messier/IC/NGC labels on solved image."""

import pathlib
import warnings

import matplotlib
matplotlib.use("Agg")

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from PIL import Image


def create_annotated_image(
    image_path: str | pathlib.Path,
    wcs: WCS,
    catalog_objects: dict[str, list[dict]],
    constellation_lines: list[dict],
    output_path: str | pathlib.Path,
    plate_scale_arcsec: float = 1.0,
    show_ngc: bool = True,
    verbose: bool = False,
) -> pathlib.Path:
    """Render an annotated overlay on the original image.

    Constellation lines and the RA/Dec grid are projected to pixel space
    via wcs.world_to_pixel() so that SIP distortion is handled correctly.
    """
    output_path = pathlib.Path(output_path)

    # Load original image
    img = Image.open(image_path)
    img_array = np.array(img)
    height, width = img_array.shape[:2]

    # Create figure with plain axes (origin="upper" so y=0 is top)
    dpi = 100
    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    ax.imshow(img_array, origin="upper", aspect="equal")
    ax.set_xlim(-0.5, width - 0.5)
    ax.set_ylim(height - 0.5, -0.5)
    ax.axis("off")

    # Base scaling factor for annotations based on image size
    img_scale = min(width, height) / 1000.0
    base_font = max(6, min(14, 9 * img_scale))
    base_marker = max(8, min(30, 15 * img_scale))
    line_width = max(0.5, min(2.0, 1.0 * img_scale))

    # --- RA/Dec coordinate grid ---
    _draw_grid(ax, wcs, width, height, base_font)

    # --- Constellation lines (project RA/Dec to pixel via WCS) ---
    _draw_constellation_lines(ax, wcs, constellation_lines, width, height,
                              line_width, base_font)

    # --- Catalog objects (already in pixel coordinates) ---
    _draw_objects(
        ax, catalog_objects.get("messier", []),
        color="#ffdd00", base_marker=base_marker, base_font=base_font,
        plate_scale_arcsec=plate_scale_arcsec,
        label_key="name", common_name_key="common_name",
    )

    _draw_objects(
        ax, catalog_objects.get("ic", []),
        color="#00dddd", base_marker=base_marker, base_font=base_font,
        plate_scale_arcsec=plate_scale_arcsec,
        label_key="name",
    )

    if show_ngc:
        _draw_objects(
            ax, catalog_objects.get("ngc", []),
            color="#44cc44", base_marker=base_marker, base_font=base_font,
            plate_scale_arcsec=plate_scale_arcsec,
            label_key="name",
        )

    # Save at exact image resolution
    fig.savefig(output_path, dpi=dpi, facecolor="black")
    plt.close(fig)

    return output_path


# ---------------------------------------------------------------------------
#  Line clipping (Liang-Barsky algorithm)
# ---------------------------------------------------------------------------

def _clip_segment(x1, y1, x2, y2, xmin, ymin, xmax, ymax):
    """Clip a line segment to a rectangle using Liang-Barsky algorithm.

    Returns (cx1, cy1, cx2, cy2) or None if the segment is entirely outside.
    """
    dx = x2 - x1
    dy = y2 - y1
    t_min, t_max = 0.0, 1.0

    for p, q in [(-dx, x1 - xmin), (dx, xmax - x1),
                 (-dy, y1 - ymin), (dy, ymax - y1)]:
        if p == 0:
            if q < 0:
                return None  # parallel and outside
        else:
            t = q / p
            if p < 0:
                t_min = max(t_min, t)
            else:
                t_max = min(t_max, t)

    if t_min > t_max:
        return None

    cx1 = x1 + t_min * dx
    cy1 = y1 + t_min * dy
    cx2 = x1 + t_max * dx
    cy2 = y1 + t_max * dy
    return cx1, cy1, cx2, cy2


# ---------------------------------------------------------------------------
#  RA/Dec grid
# ---------------------------------------------------------------------------

def _nice_grid_spacing(fov_deg: float) -> float:
    """Pick an aesthetically pleasing grid spacing for a given FOV."""
    candidates = [0.1, 0.25, 0.5, 1, 2, 5, 10, 15, 30, 45]
    for s in candidates:
        if fov_deg / s <= 10:
            return s
    return 45.0


def _draw_grid(ax, wcs: WCS, width: int, height: int, base_font: float) -> None:
    """Draw an RA/Dec coordinate grid by projecting grid lines to pixel space."""
    corners_x = [0, width - 1, 0, width - 1]
    corners_y = [0, 0, height - 1, height - 1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        corner_sky = wcs.pixel_to_world(corners_x, corners_y)

    ra_vals = corner_sky.ra.deg
    dec_vals = corner_sky.dec.deg

    # Handle RA wrap-around
    ra_range = ra_vals.max() - ra_vals.min()
    if ra_range > 180:
        ra_shifted = np.where(ra_vals > 180, ra_vals - 360, ra_vals)
        ra_min, ra_max = ra_shifted.min(), ra_shifted.max()
    else:
        ra_min, ra_max = ra_vals.min(), ra_vals.max()

    dec_min, dec_max = dec_vals.min(), dec_vals.max()

    # Margin so grid extends slightly beyond the visible corners
    ra_margin = (ra_max - ra_min) * 0.15
    dec_margin = (dec_max - dec_min) * 0.15
    ra_min -= ra_margin
    ra_max += ra_margin
    dec_min = max(-90, dec_min - dec_margin)
    dec_max = min(90, dec_max + dec_margin)

    ra_spacing = _nice_grid_spacing(ra_max - ra_min)
    dec_spacing = _nice_grid_spacing(dec_max - dec_min)

    ra_start = np.floor(ra_min / ra_spacing) * ra_spacing
    ra_end = np.ceil(ra_max / ra_spacing) * ra_spacing
    dec_start = np.floor(dec_min / dec_spacing) * dec_spacing
    dec_end = np.ceil(dec_max / dec_spacing) * dec_spacing

    n_samples = 200
    grid_color = "white"
    grid_alpha = 0.3
    grid_lw = 0.5
    label_size = max(5, base_font * 0.65)

    # Lines of constant RA
    for ra in np.arange(ra_start, ra_end + ra_spacing / 2, ra_spacing):
        ra_norm = ra % 360
        decs = np.linspace(dec_min, dec_max, n_samples)
        coords = SkyCoord(ra=np.full_like(decs, ra_norm) * u.deg, dec=decs * u.deg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            px, py = wcs.world_to_pixel(coords)
        mask = np.isfinite(px) & np.isfinite(py)
        mask &= (px >= -width * 0.05) & (px <= width * 1.05)
        mask &= (py >= -height * 0.05) & (py <= height * 1.05)
        if mask.sum() < 2:
            continue
        _plot_clipped(ax, px, py, mask, color=grid_color, alpha=grid_alpha,
                      linewidth=grid_lw, linestyle="--")
        _place_edge_label(ax, px, py, mask, width, height,
                          _format_ra(ra_norm), label_size, edge="bottom")

    # Lines of constant Dec
    for dec in np.arange(dec_start, dec_end + dec_spacing / 2, dec_spacing):
        dec = np.clip(dec, -90, 90)
        ras = np.linspace(ra_min, ra_max, n_samples) % 360
        coords = SkyCoord(ra=ras * u.deg, dec=np.full_like(ras, dec) * u.deg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            px, py = wcs.world_to_pixel(coords)
        mask = np.isfinite(px) & np.isfinite(py)
        mask &= (px >= -width * 0.05) & (px <= width * 1.05)
        mask &= (py >= -height * 0.05) & (py <= height * 1.05)
        if mask.sum() < 2:
            continue
        _plot_clipped(ax, px, py, mask, color=grid_color, alpha=grid_alpha,
                      linewidth=grid_lw, linestyle="--")
        _place_edge_label(ax, px, py, mask, width, height,
                          _format_dec(dec), label_size, edge="left")


def _plot_clipped(ax, px, py, mask, **kwargs) -> None:
    """Plot line segments, breaking at gaps in the mask."""
    indices = np.where(mask)[0]
    if len(indices) == 0:
        return
    breaks = np.where(np.diff(indices) > 1)[0] + 1
    segments = np.split(indices, breaks)
    for seg_idx in segments:
        if len(seg_idx) >= 2:
            ax.plot(px[seg_idx], py[seg_idx], **kwargs)


def _place_edge_label(ax, px, py, mask, width, height,
                      label: str, fontsize: float, edge: str) -> None:
    """Place a coordinate label where a grid line crosses an image edge."""
    valid_x, valid_y = px[mask], py[mask]
    if edge == "bottom":
        idx = np.argmax(valid_y)
        x, y = float(valid_x[idx]), float(valid_y[idx])
        if 0 <= x <= width and abs(y - height) < height * 0.15:
            ax.text(x, height - 2, label, color="white", alpha=0.5,
                    fontsize=fontsize, ha="center", va="bottom")
    elif edge == "left":
        idx = np.argmin(valid_x)
        x, y = float(valid_x[idx]), float(valid_y[idx])
        if abs(x) < width * 0.15 and 0 <= y <= height:
            ax.text(2, y, label, color="white", alpha=0.5,
                    fontsize=fontsize, ha="left", va="center")


def _format_ra(ra_deg: float) -> str:
    """Format RA in hours:minutes."""
    ra_deg = ra_deg % 360
    hours = ra_deg / 15.0
    h = int(hours)
    m = int((hours - h) * 60)
    return f"{h:02d}h{m:02d}m"


def _format_dec(dec_deg: float) -> str:
    """Format Dec in degrees:arcminutes."""
    sign = "+" if dec_deg >= 0 else "-"
    dec_abs = abs(dec_deg)
    d = int(dec_abs)
    m = int((dec_abs - d) * 60)
    return f"{sign}{d:d}\u00b0{m:02d}'"


# ---------------------------------------------------------------------------
#  Constellation lines
# ---------------------------------------------------------------------------

def _draw_constellation_lines(ax, wcs: WCS, constellation_lines: list[dict],
                              width: int, height: int,
                              line_width: float, base_font: float) -> None:
    """Project constellation RA/Dec segments to pixel space and draw them.

    Segments are clipped to the image rectangle using Liang-Barsky so that
    lines extending off-frame are drawn up to the image edge.
    """
    for constellation in constellation_lines:
        seg_px = []
        for seg in constellation["segments"]:
            coord1 = SkyCoord(ra=seg["ra1"] * u.deg, dec=seg["dec1"] * u.deg)
            coord2 = SkyCoord(ra=seg["ra2"] * u.deg, dec=seg["dec2"] * u.deg)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                x1, y1 = wcs.world_to_pixel(coord1)
                x2, y2 = wcs.world_to_pixel(coord2)

            x1, y1, x2, y2 = float(x1), float(y1), float(x2), float(y2)
            if any(np.isnan(v) for v in [x1, y1, x2, y2]):
                continue

            clipped = _clip_segment(x1, y1, x2, y2,
                                    -0.5, -0.5, width - 0.5, height - 0.5)
            if clipped is None:
                continue

            cx1, cy1, cx2, cy2 = clipped
            ax.plot([cx1, cx2], [cy1, cy2],
                    color="#4488ff", alpha=0.4, linewidth=line_width)
            seg_px.append((cx1, cy1, cx2, cy2))

        # Label at centroid of drawn segments
        if seg_px:
            all_x = [s[0] for s in seg_px] + [s[2] for s in seg_px]
            all_y = [s[1] for s in seg_px] + [s[3] for s in seg_px]
            cx, cy = np.mean(all_x), np.mean(all_y)
            if 0 <= cx <= width and 0 <= cy <= height:
                ax.text(cx, cy, constellation["abbr"],
                        color="#4488ff", alpha=0.6,
                        fontsize=base_font * 0.9, fontweight="bold",
                        ha="center", va="center")


# ---------------------------------------------------------------------------
#  Catalog objects
# ---------------------------------------------------------------------------

def _prominence(mag: float | None) -> float:
    """Map magnitude to a prominence factor (0.3 to 1.0)."""
    if mag is None:
        return 0.5
    mag = max(0.0, min(mag, 15.0))
    return max(0.3, 1.0 - (mag - 1.0) * (0.7 / 12.0))


def _draw_objects(
    ax,
    objects: list[dict],
    color: str,
    base_marker: float,
    base_font: float,
    plate_scale_arcsec: float = 1.0,
    label_key: str = "name",
    common_name_key: str | None = None,
) -> None:
    """Draw catalog objects with circles sized to real angular extent."""
    for obj in objects:
        x, y = obj["x"], obj["y"]
        mag = obj.get("magnitude")
        size_arcmin = obj.get("major_axis_arcmin")
        p = _prominence(mag)

        if size_arcmin and size_arcmin > 0 and plate_scale_arcsec > 0:
            radius_px = (size_arcmin * 60.0) / plate_scale_arcsec / 2.0
            radius_px = max(radius_px, base_marker * 0.3)
        else:
            radius_px = base_marker * 0.5 * (0.5 + p)

        font_size = base_font * (0.6 + 0.8 * p)
        alpha = 0.4 + 0.55 * p
        lw = 0.6 + 1.8 * p

        circle = plt.Circle(
            (x, y), radius_px,
            fill=False, edgecolor=color, linewidth=lw, alpha=alpha,
        )
        ax.add_patch(circle)

        label = obj.get(label_key, "")
        if common_name_key and obj.get(common_name_key):
            label = f"{label}\n{obj[common_name_key]}"

        label_offset = radius_px + base_font * 0.3
        ax.text(
            x + label_offset, y - label_offset,
            label,
            color=color,
            fontsize=font_size,
            fontweight="bold" if p > 0.7 else "normal",
            alpha=alpha,
            ha="left",
            va="top",
        )
