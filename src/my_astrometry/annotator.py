"""Overlay rendering: constellation lines, Messier/IC/NGC labels on solved image."""

import pathlib

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
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

    Args:
        image_path: Path to the original image.
        wcs: WCS solution.
        catalog_objects: Dict with 'messier', 'ic', 'ngc' lists from catalogs.get_catalog_objects().
        constellation_lines: List from catalogs.get_constellation_lines_in_fov().
        output_path: Where to save the annotated PNG.
        plate_scale_arcsec: Plate scale in arcsec/pixel (for sizing circles to real angular size).
        show_ngc: Whether to show NGC labels.
        verbose: Print progress.

    Returns:
        Path to the saved annotated image.
    """
    output_path = pathlib.Path(output_path)

    # Load original image
    img = Image.open(image_path)
    img_array = np.array(img)
    height, width = img_array.shape[:2]

    # Create figure at original resolution using plain axes (not WCSAxes)
    # All coordinates are already in pixel space from catalogs.py
    dpi = 100
    fig, ax = plt.subplots(
        figsize=(width / dpi, height / dpi),
        dpi=dpi,
    )
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Display the original image (origin="upper" so y=0 is top, matching pixel coords)
    ax.imshow(img_array, origin="upper", aspect="equal")
    ax.set_xlim(-0.5, width - 0.5)
    ax.set_ylim(height - 0.5, -0.5)
    ax.axis("off")

    # Base scaling factor for annotations based on image size
    img_scale = min(width, height) / 1000.0
    base_font = max(6, min(14, 9 * img_scale))
    base_marker = max(8, min(30, 15 * img_scale))
    line_width = max(0.5, min(2.0, 1.0 * img_scale))

    # Draw constellation lines
    for constellation in constellation_lines:
        for seg in constellation["segments"]:
            ax.plot(
                [seg["x1"], seg["x2"]],
                [seg["y1"], seg["y2"]],
                color="#4488ff",
                alpha=0.4,
                linewidth=line_width,
            )

    # Label constellation names (at centroid of their visible segments)
    for constellation in constellation_lines:
        segs = constellation["segments"]
        all_x = [s["x1"] for s in segs] + [s["x2"] for s in segs]
        all_y = [s["y1"] for s in segs] + [s["y2"] for s in segs]
        cx, cy = np.mean(all_x), np.mean(all_y)
        if 0 <= cx <= width and 0 <= cy <= height:
            ax.text(
                cx, cy,
                constellation["abbr"],
                color="#4488ff",
                alpha=0.6,
                fontsize=base_font * 0.9,
                fontweight="bold",
                ha="center",
                va="center",
            )

    # Draw catalog objects scaled by magnitude and angular size
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

    # Save
    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0,
        facecolor="black",
    )
    plt.close(fig)

    return output_path


def _prominence(mag: float | None) -> float:
    """Map magnitude to a prominence factor (0.3 to 1.0).

    Bright objects (mag ~1-4) get ~1.0, faint objects (mag ~12+) get ~0.3.
    Objects with no magnitude data get 0.5.
    """
    if mag is None:
        return 0.5
    mag = max(0.0, min(mag, 15.0))
    # Linear mapping: mag 1 -> 1.0, mag 13 -> 0.3
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

        # Circle radius in pixels: use real angular size if available,
        # otherwise fall back to a small marker scaled by prominence
        if size_arcmin and size_arcmin > 0 and plate_scale_arcsec > 0:
            radius_px = (size_arcmin * 60.0) / plate_scale_arcsec / 2.0
            # Enforce a minimum so tiny objects are still visible
            radius_px = max(radius_px, base_marker * 0.3)
        else:
            radius_px = base_marker * 0.5 * (0.5 + p)

        # Scale font, stroke, alpha by prominence
        font_size = base_font * (0.6 + 0.8 * p)
        alpha = 0.4 + 0.55 * p     # 0.4 for faintest, 0.95 for brightest
        lw = 0.6 + 1.8 * p         # 0.6 to 2.4

        # Draw circle marker
        circle = plt.Circle(
            (x, y), radius_px,
            fill=False, edgecolor=color, linewidth=lw, alpha=alpha,
        )
        ax.add_patch(circle)

        # Build label
        label = obj.get(label_key, "")
        if common_name_key and obj.get(common_name_key):
            label = f"{label}\n{obj[common_name_key]}"

        # Draw label offset outside the circle
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
