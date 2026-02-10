# my-astrometry

A local, offline astrometry plate solver and image annotator for astrophotography.

**my-astrometry** takes your astrophotography images, identifies the star field by matching detected stars against index catalogs (plate solving), and produces annotated images with labeled deep-sky objects (Messier, NGC, IC) and constellation stick figures overlaid on the original photo.

Everything runs locally on your machine — no internet connection needed after the initial setup.

## Features

- **Offline plate solving** powered by [astrometry.net](https://astrometry.net/) index files via the [`astrometry`](https://github.com/neuromorphicsystems/astrometry) Python wrapper
- **Automatic star detection** using [SEP](https://github.com/pjw/sep) (Source Extractor Python)
- **WCS persistence** — solve once, annotate many times. The WCS solution is saved directly into the image's EXIF metadata
- **Deep-sky catalog overlays** — Messier (all 110 objects), NGC, and IC catalogs with circles sized by real angular extent
- **Constellation lines** — stick-figure constellation overlays projected onto your image
- **RA/Dec coordinate grid** — with labeled axes along the image edges
- **Smart hints** — provide an object name (e.g., `m42`), coordinates, or plate scale to speed up solving
- **EXIF-aware** — automatically estimates plate scale from your camera's EXIF focal length data

## Installation

Requires **Python 3.10** or later.

```bash
# Clone the repository
git clone https://github.com/yourusername/my-astrometry.git
cd my-astrometry

# Install in editable mode (using uv or pip)
uv pip install -e .
# or
pip install -e .
```

### One-time data setup

Before solving images, you need to download the astrometry.net index files and the OpenNGC catalog. This downloads approximately **360 MB** of data to `~/.my_astrometry/`.

```bash
my-astrometry setup
```

## Usage

The workflow is split into two steps: **solve** (identify the sky coordinates) and **annotate** (render the overlay).

### Quick start

```bash
# 1. Plate-solve your image (saves WCS solution into image EXIF)
my-astrometry solve photo.jpg

# 2. Annotate with deep-sky objects and constellations
my-astrometry annotate photo.jpg
# -> Produces photo_annotated.png
```

---

### `setup` — Download index files and catalogs

Downloads the astrometry.net 4100 series index files and the OpenNGC catalog. Data is stored in `~/.my_astrometry/` by default.

```bash
my-astrometry setup
```

| Option | Description |
|--------|-------------|
| `--data-dir PATH` | Store data in a custom directory instead of `~/.my_astrometry/` |
| `--with-5200` | Also download the 5200 series index files (Gaia-DR2 based, healpix-partitioned). Adds ~9 GB for default scales |
| `--scales 2,3,4` | Comma-separated list of 5200 scale numbers to download (default: 2,3,4,5,6). Requires `--with-5200` |
| `--all-scales` | Download all 5200 scales 0 through 6 (~35 GB). Requires `--with-5200` |

**When do you need the 5200 series?** The 4100 series covers most use cases. The 5200 series (based on Gaia DR2) provides denser star coverage and can help with very narrow fields of view or crowded star fields.

---

### `solve` — Plate-solve an image

Detects stars in your image, matches the star pattern against index files, and determines the precise sky coordinates (WCS — World Coordinate System). The solution is saved into the image's EXIF metadata so that annotation can be done later without re-solving.

```bash
my-astrometry solve photo.jpg
```

**Output:**
```
Solving: photo.jpg
  Loading image (photo.jpg)...
  Image size: 6000 x 4000
  Detecting stars...
  Found 200 stars (0.3s)
  EXIF auto-detected scale: 1.05-1.95 arcsec/px
  Plate solving (42 index files)...
  Solved in 2.1s (total 2.8s)

Solution:
  Center RA:   05h 35m 17.2s (83.8217 deg)
  Center Dec:  -05° 23' 28.0" (-5.3911 deg)
  FOV:         1.75° x 1.17° (105.0' x 70.0')
  Plate scale: 1.05 arcsec/pixel
  Log-odds:    78.3
  Stars used:  200

WCS solution saved to EXIF in photo.jpg
```

#### Solving hints

Providing hints dramatically speeds up solving. Without hints, the solver must search the entire sky at all scales — with good hints, solving typically completes in 1-3 seconds.

| Option | Description |
|--------|-------------|
| `--object NAME` | Search near a named object. Accepts Messier names (`m42`, `M31`), NGC/IC numbers (`NGC1976`, `IC434`), or common names (`"orion nebula"`, `"andromeda galaxy"`). Automatically sets RA/Dec/radius hints |
| `--ra DEGREES` | Right ascension hint in degrees (0-360). Use together with `--dec` |
| `--dec DEGREES` | Declination hint in degrees (-90 to +90). Use together with `--ra` |
| `--radius DEGREES` | Search radius around the RA/Dec hint (default: 10 degrees when `--object` is used) |
| `--scale-low ARCSEC` | Lower bound of expected plate scale in arcsec/pixel |
| `--scale-high ARCSEC` | Upper bound of expected plate scale in arcsec/pixel |

**Examples:**

```bash
# Search near the Orion Nebula
my-astrometry solve photo.jpg --object m42

# Search near specific coordinates with a tight radius
my-astrometry solve photo.jpg --ra 83.8 --dec -5.4 --radius 5

# Constrain the plate scale (useful for known optical setups)
my-astrometry solve photo.jpg --scale-low 1.0 --scale-high 2.0

# Combine hints for fastest solving
my-astrometry solve photo.jpg --object m42 --scale-low 1.0 --scale-high 2.0
```

#### Advanced options

| Option | Description |
|--------|-------------|
| `--max-stars N` | Maximum number of stars to detect (default: 200). Increase for dense fields, decrease for speed |
| `--parity {both,normal,flip}` | Image parity. Set to `normal` for a ~2x speedup if you know your image is not mirrored. Default: `both` |
| `--no-exif-scale` | Disable automatic plate scale estimation from camera EXIF data (focal length, sensor size) |
| `--no-save-exif` | Solve the image but don't write the WCS solution to its EXIF. Useful for read-only files or when you just want the coordinates printed |
| `--data-dir PATH` | Use index files from a custom directory |
| `-v, --verbose` | Print detailed progress and debug information |

---

### `annotate` — Annotate a solved image

Reads the WCS solution from a previously solved image's EXIF metadata and renders an annotated overlay with deep-sky objects and constellation lines. The annotated image is saved as a separate file (the original is not modified).

```bash
my-astrometry annotate photo.jpg
# -> Saves photo_annotated.png
```

**The image must be solved first.** If no WCS solution is found in the EXIF, you'll see:
```
Error: No WCS solution found in photo.jpg.
Run 'my-astrometry solve photo.jpg' first.
```

#### Annotation options

| Option | Description |
|--------|-------------|
| `-o, --output PATH` | Custom output path (default: `<input>_annotated.png`) |
| `--no-ngc` | Don't label NGC objects. Useful to reduce clutter on wide-field images |
| `--mag-limit MAG` | Only show IC/NGC objects brighter than this magnitude (e.g., `--mag-limit 10`). Lower values = fewer, brighter objects |
| `--data-dir PATH` | Use catalog data from a custom directory |
| `-v, --verbose` | Print detailed progress information |

**Examples:**

```bash
# Annotate with default settings
my-astrometry annotate photo.jpg

# Save to a specific output file
my-astrometry annotate photo.jpg -o ~/output/annotated.png

# Clean annotation: only Messier and IC objects, nothing fainter than magnitude 10
my-astrometry annotate photo.jpg --no-ngc --mag-limit 10
```

#### What gets annotated

- **Messier objects** (yellow) — all 110 Messier catalog objects with common names
- **IC objects** (cyan) — Index Catalogue objects within the field of view
- **NGC objects** (green) — New General Catalogue objects (disable with `--no-ngc`)
- **Constellation lines** (blue) — stick-figure constellation outlines
- **RA/Dec grid** (white, dashed) — coordinate grid with labeled axes

Object circles are sized proportionally to their real angular extent using the image's plate scale. Brighter objects appear with stronger labels and thicker outlines.

## How it works

1. **Star Detection** — SEP (Source Extractor Python) subtracts the sky background and identifies point sources, returning pixel coordinates sorted by brightness
2. **Plate Solving** — the detected star pattern is matched against pre-computed quad index files using the astrometry.net solver, yielding a WCS (World Coordinate System) transformation
3. **WCS Persistence** — the WCS solution is serialized as JSON and stored in the image's EXIF `UserComment` tag, enabling annotation without re-solving
4. **Catalog Matching** — Messier, NGC, and IC catalog positions are projected through the WCS to determine which objects fall within the image's field of view
5. **Annotation Rendering** — matplotlib renders constellation lines, object labels, and a coordinate grid as an overlay on the original image

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| [astrometry](https://github.com/neuromorphicsystems/astrometry) | >= 4.0 | Python wrapper for the astrometry.net plate solver |
| [sep-pjw](https://github.com/pjw/sep) | >= 1.2 | Source Extractor for star detection |
| [astropy](https://www.astropy.org/) | >= 5.0 | WCS handling, FITS headers, sky coordinates |
| [numpy](https://numpy.org/) | >= 1.24 | Array operations |
| [Pillow](https://python-pillow.org/) | >= 9.0 | Image I/O and EXIF handling |
| [matplotlib](https://matplotlib.org/) | >= 3.6 | Annotation rendering |

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
