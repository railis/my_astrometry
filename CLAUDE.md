# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Environment

- Python 3.12 on WSL2 Ubuntu
- Use `uv` (at `~/.local/bin/uv`) instead of pip — no system pip available
- Virtual env: `.venv/` managed via `uv pip install -e .`
- Install/reinstall after code changes: `uv pip install -e .`

## Commands

```bash
# Install dependencies and package in editable mode
uv pip install -e .

# Run the CLI
my-astrometry setup              # One-time: downloads ~360MB of index files + OpenNGC catalog to ~/.my_astrometry
my-astrometry solve image.jpg    # Plate-solve and annotate an image
my-astrometry solve image.jpg --no-annotate  # Solve only, print coordinates
```

No test suite, linter, or formatter is configured.

## Architecture

Pipeline: **Image → Star Detection → Plate Solving → Coordinate Extraction → Catalog Matching → Annotation**

### Data flow through `cli.py _cmd_solve()`

1. `solver.solve()` — loads image via Pillow, detects stars with `sep` (Source Extractor), runs `astrometry.Solver` to get WCS solution
2. `coordinates.get_field_info()` — extracts center RA/Dec, FOV, plate scale from WCS
3. `catalogs.get_catalog_objects()` — loads Messier (bundled CSV) + OpenNGC (downloaded CSV), filters to objects in FOV using `wcs.world_to_pixel()`, returns pixel coords
4. `catalogs.get_constellation_lines_in_fov()` — loads bundled constellation JSON, projects line segments to pixel coords
5. `annotator.create_annotated_image()` — renders overlay on original image using matplotlib (plain axes, not WCSAxes)

### Key design decisions

- All catalog/constellation coordinates are converted to **pixel space** in `catalogs.py` before reaching the annotator. The annotator works entirely in pixel coordinates.
- Messier catalog is bundled (`data/messier.csv`) for reliability; NGC/IC come from downloaded OpenNGC CSV (semicolon-delimited, parsed manually).
- Constellation data: `data/constellations.json` (line segments as RA/Dec pairs) + `data/constellation_names.json` (IAU abbreviation → full name).
- Index file pre-filtering in `solver.py:_filter_index_files()` reduces solving time by matching quad angular scales to the expected plate scale range.
- The `astrometry` package (neuromorphicsystems/astrometry) wraps astrometry.net's C solver. Use `sep-pjw` (NOT `sep`) for star extraction.

### Module responsibilities

- **`setup_data.py`**: Downloads 4100 series index files via `astrometry.series_4100.index_files()` and OpenNGC CSV. Data stored in `~/.my_astrometry/`.
- **`solver.py`**: Star detection (`sep.Background` → `sep.extract`), plate solving (`astrometry.Solver`). Returns `(WCS, metadata_dict)`.
- **`coordinates.py`**: Pure WCS → human-readable coordinate extraction (center, corners, FOV, plate scale).
- **`catalogs.py`**: Loads/filters Messier, OpenNGC, constellation data. `lookup_object()` resolves names (M42, "orion nebula", NGC1976) to RA/Dec. `filter_objects_in_fov()` does bulk WCS projection.
- **`annotator.py`**: Matplotlib rendering. Object circles sized by real angular extent via plate scale. Prominence (alpha, stroke, font) scaled by magnitude.
