"""CLI entry point for my-astrometry."""

import argparse
import pathlib
import sys


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="my-astrometry",
        description="Local offline astrometry plate solver with annotation.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- setup command ---
    setup_parser = subparsers.add_parser(
        "setup",
        help="Download index files and catalog data (~360 MB, one-time).",
    )
    setup_parser.add_argument(
        "--data-dir",
        type=pathlib.Path,
        default=None,
        help="Directory to store downloaded data (default: ~/.my_astrometry).",
    )
    setup_parser.add_argument(
        "--with-5200",
        action="store_true",
        help="Also download 5200 series (Gaia-DR2, healpix-partitioned, ~9 GB for default scales).",
    )
    setup_parser.add_argument(
        "--scales",
        type=str,
        default=None,
        help="Comma-separated 5200 scale numbers to download (default: 2,3,4,5,6). Requires --with-5200.",
    )
    setup_parser.add_argument(
        "--all-scales",
        action="store_true",
        help="Download all 5200 scales 0-6 (~35 GB). Requires --with-5200.",
    )

    # --- solve command ---
    solve_parser = subparsers.add_parser(
        "solve",
        help="Plate-solve an image and save WCS solution to EXIF.",
    )
    solve_parser.add_argument(
        "image",
        type=pathlib.Path,
        help="Path to the input image (JPG/PNG/TIFF).",
    )
    solve_parser.add_argument(
        "--data-dir",
        type=pathlib.Path,
        default=None,
        help="Directory containing downloaded data (default: ~/.my_astrometry).",
    )
    solve_parser.add_argument(
        "--object",
        type=str,
        default=None,
        help="Object name to center search on (e.g., m42, NGC1976, 'orion nebula').",
    )
    solve_parser.add_argument(
        "--scale-low",
        type=float,
        default=None,
        help="Lower bound of plate scale in arcsec/pixel.",
    )
    solve_parser.add_argument(
        "--scale-high",
        type=float,
        default=None,
        help="Upper bound of plate scale in arcsec/pixel.",
    )
    solve_parser.add_argument(
        "--ra",
        type=float,
        default=None,
        help="RA hint in degrees (0-360).",
    )
    solve_parser.add_argument(
        "--dec",
        type=float,
        default=None,
        help="Dec hint in degrees (-90 to 90).",
    )
    solve_parser.add_argument(
        "--radius",
        type=float,
        default=None,
        help="Search radius in degrees (requires --ra and --dec).",
    )
    solve_parser.add_argument(
        "--max-stars",
        type=int,
        default=200,
        help="Maximum number of stars to detect (default: 200).",
    )
    solve_parser.add_argument(
        "--parity",
        choices=["both", "normal", "flip"],
        default="both",
        help="Image parity: 'normal' skips flip check (~2x speedup), 'both' tries both (default).",
    )
    solve_parser.add_argument(
        "--no-exif-scale",
        action="store_true",
        help="Disable automatic plate scale estimation from EXIF data.",
    )
    solve_parser.add_argument(
        "--no-save-exif",
        action="store_true",
        help="Don't save WCS solution to image EXIF.",
    )
    solve_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress information.",
    )

    # --- annotate command ---
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate a solved image with catalog objects and constellations.",
    )
    annotate_parser.add_argument(
        "image",
        type=pathlib.Path,
        help="Path to a plate-solved image (with WCS in EXIF).",
    )
    annotate_parser.add_argument(
        "-o", "--output",
        type=pathlib.Path,
        default=None,
        help="Output path for annotated image (default: <input>_annotated.png).",
    )
    annotate_parser.add_argument(
        "--data-dir",
        type=pathlib.Path,
        default=None,
        help="Directory containing downloaded data (default: ~/.my_astrometry).",
    )
    annotate_parser.add_argument(
        "--no-ngc",
        action="store_true",
        help="Don't label NGC objects (reduces clutter).",
    )
    annotate_parser.add_argument(
        "--mag-limit",
        type=float,
        default=None,
        help="Magnitude limit for IC/NGC objects (e.g., 12).",
    )
    annotate_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress information.",
    )

    args = parser.parse_args()

    if args.command == "setup":
        _cmd_setup(args)
    elif args.command == "solve":
        _cmd_solve(args)
    elif args.command == "annotate":
        _cmd_annotate(args)


def _cmd_setup(args: argparse.Namespace) -> None:
    from my_astrometry.setup_data import run_setup

    scales_5200 = None
    if args.with_5200:
        if args.all_scales:
            scales_5200 = {0, 1, 2, 3, 4, 5, 6}
        elif args.scales:
            scales_5200 = {int(s.strip()) for s in args.scales.split(",")}
        # else: None -> uses DEFAULT_5200_SCALES in setup_data.py

    run_setup(
        data_dir=args.data_dir,
        include_5200=args.with_5200,
        scales_5200=scales_5200,
    )


def _cmd_solve(args: argparse.Namespace) -> None:
    from my_astrometry.setup_data import check_setup, get_index_dir, get_openngc_path
    from my_astrometry.solver import solve
    from my_astrometry.coordinates import get_field_info, print_field_info

    # Check that setup has been run
    if not check_setup(args.data_dir):
        print("Error: Data not found. Run 'my-astrometry setup' first.", file=sys.stderr)
        sys.exit(1)

    # Validate image exists
    if not args.image.exists():
        print(f"Error: Image not found: {args.image}", file=sys.stderr)
        sys.exit(1)

    # Resolve --object to RA/Dec hints
    ra_hint = args.ra
    dec_hint = args.dec
    radius_hint = args.radius

    if args.object:
        from my_astrometry.catalogs import lookup_object

        openngc_path = get_openngc_path(args.data_dir)
        result = lookup_object(args.object, openngc_path=openngc_path)
        if result is None:
            print(f"Error: Unknown object '{args.object}'.", file=sys.stderr)
            sys.exit(1)
        obj_ra, obj_dec, obj_display = result
        print(f"Object: {obj_display} (RA={obj_ra:.4f}, Dec={obj_dec:.4f})")
        if ra_hint is None:
            ra_hint = obj_ra
        if dec_hint is None:
            dec_hint = obj_dec
        if radius_hint is None:
            radius_hint = 10.0

    # Solve
    print(f"Solving: {args.image}")
    try:
        wcs, metadata = solve(
            image_path=args.image,
            index_dir=get_index_dir(args.data_dir),
            scale_low=args.scale_low,
            scale_high=args.scale_high,
            ra_hint=ra_hint,
            dec_hint=dec_hint,
            radius_hint=radius_hint,
            max_stars=args.max_stars,
            verbose=args.verbose,
            parity=args.parity,
            use_exif=not args.no_exif_scale,
        )
    except RuntimeError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Print results
    width = metadata["image_width"]
    height = metadata["image_height"]
    field_info = get_field_info(wcs, width, height)

    print()
    print("Solution:")
    print_field_info(field_info)
    print(f"  Log-odds:    {metadata['logodds']:.1f}")
    print(f"  Stars used:  {metadata['num_stars_detected']}")

    # Save WCS to EXIF
    if not args.no_save_exif:
        from my_astrometry.wcs_exif import save_solution

        save_solution(args.image, wcs, metadata)
        print(f"\nWCS solution saved to EXIF in {args.image}")
    else:
        print("\nWCS solution NOT saved (--no-save-exif).")


def _cmd_annotate(args: argparse.Namespace) -> None:
    from my_astrometry.setup_data import check_setup, get_openngc_path
    from my_astrometry.wcs_exif import load_solution
    from my_astrometry.catalogs import get_catalog_objects, get_constellation_lines_in_fov
    from my_astrometry.annotator import create_annotated_image

    # Check that setup has been run
    if not check_setup(args.data_dir):
        print("Error: Data not found. Run 'my-astrometry setup' first.", file=sys.stderr)
        sys.exit(1)

    # Validate image exists
    if not args.image.exists():
        print(f"Error: Image not found: {args.image}", file=sys.stderr)
        sys.exit(1)

    # Load WCS from EXIF
    result = load_solution(args.image)
    if result is None:
        print(
            f"Error: No WCS solution found in {args.image}.\n"
            f"Run 'my-astrometry solve {args.image}' first.",
            file=sys.stderr,
        )
        sys.exit(1)

    wcs, metadata = result
    width = metadata["image_width"]
    height = metadata["image_height"]

    print(f"Annotating: {args.image}")
    print(f"  WCS center: RA={metadata['center_ra_deg']:.4f}, "
          f"Dec={metadata['center_dec_deg']:.4f}")

    openngc_path = get_openngc_path(args.data_dir)

    print("  Loading catalogs...")
    catalog_objects = get_catalog_objects(
        wcs, width, height,
        openngc_path=openngc_path,
        mag_limit=args.mag_limit,
        include_ngc=not args.no_ngc,
    )

    constellation_lines = get_constellation_lines_in_fov(wcs, width, height)

    # Count objects
    n_messier = len(catalog_objects["messier"])
    n_ic = len(catalog_objects["ic"])
    n_ngc = len(catalog_objects["ngc"])
    n_const = len(constellation_lines)

    print(f"  Found {n_messier} Messier, {n_ic} IC, {n_ngc} NGC objects in FOV")
    print(f"  Found {n_const} constellations in FOV")

    # List Messier objects by name
    if n_messier > 0:
        labels = []
        for obj in catalog_objects["messier"]:
            label = obj["name"]
            if obj.get("common_name"):
                label += f" ({obj['common_name']})"
            labels.append(label)
        print(f"  Messier: {', '.join(labels)}")

    # Determine output path
    output_path = args.output
    if output_path is None:
        output_path = args.image.with_name(args.image.stem + "_annotated.png")

    print(f"  Rendering overlay...")

    create_annotated_image(
        image_path=args.image,
        wcs=wcs,
        catalog_objects=catalog_objects,
        constellation_lines=constellation_lines,
        output_path=output_path,
        plate_scale_arcsec=metadata["scale_arcsec_per_pixel"],
        show_ngc=not args.no_ngc,
        verbose=args.verbose,
    )

    print(f"  Saved: {output_path}")


if __name__ == "__main__":
    main()
