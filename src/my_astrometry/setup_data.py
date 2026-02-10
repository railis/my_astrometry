"""One-time download of index files and catalog data."""

import pathlib
import urllib.request

DEFAULT_DATA_DIR = pathlib.Path.home() / ".my_astrometry"

OPENNGC_URL = (
    "https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files/NGC.csv"
)


def get_data_dir(data_dir: pathlib.Path | None = None) -> pathlib.Path:
    return data_dir or DEFAULT_DATA_DIR


def get_index_dir(data_dir: pathlib.Path | None = None) -> pathlib.Path:
    return get_data_dir(data_dir) / "index_files"


def get_openngc_path(data_dir: pathlib.Path | None = None) -> pathlib.Path:
    return get_data_dir(data_dir) / "openngc" / "NGC.csv"


def check_setup(data_dir: pathlib.Path | None = None) -> bool:
    """Return True if data has been downloaded."""
    index_dir = get_index_dir(data_dir)
    openngc_path = get_openngc_path(data_dir)
    return index_dir.exists() and any(index_dir.iterdir()) and openngc_path.exists()


DEFAULT_5200_SCALES = {2, 3, 4, 5, 6}  # [4-22] arcmin, ~9.3 GB


def download_index_files(data_dir: pathlib.Path | None = None) -> list[pathlib.Path]:
    """Download astrometry.net 4100 series index files."""
    import astrometry

    index_dir = get_index_dir(data_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading 4100 series index files to {index_dir}...")
    print("This is ~355 MB and may take a few minutes.")
    index_files = astrometry.series_4100.index_files(
        cache_directory=index_dir,
    )
    print(f"Downloaded {len(index_files)} index files.")
    return index_files


def download_5200_index_files(
    data_dir: pathlib.Path | None = None,
    scales: set[int] | None = None,
) -> list[pathlib.Path]:
    """Download astrometry.net 5200 series (LIGHT) index files.

    The 5200 series is built from Tycho-2 + Gaia-DR2, healpix-partitioned
    (48 sky tiles per scale), covering fine scales [2-22] arcmin.
    """
    import astrometry

    scales = scales if scales is not None else DEFAULT_5200_SCALES
    index_dir = get_index_dir(data_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    total_size = astrometry.series_5200.size(scales=scales)
    size_str = f"{total_size / 1e9:.1f} GB"
    print(f"Downloading 5200 series (scales {sorted(scales)}) to {index_dir}...")
    print(f"This is ~{size_str} and may take a while.")
    index_files = astrometry.series_5200.index_files(
        cache_directory=index_dir,
        scales=scales,
    )
    print(f"Downloaded {len(index_files)} index files.")
    return index_files


def download_openngc(data_dir: pathlib.Path | None = None) -> pathlib.Path:
    """Download OpenNGC catalog CSV."""
    openngc_path = get_openngc_path(data_dir)
    openngc_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading OpenNGC catalog to {openngc_path}...")
    urllib.request.urlretrieve(OPENNGC_URL, openngc_path)
    size_mb = openngc_path.stat().st_size / (1024 * 1024)
    print(f"Downloaded OpenNGC catalog ({size_mb:.1f} MB).")
    return openngc_path


def run_setup(
    data_dir: pathlib.Path | None = None,
    include_5200: bool = False,
    scales_5200: set[int] | None = None,
) -> None:
    """Run full setup: download index files and catalogs."""
    data_dir_resolved = get_data_dir(data_dir)
    print(f"Setting up data in {data_dir_resolved}")
    print()
    download_index_files(data_dir)
    if include_5200:
        print()
        download_5200_index_files(data_dir, scales=scales_5200)
    print()
    download_openngc(data_dir)
    print()
    print("Setup complete! You can now use 'my-astrometry solve'.")
