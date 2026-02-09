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


def download_openngc(data_dir: pathlib.Path | None = None) -> pathlib.Path:
    """Download OpenNGC catalog CSV."""
    openngc_path = get_openngc_path(data_dir)
    openngc_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading OpenNGC catalog to {openngc_path}...")
    urllib.request.urlretrieve(OPENNGC_URL, openngc_path)
    size_mb = openngc_path.stat().st_size / (1024 * 1024)
    print(f"Downloaded OpenNGC catalog ({size_mb:.1f} MB).")
    return openngc_path


def run_setup(data_dir: pathlib.Path | None = None) -> None:
    """Run full setup: download index files and catalogs."""
    data_dir_resolved = get_data_dir(data_dir)
    print(f"Setting up data in {data_dir_resolved}")
    print()
    download_index_files(data_dir)
    print()
    download_openngc(data_dir)
    print()
    print("Setup complete! You can now use 'my-astrometry solve'.")
