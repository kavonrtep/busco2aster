"""Install an official IQ-TREE 3 binary into work/tools/iqtree3/current."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import tarfile
import zipfile
from pathlib import Path


LATEST_RELEASE_API = "https://api.github.com/repos/iqtree/iqtree3/releases/latest"
PLATFORM_SUFFIXES = {
    ("Linux", "x86_64"): ("Linux-intel", ".tar.gz"),
    ("Linux", "amd64"): ("Linux-intel", ".tar.gz"),
    ("Linux", "aarch64"): ("Linux-arm", ".tar.gz"),
    ("Linux", "arm64"): ("Linux-arm", ".tar.gz"),
    ("Darwin", "x86_64"): ("macOS-intel", ".zip"),
    ("Darwin", "arm64"): ("macOS-arm", ".zip"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--install-root",
        default="work/tools/iqtree3",
        help="Root directory for downloaded releases and the current symlink.",
    )
    return parser.parse_args()


def load_latest_release() -> dict:
    result = subprocess.run(
        ["curl", "-fsSL", LATEST_RELEASE_API],
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout)


def select_asset(release: dict) -> dict:
    key = (platform.system(), platform.machine())
    if key not in PLATFORM_SUFFIXES:
        raise RuntimeError(f"Unsupported platform for automatic IQ-TREE 3 install: {key}")
    suffix, extension = PLATFORM_SUFFIXES[key]
    expected_name = f"iqtree-{release['tag_name'].removeprefix('v')}-{suffix}{extension}"
    for asset in release.get("assets", []):
        if asset.get("name") == expected_name:
            return asset
    raise RuntimeError(f"Could not find matching IQ-TREE asset {expected_name!r} in release metadata.")


def download_asset(asset: dict, archive_path: Path) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    if archive_path.is_file():
        return
    subprocess.run(
        ["curl", "-fsSL", "-o", archive_path.as_posix(), asset["browser_download_url"]],
        check=True,
    )


def verify_asset_digest(asset: dict, archive_path: Path) -> None:
    digest = asset.get("digest", "")
    if not digest.startswith("sha256:"):
        return
    expected = digest.split(":", 1)[1]
    actual = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError(
            f"Downloaded IQ-TREE archive checksum mismatch for {archive_path}: "
            f"expected {expected}, observed {actual}"
        )


def extract_archive(archive_path: Path, extract_root: Path, extracted_name: str) -> Path:
    extract_root.mkdir(parents=True, exist_ok=True)
    if archive_path.suffix == ".zip":
        with zipfile.ZipFile(archive_path) as handle:
            handle.extractall(extract_root)
    else:
        with tarfile.open(archive_path, "r:gz") as handle:
            handle.extractall(extract_root)
    extracted_dir = extract_root / extracted_name
    if not extracted_dir.is_dir():
        raise RuntimeError(f"Expected extracted IQ-TREE directory at {extracted_dir}")
    return extracted_dir


def ensure_current_symlink(install_root: Path, extracted_dir: Path) -> Path:
    current_link = install_root / "current"
    if current_link.exists() or current_link.is_symlink():
        current_link.unlink()
    current_link.symlink_to(extracted_dir.resolve(), target_is_directory=True)
    executable = current_link / "bin" / "iqtree3"
    if not executable.is_file():
        raise RuntimeError(f"Installed IQ-TREE executable not found at {executable}")
    return executable


def main() -> int:
    args = parse_args()
    install_root = Path(args.install_root)
    release = load_latest_release()
    asset = select_asset(release)

    archive_path = install_root / "releases" / asset["name"]
    download_asset(asset, archive_path)
    verify_asset_digest(asset, archive_path)

    extracted_name = asset["name"].removesuffix(".tar.gz").removesuffix(".zip")
    extracted_dir = install_root / "releases" / extracted_name
    if not extracted_dir.is_dir():
        extracted_dir = extract_archive(archive_path, install_root / "releases", extracted_name)

    executable = ensure_current_symlink(install_root, extracted_dir)
    print(executable)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
