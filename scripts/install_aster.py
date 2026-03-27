"""Install official ASTER branch binaries into work/tools/aster/current."""

from __future__ import annotations

import argparse
import platform
import subprocess
import zipfile
from pathlib import Path


PLATFORM_BRANCHES = {
    ("Linux", "x86_64"): ("Linux", "ASTER-Linux"),
    ("Linux", "amd64"): ("Linux", "ASTER-Linux"),
    ("Darwin", "arm64"): ("MacOS", "ASTER-MacOS"),
    ("Darwin", "aarch64"): ("MacOS", "ASTER-MacOS"),
    ("Darwin", "x86_64"): ("MacOSx86", "ASTER-MacOSx86"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--install-root",
        default="work/tools/aster",
        help="Root directory for downloaded ASTER archives and the current symlink.",
    )
    return parser.parse_args()


def select_branch() -> tuple[str, str]:
    key = (platform.system(), platform.machine())
    if key not in PLATFORM_BRANCHES:
        raise RuntimeError(f"Unsupported platform for automatic ASTER install: {key}")
    return PLATFORM_BRANCHES[key]


def download_archive(branch: str, archive_path: Path) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    if archive_path.is_file():
        return

    url = f"https://github.com/chaoszhang/ASTER/archive/refs/heads/{branch}.zip"
    subprocess.run(["curl", "-fsSL", "-o", archive_path.as_posix(), url], check=True)


def extract_archive(archive_path: Path, extract_root: Path, extracted_name: str) -> Path:
    extracted_dir = extract_root / extracted_name
    if extracted_dir.is_dir():
        return extracted_dir

    extract_root.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as handle:
        handle.extractall(extract_root)

    if not extracted_dir.is_dir():
        raise RuntimeError(f"Expected extracted ASTER directory at {extracted_dir}")
    return extracted_dir


def build_required_binaries(extracted_dir: Path) -> None:
    subprocess.run(["make", "wastral", "astral"], cwd=extracted_dir, check=True)


def ensure_current_symlink(install_root: Path, extracted_dir: Path) -> tuple[Path, Path]:
    current_link = install_root / "current"
    if current_link.exists() or current_link.is_symlink():
        current_link.unlink()
    current_link.symlink_to(extracted_dir.resolve(), target_is_directory=True)

    wastral = current_link / "bin" / "wastral"
    astral4 = current_link / "bin" / "astral4"
    for executable in (wastral, astral4):
        if not executable.is_file():
            raise RuntimeError(f"Installed ASTER executable not found at {executable}")
        executable.chmod(0o755)
    return wastral, astral4


def main() -> int:
    args = parse_args()
    install_root = Path(args.install_root)
    branch, extracted_name = select_branch()
    archive_path = install_root / "downloads" / f"{extracted_name}.zip"

    download_archive(branch, archive_path)
    extracted_dir = extract_archive(archive_path, install_root / "releases", extracted_name)
    build_required_binaries(extracted_dir)
    wastral, _ = ensure_current_symlink(install_root, extracted_dir)
    print(wastral)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
