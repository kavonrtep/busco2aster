#!/usr/bin/env python3
"""Run the busco2aster workflow locally or inside an Apptainer/Singularity image."""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml


PIPELINE_DIR = Path("/opt/pipeline")
SNAKEFILE = PIPELINE_DIR / "Snakefile"
CONDA_ENVS_PATH = Path("/opt/conda/envs")
IQTREE_EXECUTABLE = Path("/opt/tools/iqtree3/current/bin/iqtree3")
WASTRAL_EXECUTABLE = Path("/opt/tools/aster/current/bin/wastral")
ASTRAL4_EXECUTABLE = Path("/opt/tools/aster/current/bin/astral4")
REQUIRED_SAMPLE_COLUMNS = ("sample_id", "taxon_id", "assembly_fasta")


@dataclass(frozen=True)
class InputPathSpec:
    path: Path
    raw_text: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-c",
        "--config",
        default="config/config.yaml",
        help="Workflow config path. Relative paths are resolved against --repo-root.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Maximum Snakemake core count.",
    )
    parser.add_argument(
        "--target",
        help="Optional Snakemake target file or rule. Defaults to the workflow 'all' target.",
    )
    parser.add_argument(
        "--directory",
        help=(
            "Writable Snakemake working directory. Defaults to --repo-root. "
            "All results/, work/, and .snakemake/ state will be created here."
        ),
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help=(
            "Repository root used to resolve config-relative sample manifests and "
            "assembly paths. Relative paths are resolved against the current directory."
        ),
    )
    parser.add_argument(
        "-S",
        "--snakemake-args",
        default="",
        help=(
            'Additional Snakemake arguments. Use --snakemake-args="--dry-run" '
            'or -S=--dry-run for option-like values.'
        ),
    )
    return parser.parse_args()


def _is_container() -> bool:
    return SNAKEFILE.is_file()


def _local_pipeline_dir() -> Path:
    return Path(__file__).resolve().parent


def _pipeline_dir() -> Path:
    return PIPELINE_DIR if _is_container() else _local_pipeline_dir()


def _snakefile_path() -> Path:
    return SNAKEFILE if _is_container() else _local_pipeline_dir() / "Snakefile"


def resolve_path(path_text: str, base_dir: Path) -> Path:
    path = Path(path_text).expanduser()
    return path if path.is_absolute() else (base_dir / path).resolve()


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Config file {path} does not contain a YAML mapping.")
    return data


def load_samples_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        missing = [column for column in REQUIRED_SAMPLE_COLUMNS if column not in fieldnames]
        if missing:
            raise ValueError(
                f"Sample manifest {path} is missing required columns: {', '.join(missing)}"
            )
        return [dict(row) for row in reader]


def collect_input_paths(
    config_path: Path,
    config_obj: dict,
    repo_root: Path,
) -> tuple[Path, dict[str, InputPathSpec]]:
    input_paths: dict[str, InputPathSpec] = {"config": InputPathSpec(config_path)}
    samples_path = resolve_path(str(config_obj.get("samples", "config/samples.tsv")), repo_root)
    input_paths["samples"] = InputPathSpec(samples_path)

    if samples_path.is_file():
        for row in load_samples_manifest(samples_path):
            sample_id = row["sample_id"].strip() or "unknown"
            raw_assembly_path = row["assembly_fasta"]
            assembly_path = resolve_path(raw_assembly_path.strip(), repo_root)
            input_paths[f"assembly:{sample_id}"] = InputPathSpec(
                path=assembly_path,
                raw_text=raw_assembly_path,
            )

    return samples_path, input_paths


def inferred_mount_root(path: Path) -> Path:
    parts = path.parts
    if len(parts) >= 3:
        return Path(parts[0], parts[1], parts[2])
    if len(parts) >= 2:
        return Path(parts[0], parts[1])
    return Path(path.anchor or "/")


def nearest_existing_ancestor(path: Path) -> Path | None:
    current = path
    while True:
        if current.exists():
            return current
        if current == current.parent:
            return None
        current = current.parent


def resolved_symlink_target(path: Path) -> Path | None:
    if not path.is_symlink():
        return None
    target = path.readlink()
    if target.is_absolute():
        return target
    return (path.parent / target).resolve(strict=False)


def missing_path_diagnosis(path: Path, raw_text: str | None = None) -> str:
    if raw_text is not None and raw_text != raw_text.strip():
        return "manifest value contains leading or trailing whitespace"

    symlink_target = resolved_symlink_target(path)
    if symlink_target is not None:
        target_mount_root = inferred_mount_root(symlink_target)
        if not symlink_target.exists():
            if not target_mount_root.exists():
                return (
                    f"symlink target {symlink_target} is not visible inside the container "
                    f"(missing mount root {target_mount_root})"
                )
            return (
                f"symlink target {symlink_target} is missing or inaccessible inside a visible tree"
            )

    if not path.is_absolute():
        return "relative path could not be resolved to an existing file"

    mount_root = inferred_mount_root(path)

    if not mount_root.exists():
        return f"mount root {mount_root} is not visible inside the container"

    existing_ancestor = nearest_existing_ancestor(path)
    if existing_ancestor is None:
        return "path is not accessible inside the container"

    if existing_ancestor == path.parent:
        return f"path is inside a visible tree, but the file does not exist under {path.parent}"

    return (
        "path is inside a visible tree, but an intermediate directory is missing under "
        f"{existing_ancestor}"
    )


def suggested_bind_dir(path: Path) -> Path:
    if path.exists():
        return path.parent if path.suffix else path

    symlink_target = resolved_symlink_target(path)
    if symlink_target is not None:
        if symlink_target.exists():
            return symlink_target.parent if symlink_target.suffix else symlink_target
        target_mount_root = inferred_mount_root(symlink_target)
        if not target_mount_root.exists():
            return target_mount_root
        existing_target_ancestor = nearest_existing_ancestor(symlink_target)
        if existing_target_ancestor is not None:
            return (
                existing_target_ancestor
                if existing_target_ancestor.is_dir()
                else existing_target_ancestor.parent
            )

    if path.is_absolute():
        mount_root = inferred_mount_root(path)
        if not mount_root.exists():
            return mount_root

    existing_ancestor = nearest_existing_ancestor(path)
    if existing_ancestor is not None:
        return existing_ancestor if existing_ancestor.is_dir() else existing_ancestor.parent

    return path.parent if path.suffix else path


def format_bind_suggestion(
    *,
    image_name: str,
    bind_dirs: list[Path],
    config_path: Path,
    repo_root: Path,
    workdir: Path,
    threads: int,
) -> str:
    bind_flags = " ".join(f"-B {path}" for path in bind_dirs)
    command = (
        f"apptainer run {bind_flags} {image_name} "
        f"-c {config_path} --repo-root {repo_root} --directory {workdir} -t {threads}"
    )
    return command


def show_bind_suggestions(
    *,
    input_paths: dict[str, InputPathSpec],
    config_path: Path,
    repo_root: Path,
    workdir: Path,
    threads: int,
) -> None:
    bind_dirs = {repo_root, workdir, config_path.parent}
    for spec in input_paths.values():
        path = spec.path
        bind_dirs.add(suggested_bind_dir(path))
    ordered_dirs = sorted(path.resolve() for path in bind_dirs)

    print("\nSuggested bind mounts:", file=sys.stderr)
    print(
        "  " + format_bind_suggestion(
            image_name="busco2aster.sif",
            bind_dirs=ordered_dirs,
            config_path=config_path,
            repo_root=repo_root,
            workdir=workdir,
            threads=threads,
        ),
        file=sys.stderr,
    )
    print(
        "  " + format_bind_suggestion(
            image_name="busco2aster.sif",
            bind_dirs=ordered_dirs,
            config_path=config_path,
            repo_root=repo_root,
            workdir=workdir,
            threads=threads,
        ).replace("apptainer run", "singularity run"),
        file=sys.stderr,
    )


def prepare_runtime_workdir(workdir: Path, pipeline_dir: Path) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    for name in ("config", "scripts", "reports"):
        destination = workdir / name
        source = pipeline_dir / name
        if destination.exists() or destination.is_symlink():
            continue
        destination.symlink_to(source.resolve(), target_is_directory=True)


def build_snakemake_command(
    *,
    config_path: Path,
    repo_root: Path,
    workdir: Path,
    threads: int,
    target: str | None,
    snakemake_args: str,
) -> list[str]:
    command = [
        "snakemake",
        "--snakefile",
        _snakefile_path().as_posix(),
        "--directory",
        workdir.as_posix(),
        "--configfile",
        config_path.as_posix(),
        "--cores",
        str(threads),
        "--use-conda",
        "--conda-frontend",
        "conda",
        "--show-failed-logs",
        "--keep-incomplete",
    ]

    if _is_container():
        command.extend(["--conda-prefix", CONDA_ENVS_PATH.as_posix()])

    command.append(target or "all")

    command.extend(
        [
        "--config",
        f"repo_root={repo_root.as_posix()}",
        ]
    )

    if _is_container():
        command.extend(
            [
                f"iqtree_executable={IQTREE_EXECUTABLE.as_posix()}",
                f"wastral_executable={WASTRAL_EXECUTABLE.as_posix()}",
                f"astral4_executable={ASTRAL4_EXECUTABLE.as_posix()}",
            ]
        )

    if snakemake_args:
        command.extend(shlex.split(snakemake_args))

    return command


def validate_runtime_inputs(
    *,
    config_path: Path,
    input_paths: dict[str, InputPathSpec],
    repo_root: Path,
    workdir: Path,
    threads: int,
) -> bool:
    missing = [(key, spec) for key, spec in input_paths.items() if not spec.path.exists()]
    if not missing:
        return True

    print("ERROR: Required workflow inputs are not accessible:", file=sys.stderr)
    for key, spec in missing:
        print(f"  {key}: {spec.path}", file=sys.stderr)
        print(f"    diagnosis: {missing_path_diagnosis(spec.path, spec.raw_text)}", file=sys.stderr)
    show_bind_suggestions(
        input_paths=input_paths,
        config_path=config_path,
        repo_root=repo_root,
        workdir=workdir,
        threads=threads,
    )
    return False


def run() -> int:
    args = parse_args()
    if args.threads < 1:
        print("ERROR: --threads must be at least 1.", file=sys.stderr)
        return 2

    repo_root = resolve_path(args.repo_root, Path.cwd())
    workdir = resolve_path(args.directory, repo_root) if args.directory else repo_root
    config_path = resolve_path(args.config, repo_root)

    if not config_path.is_file():
        print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
        show_bind_suggestions(
            input_paths={"config": InputPathSpec(config_path)},
            config_path=config_path,
            repo_root=repo_root,
            workdir=workdir,
            threads=args.threads,
        )
        return 1

    try:
        config_obj = load_yaml(config_path)
        _, input_paths = collect_input_paths(config_path, config_obj, repo_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if not validate_runtime_inputs(
        config_path=config_path,
        input_paths=input_paths,
        repo_root=repo_root,
        workdir=workdir,
        threads=args.threads,
    ):
        return 1

    prepare_runtime_workdir(workdir, _pipeline_dir())
    cache_dir = workdir / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    command = build_snakemake_command(
        config_path=config_path,
        repo_root=repo_root,
        workdir=workdir,
        threads=args.threads,
        target=args.target,
        snakemake_args=args.snakemake_args,
    )
    print("Running:", " ".join(shlex.quote(token) for token in command))

    env = os.environ.copy()
    pipeline_dir = _pipeline_dir()
    existing_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        f"{pipeline_dir.as_posix()}:{existing_pythonpath}"
        if existing_pythonpath
        else pipeline_dir.as_posix()
    )
    env["XDG_CACHE_HOME"] = cache_dir.as_posix()

    try:
        subprocess.run(command, check=True, env=env)
    except subprocess.CalledProcessError as exc:
        return exc.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
