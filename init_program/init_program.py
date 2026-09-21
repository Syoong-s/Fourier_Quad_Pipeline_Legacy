#!/usr/bin/env python3
#SBATCH --job-name=f77_init
#SBATCH --partition=cpu
#SBATCH --ntasks=200
#SBATCH --exclusive
#SBATCH --ntasks-per-node=40
#SBATCH --output=init_program_%j.out
#SBATCH --error=init_program_%j.err

"""Build the current C++/F77 dataset layout from FITS/FZ archives.

The Slurm directives above deliberately live in this Python file.  Submit it
directly with ``sbatch init_program.py [options]`` or run the same entry point
with Python, mpirun, or srun.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence
import uuid
import warnings


ARCHIVE_SUFFIX = ".fits.fz"
CCDNUM_KEYWORD = "CCDNUM"
DQ_STEM_FROM = "ood"
DQ_STEM_TO = "ooi"

DEFAULT_SCIENCE_ROOT = "/lustre/home/acct-phyzj/share/DES/g"
DEFAULT_DQ_ROOT = "/lustre/home/acct-phyzj/share/DES/mask_v1/g_mask"
DEFAULT_OUTPUT_ROOT = "/lustre/home/acct-phyzj/share/DES/g_band_v1"

NON_CHIP_BASE_DIRECTORIES = (
    "science",
    "dqmask",
    "stamps",
    "expolists",
    "result",
    "stamps/dat_StarInfo",
    "stamps/svg_StarLocus",
    "stamps/fits_StarP",
    "stamps/fits_PsfSrc",
    "stamps/dat_ExpoInfo",
    "stamps/dat_StarComp",
    "stamps/dat_Rescale",
    "stamps/dat_Pcs",
    "stamps/dat_StarCompV2",
    "astrometry/Head",
    "astrometry/dat_Chk",
)

CHIP_PRODUCT_DIRECTORIES = (
    "stamps/Norm",
    "stamps/cat_Orig",
    "stamps/dat_StarCanInfo",
    "stamps/fits_StarCan",
    "stamps/fits_StarCanN",
    "stamps/fits_StarCanP",
    "stamps/dat_SrcInfo",
    "stamps/fits_Src",
    "stamps/fits_Noise",
    "stamps/fits_SrcP",
    "stamps/dat_PsfFit",
    "stamps/fits_PsfLocal",
    "stamps/dat_Shear",
    "stamps/dat_StarXY",
    "stamps/fits_PsfResi",
    "astrometry/dat_Astro",
)

COMPRESSED_HEADER_KEYS = {
    "XTENSION",
    "PCOUNT",
    "GCOUNT",
    "TFIELDS",
    "ZIMAGE",
    "ZTENSION",
    "ZBITPIX",
    "ZNAXIS",
    "ZNAXIS1",
    "ZNAXIS2",
    "ZTILE1",
    "ZTILE2",
    "ZCMPTYPE",
    "ZMASKCMP",
    "ZQUANTIZ",
    "ZDITHER0",
    "ZSIMPLE",
    "ZEXTEND",
}


# ==========================================
# Class: Store normalized initializer configuration
# Method: Keep runtime paths and policies immutable after validation
# ==========================================
@dataclass(frozen=True)
class Config:
    science_root: Path
    dq_root: Path
    output_root: Path
    target: str
    prefix: str
    contains: tuple[str, ...]
    existing: str
    f77_max_path: int


# ==========================================
# Class: Describe one archive extraction task
# Method: Pair a stable task index with kind and immutable source path
# ==========================================
@dataclass(frozen=True)
class Task:
    index: int
    kind: str
    source: str


# ==========================================
# Class: Describe one planned output image
# Method: Retain its source HDU, expected identity, and staged/final paths
# ==========================================
@dataclass(frozen=True)
class PlannedImage:
    hdu_index: int
    output_number: int
    check_ccdnum: bool
    final_path: str
    staged_path: str


# ==========================================
# Class: Return a serializable per-archive extraction result
# Method: Carry successful paths, recovery state, skips, and diagnostics
# ==========================================
@dataclass
class ExtractionResult:
    success: bool = False
    resumed: bool = False
    output_paths: list[str] | None = None
    error: str = ""
    skipped_hdus: int = 0

    # ==========================================
    # Function: Normalize the optional path collection after construction
    # Method: Ensure every result owns a mutable output-path list
    # ==========================================
    def __post_init__(self) -> None:
        if self.output_paths is None:
            self.output_paths = []


# ==========================================
# Function: Build the command-line interface
# Method: Mirror the live C++ initializer while adding one safe Slurm control
# ==========================================
def parse_arguments(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Initialize the exposure-first C++/F77 data tree directly from "
            "read-only science and DQ FITS/FZ repositories."
        )
    )
    parser.add_argument("--science-root", default=DEFAULT_SCIENCE_ROOT)
    parser.add_argument("--dq-root", default=DEFAULT_DQ_ROOT)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--target", default="gband")
    parser.add_argument("--prefix", default="c4d_")
    parser.add_argument(
        "--contains",
        action="append",
        default=None,
        metavar="TEXT",
        help="repeatable basename token; an archive matching any token is selected",
    )
    parser.add_argument(
        "--existing",
        choices=("fail", "resume", "overwrite"),
        default="fail",
    )
    parser.add_argument(
        "--f77-max-path",
        type=int,
        default=512,
        help="maximum generated F77 path length; zero disables the length check",
    )
    parser.add_argument(
        "--no-srun",
        action="store_true",
        help="run serially in the sbatch process instead of launching its allocation",
    )
    parser.add_argument("--mpi-worker", action="store_true", help=argparse.SUPPRESS)
    return parser.parse_args(argv)


# ==========================================
# Function: Normalize and validate the initializer configuration
# Method: Resolve paths and reject unsafe names, limits, tokens, and whitespace
# ==========================================
def normalize_config(arguments: argparse.Namespace) -> Config:
    target = arguments.target
    if not target or target in {".", ".."} or "/" in target or "\\" in target:
        raise ValueError("--target must be one non-empty directory name")
    if not arguments.prefix:
        raise ValueError("--prefix must not be empty")
    if arguments.f77_max_path < 0:
        raise ValueError("--f77-max-path must be zero or positive")

    tokens = tuple(arguments.contains if arguments.contains is not None else ("v1",))
    if any(not token for token in tokens):
        raise ValueError("--contains values must not be empty")

    config = Config(
        science_root=Path(arguments.science_root).expanduser().resolve(),
        dq_root=Path(arguments.dq_root).expanduser().resolve(),
        output_root=Path(arguments.output_root).expanduser().resolve(),
        target=target,
        prefix=arguments.prefix,
        contains=tokens,
        existing=arguments.existing,
        f77_max_path=arguments.f77_max_path,
    )
    validate_pipeline_path(config.output_root / config.target, 0)
    return config


# ==========================================
# Function: Detect a top-level Slurm batch invocation
# Method: Distinguish the batch shell from an srun step or explicit worker
# ==========================================
def should_launch_srun(arguments: argparse.Namespace) -> bool:
    if arguments.no_srun or arguments.mpi_worker or "SLURM_JOB_ID" not in os.environ:
        return False
    step_id = os.environ.get("SLURM_STEP_ID")
    return step_id is None or step_id in {"batch", "extern"}


# ==========================================
# Function: Launch all allocated Slurm tasks into this same Python file
# Method: Use an argument vector and hidden worker flag without shell expansion
# ==========================================
def launch_srun(original_argv: Sequence[str]) -> int:
    srun = shutil.which("srun")
    if srun is None:
        raise RuntimeError("SLURM_JOB_ID is set but srun is not available")
    task_count = os.environ.get("SLURM_NTASKS", "1")
    command = [
        srun,
        "--ntasks",
        task_count,
        "--kill-on-bad-exit=1",
        sys.executable,
        str(Path(__file__).resolve()),
        *original_argv,
        "--mpi-worker",
    ]
    print(f"Launching {task_count} Slurm tasks with {sys.executable}", flush=True)
    return subprocess.run(command, check=False).returncode


# ==========================================
# Function: Load the MPI communicator after Slurm launch decisions
# Method: Keep help and coordinator startup independent of heavy dependencies
# ==========================================
def get_communicator() -> Any:
    try:
        from mpi4py import MPI
    except ImportError as exception:
        raise RuntimeError(
            "mpi4py is required; install it in the selected Python environment"
        ) from exception
    return MPI.COMM_WORLD


# ==========================================
# Function: Test whether one resolved path is inside another
# Method: Use pathlib ancestry rather than fragile string-prefix comparison
# ==========================================
def path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


# ==========================================
# Function: Validate one generated path for legacy F77 consumption
# Method: Reject whitespace and enforce the configured fixed-character bound
# ==========================================
def validate_pipeline_path(path: Path, maximum_length: int) -> None:
    text = str(path)
    if any(character.isspace() for character in text):
        raise ValueError(f"pipeline path contains whitespace: {text}")
    if maximum_length > 0 and len(text) > maximum_length:
        raise ValueError(f"pipeline path exceeds the configured F77 limit: {text}")


# ==========================================
# Function: Return an archive exposure stem
# Method: Remove only the exact configured FITS/FZ suffix
# ==========================================
def archive_stem(path: Path) -> str:
    if not path.name.endswith(ARCHIVE_SUFFIX) or len(path.name) <= len(ARCHIVE_SUFFIX):
        raise ValueError(f"archive does not end in {ARCHIVE_SUFFIX}: {path}")
    return path.name[: -len(ARCHIVE_SUFFIX)]


# ==========================================
# Function: Return the DQ exposure stem consumed by the pipeline
# Method: Remove the archive suffix and replace every ood token with ooi
# ==========================================
def dq_output_stem(path: Path) -> str:
    return archive_stem(path).replace(DQ_STEM_FROM, DQ_STEM_TO)


# ==========================================
# Function: Discover matching archives below one source root
# Method: Apply suffix, prefix, and OR-token filters then sort absolute paths
# ==========================================
def discover_archives(root: Path, prefix: str, tokens: Sequence[str]) -> list[Path]:
    if not root.is_dir():
        raise ValueError(f"archive root is not a directory: {root}")
    matches = []
    for path in root.rglob(f"*{ARCHIVE_SUFFIX}"):
        if not path.is_file() or not path.name.startswith(prefix):
            continue
        if tokens and not any(token in path.name for token in tokens):
            continue
        matches.append(path.resolve())
    return sorted(matches)


# ==========================================
# Function: Reject archives that collide on one output exposure stem
# Method: Compare science stems directly and mapped DQ stems independently
# ==========================================
def validate_unique_stems(paths: Iterable[Path], kind: str) -> None:
    seen: dict[str, Path] = {}
    for path in paths:
        stem = archive_stem(path) if kind == "science" else dq_output_stem(path)
        if stem in seen:
            raise ValueError(
                f"duplicate {kind} output exposure stem {stem}: {seen[stem]} and {path}"
            )
        seen[stem] = path


# ==========================================
# Function: Create the complete shared pipeline base-directory contract
# Method: Materialize every live C++ Standard/Lite base directory idempotently
# ==========================================
def create_pipeline_directories(target_root: Path) -> None:
    for relative_path in (*NON_CHIP_BASE_DIRECTORIES, *CHIP_PRODUCT_DIRECTORIES):
        (target_root / relative_path).mkdir(parents=True, exist_ok=True)


# ==========================================
# Function: Create every chip-product directory for one exposure
# Method: Apply the shared live layout after deterministic list publication
# ==========================================
def create_exposure_directories(target_root: Path, exposures: Iterable[str]) -> None:
    for exposure in sorted(set(exposures)):
        for relative_path in CHIP_PRODUCT_DIRECTORIES:
            (target_root / relative_path / exposure).mkdir(parents=True, exist_ok=True)


# ==========================================
# Function: Remove compressed-table structural cards from one FITS header
# Method: Preserve science metadata while allowing a normal primary image write
# ==========================================
def clean_header(header: Any) -> Any:
    cleaned = header.copy()
    for key in list(cleaned.keys()):
        if (
            key in COMPRESSED_HEADER_KEYS
            or key.startswith("ZNAME")
            or key.startswith("ZVAL")
            or key.startswith("TFORM")
            or key.startswith("TTYPE")
        ):
            del cleaned[key]
    return cleaned


# ==========================================
# Function: Plan every extractable image HDU in one archive
# Method: Use occurrence numbering for science and physical CCDNUM for DQ
# ==========================================
def plan_images(
    images: Any,
    source: Path,
    kind: str,
    final_directory: Path,
    staging_directory: Path,
) -> list[PlannedImage]:
    output_stem = archive_stem(source) if kind == "science" else dq_output_stem(source)
    science_index = 0
    output_names: set[str] = set()
    plans = []
    for hdu_index, hdu in enumerate(images):
        if int(hdu.header.get("NAXIS", 0)) != 2:
            continue
        if kind == "science":
            science_index += 1
            output_number = science_index
            check_ccdnum = False
        else:
            try:
                output_number = int(hdu.header[CCDNUM_KEYWORD])
            except (KeyError, TypeError, ValueError):
                continue
            check_ccdnum = True
        filename = f"{output_stem}_{output_number}.fits"
        if filename in output_names:
            continue
        output_names.add(filename)
        plans.append(
            PlannedImage(
                hdu_index=hdu_index,
                output_number=output_number,
                check_ccdnum=check_ccdnum,
                final_path=str(final_directory / filename),
                staged_path=str(staging_directory / filename),
            )
        )
    if not plans:
        raise RuntimeError("archive contains no extractable two-dimensional image HDUs")
    return plans


# ==========================================
# Function: Validate one existing uncompressed output image
# Method: Require a two-dimensional primary image and the planned DQ CCDNUM
# ==========================================
def validate_existing_output(plan: PlannedImage) -> tuple[bool, str]:
    from astropy.io import fits

    path = Path(plan.final_path)
    try:
        with fits.open(path, mode="readonly", memmap=True) as images:
            primary = images[0]
            if int(primary.header.get("NAXIS", 0)) != 2 or "ZIMAGE" in primary.header:
                return False, f"existing output is not an uncompressed 2D FITS image: {path}"
            if plan.check_ccdnum:
                ccdnum = int(primary.header[CCDNUM_KEYWORD])
                if ccdnum != plan.output_number:
                    return False, f"existing DQ output has the wrong CCDNUM: {path}"
    except Exception as exception:
        return False, f"cannot validate existing output {path}: {exception}"
    return True, ""


# ==========================================
# Function: Write one logical FITS image as an uncompressed primary image
# Method: Clean compression cards and publish first into the task staging tree
# ==========================================
def write_staged_image(hdu: Any, staged_path: Path) -> None:
    from astropy.io import fits
    from astropy.utils.exceptions import AstropyWarning

    staged_path.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        fits.writeto(
            staged_path,
            hdu.data,
            clean_header(hdu.header),
            overwrite=True,
            output_verify="fix",
        )


# ==========================================
# Function: Extract one archive with staged atomic output publication
# Method: Apply fail/resume/overwrite policy and continue past bad HDUs
# ==========================================
def extract_archive(
    task: Task,
    target_root: Path,
    staging_root: Path,
    existing_policy: str,
    rank: int,
) -> ExtractionResult:
    from astropy.io import fits

    source = Path(task.source)
    exposure = archive_stem(source) if task.kind == "science" else dq_output_stem(source)
    final_directory = target_root / task.kind / exposure
    task_staging = staging_root / f"rank_{rank}" / f"task_{task.index}"
    result = ExtractionResult()
    try:
        with fits.open(source, mode="readonly", memmap=True, lazy_load_hdus=True) as images:
            plans = plan_images(images, source, task.kind, final_directory, task_staging)
            outputs_exist = any(Path(plan.final_path).exists() for plan in plans)
            if existing_policy == "fail" and outputs_exist:
                result.error = f"output already exists: {final_directory}"
                return result

            to_extract = []
            for plan in plans:
                final_path = Path(plan.final_path)
                if final_path.exists() and existing_policy == "resume":
                    valid, _ = validate_existing_output(plan)
                    if valid:
                        result.output_paths.append(str(final_path.resolve()))
                        continue
                to_extract.append(plan)

            if not to_extract:
                result.success = True
                result.resumed = True
                return result

            successful_plans = []
            for plan in to_extract:
                try:
                    write_staged_image(images[plan.hdu_index], Path(plan.staged_path))
                    successful_plans.append(plan)
                except Exception as exception:
                    result.skipped_hdus += 1
                    Path(plan.staged_path).unlink(missing_ok=True)
                    print(
                        f"[rank {rank}] skipped {source} HDU {plan.hdu_index + 1}: {exception}",
                        file=sys.stderr,
                        flush=True,
                    )

        if not successful_plans:
            result.error = "all planned HDUs failed extraction"
            return result

        final_directory.mkdir(parents=True, exist_ok=True)
        for plan in successful_plans:
            final_path = Path(plan.final_path)
            os.replace(plan.staged_path, final_path)
            result.output_paths.append(str(final_path.resolve()))
        result.output_paths.sort()
        result.success = True
        return result
    except Exception as exception:
        result.error = str(exception)
        return result
    finally:
        shutil.rmtree(task_staging, ignore_errors=True)


# ==========================================
# Function: Publish one small text or JSON file atomically
# Method: Flush a sibling temporary file then replace the destination
# ==========================================
def write_atomic(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write(content)
            stream.flush()
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


# ==========================================
# Function: Publish one deterministic science-chip list per exposure
# Method: Sort successful output paths and write only complete task results
# ==========================================
def publish_exposure_lists(
    config: Config,
    target_root: Path,
    tasks: Sequence[Task],
    results: Sequence[ExtractionResult],
) -> None:
    for task, result in zip(tasks, results, strict=True):
        if task.kind != "science" or not result.success or not result.output_paths:
            continue
        source = Path(task.source)
        list_path = (target_root / "expolists" / f"{archive_stem(source)}.list").resolve()
        validate_pipeline_path(list_path, config.f77_max_path)
        paths = sorted(Path(path).resolve() for path in result.output_paths)
        for path in paths:
            validate_pipeline_path(path, config.f77_max_path)
        write_atomic(list_path, "".join(f"{path}\n" for path in paths))


# ==========================================
# Function: Publish the top exposure list and flat compatibility FITS list
# Method: Scan sorted expolists without rotating the first record to the end
# ==========================================
def publish_pipeline_lists(config: Config, target_root: Path) -> list[str]:
    exposure_lists = sorted((target_root / "expolists").glob("*.list"))
    top_lines = []
    fits_lines = []
    exposures = []
    for list_path in exposure_lists:
        absolute_list = list_path.resolve()
        validate_pipeline_path(absolute_list, config.f77_max_path)
        chip_paths = [
            line.strip()
            for line in list_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        for chip_path in chip_paths:
            validate_pipeline_path(Path(chip_path), config.f77_max_path)
            fits_lines.append(f"{chip_path}\n")
        top_lines.append(f'"{absolute_list}"     {len(chip_paths)}\n')
        exposures.append(list_path.stem)

    top_path = (config.output_root / f"expo_{config.target}.list").resolve()
    fits_path = (config.output_root / f"fits_{config.target}.list").resolve()
    validate_pipeline_path(top_path, config.f77_max_path)
    validate_pipeline_path(fits_path, config.f77_max_path)
    write_atomic(top_path, "".join(top_lines))
    write_atomic(fits_path, "".join(fits_lines))
    return exposures


# ==========================================
# Function: Build a durable initialization manifest
# Method: Record live naming/order policies, counts, failures, and provenance
# ==========================================
def build_manifest(
    config: Config,
    tasks: Sequence[Task],
    results: Sequence[ExtractionResult],
    lists_published: bool,
    exposure_directories_created: bool,
    final_error: str,
) -> dict[str, Any]:
    task_results = list(zip(tasks, results, strict=True))
    science_results = [result for task, result in task_results if task.kind == "science"]
    dq_results = [result for task, result in task_results if task.kind == "dqmask"]
    failed_sources = [task.source for task, result in task_results if not result.success]
    partial_sources = [
        task.source
        for task, result in task_results
        if result.success and result.skipped_hdus > 0
    ]
    status = "success"
    if not lists_published or not exposure_directories_created:
        status = "failed"
    elif failed_sources or partial_sources:
        status = "partial"
    return {
        "schema_version": 2,
        "status": status,
        "implementation": "init_program.py",
        "direct_source_read": True,
        "copy_staging": False,
        "exposure_order": "corrected_lexical_no_rotation",
        "science_numbering": "two_dimensional_hdu_occurrence",
        "dq_numbering": CCDNUM_KEYWORD,
        "science_root": str(config.science_root),
        "dq_root": str(config.dq_root),
        "output_root": str(config.output_root),
        "target": config.target,
        "filename_prefix": config.prefix,
        "filename_tokens": list(config.contains),
        "existing_policy": config.existing,
        "f77_max_path": config.f77_max_path,
        "science_sources": len(science_results),
        "dq_sources": len(dq_results),
        "science_images": sum(
            len(result.output_paths) for result in science_results if result.success
        ),
        "dq_images": sum(
            len(result.output_paths) for result in dq_results if result.success
        ),
        "resumed_sources": sum(result.resumed for result in results),
        "lists_published": lists_published,
        "exposure_directories_created": exposure_directories_created,
        "error": final_error,
        "failed_sources": failed_sources,
        "skipped_hdus": sum(result.skipped_hdus for result in results),
        "partial_sources": partial_sources,
    }


# ==========================================
# Function: Run discovery, distributed extraction, and metadata publication
# Method: Broadcast one rank-zero plan and assign tasks by stable rank stride
# ==========================================
def run_initializer(config: Config, communicator: Any) -> int:
    rank = communicator.Get_rank()
    process_count = communicator.Get_size()
    target_root = config.output_root / config.target
    setup_payload: dict[str, Any] | None = None

    if rank == 0:
        try:
            if path_is_within(target_root, config.science_root) or path_is_within(
                target_root, config.dq_root
            ):
                raise ValueError("target output must not be inside either source repository")
            science_sources = discover_archives(
                config.science_root, config.prefix, config.contains
            )
            dq_sources = discover_archives(config.dq_root, config.prefix, config.contains)
            if not science_sources:
                raise ValueError(f"no matching science {ARCHIVE_SUFFIX} archives were found")
            if not dq_sources:
                raise ValueError(f"no matching DQ {ARCHIVE_SUFFIX} archives were found")
            validate_unique_stems(science_sources, "science")
            validate_unique_stems(dq_sources, "dqmask")
            create_pipeline_directories(target_root)
            run_token = f"run_{time.time_ns()}_{os.getpid()}"
            staging_root = target_root / ".fq_init_tmp" / run_token
            staging_root.mkdir(parents=True, exist_ok=False)
            tasks = [
                Task(index=index, kind="science", source=str(path))
                for index, path in enumerate(science_sources)
            ]
            dq_start_index = len(tasks)
            tasks.extend(
                Task(index=dq_start_index + offset, kind="dqmask", source=str(path))
                for offset, path in enumerate(dq_sources)
            )
            setup_payload = {
                "ok": True,
                "tasks": [asdict(task) for task in tasks],
                "staging_root": str(staging_root),
            }
        except Exception as exception:
            setup_payload = {"ok": False, "error": str(exception)}

    setup_payload = communicator.bcast(setup_payload, root=0)
    if not setup_payload["ok"]:
        if rank == 0:
            print(f"Initializer setup failed: {setup_payload['error']}", file=sys.stderr)
        return 1

    tasks = [Task(**item) for item in setup_payload["tasks"]]
    staging_root = Path(setup_payload["staging_root"])
    local_results = []
    for task_index in range(rank, len(tasks), process_count):
        task = tasks[task_index]
        print(
            f"[rank {rank}] {task_index + 1}/{len(tasks)} {task.kind}: {task.source}",
            flush=True,
        )
        local_results.append(
            (
                task.index,
                asdict(
                    extract_archive(
                        task, target_root, staging_root, config.existing, rank
                    )
                ),
            )
        )

    gathered = communicator.gather(local_results, root=0)
    communicator.Barrier()
    final_status = 0
    if rank == 0:
        results = [ExtractionResult(error="task was not executed") for _ in tasks]
        for rank_results in gathered:
            for task_index, result_data in rank_results:
                results[task_index] = ExtractionResult(**result_data)

        lists_published = False
        exposure_directories_created = False
        final_errors = []
        failed_results = [result for result in results if not result.success]
        if failed_results:
            final_errors.append(f"{len(failed_results)} source archive(s) failed")
            final_status = 1
        try:
            publish_exposure_lists(config, target_root, tasks, results)
            exposures = publish_pipeline_lists(config, target_root)
            lists_published = True
            create_exposure_directories(target_root, exposures)
            exposure_directories_created = True
        except Exception as exception:
            final_status = 1
            final_errors.append(str(exception))

        final_error = "; ".join(final_errors)
        manifest_path = config.output_root / f"init_{config.target}_manifest.json"
        try:
            manifest = build_manifest(
                config,
                tasks,
                results,
                lists_published,
                exposure_directories_created,
                final_error,
            )
            manifest_text = json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
            write_atomic(manifest_path, manifest_text)
        except Exception as exception:
            final_status = 1
            final_error = f"{final_error}; {exception}" if final_error else str(exception)

        shutil.rmtree(staging_root, ignore_errors=True)
        staging_parent = target_root / ".fq_init_tmp"
        try:
            staging_parent.rmdir()
        except OSError:
            pass

        if final_status == 0:
            science_count = sum(task.kind == "science" for task in tasks)
            dq_count = sum(task.kind == "dqmask" for task in tasks)
            print(
                f"Initialization complete: {science_count} science archives, "
                f"{dq_count} DQ archives; exposure list: "
                f"{config.output_root / ('expo_' + config.target + '.list')}",
                flush=True,
            )
        else:
            print(f"Initialization incomplete: {final_error}", file=sys.stderr, flush=True)

    final_status = communicator.bcast(final_status, root=0)
    communicator.Barrier()
    return final_status


# ==========================================
# Function: Coordinate CLI, Slurm relaunch, MPI, and initializer execution
# Method: Keep every supported execution mode behind one Python entry point
# ==========================================
def main(argv: Sequence[str] | None = None) -> int:
    original_argv = list(sys.argv[1:] if argv is None else argv)
    try:
        arguments = parse_arguments(original_argv)
        if should_launch_srun(arguments):
            return launch_srun(original_argv)
        config = normalize_config(arguments)
        communicator = get_communicator()
        return run_initializer(config, communicator)
    except Exception as exception:
        print(f"Initializer error: {exception}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
