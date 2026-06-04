#!/usr/bin/env python3
"""
02_run_orthofinder.py — Run OrthoFinder3 on passed proteomes

Auto-detects the OrthoFinder executable, runs it on the input directory,
captures stdout/stderr to a log file, and records the exact command used.
"""

import argparse
import logging
import os
import shutil
import signal
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def setup_logger(log_path: str) -> logging.Logger:
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def log_header(logger: logging.Logger, args: argparse.Namespace, start_time: datetime) -> None:
    logger.info("=" * 54)
    logger.info("Script:      02_run_orthofinder.py")
    logger.info(f"Started:     {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Parameters:  --input_dir {args.input_dir}")
    logger.info(f"             --output_dir {args.output_dir}")
    logger.info(f"             --threads {args.threads}")
    logger.info(f"             --search_prog {args.search_prog}  (-S)")
    logger.info(f"             --gene_tree_method {args.gene_tree_method}  (-M)")
    if args.gene_tree_method == "msa":
        logger.info(f"             --msa_prog {args.msa_prog}  (-A)")
    if args.species_tree:
        logger.info(f"             --species_tree {args.species_tree}  (-s)")
    if args.timeout > 0:
        logger.info(f"             --timeout {args.timeout}s")
    if args.orthofinder:
        logger.info(f"             --orthofinder {args.orthofinder}")
    if args.extra_args:
        logger.info(f"             --extra_args {args.extra_args}")
    logger.info("=" * 54)


def log_footer(logger: logging.Logger, start_time: datetime, status: str) -> None:
    end_time = datetime.now()
    duration = end_time - start_time
    minutes, seconds = divmod(int(duration.total_seconds()), 60)
    logger.info("=" * 54)
    logger.info(f"Finished:    {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Duration:    {minutes}m {seconds}s")
    logger.info(f"Exit status: {status}")
    logger.info("=" * 54)


def find_orthofinder() -> str | None:
    candidates = [
        "orthofinder3", "orthofinder",
        "/usr/local/bin/orthofinder3", "/usr/local/bin/orthofinder",
        "/opt/homebrew/bin/orthofinder3", "/opt/homebrew/bin/orthofinder",
        os.path.expanduser("~/miniconda3/bin/orthofinder"),
        os.path.expanduser("~/anaconda3/bin/orthofinder"),
    ]
    for candidate in candidates:
        if shutil.which(candidate) or (os.path.isfile(candidate) and os.access(candidate, os.X_OK)):
            return candidate
    return None


def get_version(executable: str) -> str:
    try:
        result = subprocess.run(
            [executable, "--version"],
            capture_output=True, text=True, timeout=30
        )
        version_line = (result.stdout or result.stderr or "").strip().split("\n")[0]
        return version_line if version_line else "unknown"
    except Exception:
        return "unknown"


def find_results_dir(output_dir: Path) -> Path | None:
    matches = sorted(output_dir.glob("Results_*"))
    return matches[-1] if matches else None


def is_completed_run(results_dir: Path) -> bool:
    """Check if a Results_* directory contains a completed OrthoFinder run."""
    gc_file = results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    return gc_file.exists()


def verify_tool(name: str, logger: logging.Logger) -> bool:
    """Return True if the tool is found in PATH, log an error and return False otherwise."""
    if shutil.which(name):
        return True
    logger.error(
        f"Required tool not found in PATH: {name}\n"
        f"Install via conda: conda install -c bioconda {name}"
    )
    return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run OrthoFinder3 on passed proteomes."
    )
    parser.add_argument("--input_dir", default="qc_results/passed_proteomes",
                        help="Directory of .faa files to analyse (default: qc_results/passed_proteomes)")
    parser.add_argument("--output_dir", default="orthofinder_results",
                        help="Where OrthoFinder writes results (default: orthofinder_results)")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="Number of CPU threads (default: all available)")
    parser.add_argument(
        "--search_prog", default="diamond",
        choices=["diamond", "diamond_ultra_sens", "blast", "blastp", "mmseqs", "blastn"],
        help=(
            "Sequence search program passed to OrthoFinder -S (default: diamond). "
            "diamond=fast; diamond_ultra_sens=sensitive; blast/blastp=classic NCBI BLAST; "
            "mmseqs=ultra-fast (requires MMseqs2 installed)"
        ),
    )
    parser.add_argument(
        "--gene_tree_method", default="msa",
        choices=["dendroblast", "msa"],
        help=(
            "Gene tree inference method passed to OrthoFinder -M (default: msa). "
            "msa=alignment-based, recommended in OrthoFinder3; "
            "dendroblast=fast distance-based (legacy)"
        ),
    )
    parser.add_argument(
        "--msa_prog", default="mafft",
        choices=["mafft", "muscle", "famsa"],
        help=(
            "MSA program for gene trees, passed to OrthoFinder -A "
            "(only used when --gene_tree_method msa; default: mafft). "
            "famsa=fastest; mafft=balanced; muscle=alternative"
        ),
    )
    parser.add_argument(
        "--species_tree", default="",
        help=(
            "Path to a user-supplied ROOTED species tree in Newick format, "
            "passed to OrthoFinder -s. "
            "If provided, OrthoFinder uses this tree instead of inferring one. "
            "Leave blank to let OrthoFinder infer the species tree (default)."
        ),
    )
    parser.add_argument(
        "--timeout", type=int, default=86400,
        help=(
            "Maximum seconds to wait for OrthoFinder before killing it (default: 86400 = 24h). "
            "Set to 0 to disable the timeout."
        ),
    )
    parser.add_argument("--orthofinder", default="",
                        help="Path to OrthoFinder executable (auto-detected if omitted)")
    parser.add_argument("--extra_args", default="",
                        help="Additional OrthoFinder arguments as a quoted string")
    return parser.parse_args()


def main() -> None:
    start_time = datetime.now()
    args = parse_args()

    out_dir = Path(args.output_dir)

    # Log to cwd so it exists before (and survives removal of) out_dir
    log_path = "orthofinder_step.log"
    logger = setup_logger(log_path)
    log_header(logger, args, start_time)

    # ── Cap threads to available CPU count ──────────────────────
    max_cores = os.cpu_count() or 4
    if args.threads > max_cores:
        logger.warning(
            f"Requested threads ({args.threads}) exceeds available cores ({max_cores}). "
            f"Capping to {max_cores} to prevent hangs."
        )
        args.threads = max_cores

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        logger.error(f"Input directory does not exist: {input_dir}")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    faa_files = list(input_dir.glob("*.faa"))
    if not faa_files:
        logger.error(f"No .faa files found in {input_dir}")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)
    logger.info(f"Found {len(faa_files)} .faa files in {input_dir}")

    # ── Verify required external tools before launching ─────────
    # Search tool
    search_exe_map = {
        "diamond": "diamond",
        "diamond_ultra_sens": "diamond",
        "blast": "blastp",
        "blastp": "blastp",
        "mmseqs": "mmseqs",
    }
    search_tool = search_exe_map.get(args.search_prog, args.search_prog)
    if not verify_tool(search_tool, logger):
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    # MSA tool (only needed for msa gene tree method)
    if args.gene_tree_method == "msa":
        if not verify_tool(args.msa_prog, logger):
            log_footer(logger, start_time, "FAILURE")
            sys.exit(1)

    # Skip if a completed run already exists inside out_dir
    if out_dir.exists():
        existing_results = find_results_dir(out_dir)
        if existing_results and is_completed_run(existing_results):
            logger.info(f"Completed OrthoFinder run already exists: {existing_results}")
            logger.info("Skipping OrthoFinder run. Delete the Results_* folder to re-run.")
            print(f"\nOrthoFinder results directory: {existing_results.resolve()}")
            log_footer(logger, start_time, "SKIPPED (existing results)")
            return
        # Incomplete or empty directory — OrthoFinder refuses to write to an
        # existing directory, so remove it and let OrthoFinder re-create it.
        logger.info(
            f"Removing incomplete output directory (no completed Results_* found): {out_dir}"
        )
        shutil.rmtree(out_dir)

    of_exec = args.orthofinder if args.orthofinder else find_orthofinder()
    if not of_exec:
        logger.error(
            "OrthoFinder executable not found. Install it via conda (bioconda channel) "
            "or specify --orthofinder /path/to/orthofinder"
        )
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    version_str = get_version(of_exec)
    logger.info(f"OrthoFinder executable: {of_exec}")
    logger.info(f"OrthoFinder version:    {version_str}")

    # Validate species tree path if supplied
    if args.species_tree:
        stree_path = Path(args.species_tree)
        if not stree_path.exists():
            logger.error(f"Species tree file not found: {stree_path}")
            log_footer(logger, start_time, "FAILURE")
            sys.exit(1)
        if stree_path.stat().st_size == 0:
            logger.error(f"Species tree file is empty: {stree_path}")
            log_footer(logger, start_time, "FAILURE")
            sys.exit(1)
        logger.info(f"User-supplied species tree: {stree_path.resolve()}")

    cmd = [
        of_exec,
        "-f", str(input_dir.resolve()),
        "-o", str(out_dir.resolve()),
        "-t", str(args.threads),
        "-a", str(args.threads),
        "-S", args.search_prog,
        "-M", args.gene_tree_method,
    ]
    if args.gene_tree_method == "msa":
        cmd += ["-A", args.msa_prog]
    if args.species_tree:
        cmd += ["-s", str(Path(args.species_tree).resolve())]
    if args.extra_args:
        import shlex
        cmd.extend(shlex.split(args.extra_args))

    cmd_str = " ".join(cmd)
    logger.info(f"OrthoFinder command (full): {cmd_str}")
    if args.timeout > 0:
        logger.info(f"Timeout: {args.timeout}s ({args.timeout // 3600}h {(args.timeout % 3600) // 60}m)")

    run_log_path = Path("orthofinder_stdout.log")
    try:
        with open(run_log_path, "w") as run_log:
            run_log.write(f"OrthoFinder run started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            run_log.write(f"Command: {cmd_str}\n")
            if args.timeout > 0:
                run_log.write(f"Timeout: {args.timeout}s\n")
            run_log.write("=" * 60 + "\n")
            run_log.flush()

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            timed_out = False
            deadline = None
            if args.timeout > 0:
                import threading

                def _killer():
                    nonlocal timed_out
                    timed_out = True
                    logger.error(
                        f"OrthoFinder timed out after {args.timeout}s "
                        f"— killing process (PID {process.pid})"
                    )
                    try:
                        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                    except Exception:
                        process.kill()

                deadline = threading.Timer(args.timeout, _killer)
                deadline.start()

            try:
                for line in process.stdout:
                    run_log.write(line)
                    run_log.flush()
                    sys.stdout.write(line)
                    sys.stdout.flush()

                process.wait()
            finally:
                if deadline is not None:
                    deadline.cancel()

            end_time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            run_log.write("=" * 60 + "\n")
            run_log.write(f"OrthoFinder run finished: {end_time_str}\n")
            run_log.write(f"Return code: {process.returncode}\n")
            if timed_out:
                run_log.write("TIMED OUT\n")

        if timed_out:
            logger.error(
                f"OrthoFinder was killed after exceeding the {args.timeout}s timeout. "
                f"Check {run_log_path} for partial output. "
                f"Increase OF_TIMEOUT in pipeline_runner.sh if the dataset is large."
            )
            log_footer(logger, start_time, "FAILURE (timeout)")
            sys.exit(1)

        if process.returncode != 0:
            logger.error(f"OrthoFinder exited with return code {process.returncode}")
            logger.error(f"Check {run_log_path} for details")
            log_footer(logger, start_time, "FAILURE")
            sys.exit(1)

    except FileNotFoundError:
        logger.error(f"OrthoFinder executable not found at: {of_exec}")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)
    except Exception as e:
        logger.error(f"OrthoFinder run failed: {e}")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    results_dir = find_results_dir(out_dir)
    if not results_dir:
        logger.error("OrthoFinder completed but Results_* directory not found.")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    # out_dir now exists (created by OrthoFinder) — write bookkeeping files into it
    cmd_file = out_dir / "orthofinder_command.txt"
    with open(cmd_file, "w") as f:
        f.write(f"Executable: {of_exec}\n")
        f.write(f"Version:    {version_str}\n")
        f.write(f"Run at:     {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Command:    {cmd_str}\n")
        if args.timeout > 0:
            f.write(f"Timeout:    {args.timeout}s\n")
    logger.info(f"Command written to {cmd_file}")

    # Move the log files from cwd into out_dir
    for fname in ("orthofinder_step.log", "orthofinder_stdout.log"):
        src = Path(fname)
        if src.exists():
            shutil.move(str(src), str(out_dir / fname))

    logger.info(f"OrthoFinder completed successfully.")
    logger.info(f"Results directory: {results_dir.resolve()}")
    print(f"\nOrthoFinder results directory: {results_dir.resolve()}")

    log_footer(logger, start_time, "SUCCESS")


if __name__ == "__main__":
    main()
