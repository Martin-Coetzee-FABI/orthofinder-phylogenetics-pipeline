#!/usr/bin/env python3
"""
05_run_iqtree.py — Run IQ-TREE3 species tree

IQ-TREE3 reads trimmed alignment files directly from a directory — no
concatenated supermatrix file is needed.

Two modes (--mode):
  partition  (default) : IQ-TREE3 uses -p <dir>, ModelFinder selects the
                         best model per gene, models may be merged.
  single               : IQ-TREE3 uses -s <dir>, one model for the whole
                         concatenated alignment.

Both modes apply ultrafast bootstrap (-B) AND SH-aLRT (-alrt) branch support.
"""

import argparse
import logging
import os
import re
import shutil
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


def log_header(logger: logging.Logger, args: argparse.Namespace, start_time: datetime,
               iqtree_ver: str) -> None:
    logger.info("=" * 54)
    logger.info("Script:      05_run_iqtree.py")
    logger.info(f"Started:     {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Parameters:  --trimmed_dir {args.trimmed_dir}")
    logger.info(f"             --output_dir {args.output_dir}")
    logger.info(f"             --mode {args.mode}")
    logger.info(f"             --model {args.model}")
    logger.info(f"             --bootstrap {args.bootstrap}")
    logger.info(f"             --alrt {args.alrt}")
    logger.info(f"             --threads {args.threads}")
    logger.info(f"             --prefix {args.prefix}")
    logger.info(f"Tool versions:")
    logger.info(f"  IQ-TREE:   {iqtree_ver}")
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


def find_iqtree() -> str | None:
    candidates = [
        "iqtree3", "iqtree",
        "/usr/local/bin/iqtree3", "/opt/homebrew/bin/iqtree3",
        "/usr/local/bin/iqtree", "/opt/homebrew/bin/iqtree",
    ]
    for c in candidates:
        if shutil.which(c) or (os.path.isfile(c) and os.access(c, os.X_OK)):
            return c
    return None


def get_version(exe: str) -> str:
    try:
        result = subprocess.run(
            [exe, "--version"], capture_output=True, text=True, timeout=15
        )
        out = (result.stdout or result.stderr or "").strip()
        return out.split("\n")[0] if out else "unknown"
    except Exception:
        return "unknown"


def make_unbuffered_cmd(cmd: list[str]) -> tuple[list[str], str]:
    """
    Wrap cmd so the child process uses line-buffered stdout.
    IQ-TREE buffers output when writing to a pipe (not a terminal).
    stdbuf/unbuffer forces it to flush each line immediately.
    Returns (wrapped_cmd, method_used).
    """
    for wrapper in ["stdbuf", "gstdbuf"]:
        exe = shutil.which(wrapper)
        if exe:
            return [exe, "-oL", "-eL"] + cmd, f"{wrapper} -oL -eL"
    unbuffer = shutil.which("unbuffer")
    if unbuffer:
        return [unbuffer] + cmd, "unbuffer"
    return cmd, "none (output may be delayed until IQ-TREE internal buffer fills)"


def parse_iqtree_report(iqtree_report: Path, logger: logging.Logger) -> dict:
    """Extract model, log-likelihood, and parsimony-informative sites from .iqtree file."""
    stats = {
        "model": "unknown",
        "log_likelihood": "unknown",
        "parsimony_sites": "unknown",
    }
    if not iqtree_report.exists():
        return stats
    try:
        text = iqtree_report.read_text()

        m = re.search(r"Best-fit model.*?:\s*(\S+)", text)
        if m:
            stats["model"] = m.group(1)

        m = re.search(r"Log-likelihood of the tree:\s*([\-\d\.]+)", text)
        if m:
            stats["log_likelihood"] = m.group(1)

        m = re.search(r"Number of parsimony.informative sites:\s*(\d+)", text, re.IGNORECASE)
        if m:
            stats["parsimony_sites"] = m.group(1)
    except Exception as e:
        logger.warning(f"Could not parse IQ-TREE report: {e}")
    return stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run IQ-TREE3 species tree from a directory of trimmed alignments. "
            "IQ-TREE3 reads all alignment files in the directory directly."
        )
    )
    parser.add_argument("--trimmed_dir", default="alignment_results/trimmed",
                        help="Directory of trimmed .faa alignment files from Step 4 "
                             "(default: alignment_results/trimmed)")
    parser.add_argument("--output_dir", default="iqtree_results",
                        help="IQ-TREE output directory (default: iqtree_results)")
    parser.add_argument("--mode", default="partition", choices=["single", "partition"],
                        help="partition: per-gene model via -p (recommended); "
                             "single: one model for all genes via -s (default: partition)")
    parser.add_argument("--model", default="LG+F+G4",
                        help=(
                            "Substitution model passed to IQ-TREE -m (applies to BOTH modes). "
                            "Examples: LG+F+G4 (fast, fixed model — recommended default), "
                            "MFP (ModelFinder per partition, slow), "
                            "MFP+MERGE (ModelFinder + merge partitions, slowest). "
                            "(default: LG+F+G4)"
                        ))
    parser.add_argument("--bootstrap", type=int, default=1000,
                        help="Ultrafast bootstrap replicates -B (default: 1000)")
    parser.add_argument("--alrt", type=int, default=1000,
                        help="SH-aLRT replicates -alrt (default: 1000)")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="Threads (default: all available)")
    parser.add_argument("--iqtree", default="",
                        help="Path to IQ-TREE3 executable (auto-detect)")
    parser.add_argument("--prefix", default="species_tree",
                        help="Output prefix for IQ-TREE files (default: species_tree)")
    return parser.parse_args()


def main() -> None:
    start_time = datetime.now()
    args = parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    iqtree_exe = args.iqtree or find_iqtree()
    if not iqtree_exe:
        print(
            "ERROR: IQ-TREE3 not found. Install via conda: "
            "conda install -c bioconda iqtree\n"
            "Or specify --iqtree /path/to/iqtree3",
            file=sys.stderr
        )
        sys.exit(1)

    iqtree_ver = get_version(iqtree_exe)
    log_path = str(out_dir / "iqtree_step.log")
    logger = setup_logger(log_path)
    log_header(logger, args, start_time, iqtree_ver)

    trimmed_dir = Path(args.trimmed_dir)
    if not trimmed_dir.exists():
        logger.error(
            f"Trimmed alignment directory not found: {trimmed_dir}\n"
            "Run Step 4 (04_align_trim_concat.py) first."
        )
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    aln_files = sorted(trimmed_dir.glob("*.faa")) + sorted(trimmed_dir.glob("*.fas"))
    if not aln_files:
        logger.error(
            f"No alignment files (.faa or .fas) found in {trimmed_dir}\n"
            "Run Step 4 first or check --trimmed_dir path."
        )
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    logger.info(f"Found {len(aln_files)} trimmed alignment files in {trimmed_dir}")

    prefix_path = out_dir / args.prefix

    if args.mode == "single":
        cmd = [
            iqtree_exe,
            "-s", str(trimmed_dir.resolve()),
            "-m", args.model,
            "-B", str(args.bootstrap),
            "-alrt", str(args.alrt),
            "-T", str(args.threads),
            "--prefix", str(prefix_path),
            "--redo",
        ]
        logger.info(
            f"Mode: single — model: {args.model}, "
            f"-B {args.bootstrap}, -alrt {args.alrt}"
        )
    else:
        cmd = [
            iqtree_exe,
            "-p", str(trimmed_dir.resolve()),
            "-m", args.model,
            "-B", str(args.bootstrap),
            "-alrt", str(args.alrt),
            "-T", str(args.threads),
            "--prefix", str(prefix_path),
            "--redo",
        ]
        logger.info(
            f"Mode: partition — model: {args.model}, "
            f"-B {args.bootstrap}, -alrt {args.alrt}"
        )

    cmd_str = " ".join(cmd)
    cmd_file = out_dir / "iqtree_command.txt"
    with open(cmd_file, "w") as f:
        f.write(f"IQ-TREE version: {iqtree_ver}\n")
        f.write(f"Mode:    {args.mode}\n")
        f.write(f"Run at:  {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Command: {cmd_str}\n")
    logger.info(f"Command written to {cmd_file}")

    cmd_run, stream_method = make_unbuffered_cmd(cmd)
    logger.info(f"Live output method: {stream_method}")
    logger.info(f"Running: {cmd_str}")

    print("")
    print("=" * 62)
    print("  IQ-TREE3 output (live):")
    print("=" * 62)
    sys.stdout.flush()

    iqtree_native_log = str(out_dir / f"{args.prefix}.log")
    try:
        process = subprocess.Popen(
            cmd_run,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        with open(iqtree_native_log, "w") as native_log:
            for line in process.stdout:
                native_log.write(line)
                native_log.flush()
                sys.stdout.write(line)
                sys.stdout.flush()

        process.wait()

        print("=" * 62)
        print("")
        sys.stdout.flush()

        if process.returncode != 0:
            logger.error(f"IQ-TREE exited with return code {process.returncode}")
            logger.error(f"Check {iqtree_native_log} for details")
            log_footer(logger, start_time, "FAILURE")
            sys.exit(1)

    except FileNotFoundError:
        logger.error(f"IQ-TREE executable not found: {iqtree_exe}")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)
    except Exception as e:
        logger.error(f"IQ-TREE run failed: {e}")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    iqtree_report = out_dir / f"{args.prefix}.iqtree"
    stats = parse_iqtree_report(iqtree_report, logger)

    logger.info("IQ-TREE completed successfully.")
    logger.info(f"  Best-fit model(s):           {stats['model']}")
    logger.info(f"  Log-likelihood:              {stats['log_likelihood']}")
    logger.info(f"  Parsimony-informative sites: {stats['parsimony_sites']}")

    treefile = out_dir / f"{args.prefix}.treefile"
    if treefile.exists():
        logger.info(f"  Species tree (Newick):       {treefile}")
        print(f"\nSpecies tree: {treefile.resolve()}")
        print(
            "Branch support: ultrafast bootstrap (values on nodes) + SH-aLRT\n"
            "Visualise in FigTree (https://tree.bio.ed.ac.uk/software/figtree/) "
            "or iTOL (https://itol.embl.de/)"
        )
    else:
        logger.warning(f"Expected treefile not found: {treefile}")

    log_footer(logger, start_time, "SUCCESS")


if __name__ == "__main__":
    main()
