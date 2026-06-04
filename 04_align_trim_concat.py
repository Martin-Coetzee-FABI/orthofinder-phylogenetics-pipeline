#!/usr/bin/env python3
"""
04_align_trim_concat.py — Align and trim per-ortholog FASTAs

MAFFT (align) → trimAl (trim)

The trimmed alignments in trimmed/ are the direct input to IQ-TREE3 (Step 5),
which reads the directory without needing a concatenated supermatrix file.

A fasconcat_ready/ folder is also produced: it contains the same trimmed
alignments renamed to .fas extension so the user can optionally run
FASconCAT-G manually in that directory.
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False


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
               mafft_ver: str, trimal_ver: str) -> None:
    logger.info("=" * 54)
    logger.info("Script:      04_align_trim_concat.py")
    logger.info(f"Started:     {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Parameters:  --gene_fasta_dir {args.gene_fasta_dir}")
    logger.info(f"             --output_dir {args.output_dir}")
    logger.info(f"             --threads {args.threads}")
    logger.info(f"             --trimal_mode {args.trimal_mode}")
    logger.info(f"Tool versions:")
    logger.info(f"  MAFFT:  {mafft_ver}")
    logger.info(f"  trimAl: {trimal_ver}")
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


def find_executable(names: list[str]) -> str | None:
    for name in names:
        found = shutil.which(name)
        if found:
            return found
        if os.path.isfile(name) and os.access(name, os.X_OK):
            return name
    return None


def get_version(exe: str, flag: str = "--version") -> str:
    try:
        result = subprocess.run(
            [exe, flag], capture_output=True, text=True, timeout=15
        )
        out = (result.stdout or result.stderr or "").strip()
        return out.split("\n")[0] if out else "unknown"
    except Exception:
        return "unknown"


def run_command(
    cmd: list[str],
    logger: logging.Logger,
    cmd_log_path: Path,
    timeout: int = 3600,
) -> int:
    """Run a command, log it with timestamp and exit code."""
    cmd_str = " ".join(cmd)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout
        )
        rc = result.returncode
        with open(cmd_log_path, "a") as f:
            f.write(f"[{ts}] CMD: {cmd_str}\n")
            f.write(f"[{ts}] RC:  {rc}\n")
            if result.stderr.strip():
                f.write(f"[{ts}] STDERR: {result.stderr.strip()[:500]}\n")
            f.write("\n")
        if rc != 0:
            logger.error(f"Command failed (rc={rc}): {cmd_str}")
            logger.error(f"stderr: {result.stderr.strip()[:300]}")
        return rc
    except subprocess.TimeoutExpired:
        with open(cmd_log_path, "a") as f:
            f.write(f"[{ts}] CMD: {cmd_str}\n")
            f.write(f"[{ts}] RC:  TIMEOUT\n\n")
        logger.error(f"Command timed out: {cmd_str}")
        return -1
    except Exception as e:
        with open(cmd_log_path, "a") as f:
            f.write(f"[{ts}] CMD: {cmd_str}\n")
            f.write(f"[{ts}] RC:  ERROR ({e})\n\n")
        logger.error(f"Command error: {e}")
        return -1


def alignment_length(fasta_path: Path) -> int:
    """Return alignment length (assumes all sequences equal length)."""
    from Bio import SeqIO
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        return len(str(record.seq))
    return 0


def trimal_flag(mode: str) -> str:
    mapping = {
        "auto": "-automated1",
        "gappyout": "-gappyout",
        "strict": "-strict",
        "strictplus": "-strictplus",
    }
    return mapping.get(mode, "-automated1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align (MAFFT) and trim (trimAl) per-ortholog FASTAs for IQ-TREE3."
    )
    parser.add_argument("--gene_fasta_dir", default="ortholog_results/gene_fastas",
                        help="Directory with per-ortholog .faa files (default: ortholog_results/gene_fastas)")
    parser.add_argument("--output_dir", default="alignment_results",
                        help="Output directory (default: alignment_results)")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="Threads for MAFFT (default: all available)")
    parser.add_argument("--mafft", default="",
                        help="Path to MAFFT executable (auto-detect)")
    parser.add_argument("--trimal", default="",
                        help="Path to trimAl executable (auto-detect)")
    parser.add_argument("--trimal_mode", default="auto",
                        choices=["auto", "gappyout", "strict", "strictplus"],
                        help="trimAl trimming strategy (default: auto)")
    return parser.parse_args()


def main() -> None:
    start_time = datetime.now()
    args = parse_args()

    out_dir = Path(args.output_dir)
    aligned_dir = out_dir / "aligned"
    trimmed_dir = out_dir / "trimmed"
    fasconcat_ready_dir = out_dir / "fasconcat_ready"
    for d in [out_dir, aligned_dir, trimmed_dir, fasconcat_ready_dir]:
        d.mkdir(parents=True, exist_ok=True)

    log_path = str(out_dir / "align_trim_step.log")
    cmd_log_path = out_dir / "commands.log"
    cmd_log_path.touch()

    mafft_exe = args.mafft or find_executable(["mafft"])
    trimal_exe = args.trimal or find_executable(["trimal", "trimAl"])

    mafft_ver = get_version(mafft_exe, "--version") if mafft_exe else "not found"
    trimal_ver = get_version(trimal_exe, "--version") if trimal_exe else "not found"

    logger = setup_logger(log_path)
    log_header(logger, args, start_time, mafft_ver, trimal_ver)

    if not mafft_exe:
        logger.error("MAFFT not found. Install via conda: conda install -c bioconda mafft")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)
    if not trimal_exe:
        logger.error("trimAl not found. Install via conda: conda install -c bioconda trimal")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    gene_fasta_dir = Path(args.gene_fasta_dir)
    faa_files = sorted(gene_fasta_dir.glob("*.faa"))
    if not faa_files:
        logger.error(f"No .faa files found in {gene_fasta_dir}")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    n_genes = len(faa_files)
    logger.info(f"Found {n_genes} orthologs to align and trim")

    n_aligned = 0
    n_trimmed = 0
    n_failed = 0
    stats_rows = []

    if TQDM_AVAILABLE:
        pbar = tqdm(total=n_genes, desc="Align+trim", unit="gene",
                    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]")

    for i, faa_path in enumerate(faa_files, 1):
        og_id = faa_path.stem

        if TQDM_AVAILABLE:
            pbar.set_postfix_str(og_id, refresh=True)
        else:
            print(f"  [{i:>{len(str(n_genes))}}/{n_genes}] {og_id}", flush=True)

        aln_path = aligned_dir / f"{og_id}.aln.faa"
        trimmed_path = trimmed_dir / f"{og_id}.trimmed.faa"

        # --- MAFFT ---
        mafft_cmd = [
            mafft_exe, "--auto",
            "--thread", str(args.threads),
            "--quiet",
            str(faa_path),
        ]
        cmd_str = " ".join(mafft_cmd)
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        try:
            with open(aln_path, "w") as aln_out:
                result = subprocess.run(
                    mafft_cmd, stdout=aln_out, stderr=subprocess.PIPE, text=True
                )
            rc = result.returncode
            with open(cmd_log_path, "a") as cl:
                cl.write(f"[{ts}] CMD: {cmd_str} > {aln_path}\n")
                cl.write(f"[{ts}] RC:  {rc}\n\n")
        except Exception as e:
            logger.error(f"MAFFT error for {og_id}: {e}")
            n_failed += 1
            if TQDM_AVAILABLE:
                pbar.update(1)
            continue

        if rc != 0:
            logger.error(f"MAFFT failed for {og_id} (rc={rc}): {result.stderr.strip()[:200]}")
            n_failed += 1
            stats_rows.append({
                "orthogroup": og_id,
                "raw_length_aa": 0,
                "trimmed_length_aa": 0,
                "pct_retained": 0.0,
                "status": "failed_alignment",
            })
            if TQDM_AVAILABLE:
                pbar.update(1)
            continue

        raw_len = alignment_length(aln_path)
        n_aligned += 1

        # --- trimAl ---
        trimal_cmd = [
            trimal_exe,
            "-in", str(aln_path),
            "-out", str(trimmed_path),
            trimal_flag(args.trimal_mode),
        ]
        rc = run_command(trimal_cmd, logger, cmd_log_path)
        if rc != 0:
            logger.error(f"trimAl failed for {og_id}")
            n_failed += 1
            stats_rows.append({
                "orthogroup": og_id,
                "raw_length_aa": raw_len,
                "trimmed_length_aa": 0,
                "pct_retained": 0.0,
                "status": "failed_trimming",
            })
            if TQDM_AVAILABLE:
                pbar.update(1)
            continue

        trimmed_len = alignment_length(trimmed_path)
        pct_retained = 100.0 * trimmed_len / raw_len if raw_len > 0 else 0.0
        n_trimmed += 1

        if pct_retained < 50.0:
            msg = (
                f"{og_id}: trimming removed >50% of columns "
                f"({raw_len} → {trimmed_len} aa, {pct_retained:.1f}% retained) — "
                "consider checking this alignment"
            )
            if TQDM_AVAILABLE:
                tqdm.write(f"  WARNING: {msg}")
            logger.warning(msg)

        stats_rows.append({
            "orthogroup": og_id,
            "raw_length_aa": raw_len,
            "trimmed_length_aa": trimmed_len,
            "pct_retained": round(pct_retained, 1),
            "status": "success",
        })

        if TQDM_AVAILABLE:
            pbar.update(1)

    if TQDM_AVAILABLE:
        pbar.close()

    logger.info(
        f"Alignment complete: {n_aligned} aligned, {n_trimmed} trimmed, {n_failed} failed"
    )

    import pandas as pd
    pd.DataFrame(stats_rows).to_csv(out_dir / "alignment_stats.tsv", sep="\t", index=False)
    logger.info(f"Alignment stats written to {out_dir / 'alignment_stats.tsv'}")

    # --- Copy trimmed files to fasconcat_ready/ with .fas extension ---
    # FASconCAT-G only reads .fas files. Copy here so user can run it manually.
    trimmed_files = sorted(trimmed_dir.glob("*.trimmed.faa"))
    n_copied = 0
    for tf in trimmed_files:
        og_stem = tf.stem.replace(".trimmed", "")
        dst = fasconcat_ready_dir / f"{og_stem}.fas"
        shutil.copy2(tf, dst)
        n_copied += 1

    logger.info(
        f"Copied {n_copied} trimmed alignments to {fasconcat_ready_dir} "
        f"(.fas extension, ready for FASconCAT-G if needed)"
    )
    logger.info(
        "To run FASconCAT-G manually:\n"
        f"  cd {fasconcat_ready_dir.resolve()}\n"
        "  perl /path/to/FASconCAT-G_v1.05.pl -s -p"
    )

    logger.info("")
    logger.info("Step 4 complete. Next step:")
    logger.info(f"  IQ-TREE3 input directory (trimmed alignments): {trimmed_dir.resolve()}")
    logger.info(
        f"  Run Step 5:  python 05_run_iqtree.py "
        f"--trimmed_dir {trimmed_dir}"
    )

    log_footer(logger, start_time, "SUCCESS")


if __name__ == "__main__":
    main()
