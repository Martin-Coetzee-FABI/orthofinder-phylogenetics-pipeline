#!/usr/bin/env python3
"""
01_qc_proteomes.py — Proteome quality control before OrthoFinder

Checks all .faa files in a directory and flags low-quality genomes.
Outputs: QC summary TSV, human-readable report, symlinks to passing genomes,
         list of failed genomes, and a timestamped log.
"""

import argparse
import logging
import os
import shutil
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

from Bio import SeqIO


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
    logger.info(f"Script:      01_qc_proteomes.py")
    logger.info(f"Started:     {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Parameters:  --input_dir {args.input_dir}")
    logger.info(f"             --min_proteins {args.min_proteins}")
    logger.info(f"             --min_median_len {args.min_median_len}")
    logger.info(f"             --max_short_pct {args.max_short_pct}")
    logger.info(f"             --output_dir {args.output_dir}")
    if args.exclude:
        logger.info(f"             --exclude {args.exclude}")
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


def qc_faa_file(faa_path: Path) -> dict:
    """Compute QC metrics for a single .faa file."""
    lengths = []
    ids_seen = []

    try:
        for record in SeqIO.parse(str(faa_path), "fasta"):
            ids_seen.append(record.id)
            lengths.append(len(record.seq))
    except Exception as e:
        return {"error": str(e)}

    if not lengths:
        return {"error": "empty file or no sequences parsed"}

    ids_seen_counter = Counter(ids_seen)
    n_dup = sum(v - 1 for v in ids_seen_counter.values() if v > 1)
    pct_dup = 100.0 * n_dup / len(ids_seen) if ids_seen else 0.0

    sorted_lens = sorted(lengths)
    n = len(sorted_lens)
    if n % 2 == 1:
        median_len = sorted_lens[n // 2]
    else:
        median_len = (sorted_lens[n // 2 - 1] + sorted_lens[n // 2]) / 2.0

    n_short = sum(1 for l in lengths if l < 50)
    pct_short = 100.0 * n_short / n if n else 0.0

    return {
        "n_proteins": n,
        "median_len": median_len,
        "pct_short": pct_short,
        "pct_dup": pct_dup,
        "n_dup": n_dup,
        "error": None,
    }


def determine_verdict(metrics: dict, args: argparse.Namespace, manually_excluded: set) -> tuple[str, str]:
    """Return (verdict, reason) for a genome."""
    genome_name = metrics.get("name", "")

    if genome_name in manually_excluded:
        return "FAIL", "manually excluded via --exclude"

    if metrics.get("error"):
        return "FAIL", f"file error: {metrics['error']}"

    flags = []
    if metrics["n_proteins"] < args.min_proteins:
        flags.append(f"protein count {metrics['n_proteins']} < {args.min_proteins}")
    if metrics["median_len"] < args.min_median_len:
        flags.append(f"median length {metrics['median_len']:.1f} aa < {args.min_median_len}")
    if metrics["pct_short"] > args.max_short_pct:
        flags.append(f"{metrics['pct_short']:.1f}% proteins <50 aa > {args.max_short_pct}%")
    if metrics["pct_dup"] > 0:
        flags.append(f"duplicate IDs: {metrics['n_dup']} duplicates ({metrics['pct_dup']:.2f}%)")

    if metrics["pct_dup"] > 0 or len(flags) >= 2:
        return "FAIL", "; ".join(flags)
    elif len(flags) == 1:
        return "WARN", flags[0]
    else:
        return "PASS", ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="QC check all .faa proteome files before OrthoFinder."
    )
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing .faa files")
    parser.add_argument("--min_proteins", type=int, default=500,
                        help="Minimum protein count threshold (default: 500)")
    parser.add_argument("--min_median_len", type=float, default=100.0,
                        help="Minimum median protein length in aa (default: 100)")
    parser.add_argument("--max_short_pct", type=float, default=20.0,
                        help="Maximum %% proteins <50 aa (default: 20.0)")
    parser.add_argument("--exclude", default="",
                        help="Comma-separated genome names to exclude regardless of QC")
    parser.add_argument("--output_dir", default="qc_results",
                        help="Output directory (default: qc_results)")
    parser.add_argument("--copy", action="store_true",
                        help="Copy files to passed_proteomes/ instead of symlinking")
    return parser.parse_args()


def main() -> None:
    start_time = datetime.now()
    args = parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    passed_dir = out_dir / "passed_proteomes"
    passed_dir.mkdir(parents=True, exist_ok=True)

    log_path = str(out_dir / "qc_step.log")
    logger = setup_logger(log_path)
    log_header(logger, args, start_time)

    manually_excluded = {s.strip() for s in args.exclude.split(",") if s.strip()}
    if manually_excluded:
        logger.info(f"Manually excluded genomes: {', '.join(sorted(manually_excluded))}")

    input_dir = Path(args.input_dir)
    faa_files = sorted(input_dir.glob("*.faa"))

    if not faa_files:
        logger.error(f"No .faa files found in {input_dir}")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    logger.info(f"Found {len(faa_files)} .faa files in {input_dir}")
    logger.info(f"QC thresholds: min_proteins={args.min_proteins}, "
                f"min_median_len={args.min_median_len}, "
                f"max_short_pct={args.max_short_pct}%")

    n_total = len(faa_files)
    results = []
    for i, faa_path in enumerate(faa_files, 1):
        genome_name = faa_path.stem
        print(f"  [{i:>{len(str(n_total))}}/{n_total}] Checking: {genome_name}", flush=True)
        metrics = qc_faa_file(faa_path)
        metrics["name"] = genome_name
        metrics["path"] = str(faa_path.resolve())
        verdict, reason = determine_verdict(metrics, args, manually_excluded)
        metrics["verdict"] = verdict
        metrics["reason"] = reason
        results.append(metrics)
        logger.debug(f"  {genome_name}: {verdict} — {reason if reason else 'all metrics OK'}")

    results.sort(key=lambda r: r.get("n_proteins", 0) if not r.get("error") else -1)

    tsv_path = out_dir / "proteome_qc_summary.tsv"
    with open(tsv_path, "w") as f:
        header = "\t".join([
            "genome", "n_proteins", "median_len_aa",
            "pct_short", "pct_dup", "n_dup", "verdict", "reason"
        ])
        f.write(header + "\n")
        for r in results:
            if r.get("error") and r["verdict"] == "FAIL":
                row = "\t".join([
                    r["name"], "NA", "NA", "NA", "NA", "NA",
                    r["verdict"], r.get("reason", r["error"])
                ])
            else:
                row = "\t".join([
                    r["name"],
                    str(r.get("n_proteins", "NA")),
                    f"{r.get('median_len', 0):.1f}",
                    f"{r.get('pct_short', 0):.2f}",
                    f"{r.get('pct_dup', 0):.2f}",
                    str(r.get("n_dup", "NA")),
                    r["verdict"],
                    r.get("reason", ""),
                ])
            f.write(row + "\n")
    logger.info(f"QC summary written to {tsv_path}")

    report_lines = []
    report_lines.append("=" * 72)
    report_lines.append("  PROTEOME QC REPORT")
    report_lines.append(f"  Run at: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"  Input:  {input_dir}")
    report_lines.append("=" * 72)
    report_lines.append(
        f"{'GENOME':<35} {'PROTEINS':>10} {'MED_LEN':>9} {'%SHORT':>7} {'%DUP':>6}  VERDICT"
    )
    report_lines.append("-" * 72)
    for r in results:
        if r.get("error") and r["verdict"] == "FAIL":
            line = f"{r['name']:<35} {'ERROR':>10} {'':>9} {'':>7} {'':>6}  {r['verdict']}"
        else:
            line = (
                f"{r['name']:<35} "
                f"{r.get('n_proteins', 0):>10} "
                f"{r.get('median_len', 0):>9.1f} "
                f"{r.get('pct_short', 0):>7.2f} "
                f"{r.get('pct_dup', 0):>6.2f}  "
                f"{r['verdict']}"
            )
            if r.get("reason"):
                line += f"  [{r['reason']}]"
        report_lines.append(line)

    n_pass = sum(1 for r in results if r["verdict"] == "PASS")
    n_warn = sum(1 for r in results if r["verdict"] == "WARN")
    n_fail = sum(1 for r in results if r["verdict"] == "FAIL")
    report_lines.append("-" * 72)
    report_lines.append(f"  PASS: {n_pass}   WARN: {n_warn}   FAIL: {n_fail}   TOTAL: {len(results)}")
    report_lines.append("=" * 72)

    report_text = "\n".join(report_lines)
    print(report_text)

    report_path = out_dir / "proteome_qc_report.txt"
    with open(report_path, "w") as f:
        f.write(report_text + "\n")
    logger.info(f"Human-readable report written to {report_path}")

    passed_genomes = [r for r in results if r["verdict"] in ("PASS", "WARN")]
    failed_genomes = [r for r in results if r["verdict"] == "FAIL"]

    failed_path = out_dir / "failed_genomes.txt"
    with open(failed_path, "w") as f:
        for r in failed_genomes:
            f.write(f"{r['name']}\t{r.get('reason', '')}\n")
    logger.info(f"Failed genomes written to {failed_path}")

    for r in passed_genomes:
        src = Path(r["path"])
        dst = passed_dir / src.name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        if args.copy:
            shutil.copy2(src, dst)
            logger.debug(f"Copied {src.name} to passed_proteomes/")
        else:
            try:
                dst.symlink_to(src)
                logger.debug(f"Symlinked {src.name} to passed_proteomes/")
            except OSError:
                shutil.copy2(src, dst)
                logger.warning(f"Symlink failed for {src.name}, copied instead")

    logger.info(f"Passing genomes ({len(passed_genomes)}) written to {passed_dir}")

    for r in results:
        verdict = r["verdict"]
        name = r["name"]
        reason = r.get("reason", "")
        if verdict == "PASS":
            logger.info(f"  PASS  {name}")
        elif verdict == "WARN":
            logger.warning(f"  WARN  {name} — {reason}")
        else:
            logger.warning(f"  FAIL  {name} — {reason}")

    logger.info(
        f"QC complete: {n_pass} PASS, {n_warn} WARN (included), {n_fail} FAIL (excluded) "
        f"out of {len(results)} genomes"
    )

    if len(passed_genomes) == 0:
        logger.error("Zero genomes passed QC — cannot proceed with OrthoFinder.")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    log_footer(logger, start_time, "SUCCESS")


if __name__ == "__main__":
    main()
