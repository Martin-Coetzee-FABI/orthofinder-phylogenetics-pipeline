#!/usr/bin/env python3
"""
pipeline_runner_gui.py — Detached pipeline runner for the Streamlit GUI
=======================================================================
This script is launched as a DETACHED subprocess by app.py.
It runs the six pipeline steps in sequence and writes all output to
pipeline_run.log.  It should NOT be run directly by the user.

Usage (internal, called by app.py):
    python pipeline_runner_gui.py <config.json>
"""

import glob
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

PIPELINE_LOG = Path("pipeline_run.log")
PID_FILE     = Path(".pipeline_gui.pid")


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log_write(log_f, msg: str):
    line = f"[{ts()}] {msg}"
    log_f.write(line + "\n")
    log_f.flush()


def run_step(label: str, args: list[str], log_f) -> int:
    """Run one pipeline step, streaming its output to the log file.
    Returns the subprocess exit code."""
    sep = "=" * 56
    for line in ["", sep, f"  {label}", f"  Started: {ts()}", sep, ""]:
        log_f.write(line + "\n")
    log_f.flush()

    proc = subprocess.Popen(
        [sys.executable] + args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    # Update PID file to the currently-running step subprocess so that
    # the GUI's Stop button can kill it.
    PID_FILE.write_text(str(proc.pid))

    for line in proc.stdout:
        log_f.write(line)
        log_f.flush()

    proc.wait()

    log_f.write(f"\n  Finished: {ts()}  |  exit code: {proc.returncode}\n")
    log_f.flush()

    return proc.returncode


def main():
    if len(sys.argv) < 2:
        print("Usage: pipeline_runner_gui.py <config.json>", file=sys.stderr)
        sys.exit(1)

    cfg_path = Path(sys.argv[1])
    if not cfg_path.exists():
        print(f"Config file not found: {cfg_path}", file=sys.stderr)
        sys.exit(1)

    cfg = json.loads(cfg_path.read_text())

    # Step range from GUI settings (default: all steps)
    step_start = int(cfg.get("step_start", 1))
    step_end   = int(cfg.get("step_end",   5))

    # Write our own PID so the GUI can track / kill us
    PID_FILE.write_text(str(os.getpid()))

    with open(PIPELINE_LOG, "a") as log_f:

        log_f.write(f"\n{'#'*58}\n")
        log_f.write(f"#  OrthoFinder3 Pipeline — GUI run\n")
        log_f.write(f"#  Started : {ts()}\n")
        log_f.write(f"#  Config  : {cfg_path.resolve()}\n")
        log_f.write(f"#  Dir     : {Path('.').resolve()}\n")
        log_f.write(f"#  Steps   : {step_start}–{step_end}  (PDF report always runs)\n")
        log_f.write(f"{'#'*58}\n\n")
        log_f.write(f"  Settings:\n")
        for k, v in cfg.items():
            log_f.write(f"    {k:<22}: {v}\n")
        log_f.write("\n")
        log_f.flush()

        # Build re-usable arg fragments
        species_tree_args = (
            ["--species_tree", cfg["species_tree"]]
            if cfg.get("species_tree") else []
        )
        timeout_args = ["--timeout", str(cfg.get("timeout", 86400))]

        # ── Step 1 — QC ──────────────────────────────────────────────────
        if step_start > 1:
            log_write(log_f, f"[SKIP] Step 1 — QC proteomes "
                             f"(outside selected range: steps {step_start}–{step_end})")
        elif Path("qc_results/.step_complete").exists():
            log_write(log_f, "[SKIP] Step 1 — QC proteomes (already complete)")
        else:
            rc = run_step("Step 1 — QC proteomes", [
                "01_qc_proteomes.py",
                "--input_dir",    cfg["input_dir"],
                "--min_proteins", str(cfg["min_proteins"]),
                "--output_dir",   "qc_results",
            ], log_f)
            if rc != 0:
                log_write(log_f, f"[FATAL] Step 1 failed (rc={rc}). Pipeline aborted.")
                PID_FILE.unlink(missing_ok=True)
                sys.exit(rc)
            Path("qc_results/.step_complete").touch()
            log_write(log_f, "Step 1 complete ✓")

        # ── Step 2 — OrthoFinder ─────────────────────────────────────────
        if step_start > 2 or step_end < 2:
            log_write(log_f, f"[SKIP] Step 2 — OrthoFinder "
                             f"(outside selected range: steps {step_start}–{step_end})")
        else:
            rc = run_step("Step 2 — OrthoFinder3", [
                "02_run_orthofinder.py",
                "--input_dir",        "qc_results/passed_proteomes",
                "--output_dir",       "orthofinder_results",
                "--threads",          str(cfg["threads"]),
                "--search_prog",      cfg["search_prog"],
                "--gene_tree_method", cfg["gene_tree_method"],
                "--msa_prog",         cfg["msa_prog"],
            ] + species_tree_args + timeout_args, log_f)
            if rc != 0:
                log_write(log_f, f"[FATAL] Step 2 failed (rc={rc}). Pipeline aborted.")
                PID_FILE.unlink(missing_ok=True)
                sys.exit(rc)
            log_write(log_f, "Step 2 complete ✓")

        # Locate the OrthoFinder Results_* directory (needed by Step 3)
        of_dirs = sorted(glob.glob("orthofinder_results/Results_*"))
        of_dir  = of_dirs[-1] if of_dirs else ""
        if of_dir:
            log_write(log_f, f"OrthoFinder results: {of_dir}")

        # ── Step 3 — Extract orthologs ────────────────────────────────────
        if step_start > 3 or step_end < 3:
            log_write(log_f, f"[SKIP] Step 3 — Extract orthologs "
                             f"(outside selected range: steps {step_start}–{step_end})")
        elif Path("ortholog_results/.step_complete").exists():
            log_write(log_f, "[SKIP] Step 3 — Extract orthologs (already complete)")
        else:
            if not of_dir:
                log_write(log_f, "[FATAL] Step 3 requires OrthoFinder results. "
                                 "Run Step 2 first.")
                PID_FILE.unlink(missing_ok=True)
                sys.exit(1)
            rc = run_step("Step 3 — Extract single-copy orthologs", [
                "03_extract_orthologs.py",
                "--orthofinder_dir", of_dir,
                "--faa_dir",         "qc_results/passed_proteomes",
                "--min_presence",    str(cfg["min_presence"]),
                "--output_dir",      "ortholog_results",
            ], log_f)
            if rc != 0:
                log_write(log_f, f"[FATAL] Step 3 failed (rc={rc}). Pipeline aborted.")
                PID_FILE.unlink(missing_ok=True)
                sys.exit(rc)
            Path("ortholog_results/.step_complete").touch()
            log_write(log_f, "Step 3 complete ✓")

        # ── Step 4 — Align & trim ─────────────────────────────────────────
        if step_start > 4 or step_end < 4:
            log_write(log_f, f"[SKIP] Step 4 — Align & trim "
                             f"(outside selected range: steps {step_start}–{step_end})")
        elif Path("alignment_results/.step_complete").exists():
            log_write(log_f, "[SKIP] Step 4 — Align & trim (already complete)")
        else:
            if not Path("ortholog_results/gene_fastas").is_dir():
                log_write(log_f, "[FATAL] Step 4 requires ortholog_results/gene_fastas/. "
                                 "Run Step 3 first.")
                PID_FILE.unlink(missing_ok=True)
                sys.exit(1)
            rc = run_step("Step 4 — Align (MAFFT) & trim (trimAl)", [
                "04_align_trim_concat.py",
                "--gene_fasta_dir", "ortholog_results/gene_fastas",
                "--output_dir",     "alignment_results",
                "--threads",        str(cfg["threads"]),
                "--trimal_mode",    cfg["trimal_mode"],
            ], log_f)
            if rc != 0:
                log_write(log_f, f"[FATAL] Step 4 failed (rc={rc}). Pipeline aborted.")
                PID_FILE.unlink(missing_ok=True)
                sys.exit(rc)
            Path("alignment_results/.step_complete").touch()
            log_write(log_f, "Step 4 complete ✓")

        # ── Step 5 — IQ-TREE ─────────────────────────────────────────────
        if step_start > 5 or step_end < 5:
            log_write(log_f, f"[SKIP] Step 5 — IQ-TREE3 "
                             f"(outside selected range: steps {step_start}–{step_end})")
        elif Path("iqtree_results/.step_complete").exists():
            log_write(log_f, "[SKIP] Step 5 — IQ-TREE3 (already complete)")
        else:
            if not Path("alignment_results/trimmed").is_dir():
                log_write(log_f, "[FATAL] Step 5 requires alignment_results/trimmed/. "
                                 "Run Step 4 first.")
                PID_FILE.unlink(missing_ok=True)
                sys.exit(1)
            rc = run_step("Step 5 — IQ-TREE3 species tree", [
                "05_run_iqtree.py",
                "--trimmed_dir", "alignment_results/trimmed",
                "--output_dir",  "iqtree_results",
                "--mode",        cfg["iqtree_mode"],
                "--model",       cfg["iqtree_model"],
                "--bootstrap",   str(cfg["bootstrap"]),
                "--alrt",        str(cfg["alrt"]),
                "--threads",     str(cfg["threads"]),
            ], log_f)
            if rc != 0:
                log_write(log_f, f"[FATAL] Step 5 failed (rc={rc}). Pipeline aborted.")
                PID_FILE.unlink(missing_ok=True)
                sys.exit(rc)
            Path("iqtree_results/.step_complete").touch()
            log_write(log_f, "Step 5 complete ✓")

        # ── Write run_settings.json for the PDF report ────────────────────
        try:
            run_settings = {
                "input_dir":       cfg["input_dir"],
                "min_proteins":    str(cfg["min_proteins"]),
                "min_presence":    str(cfg["min_presence"]),
                "threads":         str(cfg["threads"]),
                "of_search":       cfg["search_prog"],
                "of_gene_tree":    cfg["gene_tree_method"],
                "of_msa":          cfg["msa_prog"],
                "of_species_tree": cfg.get("species_tree", ""),
                "of_timeout":      str(cfg.get("timeout", 86400)),
                "trimal_mode":     cfg["trimal_mode"],
                "iqtree_mode":     cfg["iqtree_mode"],
                "iqtree_model":    cfg["iqtree_model"],
                "bootstrap":       str(cfg["bootstrap"]),
                "alrt":            str(cfg["alrt"]),
                "step_start":      str(step_start),
                "step_end":        str(step_end),
            }
            with open("run_settings.json", "w") as sf:
                json.dump(run_settings, sf, indent=2)
            log_write(log_f, "run_settings.json written for PDF report")
        except Exception as e:
            log_write(log_f, f"[WARNING] Could not write run_settings.json: {e}")

        # ── Step 6 — PDF report ───────────────────────────────────────────
        rc = run_step("Step 6 — Generate PDF report", [
            "06_generate_report.py",
            "--output_dir",    "report",
            "--settings_file", "run_settings.json",
        ], log_f)
        if rc != 0:
            # Non-fatal: the tree is still available
            log_write(log_f, f"[WARNING] Step 6 failed (rc={rc}). "
                             "PDF report not generated, but species tree is available.")
        else:
            log_write(log_f, "Step 6 complete ✓")

        # ── Done ──────────────────────────────────────────────────────────
        log_f.write(f"\n{'#'*58}\n")
        log_f.write(f"#  Pipeline complete!\n")
        log_f.write(f"#  Finished : {ts()}\n")
        log_f.write(f"#\n")
        log_f.write(f"#  Species tree : iqtree_results/species_tree.treefile\n")
        log_f.write(f"#  PDF report   : report/pipeline_summary.pdf\n")
        log_f.write(f"{'#'*58}\n\n")
        log_f.flush()

    PID_FILE.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
