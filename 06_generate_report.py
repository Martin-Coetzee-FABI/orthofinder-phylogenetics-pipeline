#!/usr/bin/env python3
"""
06_generate_report.py — Generate a PDF summary of all pipeline results.

Reads outputs from Steps 1–5 and produces a multi-page PDF in report/.
Any step whose output is missing is noted but does not prevent the report.

Pages produced:
  1. Title & at-a-glance overview
  2. Pipeline: Steps & Settings Used
  3. QC results — table of all genomes with verdicts
  4. Ortholog selection — counts, completeness, presence matrix heat-map
  5. Functional overview — broad categories of selected orthologs
  6. Alignment statistics — length before/after trimming, per-gene distribution
  7. Species tree — rendered phylogeny + raw Newick string
  8. References
"""

import argparse
import json
import logging
import re
import sys
import textwrap
from collections import Counter
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                           # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np

try:
    from Bio import Phylo, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

# ── Colour palette ────────────────────────────────────────────────────────────
C_HEADER  = "#2C3E50"   # dark blue-grey  (page header band)
C_ACCENT  = "#2980B9"   # blue            (chart bars, borders)
C_PASS    = "#27AE60"   # green
C_WARN    = "#F39C12"   # amber
C_FAIL    = "#E74C3C"   # red
C_LIGHT   = "#ECF0F1"   # near-white      (odd table rows)
C_MID     = "#BDC3C7"   # grey            (table borders)
C_TEXT    = "#2C3E50"   # near-black      (body text)

# ── Functional gene categories (keyword → category) ───────────────────────────
FUNC_CATEGORIES = [
    ("Ribosome & Translation",
     ["ribosom", "trna", "aminoacyl", "translation elongation",
      "translation initiation", "peptidyl transferase", "peptide chain release",
      "ribosome recycling"]),
    ("DNA Replication & Repair",
     ["dna replic", "dna repair", "dna helicase", "topoisomerase", "gyrase",
      "primase", "single-strand", "mismatch repair", "nucleotide excision",
      "dnaa", "dnab", "dnac", "dnag"]),
    ("Transcription",
     ["rna polymerase", "transcription factor", "sigma factor", "anti-sigma",
      "rna-binding", "transcriptional regulator"]),
    ("Cell Division & Wall",
     ["cell division", "fts", "mur", "murein", "peptidoglycan", "divisome",
      "septum"]),
    ("Energy Metabolism",
     ["atp synthase", "atpase", "nadh", "oxidoreductase", "cytochrome",
      "electron transport", "proton pump"]),
    ("Metabolic Enzymes",
     ["transferase", "synthase", "synthetase", "kinase", "dehydrogenase",
      "reductase", "isomerase", "lyase", "hydrolase", "phosphatase",
      "carboxylase"]),
    ("Transport & Secretion",
     ["transport", "permease", "abc transporter", "porin", "efflux",
      "secretion", "type ii", "type iii", "type iv"]),
    ("Chaperones & Stress",
     ["chaperone", "heat shock", "groel", "groe", "dnak", "dnaj",
      "stress", "protease clp", "lon protease"]),
    ("Signal Transduction",
     ["two-component", "histidine kinase", "response regulator",
      "chemotaxis", "signal transduction"]),
    ("Hypothetical / Unknown",
     ["hypothetical", "uncharacterized", "unknown function", "duf",
      "putative"]),
]


# ── Logging ───────────────────────────────────────────────────────────────────

def setup_logger(log_path: str) -> logging.Logger:
    logger = logging.getLogger("report")
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


def log_header(logger, args, start_time):
    logger.info("=" * 54)
    logger.info("Script:      06_generate_report.py")
    logger.info(f"Started:     {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Parameters:  --output_dir {args.output_dir}")
    logger.info("=" * 54)


def log_footer(logger, start_time, status):
    end = datetime.now()
    d = end - start_time
    m, s = divmod(int(d.total_seconds()), 60)
    logger.info("=" * 54)
    logger.info(f"Finished:    {end.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Duration:    {m}m {s}s")
    logger.info(f"Exit status: {status}")
    logger.info("=" * 54)


# ── Data loading ──────────────────────────────────────────────────────────────

def load_qc_summary(qc_dir: Path, logger) -> pd.DataFrame | None:
    path = qc_dir / "proteome_qc_summary.tsv"
    if not path.exists():
        logger.warning(f"QC summary not found: {path}")
        return None
    df = pd.read_csv(path, sep="\t")
    return df


def load_ortholog_data(ortholog_dir: Path, logger) -> tuple[dict, pd.DataFrame | None, list]:
    summary_path = ortholog_dir / "tables" / "ortholog_summary.tsv"
    included_path = ortholog_dir / "tables" / "orthologs_included.tsv"
    og_list_path  = ortholog_dir / "included_orthologs.txt"

    summary = {}
    if summary_path.exists():
        df = pd.read_csv(summary_path, sep="\t")
        summary = df.iloc[0].to_dict() if not df.empty else {}
    else:
        logger.warning(f"Ortholog summary not found: {summary_path}")

    included_df = None
    if included_path.exists():
        included_df = pd.read_csv(included_path, sep="\t")
    else:
        logger.warning(f"Orthologs included table not found: {included_path}")

    included_ogs = []
    if og_list_path.exists():
        included_ogs = og_list_path.read_text().strip().splitlines()
    else:
        logger.warning(f"included_orthologs.txt not found: {og_list_path}")

    return summary, included_df, included_ogs


def load_alignment_stats(alignment_dir: Path, logger) -> pd.DataFrame | None:
    path = alignment_dir / "alignment_stats.tsv"
    if not path.exists():
        logger.warning(f"Alignment stats not found: {path}")
        return None
    return pd.read_csv(path, sep="\t")


def find_orthofinder_dir(orthofinder_results: Path, logger) -> Path | None:
    hits = sorted(orthofinder_results.glob("Results_*"))
    if hits:
        return hits[-1]
    logger.warning(f"No OrthoFinder Results_* directory found under {orthofinder_results}")
    return None


def extract_description(header_str: str) -> str:
    """Pull a plain product name from a Prokka or PGAP FASTA header."""
    parts = header_str.split(None, 1)
    if len(parts) < 2:
        return "unknown"
    desc = parts[1]
    # PGAP: [protein=...] or [gene=...]
    m = re.search(r'\[protein=([^\]]+)\]', desc)
    if m:
        return m.group(1).strip()
    # Strip any bracket annotations common in NCBI headers
    desc = re.sub(r'\[.*?\]', '', desc).strip()
    return desc if desc else "unknown"


def load_gene_annotations(
    faa_dir: Path,
    orthofinder_dir: Path | None,
    included_ogs: list[str],
    logger,
) -> dict[str, str]:
    """Return {og_id: product_description} for each included ortholog."""
    if not BIO_AVAILABLE:
        logger.warning("BioPython not available — skipping gene annotation loading")
        return {}
    if not faa_dir.exists():
        logger.warning(f"FAA directory not found: {faa_dir}")
        return {}

    # Build locus_tag → description
    logger.info("Loading gene product descriptions from .faa headers ...")
    locus_to_desc: dict[str, str] = {}
    for faa_path in sorted(faa_dir.glob("*.faa")):
        for rec in SeqIO.parse(str(faa_path), "fasta"):
            locus_to_desc[rec.id] = extract_description(rec.description)

    if not orthofinder_dir:
        return {}

    og_txt = orthofinder_dir / "Orthogroups" / "Orthogroups.txt"
    if not og_txt.exists():
        logger.warning(f"Orthogroups.txt not found: {og_txt}")
        return {}

    og_to_tags: dict[str, list[str]] = {}
    with open(og_txt) as f:
        for line in f:
            parts = line.strip().split(":", 1)
            if len(parts) == 2:
                og_to_tags[parts[0].strip()] = parts[1].split()

    og_to_desc: dict[str, str] = {}
    for og_id in included_ogs:
        tags = og_to_tags.get(og_id, [])
        chosen = "unknown"
        for tag in tags:
            d = locus_to_desc.get(tag, "")
            if d and "hypothetical" not in d.lower():
                chosen = d
                break
        if chosen == "unknown" and tags:
            chosen = locus_to_desc.get(tags[0], "unknown")
        og_to_desc[og_id] = chosen

    return og_to_desc


def categorize_genes(og_to_desc: dict[str, str]) -> Counter:
    """Assign each OG to a functional category; return counts."""
    counts: Counter = Counter()
    for og_id, desc in og_to_desc.items():
        d = desc.lower()
        matched = False
        for cat_name, keywords in FUNC_CATEGORIES:
            if any(kw in d for kw in keywords):
                counts[cat_name] += 1
                matched = True
                break
        if not matched:
            counts["Other / Miscellaneous"] += 1
    return counts


def load_run_settings(settings_file: Path, logger) -> dict:
    """Load pipeline run settings from JSON file or fall back to parsing command files."""
    defaults = {
        "input_dir": "N/A",
        "min_proteins": "N/A",
        "min_presence": "N/A",
        "threads": "N/A",
        "of_search": "N/A",
        "of_gene_tree": "N/A",
        "of_msa": "N/A",
        "of_species_tree": None,
        "trimal_mode": "N/A",
        "iqtree_mode": "N/A",
        "iqtree_model": "N/A",
        "bootstrap": "N/A",
        "alrt": "N/A",
        "step_start": 1,
        "step_end": 5,
    }

    if settings_file.exists():
        try:
            with open(settings_file) as f:
                data = json.load(f)
            merged = {**defaults, **data}
            logger.info(f"Loaded run settings from {settings_file}")
            return merged
        except Exception as e:
            logger.warning(f"Could not parse {settings_file}: {e}")

    # Fall back: try to parse command text files
    settings = dict(defaults)

    iqtree_cmd_path = Path("iqtree_results/iqtree_command.txt")
    if iqtree_cmd_path.exists():
        try:
            text = iqtree_cmd_path.read_text()
            for line in text.splitlines():
                if line.strip().startswith("Command:"):
                    cmd = line.split("Command:", 1)[1].strip()
                    # Parse -m model
                    m = re.search(r'-m\s+(\S+)', cmd)
                    if m:
                        settings["iqtree_model"] = m.group(1)
                    # Parse -B bootstrap
                    m = re.search(r'-B\s+(\d+)', cmd)
                    if m:
                        settings["bootstrap"] = m.group(1)
                    # Parse -alrt
                    m = re.search(r'-alrt\s+(\d+)', cmd)
                    if m:
                        settings["alrt"] = m.group(1)
                    # Parse -T threads
                    m = re.search(r'-T\s+(\S+)', cmd)
                    if m:
                        settings["threads"] = m.group(1)
                    break
        except Exception as e:
            logger.warning(f"Could not parse iqtree_command.txt: {e}")

    of_cmd_path = Path("orthofinder_results/orthofinder_command.txt")
    if of_cmd_path.exists():
        try:
            text = of_cmd_path.read_text()
            for line in text.splitlines():
                if line.strip().startswith("Command:"):
                    cmd = line.split("Command:", 1)[1].strip()
                    m = re.search(r'-S\s+(\S+)', cmd)
                    if m:
                        settings["of_search"] = m.group(1)
                    m = re.search(r'-M\s+(\S+)', cmd)
                    if m:
                        settings["of_gene_tree"] = m.group(1)
                    m = re.search(r'-A\s+(\S+)', cmd)
                    if m:
                        settings["of_msa"] = m.group(1)
                    break
        except Exception as e:
            logger.warning(f"Could not parse orthofinder_command.txt: {e}")

    return settings


# ── PDF page helpers ──────────────────────────────────────────────────────────

PAGE_W, PAGE_H = 8.5, 11.0   # letter portrait, inches


def new_page(title: str, run_date: str) -> tuple[plt.Figure, plt.Axes]:
    """Create a new figure with a coloured header band; return (fig, content_ax)."""
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))

    # Header band (top 0.7 inches)
    hdr = fig.add_axes([0, 0.91, 1, 0.09])
    hdr.set_facecolor(C_HEADER)
    hdr.set_xlim(0, 1); hdr.set_ylim(0, 1)
    hdr.axis("off")
    hdr.text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
             color="white", fontsize=10, fontweight="bold", va="center")
    hdr.text(0.98, 0.55, run_date,
             color="#BDC3C7", fontsize=8, va="center", ha="right")

    # Thin accent line below header
    accent = fig.add_axes([0, 0.905, 1, 0.007])
    accent.set_facecolor(C_ACCENT); accent.axis("off")

    # Page title
    title_ax = fig.add_axes([0.03, 0.865, 0.94, 0.04])
    title_ax.axis("off")
    title_ax.text(0, 0.5, title, color=C_HEADER, fontsize=14,
                  fontweight="bold", va="center")

    # Thin separator
    sep = fig.add_axes([0.03, 0.858, 0.94, 0.003])
    sep.set_facecolor(C_MID); sep.axis("off")

    # Main content area
    content = fig.add_axes([0.04, 0.04, 0.92, 0.81])
    content.axis("off")
    return fig, content


def df_to_table(ax: plt.Axes, df: pd.DataFrame,
                col_widths: list[float] | None = None,
                row_colors: bool = True,
                fontsize: int = 8) -> None:
    """Render a DataFrame as a styled table on the given axes."""
    n_rows, n_cols = df.shape

    if col_widths is None:
        col_widths = [1.0 / n_cols] * n_cols

    cell_text = df.astype(str).values.tolist()
    col_labels = list(df.columns)

    row_colour_list = (
        [C_LIGHT if i % 2 == 0 else "white" for i in range(n_rows)]
        if row_colors else ["white"] * n_rows
    )

    tbl = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellLoc="left",
        loc="upper left",
        colWidths=col_widths,
        bbox=[0, 0, 1, 1],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)

    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(C_MID)
        if r == 0:
            cell.set_facecolor(C_HEADER)
            cell.set_text_props(color="white", fontweight="bold")
        else:
            cell.set_facecolor(row_colour_list[r - 1])
            cell.set_text_props(color=C_TEXT)


# ── Page drawing functions ─────────────────────────────────────────────────────

def page_overview(pdf: PdfPages, run_date: str,
                  qc_df: pd.DataFrame | None,
                  summary: dict,
                  aln_df: pd.DataFrame | None,
                  logger) -> None:

    fig, ax = new_page("1 · Pipeline At-a-Glance", run_date)

    stats = []

    if qc_df is not None:
        n_total = len(qc_df)
        n_pass  = (qc_df["verdict"] == "PASS").sum()
        n_warn  = (qc_df["verdict"] == "WARN").sum()
        n_fail  = (qc_df["verdict"] == "FAIL").sum()
        stats += [
            ("Genomes analysed",              str(n_total)),
            ("Genomes passed QC (PASS)",       str(n_pass)),
            ("Genomes with warnings (WARN)",   str(n_warn)),
            ("Genomes failed QC (FAIL)",       str(n_fail)),
        ]

    if summary:
        stats += [
            ("Orthogroups assessed",
             f"{int(summary.get('total_orthogroups_assessed', 0)):,}"),
            ("Orthologs passing strict filter (100% presence)",
             f"{int(summary.get('orthogroups_passing_strict', 0)):,}"),
            ("Orthologs selected for analysis",
             f"{int(summary.get('orthogroups_selected', 0)):,}"),
            ("Species in analysis",
             f"{int(summary.get('total_species', 0)):,}"),
        ]

    if aln_df is not None:
        success = aln_df[aln_df["status"] == "success"]
        total_raw     = int(success["raw_length_aa"].sum())
        total_trimmed = int(success["trimmed_length_aa"].sum())
        pct_kept      = 100.0 * total_trimmed / total_raw if total_raw > 0 else 0
        stats += [
            ("Genes aligned & trimmed",  f"{len(success):,}"),
            ("Total alignment length (raw)",     f"{total_raw:,} aa"),
            ("Total alignment length (trimmed)", f"{total_trimmed:,} aa"),
            ("Columns retained after trimming",  f"{pct_kept:.1f}%"),
        ]

    if not stats:
        ax.text(0.5, 0.5, "No pipeline output found.\nRun the pipeline first.",
                ha="center", va="center", fontsize=12, color=C_FAIL,
                transform=ax.transAxes)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig); return

    # Draw as two-column key-value table
    fig2 = plt.figure(figsize=(PAGE_W, PAGE_H))
    hdr = fig2.add_axes([0, 0.91, 1, 0.09])
    hdr.set_facecolor(C_HEADER); hdr.set_xlim(0,1); hdr.set_ylim(0,1); hdr.axis("off")
    hdr.text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
             color="white", fontsize=10, fontweight="bold", va="center")
    hdr.text(0.98, 0.55, run_date, color="#BDC3C7", fontsize=8, va="center", ha="right")
    acc = fig2.add_axes([0, 0.905, 1, 0.007]); acc.set_facecolor(C_ACCENT); acc.axis("off")
    tit = fig2.add_axes([0.03, 0.865, 0.94, 0.04]); tit.axis("off")
    tit.text(0, 0.5, "1 · Pipeline At-a-Glance", color=C_HEADER, fontsize=14, fontweight="bold", va="center")
    sep = fig2.add_axes([0.03, 0.858, 0.94, 0.003]); sep.set_facecolor(C_MID); sep.axis("off")
    plt.close(fig)

    # Stats box on fig2
    box = fig2.add_axes([0.08, 0.30, 0.84, 0.52])
    box.set_facecolor(C_LIGHT); box.set_xlim(0,1); box.set_ylim(0, len(stats) + 0.5)
    box.axis("off")

    for i, (label, value) in enumerate(reversed(stats)):
        y = i + 0.5
        row_bg = C_LIGHT if i % 2 == 0 else "white"
        box.add_patch(mpatches.FancyBboxPatch(
            (0, i), 1, 1, boxstyle="square,pad=0",
            facecolor=row_bg, edgecolor=C_MID, linewidth=0.5
        ))
        box.text(0.02, y, label, va="center", ha="left",
                 fontsize=10, color=C_TEXT)
        box.text(0.98, y, value, va="center", ha="right",
                 fontsize=10, fontweight="bold", color=C_ACCENT)

    box.add_patch(mpatches.FancyBboxPatch(
        (0, 0), 1, len(stats), boxstyle="square,pad=0",
        facecolor="none", edgecolor=C_ACCENT, linewidth=1.5
    ))

    pdf.savefig(fig2, bbox_inches="tight")
    plt.close(fig2)
    logger.info("  Page 1 (Overview) done")


def page_pipeline_flow(pdf: PdfPages, run_date: str, settings: dict, logger) -> None:
    """Page 2: vertical pipeline flow diagram with step boxes and settings."""
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))

    def _band(rect, fc):
        a = fig.add_axes(rect); a.set_facecolor(fc); a.axis("off"); return a

    h = _band([0, 0.91, 1, 0.09], C_HEADER)
    h.text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
           color="white", fontsize=10, fontweight="bold", va="center")
    h.text(0.98, 0.55, run_date, color="#BDC3C7", fontsize=8, va="center", ha="right")
    _band([0, 0.905, 1, 0.007], C_ACCENT)
    t = fig.add_axes([0.03, 0.865, 0.94, 0.04]); t.axis("off")
    t.text(0, 0.5, "2 · Pipeline: Steps & Settings Used",
           color=C_HEADER, fontsize=14, fontweight="bold", va="center")
    _band([0.03, 0.858, 0.94, 0.003], C_MID)

    # Main drawing area
    draw_ax = fig.add_axes([0.04, 0.06, 0.92, 0.79])
    draw_ax.set_xlim(0, 1)
    draw_ax.set_ylim(0, 1)
    draw_ax.axis("off")

    step_start = int(settings.get("step_start", 1))
    step_end   = int(settings.get("step_end", 5))

    of_species_tree_val = settings.get("of_species_tree")
    if of_species_tree_val:
        of_species_tree_str = str(of_species_tree_val)
    else:
        of_species_tree_str = "(infer from data)"

    steps = [
        {
            "num": 1,
            "name": "QC Filtering",
            "tool": "Custom QC (Python)",
            "color": C_PASS,
            "settings": [
                f"input_dir: {settings.get('input_dir', 'N/A')}",
                f"min_proteins: {settings.get('min_proteins', 'N/A')}",
            ],
            "output": "qc_results/passed_proteomes/",
        },
        {
            "num": 2,
            "name": "OrthoFinder",
            "tool": "OrthoFinder",
            "color": C_ACCENT,
            "settings": [
                f"-S {settings.get('of_search', 'N/A')}",
                f"-M {settings.get('of_gene_tree', 'N/A')}",
                f"-A {settings.get('of_msa', 'N/A')}",
                f"-s {of_species_tree_str}",
            ],
            "output": "orthofinder_results/Results_*/",
        },
        {
            "num": 3,
            "name": "Extract Orthologs",
            "tool": "Custom Extract (Python)",
            "color": C_ACCENT,
            "settings": [
                f"min_presence: {settings.get('min_presence', 'N/A')}",
            ],
            "output": "ortholog_results/gene_fastas/",
        },
        {
            "num": 4,
            "name": "Align & Trim",
            "tool": "MAFFT + trimAl",
            "color": C_ACCENT,
            "settings": [
                "MAFFT (--auto)",
                f"trimAl mode: {settings.get('trimal_mode', 'N/A')}",
            ],
            "output": "alignment_results/trimmed/",
        },
        {
            "num": 5,
            "name": "Phylogenetic Inference",
            "tool": "IQ-TREE",
            "color": "#8E44AD",
            "settings": [
                f"mode: {settings.get('iqtree_mode', 'N/A')}",
                f"model: {settings.get('iqtree_model', 'N/A')}",
                f"-B {settings.get('bootstrap', 'N/A')}",
                f"-alrt {settings.get('alrt', 'N/A')}",
            ],
            "output": "iqtree_results/species_tree.treefile",
        },
    ]

    n_steps = len(steps)
    # Vertical layout: spread steps evenly in the drawing area
    # Reserve bottom for footnote and possible warning
    top_y    = 0.97
    bottom_y = 0.08
    usable   = top_y - bottom_y
    slot_h   = usable / n_steps

    box_left   = 0.04
    box_width  = 0.44
    box_height = slot_h * 0.68

    arrow_x = box_left + box_width / 2

    for idx, step in enumerate(steps):
        # Vertical centre of this slot (from top downward)
        slot_centre = top_y - (idx + 0.5) * slot_h
        box_bottom  = slot_centre - box_height / 2
        box_top     = slot_centre + box_height / 2

        # Step box
        rect = mpatches.FancyBboxPatch(
            (box_left, box_bottom), box_width, box_height,
            boxstyle="round,pad=0.01",
            facecolor=step["color"], edgecolor="white", linewidth=1.5,
            transform=draw_ax.transAxes, clip_on=False,
            zorder=2,
        )
        draw_ax.add_patch(rect)

        # Step number circle
        circle_cx = box_left + 0.055
        circle_cy = slot_centre
        circle = plt.Circle(
            (circle_cx, circle_cy), 0.028,
            color="white", zorder=3, transform=draw_ax.transAxes,
            clip_on=False,
        )
        draw_ax.add_patch(circle)
        draw_ax.text(circle_cx, circle_cy, str(step["num"]),
                     ha="center", va="center", fontsize=9,
                     fontweight="bold", color=step["color"], zorder=4,
                     transform=draw_ax.transAxes)

        # Step name and tool
        text_x = box_left + 0.105
        draw_ax.text(text_x, slot_centre + 0.012, step["name"],
                     ha="left", va="center", fontsize=9,
                     fontweight="bold", color="white", zorder=4,
                     transform=draw_ax.transAxes)
        draw_ax.text(text_x, slot_centre - 0.018, step["tool"],
                     ha="left", va="center", fontsize=7.5,
                     color="#DDEEFF", style="italic", zorder=4,
                     transform=draw_ax.transAxes)

        # Settings panel on the right
        settings_x = box_left + box_width + 0.03
        settings_top = box_top - 0.005

        draw_ax.text(settings_x, settings_top, "Settings:",
                     ha="left", va="top", fontsize=7.5,
                     fontweight="bold", color=C_HEADER,
                     transform=draw_ax.transAxes)

        for si, s_line in enumerate(step["settings"]):
            draw_ax.text(settings_x + 0.01,
                         settings_top - 0.022 - si * 0.020,
                         f"• {s_line}",
                         ha="left", va="top", fontsize=7,
                         color=C_TEXT,
                         transform=draw_ax.transAxes)

        # Output directory in grey italic
        out_y = settings_top - 0.022 - len(step["settings"]) * 0.020 - 0.006
        draw_ax.text(settings_x, out_y, step["output"],
                     ha="left", va="top", fontsize=6.5,
                     color="#888888", style="italic",
                     transform=draw_ax.transAxes)

        # Downward arrow between steps
        if idx < n_steps - 1:
            next_slot_centre = top_y - (idx + 1.5) * slot_h
            arrow_y_start = slot_centre - box_height / 2 - 0.005
            arrow_y_end   = next_slot_centre + box_height / 2 + 0.005
            draw_ax.annotate(
                "", xy=(arrow_x, arrow_y_end),
                xytext=(arrow_x, arrow_y_start),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(
                    arrowstyle="-|>", color=C_MID,
                    lw=1.8, mutation_scale=14,
                ),
                zorder=1,
            )

    # Warning note if partial run
    if step_start > 1 or step_end < 5:
        draw_ax.text(0.5, 0.03,
                     f"Note: Pipeline ran steps {step_start}–{step_end} only.",
                     ha="center", va="bottom", fontsize=8,
                     color=C_WARN, fontweight="bold",
                     transform=draw_ax.transAxes)

    # Footnote at the bottom
    fig.text(0.5, 0.025,
             "Branch support: SH-aLRT / ultrafast bootstrap (UFBoot). "
             "Node label format: aLRT/UFBoot",
             ha="center", fontsize=7, color="#7F8C8D", style="italic")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Page 2 (Pipeline flow) done")


def page_qc(pdf: PdfPages, run_date: str, qc_df: pd.DataFrame, logger) -> None:
    verdicts = qc_df["verdict"].value_counts()

    # ---- Pie chart page ----
    fig, ax = new_page("3 · QC Results", run_date)
    plt.close(fig)

    fig2 = plt.figure(figsize=(PAGE_W, PAGE_H))
    for band, kw in [
        ([0, 0.91, 1, 0.09], dict(fc=C_HEADER)),
        ([0, 0.905, 1, 0.007], dict(fc=C_ACCENT)),
    ]:
        a = fig2.add_axes(band); a.set_facecolor(kw["fc"]); a.axis("off")
    fig2.axes[0].text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
                      color="white", fontsize=10, fontweight="bold", va="center")
    fig2.axes[0].text(0.98, 0.55, run_date, color="#BDC3C7", fontsize=8,
                      va="center", ha="right")
    tit2 = fig2.add_axes([0.03, 0.865, 0.94, 0.04]); tit2.axis("off")
    tit2.text(0, 0.5, "3 · QC Results", color=C_HEADER, fontsize=14,
              fontweight="bold", va="center")
    sep2 = fig2.add_axes([0.03, 0.858, 0.94, 0.003]); sep2.set_facecolor(C_MID); sep2.axis("off")

    # Bar chart of verdicts
    bar_ax = fig2.add_axes([0.08, 0.62, 0.40, 0.22])
    cats   = ["PASS", "WARN", "FAIL"]
    colors = [C_PASS, C_WARN, C_FAIL]
    values = [verdicts.get(c, 0) for c in cats]
    bars = bar_ax.bar(cats, values, color=colors, edgecolor="white", width=0.5)
    for bar, val in zip(bars, values):
        bar_ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                    str(val), ha="center", va="bottom", fontweight="bold", fontsize=11)
    bar_ax.set_ylabel("Number of genomes", fontsize=9)
    bar_ax.set_title("Genome QC Verdicts", fontsize=11, fontweight="bold", color=C_HEADER)
    bar_ax.spines[["top", "right"]].set_visible(False)
    bar_ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # Summary text
    txt_ax = fig2.add_axes([0.55, 0.62, 0.38, 0.22])
    txt_ax.axis("off")
    total = len(qc_df)
    n_inc = verdicts.get("PASS", 0) + verdicts.get("WARN", 0)
    txt_ax.text(0.05, 0.85, f"Total genomes:      {total}", fontsize=10, color=C_TEXT)
    txt_ax.text(0.05, 0.70, f"Included (PASS):    {verdicts.get('PASS',0)}", fontsize=10, color=C_PASS, fontweight="bold")
    txt_ax.text(0.05, 0.55, f"Included (WARN):    {verdicts.get('WARN',0)}", fontsize=10, color=C_WARN, fontweight="bold")
    txt_ax.text(0.05, 0.40, f"Excluded (FAIL):    {verdicts.get('FAIL',0)}", fontsize=10, color=C_FAIL, fontweight="bold")
    txt_ax.text(0.05, 0.20, "WARN genomes have one flagged\nmetric but are included.",
                fontsize=8, color="#7F8C8D", style="italic")

    # Full QC table
    tbl_ax = fig2.add_axes([0.04, 0.04, 0.92, 0.54])
    tbl_ax.axis("off")

    disp_cols = ["genome", "n_proteins", "median_len_aa", "pct_short", "pct_dup", "verdict", "reason"]
    existing  = [c for c in disp_cols if c in qc_df.columns]
    disp_df   = qc_df[existing].copy()
    disp_df.columns = [c.replace("_", " ").title() for c in existing]

    n = len(disp_df)
    fontsize = max(5, min(8, 200 // max(n, 1)))

    cell_text = disp_df.astype(str).values.tolist()
    verdict_col = [c.upper() for c in list(disp_df.get("Verdict", pd.Series()).values)] if "Verdict" in disp_df.columns else []
    color_map = {"PASS": C_PASS, "WARN": C_WARN, "FAIL": C_FAIL}

    tbl = tbl_ax.table(
        cellText=cell_text,
        colLabels=list(disp_df.columns),
        cellLoc="left",
        loc="upper left",
        bbox=[0, 0, 1, 1],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(C_MID)
        if r == 0:
            cell.set_facecolor(C_HEADER)
            cell.set_text_props(color="white", fontweight="bold")
        else:
            row_bg = C_LIGHT if r % 2 == 1 else "white"
            cell.set_facecolor(row_bg)
            row_verdict = disp_df.iloc[r - 1].get("Verdict", "") if "Verdict" in disp_df.columns else ""
            if c == list(disp_df.columns).index("Verdict") if "Verdict" in disp_df.columns else -1:
                col = color_map.get(str(row_verdict).upper(), C_TEXT)
                cell.set_text_props(color=col, fontweight="bold")
            else:
                cell.set_text_props(color=C_TEXT)

    pdf.savefig(fig2, bbox_inches="tight")
    plt.close(fig2)
    logger.info("  Page 3 (QC) done")


def page_orthologs(pdf: PdfPages, run_date: str,
                   summary: dict, included_df: pd.DataFrame | None,
                   logger) -> None:
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))

    def _band(rect, fc):
        a = fig.add_axes(rect); a.set_facecolor(fc); a.axis("off")
        return a

    h = _band([0, 0.91, 1, 0.09], C_HEADER)
    h.text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
           color="white", fontsize=10, fontweight="bold", va="center")
    h.text(0.98, 0.55, run_date, color="#BDC3C7", fontsize=8, va="center", ha="right")
    _band([0, 0.905, 1, 0.007], C_ACCENT)
    t = fig.add_axes([0.03, 0.865, 0.94, 0.04]); t.axis("off")
    t.text(0, 0.5, "4 · Ortholog Selection", color=C_HEADER, fontsize=14,
           fontweight="bold", va="center")
    _band([0.03, 0.858, 0.94, 0.003], C_MID)

    if not summary and included_df is None:
        ax = fig.add_axes([0.1, 0.4, 0.8, 0.4]); ax.axis("off")
        ax.text(0.5, 0.5, "Ortholog data not available.", ha="center",
                va="center", color=C_FAIL, fontsize=12)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig); return

    # Summary stats box
    stats_pairs = [
        ("Species in analysis",           summary.get("total_species", "N/A")),
        ("Total orthogroups assessed",     f"{int(summary.get('total_orthogroups_assessed', 0)):,}"),
        ("Strict filter (100% presence)", f"{int(summary.get('orthogroups_passing_strict', 0)):,}"),
        ("Orthologs selected",            f"{int(summary.get('orthogroups_selected', 0)):,}"),
        ("Total aa in selected set",       f"{int(summary.get('total_aa_in_concat', 0)):,}"),
    ]
    sp = fig.add_axes([0.05, 0.68, 0.90, 0.17])
    sp.set_facecolor(C_LIGHT); sp.set_xlim(0,1); sp.set_ylim(0, len(stats_pairs))
    sp.axis("off")
    for i, (lbl, val) in enumerate(reversed(stats_pairs)):
        y = i + 0.5
        sp.add_patch(mpatches.Rectangle((0, i), 1, 1,
            facecolor=C_LIGHT if i % 2 == 0 else "white",
            edgecolor=C_MID, linewidth=0.5))
        sp.text(0.015, y, lbl, va="center", fontsize=9, color=C_TEXT)
        sp.text(0.985, y, str(val), va="center", ha="right", fontsize=9,
                fontweight="bold", color=C_ACCENT)
    sp.add_patch(mpatches.Rectangle((0,0), 1, len(stats_pairs),
        facecolor="none", edgecolor=C_ACCENT, linewidth=1.2))

    # Bar chart: per-gene completeness distribution
    if included_df is not None and "pct_species_present" in included_df.columns:
        bax = fig.add_axes([0.08, 0.34, 0.84, 0.28])
        pcts = included_df["pct_species_present"].dropna()
        bins = list(range(50, 105, 5))
        bax.hist(pcts, bins=bins, color=C_ACCENT, edgecolor="white", rwidth=0.85)
        bax.set_xlabel("% species with ortholog", fontsize=9)
        bax.set_ylabel("Number of orthologs", fontsize=9)
        bax.set_title("Per-Ortholog Species Completeness", fontsize=10,
                       fontweight="bold", color=C_HEADER)
        bax.spines[["top", "right"]].set_visible(False)
        bax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # Length stats table
    if included_df is not None and not included_df.empty:
        tbl_cols = ["orthogroup_id", "n_species_present", "pct_species_present",
                    "mean_seq_length_aa", "min_seq_length_aa", "max_seq_length_aa"]
        show = included_df[[c for c in tbl_cols if c in included_df.columns]].head(20)
        show.columns = [c.replace("_", " ").replace("aa", "(aa)").title()
                        for c in show.columns]
        show = show.round(1).astype(str)
        tbl_ax = fig.add_axes([0.04, 0.04, 0.92, 0.27])
        tbl_ax.axis("off")
        n = len(show); fs = max(5, min(7, 150 // max(n, 1)))
        tbl = tbl_ax.table(cellText=show.values.tolist(),
                           colLabels=list(show.columns),
                           cellLoc="left", loc="upper left",
                           bbox=[0, 0, 1, 1])
        tbl.auto_set_font_size(False); tbl.set_fontsize(fs)
        for (r, c), cell in tbl.get_celld().items():
            cell.set_edgecolor(C_MID)
            if r == 0:
                cell.set_facecolor(C_HEADER)
                cell.set_text_props(color="white", fontweight="bold")
            else:
                cell.set_facecolor(C_LIGHT if r % 2 == 1 else "white")
                cell.set_text_props(color=C_TEXT)
        if len(included_df) > 20:
            fig.text(0.5, 0.025, f"(showing first 20 of {len(included_df)} orthologs)",
                     ha="center", fontsize=7, color="#7F8C8D", style="italic")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Page 4 (Ortholog selection) done")


def page_functional(pdf: PdfPages, run_date: str,
                    og_to_desc: dict[str, str], logger) -> None:
    if not og_to_desc:
        logger.warning("  Skipping functional page — no gene annotations loaded")
        return

    counts = categorize_genes(og_to_desc)
    total  = sum(counts.values())

    fig = plt.figure(figsize=(PAGE_W, PAGE_H))

    def _band(rect, fc):
        a = fig.add_axes(rect); a.set_facecolor(fc); a.axis("off"); return a

    h = _band([0, 0.91, 1, 0.09], C_HEADER)
    h.text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
           color="white", fontsize=10, fontweight="bold", va="center")
    h.text(0.98, 0.55, run_date, color="#BDC3C7", fontsize=8, va="center", ha="right")
    _band([0, 0.905, 1, 0.007], C_ACCENT)
    t = fig.add_axes([0.03, 0.865, 0.94, 0.04]); t.axis("off")
    t.text(0, 0.5, "5 · Functional Overview of Selected Orthologs",
           color=C_HEADER, fontsize=14, fontweight="bold", va="center")
    _band([0.03, 0.858, 0.94, 0.003], C_MID)

    # Horizontal bar chart
    cats_ordered = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    labels = [c for c, _ in cats_ordered]
    values = [v for _, v in cats_ordered]

    palette = plt.cm.Blues(np.linspace(0.4, 0.85, len(labels)))

    bar_ax = fig.add_axes([0.36, 0.20, 0.58, 0.62])
    bars = bar_ax.barh(range(len(labels)), values, color=palette,
                       edgecolor="white", height=0.65)
    for bar, val in zip(bars, values):
        pct = 100 * val / total if total else 0
        bar_ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
                    f"{val}  ({pct:.1f}%)", va="center", fontsize=8, color=C_TEXT)
    bar_ax.set_yticks(range(len(labels)))
    bar_ax.set_yticklabels(labels, fontsize=8)
    bar_ax.set_xlabel("Number of orthologs", fontsize=9)
    bar_ax.set_title(f"Functional Categories  (n = {total} orthologs)",
                      fontsize=10, fontweight="bold", color=C_HEADER)
    bar_ax.spines[["top", "right"]].set_visible(False)
    bar_ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    bar_ax.set_xlim(0, max(values) * 1.35 if values else 1)

    # Note
    note_ax = fig.add_axes([0.04, 0.10, 0.28, 0.72])
    note_ax.axis("off")
    note_ax.text(0, 1.0,
        "Categories are assigned by\n"
        "keyword matching on gene\n"
        "product descriptions from\n"
        "the input .faa headers.\n\n"
        "These are broad, indicative\n"
        "groupings — not COG/KEGG\n"
        "annotations.\n\n"
        "Genes described as\n"
        "'hypothetical protein' are\n"
        "placed in Hypothetical /\n"
        "Unknown.\n\n"
        "Most selected orthologs\n"
        "should fall into core\n"
        "housekeeping categories\n"
        "(ribosomal proteins, DNA\n"
        "replication, etc.), which\n"
        "is the expected outcome\n"
        "for single-copy universal\n"
        "orthologs.",
        va="top", ha="left", fontsize=7.5, color="#555555",
        style="italic", transform=note_ax.transAxes,
        wrap=True
    )

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Page 5 (Functional overview) done")


def page_alignment(pdf: PdfPages, run_date: str,
                   aln_df: pd.DataFrame, logger) -> None:
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))

    def _band(rect, fc):
        a = fig.add_axes(rect); a.set_facecolor(fc); a.axis("off"); return a

    h = _band([0, 0.91, 1, 0.09], C_HEADER)
    h.text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
           color="white", fontsize=10, fontweight="bold", va="center")
    h.text(0.98, 0.55, run_date, color="#BDC3C7", fontsize=8, va="center", ha="right")
    _band([0, 0.905, 1, 0.007], C_ACCENT)
    t = fig.add_axes([0.03, 0.865, 0.94, 0.04]); t.axis("off")
    t.text(0, 0.5, "6 · Alignment Statistics", color=C_HEADER, fontsize=14,
           fontweight="bold", va="center")
    _band([0.03, 0.858, 0.94, 0.003], C_MID)

    success = aln_df[aln_df["status"] == "success"]
    failed  = aln_df[aln_df["status"] != "success"]

    total_raw     = int(success["raw_length_aa"].sum()) if not success.empty else 0
    total_trimmed = int(success["trimmed_length_aa"].sum()) if not success.empty else 0
    pct_ret       = 100.0 * total_trimmed / total_raw if total_raw else 0
    mean_ret      = success["pct_retained"].mean() if not success.empty else 0

    stats = [
        ("Orthologs successfully aligned",  f"{len(success):,}"),
        ("Orthologs failed alignment/trim",  f"{len(failed):,}"),
        ("Total alignment length (raw)",      f"{total_raw:,} aa"),
        ("Total alignment length (trimmed)",  f"{total_trimmed:,} aa"),
        ("Columns removed by trimming",       f"{total_raw - total_trimmed:,} aa  ({100-pct_ret:.1f}%)"),
        ("Mean column retention per gene",    f"{mean_ret:.1f}%"),
    ]

    sp = fig.add_axes([0.05, 0.66, 0.90, 0.18])
    sp.set_facecolor(C_LIGHT); sp.set_xlim(0,1); sp.set_ylim(0, len(stats))
    sp.axis("off")
    for i, (lbl, val) in enumerate(reversed(stats)):
        y = i + 0.5
        sp.add_patch(mpatches.Rectangle((0, i), 1, 1,
            facecolor=C_LIGHT if i % 2 == 0 else "white",
            edgecolor=C_MID, linewidth=0.5))
        sp.text(0.015, y, lbl, va="center", fontsize=9, color=C_TEXT)
        sp.text(0.985, y, val, va="center", ha="right", fontsize=9,
                fontweight="bold", color=C_ACCENT)
    sp.add_patch(mpatches.Rectangle((0,0), 1, len(stats),
        facecolor="none", edgecolor=C_ACCENT, linewidth=1.2))

    # Before/after bar comparison
    if total_raw > 0:
        comp_ax = fig.add_axes([0.08, 0.40, 0.38, 0.22])
        comp_ax.bar(["Raw", "Trimmed"], [total_raw, total_trimmed],
                    color=[C_WARN, C_PASS], edgecolor="white", width=0.4)
        comp_ax.set_ylabel("Total alignment length (aa)", fontsize=8)
        comp_ax.set_title("Before vs After Trimming", fontsize=9,
                           fontweight="bold", color=C_HEADER)
        for val, x in zip([total_raw, total_trimmed], [0, 1]):
            comp_ax.text(x, val * 1.01, f"{val:,}", ha="center",
                         va="bottom", fontsize=7, fontweight="bold")
        comp_ax.spines[["top", "right"]].set_visible(False)
        comp_ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    # Histogram of per-gene column retention
    if not success.empty and "pct_retained" in success.columns:
        hist_ax = fig.add_axes([0.58, 0.40, 0.38, 0.22])
        hist_ax.hist(success["pct_retained"].dropna(), bins=20,
                     color=C_ACCENT, edgecolor="white", rwidth=0.9)
        hist_ax.axvline(50, color=C_FAIL, linestyle="--", linewidth=1,
                        label=">50% removed (warning)")
        hist_ax.set_xlabel("% columns retained", fontsize=8)
        hist_ax.set_ylabel("Number of orthologs", fontsize=8)
        hist_ax.set_title("Per-Gene Column Retention", fontsize=9,
                           fontweight="bold", color=C_HEADER)
        hist_ax.legend(fontsize=6)
        hist_ax.spines[["top", "right"]].set_visible(False)

    # Bottom: per-gene table (top/worst retained)
    tbl_ax = fig.add_axes([0.04, 0.04, 0.92, 0.33])
    tbl_ax.axis("off")

    show = success.nsmallest(20, "pct_retained")[
        ["orthogroup", "raw_length_aa", "trimmed_length_aa", "pct_retained"]
    ].copy() if not success.empty else pd.DataFrame()

    if not show.empty:
        show.columns = ["Orthogroup", "Raw (aa)", "Trimmed (aa)", "% Retained"]
        show = show.round(1).astype(str)
        fs = max(5, min(7, 200 // max(len(show), 1)))
        tbl = tbl_ax.table(cellText=show.values.tolist(),
                           colLabels=list(show.columns),
                           cellLoc="left", loc="upper left",
                           bbox=[0, 0, 1, 1])
        tbl.auto_set_font_size(False); tbl.set_fontsize(fs)
        for (r, c), cell in tbl.get_celld().items():
            cell.set_edgecolor(C_MID)
            if r == 0:
                cell.set_facecolor(C_HEADER)
                cell.set_text_props(color="white", fontweight="bold")
            else:
                cell.set_facecolor(C_LIGHT if r % 2 == 1 else "white")
                cell.set_text_props(color=C_TEXT)
        tbl_ax.set_title("Orthologs with lowest column retention (top 20)",
                          fontsize=8, color=C_TEXT, pad=2)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Page 6 (Alignment stats) done")


def page_tree(pdf: PdfPages, run_date: str, treefile: Path, logger) -> None:
    if not BIO_AVAILABLE:
        logger.warning("  Skipping tree page — BioPython not available")
        return
    if not treefile.exists():
        logger.warning(f"  Skipping tree page — treefile not found: {treefile}")
        return

    try:
        tree = Phylo.read(str(treefile), "newick")
    except Exception as e:
        logger.warning(f"  Could not parse treefile: {e}")
        return

    terminals = tree.get_terminals()
    n_tips = len(terminals)

    # Scale figure height to number of tips (min 6, max 60 inches)
    fig_h = max(6.0, min(60.0, n_tips * 0.35))
    fig = plt.figure(figsize=(PAGE_W, fig_h))

    # Header band scaled to figure
    hdr_frac = 0.7 / fig_h
    h = fig.add_axes([0, 1 - hdr_frac, 1, hdr_frac])
    h.set_facecolor(C_HEADER); h.set_xlim(0,1); h.set_ylim(0,1); h.axis("off")
    h.text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
           color="white", fontsize=10, fontweight="bold", va="center")
    h.text(0.98, 0.55, run_date, color="#BDC3C7", fontsize=8, va="center", ha="right")
    acc = fig.add_axes([0, 1 - hdr_frac - 0.005/fig_h*11, 1, 0.005*11/fig_h])
    acc.set_facecolor(C_ACCENT); acc.axis("off")

    title_frac = 0.4 / fig_h
    t = fig.add_axes([0.03, 1 - hdr_frac - title_frac - 0.01, 0.94, title_frac])
    t.axis("off")
    t.text(0, 0.5, "7 · Species Tree (IQ-TREE3)", color=C_HEADER, fontsize=14,
           fontweight="bold", va="center")

    # Newick string — wrap before determining height
    newick_str = treefile.read_text().strip()
    wrapped_nwk = textwrap.fill(newick_str, width=95,
                                break_on_hyphens=False, break_long_words=True)

    # Tree axes
    tree_top    = 1 - hdr_frac - title_frac - 0.05
    newick_frac = max(0.8, len(wrapped_nwk.splitlines()) * 0.15) / fig_h
    tree_bottom = newick_frac + 0.03
    tree_ax = fig.add_axes([0.04, tree_bottom, 0.92, tree_top - tree_bottom])

    Phylo.draw(tree, axes=tree_ax, do_show=False)
    tree_ax.set_title(
        f"{n_tips} taxa  ·  "
        "node labels = SH-aLRT / UFBoot  ·  "
        "branch lengths = substitutions per site",
        fontsize=7, color="#555555", pad=4
    )
    tree_ax.spines[["top", "right"]].set_visible(False)

    # Newick string at bottom
    newick_ax = fig.add_axes([0.04, 0.005, 0.92, newick_frac - 0.015])
    newick_ax.axis("off")
    newick_ax.text(0, 1.0, "Newick string:", fontsize=7, fontweight="bold",
                   color=C_HEADER, va="top")
    newick_ax.text(0, 0.85, wrapped_nwk, fontsize=5.0, color="#444444",
                   va="top", family="monospace",
                   transform=newick_ax.transAxes)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Page 7 (Species tree) done")


def page_references(pdf: PdfPages, run_date: str, logger) -> None:
    """Page 8: References for all tools used in the pipeline."""
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))

    def _band(rect, fc):
        a = fig.add_axes(rect); a.set_facecolor(fc); a.axis("off"); return a

    h = _band([0, 0.91, 1, 0.09], C_HEADER)
    h.text(0.02, 0.55, "OrthoFinder3 · Species Tree Pipeline — Summary Report",
           color="white", fontsize=10, fontweight="bold", va="center")
    h.text(0.98, 0.55, run_date, color="#BDC3C7", fontsize=8, va="center", ha="right")
    _band([0, 0.905, 1, 0.007], C_ACCENT)
    t = fig.add_axes([0.03, 0.865, 0.94, 0.04]); t.axis("off")
    t.text(0, 0.5, "8 · References", color=C_HEADER, fontsize=14,
           fontweight="bold", va="center")
    _band([0.03, 0.858, 0.94, 0.003], C_MID)

    ref_ax = fig.add_axes([0.05, 0.08, 0.90, 0.77])
    ref_ax.axis("off")
    ref_ax.set_xlim(0, 1)
    ref_ax.set_ylim(0, 1)

    note_text = (
        "Please cite the following tools when reporting results from this "
        "pipeline in a publication."
    )
    ref_ax.text(0, 1.0, note_text, ha="left", va="top", fontsize=8.5,
                color=C_TEXT, style="italic", transform=ref_ax.transAxes)

    references = [
        ("[1]",
         "Emms DM & Kelly S (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics.",
         "Genome Biology 20:238.",
         "https://doi.org/10.1186/s13059-019-1832-y"),
        ("[2]",
         "Buchfink B, Xie C & Huson DH (2015) Fast and sensitive protein alignment using DIAMOND.",
         "Nature Methods 12:59-60.",
         "https://doi.org/10.1038/nmeth.3176"),
        ("[3]",
         "Katoh K & Standley DM (2013) MAFFT Multiple Sequence Alignment Software Version 7: "
         "Improvements in Performance and Usability.",
         "Molecular Biology and Evolution 30(4):772-780.",
         "https://doi.org/10.1093/molbev/mst010"),
        ("[4]",
         "Capella-Gutierrez S, Silla-Martinez JM & Gabaldon T (2009) trimAl: a tool for automated "
         "alignment trimming in large-scale phylogenetic analyses.",
         "Bioinformatics 25(15):1972-1973.",
         "https://doi.org/10.1093/bioinformatics/btp348"),
        ("[5]",
         "Minh BQ et al. (2020) IQ-TREE 2: New Models and Methods for Phylogenetic Inference.",
         "Molecular Biology and Evolution 37(5):1530-1534.",
         "https://doi.org/10.1093/molbev/msaa015"),
        ("[6]",
         "Kalyaanamoorthy S et al. (2017) ModelFinder: fast model selection for accurate "
         "phylogenetic estimates.",
         "Nature Methods 14:587-589.  [cite when using MFP or MFP+MERGE model]",
         "https://doi.org/10.1038/nmeth.4285"),
        ("[7]",
         "Hoang DT et al. (2018) UFBoot2: Improving the Ultrafast Bootstrap Approximation.",
         "Molecular Biology and Evolution 35(2):518-522.",
         "https://doi.org/10.1093/molbev/msx281"),
        ("[8]",
         "Cock PJA et al. (2009) Biopython: freely available Python tools for computational "
         "molecular biology and bioinformatics.",
         "Bioinformatics 25(11):1422-1423.",
         "https://doi.org/10.1093/bioinformatics/btp163"),
    ]

    # Layout: start below the note, each entry ~0.082 of axes height
    entry_h  = 0.082
    start_y  = 0.92
    indent   = 0.045

    for i, (num, authors_title, journal, doi) in enumerate(references):
        y_top = start_y - i * entry_h

        # Reference number
        ref_ax.text(0, y_top, num, ha="left", va="top", fontsize=8,
                    fontweight="bold", color=C_ACCENT,
                    transform=ref_ax.transAxes)

        # Authors + title line
        ref_ax.text(indent, y_top, authors_title,
                    ha="left", va="top", fontsize=7.5, color=C_TEXT,
                    transform=ref_ax.transAxes)

        # Journal line
        ref_ax.text(indent, y_top - 0.027, journal,
                    ha="left", va="top", fontsize=7.5, color=C_TEXT,
                    style="italic", transform=ref_ax.transAxes)

        # DOI line
        ref_ax.text(indent, y_top - 0.052, doi,
                    ha="left", va="top", fontsize=7, color="#2471A3",
                    transform=ref_ax.transAxes)

    # Suggested Methods & Materials text
    mm_y = start_y - len(references) * entry_h - 0.04
    ref_ax.text(0, mm_y, "Suggested Materials & Methods text:",
                ha="left", va="top", fontsize=9,
                fontweight="bold", color=C_HEADER,
                transform=ref_ax.transAxes)

    mm_para = (
        "Proteome sequences were quality-filtered using custom QC criteria before ortholog "
        "identification with OrthoFinder [1] using DIAMOND [2] for sequence search. "
        "Single-copy orthologs present in the required fraction of taxa were extracted, "
        "aligned with MAFFT [3] (--auto), and trimmed with trimAl [4]. "
        "Concatenated alignments were used to infer a maximum-likelihood species tree with "
        "IQ-TREE 2 [5] using ModelFinder [6] for model selection, ultrafast bootstrap (UFBoot2) [7], "
        "and SH-aLRT branch support; Biopython [8] was used for result parsing."
    )
    wrapped_mm = textwrap.fill(mm_para, width=110,
                               break_on_hyphens=False, break_long_words=False)
    ref_ax.text(0, mm_y - 0.040, wrapped_mm,
                ha="left", va="top", fontsize=7.5, color=C_TEXT,
                transform=ref_ax.transAxes)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Page 8 (References) done")


# ── CLI & main ────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate a PDF summary report of all pipeline results."
    )
    p.add_argument("--output_dir",      default="report",
                   help="Directory to write the PDF (default: report)")
    p.add_argument("--qc_dir",          default="qc_results")
    p.add_argument("--ortholog_dir",    default="ortholog_results")
    p.add_argument("--alignment_dir",   default="alignment_results")
    p.add_argument("--iqtree_dir",      default="iqtree_results")
    p.add_argument("--orthofinder_results", default="orthofinder_results",
                   help="Parent directory of OrthoFinder Results_*/ (default: orthofinder_results)")
    p.add_argument("--faa_dir",         default="qc_results/passed_proteomes",
                   help="Directory with original .faa files for gene annotation lookup")
    p.add_argument("--prefix",          default="species_tree",
                   help="IQ-TREE output prefix (default: species_tree)")
    p.add_argument("--output_file",     default="pipeline_summary.pdf",
                   help="PDF filename (default: pipeline_summary.pdf)")
    p.add_argument("--settings_file",   default="run_settings.json",
                   help="JSON file with pipeline run settings (written automatically by pipeline_runner.sh)")
    return p.parse_args()


def main() -> None:
    start_time = datetime.now()
    args = parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = str(out_dir / "report.log")
    logger   = setup_logger(log_path)
    log_header(logger, args, start_time)

    run_date = start_time.strftime("%Y-%m-%d %H:%M")
    pdf_path = out_dir / args.output_file

    logger.info("Loading pipeline output data ...")

    qc_df        = load_qc_summary(Path(args.qc_dir), logger)
    summary, included_df, included_ogs = load_ortholog_data(
        Path(args.ortholog_dir), logger
    )
    aln_df       = load_alignment_stats(Path(args.alignment_dir), logger)
    of_dir       = find_orthofinder_dir(Path(args.orthofinder_results), logger)
    treefile     = Path(args.iqtree_dir) / f"{args.prefix}.treefile"
    settings     = load_run_settings(Path(args.settings_file), logger)

    og_to_desc = load_gene_annotations(
        Path(args.faa_dir), of_dir, included_ogs, logger
    )

    logger.info(f"Writing PDF: {pdf_path}")

    with PdfPages(str(pdf_path)) as pdf:
        # PDF metadata
        d = pdf.infodict()
        d["Title"]   = "OrthoFinder3 Species Tree Pipeline — Summary Report"
        d["Author"]  = "06_generate_report.py"
        d["Subject"] = "Phylogenomics pipeline summary"
        d["CreationDate"] = start_time

        logger.info("Generating pages ...")
        page_overview(pdf, run_date, qc_df, summary, aln_df, logger)
        logger.info("  Page 1 (Overview) written")

        page_pipeline_flow(pdf, run_date, settings, logger)

        if qc_df is not None:
            page_qc(pdf, run_date, qc_df, logger)
        else:
            logger.warning("  Skipping QC page — data not available")

        page_orthologs(pdf, run_date, summary, included_df, logger)
        page_functional(pdf, run_date, og_to_desc, logger)

        if aln_df is not None:
            page_alignment(pdf, run_date, aln_df, logger)
        else:
            logger.warning("  Skipping alignment page — data not available")

        page_tree(pdf, run_date, treefile, logger)
        page_references(pdf, run_date, logger)

    logger.info(f"Report written to: {pdf_path.resolve()}")
    print(f"\nReport: {pdf_path.resolve()}")

    log_footer(logger, start_time, "SUCCESS")


if __name__ == "__main__":
    main()
