#!/usr/bin/env python3
"""
03_extract_orthologs.py — Extract and reformat single-copy orthologs

Reads OrthoFinder results, selects single-copy orthologs based on presence
threshold, and produces:
  - Per-gene FASTA files (gene_fastas/)
  - Concatenated archive FASTA (species_tree_input/)
  - Summary tables (tables/)
  - List of included orthologs (included_orthologs.txt)
"""

import argparse
import logging
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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


def log_header(logger: logging.Logger, args: argparse.Namespace, start_time: datetime) -> None:
    try:
        import Bio
        biopython_ver = Bio.__version__
    except Exception:
        biopython_ver = "unknown"
    try:
        pandas_ver = pd.__version__
    except Exception:
        pandas_ver = "unknown"

    logger.info("=" * 54)
    logger.info("Script:      03_extract_orthologs.py")
    logger.info(f"Started:     {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Parameters:  --orthofinder_dir {args.orthofinder_dir}")
    logger.info(f"             --faa_dir {args.faa_dir}")
    logger.info(f"             --min_presence {args.min_presence}")
    logger.info(f"             --output_dir {args.output_dir}")
    logger.info(f"Tool versions: biopython {biopython_ver}, pandas {pandas_ver}")
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


def sanitise_name(name: str) -> str:
    """Sanitise a species name for use in FASTA headers and IQ-TREE."""
    name = name.replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9_\-]", "", name)
    return name


def build_sequence_index(
    faa_dir: Path, logger: logging.Logger
) -> tuple[dict[str, SeqRecord], dict[str, str], dict[str, str]]:
    """
    Build three mappings:
      locus_tag  -> SeqRecord  (for sequence retrieval)
      locus_tag  -> species    (for species attribution)
      file_stem  -> species    (for gc_df column mapping)

    By reading each .faa file and noting the source filename, we get a direct,
    unambiguous locus-tag-to-species link without relying on prefix matching.
    """
    logger.info(f"Building sequence index from {faa_dir} ...")
    locus_to_seq: dict[str, SeqRecord] = {}
    locus_to_species: dict[str, str] = {}
    stem_to_species: dict[str, str] = {}

    faa_files = sorted(faa_dir.glob("*.faa"))
    if not faa_files:
        logger.error(f"No .faa files found in {faa_dir}")
        return locus_to_seq, locus_to_species, stem_to_species

    file_iter = (
        tqdm(faa_files, desc="Indexing genomes", unit="genome")
        if TQDM_AVAILABLE else faa_files
    )
    for faa_path in file_iter:
        species = sanitise_name(faa_path.stem)
        stem_to_species[faa_path.stem] = species
        count = 0
        for record in SeqIO.parse(str(faa_path), "fasta"):
            locus_to_seq[record.id] = record
            locus_to_species[record.id] = species
            count += 1
        logger.debug(f"  Indexed {count} sequences from {faa_path.name} → species '{species}'")

    logger.info(
        f"Indexed {len(locus_to_seq)} total sequences "
        f"from {len(faa_files)} genomes"
    )
    return locus_to_seq, locus_to_species, stem_to_species


def load_gene_count_matrix(orthofinder_dir: Path, logger: logging.Logger) -> pd.DataFrame:
    gc_path = orthofinder_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    if not gc_path.exists():
        logger.error(f"Gene count matrix not found: {gc_path}")
        sys.exit(1)
    df = pd.read_csv(gc_path, sep="\t", index_col=0)
    if "Total" in df.columns:
        df = df.drop(columns=["Total"])
    logger.info(
        f"Loaded gene count matrix: "
        f"{df.shape[0]} orthogroups × {df.shape[1]} species"
    )
    return df


def load_orthogroups_txt(
    orthofinder_dir: Path, logger: logging.Logger
) -> dict[str, list[str]]:
    """
    Parse Orthogroups/Orthogroups.txt into:
      orthogroup_id -> [locus_tag, ...]
    """
    og_path = orthofinder_dir / "Orthogroups" / "Orthogroups.txt"
    if not og_path.exists():
        logger.error(f"Orthogroups.txt not found: {og_path}")
        sys.exit(1)

    og_dict: dict[str, list[str]] = {}
    with open(og_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(":", 1)
            if len(parts) != 2:
                continue
            og_id = parts[0].strip()
            tags = [t.strip() for t in parts[1].split() if t.strip()]
            og_dict[og_id] = tags

    logger.info(f"Loaded {len(og_dict)} orthogroups from Orthogroups.txt")
    return og_dict


def select_orthologs(
    gc_df: pd.DataFrame,
    min_presence: float,
    logger: logging.Logger,
) -> tuple[list[str], list[str]]:
    """
    Return (strict_list, relaxed_list) of orthogroup IDs.
    strict:  present in ALL species, exactly 1 copy each
    relaxed: present in >= min_presence fraction, exactly 1 copy where present
    """
    n_species = gc_df.shape[1]
    min_count = n_species if min_presence >= 1.0 else int(n_species * min_presence)

    strict_list: list[str] = []
    relaxed_list: list[str] = []

    for og_id in gc_df.index:
        counts = gc_df.loc[og_id].values
        n_present = int((counts > 0).sum())
        n_multi = int((counts > 1).sum())

        if n_multi > 0:
            continue

        if n_present == n_species:
            strict_list.append(og_id)

        if n_present >= min_count:
            relaxed_list.append(og_id)

    logger.info(
        f"Strict filter (all {n_species} species, 1 copy): "
        f"{len(strict_list)} orthogroups"
    )
    logger.info(
        f"Relaxed filter (≥{min_presence*100:.0f}% = "
        f"{min_count}/{n_species} species, no multi-copy): "
        f"{len(relaxed_list)} orthogroups"
    )
    return strict_list, relaxed_list


def build_og_seq_map(
    selected: list[str],
    og_txt: dict[str, list[str]],
    gc_df: pd.DataFrame,
    locus_to_seq: dict[str, SeqRecord],
    locus_to_species: dict[str, str],
    stem_to_species: dict[str, str],
    logger: logging.Logger,
) -> dict[str, dict[str, SeqRecord]]:
    """
    For each selected orthogroup, build {species_name: SeqRecord}.
    Uses the direct locus_to_species mapping (no prefix guessing).
    gc_df column names are file stems; map them to species names first.
    """
    col_to_species: dict[str, str] = {}
    for col in gc_df.columns:
        col_to_species[col] = stem_to_species.get(col, sanitise_name(col))

    og_seq_map: dict[str, dict[str, SeqRecord]] = {}

    for og_id in selected:
        tags = og_txt.get(og_id, [])
        species_seqs: dict[str, SeqRecord] = {}

        for tag in tags:
            species = locus_to_species.get(tag)
            if species is None:
                logger.warning(
                    f"{og_id}: locus tag '{tag}' not found in sequence index — skipping"
                )
                continue
            record = locus_to_seq[tag]
            species_seqs[species] = SeqRecord(
                seq=record.seq, id=species, description=""
            )

        og_seq_map[og_id] = species_seqs

    return og_seq_map


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract single-copy orthologs from OrthoFinder results."
    )
    parser.add_argument("--orthofinder_dir", required=True,
                        help="Path to OrthoFinder Results_*/ folder")
    parser.add_argument("--faa_dir", required=True,
                        help="Directory with original .faa files")
    parser.add_argument("--min_presence", type=float, default=1.0,
                        help="Min fraction of species with the gene (default: 1.0 = strict)")
    parser.add_argument("--output_dir", default="ortholog_results",
                        help="Output directory (default: ortholog_results)")
    return parser.parse_args()


def main() -> None:
    start_time = datetime.now()
    args = parse_args()

    out_dir = Path(args.output_dir)
    gene_fasta_dir = out_dir / "gene_fastas"
    species_tree_dir = out_dir / "species_tree_input"
    tables_dir = out_dir / "tables"
    for d in [out_dir, gene_fasta_dir, species_tree_dir, tables_dir]:
        d.mkdir(parents=True, exist_ok=True)

    log_path = str(out_dir / "extract_step.log")
    logger = setup_logger(log_path)
    log_header(logger, args, start_time)

    orthofinder_dir = Path(args.orthofinder_dir)
    faa_dir = Path(args.faa_dir)

    for p, label in [(orthofinder_dir, "--orthofinder_dir"), (faa_dir, "--faa_dir")]:
        if not p.exists():
            logger.error(f"{label} path does not exist: {p}")
            log_footer(logger, start_time, "FAILURE")
            sys.exit(1)

    locus_to_seq, locus_to_species, stem_to_species = build_sequence_index(faa_dir, logger)
    if not locus_to_seq:
        logger.error("No sequences indexed — check --faa_dir")
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    gc_df = load_gene_count_matrix(orthofinder_dir, logger)
    og_txt = load_orthogroups_txt(orthofinder_dir, logger)

    col_to_species: dict[str, str] = {
        col: stem_to_species.get(col, sanitise_name(col))
        for col in gc_df.columns
    }
    all_species_ordered = [col_to_species[c] for c in sorted(gc_df.columns)]
    n_species = gc_df.shape[1]

    strict_list, relaxed_list = select_orthologs(gc_df, args.min_presence, logger)

    selected = strict_list if args.min_presence >= 1.0 else relaxed_list
    mode_label = "strict" if args.min_presence >= 1.0 else f"relaxed ({args.min_presence*100:.0f}%)"
    logger.info(f"Using {mode_label} mode: {len(selected)} orthologs selected")

    if not selected:
        logger.error(
            "No orthologs passed the selection filter. "
            "Try lowering --min_presence (e.g. 0.95) or check OrthoFinder results."
        )
        log_footer(logger, start_time, "FAILURE")
        sys.exit(1)

    og_seq_map = build_og_seq_map(
        selected, og_txt, gc_df, locus_to_seq, locus_to_species, stem_to_species, logger
    )

    n_selected = len(selected)
    logger.info(f"Writing {n_selected} per-ortholog FASTA files to {gene_fasta_dir} ...")

    inclusion_rows = []
    presence_matrix: dict[str, dict[str, str]] = {sp: {} for sp in all_species_ordered}
    gap_insertions = 0
    included_orthologs: list[str] = []

    if TQDM_AVAILABLE:
        og_iter = tqdm(selected, desc="Extracting orthologs", unit="OG",
                       bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]")
    else:
        og_iter = selected

    for i_og, og_id in enumerate(og_iter, 1):
        if not TQDM_AVAILABLE:
            print(f"  [{i_og:>{len(str(n_selected))}}/{n_selected}] {og_id}", flush=True)
        seqs = og_seq_map.get(og_id, {})
        gc_row = gc_df.loc[og_id] if og_id in gc_df.index else None

        present_lengths = [len(str(r.seq)) for r in seqs.values()]
        mean_len = (
            int(round(sum(present_lengths) / len(present_lengths)))
            if present_lengths else 0
        )

        records_out: list[SeqRecord] = []
        for species in all_species_ordered:
            if species in seqs:
                records_out.append(seqs[species])
                presence_matrix[species][og_id] = "1"
            else:
                is_multicopy = False
                if gc_row is not None:
                    for col, sp in col_to_species.items():
                        if sp == species and gc_row[col] > 1:
                            is_multicopy = True
                            break
                if is_multicopy:
                    presence_matrix[species][og_id] = "M"
                else:
                    gap_rec = SeqRecord(Seq("-" * mean_len), id=species, description="")
                    records_out.append(gap_rec)
                    presence_matrix[species][og_id] = "0"
                    logger.warning(
                        f"Gap inserted: species='{species}' orthogroup='{og_id}' "
                        f"(gap length={mean_len} aa)"
                    )
                    gap_insertions += 1

        SeqIO.write(records_out, str(gene_fasta_dir / f"{og_id}.faa"), "fasta")

        n_present = len(seqs)
        pct_present = 100.0 * n_present / n_species
        inclusion_rows.append({
            "orthogroup_id": og_id,
            "n_species_present": n_present,
            "pct_species_present": round(pct_present, 1),
            "mean_seq_length_aa": mean_len,
            "min_seq_length_aa": min(present_lengths) if present_lengths else 0,
            "max_seq_length_aa": max(present_lengths) if present_lengths else 0,
            "included": "yes",
        })
        included_orthologs.append(og_id)

    logger.info(
        f"Wrote {len(included_orthologs)} per-ortholog FASTA files"
    )
    if gap_insertions > 0:
        logger.warning(
            f"Total gap sequences inserted (relaxed mode missing species): {gap_insertions}"
        )

    logger.info(f"Writing concatenated archive FASTA to {species_tree_dir} ...")
    concat_records: list[SeqRecord] = []
    for species in all_species_ordered:
        parts = []
        for og_id in included_orthologs:
            og_path = gene_fasta_dir / f"{og_id}.faa"
            seq_found = ""
            for record in SeqIO.parse(str(og_path), "fasta"):
                if record.id == species:
                    seq_found = str(record.seq)
                    break
            parts.append(seq_found)
        concat_records.append(
            SeqRecord(Seq("".join(parts)), id=species, description="")
        )
    concat_path = species_tree_dir / "all_species_all_orthologs.faa"
    SeqIO.write(concat_records, str(concat_path), "fasta")
    logger.info(f"Archive FASTA written: {concat_path}")

    included_path = out_dir / "included_orthologs.txt"
    with open(included_path, "w") as f:
        for og_id in included_orthologs:
            f.write(og_id + "\n")

    logger.info(f"Writing summary tables to {tables_dir} ...")

    pd.DataFrame(inclusion_rows).to_csv(
        tables_dir / "orthologs_included.tsv", sep="\t", index=False
    )

    presence_rows = []
    for species in all_species_ordered:
        row = {"species": species}
        for og_id in included_orthologs:
            row[og_id] = presence_matrix[species].get(og_id, "0")
        presence_rows.append(row)
    pd.DataFrame(presence_rows).set_index("species").to_csv(
        tables_dir / "species_gene_presence.tsv", sep="\t"
    )

    total_aa = sum(len(str(r.seq)) for r in concat_records)
    pd.DataFrame([{
        "total_species": n_species,
        "total_orthogroups_assessed": len(gc_df),
        "orthogroups_passing_strict": len(strict_list),
        f"orthogroups_passing_relaxed_{int(args.min_presence*100)}pct": len(relaxed_list),
        "orthogroups_selected": len(included_orthologs),
        "total_aa_in_concat": total_aa,
    }]).to_csv(tables_dir / "ortholog_summary.tsv", sep="\t", index=False)

    logger.info("Summary tables written.")
    logger.info(
        f"FINAL TALLY: strict={len(strict_list)}, "
        f"relaxed@{args.min_presence*100:.0f}%={len(relaxed_list)}, "
        f"selected ({mode_label})={len(included_orthologs)}, "
        f"species={n_species}, "
        f"gap_insertions={gap_insertions}"
    )

    log_footer(logger, start_time, "SUCCESS")


if __name__ == "__main__":
    main()
