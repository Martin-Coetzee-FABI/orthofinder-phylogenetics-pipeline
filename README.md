# OrthoFinder3 Single-Copy Ortholog & Species Tree Pipeline

A modular, reproducible pipeline for identifying single-copy orthologs across
bacterial/archaeal proteomes and inferring a concatenated species tree.

---

## Overview

This pipeline takes annotated proteomes (`.faa` files, one per genome) and
produces a species phylogeny via:

1. **QC** — flag and exclude low-quality proteomes before analysis
2. **OrthoFinder3** — cluster proteins into orthogroups using DIAMOND
3. **Ortholog extraction** — select single-copy orthologs (strict or relaxed)
4. **Alignment & trimming** — MAFFT (align) → trimAl (trim)
5. **Species tree** — IQ-TREE3 reads trimmed alignments directly from the directory;
   ultrafast bootstrap + SH-aLRT branch support
6. **PDF report** — multi-page summary of QC results, ortholog selection, alignment
   statistics, functional gene categories, and the species tree

All sequences are **amino acid (protein)**. The pipeline runs on Mac or Linux.

---

## Pipeline diagram

```
proteomes/
  Strain_A.faa
  Strain_B.faa          Input: one .faa per genome
  ...
       │
       ▼
01_qc_proteomes.py      QC check → passed_proteomes/
       │
       ▼
02_run_orthofinder.py   OrthoFinder3 → orthofinder_results/Results_*/
       │
       ▼
03_extract_orthologs.py Select single-copy orthologs
       │                → gene_fastas/OG*.faa  (one per gene)
       │                → species_tree_input/   (archive concat)
       │
       ▼
04_align_trim_concat.py MAFFT → trimAl
       │                → aligned/  trimmed/  fasconcat_ready/
       │
       ▼
05_run_iqtree.py        IQ-TREE3 species tree
                        → iqtree_results/species_tree.treefile
       │
       ▼
06_generate_report.py   Multi-page PDF summary
                        → report/pipeline_summary.pdf
```

---

## Dependencies

### Install via conda (recommended)

```bash
conda env create -f environment.yml
conda activate orthofinder_pipeline
```

This installs Python 3.11, OrthoFinder ≥ 3.0, DIAMOND, MAFFT, trimAl,
IQ-TREE ≥ 3.0, BioPython, pandas, numpy, tqdm, matplotlib, and perl.

### FASconCAT-G (optional, manual install)

FASconCAT-G is **not required** by the pipeline. IQ-TREE3 reads the trimmed
alignment files directly from a directory without needing a pre-concatenated
supermatrix.

However, if you want a concatenated supermatrix file (e.g. for use with other
tools, or for manual inspection), Step 4 produces a `fasconcat_ready/` folder
containing all trimmed alignments renamed to `.fas` extension. Run FASconCAT-G
inside that folder:

```bash
# Install FASconCAT-G
wget https://github.com/PatrickKueck/FASconCAT-G/raw/master/FASconCAT-G_v1.05.pl
chmod +x FASconCAT-G_v1.05.pl

# Run it inside the fasconcat_ready folder
cd alignment_results/fasconcat_ready
perl /path/to/FASconCAT-G_v1.05.pl -s -p
# Produces: FcC_supermatrix.fas and FcC_supermatrix_partition.txt
```

### Tool versions tested

| Tool          | Version  |
|---------------|----------|
| Python        | 3.11     |
| OrthoFinder   | ≥ 3.0    |
| DIAMOND       | ≥ 2.1    |
| MAFFT         | ≥ 7.5    |
| trimAl        | ≥ 1.4    |
| IQ-TREE       | ≥ 3.0    |
| BioPython     | ≥ 1.81   |
| matplotlib    | ≥ 3.7    |
| FASconCAT-G   | 1.05     |

---

## Input format

One `.faa` file per genome in a single directory. Files must be named after
the strain/species (the filename stem becomes the species identifier):

```
proteomes/
  Escherichia_coli_K12.faa
  Salmonella_enterica_LT2.faa
  Klebsiella_pneumoniae_NTUH.faa
```

Headers within each `.faa` file can be standard Prokka or PGAP locus tags —
the pipeline renames all headers to the species name automatically.

---

## Test dataset

A small, anonymised test dataset is included in the `proteomes/` directory so
you can verify that the pipeline runs correctly before using your own data.

### What is included

| File | Role | Notes |
|------|------|-------|
| `Genome_out.faa` | **Outgroup** | Use as guide-tree root |
| `Genome_01.faa` | Ingroup | Divergent singleton species A |
| `Genome_02.faa` | Ingroup | Divergent singleton species B |
| `Genome_03.faa` | Ingroup | Divergent singleton species C |
| `Genome_04.faa` | Ingroup | Divergent singleton species D |
| `Genome_05.faa` | Ingroup | Divergent singleton species E |
| `Genome_06.faa` | Ingroup | Species complex — **clade 1** representative |
| `Genome_07.faa` | Ingroup | Species complex — **clade 2** representative |
| `Genome_08.faa` | Ingroup | Species complex — **clade 3** representative |

Nine genomes in total (~2.5–3 MB each, ~25 MB combined). The dataset is
designed to:

- Include an **outgroup** genome (`Genome_out`) that is phylogenetically
  distant from the ingroup — a realistic scenario for rooting.
- Represent **three distinct clades** within one species complex
  (`Genome_06`–`08`), testing the pipeline's ability to resolve closely
  related genomes.
- Include several **divergent singleton species** (`Genome_01`–`05`),
  exercising the QC and ortholog-selection steps across a range of
  evolutionary distances.
- Be small enough to **complete in 15–30 minutes** on a modern laptop or
  desktop.

The genome sequences are real annotated bacterial proteomes (PROKKA format).
Only the file names have been anonymised; all internal locus tags and protein
sequences are unchanged.

### Guide tree

A rooted Newick guide tree for OrthoFinder is provided as
`test_guide_tree.tre`. It was derived from the published phylogeny and
pruned to the nine test genomes. Branch lengths are preserved; internal
node support values have been removed.

```
(((((Genome_01:0.012,Genome_02:0.013):0.013,
    ((Genome_06:0.008,Genome_07:0.006):0.002,
      Genome_08:0.008):0.013):0.006,
    (Genome_03:0.022,Genome_04:0.025):0.009):0.036,
  Genome_05:0.059):0.263,
Genome_out:0.000);
```

`Genome_out` sits on the long outgroup branch (0.263 substitutions/site),
confirming the correct rooting.

### Running the test

#### Option A — GUI

```bash
conda activate orthofinder_pipeline
streamlit run app.py
```

In the sidebar: the `proteomes/` directory and `test_guide_tree.tre` are
pre-filled. Set **Threads** to 4–8, then click **▶ Run Pipeline**.

#### Option B — Bash script

```bash
conda activate orthofinder_pipeline
cd /path/to/pipeline

bash pipeline_runner.sh
```

When the settings screen appears, verify:

| Item | Value |
|------|-------|
| **[1] Input dir** | `proteomes` |
| **[8] Species tree** | `test_guide_tree.tre` |
| **[12] Model** | `LG+F+G4` |
| **[4] Threads** | 4–8 |

Press **Enter** to start. All six steps should complete without errors.

### Expected test output

```
qc_results/
  proteome_qc_summary.tsv    — all 9 genomes should PASS
  passed_proteomes/          — 9 symlinks

ortholog_results/
  gene_fastas/               — ~1 000–2 000 single-copy ortholog FASTAs
  tables/ortholog_summary.tsv

alignment_results/
  trimmed/                   — aligned + trimmed orthologs

iqtree_results/
  species_tree.treefile      — 9-taxon tree; Genome_out should be outgroup

report/
  pipeline_summary.pdf       — full summary with flow diagram and references
```

### Replacing the test data with your own

Once the test passes, simply replace the contents of `proteomes/` with your
own `.faa` files (one per genome, named after the strain/species). Delete or
update `test_guide_tree.tre` if you want to use a different guide tree, or
leave `OF_SPECIES_TREE` blank to let OrthoFinder infer the species tree.

---

## Two ways to run the pipeline

The pipeline can be run in two ways. Both use the same underlying Python
scripts and produce identical output. Choose whichever suits your environment.

---

> **Windows is not supported.** OrthoFinder, DIAMOND, MAFFT, trimAl, and
> IQ-TREE are Linux/macOS tools. Windows users should use
> [WSL 2](https://learn.microsoft.com/en-us/windows/wsl/install) and run
> the pipeline inside a WSL Ubuntu terminal.

### Option A — Browser GUI (`app.py`)

Best for: **local workstations** (Mac or Linux desktop) where a browser is
available. Suitable for users who are not comfortable with the command line.

> ⚠️ **Activate your conda environment *before* launching the GUI.**
> If you launch `streamlit run app.py` from a shell where the environment is
> not active, all tools will appear as missing in the pre-flight check even
> though they are correctly installed.

```bash
# 1. Activate the environment that contains orthofinder, diamond, mafft, etc.
conda activate orthofinder3     # replace with your actual environment name

# 2. One-time: install GUI dependencies into that same environment
conda install -c conda-forge "streamlit>=1.37" psutil

# 3. Launch from the pipeline directory
streamlit run app.py
```

Streamlit opens the GUI automatically at `http://localhost:8501`.

**What the GUI provides:**
- Point-and-click settings panel (all 15+ parameters, with descriptions and sensible defaults)
- **Step selection** — choose which steps to run (e.g. start from Step 4 if alignment already exists)
- Pre-flight tool check before running, with hints for missing tools
- ▶ Run / ⏹ Stop controls
- Live log that auto-refreshes every 3 seconds while the pipeline runs
- Step progress indicators (⬜ pending → 🔄 running → ✅ done → ❌ failed / ⏭️ skipped)
- Results tab: QC table, ortholog stats, alignment stats, rendered species tree
- Download buttons for the species tree (.treefile / .nwk) and PDF report

The pipeline process is **detached** from the browser — closing the tab does
not stop the pipeline. Re-opening `app.py` reconnects to the running job by
reading the log and marker files.

> **Note:** The GUI will not work on a headless Linux server accessed over SSH
> without X forwarding or a browser. Use Option B in that case.

---

### Option B — Bash wrapper (`pipeline_runner.sh`)

Best for: **Linux servers, HPC clusters, SSH sessions**, or any situation
where a browser is unavailable or undesirable. Also suitable for scripted /
automated runs.

> ⚠️ **Activate your conda environment *before* running the script.**
> If the environment is not active, the tool checks will fail immediately.

```bash
# 1. Activate the environment first
conda activate orthofinder3     # replace with your actual environment name

# 2. Interactive mode — confirms settings and prompts between each step
bash pipeline_runner.sh

# 3. Automatic mode — runs all steps without prompts (good for nohup / screen)
bash pipeline_runner.sh --auto
```

Edit the `USER CONFIGURATION` block at the top of `pipeline_runner.sh` to
set all parameters before running. The script prints a full settings summary
and prompts for the number of CPU threads before any analysis begins.

---

### Shared behaviour (both options)

Both options call the same six Python scripts in the same order, write to the
same output directories, and respect the same `.step_complete` marker files.
A run started with one option can be inspected or continued with the other.

---

## Quick start

### First run — use the included test dataset

The repository ships with 9 anonymised test proteomes in `proteomes/` and a
matching rooted guide tree (`test_guide_tree.tre`).  Running the pipeline on
this dataset confirms your installation is working before you commit to a
full analysis.

```bash
# 1. Install the conda environment (one-time)
conda env create -f environment.yml
conda activate orthofinder_pipeline

# 2. Run the test — the proteomes/ folder is pre-populated
bash pipeline_runner.sh
# When the settings screen opens:
#   [8]  Species tree → test_guide_tree.tre   (auto-detected)
#   [12] Model        → LG+F+G4
#   [4]  Threads      → 4–8
# Press Enter to start.  Total runtime: ~15–30 min on a modern laptop.
```

A successful test produces `iqtree_results/species_tree.treefile` (9 taxa,
`Genome_out` as the outgroup) and `report/pipeline_summary.pdf`.

### Your own data

```bash
# 1. Activate the conda environment
conda activate orthofinder_pipeline   # or your actual environment name

# 2. Place your .faa files in the proteomes/ directory (replace test files)

# ── Option A: Browser GUI (local Mac/Linux workstation) ───────────────────
streamlit run app.py
# Opens http://localhost:8501 — configure settings in the sidebar, then click ▶ Run

# ── Option B: Bash script (Linux server / SSH / automated) ────────────────
bash pipeline_runner.sh           # interactive — review settings, confirm each step
bash pipeline_runner.sh --auto    # automatic  — runs all steps without prompts

# 3. View the species tree
# Open iqtree_results/species_tree.treefile in FigTree or upload to iTOL

# 4. View the summary report
open report/pipeline_summary.pdf      # macOS
xdg-open report/pipeline_summary.pdf  # Linux
```

---

## Running the pipeline — `pipeline_runner.sh`

The master script runs any combination of the six pipeline steps, logs
everything to `pipeline_run.log`, and presents a single settings screen
before the first step executes.

### Settings screen

Every run (interactive or automatic) opens a numbered settings table where
you can review and edit every parameter before committing:

```
  ╔══════════════════════════════════════════════════════╗
  ║  Pipeline Settings                                   ║
  ╚══════════════════════════════════════════════════════╝

   General
   ─────────────────────────────────────────────────────
    [ 1]  Input proteomes dir   :  proteomes
    [ 2]  Min proteins (QC)     :  500
    [ 3]  Min ortholog presence :  1.0
    [ 4]  Threads               :  8  (176 available)

   OrthoFinder
   ─────────────────────────────────────────────────────
    [ 5]  Search program  (-S)  :  diamond
    [ 6]  Gene tree       (-M)  :  msa
    [ 7]  MSA program     (-A)  :  mafft
    [ 8]  Species tree    (-s)  :  mlsa_prokka.tre
    [ 9]  Timeout               :  86400 s

   Alignment & Trimming
   ─────────────────────────────────────────────────────
    [10]  trimAl mode           :  auto

   IQ-TREE
   ─────────────────────────────────────────────────────
    [11]  Mode                  :  partition
    [12]  Model                 :  LG+F+G4
    [13]  Bootstrap reps        :  1000
    [14]  SH-aLRT reps          :  1000

   Steps to run  (PDF report always runs at the end)
   ─────────────────────────────────────────────────────
    [16]  Start from step      :  1  (QC proteomes)
    [17]  Stop after step      :  5  (IQ-TREE species tree)

   Pipeline
   ─────────────────────────────────────────────────────
    [15]  Run mode              :  interactive

   ─────────────────────────────────────────────────────
   Enter number to edit a setting, or [Enter] to start.
   ─────────────────────────────────────────────────────

   >
```

Type any number to edit that setting. The table redraws immediately
after each change so you can see the effect. Press **Enter** on an empty
line to confirm all settings and start the pipeline.

### Step selection — items [16] and [17]

The most powerful feature of the settings screen: you can run **any
contiguous range of steps** rather than always starting from the
beginning. When you type `16` or `17`, the pipeline shows a brief
description of every step before asking for your selection:

```
    [1]  QC proteomes
         Checks every .faa file for protein count, length, and duplicates.
         Flags WARN genomes (one issue) and excludes FAIL genomes.
         Needs  : proteomes/ directory with .faa files
         Outputs: qc_results/passed_proteomes/

    [2]  OrthoFinder
         Runs DIAMOND all-vs-all search and clusters proteins into orthogroups.
         Needs  : qc_results/passed_proteomes/
         Outputs: orthofinder_results/Results_*/

    [3]  Extract single-copy orthologs
         Selects orthogroups present in all (or most) species with 1 copy each.
         Writes one .faa file per ortholog with species names as headers.
         Needs  : orthofinder_results/Results_*/
         Outputs: ortholog_results/gene_fastas/

    [4]  Align & trim (MAFFT + trimAl)
         Aligns each ortholog with MAFFT and removes gappy columns with trimAl.
         Needs  : ortholog_results/gene_fastas/
         Outputs: alignment_results/trimmed/

    [5]  IQ-TREE species tree
         Infers the species phylogeny from the trimmed alignments.
         Needs  : alignment_results/trimmed/
         Outputs: iqtree_results/species_tree.treefile

    [6]  PDF report  (always runs at the end)
         Summarises QC, orthologs, alignment stats, and the species tree.
         Uses whatever output directories are already present.
         Outputs: report/pipeline_summary.pdf
```

**Common partial-run scenarios:**

| Situation | Start | Stop |
|-----------|-------|------|
| Full pipeline (default) | 1 | 5 |
| OrthoFinder already done; redo alignment + tree | 4 | 5 |
| Only want orthogroup extraction, not the tree | 1 | 3 |
| Trimmed alignments exist; redo IQ-TREE only | 5 | 5 |
| Full pipeline minus the final tree | 1 | 4 |

The PDF report (Step 6) **always runs** at the end of whatever range is
selected, using whichever result directories are present.

If the required input directory for a step is missing (e.g. you start at
Step 4 but `ortholog_results/gene_fastas/` does not exist), the script
exits immediately with a clear error message telling you which earlier
step to run first.

### Run mode — item [15]

| Mode | Behaviour |
|------|-----------|
| **interactive** (default) | Pauses between steps with `[Enter] Continue / [q] Quit` prompt. Lets you inspect each step's output before proceeding. |
| **automatic** | Runs all selected steps straight through without stopping. Ideal for overnight jobs or server runs. |

You can also force automatic mode from the command line:

```bash
bash pipeline_runner.sh --auto
```

Non-interactive sessions (SSH without a terminal, cron, nohup) detect
the missing terminal automatically and switch to automatic mode.

### Resuming a stopped run

Each step writes a `.step_complete` marker file when it finishes
(`qc_results/.step_complete`, `ortholog_results/.step_complete`, etc.).
If you re-run the pipeline with the same settings, steps that have
already completed are skipped with a `[SKIP]` message. You never have to
re-run a long OrthoFinder search just to redo the IQ-TREE step.

To force a step to re-run, delete its marker file:
```bash
rm iqtree_results/.step_complete   # forces Step 5 to re-run
```

### Configuration defaults

Edit the **USER CONFIGURATION** block at the top of `pipeline_runner.sh`
to change the defaults that appear in the settings screen:

**General**

| Variable        | Default     | Description                                                         |
|-----------------|-------------|---------------------------------------------------------------------|
| `INPUT_FAA_DIR` | `proteomes` | Directory containing `.faa` files                                   |
| `MIN_PROTEINS`  | `500`       | QC minimum protein count per genome                                 |
| `MIN_PRESENCE`  | `1.0`       | Ortholog selection (1.0 = strict, 0.95 = relaxed)                   |
| `THREADS`       | `8`         | CPU threads; always shown in the settings screen before the run     |
| `STEP_START`    | `1`         | Default first step (1–5)                                            |
| `STEP_END`      | `5`         | Default last step (1–5); PDF report always runs after               |

**OrthoFinder**

| Variable          | Default       | Description                                                                    |
|-------------------|---------------|--------------------------------------------------------------------------------|
| `OF_SEARCH`       | `diamond`     | Search program (`-S`): `diamond`, `diamond_ultra_sens`, `blast`, `mmseqs`      |
| `OF_GENE_TREE`    | `msa`         | Gene tree method (`-M`): `msa` (recommended) or `dendroblast` (faster)        |
| `OF_MSA`          | `mafft`       | MSA program (`-A`): `mafft` / `famsa` / `muscle`. Only used with `msa` method |
| `OF_SPECIES_TREE` | auto-detected | Rooted Newick tree (`-s`). Auto-detects first `.tre`/`.tree` in working dir.  |
| `OF_TIMEOUT`      | `86400`       | Kill OrthoFinder after this many seconds (0 = no limit)                        |

**Alignment & Trimming**

| Variable      | Default | Description                                               |
|---------------|---------|-----------------------------------------------------------|
| `TRIMAL_MODE` | `auto`  | trimAl strategy: `auto` / `gappyout` / `strict` / `strictplus` |

**IQ-TREE**

| Variable       | Default    | Description                                                                          |
|----------------|------------|--------------------------------------------------------------------------------------|
| `IQTREE_MODE`  | `partition`| `partition` (-p, per-gene model) or `single` (-s, one model for all genes)          |
| `IQTREE_MODEL` | `LG+F+G4`  | Substitution model for **both** modes. `LG+F+G4` = fast; `MFP` = ModelFinder per gene (slow); `MFP+MERGE` = ModelFinder + merge partitions (slowest) |
| `BOOTSTRAP`    | `1000`     | Ultrafast bootstrap replicates (`-B`)                                                |
| `ALRT`         | `1000`     | SH-aLRT replicates (`-alrt`); set to `0` to disable                                 |

---

## Step-by-step guide

### Step 1 — `01_qc_proteomes.py`

Checks every `.faa` file and assigns PASS / WARN / FAIL:

| Metric                 | Threshold      | Verdict         |
|------------------------|----------------|-----------------|
| Protein count          | < 500          | Flag            |
| Median protein length  | < 100 aa       | Flag            |
| % proteins < 50 aa     | > 20%          | Flag            |
| Duplicate sequence IDs | any            | FAIL (always)   |
| Two or more flags      | —              | FAIL            |
| Exactly one flag       | —              | WARN (included) |
| No flags               | —              | PASS            |

WARN genomes are **included** in the analysis (only one marginal metric).
FAIL genomes are **excluded** and listed in `failed_genomes.txt`.

A per-genome progress counter is printed to the screen as each file is
checked: `[ 3/45] Checking: Strain_C`.

```bash
python 01_qc_proteomes.py \
  --input_dir proteomes \
  --min_proteins 500 \
  --min_median_len 100 \
  --max_short_pct 20.0 \
  --output_dir qc_results
```

**Parameters:**

| Parameter        | Default         | Description                                      |
|------------------|-----------------|--------------------------------------------------|
| `--input_dir`    | required        | Directory containing .faa files                  |
| `--min_proteins` | 500             | Minimum protein count threshold                  |
| `--min_median_len` | 100           | Minimum median protein length (aa)               |
| `--max_short_pct` | 20.0          | Maximum % proteins < 50 aa                       |
| `--exclude`      | ""              | Comma-separated genome names to force-exclude    |
| `--output_dir`   | qc_results      | Output directory                                 |
| `--copy`         | False (symlink) | Copy files instead of symlinking                 |

---

### Step 2 — `02_run_orthofinder.py`

Runs OrthoFinder3 on the passed proteomes. Auto-detects the OrthoFinder
executable from common install locations. Before launching, the script
verifies that the configured sequence search tool and (if using `msa` gene
tree inference) the MSA tool are present in PATH — it exits immediately with
a clear error message if either is missing.

```bash
python 02_run_orthofinder.py \
  --input_dir qc_results/passed_proteomes \
  --output_dir orthofinder_results \
  --threads 8 \
  --search_prog diamond \
  --gene_tree_method msa \
  --msa_prog mafft \
  --species_tree mlsa_prokka.tre \
  --timeout 86400
```

**Parameters:**

| Parameter              | Default                         | Description                                                                  |
|------------------------|---------------------------------|------------------------------------------------------------------------------|
| `--input_dir`          | qc_results/passed_proteomes     | Directory of .faa files                                                      |
| `--output_dir`         | orthofinder_results             | Where OrthoFinder writes results                                             |
| `--threads`            | all available                   | CPU threads; capped to `os.cpu_count()` automatically                        |
| `--search_prog`        | `diamond`                       | Sequence search program (`-S`): `diamond`, `diamond_ultra_sens`, `blast`, `mmseqs` |
| `--gene_tree_method`   | `msa`                           | Gene tree method (`-M`): `msa` or `dendroblast`                              |
| `--msa_prog`           | `mafft`                         | MSA program (`-A`): `mafft`, `famsa`, `muscle`. Only used with `--gene_tree_method msa` |
| `--species_tree`       | ""                              | Path to a user-supplied rooted Newick species tree (`-s`)                    |
| `--timeout`            | `86400`                         | Kill OrthoFinder after this many seconds; `0` = no timeout                   |
| `--orthofinder`        | auto-detected                   | Path to OrthoFinder executable                                               |
| `--extra_args`         | ""                              | Additional OrthoFinder arguments as a quoted string                          |

If a completed `Results_*` directory already exists inside `--output_dir`,
the step is skipped and the existing results are used. An incomplete
`Results_*` directory (e.g. from a killed run) is removed automatically
before OrthoFinder is re-launched.

---

### Step 3 — `03_extract_orthologs.py`

Reads the OrthoFinder gene count matrix and selects single-copy orthologs.
Writes per-gene FASTAs, an archive concatenation, and summary tables.

```bash
python 03_extract_orthologs.py \
  --orthofinder_dir orthofinder_results/Results_Jun01 \
  --faa_dir qc_results/passed_proteomes \
  --min_presence 1.0 \
  --output_dir ortholog_results
```

**Parameters:**

| Parameter           | Default          | Description                                     |
|---------------------|------------------|-------------------------------------------------|
| `--orthofinder_dir` | required         | Path to OrthoFinder `Results_*/` folder         |
| `--faa_dir`         | required         | Directory with original .faa files              |
| `--min_presence`    | 1.0              | Min fraction of species with the gene           |
| `--output_dir`      | ortholog_results | Output directory                                |

See **Ortholog selection** section for guidance on strict vs relaxed mode.

Progress bars (via `tqdm`) are shown for genome indexing and ortholog
extraction. If `tqdm` is not installed, a plain `[X/N]` counter is printed
instead.

---

### Step 4 — `04_align_trim_concat.py`

Aligns each ortholog with MAFFT and trims with trimAl. The trimmed alignment
files in `trimmed/` are the direct input to IQ-TREE3 in Step 5 — no
concatenated supermatrix is needed.

A `fasconcat_ready/` folder is also produced, containing the same trimmed
alignments renamed to `.fas` extension, ready for FASconCAT-G if you want
to produce a supermatrix manually (see **FASconCAT-G** under Dependencies).

```bash
python 04_align_trim_concat.py \
  --gene_fasta_dir ortholog_results/gene_fastas \
  --output_dir alignment_results \
  --threads 8 \
  --trimal_mode auto
```

**Parameters:**

| Parameter          | Default                      | Description                             |
|--------------------|------------------------------|-----------------------------------------|
| `--gene_fasta_dir` | ortholog_results/gene_fastas | Per-ortholog .faa files from Step 3     |
| `--output_dir`     | alignment_results            | Output directory                        |
| `--threads`        | all available                | Threads for MAFFT                       |
| `--mafft`          | auto-detected                | Path to MAFFT executable                |
| `--trimal`         | auto-detected                | Path to trimAl executable               |
| `--trimal_mode`    | auto                         | auto / gappyout / strict / strictplus   |

**trimAl modes:**

| Mode         | Behaviour                                                       |
|--------------|-----------------------------------------------------------------|
| `auto`       | `-automated1`: smart trimming, avoids over-trimming (default)   |
| `gappyout`   | Remove columns with high gap frequency                          |
| `strict`     | More aggressive gap-based trimming                              |
| `strictplus` | Even more aggressive                                            |

A `tqdm` progress bar tracks alignment and trimming across all orthologs
(`[42/1250]`). Genes that lose >50% of columns after trimming are flagged
as warnings above the progress bar.

---

### Step 5 — `05_run_iqtree.py`

Runs IQ-TREE3 directly on the directory of trimmed alignment files from
Step 4. No concatenated supermatrix file is needed — IQ-TREE3 reads and
processes all alignment files in the directory itself.

Branch support is always computed with **both** ultrafast bootstrap (`-B`)
and SH-aLRT (`-alrt`). The `.treefile` will have two support values per
node, e.g. `95/98`, where the first is SH-aLRT and the second is UFBoot.

IQ-TREE output is streamed live to the screen during the run, framed by a
clear banner:

```
==============================================================
  IQ-TREE3 output (live):
==============================================================
IQ-TREE multicore version 3.x ...
Scanning alignment files in directory alignment_results/trimmed/ ...
...
==============================================================
```

To ensure output appears line-by-line rather than in delayed bursts,
the script wraps the IQ-TREE command with `stdbuf -oL` (Linux / macOS
with GNU coreutils) or `unbuffer` (from the `expect` package) when either
is available. The method used is recorded in `iqtree_step.log`. If neither
tool is found, output still appears but may be delayed until IQ-TREE's
internal buffer fills.

```bash
python 05_run_iqtree.py \
  --trimmed_dir alignment_results/trimmed \
  --output_dir iqtree_results \
  --mode partition \
  --bootstrap 1000 \
  --alrt 1000 \
  --threads 8
```

**Parameters:**

| Parameter        | Default                        | Description                                    |
|------------------|--------------------------------|------------------------------------------------|
| `--trimmed_dir`  | alignment_results/trimmed      | Directory of trimmed alignments from Step 4    |
| `--output_dir`   | iqtree_results                 | Output directory                               |
| `--mode`         | partition                      | `partition`: `-p dir/` (per-gene ModelFinder, recommended); `single`: `-s dir/` (one model) |
| `--model`        | LG+F+G4                        | Substitution model passed to IQ-TREE `-m`. Applies to **both** modes. `LG+F+G4` = fast fixed model (recommended); `MFP` = ModelFinder per gene (slow); `MFP+MERGE` = ModelFinder + merge partitions (slowest) |
| `--bootstrap`    | 1000                           | Ultrafast bootstrap replicates (`-B`); set via `BOOTSTRAP` in `pipeline_runner.sh` |
| `--alrt`         | 1000                           | SH-aLRT replicates (`-alrt`); set independently via `ALRT` in `pipeline_runner.sh`. Set to `0` to disable. |
| `--threads`      | all available                  | CPU threads                                    |
| `--iqtree`       | auto-detected                  | Path to IQ-TREE3 executable                    |
| `--prefix`       | species_tree                   | Output file prefix                             |

---

### Step 6 — `06_generate_report.py`

Generates a multi-page PDF summary of the entire pipeline run and writes it to
`report/pipeline_summary.pdf`. The report is self-contained — it collects data
from all previous output directories automatically.

```bash
python 06_generate_report.py \
  --output_dir report
```

**Parameters:**

| Parameter             | Default                           | Description                                                |
|-----------------------|-----------------------------------|------------------------------------------------------------|
| `--output_dir`        | report                            | Directory where the PDF is written                         |
| `--qc_dir`            | qc_results                        | Step 1 output directory                                    |
| `--ortholog_dir`      | ortholog_results                  | Step 3 output directory                                    |
| `--alignment_dir`     | alignment_results                 | Step 4 output directory                                    |
| `--iqtree_dir`        | iqtree_results                    | Step 5 output directory                                    |
| `--orthofinder_results` | (auto-detected)                 | Path to OrthoFinder `Results_*/` folder                    |
| `--faa_dir`           | qc_results/passed_proteomes       | Directory of passed .faa files (used for gene annotations) |
| `--prefix`            | species_tree                      | IQ-TREE output file prefix                                 |
| `--output_file`       | pipeline_summary.pdf              | PDF filename within `--output_dir`                         |

**PDF pages:**

| Page | Content                                                                             |
|------|-------------------------------------------------------------------------------------|
| 1    | **At-a-glance stats** — genomes analysed/passed/failed, orthologs selected, total alignment length before and after trimming |
| 2    | **QC results** — verdict bar chart (PASS/WARN/FAIL counts) + full genome table with colour-coded verdicts |
| 3    | **Ortholog selection** — stats box, species completeness histogram, table of included orthologs with functional annotations |
| 4    | **Functional categories** — horizontal bar chart grouping orthologs into broad categories (Ribosome & Translation, DNA Replication, Transcription, etc.) with a note that categories are indicative keyword matches |
| 5    | **Alignment statistics** — before/after trimming bar chart, per-gene column retention histogram, table of the 10 most aggressively trimmed genes |
| 6    | **Species tree** — rendered tree (Bio.Phylo) with branch support values, plus the raw Newick string for copy-paste |

All pages degrade gracefully: if an output directory or file from a previous step
is missing, that page is either skipped or shows a "data not found" message, so
the report can be generated even from a partially completed run.

**Requirements:** `matplotlib` (included in `environment.yml`), `biopython`, `pandas`.

---

## Output file reference

```
qc_results/
  proteome_qc_summary.tsv       One row per genome: all QC metrics + verdict
  proteome_qc_report.txt        Human-readable QC table
  passed_proteomes/             Symlinks to PASS/WARN genomes (input to Step 2)
  failed_genomes.txt            Names of FAIL genomes with reasons
  qc_step.log                   Timestamped QC decisions and thresholds

orthofinder_results/
  Results_[date]/               Raw OrthoFinder output (untouched)
  orthofinder_run.log           Full OrthoFinder stdout/stderr
  orthofinder_stdout.log        Live OrthoFinder output stream
  orthofinder_command.txt       Exact command + version string

ortholog_results/
  gene_fastas/OG*.faa           One FASTA per ortholog; headers = species name
  species_tree_input/
    all_species_all_orthologs.faa  Archive: all orthologs concatenated per species
  included_orthologs.txt        Plain list of selected orthogroup IDs
  extract_step.log              Filtering decisions, gap insertions, final tally
  tables/
    orthologs_included.tsv      Per-ortholog stats (length, % presence)
    species_gene_presence.tsv   Species × orthogroup matrix (1/0/M)
    ortholog_summary.tsv        Run summary: counts, total aa

alignment_results/
  aligned/OG*.aln.faa           MAFFT output (unfiltered alignments)
  trimmed/OG*.trimmed.faa       trimAl output → direct input to IQ-TREE3
  fasconcat_ready/OG*.fas       Same trimmed files renamed .fas → run FASconCAT-G here
  alignment_stats.tsv           Per-gene: raw length, trimmed length, % retained
  commands.log                  Every MAFFT/trimAl command with timestamp + exit code
  align_trim_step.log           Tool versions, per-gene status, trimming warnings

iqtree_results/
  species_tree.treefile         Final species tree (Newick)
  species_tree.iqtree           Full IQ-TREE report (models, likelihood, etc.)
  species_tree.log              IQ-TREE native run log
  species_tree.contree          Consensus tree with bootstrap values on nodes
  iqtree_command.txt            Exact IQ-TREE command used
  iqtree_step.log               Version, mode, models selected, likelihood score

report/
  pipeline_summary.pdf          Multi-page PDF summary report
  report_step.log               PDF generation log (data sources, warnings)

pipeline_run.log                Master log: all steps, versions, timestamps
```

---

## Log files guide

| Log file                      | When to consult                                               |
|-------------------------------|---------------------------------------------------------------|
| `qc_results/qc_step.log`      | Any genome flagged; understanding why a genome was excluded   |
| `orthofinder_results/orthofinder_run.log` | OrthoFinder fails or produces unexpected results |
| `ortholog_results/extract_step.log` | Zero or very few orthologs; gap insertion warnings    |
| `alignment_results/commands.log` | MAFFT or trimAl fails on a specific gene               |
| `alignment_results/align_trim_step.log` | Genes losing >50% of columns after trimming       |
| `iqtree_results/iqtree_step.log` | IQ-TREE model selection, likelihood, convergence        |
| `report/report_step.log`      | PDF generation; which data files were found/missing           |
| `pipeline_run.log`            | Full run overview; timestamps for each step                   |

**Warnings to watch for:**

- `align_trim_step.log`: "trimming removed >50% of columns" — the alignment for
  that gene may be poor. Consider inspecting it manually.
- `extract_step.log`: "Gap inserted" — a species is missing from an ortholog in
  relaxed mode; a gap sequence has been inserted.
- `qc_step.log`: WARN genomes — included but flag one QC metric. Monitor their
  effect on orthogroup completeness.

---

## Ortholog selection

### Strict mode (`--min_presence 1.0`, default)

Every selected orthogroup must be present in **all N species** with **exactly 1
copy per species**. This produces the most reliable phylogenetic signal but may
drastically reduce the number of usable orthologs if some genomes are
incomplete or distant.

**Use when:**
- All genomes are high-quality complete assemblies
- You want maximum phylogenetic signal per site
- You are comfortable with a smaller but clean gene set

### Relaxed mode (e.g. `--min_presence 0.95`)

Orthogroups present in ≥ 95% of species, with exactly 1 copy where present.
Missing species receive a gap sequence (`-` × mean gene length) so the
concatenated alignment stays rectangular. Species with multi-copy genes in an
orthogroup are always excluded from that orthogroup.

**Use when:**
- Some genomes are incomplete (draft assemblies)
- A strict filter leaves fewer than ~50 orthologs
- You have a large dataset where 100% universality is unrealistic

**Practical guidance:**
- Start with strict mode. If you get < 50 orthologs, try 0.99, then 0.95.
- Check `tables/species_gene_presence.tsv` to identify which species cause
  the most missing data — these may be distant relatives or poor assemblies.

---

## IQ-TREE model choice

### How IQ-TREE3 reads the alignment directory

IQ-TREE3 can read a directory of alignment files directly without requiring a
pre-concatenated supermatrix. The behaviour depends on the flag used:

| Flag      | Behaviour                                                          |
|-----------|--------------------------------------------------------------------|
| `-p dir/` | Partition analysis — each file is one partition, best model per gene |
| `-s dir/` | Supermatrix analysis — all files concatenated, one model for all   |

The pipeline passes `alignment_results/trimmed/` to IQ-TREE3.

### Partition mode (recommended, `--mode partition`)

IQ-TREE runs ModelFinder on each gene partition and selects the best-fitting
substitution model per gene (`-m MFP+MERGE`). Models may be merged if they
are statistically equivalent (reducing over-parameterisation). This accounts
for rate variation across genes and is the standard approach for multi-gene
concatenated analyses.

### Single model mode (`--mode single`)

One model (`LG+F+G4` by default) is applied to all genes treated as one
concatenated alignment.

`LG+F+G4` explained:
- **LG**: Le & Gascuel (2008) general amino acid substitution matrix
- **+F**: empirical amino acid frequencies estimated from the data
- **+G4**: Gamma distribution with 4 rate categories (accounts for among-site rate variation)

**Use single mode when:**
- You want a fast run to check tree topology before a full partition run
- You are analysing closely related genomes with very similar gene compositions

### Branch support: SH-aLRT + ultrafast bootstrap

Both metrics are always computed:

| Metric              | Flag     | Interpretation                                      |
|---------------------|----------|-----------------------------------------------------|
| SH-aLRT             | `-alrt`  | ≥ 80% considered significant                        |
| Ultrafast bootstrap | `-B`     | ≥ 95% considered significant                        |

In FigTree, node labels will show two values separated by `/` (SH-aLRT/UFBoot).
A node is considered well-supported if **both** values exceed their thresholds.
Using both metrics reduces the risk of falsely high bootstrap values that can
occur with ultrafast bootstrap alone on certain topologies.

---

## Downstream use — individual gene trees

The `ortholog_results/gene_fastas/` directory contains one FASTA per
ortholog, with species names as headers. These are suitable for individual
gene tree inference, which can be used to:

- Detect horizontal gene transfer (genes with highly discordant topologies)
- Perform coalescent-based species tree inference (e.g. ASTRAL)
- Assess phylogenetic signal per gene

Example — build a gene tree for OG0000001:
```bash
mafft --auto OG0000001.faa > OG0000001.aln.faa
trimal -in OG0000001.aln.faa -out OG0000001.trimmed.faa -automated1
iqtree3 -s OG0000001.trimmed.faa -m MFP -B 1000 -T 4
```

---

## Re-running after a partial failure

Each pipeline step writes a small marker file when it completes successfully.
If the pipeline is interrupted (e.g. OrthoFinder killed, server rebooted,
session timeout), simply re-run `bash pipeline_runner.sh` — completed steps
are detected automatically and skipped, so only the failed step and everything
after it will re-run.

| Step | Marker file                          |
|------|--------------------------------------|
| 1    | `qc_results/.step_complete`          |
| 2    | `orthofinder_results/Results_*/Orthogroups/Orthogroups.GeneCount.tsv` (OrthoFinder's own output) |
| 3    | `ortholog_results/.step_complete`    |
| 4    | `alignment_results/.step_complete`   |
| 5    | `iqtree_results/.step_complete`      |

To **force a step to re-run**, delete its marker:

```bash
# Re-run ortholog extraction onwards
rm ortholog_results/.step_complete

# Re-run OrthoFinder (also removes the incomplete Results_* directory automatically)
rm -rf orthofinder_results/
```

Step 2 (OrthoFinder) never re-runs if a completed `Results_*` directory
already exists inside `orthofinder_results/`. To force a full re-run of
OrthoFinder, remove that directory:

```bash
rm -rf orthofinder_results/Results_*/
```

---

## Troubleshooting

### Zero orthologs found (strict mode)

**Symptom:** `extract_step.log` reports 0 orthologs passing strict filter.

**Causes and fixes:**
1. Some genomes are incomplete — check `qc_results/proteome_qc_summary.tsv`
2. Some taxa are very distantly related — try `--min_presence 0.95`
3. OrthoFinder run was incomplete — check `orthofinder_run.log`

### Running FASconCAT-G to produce a supermatrix

FASconCAT-G is not required by the pipeline. IQ-TREE3 reads the trimmed
alignments directly. However, if you need a concatenated supermatrix (e.g.
for use with RAxML or other tools), Step 4 has already prepared the
`alignment_results/fasconcat_ready/` folder for you:

```bash
cd alignment_results/fasconcat_ready
perl /path/to/FASconCAT-G_v1.05.pl -s -p
# Produces FcC_supermatrix.fas and FcC_supermatrix_partition.txt
```

To install FASconCAT-G:
```bash
wget https://github.com/PatrickKueck/FASconCAT-G/raw/master/FASconCAT-G_v1.05.pl
chmod +x FASconCAT-G_v1.05.pl
```

### IQ-TREE runs out of memory (partition mode)

**Symptom:** IQ-TREE crashes with a memory error on large datasets.

**Fixes:**
- Reduce threads: `--threads 4`
- Switch to single model: `--mode single`
- Increase `--min_presence` to reduce the number of genes
- Use a server with more RAM

### trimAl removes too many columns

**Symptom:** `align_trim_step.log` warns ">50% columns removed" for many genes.

**Causes and fixes:**
1. Sequences are highly divergent — try `--trimal_mode gappyout` (less aggressive)
2. A few sequences in the alignment are very divergent — inspect those FASTAs
3. The ortholog may contain paralogs not caught by OrthoFinder — inspect manually

### OrthoFinder fails with "diamond not found" (or other missing tool)

The pipeline checks for all required tools before any step runs and exits
immediately with a specific error message and install command. Common fixes:

```bash
conda install -c bioconda diamond       # for diamond / diamond_ultra_sens
conda install -c bioconda blast         # for blastp
conda install -c bioconda mmseqs2       # for mmseqs
conda install -c bioconda mafft         # for mafft
conda install -c bioconda trimal        # for trimal
conda install -c bioconda iqtree        # for iqtree3 / iqtree
```

If tools are installed but not found, ensure the correct conda environment is
activated (`conda activate orthofinder_pipeline`) before running the pipeline.

### OrthoFinder times out

**Symptom:** The pipeline exits with "OrthoFinder was killed after exceeding
the Ns timeout."

**Cause:** The search (especially BLAST all-vs-all with many proteomes) is
taking longer than the configured timeout.

**Fixes:**
1. Use `diamond` instead of `blast` (`OF_SEARCH="diamond"` in the config).
   DIAMOND is typically 500–1000× faster than BLAST for large datasets.
2. Increase the timeout: set `OF_TIMEOUT=172800` (48 h) or `OF_TIMEOUT=0`
   (no timeout) in `pipeline_runner.sh`.
3. Reduce `THREADS` — you are prompted for a value at the start of each run, or set it directly in the config block.

### IQ-TREE output is not appearing live (delayed in large bursts)

**Symptom:** Nothing prints during IQ-TREE run, then a large block of text
appears at once.

**Cause:** `stdbuf` and `unbuffer` are not installed, so IQ-TREE uses full
block buffering when writing to a pipe.

**Fix — option 1:** Install GNU coreutils (provides `stdbuf`):
```bash
# macOS
brew install coreutils

# Linux (Debian/Ubuntu)
sudo apt-get install coreutils
```

**Fix — option 2:** Install `expect` (provides `unbuffer`):
```bash
conda install -c conda-forge expect
```

After installing either tool, re-run Step 5 — live output will work
automatically.

---

## Citation

If you use this pipeline, please cite the tools it uses:

- **OrthoFinder**: Emms & Kelly (2019) *Genome Biology* 20:238
- **DIAMOND**: Buchfink et al. (2015) *Nature Methods* 12:59–60
- **MAFFT**: Katoh & Standley (2013) *Molecular Biology and Evolution* 30:772–780
- **trimAl**: Capella-Gutiérrez et al. (2009) *Bioinformatics* 25:1972–1973
- **FASconCAT-G**: Kück & Longo (2014) *Frontiers in Zoology* 11:81
- **IQ-TREE**: Minh et al. (2020) *Molecular Biology and Evolution* 37:1530–1534

---

## AI Declaration

This pipeline was developed with the assistance of
[Claude Code](https://claude.ai/code) (Anthropic), an AI-powered coding
assistant. Claude Code was used throughout the development process for:

- Designing and implementing the modular Python pipeline scripts
  (`01_qc_proteomes.py` through `06_generate_report.py`)
- Writing the interactive Bash wrapper (`pipeline_runner.sh`), including
  the unified settings screen, step-range selection, and automatic/interactive
  run modes
- Developing the Streamlit browser GUI (`app.py` and
  `pipeline_runner_gui.py`)
- Generating the multi-page PDF summary report with pipeline flow diagram
  and reference list (`06_generate_report.py`)
- Security review and remediation (path-traversal vulnerability in the GUI)
- Writing and maintaining this documentation

All code was reviewed, tested, and validated by the author on real
bacterial proteome datasets. The underlying bioinformatics tools
(OrthoFinder, DIAMOND, MAFFT, trimAl, IQ-TREE) and the scientific
decisions regarding their application are the responsibility of the author.

> Anthropic. (2025). *Claude Code* (Version claude-sonnet-4-5).
> https://claude.ai/code
