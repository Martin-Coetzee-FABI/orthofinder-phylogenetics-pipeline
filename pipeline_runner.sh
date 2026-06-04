#!/usr/bin/env bash
# ============================================================
# OrthoFinder3 Single-Copy Ortholog & Species Tree Pipeline
# ============================================================
# Edit the USER CONFIGURATION block, then run:
#   bash pipeline_runner.sh
#
# Flags:
#   bash pipeline_runner.sh --auto    # skip all confirmation prompts

# --- USER CONFIGURATION ---
INPUT_FAA_DIR="proteomes"            # Directory with .faa files
MIN_PROTEINS=500                     # QC: minimum proteins per genome
MIN_PRESENCE=1.0                     # Ortholog selection: 1.0=strict, 0.95=relaxed
THREADS=8                            # CPU threads for MAFFT, OrthoFinder, IQ-TREE
                                     # The user is prompted to confirm or change this
                                     # before each run; override by setting a fixed
                                     # number here to skip the prompt in --auto mode.

# OrthoFinder settings
OF_SEARCH="diamond"                  # Sequence search program (-S):
                                     #   diamond            fast (recommended)
                                     #   diamond_ultra_sens sensitive DIAMOND
                                     #   blast / blastp     classic NCBI BLAST
                                     #   mmseqs             ultra-fast (requires MMseqs2)
OF_GENE_TREE="msa"                   # Gene tree inference method (-M):
                                     #   msa          alignment-based (recommended for OF3)
                                     #   dendroblast  distance-based, fast (legacy)
OF_MSA="mafft"                       # MSA program for gene trees (-A):
                                     #   only used when OF_GENE_TREE=msa
                                     #   mafft / muscle / famsa
# Auto-detect a rooted species tree in the working directory.
# The first .tre or .tree file found is used as the default.
# The user is always prompted to confirm, change, or clear this before each run.
OF_SPECIES_TREE=""
for _f in *.tre *.tree; do
    if [ -f "$_f" ]; then
        OF_SPECIES_TREE="$_f"
        break
    fi
done
unset _f
OF_TIMEOUT=86400                     # Max seconds to wait for OrthoFinder (default: 24h)
                                     #   0 = no timeout

TRIMAL_MODE="auto"                   # trimAl mode: auto, gappyout, strict, strictplus
IQTREE_MODE="partition"              # IQ-TREE mode: single or partition
                                     #   partition: -p dir/ (per-gene ModelFinder, recommended)
                                     #   single:    -s dir/ (one model for all genes)
IQTREE_MODEL="LG+F+G4"              # Substitution model — applies to BOTH modes:
                                     #   LG+F+G4     fast, fixed model (recommended)
                                     #   MFP         ModelFinder per gene (slow)
                                     #   MFP+MERGE   ModelFinder + merge partitions (slowest)
BOOTSTRAP=1000                       # Ultrafast bootstrap replicates (-B)
ALRT=1000                            # SH-aLRT replicates (-alrt); set to 0 to disable

STEP_START=1                         # First step to run  (1-5; see step list below)
STEP_END=5                           # Last step to run   (1-5; PDF report always runs)
                                     #   1  QC proteomes
                                     #   2  OrthoFinder
                                     #   3  Extract single-copy orthologs
                                     #   4  Align & trim (MAFFT + trimAl)
                                     #   5  IQ-TREE species tree

# ============================================================
# Parse --auto flag
# ============================================================
AUTO_MODE="false"
for arg in "$@"; do
    if [ "$arg" = "--auto" ] || [ "$arg" = "-y" ]; then
        AUTO_MODE="true"
    fi
done

# ============================================================
# Helper: ask user to continue (skipped in auto mode)
# Reads from /dev/tty so it works even when stdout is redirected
# ============================================================
ask_continue() {
    local next_step="$1"

    if [ "$AUTO_MODE" = "true" ]; then
        echo ""
        echo "  [AUTO] Continuing automatically to: $next_step"
        echo ""
        return 0
    fi

    echo ""
    echo "  ┌──────────────────────────────────────────────────────┐"
    printf "  │  Next:  %s│\n" "$(_pad 45 "$next_step")"
    echo "  │                                                      │"
    echo "  │  [Enter]  Continue      [q]  Quit pipeline           │"
    echo "  └──────────────────────────────────────────────────────┘"
    printf "  > "

    local answer
    read -r answer < /dev/tty

    if [ "$answer" = "q" ] || [ "$answer" = "Q" ] || [ "$answer" = "quit" ]; then
        echo ""
        echo "  Pipeline paused by user. Re-run to resume from this step."
        echo "  (Previously completed steps will be detected and skipped.)"
        echo ""
        exit 0
    fi
    echo ""
}

# ============================================================
# Helper: print a step banner
# ============================================================
step_banner() {
    local step_num="$1"
    local step_name="$2"
    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    printf "║  Step %s — %-43s║\n" "$step_num" "$step_name"
    printf "║  Started: %-43s║\n" "$(date '+%Y-%m-%d %H:%M:%S')"
    echo "╚══════════════════════════════════════════════════════╝"
    echo ""
}

step_done() {
    local step_num="$1"
    echo ""
    echo "  ✓  Step $step_num done: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
}

# ============================================================
# Helper: left-pad a string to WIDTH *display* columns.
# printf's %-Ns counts bytes, not characters — multibyte UTF-8
# characters (e.g. em-dash U+2014) display as 1 column but
# take 3 bytes, causing right borders to shift.  This helper
# uses wc -m (character count) so padding is always correct.
# ============================================================
_pad() {
    local width="$1" str="$2"
    local dlen pad
    dlen=$(printf "%s" "$str" | wc -m | tr -d ' ')
    pad=$(( width - dlen ))
    [ "$pad" -lt 0 ] && pad=0
    printf "%s%*s" "$str" "$pad" ""
}

# ============================================================
# Helper: print all pipeline settings as a numbered table.
# Called by configure_all_settings; also re-used to echo the
# confirmed settings into the log after tee is set up.
# ============================================================
_show_settings() {
    local cores
    cores=$(nproc 2>/dev/null \
            || python -c "import os; print(os.cpu_count() or 4)" 2>/dev/null \
            || echo "?")
    local stree_val="${OF_SPECIES_TREE:-(infer from data)}"
    local run_mode_label
    [ "$AUTO_MODE" = "true" ] && run_mode_label="automatic" || run_mode_label="interactive"

    echo ""
    echo "  ╔══════════════════════════════════════════════════════╗"
    echo "  ║  Pipeline Settings                                   ║"
    echo "  ╚══════════════════════════════════════════════════════╝"
    echo ""
    echo "   General"
    echo "   ─────────────────────────────────────────────────────"
    printf "    [ 1]  Input proteomes dir   :  %s\n"         "$INPUT_FAA_DIR"
    printf "    [ 2]  Min proteins (QC)     :  %s\n"         "$MIN_PROTEINS"
    printf "    [ 3]  Min ortholog presence :  %s\n"         "$MIN_PRESENCE"
    printf "    [ 4]  Threads               :  %s  (%s available)\n" "$THREADS" "$cores"
    echo ""
    echo "   OrthoFinder"
    echo "   ─────────────────────────────────────────────────────"
    printf "    [ 5]  Search program  (-S)  :  %s\n"         "$OF_SEARCH"
    printf "    [ 6]  Gene tree       (-M)  :  %s\n"         "$OF_GENE_TREE"
    printf "    [ 7]  MSA program     (-A)  :  %s\n"         "$OF_MSA"
    printf "    [ 8]  Species tree    (-s)  :  %s\n"         "$stree_val"
    printf "    [ 9]  Timeout               :  %s s\n"       "$OF_TIMEOUT"
    echo ""
    echo "   Alignment & Trimming"
    echo "   ─────────────────────────────────────────────────────"
    printf "    [10]  trimAl mode           :  %s\n"         "$TRIMAL_MODE"
    echo ""
    echo "   IQ-TREE"
    echo "   ─────────────────────────────────────────────────────"
    printf "    [11]  Mode                  :  %s\n"         "$IQTREE_MODE"
    printf "    [12]  Model                 :  %s\n"         "$IQTREE_MODEL"
    printf "    [13]  Bootstrap reps        :  %s\n"         "$BOOTSTRAP"
    printf "    [14]  SH-aLRT reps          :  %s\n"         "$ALRT"
    echo ""
    echo "   Steps to run  (PDF report always runs at the end)"
    echo "   ─────────────────────────────────────────────────────"
    printf "    [16]  Start from step      :  %s  (%s)\n"  "$STEP_START" "$(_step_name "$STEP_START")"
    printf "    [17]  Stop after step      :  %s  (%s)\n"  "$STEP_END"   "$(_step_name "$STEP_END")"
    echo ""
    echo "   Pipeline"
    echo "   ─────────────────────────────────────────────────────"
    printf "    [15]  Run mode              :  %s\n"         "$run_mode_label"
    echo ""
    echo "   ─────────────────────────────────────────────────────"
    echo "   Enter number to edit a setting, or [Enter] to start."
    echo "   ─────────────────────────────────────────────────────"
    echo ""
}

# ============================================================
# Helper: unified settings editor.
# Displays all settings as a numbered table (via _show_settings)
# and lets the user edit any item by typing its number.
# Loops until the user presses Enter on an empty line.
# In --auto mode or non-interactive sessions the table is
# printed once and the function returns immediately.
# ============================================================
configure_all_settings() {

    if [ "$AUTO_MODE" = "true" ] || [ ! -t 0 ]; then
        _show_settings
        echo "   [AUTO] Using these settings — starting pipeline."
        echo ""
        return 0
    fi

    local choice="" v=""

    while true; do
        _show_settings
        printf "   > "
        choice=""
        read -r choice < /dev/tty || choice=""

        # Empty input → confirmed, start the pipeline
        [ -z "$choice" ] && break

        v=""
        case "$choice" in

          1)  printf "   Input proteomes dir [%s]: " "$INPUT_FAA_DIR"
              read -r v < /dev/tty || v=""
              [ -n "$v" ] && INPUT_FAA_DIR="$v" ;;

          2)  printf "   Min proteins (QC) [%s]: " "$MIN_PROTEINS"
              read -r v < /dev/tty || v=""
              [[ "$v" =~ ^[0-9]+$ ]] && MIN_PROTEINS="$v" ;;

          3)  printf "   Min presence [%s]  (1.0=strict, 0.95=relaxed): " "$MIN_PRESENCE"
              read -r v < /dev/tty || v=""
              [ -n "$v" ] && MIN_PRESENCE="$v" ;;

          4)  local cores
              cores=$(nproc 2>/dev/null || echo "?")
              printf "   Threads [%s] (%s available): " "$THREADS" "$cores"
              read -r v < /dev/tty || v=""
              [[ "$v" =~ ^[1-9][0-9]*$ ]] && THREADS="$v" ;;

          5)  echo "    [1] diamond  [2] diamond_ultra_sens  [3] blast  [4] mmseqs"
              printf "   Search program [%s] — select 1-4 or Enter to keep: " "$OF_SEARCH"
              read -r v < /dev/tty || v=""
              case "$v" in
                  1) OF_SEARCH="diamond" ;;
                  2) OF_SEARCH="diamond_ultra_sens" ;;
                  3) OF_SEARCH="blast" ;;
                  4) OF_SEARCH="mmseqs" ;;
              esac ;;

          6)  echo "    [1] msa (alignment-based, recommended)  [2] dendroblast (fast, legacy)"
              printf "   Gene tree method [%s] — select 1-2 or Enter to keep: " "$OF_GENE_TREE"
              read -r v < /dev/tty || v=""
              case "$v" in
                  1) OF_GENE_TREE="msa" ;;
                  2) OF_GENE_TREE="dendroblast" ;;
              esac ;;

          7)  echo "    [1] mafft (recommended)  [2] famsa (fastest)  [3] muscle"
              printf "   MSA program [%s] — select 1-3 or Enter to keep: " "$OF_MSA"
              read -r v < /dev/tty || v=""
              case "$v" in
                  1) OF_MSA="mafft" ;;
                  2) OF_MSA="famsa" ;;
                  3) OF_MSA="muscle" ;;
              esac ;;

          8)  printf "   Current: %s\n" "${OF_SPECIES_TREE:-(infer from data)}"
              printf "   Path to ROOTED Newick file, 'none' to clear, or Enter to keep: "
              read -r v < /dev/tty || v=""
              if [ "$v" = "none" ] || [ "$v" = "None" ]; then
                  OF_SPECIES_TREE=""
                  echo "   Cleared — OrthoFinder will infer the species tree."
              elif [ -n "$v" ]; then
                  if [ -f "$v" ]; then
                      OF_SPECIES_TREE="$v"
                  else
                      echo "   WARNING: file not found: $v — unchanged."
                  fi
              fi ;;

          9)  printf "   Timeout seconds [%s] (0 = no limit): " "$OF_TIMEOUT"
              read -r v < /dev/tty || v=""
              [[ "$v" =~ ^[0-9]+$ ]] && OF_TIMEOUT="$v" ;;

         10)  echo "    [1] auto (recommended)  [2] gappyout  [3] strict  [4] strictplus"
              printf "   trimAl mode [%s] — select 1-4 or Enter to keep: " "$TRIMAL_MODE"
              read -r v < /dev/tty || v=""
              case "$v" in
                  1) TRIMAL_MODE="auto" ;;
                  2) TRIMAL_MODE="gappyout" ;;
                  3) TRIMAL_MODE="strict" ;;
                  4) TRIMAL_MODE="strictplus" ;;
              esac ;;

         11)  echo "    [1] partition  (-p, per-gene model, recommended)"
              echo "    [2] single     (-s, one model for all genes, faster)"
              printf "   IQ-TREE mode [%s] — select 1-2 or Enter to keep: " "$IQTREE_MODE"
              read -r v < /dev/tty || v=""
              case "$v" in
                  1) IQTREE_MODE="partition" ;;
                  2) IQTREE_MODE="single" ;;
              esac ;;

         12)  echo "    Examples: LG+F+G4 (fast, recommended)  MFP (ModelFinder)  MFP+MERGE (slowest)"
              printf "   Model [%s]: " "$IQTREE_MODEL"
              read -r v < /dev/tty || v=""
              [ -n "$v" ] && IQTREE_MODEL="$v" ;;

         13)  printf "   Bootstrap reps [%s]: " "$BOOTSTRAP"
              read -r v < /dev/tty || v=""
              [[ "$v" =~ ^[0-9]+$ ]] && BOOTSTRAP="$v" ;;

         14)  printf "   SH-aLRT reps [%s] (0 to disable): " "$ALRT"
              read -r v < /dev/tty || v=""
              [[ "$v" =~ ^[0-9]+$ ]] && ALRT="$v" ;;

         15)  echo "    [1] interactive (confirm before each step)"
              echo "    [2] automatic   (run all steps without prompts)"
              printf "   Run mode [%s] — select 1-2 or Enter to keep: " \
                     "$([ "$AUTO_MODE" = "true" ] && echo "automatic" || echo "interactive")"
              read -r v < /dev/tty || v=""
              case "$v" in
                  1) AUTO_MODE="false" ;;
                  2) AUTO_MODE="true" ;;
              esac ;;

         16)  _show_step_descriptions
              printf "   Start from step [%s] — enter 1-5 or Enter to keep: " "$STEP_START"
              read -r v < /dev/tty || v=""
              if [[ "$v" =~ ^[1-5]$ ]]; then
                  STEP_START="$v"
                  [ "$STEP_START" -gt "$STEP_END" ] && STEP_END="$STEP_START"
              elif [ -n "$v" ]; then
                  echo "   Invalid — enter a number between 1 and 5."
              fi ;;

         17)  _show_step_descriptions
              printf "   Stop after step [%s] — enter 1-5 or Enter to keep: " "$STEP_END"
              read -r v < /dev/tty || v=""
              if [[ "$v" =~ ^[1-5]$ ]]; then
                  STEP_END="$v"
                  [ "$STEP_END" -lt "$STEP_START" ] && STEP_START="$STEP_END"
              elif [ -n "$v" ]; then
                  echo "   Invalid — enter a number between 1 and 5."
              fi ;;

          *)  echo "   Unknown option '$choice' — enter 1-17 or press Enter to start." ;;

        esac
    done

    echo "   Settings confirmed — starting pipeline."
    echo ""
}

# ============================================================
# Helper: return the display name for a pipeline step number.
# ============================================================
_step_name() {
    case "$1" in
        1) echo "QC proteomes" ;;
        2) echo "OrthoFinder" ;;
        3) echo "Extract single-copy orthologs" ;;
        4) echo "Align & trim (MAFFT + trimAl)" ;;
        5) echo "IQ-TREE species tree" ;;
        6) echo "PDF report" ;;
        *) echo "Unknown step" ;;
    esac
}

# ============================================================
# Helper: print a short description of every pipeline step,
# including what input it expects and what it produces.
# Shown when the user edits the start/stop step settings.
# ============================================================
_show_step_descriptions() {
    echo ""
    echo "   ╔══════════════════════════════════════════════════════════════════════╗"
    echo "   ║  Pipeline steps                                                      ║"
    echo "   ╚══════════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "    [1]  QC proteomes"
    echo "         Checks every .faa file for protein count, length, and duplicates."
    echo "         Flags WARN genomes (one issue) and excludes FAIL genomes."
    echo "         Needs  : proteomes/ directory with .faa files"
    echo "         Outputs: qc_results/passed_proteomes/"
    echo ""
    echo "    [2]  OrthoFinder"
    echo "         Runs DIAMOND all-vs-all search and clusters proteins into orthogroups."
    echo "         Needs  : qc_results/passed_proteomes/"
    echo "         Outputs: orthofinder_results/Results_*/"
    echo ""
    echo "    [3]  Extract single-copy orthologs"
    echo "         Selects orthogroups present in all (or most) species with 1 copy each."
    echo "         Writes one .faa file per ortholog with species names as headers."
    echo "         Needs  : orthofinder_results/Results_*/"
    echo "         Outputs: ortholog_results/gene_fastas/"
    echo ""
    echo "    [4]  Align & trim (MAFFT + trimAl)"
    echo "         Aligns each ortholog with MAFFT and removes gappy columns with trimAl."
    echo "         Needs  : ortholog_results/gene_fastas/"
    echo "         Outputs: alignment_results/trimmed/"
    echo ""
    echo "    [5]  IQ-TREE species tree"
    echo "         Infers the species phylogeny from the trimmed alignments."
    echo "         Needs  : alignment_results/trimmed/"
    echo "         Outputs: iqtree_results/species_tree.treefile"
    echo ""
    echo "    [6]  PDF report  (always runs at the end)"
    echo "         Summarises QC, orthologs, alignment stats, and the species tree."
    echo "         Uses whatever output directories are already present."
    echo "         Outputs: report/pipeline_summary.pdf"
    echo ""
    echo "   ──────────────────────────────────────────────────────────────────────"
    echo ""
}

# ============================================================
# Helper: check that a command exists and is executable;
# exit immediately with an installation hint if not found.
# Usage: require_tool <cmd> [<install_hint>]
# ============================================================
require_tool() {
    local cmd="$1"
    local hint="${2:-install via conda (bioconda channel)}"
    if ! command -v "$cmd" > /dev/null 2>&1; then
        echo ""
        echo "  ERROR: Required tool not found in PATH: $cmd"
        echo "         Hint: $hint"
        echo ""
        exit 1
    fi
}

PIPELINE_LOG="pipeline_run.log"

# ============================================================
# Welcome banner and settings editor — runs BEFORE the tee
# redirect so interactive prompts are never doubled on screen.
# ============================================================
echo ""
echo "╔══════════════════════════════════════════════════════╗"
echo "║     OrthoFinder3 Species Tree Pipeline               ║"
echo "╚══════════════════════════════════════════════════════╝"
echo ""
echo "  Working directory : $(pwd)"
echo "  Master log        : $PIPELINE_LOG"

# Non-interactive sessions (no terminal) run automatically
if [ ! -t 0 ] && [ "$AUTO_MODE" != "true" ]; then
    AUTO_MODE="true"
    echo "  Non-interactive session detected — running in automatic mode."
fi

# Unified settings review: show all 15 items, user edits by number
configure_all_settings

echo ""

# ============================================================
# Start logging — everything from here goes to screen + log
# ============================================================
exec > >(tee -a "$PIPELINE_LOG") 2>&1

set -euo pipefail

# ============================================================
# Re-echo confirmed settings into the log (they were shown
# before the tee redirect and would otherwise be missing).
# ============================================================
echo "  Settings used for this run:"
_show_settings
echo ""

# ============================================================
# Write run_settings.json — read by 06_generate_report.py to
# populate the pipeline flow diagram in the PDF report.
# ============================================================
cat > run_settings.json << SETTINGS_EOF
{
  "input_dir":       "$INPUT_FAA_DIR",
  "min_proteins":    $MIN_PROTEINS,
  "min_presence":    "$MIN_PRESENCE",
  "threads":         $THREADS,
  "of_search":       "$OF_SEARCH",
  "of_gene_tree":    "$OF_GENE_TREE",
  "of_msa":          "$OF_MSA",
  "of_species_tree": "$OF_SPECIES_TREE",
  "of_timeout":      $OF_TIMEOUT,
  "trimal_mode":     "$TRIMAL_MODE",
  "iqtree_mode":     "$IQTREE_MODE",
  "iqtree_model":    "$IQTREE_MODEL",
  "bootstrap":       $BOOTSTRAP,
  "alrt":            $ALRT,
  "step_start":      $STEP_START,
  "step_end":        $STEP_END
}
SETTINGS_EOF
echo "  Run settings written to run_settings.json"
echo ""

# ============================================================
# Check that a conda environment is active
# Tools installed via conda will not be on PATH until the
# environment is activated with: conda activate <env_name>
# ============================================================
_conda_env="${CONDA_DEFAULT_ENV:-}"
if [ -z "$_conda_env" ] || [ "$_conda_env" = "base" ]; then
    echo ""
    echo "  ╔══════════════════════════════════════════════════════╗"
    echo "  ║  WARNING: no conda environment is active             ║"
    echo "  ║                                                      ║"
    echo "  ║  Tools installed in a conda environment will not     ║"
    echo "  ║  be found until that environment is activated.       ║"
    echo "  ║                                                      ║"
    echo "  ║  Before running this script, activate your env:      ║"
    echo "  ║    conda activate orthofinder3                       ║"
    echo "  ║  (substitute the name of your actual environment)   ║"
    echo "  ╚══════════════════════════════════════════════════════╝"
    echo ""
    if [ "$AUTO_MODE" != "true" ] && [ -t 0 ]; then
        printf "  Press Enter to continue anyway, or Ctrl+C to abort: "
        read -r _ < /dev/tty || true
        echo ""
    fi
else
    echo "  Active conda environment: $_conda_env"
    echo ""
fi
unset _conda_env

# ============================================================
# Pre-flight: tool version checks and hard exits for missing tools
# ============================================================
echo "  Tool versions:"
printf "    Python     : "; python --version 2>&1 || true
printf "    OrthoFinder: "; orthofinder --version 2>/dev/null || echo "not found in PATH"
printf "    MAFFT      : "; mafft --version 2>&1 | head -1 || true
printf "    trimAl     : "; trimal -version 2>&1 | head -1 || trimal --version 2>&1 | head -1 || echo "not found in PATH"
printf "    IQ-TREE    : "; (iqtree3 --version 2>/dev/null | head -1) || (iqtree --version 2>/dev/null | head -1) || echo "not found in PATH"
echo ""

# Hard stop if required tools are missing
require_tool "orthofinder"  "conda install -c bioconda orthofinder"
require_tool "mafft"        "conda install -c bioconda mafft"
require_tool "trimal"       "conda install -c bioconda trimal"

# Verify IQ-TREE (either iqtree3 or iqtree must exist)
if ! command -v "iqtree3" > /dev/null 2>&1 && ! command -v "iqtree" > /dev/null 2>&1; then
    echo ""
    echo "  ERROR: Required tool not found in PATH: iqtree3 (or iqtree)"
    echo "         Hint: conda install -c bioconda iqtree"
    echo ""
    exit 1
fi

# Verify the sequence search tool
case "$OF_SEARCH" in
    diamond|diamond_ultra_sens)
        require_tool "diamond" "conda install -c bioconda diamond" ;;
    blast|blastp)
        require_tool "blastp"  "conda install -c bioconda blast" ;;
    mmseqs)
        require_tool "mmseqs"  "conda install -c bioconda mmseqs2" ;;
esac

# Verify alignment tool if using MSA method
if [ "$OF_GENE_TREE" = "msa" ]; then
    require_tool "$OF_MSA" "conda install -c bioconda $OF_MSA"
fi

# Validate species tree file if supplied
if [ -n "$OF_SPECIES_TREE" ]; then
    if [ ! -f "$OF_SPECIES_TREE" ]; then
        echo ""
        echo "  ERROR: Species tree file not found: $OF_SPECIES_TREE"
        echo "         Check OF_SPECIES_TREE in the USER CONFIGURATION block."
        echo ""
        exit 1
    fi
    if [ ! -s "$OF_SPECIES_TREE" ]; then
        echo ""
        echo "  ERROR: Species tree file is empty: $OF_SPECIES_TREE"
        echo ""
        exit 1
    fi
    echo "  Species tree validated: $OF_SPECIES_TREE"
fi

# ============================================================
# Pre-flight: disk space check
# OrthoFinder all-vs-all can need ~10x the input .faa size.
# ============================================================
INPUT_SIZE_KB=$(du -sk "$INPUT_FAA_DIR" 2>/dev/null | cut -f1 || echo 0)
REQUIRED_KB=$(( INPUT_SIZE_KB * 10 ))
AVAIL_KB=$(df -k . 2>/dev/null | awk 'NR==2 {print $4}' || echo 999999999)
if [ "$REQUIRED_KB" -gt "$AVAIL_KB" ]; then
    echo ""
    echo "  WARNING: Estimated disk space needed: ~$(( REQUIRED_KB / 1024 )) MB"
    echo "           Available disk space:          $(( AVAIL_KB / 1024 )) MB"
    echo "           The run may fail due to insufficient disk space."
    echo "           Press Ctrl+C to abort, or Enter to continue anyway."
    if [ "$AUTO_MODE" != "true" ] && [ -t 0 ]; then
        read -r _ < /dev/tty || true
    fi
else
    echo "  Disk space check: OK (~$(( REQUIRED_KB / 1024 )) MB needed, $(( AVAIL_KB / 1024 )) MB available)"
fi
echo ""

# ============================================================
# Safety cap: if THREADS still exceeds available cores
# (can happen in --auto mode with a stale config value)
# ============================================================
MAX_CORES=$(nproc 2>/dev/null || python -c "import os; print(os.cpu_count() or 4)" 2>/dev/null || echo 4)
if [ "$THREADS" -gt "$MAX_CORES" ]; then
    echo "  WARNING: THREADS=$THREADS exceeds available cores ($MAX_CORES)."
    echo "           Capping to $MAX_CORES. Edit THREADS in pipeline_runner.sh to suppress this."
    THREADS=$MAX_CORES
fi
echo "  Using $THREADS thread(s) (of $MAX_CORES available)."

# ============================================================
# STEP 1 — QC proteomes
# ============================================================
STEP1_MARKER="qc_results/.step_complete"

if [ "$STEP_START" -gt 1 ]; then
    echo "  [SKIP] Step 1 — QC proteomes (not in selected range: steps $STEP_START–$STEP_END)"
elif [ -f "$STEP1_MARKER" ]; then
    echo "  [SKIP] Step 1 — QC proteomes (already complete)"
else
    step_banner "1" "QC proteomes"
    python 01_qc_proteomes.py \
        --input_dir "$INPUT_FAA_DIR" \
        --min_proteins "$MIN_PROTEINS" \
        --output_dir qc_results
    touch "$STEP1_MARKER"
    step_done "1"
fi

[ "$STEP_END" -ge 2 ] && ask_continue "Step 2 — Run OrthoFinder3"

# ============================================================
# STEP 2 — Run OrthoFinder3
# ============================================================
ORTHOFINDER_DIR=""

if [ "$STEP_START" -gt 2 ] || [ "$STEP_END" -lt 2 ]; then
    echo "  [SKIP] Step 2 — OrthoFinder (not in selected range: steps $STEP_START–$STEP_END)"
else
    step_banner "2" "Run OrthoFinder3"

    # Build optional species-tree argument
    SPECIES_TREE_ARG=""
    if [ -n "$OF_SPECIES_TREE" ]; then
        SPECIES_TREE_ABS="$(cd "$(dirname "$OF_SPECIES_TREE")" && pwd)/$(basename "$OF_SPECIES_TREE")"
        SPECIES_TREE_ARG="--species_tree $SPECIES_TREE_ABS"
        echo "  Species tree (-s): $SPECIES_TREE_ABS"
    fi

    # Build optional timeout argument
    TIMEOUT_ARG=""
    if [ "$OF_TIMEOUT" -gt 0 ] 2>/dev/null; then
        TIMEOUT_ARG="--timeout $OF_TIMEOUT"
    fi

    python 02_run_orthofinder.py \
        --input_dir qc_results/passed_proteomes \
        --output_dir orthofinder_results \
        --threads "$THREADS" \
        --search_prog "$OF_SEARCH" \
        --gene_tree_method "$OF_GENE_TREE" \
        --msa_prog "$OF_MSA" \
        $TIMEOUT_ARG \
        $SPECIES_TREE_ARG

    ORTHOFINDER_DIR=$(ls -d orthofinder_results/Results_* 2>/dev/null | tail -1)
    if [ -z "$ORTHOFINDER_DIR" ]; then
        echo "ERROR: Could not find OrthoFinder Results_* directory" >&2
        exit 1
    fi
    echo "  OrthoFinder results directory: $ORTHOFINDER_DIR"
    step_done "2"
fi

# Always try to locate Results_* in case Step 2 was skipped
if [ -z "$ORTHOFINDER_DIR" ]; then
    ORTHOFINDER_DIR=$(ls -d orthofinder_results/Results_* 2>/dev/null | tail -1)
fi

[ "$STEP_END" -ge 3 ] && ask_continue "Step 3 — Extract single-copy orthologs"

# ============================================================
# STEP 3 — Extract orthologs
# ============================================================
STEP3_MARKER="ortholog_results/.step_complete"

if [ "$STEP_START" -gt 3 ] || [ "$STEP_END" -lt 3 ]; then
    echo "  [SKIP] Step 3 — Extract orthologs (not in selected range: steps $STEP_START–$STEP_END)"
elif [ -f "$STEP3_MARKER" ]; then
    echo "  [SKIP] Step 3 — Extract orthologs (already complete)"
else
    # Dependency check: need OrthoFinder results
    if [ -z "$ORTHOFINDER_DIR" ] || [ ! -d "$ORTHOFINDER_DIR" ]; then
        echo "  ERROR: Step 3 requires OrthoFinder results in orthofinder_results/Results_*/"
        echo "         Run Step 2 first, or set STEP_START=1 to run from the beginning."
        exit 1
    fi
    step_banner "3" "Extract single-copy orthologs"
    python 03_extract_orthologs.py \
        --orthofinder_dir "$ORTHOFINDER_DIR" \
        --faa_dir qc_results/passed_proteomes \
        --min_presence "$MIN_PRESENCE" \
        --output_dir ortholog_results
    touch "$STEP3_MARKER"
    step_done "3"
fi

[ "$STEP_END" -ge 4 ] && ask_continue "Step 4 — Align and trim (MAFFT + trimAl)"

# ============================================================
# STEP 4 — Align and trim
# ============================================================
STEP4_MARKER="alignment_results/.step_complete"

if [ "$STEP_START" -gt 4 ] || [ "$STEP_END" -lt 4 ]; then
    echo "  [SKIP] Step 4 — Align & trim (not in selected range: steps $STEP_START–$STEP_END)"
elif [ -f "$STEP4_MARKER" ]; then
    echo "  [SKIP] Step 4 — Align & trim (already complete)"
else
    # Dependency check: need gene FASTAs from Step 3
    if [ ! -d "ortholog_results/gene_fastas" ] || \
       [ -z "$(ls -A ortholog_results/gene_fastas/*.faa 2>/dev/null)" ]; then
        echo "  ERROR: Step 4 requires ortholog_results/gene_fastas/*.faa"
        echo "         Run Step 3 first, or set STEP_START to 3 or lower."
        exit 1
    fi
    step_banner "4" "Align and trim (MAFFT + trimAl)"
    python 04_align_trim_concat.py \
        --gene_fasta_dir ortholog_results/gene_fastas \
        --output_dir alignment_results \
        --threads "$THREADS" \
        --trimal_mode "$TRIMAL_MODE"
    touch "$STEP4_MARKER"
    step_done "4"
fi

[ "$STEP_END" -ge 5 ] && ask_continue "Step 5 — IQ-TREE3 species tree"

# ============================================================
# STEP 5 — IQ-TREE species tree
# ============================================================
STEP5_MARKER="iqtree_results/.step_complete"

if [ "$STEP_START" -gt 5 ] || [ "$STEP_END" -lt 5 ]; then
    echo "  [SKIP] Step 5 — IQ-TREE (not in selected range: steps $STEP_START–$STEP_END)"
elif [ -f "$STEP5_MARKER" ]; then
    echo "  [SKIP] Step 5 — IQ-TREE (already complete)"
else
    # Dependency check: need trimmed alignments from Step 4
    if [ ! -d "alignment_results/trimmed" ] || \
       [ -z "$(ls -A alignment_results/trimmed/*.faa 2>/dev/null)" ]; then
        echo "  ERROR: Step 5 requires alignment_results/trimmed/*.faa"
        echo "         Run Step 4 first, or set STEP_START to 4 or lower."
        exit 1
    fi
    step_banner "5" "IQ-TREE3 species tree"
    python 05_run_iqtree.py \
        --trimmed_dir alignment_results/trimmed \
        --output_dir iqtree_results \
        --mode "$IQTREE_MODE" \
        --model "$IQTREE_MODEL" \
        --bootstrap "$BOOTSTRAP" \
        --alrt "$ALRT" \
        --threads "$THREADS"
    touch "$STEP5_MARKER"
    step_done "5"
fi

ask_continue "Step 6 — Generate PDF summary report"

# ============================================================
# STEP 6 — Generate PDF summary report
# ============================================================
step_banner "6" "Generate PDF summary report"

python 06_generate_report.py \
    --output_dir report \
    --settings_file run_settings.json

step_done "6"

# ============================================================
# Pipeline complete
# ============================================================
echo ""
echo "╔══════════════════════════════════════════════════════╗"
echo "║  Pipeline complete!                                  ║"
printf "║  Finished: %-42s║\n" "$(date '+%Y-%m-%d %H:%M:%S')"
echo "╠══════════════════════════════════════════════════════╣"
echo "║  Species tree → iqtree_results/species_tree.treefile ║"
echo "║  PDF report     → report/pipeline_summary.pdf        ║"
echo "║  Gene FASTAs    → ortholog_results/gene_fastas/      ║"
echo "║  Summary tables → ortholog_results/tables/           ║"
printf "║  Master log     → %-35s║\n" "$PIPELINE_LOG"
echo "╚══════════════════════════════════════════════════════╝"
echo ""
