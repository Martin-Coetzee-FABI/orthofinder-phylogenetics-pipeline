#!/usr/bin/env python3
"""
app.py — Streamlit GUI for the OrthoFinder3 Species Tree Pipeline
=================================================================
Launch from the pipeline directory with:
    streamlit run app.py

Requirements (add to conda env):
    conda install -c conda-forge streamlit psutil
    # or: pip install "streamlit>=1.37" psutil
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import streamlit as st

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="OrthoFinder3 Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
/* Dark terminal-style log box */
.log-box {
    background-color: #0d1117;
    color: #c9d1d9;
    font-family: "SFMono-Regular", "Cascadia Code", "Consolas", "Courier New", monospace;
    font-size: 12.5px;
    line-height: 1.55;
    padding: 16px 20px;
    border-radius: 8px;
    height: 460px;
    overflow-y: auto;
    white-space: pre-wrap;
    word-break: break-all;
    border: 1px solid #30363d;
    margin-top: 4px;
}
/* Metric card row spacing */
div[data-testid="metric-container"] { padding: 10px 0; }
/* Make sidebar labels a touch lighter */
.css-1d391kg { font-size: 13px; }
</style>
""", unsafe_allow_html=True)

# ── Paths & constants ─────────────────────────────────────────────────────────
PIPELINE_LOG  = Path("pipeline_run.log")
PID_FILE      = Path(".pipeline_gui.pid")
CONFIG_FILE   = Path(".pipeline_gui_config.json")
RUNNER_SCRIPT = Path("pipeline_runner_gui.py")

# Resolved once at startup — used as the containment boundary for user paths.
_CWD = Path(".").resolve()

# Accepted file extensions for Newick tree files.
_TREE_SUFFIXES = frozenset({".tre", ".tree", ".nwk", ".newick", ".txt", ".treefile"})


def _sanitise_tree_path(raw: str) -> "Path | None":
    """Canonicalise and validate a user-supplied Newick tree-file path.

    Security notes
    --------------
    * The resolved path must lie **within the current working directory**
      (``_CWD``).  Any attempt to escape via ``../``, absolute paths to
      system files, or symbolic links that point outside ``_CWD`` is
      rejected, preventing path-traversal attacks.
    * Only file extensions from ``_TREE_SUFFIXES`` are accepted; other
      extensions (e.g. ``.sh``, ``.py``) are rejected outright.
    * This is the **single point** at which the untrusted string from
      ``st.text_input`` is converted to a ``Path``; all subsequent
      filesystem operations on tree files use the returned object, never
      the raw string.

    Returns the resolved, validated ``Path``, or ``None`` if validation
    fails for any reason.
    """
    raw = (raw or "").strip()
    if not raw:
        return None

    candidate = Path(raw)
    # Make relative paths relative to the working directory, not the OS cwd
    # (which may differ if the script was launched from a different shell dir).
    if not candidate.is_absolute():
        candidate = _CWD / candidate

    try:
        resolved = candidate.resolve()          # removes .., resolves symlinks
        resolved.relative_to(_CWD)             # raises ValueError if outside CWD
    except (ValueError, OSError, RuntimeError):
        return None                             # path traversal attempt or OS error

    if resolved.suffix.lower() not in _TREE_SUFFIXES:
        return None                             # extension not in allowlist

    return resolved

STEP_MARKERS = {
    1: Path("qc_results/.step_complete"),
    3: Path("ortholog_results/.step_complete"),
    4: Path("alignment_results/.step_complete"),
    5: Path("iqtree_results/.step_complete"),
}

# Step descriptions shown in selection widgets and the Help tab
STEP_DESCRIPTIONS = {
    1: "1 · QC proteomes — check .faa files, flag/exclude low-quality genomes",
    2: "2 · OrthoFinder — DIAMOND all-vs-all; orthogroup clustering",
    3: "3 · Extract orthologs — select single-copy genes; write per-gene FASTAs",
    4: "4 · Align & trim — MAFFT alignment + trimAl column trimming",
    5: "5 · IQ-TREE — infer species tree from trimmed alignments",
}

STEP_INPUTS = {
    1: "proteomes/  (.faa files)",
    2: "qc_results/passed_proteomes/",
    3: "orthofinder_results/Results_*/",
    4: "ortholog_results/gene_fastas/",
    5: "alignment_results/trimmed/",
}

STEPS = [
    (1, "QC"),
    (2, "OrthoFinder"),
    (3, "Orthologs"),
    (4, "Align & trim"),
    (5, "IQ-TREE3"),
    (6, "Report"),
]

STEP_ICONS = {
    "done":    "✅",
    "running": "🔄",
    "pending": "⬜",
    "failed":  "❌",
    "skipped": "⏭️",
}

# ── Session state defaults ────────────────────────────────────────────────────
_defaults = {"log_offset": 0}
for _k, _v in _defaults.items():
    if _k not in st.session_state:
        st.session_state[_k] = _v

# ═════════════════════════════════════════════════════════════════════════════
# Helper functions
# ═════════════════════════════════════════════════════════════════════════════

def available_cores() -> int:
    return os.cpu_count() or 4


def detect_tree_file() -> str:
    """Return the first .tre or .tree file found in the working directory."""
    for pat in ("*.tre", "*.tree"):
        hits = sorted(Path(".").glob(pat))
        if hits:
            return str(hits[0])
    return ""


def tool_ok(name: str) -> bool:
    return shutil.which(name) is not None


def active_conda_env() -> str:
    """Return the name of the currently active conda environment, or ''."""
    return os.environ.get("CONDA_DEFAULT_ENV", "")


@st.cache_data(ttl=60)
def find_tool_in_conda(tool: str) -> str:
    """
    Search all conda environments for a tool that is not in the active PATH.
    Returns the environment name if found, or '' if not found anywhere.
    Cached for 60 s so repeated calls during a session are fast.
    """
    conda_exe = shutil.which("conda") or shutil.which("mamba")
    if not conda_exe:
        return ""
    try:
        result = subprocess.run(
            [conda_exe, "env", "list", "--json"],
            capture_output=True, text=True, timeout=15,
        )
        envs = json.loads(result.stdout).get("envs", [])
        for env_path in envs:
            # bin/ on Mac/Linux; Scripts/ on Windows
            for bin_dir in ("bin", "Scripts"):
                p = Path(env_path) / bin_dir / tool
                if p.exists() and os.access(str(p), os.X_OK):
                    return Path(env_path).name
    except Exception:
        pass
    return ""


def pipeline_is_running() -> bool:
    """Derive running state from the PID file — always accurate after crashes."""
    if not PID_FILE.exists():
        return False
    try:
        pid = int(PID_FILE.read_text().strip())
        os.kill(pid, 0)   # signal 0 = check existence only
        return True
    except (ValueError, ProcessLookupError, PermissionError, OSError):
        PID_FILE.unlink(missing_ok=True)
        return False


def get_step_status() -> dict[int, str]:
    """Derive step completion from marker files and OrthoFinder output."""
    s = {i: "pending" for i in range(1, 7)}

    if STEP_MARKERS[1].exists():
        s[1] = "done"

    of_results = Path("orthofinder_results")
    if of_results.exists():
        of_dirs = sorted(of_results.glob("Results_*"))
        gc = of_dirs[-1] / "Orthogroups" / "Orthogroups.GeneCount.tsv" if of_dirs else None
        if gc and gc.exists():
            s[2] = "done"

    for n in [3, 4, 5]:
        if STEP_MARKERS[n].exists():
            s[n] = "done"

    if Path("report/pipeline_summary.pdf").exists():
        s[6] = "done"

    if pipeline_is_running():
        for i in range(1, 7):
            if s[i] == "pending":
                s[i] = "running"
                break

    return s


def read_log_tail() -> str:
    """Read the last ~60 KB of the pipeline log efficiently."""
    if not PIPELINE_LOG.exists():
        return "(Log file not yet created. Click ▶ Run Pipeline to start.)"
    try:
        with open(PIPELINE_LOG, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 61_440))   # last ~60 KB
            raw = f.read()
        return raw.decode("utf-8", errors="replace")
    except Exception as e:
        return f"(Error reading log: {e})"


def stop_pipeline():
    """Kill the full pipeline process tree cross-platform."""
    if not PID_FILE.exists():
        return
    try:
        pid = int(PID_FILE.read_text().strip())
    except (ValueError, OSError):
        PID_FILE.unlink(missing_ok=True)
        return

    try:
        import psutil
        try:
            parent = psutil.Process(pid)
            children = parent.children(recursive=True)
            targets = [parent] + children
        except psutil.NoSuchProcess:
            targets = []

        for p in targets:
            try:
                p.terminate()
            except psutil.NoSuchProcess:
                pass

        _, alive = psutil.wait_procs(targets, timeout=5)
        for p in alive:
            try:
                p.kill()
            except psutil.NoSuchProcess:
                pass

    except ImportError:
        # psutil not available — use OS-specific fallback
        try:
            if os.name == "nt":
                subprocess.run(
                    ["taskkill", "/F", "/T", "/PID", str(pid)],
                    capture_output=True,
                )
            else:
                os.killpg(os.getpgid(pid), 15)  # SIGTERM to process group
        except Exception:
            try:
                os.kill(pid, 15)
            except Exception:
                pass

    PID_FILE.unlink(missing_ok=True)

    with open(PIPELINE_LOG, "a") as f:
        f.write(f"\n[GUI] Pipeline stopped by user — {datetime.now():%Y-%m-%d %H:%M:%S}\n")


def launch_pipeline(cfg: dict):
    """Write config and start the runner as a detached subprocess."""
    if not RUNNER_SCRIPT.exists():
        st.error(f"Runner script not found: {RUNNER_SCRIPT}\n"
                 "Make sure pipeline_runner_gui.py is in the same folder as app.py.")
        return

    CONFIG_FILE.write_text(json.dumps(cfg, indent=2))

    # Reset log-tail offset so we show new output from this run
    if PIPELINE_LOG.exists():
        st.session_state.log_offset = PIPELINE_LOG.stat().st_size
    else:
        st.session_state.log_offset = 0

    popen_kwargs: dict = {
        "stdout": subprocess.DEVNULL,
        "stderr": subprocess.DEVNULL,
    }
    if os.name == "nt":
        popen_kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
    else:
        popen_kwargs["start_new_session"] = True

    proc = subprocess.Popen(
        [sys.executable, str(RUNNER_SCRIPT), str(CONFIG_FILE)],
        **popen_kwargs,
    )
    # Write PID immediately as a backup (runner also writes it on startup)
    PID_FILE.write_text(str(proc.pid))


@st.cache_data(ttl=60)
def check_tools(
    search_prog: str, gene_tree_method: str, msa_prog: str
) -> list[tuple[str, bool, str, str]]:
    """
    Return list of (tool, in_path, install_hint, conda_env_if_found_elsewhere).
    - in_path         : True if the tool is on the current PATH
    - install_hint    : conda install command shown when the tool is missing everywhere
    - conda_env_if_found_elsewhere : non-empty string if the tool exists in a conda
                        environment that is not currently active; empty otherwise.
    Cached for 60 s so repeated page reruns don't re-run conda env list every time.
    """
    search_exe = {"diamond": "diamond", "diamond_ultra_sens": "diamond",
                  "blast": "blastp", "blastp": "blastp",
                  "mmseqs": "mmseqs"}.get(search_prog, search_prog)

    candidates = [
        ("orthofinder", "conda install -c bioconda orthofinder"),
        (search_exe,    f"conda install -c bioconda {search_exe}"),
        ("mafft",       "conda install -c bioconda mafft"),
        ("trimal",      "conda install -c bioconda trimal"),
    ]
    if gene_tree_method == "msa" and msa_prog != "mafft":
        candidates.append((msa_prog, f"conda install -c bioconda {msa_prog}"))
    iqtree_name = "iqtree3" if tool_ok("iqtree3") else "iqtree"
    candidates.append((iqtree_name, "conda install -c bioconda iqtree"))

    results = []
    for tool, hint in candidates:
        in_path = tool_ok(tool)
        found_in_env = "" if in_path else find_tool_in_conda(tool)
        results.append((tool, in_path, hint, found_in_env))
    return results


def load_tsv(path: str):
    """Load a TSV to DataFrame, or return None gracefully."""
    p = Path(path)
    if not p.exists():
        return None
    try:
        import pandas as pd
        return pd.read_csv(p, sep="\t")
    except Exception:
        return None


def _extract_iqtree_stat(text: str, pattern: str, flags=0) -> str:
    m = re.search(pattern, text, flags)
    return m.group(1) if m else "—"


# ═════════════════════════════════════════════════════════════════════════════
# Sidebar — Configuration
# ═════════════════════════════════════════════════════════════════════════════

with st.sidebar:
    st.markdown("## 🧬 OrthoFinder3 Pipeline")
    st.caption(f"📂 `{Path('.').resolve()}`")
    st.divider()

    # ── Input ─────────────────────────────────────────────────────────────
    st.markdown("#### 📁 Input")
    input_dir = st.text_input(
        "Proteomes directory", value="proteomes",
        help="Folder with one .faa file per genome (amino acid sequences).",
    )
    faa_files = sorted(Path(input_dir).glob("*.faa")) \
        if Path(input_dir).is_dir() else []
    if faa_files:
        st.caption(f"✅ {len(faa_files)} .faa file{'s' if len(faa_files) != 1 else ''} found")
    elif Path(input_dir).is_dir():
        st.caption("⚠️ No .faa files found in this directory")
    else:
        st.caption("⚠️ Directory not found")

    st.divider()

    # ── QC ────────────────────────────────────────────────────────────────
    st.markdown("#### 🔬 QC")
    min_proteins = st.number_input(
        "Min proteins per genome", value=500,
        min_value=1, max_value=50_000, step=50,
        help="Genomes with fewer proteins are excluded.",
    )
    min_presence = st.slider(
        "Ortholog presence threshold", 0.80, 1.00, 1.00, 0.01,
        help="Fraction of species that must have the gene. "
             "1.0 = all species required (strict). 0.95 = 95% required.",
    )
    if min_presence == 1.0:
        st.caption("🔒 Strict — gene must be present in **all** species")
    else:
        st.caption(f"🔓 Relaxed — gene present in ≥ {min_presence*100:.0f}% of species")

    st.divider()

    # ── OrthoFinder ───────────────────────────────────────────────────────
    st.markdown("#### ⚙️ OrthoFinder")
    search_prog = st.selectbox(
        "Search method",
        ["diamond", "diamond_ultra_sens", "blast", "mmseqs"], index=0,
        help="DIAMOND is fast and recommended. BLAST is slower but classic.",
    )
    gene_tree_method = st.selectbox(
        "Gene tree method", ["msa", "dendroblast"], index=0,
        help="MSA = alignment-based (OrthoFinder3 default, more accurate). "
             "Dendroblast = fast distance-based (legacy).",
    )
    msa_prog = st.selectbox(
        "MSA program", ["mafft", "famsa", "muscle"], index=0,
        disabled=(gene_tree_method != "msa"),
        help="Alignment tool for gene trees. Only used when gene tree method = MSA.",
    )

    default_tree = detect_tree_file()
    species_tree = st.text_input(
        "Guide species tree (optional)", value=default_tree,
        help="Rooted Newick file (.tre / .tree). "
             "Auto-detected from the working directory. "
             "Must be located inside the pipeline's working directory. "
             "Leave blank to let OrthoFinder infer the tree.",
    )
    # Sanitise once; all subsequent code uses this validated Path (or None).
    _validated_tree: "Path | None" = _sanitise_tree_path(species_tree)
    if species_tree:
        if _validated_tree is None:
            st.caption(
                "⚠️ Invalid path — use a .tre / .tree / .nwk file "
                "inside the working directory"
            )
        elif _validated_tree.exists():
            st.caption(f"✅ {_validated_tree.name}")
        else:
            st.caption("⚠️ File not found")

    of_timeout = st.number_input(
        "OrthoFinder timeout (seconds)", value=86_400, min_value=0, step=3_600,
        help="Kill OrthoFinder if it exceeds this time. 0 = no limit. Default = 24 h.",
    )

    st.divider()

    # ── Threads ───────────────────────────────────────────────────────────
    st.markdown("#### 🧵 Threads")
    max_cores = available_cores()
    threads = st.slider(
        "CPU threads", min_value=1, max_value=max_cores, value=min(8, max_cores),
        help=f"This machine has {max_cores} cores. "
             "Leave some free for other users on a shared server.",
    )
    st.caption(f"{max_cores} cores available on this machine")

    st.divider()

    # ── Alignment & Tree ──────────────────────────────────────────────────
    st.markdown("#### ✂️ Alignment & Tree")
    trimal_mode = st.selectbox(
        "trimAl mode", ["auto", "gappyout", "strict", "strictplus"], index=0,
        help="auto (-automated1) is recommended for most datasets.",
    )
    iqtree_mode = st.selectbox(
        "IQ-TREE mode", ["partition", "single"], index=0,
        help="Partition = per-gene ModelFinder (recommended, slower). "
             "Single = one model for all genes (faster).",
    )
    iqtree_model = st.text_input(
        "Model", value="LG+F+G4",
        help=(
            "Substitution model passed to IQ-TREE -m. Applies to **both** modes.\n\n"
            "- `LG+F+G4` — fast, fixed model (recommended default)\n"
            "- `MFP` — ModelFinder selects best model per gene (slow)\n"
            "- `MFP+MERGE` — ModelFinder + merge similar partitions (slowest)\n\n"
            "LG+F+G4 is suitable for most bacterial/archaeal proteome analyses."
        ),
    )
    bootstrap = st.number_input(
        "Bootstrap replicates (-B)", value=1000, min_value=100, step=100,
        help="Ultrafast bootstrap. 1000 is standard for publication.",
    )
    alrt = st.number_input(
        "SH-aLRT replicates (-alrt)", value=1000, min_value=0, step=100,
        help="Set to 0 to disable SH-aLRT. 1000 is recommended.",
    )

    st.divider()

    # ── Step selection ────────────────────────────────────────────────────
    st.markdown("#### 🔢 Steps to run")
    st.caption("The PDF report (Step 6) always runs at the end.")

    step_start = st.selectbox(
        "Start from step",
        options=list(STEP_DESCRIPTIONS.keys()),
        format_func=lambda x: STEP_DESCRIPTIONS[x],
        index=0,
        help=(
            "Skip all steps before this one — the required input must already exist.\n\n"
            + "\n".join(
                f"Step {k}: needs  {v}"
                for k, v in STEP_INPUTS.items()
            )
        ),
    )

    valid_end_options = [s for s in STEP_DESCRIPTIONS if s >= step_start]
    step_end = st.selectbox(
        "Stop after step",
        options=valid_end_options,
        format_func=lambda x: STEP_DESCRIPTIONS[x],
        index=len(valid_end_options) - 1,
        help=(
            "Stop the pipeline after this step. "
            "The PDF report always runs after the final selected step."
        ),
    )

    if step_start > 1 or step_end < 5:
        st.caption(
            f"⚡ Running steps {step_start}–{step_end} only. "
            "Ensure the required input directories already exist."
        )

# ── Assemble config dict ──────────────────────────────────────────────────────
cfg = {
    "input_dir":        input_dir,
    # Use the sanitised, canonical path string — never the raw user input.
    "species_tree":     str(_validated_tree) if _validated_tree else "",
    "min_proteins":     int(min_proteins),
    "min_presence":     float(min_presence),
    "search_prog":      search_prog,
    "gene_tree_method": gene_tree_method,
    "msa_prog":         msa_prog,
    "timeout":          int(of_timeout),
    "threads":          int(threads),
    "trimal_mode":      trimal_mode,
    "iqtree_mode":      iqtree_mode,
    "iqtree_model":     iqtree_model,
    "bootstrap":        int(bootstrap),
    "alrt":             int(alrt),
    "step_start":       int(step_start),
    "step_end":         int(step_end),
}

# ═════════════════════════════════════════════════════════════════════════════
# Main area
# ═════════════════════════════════════════════════════════════════════════════

st.title("🧬 OrthoFinder3 Species Tree Pipeline")
tab_run, tab_results, tab_help = st.tabs(["🚀 Run Pipeline", "📊 Results", "❓ Help"])

# ════════════════════════════════════════════
# Tab 1 — Run Pipeline
# ════════════════════════════════════════════
with tab_run:

    # ── Step progress indicators ───────────────────────────────────────────
    step_status = get_step_status()
    step_cols = st.columns(len(STEPS))
    for col, (num, name) in zip(step_cols, STEPS):
        s = step_status[num]
        # Dim steps outside the selected range
        if num != 6 and (num < step_start or num > step_end):
            col.markdown(f"**⏭️ {name}**")
            col.caption(":gray[Skipped]")
            continue
        icon = STEP_ICONS[s]
        col.markdown(f"**{icon} {name}**")
        colour = {"done": ":green", "running": ":blue",
                  "failed": ":red", "pending": ":gray"}.get(s, ":gray")
        col.caption(f":{colour[1:]}[{s.capitalize()}]")

    st.divider()

    # ── Conda environment banner ───────────────────────────────────────────
    running = pipeline_is_running()
    active_env = active_conda_env()

    if not active_env or active_env == "base":
        st.warning(
            "**No conda environment is active** (or you are in `base`).  \n"
            "Tools installed in a conda environment will appear as missing below "
            "until that environment is activated.  \n"
            "**Fix:** close the terminal, activate your environment, then relaunch:\n"
            "```\nconda activate orthofinder3\nstreamlit run app.py\n```\n"
            "*(Replace `orthofinder3` with the name of your environment.)*",
            icon="⚠️",
        )
    else:
        st.info(f"Active conda environment: **`{active_env}`**", icon="🐍")

    # ── Pre-flight check ───────────────────────────────────────────────────
    with st.expander("🔍 Pre-flight checks", expanded=(not running)):
        tool_checks = check_tools(search_prog, gene_tree_method, msa_prog)
        all_tools_ok = all(in_path or found_env
                           for _, in_path, _, found_env in tool_checks)

        for tool, in_path, hint, found_env in tool_checks:
            if in_path:
                st.markdown(f"✅ **{tool}**")
            elif found_env:
                # Tool exists in a conda env but that env isn't active
                st.markdown(
                    f"⚠️ **{tool}** — found in conda env **`{found_env}`** "
                    f"but not on the current PATH.  \n"
                    f"&nbsp;&nbsp;&nbsp;&nbsp;Activate it first: "
                    f"`conda activate {found_env}` then relaunch the GUI."
                )
            else:
                # Not found anywhere — genuinely needs installing
                st.markdown(
                    f"❌ **{tool}** — not found in PATH or any conda environment.  \n"
                    f"&nbsp;&nbsp;&nbsp;&nbsp;Install: `{hint}`"
                )

        env_ok = True
        if species_tree:
            if _validated_tree is None:
                env_ok = False
                st.error(
                    "Species tree path is invalid. "
                    "The file must be a .tre / .tree / .nwk file inside the "
                    "pipeline working directory (no path traversal allowed)."
                )
            elif not _validated_tree.exists():
                env_ok = False
                st.error(f"Species tree file not found: `{_validated_tree.name}`")
        if not faa_files:
            env_ok = False
            st.error(f"No .faa files found in `{input_dir}`")

        if all_tools_ok and env_ok:
            st.success("All checks passed — ready to run ✓")
        elif not all_tools_ok and not active_env:
            st.info("Activate your conda environment and relaunch to re-check.")

    st.markdown("")

    # ── Run / Stop / Clear controls ────────────────────────────────────────
    c1, c2, c3 = st.columns([3, 1, 1])

    with c1:
        run_label = (
            f"▶  Run Pipeline  (steps {step_start}–{step_end})"
            if (step_start > 1 or step_end < 5)
            else "▶  Run Pipeline"
        )
        if st.button(
            run_label, type="primary",
            disabled=(running or not faa_files),
            use_container_width=True,
            help=(
                f"Runs steps {step_start}–{step_end} "
                f"({'full pipeline' if step_start == 1 and step_end == 5 else 'partial run'}). "
                "The PDF report always runs at the end. "
                "Previously completed steps are skipped automatically."
            ),
        ):
            launch_pipeline(cfg)
            time.sleep(0.6)
            st.rerun()

    with c2:
        if st.button(
            "⏹  Stop", type="secondary",
            disabled=not running,
            use_container_width=True,
            help="Terminate the running pipeline and all its subprocesses.",
        ):
            stop_pipeline()
            st.toast("Pipeline stopped.", icon="⏹")
            time.sleep(0.5)
            st.rerun()

    with c3:
        if st.button(
            "🗑  Clear log",
            disabled=running,
            use_container_width=True,
            help="Delete the pipeline log file.",
        ):
            if PIPELINE_LOG.exists():
                PIPELINE_LOG.unlink()
            st.session_state.log_offset = 0
            st.rerun()

    # ── Live log viewer ────────────────────────────────────────────────────
    st.markdown("#### 📋 Pipeline log")
    st.caption("Auto-refreshes every 3 seconds while the pipeline is running.")

    # st.fragment with run_every auto-refreshes this block independently,
    # without triggering a full page rerun (requires Streamlit ≥ 1.37).
    @st.fragment(run_every=3)
    def _live_log():
        log_text = read_log_tail()
        safe = (log_text
                .replace("&", "&amp;")
                .replace("<", "&lt;")
                .replace(">", "&gt;"))
        st.markdown(
            f'<div class="log-box" id="logbox">{safe}</div>'
            # Auto-scroll to bottom each refresh
            '<script>'
            'var b = window.parent.document.getElementById("logbox");'
            'if (b) b.scrollTop = b.scrollHeight;'
            '</script>',
            unsafe_allow_html=True,
        )
        # Show elapsed time while running
        if pipeline_is_running():
            st.caption(f"🔄 Running … last refreshed {datetime.now():%H:%M:%S}")
        else:
            step_s = get_step_status()
            done = sum(1 for s in step_s.values() if s == "done")
            if done == 6:
                st.success("🎉 Pipeline complete! Switch to the Results tab.")
            elif done > 0:
                st.info(f"{done}/6 steps complete. Press ▶ Run Pipeline to continue.")

    _live_log()


# ════════════════════════════════════════════
# Tab 2 — Results
# ════════════════════════════════════════════
with tab_results:

    r_qc, r_orth, r_aln, r_tree = st.tabs(
        ["🔬 QC", "🧩 Orthologs", "📐 Alignments", "🌲 Species Tree"]
    )

    # ── QC ─────────────────────────────────────────────────────────────────
    with r_qc:
        st.markdown("### Proteome QC Summary")
        qc_df = load_tsv("qc_results/proteome_qc_summary.tsv")
        if qc_df is not None:
            verdicts = qc_df["verdict"].value_counts() if "verdict" in qc_df.columns else {}
            m = st.columns(4)
            m[0].metric("Total genomes",  len(qc_df))
            m[1].metric("✅ PASS", int(verdicts.get("PASS", 0)))
            m[2].metric("⚠️ WARN", int(verdicts.get("WARN", 0)))
            m[3].metric("❌ FAIL", int(verdicts.get("FAIL", 0)))

            def _color_verdict(val):
                return {
                    "PASS": "background-color:#1a3a1a; color:#5dbb63",
                    "WARN": "background-color:#3a3010; color:#e3b341",
                    "FAIL": "background-color:#3a1010; color:#f85149",
                }.get(str(val), "")

            styled = qc_df.copy()
            if "verdict" in styled.columns:
                st.dataframe(
                    styled.style.applymap(_color_verdict, subset=["verdict"]),
                    use_container_width=True, hide_index=True,
                )
            else:
                st.dataframe(styled, use_container_width=True, hide_index=True)

            dl1, dl2 = st.columns(2)
            dl1.download_button(
                "⬇️ QC table (TSV)",
                qc_df.to_csv(sep="\t", index=False),
                "proteome_qc_summary.tsv", "text/tab-separated-values",
            )
            report_txt = Path("qc_results/proteome_qc_report.txt")
            if report_txt.exists():
                dl2.download_button(
                    "⬇️ QC report (TXT)",
                    report_txt.read_text(),
                    "proteome_qc_report.txt", "text/plain",
                )
        else:
            st.info("QC results not yet available. Run Step 1 first.")

    # ── Orthologs ────────────────────────────────────────────────────────────
    with r_orth:
        st.markdown("### Ortholog Selection")
        summary_df  = load_tsv("ortholog_results/tables/ortholog_summary.tsv")
        included_df = load_tsv("ortholog_results/tables/orthologs_included.tsv")

        if summary_df is not None and len(summary_df) > 0:
            row = summary_df.iloc[0]
            m = st.columns(4)
            m[0].metric("Orthogroups assessed",
                        int(row.get("total_orthogroups_assessed", 0)))
            m[1].metric("Selected for tree",
                        int(row.get("orthogroups_selected", 0)))
            m[2].metric("Species",
                        int(row.get("total_species", 0)))
            total_aa = int(row.get("total_aa_in_concat", 0))
            m[3].metric("Concatenated length",
                        f"{total_aa:,} aa" if total_aa else "—")

        if included_df is not None:
            st.markdown("**Included orthologs** (sorted by species coverage, best first)")
            if "pct_species_present" in included_df.columns:
                included_df = included_df.sort_values(
                    "pct_species_present", ascending=False
                )
            st.dataframe(included_df, use_container_width=True,
                         hide_index=True, height=320)
            st.download_button(
                "⬇️ Orthologs table (TSV)",
                included_df.to_csv(sep="\t", index=False),
                "orthologs_included.tsv", "text/tab-separated-values",
            )
        elif summary_df is None:
            st.info("Ortholog results not yet available. Run Steps 1–3 first.")

    # ── Alignments ────────────────────────────────────────────────────────────
    with r_aln:
        st.markdown("### Alignment Statistics")
        aln_df = load_tsv("alignment_results/alignment_stats.tsv")

        if aln_df is not None:
            ok_mask = (aln_df["status"] == "success") \
                if "status" in aln_df.columns else [True] * len(aln_df)

            m = st.columns(4)
            m[0].metric("Total alignments", len(aln_df))
            m[1].metric("Successful", int(sum(ok_mask)))
            if "pct_retained" in aln_df.columns:
                m[2].metric("Median % retained",
                            f"{aln_df['pct_retained'].median():.1f}%")
                m[3].metric("Min % retained",
                            f"{aln_df['pct_retained'].min():.1f}%")

            # Sort worst-first so problematic alignments surface at the top
            if "pct_retained" in aln_df.columns:
                aln_df = aln_df.sort_values("pct_retained")

            def _color_pct(val):
                try:
                    v = float(val)
                    if v < 30:  return "background-color:#3a1010; color:#f85149"
                    if v < 60:  return "background-color:#3a2a10; color:#e3b341"
                except Exception:
                    pass
                return ""

            styled = aln_df.style.applymap(_color_pct, subset=["pct_retained"]) \
                if "pct_retained" in aln_df.columns else aln_df
            st.dataframe(styled, use_container_width=True,
                         hide_index=True, height=320)
            st.caption("Rows with < 30% retained are shown in red; < 60% in yellow. "
                       "These alignments may have quality issues and are worth inspecting.")
            st.download_button(
                "⬇️ Alignment stats (TSV)",
                aln_df.to_csv(sep="\t", index=False),
                "alignment_stats.tsv", "text/tab-separated-values",
            )
        else:
            st.info("Alignment results not yet available. Run Steps 1–4 first.")

    # ── Species tree ──────────────────────────────────────────────────────────
    with r_tree:
        st.markdown("### Species Tree")
        treefile     = Path("iqtree_results/species_tree.treefile")
        iqtree_rpt   = Path("iqtree_results/species_tree.iqtree")
        pdf_path     = Path("report/pipeline_summary.pdf")

        if treefile.exists():
            tree_str    = treefile.read_text().strip()
            iqtree_text = iqtree_rpt.read_text() if iqtree_rpt.exists() else ""

            # Key stats from the IQ-TREE report
            model  = _extract_iqtree_stat(iqtree_text, r"Best-fit model.*?:\s*(\S+)")
            lnl    = _extract_iqtree_stat(iqtree_text,
                                          r"Log-likelihood of the tree:\s*([\-\d\.]+)")
            sites  = _extract_iqtree_stat(iqtree_text,
                                          r"Number of parsimony.informative sites:\s*(\d+)",
                                          re.IGNORECASE)
            n_taxa = tree_str.count(",") + 1

            m = st.columns(4)
            m[0].metric("Best-fit model", model)
            m[1].metric("Log-likelihood", lnl)
            m[2].metric("Parsimony-inf. sites", sites)
            m[3].metric("Taxa in tree", n_taxa)

            # Tree visualisation via Bio.Phylo + matplotlib
            try:
                from io import StringIO
                from Bio import Phylo
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt

                tree_obj = Phylo.read(StringIO(tree_str), "newick")
                n_leaves = n_taxa
                fig_h = max(5, n_leaves * 0.28)
                fig, ax = plt.subplots(figsize=(11, fig_h))
                fig.patch.set_facecolor("#0e1117")
                ax.set_facecolor("#0e1117")
                Phylo.draw(tree_obj, axes=ax, do_show=False)
                ax.set_title("Species tree — IQ-TREE3",
                             fontsize=11, color="#c9d1d9")
                for spine in ax.spines.values():
                    spine.set_edgecolor("#30363d")
                ax.tick_params(colors="#8b949e")
                ax.yaxis.label.set_color("#8b949e")
                ax.xaxis.label.set_color("#8b949e")
                plt.tight_layout()
                st.pyplot(fig, use_container_width=True)
                plt.close(fig)
            except Exception as e:
                st.info(
                    f"Graphical tree preview unavailable ({type(e).__name__}: {e}). "
                    "Download the .treefile and open it in FigTree or iTOL."
                )

            st.markdown("**Newick string:**")
            st.code(tree_str, language=None)

            # Branch support note
            st.caption(
                "Branch support: SH-aLRT / UFBoot values on each node. "
                "Well-supported branches: **SH-aLRT ≥ 80%** AND **UFBoot ≥ 95%**."
            )

            # Downloads
            dl1, dl2, dl3 = st.columns(3)
            dl1.download_button(
                "⬇️ Species tree (.treefile)",
                tree_str, "species_tree.treefile", "text/plain",
            )
            dl2.download_button(
                "⬇️ Species tree (.nwk)",
                tree_str, "species_tree.nwk", "text/plain",
            )
            if pdf_path.exists():
                dl3.download_button(
                    "⬇️ PDF report",
                    pdf_path.read_bytes(),
                    "pipeline_summary.pdf", "application/pdf",
                )

            # Full IQ-TREE report
            if iqtree_text:
                with st.expander("📄 Full IQ-TREE3 report"):
                    st.text(iqtree_text)
        else:
            st.info("Species tree not yet available. Run the full pipeline first.")


# ════════════════════════════════════════════
# Tab 3 — Help
# ════════════════════════════════════════════
with tab_help:
    col_a, col_b = st.columns([3, 2])

    with col_a:
        st.markdown("""
## Quick start

1. **Set the proteomes directory** in the left sidebar — it should contain one `.faa` file per genome.
2. **Adjust settings** as needed. Sensible defaults are pre-filled.
3. Click **🔍 Pre-flight checks** to confirm all required tools are installed.
4. Click **▶ Run Pipeline**.
5. Watch the live log on the **Run Pipeline** tab.
6. Switch to the **📊 Results** tab at any time — each sub-tab populates as its step completes.

---

## Settings reference

| Setting | Default | Notes |
|---|---|---|
| Search method | `diamond` | Fast; `diamond_ultra_sens` for sensitive searches on divergent genomes |
| Gene tree method | `msa` | OrthoFinder3 default; more accurate than `dendroblast` |
| MSA program | `mafft` | Only applies when gene tree method = MSA |
| Species tree | auto-detected | Rooted Newick file (.tre/.tree) to guide OrthoFinder. Leave blank to infer. |
| Threads | 8 | **Leave cores free on a shared server** |
| Ortholog presence | 1.0 | Lower to 0.95 if fewer than 50 orthologs are found |
| trimAl mode | `auto` | `-automated1`; rarely needs changing |
| IQ-TREE mode | `partition` | Per-gene model (`-p`); better but slower than `single` (`-s`) |
| Model | `LG+F+G4` | Applies to **both** modes. `MFP` = ModelFinder (slow). `MFP+MERGE` = slowest |
| Bootstrap | 1000 | Minimum for publication; rarely need more |
| SH-aLRT | 1000 | Set to 0 to disable |
| Start from step | 1 | Skip early steps when inputs already exist |
| Stop after step | 5 | Stop before IQ-TREE if you only want orthologs, etc. |

---

## Step descriptions

| Step | What it does | Needs | Produces |
|------|-------------|-------|---------|
| 1 · QC | Checks protein count, length, duplicates; flags/excludes low-quality genomes | `proteomes/*.faa` | `qc_results/passed_proteomes/` |
| 2 · OrthoFinder | DIAMOND all-vs-all search; clusters proteins into orthogroups | passed proteomes | `orthofinder_results/Results_*/` |
| 3 · Extract | Selects single-copy orthologs; writes one .faa per gene with species-name headers | OrthoFinder results | `ortholog_results/gene_fastas/` |
| 4 · Align & trim | MAFFT alignment + trimAl column trimming per gene | gene FASTAs | `alignment_results/trimmed/` |
| 5 · IQ-TREE | Infers species phylogeny from trimmed alignments | trimmed alignments | `iqtree_results/species_tree.treefile` |
| 6 · Report | Multi-page PDF summary (always runs at end) | any available outputs | `report/pipeline_summary.pdf` |

---

## Re-running after a failure

Completed steps are **skipped automatically** on re-run. To force a step to re-run, delete its marker file from a terminal:

```bash
rm qc_results/.step_complete          # re-run QC
rm -rf orthofinder_results/Results_*  # re-run OrthoFinder
rm ortholog_results/.step_complete    # re-run ortholog extraction
rm alignment_results/.step_complete   # re-run alignment & trimming
rm iqtree_results/.step_complete      # re-run IQ-TREE
```

---

## Installing required tools

```bash
conda install -c bioconda orthofinder diamond mafft trimal iqtree
conda install -c conda-forge streamlit psutil
```

---

## Visualising the species tree

Open `iqtree_results/species_tree.treefile` in:
- **[FigTree](https://tree.bio.ed.ac.uk/software/figtree/)** — free desktop viewer; shows bootstrap support values on nodes
- **[iTOL](https://itol.embl.de/)** — web-based; drag-and-drop the treefile

Branch support values appear as `SH-aLRT/UFBoot` on each node. A branch is
considered well-supported when **SH-aLRT ≥ 80 %** AND **UFBoot ≥ 95 %**.
        """)

    with col_b:
        st.markdown("""
## Ortholog presence threshold

| Value | Behaviour |
|---|---|
| **1.0** | Gene must be in **all** species (strict, default) |
| **0.99** | One missing species allowed per gene |
| **0.95** | Gene in ≥ 95% of species |
| **0.90** | Gene in ≥ 90% of species |

Start with 1.0. If fewer than 50 orthologs pass, try 0.95.

---

## Pipeline outputs

```
qc_results/
  proteome_qc_summary.tsv
  passed_proteomes/

orthofinder_results/
  Results_[date]/

ortholog_results/
  gene_fastas/OG*.faa
  tables/

alignment_results/
  aligned/   trimmed/

iqtree_results/
  species_tree.treefile  ← main result
  species_tree.iqtree

report/
  pipeline_summary.pdf

pipeline_run.log        ← full run log
```

---

## Citation

If you publish results from this pipeline, please cite:

- **OrthoFinder**: Emms & Kelly (2019) *Genome Biol* 20:238
- **DIAMOND**: Buchfink et al. (2015) *Nat Methods* 12:59
- **MAFFT**: Katoh & Standley (2013) *Mol Biol Evol* 30:772
- **trimAl**: Capella-Gutiérrez et al. (2009) *Bioinformatics* 25:1972
- **IQ-TREE**: Minh et al. (2020) *Mol Biol Evol* 37:1530
        """)
