#!/usr/bin/env bash

#SBATCH --job-name=minervachem-demos
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --qos=standard
#SBATCH --account=w26_pfas_g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --gpus-per-node=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --no-requeue
#SBATCH --output=logs/demo-notebooks-%j.out
#SBATCH --error=logs/demo-notebooks-%j.err

set -euo pipefail

# When run directly, read .env and submit this same script as a Slurm job.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
    JOB_SCRIPT="$SCRIPT_DIR/run_demo_notebooks.sh"
    LOG_DIR="$PROJECT_DIR/logs"

    # Slurm must be able to open these paths before the batch script starts.
    mkdir -p "$LOG_DIR"

    if [[ -f "$PROJECT_DIR/.env" ]]; then
        # shellcheck disable=SC1091
        source "$PROJECT_DIR/.env"
    fi

    if [[ -z "${MAIL_USER:-}" ]]; then
        echo "ERROR: MAIL_USER is not set in $PROJECT_DIR/.env" >&2
        exit 1
    fi

    echo "Submitting from project directory: $PROJECT_DIR"
    exec sbatch \
        --chdir="$PROJECT_DIR" \
        --mail-user="$MAIL_USER" \
        --export=ALL,MINERVACHEM_PROJECT_DIR="$PROJECT_DIR" \
        "$JOB_SCRIPT" "$@"
fi

# The project path is explicitly exported by the self-submission above.
if [[ -z "${MINERVACHEM_PROJECT_DIR:-}" ]]; then
    echo "ERROR: MINERVACHEM_PROJECT_DIR is missing." >&2
    echo "Run this script directly: ./scripts/run_demo_notebooks.sh" >&2
    exit 1
fi

PROJECT_DIR="$MINERVACHEM_PROJECT_DIR"
NOTEBOOK_DIR="$PROJECT_DIR/demos"
LOG_DIR="$PROJECT_DIR/logs"
UV_RUN=(uv run --project "$PROJECT_DIR")

mkdir -p "$LOG_DIR"

echo "Project directory: $PROJECT_DIR"
echo "Notebook directory: $NOTEBOOK_DIR"

if [[ ! -d "$NOTEBOOK_DIR" ]]; then
    echo "ERROR: notebook directory does not exist: $NOTEBOOK_DIR" >&2
    exit 1
fi

if ! command -v uv >/dev/null 2>&1; then
    echo "ERROR: uv is not available on the compute node." >&2
    exit 1
fi

# uv run automatically synchronizes the project environment and selects it for
# nbconvert and every spawned Python kernel.
export PYTHONNOUSERSITE=1

# Keep runtime and temporary files on compute-node-local storage.
RUNTIME_BASE="${SLURM_TMPDIR:-/tmp}/minervachem-${SLURM_JOB_ID}"
export JUPYTER_RUNTIME_DIR="$RUNTIME_BASE/jupyter"
export IPYTHONDIR="$RUNTIME_BASE/ipython"
export MPLCONFIGDIR="$RUNTIME_BASE/matplotlib"
export MPLBACKEND=Agg
export XDG_CACHE_HOME="$RUNTIME_BASE/cache"
export JOBLIB_TEMP_FOLDER="$RUNTIME_BASE/joblib"
mkdir -p \
    "$JUPYTER_RUNTIME_DIR" \
    "$IPYTHONDIR" \
    "$MPLCONFIGDIR" \
    "$XDG_CACHE_HOME" \
    "$JOBLIB_TEMP_FOLDER"

# Fail once with a useful diagnostic instead of failing every notebook.
if ! "${UV_RUN[@]}" python -c 'import ipykernel, nbconvert, minervachem' \
    2>"$LOG_DIR/environment-check.log"; then
    echo "ERROR: required Python packages cannot be imported on the compute node." >&2
    echo "See $LOG_DIR/environment-check.log" >&2
    exit 1
fi

mapfile -t NOTEBOOKS < <(
    find "$NOTEBOOK_DIR" \
        -type f \
        -name '*.ipynb' \
        -not -path '*/.ipynb_checkpoints/*' \
        | sort
)

if [[ ${#NOTEBOOKS[@]} -eq 0 ]]; then
    echo "No notebooks found under $NOTEBOOK_DIR" >&2
    exit 1
fi

# Inspect code cells only. This avoids classifying a notebook from stale output
# text and avoids requiring ripgrep on the compute node.
notebook_uses_parallel_jobs() {
    "${UV_RUN[@]}" python -c '
import json
import re
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    notebook = json.load(handle)

source = "\n".join(
    "".join(cell.get("source", []))
    for cell in notebook.get("cells", [])
    if cell.get("cell_type") == "code"
)

# n_jobs=1 is serial. Negative values and values greater than one are parallel.
parallel = re.search(r"\bn_jobs\s*=\s*(?!1\b)", source) is not None
sys.exit(0 if parallel else 1)
' "$1"
}

failures=0

for notebook in "${NOTEBOOKS[@]}"; do
    relative="${notebook#"$PROJECT_DIR"/}"
    safe_name="${relative////__}"
    log_file="$LOG_DIR/${safe_name%.ipynb}.log"

    if notebook_uses_parallel_jobs "$notebook"; then
        cores=64
    else
        cores=1
    fi

    echo "Running $relative with $cores CPU(s)"

    # n_jobs controls joblib workers. Keep each worker's native libraries at
    # one thread to prevent process_count x thread_count oversubscription.
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

    if srun --exclusive \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="$cores" \
        --cpu-bind=cores \
        "${UV_RUN[@]}" jupyter nbconvert \
            --to notebook \
            --execute "$notebook" \
            --inplace \
            --ExecutePreprocessor.kernel_name=python3 \
            --ExecutePreprocessor.timeout=-1 \
            >"$log_file" 2>&1; then
        echo "SUCCESS: $relative"
    else
        echo "FAILED:  $relative (see $log_file)" >&2
        failures=$((failures + 1))
    fi
done

echo "Completed ${#NOTEBOOKS[@]} notebook(s); failures: $failures"
exit "$failures"
