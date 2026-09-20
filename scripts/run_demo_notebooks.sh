#!/usr/bin/env bash
#
# Execute every Jupyter notebook under demos/ on a Slurm HPC cluster.
# Notebooks containing n_jobs are allocated 64 CPUs;
# all other notebooks are allocated 1 CPU.
#
# Submit with:
#   sbatch scripts/run_demo_notebooks.sh
#

#SBATCH --job-name=minervachem-demos
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --qos=standard
#SBATCH --account=w26_pfas_g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-task=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --no-requeue
#SBATCH --output=logs/demo-notebooks-%j.out
#SBATCH --error=logs/demo-notebooks-%j.err

set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
NOTEBOOK_DIR="$PROJECT_DIR/demos"
LOG_DIR="$PROJECT_DIR/logs"

if [[ -f "$PROJECT_DIR/.env" ]]; then
    # shellcheck disable=SC1091
    source "$PROJECT_DIR/.env"
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    if [[ -z "${MAIL_USER:-}" ]]; then
        echo "ERROR: MAIL_USER is not set in $PROJECT_DIR/.env" >&2
        exit 1
    fi

    echo "Submitting Slurm job with mail notifications for $MAIL_USER"
    exec sbatch --chdir="$PROJECT_DIR" --mail-user="$MAIL_USER" "$0" "$@"
fi

# Activate the environment used by the notebooks, if needed. For example:
# source "$PROJECT_DIR/.venv/bin/activate"
# or:
# source "$HOME/miniconda3/etc/profile.d/conda.sh"
# conda activate minervachem

mkdir -p "$LOG_DIR"

if [[ ! -d "$NOTEBOOK_DIR" ]]; then
    echo "ERROR: notebook directory does not exist: $NOTEBOOK_DIR" >&2
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

failures=0

for notebook in "${NOTEBOOKS[@]}"; do
    relative="${notebook#"$PROJECT_DIR"/}"
    safe_name="${relative//\//__}"
    log_file="$LOG_DIR/${safe_name%.ipynb}.log"

    # Allocate extra CPUs only when the notebook actually sets n_jobs.
    # Merely using a GraphletFingerprinter does not imply parallel execution.
    if rg -q -i 'n_jobs[[:space:]]*=' "$notebook"; then
        cores=64
    else
        cores=1
    fi

    echo "Running $relative with $cores CPU(s)"

    # Limit common numerical libraries to the CPUs assigned to this notebook.
    export OMP_NUM_THREADS="$cores"
    export MKL_NUM_THREADS="$cores"
    export OPENBLAS_NUM_THREADS="$cores"
    export NUMEXPR_NUM_THREADS="$cores"

    if srun --exclusive \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="$cores" \
        --cpu-bind=cores \
        jupyter nbconvert \
            --to notebook \
            --execute "$notebook" \
            --inplace \
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
