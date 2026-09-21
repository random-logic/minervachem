#!/usr/bin/env bash

# Submit exactly one demo notebook to Slurm.
# Usage:
#   ./scripts/run_demo_notebook.sh --notebook demos/2_Regression_and_visualization.ipynb
#   ./scripts/run_demo_notebook.sh -n demos/mcts/2_MCTS_with_MinervaChem.ipynb

set -euo pipefail

usage() {
    echo "Usage: $0 --notebook PATH" >&2
    exit 2
}

NOTEBOOK_ARG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--notebook)
            [[ $# -ge 2 ]] || usage
            NOTEBOOK_ARG="$2"
            shift 2
            ;;
        --notebook=*)
            NOTEBOOK_ARG="${1#--notebook=}"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage
            ;;
    esac
done

[[ -n "$NOTEBOOK_ARG" ]] || usage

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
WORKER_SCRIPT="$SCRIPT_DIR/run_demo_notebooks.sh"
LOG_DIR="$PROJECT_DIR/logs"

case "$NOTEBOOK_ARG" in
    /*) NOTEBOOK="$NOTEBOOK_ARG" ;;
    *) NOTEBOOK="$PROJECT_DIR/$NOTEBOOK_ARG" ;;
esac

if [[ ! -f "$NOTEBOOK" ]]; then
    echo "ERROR: notebook does not exist: $NOTEBOOK" >&2
    exit 1
fi

case "$NOTEBOOK" in
    "$PROJECT_DIR"/demos/*.ipynb) ;;
    *)
        echo "ERROR: notebook must be an .ipynb file under $PROJECT_DIR/demos" >&2
        exit 1
        ;;
esac

if [[ ! -f "$PROJECT_DIR/.env" ]]; then
    echo "ERROR: missing $PROJECT_DIR/.env" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$PROJECT_DIR/.env"
if [[ -z "${MAIL_USER:-}" ]]; then
    echo "ERROR: MAIL_USER is not set in $PROJECT_DIR/.env" >&2
    exit 1
fi

if ! command -v uv >/dev/null 2>&1; then
    echo "ERROR: uv is not available for job submission." >&2
    exit 1
fi

resources="$(uv run --project "$PROJECT_DIR" --all-extras python -c '
import json
import re
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    notebook = json.load(handle)

first_markdown = next(
    (
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "markdown"
    ),
    "",
)
match = re.search(
    r"<!--\s*HPC_RESOURCES:\s*cpus=(\d+)\s+mem=(\d+[KMGTP])\s*-->",
    first_markdown,
)
if match is None:
    sys.exit("missing HPC_RESOURCES marker in the first Markdown cell")

print(f"{match.group(1)} {match.group(2)}")
' "$NOTEBOOK")"
read -r cores memory <<< "$resources"

if [[ ! "$cores" =~ ^[1-9][0-9]*$ ]] || [[ ! "$memory" =~ ^[1-9][0-9]*[KMGTP]$ ]]; then
    echo "ERROR: invalid HPC_RESOURCES marker in $NOTEBOOK: $resources" >&2
    exit 1
fi

mkdir -p "$LOG_DIR"
submission="$(sbatch --parsable \
    --chdir="$PROJECT_DIR" \
    --mail-user="$MAIL_USER" \
    --mail-type=BEGIN,END,FAIL \
    --cpus-per-task="$cores" \
    --mem="$memory" \
    --export=ALL,MINERVACHEM_PROJECT_DIR="$PROJECT_DIR",MINERVACHEM_NOTEBOOK="$NOTEBOOK" \
    "$WORKER_SCRIPT")"

echo "Submitted ${NOTEBOOK#"$PROJECT_DIR"/} as job ${submission%%;*} ($cores CPUs, $memory RAM)."
