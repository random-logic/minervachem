#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
HELPER="$PROJECT_ROOT/../buhito/scripts/hpc/upload_data.py"
if [[ ! -f "$HELPER" ]]; then
    echo "Missing shared uploader: $HELPER. Keep the Buhito checkout beside MinervaChem." >&2
    exit 1
fi
exec python3 "$HELPER" --repo minervachem --source "$PROJECT_ROOT" "$@"
