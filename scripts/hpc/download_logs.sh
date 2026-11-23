#!/usr/bin/env bash
set -euo pipefail

REMOTE="se"
REMOTE_DIR="minervachem/logs"
REDOWNLOAD=false
DRY_RUN=false

usage() {
    cat <<'EOF'
Download HPC logs into this MinervaChem checkout's logs/ directory.

Usage: bash scripts/hpc/download_logs.sh [OPTIONS]

  --redownload       Transfer all files again, overwriting local copies.
  --remote HOST      SSH host or alias (default: se).
  --remote-dir PATH  Remote log directory (default: minervachem/logs,
                     relative to the remote home; absolute paths also work).
  --dry-run          Show planned transfers without writing files; connects to HPC.
  -h, --help         Show this help.

Requires rsync locally and on HPC. Existing local files are skipped even if
the remote logs have grown; use --redownload to refresh them. Local-only files
are retained. On WSL, Windows ssh.exe is used for smart-card authentication.
EOF
}

while (($#)); do
    case "$1" in
        --redownload) REDOWNLOAD=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        --remote|--remote-dir)
            if (($# < 2)) || [[ -z "$2" || "$2" == --* ]]; then
                echo "Error: $1 requires a value." >&2
                exit 2
            fi
            case "$1" in
                --remote) REMOTE="$2" ;;
                --remote-dir) REMOTE_DIR="$2" ;;
            esac
            shift 2
            ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Error: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
DEST="$PROJECT_ROOT/logs"

if grep -qi microsoft /proc/sys/kernel/osrelease 2>/dev/null; then
    SSH_CMD="ssh.exe"
    echo "WSL detected: using $SSH_CMD"
else
    SSH_CMD="ssh"
fi

command -v "$SSH_CMD" >/dev/null || {
    echo "Error: $SSH_CMD not found." >&2
    exit 1
}
command -v rsync >/dev/null || {
    echo "Error: rsync is required locally and on the HPC host." >&2
    exit 1
}

# Blocking I/O avoids EAGAIN writes with Windows ssh.exe under WSL.
OPTIONS=(--archive --omit-dir-times --human-readable --info=progress2
         --protect-args --blocking-io)
if $REDOWNLOAD; then
    OPTIONS+=(--ignore-times --whole-file)
    echo "Mode: redownload all remote files and overwrite local copies."
else
    OPTIONS+=(--ignore-existing)
    echo "Mode: download missing files only; keep existing local files."
fi
if $DRY_RUN; then
    OPTIONS+=(--dry-run --itemize-changes)
else
    mkdir -p "$DEST"
fi

echo "Source:      $REMOTE:${REMOTE_DIR%/}/"
echo "Destination: $DEST/"
echo "Connecting to $REMOTE; complete SSH/smart-card authentication if prompted."
rsync "${OPTIONS[@]}" -e "$SSH_CMD -o BatchMode=no" \
    -- "$REMOTE:${REMOTE_DIR%/}/" "$DEST/"

if $DRY_RUN; then
    echo "Dry run completed successfully; no files transferred."
else
    echo "All transfers completed successfully."
fi
