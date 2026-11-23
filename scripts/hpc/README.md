# Transfer MinervaChem data and logs

## Download logs

Run locally from this checkout (or invoke the script by its full path):

```bash
# Download missing logs from se:~/minervachem/logs into this checkout's logs/.
bash scripts/hpc/download_logs.sh

# Refresh existing logs, including logs that have grown since the last download.
bash scripts/hpc/download_logs.sh --redownload

# Preview transfers; this still connects to HPC for the file listing.
bash scripts/hpc/download_logs.sh --dry-run

# Override the SSH alias or remote directory.
bash scripts/hpc/download_logs.sh --remote se --remote-dir /users/ajowe/minervachem/logs
```

Requires rsync locally and remotely. Existing local files are skipped by default;
`--redownload` overwrites them. Local-only logs are retained and remote files are
never deleted. WSL automatically uses Windows `ssh.exe` for smart-card login;
other systems use `ssh`. Authentication happens interactively when you run it.

## Upload input data

From the current MinervaChem repository root:

```bash
# Inspect local inputs only; no connection.
bash scripts/hpc/upload_data.sh --list

# Compare against HPC without uploading (remote directories must exist).
bash scripts/hpc/upload_data.sh --dry-run

# Upload missing files to ~/minervachem on se.
bash scripts/hpc/upload_data.sh
```

The uploader prefers input files from this checkout and supplements them from
`~/minervachem` (override with `--upstream PATH`). It includes:

- Processed QM9 at `demos/qm9_processed.csv`.
- The benchmark QM9 CSV at `demos/qm9/qm9data/qm9_processed.csv`.
- Raw QM9 files under `demos/qm9_data_files/`.
- PyTDC data under `demos/data/`. Upstream root `data/` is mapped here so the
  notebooks find the cache without downloading on compute nodes.
- Additional upstream input directories `Datasets/`, `Datasets_no_split/`, and
  `Datasets_Agg/`, preserving their paths under the remote repository.

No environments, source code, credentials, generated result directories, or
notebook outputs are uploaded. Existing remote files are always skipped using
`rsync --ignore-existing`; remote-only files are retained. Missing required local
inputs are reported and give a nonzero exit status after available uploads.

All implementation lives in `scripts/`: this wrapper uses the shared helper at
`../buhito/scripts/hpc/upload_data.py`. Keep both checkouts beside one another.
Options include `--remote ch`, `--remote-dir minervachem`, `--source PATH`,
`--upstream PATH`, `--ssh-command CMD`, `--dry-run`, and `--list`.

On WSL, Windows OpenSSH is selected automatically so aliases and smart-card
settings in the Windows SSH config are honored. Run authentication yourself;
never put a PIN or password in these scripts. No remote action runs automatically.
