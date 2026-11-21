# Upload MinervaChem input data to HPC

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
