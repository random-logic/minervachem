# Demo runners

Keep the current Buhito checkout at `../buhito`. MinervaChem's `uv` configuration
uses it as an editable dependency. Both repositories and their environments must
be accessible at the same paths on the HPC compute nodes.

From the MinervaChem repository root:

```sh
./scripts/run_demo.sh --list
./scripts/run_all_demos.sh --slurm
./scripts/run_demo.sh --demo demos/1_Obtaining_a_fingerprint.ipynb --slurm
./scripts/run_demo.sh --demo demos/qm9/benchmark_graphlet_featurizers_train_test.py --slurm -- --data demos/qm9/qm9data/qm9_processed.csv --max-len 3 --n-jobs 8
./scripts/run_demo.sh --demo demos/process_qm9.py --slurm
./scripts/run_demo.sh --demo demos/mcts/create_pipeline.py --slurm -- --limit 100 --output qm9_pipeline.joblib
```

`run_all_demos.sh` discovers every Python demo and notebook under `demos/`,
recursively. It excludes checkpoints, hidden directories, and old executed copies.
With `--slurm` it submits one independent job per demo; without it demos run
locally in sequence. Scripts run from the repository root; notebooks run from
their own directory so relative data paths continue to work. Python arguments
following `--` are forwarded unchanged. All-demo runs use each demo's default
inputs, so supply the required datasets first. This is not a dependency scheduler:
run preprocessing before demos that need its output. Local failures are collected
and cause a nonzero exit after the remaining demos run. Slurm submission failures
stop further submissions; already submitted jobs remain queued.

Every notebook is executed with `nbconvert --execute --inplace`: successful results
replace its original `.ipynb`, with no `_executed`/`.executed` copy. A failed
execution leaves the original notebook unchanged. Logs go to `logs/`.

Notebook `HPC_RESOURCES` comments set CPU/memory requests. Other demos default to
1 CPU and 8G RAM. Override with `--cpus`, `--mem`, `--time`, `--partition`,
`--account`, or `--mail-user`; `HPC_PARTITION`, `HPC_ACCOUNT`, `HPC_TIME`, and
`MAIL_USER` are also supported. The scripts load `.env` when present.

The original notebook-only commands remain available:

```sh
./scripts/run_demo_notebooks.sh
./scripts/run_demo_notebook.sh --notebook demos/2_Regression_and_visualization.ipynb
```

These automatically submit Slurm jobs, use the existing cluster directives and
notebook resource markers, and require `MAIL_USER` in `.env`. They also update
notebooks in place. Installing or publishing the package does not submit jobs;
run a script explicitly.

| Demos | Required inputs |
| --- | --- |
| `1_Obtaining_a_fingerprint.ipynb` | Molecules defined inside the notebook |
| `2_Regression_and_visualization.ipynb` | QM9 processed data in `demos/` |
| `3_Advanced_fitting_and_unseen_bits.ipynb` | Datasets described in the notebook; network access may be needed |
| `5_Filtering_substructures_in_fp.ipynb` | Chemical examples and data described in the notebook |
| `mcts/*.ipynb` | Molecular search inputs and models described in each notebook |
| `process_qm9.py` | Raw QM9 archives and reference files in `demos/qm9_data_files/` |
| `qm9/benchmark_graphlet_featurizers_train_test.py` | QM9 CSV supplied with `--data` |
| `mcts/create_pipeline.py` | `demos/qm9_processed.csv`; writes a fitted pipeline |

The synthetic PADRE notebook and Reddit demos now live in Buhito and have runners
in Buhito's `scripts/` directory.
