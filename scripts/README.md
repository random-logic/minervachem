# Demo runners

Keep the current Buhito checkout at `../buhito`. MinervaChem's `uv` configuration
uses it as an editable dependency. Both repositories and their environments must
be accessible at the same paths on the HPC compute nodes.

From the MinervaChem repository root:

```sh
./scripts/run_demo.sh --list
./scripts/run_all_demos.sh
./scripts/run_demo.sh --demo demos/1_Obtaining_a_fingerprint.ipynb
./scripts/run_demo.sh --demo demos/qm9/benchmark_graphlet_featurizers_train_test.py -- --data demos/qm9/qm9data/qm9_processed.csv --max-len 3 --n-jobs 8
./scripts/run_demo.sh --demo demos/process_qm9.py
./scripts/run_demo.sh --demo demos/mcts/create_pipeline.py -- --limit 100 --output qm9_pipeline.joblib
```

`run_all_demos.sh` discovers every Python demo and notebook under `demos/`,
recursively. It excludes checkpoints, hidden directories, and old executed copies.
The shell launchers submit one independent Slurm job per selected demo by default.
`--slurm` remains accepted for compatibility. They require `sbatch`; submission
failure never falls back to local execution. Python demos and notebooks both
run through `demo_job.sh` on compute nodes. The worker refuses to execute without
a Slurm allocation. The frontend runs the submission script through `uv run`;
only the batch worker executes demo code. Scripts run from the repository root; notebooks run from
their own directory so relative data paths continue to work. Python arguments
following `--` are forwarded unchanged. All-demo runs use each demo's default
inputs, so supply the required datasets first. This is not a dependency scheduler:
run preprocessing before demos that need its output. Slurm submission failures
stop further submissions; already submitted jobs remain queued.

Every notebook is executed with `nbconvert --execute --inplace`: successful results
replace its original `.ipynb`, with no `_executed`/`.executed` copy. A failed
execution leaves the original notebook unchanged. Logs go to `logs/`.

Notebook `HPC_RESOURCES` comments in the first Markdown cell set CPU/memory
requests. The shared runner validates every notebook marker before submitting
any jobs; missing or invalid markers are errors. Joblib workers are capped at
`SLURM_CPUS_PER_TASK` and native library threads are limited to one per worker.
Python demos can set a `# HPC_RESOURCES: cpus=8 mem=32G` comment; without one,
they default to 1 CPU and 8G RAM. Override with `--cpus`, `--mem`, `--time`, `--partition`,
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

## Demo resource requests

| Demo | CPUs | RAM | Workload |
| --- | ---: | ---: | --- |
| `1_Obtaining_a_fingerprint.ipynb` | 2 | 4G | Individual molecule fingerprints |
| `2_Regression_and_visualization.ipynb` | 64 | 64G | Full QM9 featurization and regression |
| `3_Advanced_fitting_and_unseen_bits.ipynb` | 64 | 64G | Dataset featurization and model fitting |
| `5_Filtering_substructures_in_fp.ipynb` | 64 | 64G | Multiple QM9 feature/filtering runs |
| `mcts/1_MCTS_with_RDKit.ipynb` | 2 | 8G | Descriptor-based molecular search |
| `mcts/2_MCTS_with_MinervaChem.ipynb` | 64 | 64G | Molecular model training and search |
| `mcts/3_MCTS_experiments.ipynb` | 2 | 8G | Search experiments |
| `mcts/create_pipeline.py` | 8 | 32G | Full QM9 model training |
| `process_qm9.py` | 1 | 8G | Serial raw QM9 preprocessing |
| `qm9/benchmark_graphlet_featurizers_train_test.py` | 8 | 32G | QM9 graph conversion and BFS/DFS benchmarking |

These are per-job allocations. Regression, advanced fitting, and filtering use
`n_jobs=-3`, which starts up to 62 joblib workers in a 64-CPU allocation. Pipeline
training and the QM9 benchmark use the Slurm CPU allocation by default. Serial
stages do not keep every allocated CPU busy.

Local sampled workload checks on 2026-11-21 support keeping the notebook
allocations and raising the two parallel Python demos from 1 CPU/8G to 8 CPUs/32G.
See [the smoke-run report](demo_resource_smoke.md) for measurements and limits.
Full-dataset peak usage and optimal requests still require measurements on the
target HPC hardware.

The batch worker follows the original notebook runner: one node and task,
`standard` QoS, account `w26_pfas_g`, no GPUs, and no requeue. Use `--account`
or `HPC_ACCOUNT` to override the account. Runtime and cache directories use
`SLURM_TMPDIR` (or `/tmp`) inside each job. The frontend requires `uv` and `sbatch`. It runs the submission script through
`uv run` using the project Python, without importing or executing demo code.
