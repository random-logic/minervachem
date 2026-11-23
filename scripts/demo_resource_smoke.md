# Demo resource smoke runs — 2026-11-21

All 15 sampled compute runs passed. These ran locally in MinervaChem's existing
Python 3.11 environment on a 24-CPU host. The Slurm controller is unavailable from
this workspace, so these are not cluster jobs or measurements of full production
notebooks. The original notebooks and datasets were left unchanged.

The samples use real local QM9 data (130,831 molecules) and the cached
Lipophilicity AstraZeneca dataset (4,200 molecules), sampled with random seed 42.
QM9 workloads use explicit hydrogens; advanced fitting uses max_len=7 and
useHs=True, matching the notebook. Training/test splits use 80%/20%.

The workload harness used the same featurizers and regressors as the demos:

- Regression: max_len=5 fingerprints, train/test transforms, Ridge, and the
  hierarchical residual model.
- Filtering: both O/N and carbon filtering passes, retaining both feature sets,
  plus Ridge fits and predictions.
- Advanced fitting: max_len=7 fingerprints, unseen-feature transformation,
  Ridge, and hierarchical cross-validation. The smoke grid has two alphas and
  two folds, compared with 30 alphas and five folds in the full notebook.
- Pipeline/MCTS: the actual create_pipeline function followed by two search
  levels with budget 20 and beam size 3 (full notebook: five levels, budget
  10,000, beam size 10). Training and search share the measured process tree.
- RDKit MCTS and MCTS experiments: reduced LogP search only; the experiments
  notebook's full diversity and pairwise-similarity stages were not profiled.
- Benchmark: the actual benchmark main function, both BFS and DFS, max_len=3,
  train/test transforms, one repetition and no warmup (default full run: three
  repetitions). Graph conversion and result writing are included.

Model predictions were checked for finite values. Worker count was explicit;
LOKY_MAX_CPU_COUNT matched it and native numerical-library threads were limited
to one. Peak memory sums RSS over the parent and all descendants, sampled every
0.1 seconds. This includes workers and conservatively double-counts shared
pages; short peaks can still be missed. Wall times include process startup.
Plots, full notebook rendering, and raw QM9 preprocessing were not measured.

| Workload | Molecules | Workers | Wall seconds | Peak process-tree RSS (GiB) |
| --- | ---: | ---: | ---: | ---: |
| regression | 512 | 4 | 2.93 | 0.804 |
| regression | 2,048 | 4 | 7.06 | 0.923 |
| filtering | 512 | 4 | 4.03 | 0.821 |
| filtering | 2,048 | 4 | 12.20 | 0.991 |
| advanced | 32 | 4 | 3.23 | 0.923 |
| advanced | 128 | 4 | 9.68 | 0.992 |
| pipeline_mcts | 512 | 4 | 2.11 | 0.877 |
| mcts_rdkit | 128 | 1 | 1.41 | 0.284 |
| mcts_experiments | 128 | 1 | 1.31 | 0.284 |
| benchmark | 512 | 4 | 1.91 | 0.645 |
| pipeline_mcts | 8,192 | 4 | 4.13 | 1.078 |
| pipeline_mcts | 8,192 | 8 | 3.73 | 1.634 |
| benchmark | 8,192 | 4 | 9.76 | 0.967 |
| benchmark | 8,192 | 8 | 9.27 | 1.348 |
| advanced | 512 | 4 | 33.77 | 1.274 |

## Configuration decisions

- Keep the four heavy notebooks at 64 CPUs/64G. These samples provide no evidence
  that larger allocations are required. They do not certify full-run peak usage.
  Regression and filtering showed modest memory growth between 512 and 2,048
  molecules; advanced fitting stayed at 1.274 GiB on 512 molecules with four workers.
  More workers, full cross-validation, plots, and retained state add overhead.
- Raise create_pipeline.py and the QM9 benchmark from the fallback 1 CPU/8G to
  explicit 8 CPUs/32G. Their 8,192-molecule runs passed with eight workers, with
  modest speed improvements over four workers and peaks below 1.7 GiB. The
  full QM9 dataset is about 16 times larger; 32G provides additional headroom
  for graph objects, retained fingerprints, and worker overhead. Larger custom
  --max-len settings are not covered by the default benchmark smoke test.
- Keep the smaller notebook allocations and serial preprocessing allocation.
  Extra cores will not parallelize serial search or preprocessing stages.

The resource comments control actual sbatch requests. A separate mocked submission
check through uv run verified all ten discovered demos and the new 8-CPU/32G
requests. Six MinervaChem repository-script tests passed. No real cluster jobs
were submitted by these checks.
