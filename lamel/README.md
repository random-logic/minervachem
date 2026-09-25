# MinervaChem LAMEL

LAMEL is an optional meta-learning add-on for
[MinervaChem](https://github.com/lanl/minervachem).

Install published releases with:

```sh
pip install "minervachem[lamel]"
```

For development from the MinervaChem monorepo, use:

```sh
uv sync --extra lamel
```

Import its public API separately from the core package:

```python
from lamel import MetaLearner
```

## scikit-learn and MAPIE

`LAMELRegressor` uses support tasks from the LAMEL dataset while its `fit`
method learns from the target-task `X_train, y_train` you supply. All feature
matrices must be LAMEL Graphlet fingerprints made with the same dataset and
`max_subgraph_size`.

```python
from lamel import LAMELRegressor
from mapie.regression import SplitConformalRegressor

estimator = LAMELRegressor(
    database_path="support_tasks.csv",
    target_task="target_property",
    support_tasks=["property_a", "property_b"],
    max_subgraph_size=5,
)
model = SplitConformalRegressor(estimator=estimator, confidence_level=0.95)
model.fit(X_train, y_train)
model.conformalize(X_cal, y_cal)
y_pred, y_interval = model.predict_interval(X_test)
```

## Citation
If you use LAMeL in your work, please cite [our paper](https://arxiv.org/abs/2509.13527).

```latex
@article{pimonova2025metalearning,
  title={Meta-Learning Linear Models for Molecular Property Prediction},
  author={Pimonova, Yulia and Taylor, Michael G and Allen, Alice and Yang, Ping and Lubbers, Nicholas},
  journal={arXiv preprint arXiv:2509.13527},
  url={https://arxiv.org/abs/2509.13527},
  year={2025}
}
```

## Copyright and licensing

`minervachem` is released under the BSD-3 License. See LICENSE.txt for the full license.


The copyright to `minervachem` is owned by Triad National Security, LLC
and is released for open source use as project number O04631.


© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.
