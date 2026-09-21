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
