"""LAMEL: an optional meta-learning add-on for MinervaChem."""

from lamel.dataset_loader import LoadedDataset
from lamel.lamel_learner import MetaLearner

__all__ = ["LoadedDataset", "MetaLearner"]
