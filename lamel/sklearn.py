"""scikit-learn adapter for LAMEL."""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix, issparse
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_array, check_is_fitted, check_X_y

from lamel.lamel_learner import MetaLearner
from lamel.lamel_utilities import centroid_by_points, run_default_ridge_regressions


class LAMELRegressor(RegressorMixin, BaseEstimator):
    """A sklearn-compatible regressor trained with LAMEL meta-learning.

    ``database_path`` supplies a LAMEL-format dataset containing SMILES and the
    named ``support_tasks``. ``fit(X, y)`` receives target-task Graphlet
    fingerprints and labels. Its features must have the same ordering and width
    as LAMEL's fingerprints for this dataset and ``max_subgraph_size``.
    """

    def __init__(
        self, database_path, target_task, support_tasks, max_subgraph_size=5,
        random_state=16, epsilon_par=1.0, epsilon_perp=1.0,
        hierarchical=False, working_directory=None, make_log=False, a_range="a1",
    ):
        self.database_path = database_path
        self.target_task = target_task
        self.support_tasks = support_tasks
        self.max_subgraph_size = max_subgraph_size
        self.random_state = random_state
        self.epsilon_par = epsilon_par
        self.epsilon_perp = epsilon_perp
        self.hierarchical = hierarchical
        self.working_directory = working_directory
        self.make_log = make_log
        self.a_range = a_range

    def fit(self, X, y):
        """Fit LAMEL using all supplied target-task training examples."""
        X, y = check_X_y(X, y, accept_sparse="csr", y_numeric=True)
        if not issparse(X):
            X = csr_matrix(X)

        learner = MetaLearner(
            database_path=self.database_path,
            max_subgraph_size=self.max_subgraph_size,
            random_state=self.random_state,
            epsilon_par=self.epsilon_par,
            epsilon_perp=self.epsilon_perp,
            hierarchical=self.hierarchical,
            working_directory=self.working_directory,
            make_log=self.make_log,
            a_range=self.a_range,
        )
        learner.set_task_names(self.target_task, list(self.support_tasks))

        expected_features = learner.fingerprints.shape[1]
        if X.shape[1] != expected_features:
            raise ValueError(
                "X has a different number of features from LAMEL's Graphlet "
                f"fingerprints ({X.shape[1]} != {expected_features}). Generate "
                "X with this adapter's LAMEL dataset/featurization settings."
            )

        support_sets = learner.process_old_tasks(
            new_task=self.target_task, old_tasks=learner.old_task_names
        )
        support_vectors, _ = run_default_ridge_regressions(
            support_sets, a_range=learner.alpha_range
        )
        theta_bar = centroid_by_points(support_vectors)
        theta_parallel, _ = learner.solve_parallel(
            support_vectors, theta_bar, [X, y], print_alpha=False
        )
        _, theta_task, _ = learner.solve_perpendicular(
            [X, y], theta_bar, theta_parallel, print_alpha=False
        )

        self.coef_ = np.asarray(theta_task)
        self.intercept_ = 0.0
        self.n_features_in_ = X.shape[1]
        self.learner_ = learner
        return self

    def predict(self, X):
        """Predict using the meta-learned LAMEL target-task vector."""
        check_is_fitted(self, attributes=["coef_", "n_features_in_"])
        X = check_array(X, accept_sparse="csr")
        if X.shape[1] != self.n_features_in_:
            raise ValueError(
                f"X has {X.shape[1]} features, but LAMELRegressor was fitted "
                f"with {self.n_features_in_}."
            )
        return np.asarray(X @ self.coef_ + self.intercept_).ravel()
