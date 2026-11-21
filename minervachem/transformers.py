"""Provides sklearn transformer interface to minerva's fingerprinting capabiliities
"""


import pandas as pd

from joblib import Parallel, delayed

from sklearn.base import BaseEstimator, TransformerMixin

from rdkit import Chem


from buhito.transformers import GraphletTransformer, SparseFeatureEliminator

from .utils import misc


class MoleculeFeaturizer(BaseEstimator, TransformerMixin):
    """Creates RDKitMols from a pd dataframe column
    """

    def __init__(self,
                input_column='ligand',
                output_column='mol',
                n_jobs=None,
                verbose=0,
                transformer=Chem.MolFromSmiles):
        self.input_column = input_column
        self.output_column = output_column
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.transformer = transformer


    def fit(self, X, y=None):
        """Doesn't do anything, exists for compatibility"""
        self.is_fitted_ = True
        return self

    def transform(self, X):
        """Looks for smiles_column, finds unique, makes rdkit mols, joins in as output column"""

        # only convert each molecule from smiles once
        unique_mol_rows = X.drop_duplicates(subset=self.input_column)[[self.input_column]]

        n_mols = unique_mol_rows.shape[0]
        batch_size = misc.evenly_distribute_jobs(n_mols, self.n_jobs)
        print('Converting to RDKit Mol...')
        p = Parallel(n_jobs=self.n_jobs,
                    verbose=self.verbose,
                    prefer='threads',
                    batch_size=batch_size,
                    )

        # RDKit functions cannot be pickled directly, so we wrap them
        def t(s):
            return self.transformer(s)
        f = delayed(t)

        mols = p(f(mol) for mol in unique_mol_rows[self.input_column])
        unique_mol_rows[self.output_column] = mols
        X = pd.merge(X, unique_mol_rows, how='left')
        return X

    def get_feature_names(self):
        return ['mol']

class FingerprintFeaturizer(GraphletTransformer):
    """Apply a molecular fingerprinter using Buhito's shared matrix transformer."""
    def __init__(self, fingerprinter, return_dense=False, return_float=False,
                 n_jobs=None, chunk_size='even', verbose=0):
        self.fingerprinter = fingerprinter
        super().__init__(fingerprinter, return_dense=return_dense,
                         return_float=return_float, n_jobs=n_jobs,
                         chunk_size=chunk_size, verbose=verbose)

    @property
    def featurizer(self):
        return self.fingerprinter

    @featurizer.setter
    def featurizer(self, value):
        self.fingerprinter = value



    def _to_sparse(self, fps, used_bits):
        # Preserve the historical private-method DOK result for existing callers.
        return super()._to_sparse(fps, used_bits).todok()
