"""Check the chemistry boundary and the historical molecular fingerprint values."""
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from sklearn.base import clone

from buhito.transformers import GraphletTransformer
from buhito.regressors import HierarchicalResidualModel
from minervachem import regressors
from minervachem.fingerprinters import GraphletFingerprinter
from minervachem.graphlet_dags import GraphletDAG
from minervachem.transformers import FingerprintFeaturizer


@pytest.mark.parametrize('smiles', ['CCO', 'c1ccccc1', '[NH4+]'])
def test_historical_molecular_fingerprints(smiles):
    expected = json.loads((Path(__file__).parent / 'fixtures/molecular_fingerprints.json').read_text())
    fp, _ = GraphletFingerprinter(max_len=3)(Chem.MolFromSmiles(smiles))
    assert sorted([size, bit, count] for (size, bit), count in fp.items()) == expected[smiles]


def test_molecular_transformer_delegates_to_buhito():
    original = FingerprintFeaturizer(GraphletFingerprinter(), n_jobs=None)
    transformer = clone(original)
    assert isinstance(transformer, GraphletTransformer)
    mols = [Chem.MolFromSmiles(s) for s in ['CCO', 'CCC']]
    train = transformer.fit_transform(mols)
    assert np.array_equal(train.toarray(), transformer.transform(mols).toarray())
    smaller = transformer.extract_smaller_fp(train, new_max_subgraph_size=2)
    assert smaller.shape[1] == (transformer.bit_sizes_ <= 2).sum()
    replacement = GraphletFingerprinter(max_len=2)
    transformer.set_params(fingerprinter=replacement)
    assert transformer.featurizer is replacement
    assert regressors.HierarchicalResidualModel is HierarchicalResidualModel


@pytest.mark.parametrize('layer', [1, 2])
def test_molecular_projection_preserves_coefficients(layer):
    mol = Chem.MolFromSmiles('CCO')
    fp, locations = GraphletFingerprinter()(mol)
    coefficients = {bit: float(bit[0]) for bit in fp}
    expected = sum(coefficients[bit] * len(nodes) for bit, nodes in locations.items())
    dag = GraphletDAG(mol, bi=locations, bit_to_coef_map=coefficients)
    projected = dag.project_to_layer(layer)
    assert np.isclose(projected.sum(), expected)
    assert np.allclose(dag.project_to_layer(layer), projected)
