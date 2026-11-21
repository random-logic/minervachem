import networkx as nx
from minervachem.converters import smiles_to_nx


class TestConverters:
    def test_smiles_to_nx(self):
        smiles = "CCO"  # ethanol
        G, pos = smiles_to_nx(smiles)

        assert isinstance(G, nx.Graph)
        assert len(G.nodes) == 3
        assert len(G.edges) == 2

        for node in G.nodes:
            assert "atom_key" in G.nodes[node]
            assert "atom_symbol" in G.nodes[node]
        for u, v in G.edges:
            assert "bond_key" in G.edges[u, v]
