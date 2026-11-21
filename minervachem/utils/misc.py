"""Molecular annotation and compatibility imports for generic utilities."""
from buhito.utilities import evenly_distribute_jobs, uniquify_lol
from buhito.utils.misc import (
    rank,
    timeit,
    sizeof_fmt,
    bits_by_size,
    count_fragments_by_size,
    group_bitinfo_by_size,
    counts_by_bit_by_size,
    scan_magnitudes,
    find_inds,
    get_graphs_with_bit as get_mols_with_bit,
)

def mol_with_atom_index(mol, index_start=0, note=False):
    #mol = deepcopy(mol)
    for atom in mol.GetAtoms():
        ix = atom.GetIdx()+index_start
        if note:
            atom.SetProp('atomNote', str(ix))
        else:
            atom.SetAtomMapNum(ix)
    return mol
