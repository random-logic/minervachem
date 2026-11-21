from functools import partial
from io import BytesIO
import math

import cairosvg
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from PIL import Image
import sigfig
from sklearn.metrics import (r2_score,
                             mean_absolute_error as MAE,
                             mean_squared_error as MSE)
import warnings

from minervachem.utils.misc import mol_with_atom_index
from minervachem.fingerprinters import GraphletFingerprinter, MorganFingerprinter, RDKitFingerprinter


def remove_bg(im):
    arr = np.array(im)
    white = (arr.sum(2) == 765)
    alpha = np.zeros_like(white)
    alpha[~white] = 255
    alpha = Image.fromarray(alpha)
    im.putalpha(alpha)
    return


from buhito.evaluation_plotting import (
    add_identity, add_best_fit, scatter_hist, add_textbox, parity_plot, parity_plot_train_test, RMSE,
)

def draw_ss(mol, atom_ix, draw_atom_ix=False, ix_note=True):
    """Draw a substructure of a molecule given an iterable of atom indices"""
    #atoms = set([mol.GetAtomWithIdx(i) for i in atom_ix])
    rwmol = Chem.RWMol(mol)
    dos = Chem.Draw.MolDrawOptions()
    if draw_atom_ix and ix_note:
        dos.annotationFontScale=0.7
    if draw_atom_ix:
        rwmol = mol_with_atom_index(rwmol, index_start=1, note=ix_note)
    for atom in reversed(list(rwmol.GetAtoms())):
        if atom.GetIdx() not in atom_ix:
            rwmol.RemoveAtom(atom.GetIdx())
    im = Chem.Draw.MolToImage(rwmol, options=dos)
    return im

def draw_substructure(mol, atom_ix, draw_atom_ix=False, ix_note=True,
                      size=(300, 250), flatten_bonds=True, use_coordgen=True):

    """
    Render a molecular substructure defined by a set of atom indices as an SVG.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Input molecule to draw.
    atom_ix : Iterable[int]
        Atom indices (with respect to `mol`) to keep in the rendered substructure.
        Duplicates are ignored.
    draw_atom_ix : bool, default=False
        If True, annotate atoms with their indices in the substructure (after
        deletion/renumbering).
    ix_note : bool, default=True
        If True and `draw_atom_ix` is True, reduce annotation font scale to make
        indices less visually dominant.
    size : tuple[int, int], default=(300, 250)
        Output image size as (width, height) in pixels.
    flatten_bonds : bool, default=True
        If True, remove stereochemistry from the fragment prior to drawing.
    use_coordgen : bool, default=True
        If True, prefer CoordGen for 2D coordinate generation (via
        `rdDepictor.SetPreferCoordGen`).

    Returns
    -------
    str
        SVG text of the rendered substructure.
    """

    keep = sorted(set(atom_ix))
    if not keep:
        raise ValueError("No atoms selected for the substructure.")

    rwm = Chem.RWMol(mol)
    for a in sorted((i for i in range(rwm.GetNumAtoms()) if i not in keep), reverse=True):
        rwm.RemoveAtom(a)
    sub = rwm.GetMol()
    if sub.GetNumAtoms() == 0:
        raise ValueError("Fragment is empty after deletion.")

    if flatten_bonds:
        Chem.RemoveStereochemistry(sub)
    if sub.GetNumAtoms() > 1:
        sub.RemoveAllConformers()
        rdDepictor.SetPreferCoordGen(bool(use_coordgen))
        rdDepictor.Compute2DCoords(sub)

    d2d = rdMolDraw2D.MolDraw2DSVG(*size)
    opts = d2d.drawOptions()
    if draw_atom_ix:
        opts.addAtomIndices = True
        if ix_note:
            opts.annotationFontScale = 0.7

    rdMolDraw2D.PrepareAndDrawMolecule(d2d, sub)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def pil_from_svg(svg_text: str, *, scale=1.0, output_width=None, output_height=None):
    """
    Convert an SVG string into a Pillow (PIL) RGBA image.

    Parameters
    ----------
    svg_text : str
        SVG markup as a Unicode string.
    scale : float, default=1.0
        Uniform scaling factor passed to CairoSVG during rasterization.
    output_width : int | None, default=None
        Target output width in pixels. If provided, overrides the SVG's intrinsic width.
    output_height : int | None, default=None
        Target output height in pixels. If provided, overrides the SVG's intrinsic height.

    Returns
    -------
    PIL.Image.Image
        The rasterized image in RGBA mode.
    """
    png_bytes = cairosvg.svg2png(bytestring=svg_text.encode("utf-8"),
                                 scale=scale,
                                 output_width=output_width,
                                 output_height=output_height)
    return Image.open(BytesIO(png_bytes)).convert("RGBA")

def plot_fingerprint(mol,
                     fingerprinter,
                     ncol=3,
                     figsize=None,
                     decreasing=True,
                     show_count=True,
                     show_bit_ids=True,
                     show_size=True,
                     only_found=True):
    """Plot all of the substructures induced in mol by fingerprinter.

    Note that if the fingerprinter produces folded bits, this will only show the first one found.


    Args:
        mol: (rdkit.Molecule) the molecule to be fingerprinted
        fingerprinter: a minerva.fingerprinters.Fingerprinter object
        ncol: (int) the number of columns in the plot (n_rows = ceil(n_fragments / n_cols))
        figsize: passed to plt.subplots
        decreasing: (bool) sort fragments in order of decreasing size (if False, sort in increasing order)
        show_count: (bool) show the count of each fragment in the molecule
        show_bit_ids: (bool) show the ID of each bit
        show_size: (bool) show the size of the fragment (this is also included in the bit ID)
    """

    fp, bi = fingerprinter(mol)
    if only_found:
        fp = {k: v for k, v in fp.items() if v != 0}
        bi = {k: v for k, v in bi.items() if v != []}
    n_bits = len(bi.keys())

    if n_bits == 0:
        print('Found no fingerprint bits, exiting')
        return
    nrow = int(np.ceil(len(bi.keys()) / ncol))
    if figsize is None:
        figsize = (10, 5*nrow)
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize)


    # helper plotting functions for each type of fingerprinter
    _bi = {k[1]: v for k, v in bi.items()}
    if isinstance(fingerprinter, GraphletFingerprinter):
        def draw_fn(mol, bit):
            return draw_ss(mol, bi[bit][0])
        size_name = 'n'
        fingerprinter_name = 'Graphlet Fingerprint'
    elif isinstance(fingerprinter, MorganFingerprinter):
        def draw_fn(mol, bit):
            return Chem.Draw.DrawMorganBit(mol, bit[1], _bi)
        size_name = 'r'
        fingerprinter_name = 'Morgan Fingerprint'
    elif isinstance(fingerprinter, RDKitFingerprinter):
        def draw_fn(mol, bit):
            return Chem.Draw.DrawRDKitBit(mol, bit[1], _bi)
        size_name = 'l'
        fingerprinter_name = 'RDKit Fingerprint'
    else:
        raise ValueError(f'Got unsupported fingerprinter of type ({type(fingerprinter)})')

    axes = axes.ravel()
    bit_ids_by_l = sorted(bi.keys(), key=lambda k: k[0], reverse=decreasing)
    for i, bit in enumerate(bit_ids_by_l):
        plt.sca(axes[i])
        im = draw_fn(mol, bit)
        plt.imshow(im)
        plt.axis('off')

        plt_text = []
        if show_bit_ids:
            plt_text.append(f'Bit ID: {bit}')
        if show_size:
            plt_text.append(f'{size_name} = {bit[0]}')
        if show_count:
            plt_text.append(f'count = {fp[bit]}')
        plt.text(0, 10, '\n'.join(plt_text))

    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)
    #plt.suptitle(f'{fingerprinter_name} substructures for {Chem.MolToSmiles(mol)}')
    plt.tight_layout()
    return fig, axes
