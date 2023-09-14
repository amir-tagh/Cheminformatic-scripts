#!/usr/bin/env python

import copy
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import display
from qed import qed
from rdkit import Chem
import sys
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
import numpy as np
import argparse, sys, pickle, math, rdkit, matplotlib
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import *
from rdkit.Chem import AllChem as Chem
import matplotlib.pylab as plt
from rdkit.Chem.Draw import SimilarityMaps
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colorbar
from matplotlib import cm
from numpy.random import randn
import matplotlib.patheffects as PathEffects
import warnings


from rdkit.Chem import Draw
from rdkit.Chem.Fraggle import FraggleSim
from rdkit.Chem.Draw import SimilarityMaps

import sys
from optparse import OptionParser
from rdkit import Chem

from collections import defaultdict
import time

from rdkit import rdBase, Chem, DataStructs
import pandas as pd
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdFingerprintGenerator

import pandas as pd
from rdkit.Chem import Draw, Lipinski, Crippen, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import rdDepictor
from itertools import combinations


#smiles = "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
#mol = Chem.MolFromSmiles(smiles)

#img = Draw.MolToImage(mol)
#img.save('output.png')


mol1 = Chem.MolFromSmiles("Cc1ccc(CNC(=O)c2ccc3nc(c4cccc(C)c4)c(c4cccc(C)c4)nc3c2)cc1")
mol2 = Chem.MolFromSmiles("Cc1ccc(C(=O)Nc2sc(N/N=C\c3cccc(Oc4ccccc4)c3)nc2c2ccc(C)cc2)cc1")
mol3 = Chem.MolFromSmiles("Cc1ccc(NC(=O)c2ccc(Nc3ccnc(c4ccc(F)cc4)n3)cc2)cc1C")
mol4 = Chem.MolFromSmiles("Cc1ccc(c2nccc(Nc3cccc(C(=O)Nc4cccc(C)c4)c3)n2)cc1")
mol5 = Chem.MolFromSmiles("O=C(Nc1ccc(c2nc3ccccc3[nH]2)cc1)c1cc2ccccc2oc1=Nc1ccccc1")
mol6 = Chem.MolFromSmiles("CCCOc1ccc(c2nn(c3ccccc3)cc2CNC(=O)c2cc(c3ccccc3)nc3ccccc23)cc1")
mol7 = Chem.MolFromSmiles("CC1CCN(CC1N(C)C2=NC=NC3=C2C=CN3)C(=O)CC#N")
mol8 = Chem.MolFromSmiles("c1ccccc1c2ccccc2OC")

mols = [mol1,mol2,mol3,mol4,mol5,mol6,mol7,mol8]

img= Draw.MolsToGridImage(mols, molsPerRow=3,returnPNG=False)
img.save('Pic1.png')


#with open('grid.svg', 'w') as f:
#    f.write(grid.data)

def is_in_samering(idx1, idx2, bond_rings):
    for bond_ring in bond_rings:
        if idx1 in bond_ring and idx2 in bond_ring:
            return True
    return False

def getLinkerbond(mol, useScaffold=True):
    res = []
    for atom in mol.GetAtoms():
        atom.SetIntProp("orig_idx", atom.GetIdx())
    for bond in mol.GetBonds():
        bond.SetIntProp("orig_idx", bond.GetIdx())

    if useScaffold:
        mol = MurckoScaffold.GetScaffoldForMol(mol)

    ring_info = mol.GetRingInfo()
    bond_rings = ring_info.BondRings()
    ring_bonds = set()
    for ring_bond_idxs in bond_rings:
        for idx in ring_bond_idxs:
            ring_bonds.add(idx)
    all_bonds_idx = [bond.GetIdx() for bond in mol.GetBonds()]
    none_ring_bonds = set(all_bonds_idx) - ring_bonds
    for bond_idx in none_ring_bonds:
        bgn_idx = mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx()
        end_idx = mol.GetBondWithIdx(bond_idx).GetEndAtomIdx()
        if mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble() == 1.0:
            if mol.GetAtomWithIdx(bgn_idx).IsInRing()+mol.GetAtomWithIdx(end_idx).IsInRing() == 1:
                bond = mol.GetBondWithIdx(bond_idx)
                orig_idx = bond.GetIntProp("orig_idx")
                res.append(orig_idx)
            elif not is_in_samering(bgn_idx, end_idx, bond_rings) and mol.GetAtomWithIdx(bgn_idx).IsInRing()+mol.GetAtomWithIdx(end_idx).IsInRing() == 2:
                bond = mol.GetBondWithIdx(bond_idx)
                orig_idx = bond.GetIntProp("orig_idx")
                res.append(orig_idx)
    return res

for mol in mols:
    bonds = getLinkerbond(mol)
    if bonds:
        res = Chem.FragmentOnBonds(mol, bonds)
        display(res)
    else:
        display(mol)


frgs = Chem.GetMolFrags(res, asMols=True)
img = Draw.MolsToGridImage(frgs, molsPerRow=5,returnPNG=False)
img.save('Pic2.png')





