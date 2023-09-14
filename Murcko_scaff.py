#!/usr/bin/env python


import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors,Crippen

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import matplotlib.cm as cm
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans
import time


import numpy as np
import pandas as pd


from FPSim2.FPSim2lib.utils import (
    PyPopcount,
    SortResults,
)
from FPSim2.FPSim2lib import (
    TanimotoSearch,
    TverskySearch,
    SubstructureScreenout,
)
from FPSim2 import FPSim2Engine
import os
import sys
from FPSim2.io import create_db_file

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
from rdkit.Chem import rdMolDescriptors as rdmd
from rdkit.Chem.Scaffolds import MurckoScaffold

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmilesFromSmiles
import mols2grid
from tqdm.auto import tqdm
from ipywidgets import widgets
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True


from tqdm.notebook import tqdm
tqdm.pandas()

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
'''
data_1 = open(sys.argv[1], 'r')
lines_1 = data_1.readlines()
data_1.close()


for line in lines_1:
    
    q = Chem.MolFromSmiles(line)
    if q is None:
        continue
    scaff = MurckoScaffold.GetScaffoldForMol(q)
    res=Chem.MolToSmiles(scaff)
    print(res) 
    #framework = MurckoScaffold.MakeScaffoldGeneric(scaff)
    #print(Chem.MolToSmiles(framework))
'''

url = "AP_nohash_cluster_centers_3000.smi"
df = pd.read_csv(url,names=["SMILES","Name","Cluster","Center"])
print(df.head())
df['framework'] = df.SMILES.progress_apply(MurckoScaffoldSmilesFromSmiles)


df['mol'] = df.SMILES.progress_apply(Chem.MolFromSmiles)

scaffold_df = df.framework.value_counts().reset_index().copy() # cool trick
scaffold_df.columns = ["scaffold","count"]
# copy the index column for the dataframe to scaffold_idx
scaffold_df['scaffold_idx'] = scaffold_df.index

scaffold_df.to_csv('3000_clusters_scafolds.csv', sep='\t', encoding='utf-8')



