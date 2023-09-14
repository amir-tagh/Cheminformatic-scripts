#!/usr/bin/env python

###first install these packages
#conda create -n chempca_venv python=3.7
#conda activate chempca_venv
#conda install -c rdkit rdkit
#conda install pandas numpy
#pip install matplotlib
#pip install -U scikit-learn scipy matplotlib



#edit the input file
#awk -F'\t' 'NR>1{$0=$0"\t"NR-1} 1' file

#awk -F'\t' 'NR>0{$0=$0"\t"NR-1} 1' inforna20210518_molecules_smiles.smi


#convert comma sep. to tab sep.
#sed 's/\,/\t/g' file >outfile


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
'''
Plotting libraries
'''
import pandas as pd
import matplotlib.cm as cm
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

'''
What we'll need for analysis, clustering, etc.
'''
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans
from sklearn import datasets, decomposition
from sklearn.manifold import TSNE
'''
Of course the powerful RDKIT for cheminformatics <3
'''
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, Descriptors3D, Draw, rdMolDescriptors, Draw, PandasTools
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoSimMat, GetTanimotoDistMat
from rdkit.Chem.Draw import IPythonConsole

'''
Some utilities
'''
import progressbar
from math import pi
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw


m = Chem.MolFromSmiles('Cc1cccc(NC(=O)c2ccc3c(c2)C(=O)N(c2ccc(Oc4ccc(N5C(=O)c6ccc(C(=O)Nc7cccc(C)c7)cc6C5=O)cc4)cc2)C3=O)c1')
mol = MurckoScaffold.GetScaffoldForMol(m)



w = Chem.MolToSmiles(m)
print(w)









