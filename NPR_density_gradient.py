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
from scipy import stats

#config Completer.use_jedi = False
#PandasTools.RenderImagesInAllDataFrames(images=True)

#############################################
#add the index to smi file befor calculation#
#############################################


mols=Chem.SDMolSupplier(sys.argv[1],sanitize = False)

print (len(mols)) #To check how many molecules there are in the file

#bar=progressbar.ProgressBar(max_value=len(mols))
table1=pd.DataFrame()

for i,mol in enumerate(mols):
    #Chem.SanitizeMol(mol)
    
    table1.loc[i,'smiles']=Chem.MolToSmiles(mol)
    table1.loc[i,'NPR1']=rdMolDescriptors.CalcNPR1(mol)
    table1.loc[i,'NPR2']=rdMolDescriptors.CalcNPR2(mol)
    #bar.update(i+1)
pd.set_option('display.max_rows', None)
print(table1.head(251)) #Let's take a look to the table

first_column = table1.iloc[:, 1]

#plt.rcParams['axes.linewidth'] = 1.5

x=table1.iloc[:, 1].to_numpy()
print(np.prod(x.shape))

y=table1.iloc[:, 2].to_numpy()
print(y)
print(np.prod(y.shape))

bins = 100
hh, locx, locy = np.histogram2d(x, y, bins=bins)

# Sort the points by density, so that the densest points are plotted last
z = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(x,y)])
idx = z.argsort()
x2, y2, z2 = x[idx], y[idx], z[idx]



plt.figure(1,figsize=(12,8)).clf()
plt.ylim([0.45, 1])
#plt.xlabel('NPR1',fontsize=18)
#plt.ylabel('NPR2',fontsize=18)

s = plt.scatter(x2, y2, c=z2, cmap='jet', marker='.',s=30)

axes = plt.gca() #Getting the current axis

axes.spines['top'].set_visible(False)
axes.spines['right'].set_visible(False)


#change axis thickness
axes.spines['left'].set_linewidth(3.0)
axes.spines['bottom'].set_linewidth(3.0)


#change tick thickness
axes.xaxis.set_tick_params(width=2.5)
axes.yaxis.set_tick_params(width=2.5)

#axis font size
axes.xaxis.set_tick_params(labelsize=24)
axes.yaxis.set_tick_params(labelsize=24)

#plt.setp(axes.spines.values(), visible=False) 
#plt.axis('off')

#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Counts')


#plt.show()


x1, y1 = [0.5, 0], [0.5, 1]
x2, y2 = [0.5, 1], [0.5, 1]
x3, y3 = [0,1],[1,1]

plt.plot(x1, y1,x2,y2,x3,y3,c='black',ls='--',lw=1)

#plt.xlabel ('NPR1',fontsize=24,fontweight='bold',color='black',labelpad=10)

#plt.ylabel ('NPR2',fontsize=24,fontweight='bold',color='black',labelpad=10)

#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)

plt.text(0, 1.02,s='Rod',fontsize=22,horizontalalignment='center',verticalalignment='center',fontweight='bold')#1.01
plt.text(1, 1.02,s='Sphere',fontsize=22,horizontalalignment='center',verticalalignment='center',fontweight='bold')#1.01
plt.text(0.5, 0.47,s='Disc',fontsize=22,horizontalalignment='center',verticalalignment='center',fontweight='bold')#0.49
#ax.set_ylim([0.46, 1.0])

from matplotlib.lines import Line2D

#custom = [Line2D([], [], marker='o', color='lightskyblue', linestyle='None',lw=2,markersize=15),
#          Line2D([], [], marker='o', color='orange', linestyle='None',lw=2,markersize=15)]
#

#plt.legend(custom, ['Drug Bank', 'New RNA focus library'], loc='lower right',prop={'size': 24})

#plt.xticks(fontsize=60,fontname = "Arial")#fontweight='bold')
#plt.yticks(fontsize=60,fontname = "Arial")#fontweight='bold')
#plt.tick_params ('both',width=2,labelsize=18)


#plt.tight_layout()
#change axis thickness
axes.spines['left'].set_linewidth(4.0)
axes.spines['bottom'].set_linewidth(4.0)


#change tick thickness
axes.xaxis.set_tick_params(width=3.0)
axes.yaxis.set_tick_params(width=3.0)

axes.xaxis.set_tick_params(length=18.0)
axes.yaxis.set_tick_params(length=18.0)

plt.xticks(fontsize=28,weight='bold')
plt.yticks(fontsize=28,weight='bold')


plt.savefig('NPR_Main.png',dpi=400)
plt.show()



