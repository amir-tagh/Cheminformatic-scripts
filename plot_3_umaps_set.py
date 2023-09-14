#!/usr/bin/env python

# Install UMAP/hdbscan (recommend to do this before miniconda)
#pip install --quiet umap-learn hdbscan

# Install miniconda
#MINICONDA_INSTALLER_SCRIPT=Miniconda3-py37_4.9.2-Linux-x86_64.sh
#MINICONDA_PREFIX=/usr/local
#wget -q https://repo.continuum.io/miniconda/$MINICONDA_INSTALLER_SCRIPT
#chmod +x $MINICONDA_INSTALLER_SCRIPT
#./$MINICONDA_INSTALLER_SCRIPT -b -f -p $MINICONDA_PREFIX > /dev/null

# Install rdkit (should only be installed via conda)
#conda install -y --quiet -c conda-forge rdkit=2020.09.2 > /dev/null



#Import modules
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import seaborn as sns
import umap
import sys
import matplotlib.pyplot as plt
from rdkit.Chem.AtomPairs import Pairs, Torsions

#use the cv type of data type for UMAP
#UMAP built with Morgan fingerprints, 1024 data bit, radius 2

# both df1 and df2 have bivaraite normals, df1.size=200, df2.size=100

df1 = pd.read_csv('AP_NOHASH.csv')
print('df1',df1)

df2 = pd.read_csv('DRUGBANK.csv')
print('df2',df2)



# plot
# ========================================   
graph = sns.jointplot(x=df1.UMAP2, y=df1.UMAP3,
        kind="scatter",color="seagreen",
        alpha=0.1,space=0)
        #marginal_kws=dict(bins=100,alpha=0.5))
        

#graph.ax_joint.set_xlabel('UMAP1', fontweight='bold',fontsize=28)
#graph.ax_joint.set_ylabel('UMAP2', fontweight='bold')


graph.x = df2.UMAP2
graph.y = df2.UMAP3


graph.plot_marginals(plt.hist,bins=100,color="darkorange",histtype='bar',
        density=False,edgecolor = 'black',linewidth=0.1,align='left',alpha=0.3)

#graph.plot_marginals(sns.kdeplot, color="black", shade=True)

graph.plot_joint(plt.scatter, color="darkorange",s=5,marker='o',
        alpha=.2, label='Drug Bank')



plt.xticks(fontsize=26)
plt.yticks(fontsize=26)
plt.tick_params(size=8,width=2)


graph.ax_joint.spines['left'].set_linewidth(4.0)
graph.ax_joint.spines['bottom'].set_linewidth(4.0)


#turn off marginal axis
graph.ax_marg_x.set_axis_off()
graph.ax_marg_y.set_axis_off()

graph.ax_joint.set(xlabel=None)
graph.ax_joint.set(ylabel=None)

#change axis thickness
graph.ax_joint.spines['left'].set_linewidth(4.0)
graph.ax_joint.spines['bottom'].set_linewidth(4.0)
graph.ax_joint.spines['right'].set_linewidth(4.0)
graph.ax_joint.spines['top'].set_linewidth(4.0)


#change tick thickness
#graph.ax.xaxis.set_tick_params(width=3.0)
#graph.ax.yaxis.set_tick_params(width=3.0)

graph.ax_joint.tick_params(length=15.0,width=3.0)

plt.xticks(fontsize=24,weight='bold')
plt.yticks(fontsize=24,weight='bold')


plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(12.5, 10.5)

plt.savefig('UMAP-main.png',dpi=300)
plt.show()




