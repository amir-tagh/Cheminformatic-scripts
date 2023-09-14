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

#config Completer.use_jedi = False
#PandasTools.RenderImagesInAllDataFrames(images=True)

#############################################
#add the index to smi file befor calculation#
#############################################

#For this workflow I'll use FDA apporved molecules taken from ZINC15
#[https://zinc15.docking.org/catalogs/dbfda/substances/subsets/world/]

#If your file format is Mol2, maybe you can need my Mol2MolSupplier for RDKIT.
#[https://chem-workflows.com/articles/2020/03/23/building-a-multi-molecule-mol2-reader-for-rdkit-v2/]

mols=Chem.SDMolSupplier(sys.argv[1])
print (len(mols)) #To check how many molecules there are in the file

#bar=progressbar.ProgressBar(max_value=len(mols))
table=pd.DataFrame()

for i,mol in enumerate(mols):
    Chem.SanitizeMol(mol)
    table.loc[i,'smiles']=Chem.MolToSmiles(mol)
    table.loc[i,'Mol']=mol
    table.loc[i,'MolWt']=Descriptors.MolWt(mol)
    table.loc[i,'LogP']=Descriptors.MolLogP(mol)
    table.loc[i,'HBA']=Descriptors.NumHAcceptors(mol)
    table.loc[i,'HBD']=Descriptors.NumHDonors(mol)
    table.loc[i,'NumHeteroatoms']=Descriptors.NumHeteroatoms(mol)
    table.loc[i,'NumRotatableBonds']=Descriptors.NumRotatableBonds(mol)
    table.loc[i,'NumHeavyAtoms']=Descriptors.HeavyAtomCount (mol)
    table.loc[i,'NumAliphaticCarbocycles']=Descriptors.NumAliphaticCarbocycles(mol)
    table.loc[i,'NumAliphaticHeterocycles']=Descriptors.NumAliphaticHeterocycles(mol)
    table.loc[i,'NumAliphaticRings']=Descriptors.NumAliphaticRings(mol)
    table.loc[i,'NumAromaticCarbocycles']=Descriptors.NumAromaticCarbocycles(mol)
    table.loc[i,'NumAromaticHeterocycles']=Descriptors.NumAromaticHeterocycles(mol)
    table.loc[i,'NumAromaticRings']=Descriptors.NumAromaticRings(mol)
    table.loc[i,'RingCount']=Descriptors.RingCount(mol)
    table.loc[i,'FractionCSP3']=Descriptors.FractionCSP3(mol)
    
    table.loc[i,'TPSA']=Descriptors.TPSA(mol)
    table.loc[i,'NPR1']=rdMolDescriptors.CalcNPR1(mol)
    table.loc[i,'NPR2']=rdMolDescriptors.CalcNPR2(mol)
    table.loc[i,'InertialShapeFactor']=Descriptors3D.InertialShapeFactor(mol)
    table.loc[i,'RadiusOfGyration']=Descriptors3D.RadiusOfGyration(mol)
    #bar.update(i+1)
pd.set_option('display.max_rows', None)
print(table.head(251)) #Let's take a look to the table

time.sleep(5)

plt.rcParams['axes.linewidth'] = 1.5
plt.figure(figsize=(6,4))

ax=sns.scatterplot(x='NPR1',y='NPR2',data=table,s=10,linewidth=0.5,alpha=1)
x1, y1 = [0.5, 0], [0.5, 1]
x2, y2 = [0.5, 1], [0.5, 1]
x3, y3 = [0,1],[1,1]

plt.plot(x1, y1,x2,y2,x3,y3,c='gray',ls='--',lw=2)

plt.xlabel ('NPR1',fontsize=20,fontweight='bold')

plt.ylabel ('NPR2',fontsize=20,fontweight='bold')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.text(0, 1.01,s='Rod',fontsize=16,horizontalalignment='center',verticalalignment='center',fontweight='bold')
plt.text(1, 1.01,s='Sphere',fontsize=16,horizontalalignment='center',verticalalignment='center',fontweight='bold')
plt.text(0.5, 0.49,s='Disc',fontsize=16,horizontalalignment='center',verticalalignment='center',fontweight='bold')

plt.tick_params ('both',width=2,labelsize=14)
plt.tight_layout()

plt.savefig('Inforna.png',dpi=300)
plt.show()

descriptors = table[['MolWt', 'LogP','NumHeteroatoms','RingCount','FractionCSP3', 'TPSA','RadiusOfGyration']].values #The non-redundant molecular descriptors chosen for PCA
descriptors_std = StandardScaler().fit_transform(descriptors) #Important to avoid scaling problems between our different descriptors
pca = PCA()
descriptors_2d = pca.fit_transform(descriptors_std)
descriptors_pca= pd.DataFrame(descriptors_2d) # Saving PCA values to a new table
descriptors_pca.index = table.index
descriptors_pca.columns = ['PC{}'.format(i+1) for i in descriptors_pca.columns]
print(descriptors_pca.head(20)) #Displays the PCA table

print(pca.explained_variance_ratio_) #Let's plot PC1 vs PC2
print(sum(pca.explained_variance_ratio_))


# This normalization will be performed just for PC1 and PC2, but can be done for all the components.
#The normalization is to plot PCA values in 0-1 sacle and include the vectors (features to the plot)

scale1 = 1.0/(max(descriptors_pca['PC1']) - min(descriptors_pca['PC1'])) 
scale2 = 1.0/(max(descriptors_pca['PC2']) - min(descriptors_pca['PC2']))

# And we add the new values to our PCA table
descriptors_pca['PC1_normalized']=[i*scale1 for i in descriptors_pca['PC1']]
descriptors_pca['PC2_normalized']=[i*scale2 for i in descriptors_pca['PC2']]


print(descriptors_pca.head(20)) # The PCA table now has the normalized PC1 and PC2

plt.rcParams['axes.linewidth'] = 1.5
plt.figure(figsize=(6,6))

ax=sns.scatterplot(x='PC1_normalized',y='PC2_normalized',data=descriptors_pca,s=20,palette=sns.color_palette("Set2", 3),linewidth=0.2,alpha=1)

plt.xlabel ('PC1',fontsize=20,fontweight='bold')
ax.xaxis.set_label_coords(0.98, 0.45)
plt.ylabel ('PC2',fontsize=20,fontweight='bold')
ax.yaxis.set_label_coords(0.45, 0.98)

ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


lab=['MolWt', 'LogP','NumHeteroatoms','RingCount','FractionCSP3',
     'NumNHOH', 'NumNO', 'TPSA','PBF',
     'InertialShapeFactor','RadiusOfGyration'] #Feature labels

l=np.transpose(pca.components_[0:2, :]) ## We will get the components eigenvectors (main features) for PC1 and PC2

n = l.shape[0]
for i in range(n):
    plt.arrow(0, 0, l[i,0], l[i,1],color= 'k',alpha=0.5,linewidth=1.8,head_width=0.025)
    plt.text(l[i,0]*1.25, l[i,1]*1.25, lab[i], color = 'k',va = 'center', ha = 'center',fontsize=16)

circle = plt.Circle((0,0), 1, color='gray', fill=False,clip_on=True,linewidth=1.5,linestyle='--')

plt.tick_params ('both',width=2,labelsize=18)

ax.add_artist(circle)
plt.xlim(-1.2,1.2)
plt.ylim(-1.2,1.2)
plt.tight_layout()
plt.savefig('PC_circular_inforna.png',dpi=300)
plt.show()


smiles = list(table["smiles"])
smi=[Chem.MolFromSmiles(x) for x in smiles]
fps =  [MACCSkeys.GenMACCSKeys(x) for x in smi] # In this example I'll use MACCSKeys
tanimoto_sim_mat_lower_triangle=GetTanimotoSimMat(fps) #This compute a similartity matrix between all the molecules
n_mol = len(fps)
similarity_matrix = np.ones([n_mol,n_mol])
i_lower= np.tril_indices(n=n_mol,m=n_mol,k=-1)
i_upper= np.triu_indices(n=n_mol,m=n_mol,k=1)
similarity_matrix[i_lower] = tanimoto_sim_mat_lower_triangle
similarity_matrix[i_upper] = similarity_matrix.T[i_upper]
distance_matrix = np.subtract(1,similarity_matrix) #This is the similarity matrix of all vs all

TSNE_sim = TSNE(n_components=2,init='pca',random_state=90, angle = 0.3,perplexity=50).fit_transform(distance_matrix) #Remember to always tune the parameters acording your dataset!!
tsne_result = pd.DataFrame(data = TSNE_sim , columns=["TC1","TC2"]) # New table containing the tSNE results
print(tsne_result.head(5)) #A new table containing the tSNE results

plt.rcParams['axes.linewidth'] = 1.5
fig, ax = plt.subplots(figsize=(6,6))

ax=sns.scatterplot(x='TC1',y='TC2',data=tsne_result,s=15,linewidth=0.2,alpha=1)


plt.xlabel ('tSNE 1',fontsize=24,fontweight='bold')

plt.ylabel ('tSNE 2',fontsize=24,fontweight='bold')

plt.tick_params ('both',width=2,labelsize=18)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
handles, labels = ax.get_legend_handles_labels()

#ax.legend(handles=handles[1:], labels=labels[1:])

#plt.legend(loc='lower right',frameon=False,prop={'size': 22},ncol=1)

plt.tight_layout()
plt.show()

#K-means clustering
range_n_clusters = [2, 3, 4, 5, 6, 7, 8, 9, 10,12] # To explore the "best" number of cluster to clasify our molecules
for n_clusters in range_n_clusters:
    kmeans = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = kmeans.fit_predict(tsne_result[['TC1','TC2']])
    silhouette_avg = silhouette_score(tsne_result[['TC1','TC1']], cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg) #This will print the silhouette score, as higher as our data is better distributed inside the clusters

#Despite the fact that the silhouette calculations yielded very poor results in classifying the molecules,
#I'll cluster the dataset into six clusters for demonstration purposes.    

kmeans = KMeans(n_clusters=12, random_state=12) # We define the best number of clusters (6)
clusters = kmeans.fit(tsne_result[['TC1','TC2']]) #TC1vs TC2

tsne_result['Cluster'] = pd.Series(clusters.labels_, index=tsne_result.index)
print(tsne_result.head(200)) #The tSNE table now contains the numer of cluster for each element
print(type(tsne_result))
tsne_result.to_csv('t-SNE.csv')


plt.rcParams['axes.linewidth'] = 1.5
fig, ax = plt.subplots(figsize=(6,6))

ax=sns.scatterplot(x='TC1',y='TC2',data=tsne_result, hue='Cluster',s=22,palette=sns.color_palette("Set2", 12),linewidth=0.2,alpha=1)


plt.xlabel ('tSNE 1',fontsize=24,fontweight='bold')

plt.ylabel ('tSNE 2',fontsize=24,fontweight='bold')

plt.tick_params ('both',width=2,labelsize=18)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

handles, labels = ax.get_legend_handles_labels()
plt.xticks([-60,-40,-20,0,20,40,60])
ax.legend(handles=handles[1:], labels=labels[1:])

plt.legend(loc='best',frameon=False,prop={'size': 16},ncol=2)
#remove the legends
plt.legend([],[], frameon=False)
plt.tight_layout()
plt.savefig('ML_tSNE1.png',dpi=300)
plt.show()

data=pd.DataFrame() # I'll create a new table containing the normalized bRo5 values of our compounds


data['MolWt']=[i/500 for i in table['MolWt']]
data['LogP']=[i/5 for i in table['LogP']]
data['HBA']=[i/10 for i in table['HBA']]
data['HBD']=[i/5 for i in table['HBD']]
data['RotB']=[i/10 for i in table['NumRotatableBonds']]
data['TPSA']=[i/140 for i in table['TPSA']]
print(data.to_string())

categories=list(data.columns)  # This will set up the parameters for the angles of the radar plot.
N = len(categories)
print('len is',N)
values=data[categories].values[0]
print('values here are',values)
values=np.append(values,values[:1])
print(values)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]
print('angles are',angles)
Ro5_up=[1,1,1,1,1,1,1] #The upper limit for bRo5
Ro5_low=[0.5,0.1,0,0.25,0.1,0.5,0.5]  #The lower limit for bRo5

#fig=plt.figure(figsize=(6,6))

#ax = fig.add_axes([1, 1, 1, 1],projection='polar')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10),
                       subplot_kw=dict(polar=True))

plt.xticks(angles[:-1], categories,color='k',size=20,ha='center',va='top',fontweight='book')

plt.tick_params(axis='y',width=4,labelsize=12, grid_alpha=0.05)

ax.set_rlabel_position(0)

#ax.plot(angles, Ro5_up, linewidth=2, linestyle='-',color='red')
#ax.plot(angles, Ro5_low, linewidth=2, linestyle='-',color='red')

#ax.fill(angles, Ro5_up, 'red', alpha=0.2)
ax.fill(angles, Ro5_low, 'orangered', alpha=0.2)

for i in data.index[:]: #I'll just show the profile for the first 250 molecules in the table for clarity of the plot
    values=data[categories].values[i]
    values=np.append(values,values[:1])
    ax.plot(angles, values, linewidth=0.7 ,color='steelblue',alpha=0.5)
    ax.fill(angles, values, 'C2', alpha=0.025)

ax.grid(axis='y',linewidth=1.5,linestyle='dotted',alpha=0.8)
ax.grid(axis='x',linewidth=2,linestyle='-',alpha=1)

ax.autoscale()
plt.savefig('radar_chart.png',dpi=300)
plt.show()









