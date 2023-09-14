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
import pickle as pkl


#Read csv file from SARS-Cov-2 Database in ChEMBL (Downloaded in Nov. 28 2020)
#Downloaded from https://www.ebi.ac.uk/chembl/g/#browse/compounds/filter/_metadata.compound_records.src_id%3A52

#this option is to save all the table content and not just the head
#pd.set_option("display.max_rows", None, "display.max_columns", None)



#df = pd.read_csv("./sras_embl.csv", sep=';', engine='python')
#print(df)

#Define methods for fingerprints
#Computed using 2048 bits with radius 2

#Convert list of SMILES to list of fingerprints
def fingerprint_list_from_smiles_list(smiles_list, n_bits=2048):
    fingerprint_list = []
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        fingerprint_list.append(fingerprint_as_array(mol, n_bits))
    return fingerprint_list

def fingerprint_as_array(mol, n_bits=2048):
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)   #AllChem.GetMorganFingerprintAsBitVect(q,3,2048)
    array = np.zeros((1,), np.int)
    DataStructs.ConvertToNumpyArray(fingerprint, array)
    return array


#Convert list of SMILES for SARS_CoV_2 to fingerprint
#fingerprint_SARS_CoV_2_smiles_filter = fingerprint_list_from_smiles_list(sys.argv[1],2048)
fingerprint_list = []

data_1 = open(sys.argv[1], 'r')
lines_1 = data_1.readlines()
data_1.close()


for line in lines_1:
    mol = Chem.MolFromSmiles(line)
    if mol is None:
        print(line)
        continue
    else:
        
        fingerprint_list.append(fingerprint_as_array(mol, 1024))
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
        #fingerprint = Pairs.GetAtomPairFingerprint(mol)
        array = np.zeros((1,), np.int)
        DataStructs.ConvertToNumpyArray(fingerprint, array)
        #print(fingerprint_list)
        #print(array)

#Convert list of array to a single array using numpy
fingerprint_array_SARS_CoV_2_smiles_filter = np.array(fingerprint_list)    

#Compute UMAP on the SARS_CoV_2 database
umap = umap.UMAP()

umap_fingerprint_array_SARS_CoV_2_smiles_filter = umap.fit_transform(fingerprint_array_SARS_CoV_2_smiles_filter)

#Place the UMAP result into a pandas dataframe for plotting
umap_fingerprint_array_SARS_CoV_2_smiles_filter_fig = pd.DataFrame(umap_fingerprint_array_SARS_CoV_2_smiles_filter,columns=["UMAP2","UMAP3"])
#print(umap_fingerprint_array_SARS_CoV_2_smiles_filter_fig)


umap_fingerprint_array_SARS_CoV_2_smiles_filter_fig.to_pickle("BEST_D_DRUG_ALL.pkl")
output = pd.read_pickle("BEST_D_DRUG_ALL.pkl")
print(output)


#Set the figure for seaborn
sns.set(rc={'figure.figsize': (10, 10)})
sns.set(font_scale=1.5)
#sns.set_style('whitegrid')

#Visualize the UMAP using seabornd=n
fig = sns.jointplot(data=umap_fingerprint_array_SARS_CoV_2_smiles_filter_fig,x="UMAP2",y="UMAP3",kind="hex", height=10, ratio=2, color="blue")

plt.tight_layout()
plt.show()
#plt.savefig('fragments_Morgan_2048.png',dpi=400)

with open("BEST_D_DRUG_ALL.pkl", "rb") as f:
    object = pkl.load(f)
    
df = pd.DataFrame(object)
df.to_csv(r'BEST_D_DRUG_ALL.csv')






