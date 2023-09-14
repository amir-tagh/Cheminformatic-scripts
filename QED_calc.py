#!/usr/bin/env python

from qed import qed
from rdkit import Chem
import sys


# Silence non-critical RDKit warnings to minimize unnecessary outputs
from rdkit import RDLogger
#RDLogger.DisableLog('rdApp.*')
#lg = RDLogger.logger()
#lg.setLevel(RDLogger.CRITICAL)


#test = qed.test_data()
#for name in test: print(test[name], name)

'''
test='test.smi'

for name in test:
    mol = Chem.MolFromSmiles(test[name])
    p = qed.properties(mol)
    print("%6.2f\t%6.3f\t%6d\t%6d\t%6.2f\t%6d\t%6d\t%6d\t%6.3f\t%-s" % (p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],qed.default(mol),name))
print("    MW\t ALOGP\t   HBA\t   HBD\t   PSA\t  ROTB\t  AROM\tALERTS\t   QED\tNAME")
'''


data_1 = open(sys.argv[1], 'r')
lines_1 = data_1.readlines()
data_1.close()


for line in lines_1:
    mol = Chem.MolFromSmiles(line)
    if mol is None:
        continue
    #RDLogger.DisableLog('rdApp.*')
    #try:
    p = qed.properties(mol)
    #except:
        #continue
    #print("%6.2f\t%6.3f\t%6d\t%6d\t%6.2f\t%6d\t%6d\t%6d\t%6.3f\t" % (p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],qed.default(mol)))
    print("%6.2f\t%6.3f\t%6d\t%6d\t%6.2f\t%6d\t%6d\t%6d\t%6.3f\t%-s" % (p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],qed.default(mol),line))
print("    MW\t ALOGP\t   HBA\t   HBD\t   PSA\t  ROTB\t  AROM\tALERTS\t   QED\tNAME")
#print("    MW\t ALOGP\t   HBA\t   HBD\t   PSA\t  ROTB\t  AROM\tALERTS\t   QED\t")



'''
from __future__ import print_function
> import rdkit
> from rdkit import Chem
> from rdkit.Chem import AllChem
>
>
> Corefile='Porphyrin.mol'
> core_rdkit_object=Chem.MolFromMolFile(Corefile)
> core_smiles=Chem.MolToSmiles(core_rdkit_object)
> print('core_smiles',core_smiles)
> core=Chem.MolFromSmiles(core_smiles)
> print('rdkit_object',core)
> coreh=Chem.AddHs(core)
> AllChem.EmbedMolecule(coreh)
> print(Chem.MolToMolBlock(coreh))
>
>
> core_smiles_2='C1=C/C2=C/c3ccc4n3[Zn]n3/c(cc/c3=C/C3=N/C(=C\4)C=C3)=C\C1=N2'
> core=Chem.MolFromSmiles(core_smiles_2)
> print(core)
> coreh=Chem.AddHs(core)
> AllChem.EmbedMolecule(coreh)
> print(Chem.MolToMolBlock(coreh))
'''
