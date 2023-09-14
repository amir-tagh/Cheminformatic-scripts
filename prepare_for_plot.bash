#!/bin/bash

#26ExactMolWt 32HeavyAtomCount 38LabuteASA

#rm -rf *.txt

######topological descriptors#######
#1BalabanJ
#2BertzCT A topological index meant to quantify “complexity” of molecules
#3Chi0 #4Chi0n #5Chi0v #6Chi1 #7Chi1n

#30 FractionCSP3
#32 HeavyAtomCount
#38 LabuteASA
#40 MaxAbsPartialCharge
#47 MolLogP
#48 MolMR
#49 MolWt
#51 NOCount
#57 NumAromaticRings
#58 NumHAcceptors
#59 NumHDonors
#62 NumRotatableBonds
#66 NumValenceElectrons
#118 fr_ArN
#NumRadicalElectrons

#./RDkit_mol_desc_calc.py ROBIN_RNA_Binders_just_smi.smi > ROBIN_desc.txt

#sed -i "s/[][]//g" drugbank-mol-desc.dat
#sed -i "s/\,/\t/g" drugbank-mol-desc.dat

#sed -i "s/[][]//g" DEL_LIB-mol-desc.dat
#sed -i "s/\,/\t/g" DEL_LIB-mol-desc.dat

#sed -i "s/[][]//g" HITS-mol-desc.dat 
#sed -i "s/\,/\t/g" HITS-mol-desc.dat


counter=0
files=('BalabanJ' 'BertzCT' 'Chi0' 'Chi0n' 'Chi0v' 'Chi1' 'Chi1n' 'Chi1v' 'Chi2n' 'Chi2v' 'Chi3n' 'Chi3v' 'Chi4n' 'Chi4v' 'EState_VSA1' 'EState_VSA10' 'EState_VSA11' 'EState_VSA2' 'EState_VSA3' 'EState_VSA4' 'EState_VSA5' 'EState_VSA6' 'EState_VSA7' 'EState_VSA8' 'EState_VSA9' 'ExactMolWt' 'FpDensityMorgan1' 'FpDensityMorgan2' 'FpDensityMorgan3' 'FractionCSP3' 'HallKierAlpha' 'HeavyAtomCount' 'HeavyAtomMolWt' 'Ipc' 'Kappa1' 'Kappa2' 'Kappa3' 'LabuteASA' 'MaxAbsEStateIndex' 'MaxAbsPartialCharge' 'MaxEStateIndex' 'MaxPartialCharge' 'MinAbsEStateIndex' 'MinAbsPartialCharge' 'MinEStateIndex' 'MinPartialCharge' 'MolLogP' 'MolMR' 'MolWt' 'NHOHCount' 'NOCount' 'NumAliphaticCarbocycles' 'NumAliphaticHeterocycles' 'NumAliphaticRings' 'NumAromaticCarbocycles' 'NumAromaticHeterocycles' 'NumAromaticRings' 'NumHAcceptors' 'NumHDonors' 'NumHeteroatoms' 'NumRadicalElectrons' 'NumRotatableBonds' 'NumSaturatedCarbocycles' 'NumSaturatedHeterocycles' 'NumSaturatedRings' 'NumValenceElectrons' 'PEOE_VSA1' 'PEOE_VSA10' 'PEOE_VSA11' 'PEOE_VSA12' 'PEOE_VSA13' 'PEOE_VSA14' 'PEOE_VSA2' 'PEOE_VSA3' 'PEOE_VSA4' 'PEOE_VSA5' 'PEOE_VSA6' 'PEOE_VSA7' 'PEOE_VSA8' 'PEOE_VSA9' 'RingCount' 'SMR_VSA1' 'SMR_VSA10' 'SMR_VSA2' 'SMR_VSA3' 'SMR_VSA4' 'SMR_VSA5' 'SMR_VSA6' 'SMR_VSA7' 'SMR_VSA8' 'SMR_VSA9' 'SlogP_VSA1' 'SlogP_VSA10' 'SlogP_VSA11' 'SlogP_VSA12' 'SlogP_VSA2' 'SlogP_VSA3' 'SlogP_VSA4' 'SlogP_VSA5' 'SlogP_VSA6' 'SlogP_VSA7' 'SlogP_VSA8' 'SlogP_VSA9' 'TPSA' 'VSA_EState1' 'VSA_EState10' 'VSA_EState2' 'VSA_EState3' 'VSA_EState4' 'VSA_EState5' 'VSA_EState6' 'VSA_EState7' 'VSA_EState8' 'VSA_EState9' 'fr_Al_COO' 'fr_Al_OH' 'fr_Al_OH_noTert' 'fr_ArN' 'fr_Ar_COO' 'fr_Ar_N' 'fr_Ar_NH' 'fr_Ar_OH' 'fr_COO' 'fr_COO2' 'fr_C_O' 'fr_C_O_noCOO' 'fr_C_S' 'fr_HOCCN' 'fr_Imine' 'fr_NH0' 'fr_NH1' 'fr_NH2' 'fr_N_O' 'fr_Ndealkylation1' 'fr_Ndealkylation2' 'fr_Nhpyrrole' 'fr_SH' 'fr_aldehyde' 'fr_alkyl_carbamate' 'fr_alkyl_halide' 'fr_allylic_oxid' 'fr_amide' 'fr_amidine' 'fr_aniline' 'fr_aryl_methyl' 'fr_azide' 'fr_azo' 'fr_barbitur' 'fr_benzene' 'fr_benzodiazepine' 'fr_bicyclic' 'fr_diazo' 'fr_dihydropyridine' 'fr_epoxide' 'fr_ester' 'fr_ether' 'fr_furan' 'fr_guanido' 'fr_halogen' 'fr_hdrzine' 'fr_hdrzone' 'fr_imidazole' 'fr_imide' 'fr_isocyan' 'fr_isothiocyan' 'fr_ketone' 'fr_ketone_Topliss' 'fr_lactam' 'fr_lactone' 'fr_methoxy' 'fr_morpholine' 'fr_nitrile' 'fr_nitro' 'fr_nitro_arom' 'fr_nitro_arom_nonortho' 'fr_nitroso' 'fr_oxazole' 'fr_oxime' 'fr_para_hydroxylation' 'fr_phenol' 'fr_phenol_noOrthoHbond' 'fr_phos_acid' 'fr_phos_ester' 'fr_piperdine' 'fr_piperzine' 'fr_priamide' 'fr_prisulfonamd' 'fr_pyridine' 'fr_quatN' 'fr_sulfide' 'fr_sulfonamd' 'fr_sulfone' 'fr_term_acetylene' 'fr_tetrazole' 'fr_thiazole' 'fr_thiocyan' 'fr_thiophene' 'fr_unbrch_alkane' 'fr_urea' 'qed')

#files=('BalabanJ')

for f in "${files[@]}";
 do	
  counter=$((counter+1))
   echo "string is $f"
   echo "counter is $counter"

    var1=$counter

    #awk -v var=$var1  '{print $var}' drugbank-mol-desc.dat > drugbank-desc.dat
     awk -v var=$var1  '{print $var}' DEL_LIB-mol-desc.dat > DELmol-desc.dat
      awk -v var=$var1  '{print $var}' HITS-mol-desc.dat > HITS-desc.dat

      sed -i "s/BBB/$f/g" violin_polt_QED.py

       ./violin_polt_QED.py DELmol-desc.dat HITS-desc.dat

         sed -i "s/$f/BBB/g" violin_polt_QED.py
done

#rm -rf Images_RDKIT
mkdir -p Images_RDKIT_2image
mv *.png Images_RDKIT_2image



:<<'END'
for i in 7
 do
  #awk '{print $30}' DEL_HITS_MOl_Desc.dat > inforna_ExactMolWt.txt
   awk '{print $7}' APNOHASH_desc.txt > APNOHASH_mol_desc.dat
    awk  '{print $7}' ROBIN_desc.txt > ROBIN_mol_desc.dat

done

#remove compounds with MW > 600 from inforna
#awk -F'[ =]' '$1 < 500' inforna_ExactMolWt.txt > inforna_ExactMolWt_edit.txt

sed -i "s/Chi1/Chi1n/g" violin_polt_QED.py

./violin_polt_QED.py APNOHASH_mol_desc.dat ROBIN_mol_desc.dat
END
