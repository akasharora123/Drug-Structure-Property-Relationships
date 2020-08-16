#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
Code for creating the follwoing molecule 
features given the input SMILES:
   -- Morgan fp
   -- RDKit fo
   -- RDKit Descriptors
   -- MACCS Keys
   -- Customized fingerprints (ex. singlets, doublets):
      to be read from a file.
"""

import os
import numpy as np
import scipy as sp
import pandas as pd
import inspect
from collections import OrderedDict
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import RDKFingerprint
import rdkit.Chem.Descriptors as Descriptors
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import MACCSkeys

class featurizer(object):
    
   def __init__(self, df):
       self.df = df
       self.nmol = len(df.index)
       self.df['mol'] = self.df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x))


   def get_rdkit_fp(self, minP, maxP, nbits, prefix):
   
       rdkit_df = pd.DataFrame(self.df['SMILES'])
       rdkit_matrix = np.zeros((1,nbits))
       keys = [prefix+"%d"%(i) for i in range(0,nbits)]
   
       for i in range(self.nmol):
          try:
             fp = Chem.RDKFingerprint(self.df['mol'][i], minPath=minP,maxPath=maxP,fpSize=nbits) 
             rdkit_matrix = np.vstack((rdkit_matrix, fp))
   
          except:
             print('Problem calculating rdkit fp for:{}'.format(self.df['SMILES'][i]))
              
       rdkit_matrix = np.delete(rdkit_matrix, 0, axis=0)
       print('RDKit fingerprint matrix dimensions:{}'.format(rdkit_matrix.shape)) 
       rdkit_df = pd.concat([rdkit_df, pd.DataFrame(rdkit_matrix, columns=keys)], axis=1)
       
       return rdkit_df


   def get_Morgan_fp(self, radius, nbits, prefix):
   
       morgan_df = pd.DataFrame(self.df['SMILES'])
       morgan_matrix = np.zeros((1,nbits))
       keys = [prefix+"%d"%(i) for i in range(0,nbits)]
   
       for i in range(self.nmol):
          try:
             fp = Chem.AllChem.GetMorganFingerprintAsBitVect(self.df['mol'][i], radius, nBits=nbits)
             morgan_matrix = np.vstack((morgan_matrix, fp))
   
          except:
             print('Problem calculating Morgan fp for:{}'.format(self.df['SMILES'][i]))
   
       morgan_matrix = np.delete(morgan_matrix, 0, axis=0)
       print('Morgan fingerprint matrix dimensions:{}'.format(morgan_matrix.shape))
       morgan_df = pd.concat([morgan_df, pd.DataFrame(morgan_matrix, columns=keys)], axis=1)
    
       return morgan_df
   
   
   def get_rdkit_descriptors(self, prefix):
   
       getonly=['NHOHCount', 'NOCount', 'NumAliphaticCarbocycles',
                'NumAliphaticHeterocycles', 'NumAliphaticRings',
                'NumAromaticCarbocycles', 'NumAromaticHeterocycles',
                'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 
                'NumHeteroatoms', 'NumRadicalElectrons', 'NumRotatableBonds',
                'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles', 
                'NumSaturatedRings', 'NumValenceElectrons',
                'qed','TPSA', 'MolMR','BalabanJ', 'BertzCT',
                'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO',
                'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 
                'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 'fr_Imine',
                'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1',
                'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde',
                'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide',
                'fr_amidine', 'fr_aniline', 'fr_aryl_methyl', 'fr_azide', 'fr_azo',
                'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine', 'fr_bicyclic',
                'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether',
                'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone', 
                'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone',
                'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine',
                'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho', 
                'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 
                'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 
                'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd',
                'fr_pyridine', 'fr_quatN', 'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone',
                'fr_term_acetylene', 'fr_tetrazole', 'fr_thiazole', 'fr_thiocyan', 
                'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea',
                'MolWt','MolLogP']
       
       calc_props = OrderedDict(inspect.getmembers(Descriptors, inspect.isfunction))
       rdkit_spec = pd.DataFrame(self.df['SMILES'])
    
       for key in list(calc_props.keys()):
           if key.startswith('_'):
              del calc_props[key]
           
       for key,val in calc_props.items():
           if key in getonly:
              rdkit_spec[prefix+key] = self.df['mol'].apply(val)
       
       print('No. of specific descriptors from RDKIT = {:}'.format(len(calc_props)))
    
       return rdkit_spec


   def get_MACCS_fp(self, prefix):
   
       maccs_df = pd.DataFrame(self.df['SMILES'])
       fp = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles('CC')) #A dummy fp for knowing length
       maccs_matrix = np.zeros((1,len(fp)))
       keys = [prefix+"%d"%(i) for i in range(0,len(fp))]
   
       for i in range(self.nmol):
          try:
             fp = MACCSkeys.GenMACCSKeys(self.df['mol'][i])
             maccs_matrix = np.vstack((maccs_matrix, fp))
   
          except:
             print('Problem calculating MACCS key for:{}'.format(self.df['SMILES'][i]))
              
       maccs_matrix = np.delete(maccs_matrix, 0, axis=0)
       print('No. of MACCS keys:{}'.format(len(fp))) 
       maccs_df = pd.concat([maccs_df, pd.DataFrame(maccs_matrix, columns=keys)], axis=1)
       
       return maccs_df


   # Get customized features reading from file:
   # fname = vectors of features corresponding to SMILES/BigSMILES in df
   def get_custom_features(self, filename, prefix):
   
       vecs = pd.read_csv(filename, header=None)
       vecs = vecs.T
       vecs.columns = [prefix+"%d"%(i) for i in range(1,len(vecs.columns)+1)]
       vecs.insert(0,'SMILES',self.df['SMILES'])
       
       return vecs 
 

 
