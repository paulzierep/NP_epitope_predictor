from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
from output_utils import Results_To_Html, Results_To_Json, Results_To_Django
from multiprocessing_tools import df_split_apply_multip_combine

import pandas as pd
import os
import numpy as np

import matplotlib.pyplot as plt

##################################
#Get new samples
##################################

#path of the old training sets
# OLD_PATH = os.path.join('ML_data') 
OLD_PATH = os.path.join('ML_data_updated', 'epitope_update_input', '05-11-2020') 
SMILES_PATH_V1 = os.path.join(OLD_PATH, "chebi_san_assigned.csv")
df_v1 = pd.read_csv(SMILES_PATH_V1, index_col = "index")

#path of the new training sets
# UPDATE_PATH = os.path.join('ML_data_updated', 'epitope_update_input', '05-11-2020') 
UPDATE_PATH = os.path.join('ML_data_updated', 'epitope_update_input', '10-29-2020') 
SMILES_PATH_V2 = os.path.join(UPDATE_PATH, "chebi_san_assigned.csv")
df_v2 = pd.read_csv(SMILES_PATH_V2, index_col = "index")

#some molecules form the old training set can be missing in the new set, due to changes in ChEBI 
#one example is CHEBI:15433: Possilbe reson change in structure (Metallo Organic) leads to redkit parsing error

v1_idx_not_in_v2 = df_v1.index.difference(df_v2.index)
# print(df_v1.loc[v1_idx_not_in_v2,:])

# get a df only new items
v2_idx_without_v1 = df_v2.index.difference(df_v1.index)
v2_new_mols = df_v2.loc[v2_idx_without_v1,:]

#selections
v2_t = v2_new_mols.loc[(df_v2['t_cell'] == 1),:]
v2_b = v2_new_mols.loc[(df_v2['b_cell'] == 1),:]
v2_chebi = v2_new_mols.loc[((df_v2['b_cell'] == 0) & (df_v2['t_cell'] == 0)),:]
v2_bt = v2_new_mols.loc[((df_v2['b_cell'] == 1) | (df_v2['t_cell'] == 1)),:]

# print(v2_t.shape[0])
# print(v2_b.shape[0])
# exit()

##################################
#Load old clf
##################################

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(FILE_PATH, "ML_data")

predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

##################################
#Assign proba and cluster to the new samples
##################################

def assing_probas(df):

	for r_index, row in df.iterrows():

		print('*************************************************************')
		print(row['smiles'])

		results = predictor.prediction_chain(row['smiles'], only_epitopes = True, sort_order = "E")

		cluster = results["cluster_info"]["Cluster"]
		b_proba = results["classify_info_b"]['proba']
		t_proba = results["classify_info_t"]['proba']

		df.loc[r_index, 'cluster'] = cluster
		df.loc[r_index, 'b_proba'] = b_proba
		df.loc[r_index, 't_proba'] = t_proba

	return(df)

v2_new_mols_assigned = assing_probas(v2_new_mols)

# v2_new_mols_assigned = df_split_apply_multip_combine(v2_new_mols, assing_probas, None,  num_of_processes = 20)

v2_new_mols_assigned.to_csv(os.path.join(UPDATE_PATH, 'update_stats.csv'))

