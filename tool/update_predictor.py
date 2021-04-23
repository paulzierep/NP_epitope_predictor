import pandas as pd
import os

from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
from data_helper import NP_Epitope_Data_Conversion

#This is where the new data is generated, 
#The update function keeps the cluster model
#(cluster_clf.pickle) as it is.
#Also the ontology mapping of the clusters and 
#the cluster descriptions are kept. They can be updated anytime,
#using the update_onto_mapping and update_cluster_descri method of 
#NP_Epitope_Prediction.
DATA_PATH = "ML_data_updated"

#this is where the input data (IEDB csv and ChEBI sdf) is located
UPDATE_PATH = os.path.join(DATA_PATH, 'epitope_update_input', '10-29-2020') 

#the input data files
SDF_PATH = os.path.join(UPDATE_PATH, "ChEBI_lite_3star.sdf")
B_CELL = os.path.join(UPDATE_PATH, "epitope_table_b_cell_pos.csv")
T_CELL = os.path.join(UPDATE_PATH, "epitope_table_t_cell_pos.csv")

#the NP_Epitope_Data_Conversion will generate this file, which is then passed to the update_chain
SMILES_PATH = os.path.join(UPDATE_PATH, "chebi_san_assigned.csv")

print("Convert input data to NP_Epitope_Prediction training data")
converter = NP_Epitope_Data_Conversion(DATA_PATH)
converter.create_ML_data(SDF_PATH, B_CELL, T_CELL, skip_sdf = False)

print("Initiate NP_Epitope_Prediction class")
predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

print("Update the predictor")
predictor.update_chain(SMILES_PATH)

