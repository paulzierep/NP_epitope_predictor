from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
from output_utils import Results_To_Html, Results_To_Json

import pandas as pd
import os

#############################
#Run specific methods 
#############################

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(FILE_PATH, "ML_data_updated")

predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)
predictor.update_epitope_prediction_info()

#############################
#Retrain the classifier
#############################

# S_PATH = "ML_data_tests/epitope_update/11.05.2020"
# SMILES_PATH = os.path.join(S_PATH, "chebi_san_assigned.csv")

# #load data
# smiles_df = pd.read_csv(SMILES_PATH, index_col = "index")

# #print(smiles_df)
# DATA_PATH = "ML_data_tests"
# predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

# #refits the entire predictor, need reassignment of molecular groups and ontos
# # predictor.fit_chain(smiles_df)

# #update the prdictor (basically a refit, but the clustering is based on the pre-fitted clf)

# # predictor.predict_clustering_bulk()
# predictor.update_chain(smiles_df) #1595.4s

# # predictor.fit_classification(cell_type = "b_cell")
# # predictor.fit_classification(cell_type = "t_cell")
# # predictor.add_FP_imp_to_classification()

#############################
#Run a prdiction
#############################

# # DATA_PATH = "ML_data"

# # smiles = "CC(=O)NC1[C@@H](OC(C(C1O)O)CO)O"

# # smiles = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@H](C(=O)NCCCC[C@H](N)C(=O)O)[C@@H]1N[C@@H](C(=O)O)C(C)(C)S1"

# # smiles = "COC(=O)CCc1c(C)c2=CC3=[N]4C(=Cc5c(C=C)c(C)c6C=C7C(C)=C(CCC(O)=O)C8=[N]7[Mg]4(n56)n2c1=C8)C(C)=C3C=C"

# # smiles = "CCCCCCCCCCCCC\C=C\[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O[C@@H]4O[C@H](CO)[C@H](O)[C@H](O[C@@]5(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O5)[C@H](O)[C@@H](CO)O[C@@]5(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O5)[C@H](O)[C@H](O)CO)C(O)=O)C(O)=O)[C@H]4O)[C@H]3NC(C)=O)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O)[C@@H](CO)O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O)[C@H](O)CO)C(O)=O)C(O)=O)[C@H]2O)[C@H](O)[C@H]1O)NC([*])=O"

# smiles = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@H](C(=O)NCCCC[C@H](N)C(=O)O)[C@@H]1N[C@@H](C(=O)O)C(C)(C)S1"

# predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

# results = predictor.prediction_chain(smiles, only_epitopes = True, sort_order = "E")

# #generate a json output
# json_text = Results_To_Json(results)

# with open("example.json", "w") as json_file:
#  	json_file.write(json_text)

# # #generate a html output
# # #this is just an example
# # #it would probably be much better and nices to pass the 
# # #json output to a framework like Flask or Django
# html_text = Results_To_Html(results)

# with open("example.html", "w") as html_file:
#  	html_file.write(html_text)
