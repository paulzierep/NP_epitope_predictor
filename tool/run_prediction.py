from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
from output_utils import Results_To_Html, Results_To_Json

DATA_PATH = "ML_data"

smiles = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@H](C(=O)NCCCC[C@H](N)C(=O)O)[C@@H]1N[C@@H](C(=O)O)C(C)(C)S1"

predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

results = predictor.prediction_chain(smiles, only_epitopes = False, compute_k_best = 5, show_k_best = None, sort_order = "E")

json_text = Results_To_Json(results)

print(json_text)