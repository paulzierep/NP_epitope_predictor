from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
from output_utils import Results_To_Html, Results_To_Json

import pandas as pd

DATA_PATH = "ML_data"

smiles = "CC(=O)NC1[C@@H](OC(C(C1O)O)CO)O"

#smiles = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@H](C(=O)NCCCC[C@H](N)C(=O)O)[C@@H]1N[C@@H](C(=O)O)C(C)(C)S1"

#smiles = "COC(=O)CCc1c(C)c2=CC3=[N]4C(=Cc5c(C=C)c(C)c6C=C7C(C)=C(CCC(O)=O)C8=[N]7[Mg]4(n56)n2c1=C8)C(C)=C3C=C"

predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)


results = predictor.prediction_chain(smiles, only_epitopes = True, sort_order = "E")

# for k,v in results.items():
# 	print(k, type(v))
# 	if type(v) == dict:
# 		for k2,v2 in v.items():
# 			print("\t", k2, type(v2))

# exit()
#generate a json output
json_text = Results_To_Json(results)

with open("example.json", "w") as json_file:
 	json_file.write(json_text)

#generate a html output
#this is just an example
#it would probably be much better and nices to pass the 
#json output to a framework like Flask or Django
html_text = Results_To_Html(results)

with open("example.html", "w") as html_file:
 	html_file.write(html_text)
