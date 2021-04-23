from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
from output_utils import Results_To_Html, Results_To_Json, Results_To_Django

# import sys
# sys.modules['Molecule_Group_Classifier'] = Molecule_Group_Classifier
# sys.modules['Epitope_Predictor'] = Epitope_Predictor
# sys.modules['NP_Epitope_Prediction'] = NP_Epitope_Prediction

import os
import pandas as pd

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(FILE_PATH, "ML_data_updated")

predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

def prediction2django(smiles):

	results = predictor.prediction_chain(smiles, only_epitopes = True, sort_order = "E")
	results = Results_To_Django(results)

	return(results)

##################
# test
##################

# results = prediction2django('CCCCC')
# print(results)
