import pandas as pd

from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor

DATA_PATH = "ML_data"

print("initiate class")
predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

# print("Onto fold")
# predictor.get_onto_mapping_fold()

# print("Summary")
# predictor.compile_prediction_summary()

#takes approx: 2537.8s
print("FP importance update")
predictor.add_FP_imp_to_classification()

# print("Stats")
# predictor.get_cluster_stats()

# print("Update cluster descri")
# predictor.update_cluster_descri()
