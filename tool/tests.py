import unittest

import os
import pandas as pd
from pandas.util.testing import assert_frame_equal

from similarity_utils import Euclidiean_Distance_Query

#make sure to import the classes which are loaded as pickle
from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor

from output_utils import Results_To_Html, Results_To_Json

#######################
#Example for a simple function test 
#######################

class SimTests(unittest.TestCase):

	def test_Euclidiean_Distance_Query(self):

		target = pd.DataFrame([{"A":1,"B":1,"C":0},{"A":0,"B":0,"C":0},{"A":2,"B":2,"C":2}])
		query = target.iloc[[0],:]

		expected_result_df = pd.DataFrame({
		'A': {0: 1, 1: 0, 2: 2}, 
		'B': {0: 1, 1: 0, 2: 2}, 
		'C': {0: 0, 1: 0, 2: 2}, 
		'E Dist.': {0: 0.0, 1: 1.4142135623730951, 2: 2.449489742783178}
		})

		result_df = Euclidiean_Distance_Query(query,target)
		
		assert_frame_equal(result_df, expected_result_df)


#######################
#Ideally each step of fitting and prediction should be 
#tested with multiple examples (no time)
#######################

# class NP_Epitope_Prediction_Tests(unittest.TestCase):

# 	def set_up(self):

# 		TEST_DATA_PATH = "test_data"
# 		self.NP_E_Pred = NP_Epitope_Prediction(data_storage = TEST_DATA_PATH)

# 	def test_load_target(self):

# 		self.set_up()

# 		smiles_path = os.path.join(TEST_DATA_PATH, "chebi_san_assigned_for_tests.csv")
# 		smiles_df = pd.read_csv(smiles_path, index_col = "index")

# 		self.NP_E_Pred.load_target(smiles_df)

# 		#check if output is where it should be
# 		output_smiles_df = pd.read_csv(self.NP_E_Pred.target_storage, index_col = "index")

# 		assert_frame_equal(smiles_df, output_smiles_df)

	# def test_compile_pre_cluster_FPs(self):

	# 	self.set_up()

	# 	self.NP_E_Pred.compile_pre_cluster_FPs()


class NP_Epitope_Prediction_Tests_Fitted(unittest.TestCase):
	"""
	Basic test for the NP_Epitope_Prediction class, 
	if this test runs successful, the NP_Epitope_Prediction is up and 
	running and can be used. This tests only the already fitted classifier.

	Assuming that the methods:
		* fit_chain
		* update_cluster_descri
		* update_onto_mapping

	did run without errors.
	"""

	def set_up(self):

		TEST_DATA_PATH = "ML_data"
		self.NP_E_Pred = NP_Epitope_Prediction(data_storage = TEST_DATA_PATH)

	def test_json_output(self):

		self.set_up()

		test_smiles = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@H](C(=O)NCCCC[C@H](N)C(=O)O)[C@@H]1N[C@@H](C(=O)O)C(C)(C)S1"

		results = self.NP_E_Pred.prediction_chain(test_smiles, 
			only_epitopes = False, 
			compute_k_best = 5, 
			show_k_best = None, 
			sort_order = "E")

		json_text_result = Results_To_Json(results)

		with open(os.path.join("test_data", "example_output.json"),"r") as json_file:
			json_test_expected = json_file.read()

		self.assertEqual(json_text_result, json_test_expected)

if __name__ == '__main__':
    unittest.main()


############################################
#Create new example file is algo is changed:
############################################

# TEST_DATA_PATH = "ML_data"
# NP_E_Pred = NP_Epitope_Prediction(data_storage = TEST_DATA_PATH)

# test_smiles = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@H](C(=O)NCCCC[C@H](N)C(=O)O)[C@@H]1N[C@@H](C(=O)O)C(C)(C)S1"

# results = self.NP_E_Pred.prediction_chain(test_smiles, 
# 	only_epitopes = False, 
# 	compute_k_best = 5, 
# 	show_k_best = None, 
# 	sort_order = "E")

# json_text_result = Results_To_Json(results)

# with open(os.path.join("test_data", "example_output.json"),"w") as json_file:
# 	json_file.write(json_text_result)