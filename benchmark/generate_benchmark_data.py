import sys
import os

TOP_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL_PATH = os.path.join(TOP_PATH, 'tool')

sys.path.append(TOOL_PATH)

from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
from output_utils import Results_To_Html, Results_To_Json

import pandas as pd
import os

import benchmark_utils

from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier

import shutil

import numpy as np
np.random.seed(0)

#############################
# Paths
#############################

DATA_PATH = os.path.join(TOOL_PATH, "ML_data")

#############################
# Run benchmark with different FPs
#############################

# chiral_grid = {'T':True,'F':False}
# radius_grid = [0,1,2,3,4,5,6,7,8]

# for chiral in chiral_grid:
# 	for radius in radius_grid:

# 		print('*'*100)
# 		print('Compiling FPs and benchmark with specs: Chiral {0}; Radius {1}'.format(chiral, radius))

# 		#specific folder for each benchmark
# 		BENCHMARK_FOLDER_NAME = 'FP_C{0}_R{1}_MN'.format(chiral, radius)  #Fingerprint, Chiral(C): True, Radius(R): 3, Mode(M): normal
# 		BENCHMARK_DATA = os.path.join('benchmark_data',BENCHMARK_FOLDER_NAME)

# 		#copy the target to generate new FPs for each benchmark run
# 		os.makedirs(BENCHMARK_DATA, exist_ok = True)
# 		SRC_PATH = os.path.join(DATA_PATH, 'target.csv')
# 		DST_PATH = os.path.join(BENCHMARK_DATA, 'target.csv')
# 		shutil.copy(SRC_PATH, DST_PATH)

# 		#generate new FPs
# 		predictor = NP_Epitope_Prediction(data_storage = BENCHMARK_DATA)
# 		#set FP specs
# 		predictor.compile_cluster_FPs(radius = radius, useChirality = chiral_grid[chiral])

# 		#run benchmark
# 		clf = RandomForestClassifier(n_estimators = 100, random_state = 0)
# 		benchmark_utils.evaluate_all_clusters(predictor,clf,BENCHMARK_DATA)

#############################
# Run normal benchmark 
#############################

# DATA_PATH = os.path.join(TOOL_PATH, "ML_data")
# predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

# clf = RandomForestClassifier(n_estimators = 100, random_state = 0)
# data_storage = os.path.join('benchmark_data','FP_CT_R3_MN') #Fingerprint, Chiral: True, Radius: 3
# benchmark_utils.evaluate_all_clusters(predictor,clf, data_storage)

#############################
# Run benchmark with random targets 
#############################

# clf = RandomForestClassifier(n_estimators = 100, random_state = 0)
# data_storage = os.path.join('benchmark_data','FP_CT_R3_MR') #Fingerprint, Chiral: True, Radius: 3, Mode: random
# benchmark_utils.evaluate_all_clusters(predictor,clf, data_storage, mode = 'random')

#############################
# Run normal benchmark different clf
#############################

DATA_PATH = os.path.join(TOOL_PATH, "ML_data")
predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

clf = MLPClassifier(alpha=1, max_iter=1000)
data_storage = os.path.join('benchmark_data','FP_CT_R3_MN_NN') #Fingerprint, Chiral: True, Radius: 3
benchmark_utils.evaluate_all_clusters(predictor,clf, data_storage)

DATA_PATH = os.path.join(TOOL_PATH, "ML_data")
predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)

clf =  KNeighborsClassifier(3)
data_storage = os.path.join('benchmark_data','FP_CT_R3_MN_KNN') #Fingerprint, Chiral: True, Radius: 3
benchmark_utils.evaluate_all_clusters(predictor,clf, data_storage)