from scipy import spatial
import pandas as pd
import os
import numpy as np

from sklearn.metrics import jaccard_score
from sklearn.metrics import pairwise_distances

def euclidiean_distance_query(query, target, columns = None, k_best = None):

	if not k_best:
		k_best = target.shape[0]

	#only trim df if columns are provided
	if columns:

		#should be a one row df
		query = query.loc[:,columns]

		#can be multi row df
		target = target.loc[:,columns]

	#distance of each row (flatten with [0])
	distance_array = spatial.distance.cdist(query, target, metric='euclidean')[0]


	#get best indices
	top_idx = distance_array.argsort()[:k_best]

	distances = distance_array[top_idx]

	#select from df
	results = target.iloc[top_idx,:].copy()
	results["E Dist."] = distances

	#results["IDs"] = results.index

	#results.index = ["Target {0}".format(x) for x in range(k_best)]

	return(results)


def jaccard_distance_query(query, target, columns = None, k_best = None):

	if not k_best:
		k_best = target.shape[0]

	#only trim df if columns are provided
	if columns:

		#should be a one row df
		query = query.loc[:,columns]

		#can be multi row df
		target = target.loc[:,columns]

	distance_array = target.apply(lambda row: jaccard_score(row, query.iloc[0], average = 'micro'), axis = 1)

	results = target.copy()

	results["J Dist."] = distance_array

	#results.index = ["Target {0}".format(x) for x in range(k_best)]

	return(results)

DATA_PATH = "../ML_data"
FP_PATH = os.path.join(DATA_PATH, "unfolded_fps")
EXAMPLE_FP_PATH = os.path.join(FP_PATH, "6.csv")
X = pd.read_csv(EXAMPLE_FP_PATH, index_col = "index")

res = euclidiean_distance_query(X.iloc[[0],:3], X.iloc[:3,:3], )#["5685888"])
print(res)

res = jaccard_distance_query(X.iloc[[0],:3], X.iloc[:3,:3], )#["5685888"])
print(res)


#print(X.loc[ary==ary.min(),:])

#print(ary)
