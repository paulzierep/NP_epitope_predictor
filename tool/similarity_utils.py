from scipy import spatial
import pandas as pd
import os
import numpy as np

from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import Chem

def Euclidiean_Distance_Query(query, target, columns = None, k_best = None):

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

def RDKit_Mol_Similarity(query, target):

    q = Chem.MolFromSmiles(query)
    t = Chem.MolFromSmiles(target)

    fp_q = AllChem.GetMorganFingerprint(
                     q,
                     3, 
                     useFeatures = False, 
                     useChirality = True,
                     )

    fp_t = AllChem.GetMorganFingerprint(
                     t,
                     3, 
                     useFeatures = False, 
                     useChirality = True,
                     )

    sim = DataStructs.TanimotoSimilarity(fp_q,fp_t)
    #sim = DataStructs.DiceSimilarity(fp_q,fp_t)

    return(sim)

def Add_ChEBI_IP(val):
    return('<a target="_blank" href="https://www.ebi.ac.uk/chebi/searchId.do?chebiId={0}">{0}</a>'.format(val))

#def add_target_class(val):

# DATA_PATH = "ML_data"
# FP_PATH = os.path.join(DATA_PATH, "unfolded_fps")
# EXAMPLE_FP_PATH = os.path.join(FP_PATH, "0.csv")
# X = pd.read_csv(EXAMPLE_FP_PATH, index_col = "index")

# target = pd.DataFrame([{"A":1,"B":1,"C":0},{"A":0,"B":0,"C":0},{"A":2,"B":2,"C":2}])
# query = target.iloc[[0],:]

# res = Euclidiean_Distance_Query(query,target)#["5685888"])

# print(res.to_dict())

#print(X.loc[ary==ary.min(),:])

#print(ary)
