import os
import shutil
import pickle

import pandas as pd

from sklearn.cluster import KMeans

#pd.set_option('display.max_colwidth', -1)

DATA_PATH = "../ML_data"
target = os.path.join(DATA_PATH, "target.csv")
fp = os.path.join(DATA_PATH, "fp_folded.csv")

y = pd.read_csv(target, index_col = "index")

y = y.loc[(y["b_cell"] == 1),"smiles"].sample(10).to_list()

print(y)

exit()

#X = pd.read_csv(fp, index_col = "index")

#X = X.astype(float)


# y2 = pd.read_csv("target_df_pca_c_oc.csv", index_col = "acc")
#X2 = pd.read_csv("fps_oc.csv", index_col = "acc")

#print((X == X2).all())

#test_row = (X.loc[:,"12"] == X2.loc[:,"12"])

#print(test_row.loc[(test_row == False)])

#print(X.loc["CHEBI:51736", "12"])
#print(X2.loc["CHEBI:51736", "12"])


###########################
#Test the FP generation
###########################

# from rdkit import Chem
# from rdkit.Chem import AllChem, DataStructs

# def convert_to_count(val):
#     if val == 0:
#         return(0)
#     else:
#         return(1)

# def smiles_df_to_fp_df_folded(smiles_df, 
#                                 OUTPATH = None, 
#                                 count_vector = False, 
#                                 nBits = 1024, 
#                                 radius = 3, 
#                                 useFeatures = False, 
#                                 useChirality = True, 
#                                 verbose = False,
#                                 **kwargs):

#     result_df = pd.DataFrame(columns = list(range(nBits)))

#     for index, (mol_id, row) in enumerate(smiles_df.iterrows()):

#         #show process 
#         if verbose:
#             if index % 1000 == 0:
#                 print(index)

#         mol = Chem.MolFromSmiles(row["smiles"])

#         fp = AllChem.GetHashedMorganFingerprint(
#                              mol,
#                              radius, 
#                              nBits = nBits,
#                              useFeatures = useFeatures, 
#                              useChirality = useChirality,
#                              ).GetNonzeroElements()

#         result_df.loc[mol_id,:] = 0
#         result_df.loc[mol_id, fp.keys()] = fp.values()

#     #simple way to convert count to bit
#     if not count_vector:
#         result_df = result_df.applymap(convert_to_count)

#     #convert to int (saves memory)
#     result_df = result_df.applymap(int)

#     result_df.index.name = "index"

#     if OUTPATH:
#         result_df.to_csv(OUTPATH)

#     return(result_df)

# smi1 = y.loc[["CHEBI:51736"],["smiles"]]

# print(smi1)

# fp1 =  smiles_df_to_fp_df_folded(smi1)

# print(fp1.loc[:,12])

# # mol1 = Chem.MolFromSmiles(smi1)
# # smi1 = Chem.MolToSmiles(mol1)

# # print(smi1)

# smi2 = y2.loc[["CHEBI:51736"],["smiles"]]

# print(smi2)

# fp2 =  smiles_df_to_fp_df_folded(smi2)

# print(fp2.loc[:,12])

# exit()

# print("C1C2*34567C1~C3~C4~C5~C6~C7~2")

# print(y.loc["CHEBI:51736",:])
# print(y2.loc["CHEBI:51736",:])

# exit()









exit()


print(X.shape)

X.sort_index(inplace = True)
#X = X.iloc[:1000,:]

print(X.iloc[:5,:])

cluster_clf = KMeans(n_clusters=8, random_state=0, n_jobs = -1)
clusters = cluster_clf.fit_predict(X)

print(clusters)

# cluster_clf = KMeans(n_clusters=8, random_state=0, n_jobs = -1)
# clusters = cluster_clf.fit_predict(X.copy())

# print(clusters)

exit()

# exit()
#print(clusters)

y2 = pd.read_csv("target_df_pca_c_oc.csv", )
X2 = pd.read_csv("fps_oc.csv", index_col = "acc")

print(X2.shape)

X2.sort_index(inplace = True)
#X2 = X2.iloc[:1000,:]

print(X2.iloc[:5,:])

cluster_clf = KMeans(n_clusters=8, random_state=0, n_jobs = -1)
clusters2 = cluster_clf.fit_predict(X2)

print(clusters2)

print((clusters2 == clusters).all())

"""
                0    1    2    3    4    5  ...  1018  1019  1020  1021  1022  1023
acc                                         ...                                    
CHEBI:100     0.0  0.0  0.0  1.0  0.0  0.0  ...  0.0   1.0   0.0   0.0   0.0   0.0 
CHEBI:10002   0.0  0.0  0.0  0.0  0.0  0.0  ...  0.0   0.0   0.0   0.0   0.0   0.0 
CHEBI:10003   0.0  0.0  0.0  1.0  0.0  0.0  ...  0.0   1.0   0.0   0.0   1.0   0.0 
CHEBI:100147  0.0  0.0  0.0  0.0  0.0  0.0  ...  0.0   0.0   0.0   0.0   0.0   0.0 
CHEBI:100148  0.0  0.0  0.0  0.0  0.0  0.0  ...  0.0   0.0   0.0   0.0   0.0   0.0 
"""

"""
[2 2 7 ... 6 2 6]
"""