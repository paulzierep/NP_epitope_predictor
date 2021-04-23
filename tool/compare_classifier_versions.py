'''
Compares different Version of the classifier, gives stats of epitopes and samples
'''

import pandas as pd
import os
import numpy as np

import matplotlib.pyplot as plt


#path of the old training sets
# OLD_PATH = os.path.join('ML_data') 
OLD_PATH = os.path.join('ML_data_updated', 'epitope_update_input', '05-11-2020') 
SMILES_PATH_V1 = os.path.join(OLD_PATH, "chebi_san_assigned.csv")
df_v1 = pd.read_csv(SMILES_PATH_V1, index_col = "index")

#path of the new training sets
UPDATE_PATH = os.path.join('ML_data_updated', 'epitope_update_input', '10-29-2020') 
SMILES_PATH_V2 = os.path.join(UPDATE_PATH, "chebi_san_assigned.csv")
df_v2 = pd.read_csv(SMILES_PATH_V2, index_col = "index")

#some molecules form the old training set can be missing in the new set, due to changes in ChEBI 
#one example is CHEBI:15433: Possilbe reson change in structure (Metallo Organic) leads to redkit parsing error

v1_idx_not_in_v2 = df_v1.index.difference(df_v2.index)
# print(df_v1.loc[v1_idx_not_in_v2,:])

# get a df only new items
v2_idx_without_v1 = df_v2.index.difference(df_v1.index)
v2_new_mols = df_v2.loc[v2_idx_without_v1,:]

#new stats
v2_t = v2_new_mols.loc[(df_v2['t_cell'] == 1),:]
v2_b = v2_new_mols.loc[(df_v2['b_cell'] == 1),:]
v2_chebi = v2_new_mols.loc[((df_v2['b_cell'] == 0) & (df_v2['t_cell'] == 0)),:]

#old stats
v1_t = df_v1.loc[(df_v1['t_cell'] == 1),:]
v1_b = df_v1.loc[(df_v1['b_cell'] == 1),:]
v1_chebi = df_v1.loc[((df_v1['b_cell'] == 0) & (df_v1['t_cell'] == 0)),:]

####################
#plotting bar
####################

# N = 3
# ind = np.arange(N)    # the x locations for the groups

# old_stats = [v2_chebi.shape[0], v2_b.shape[0], v2_t.shape[0]]
# new_stats = [v1_chebi.shape[0], v1_b.shape[0], v1_t.shape[0]]

# print(old_stats)
# print(new_stats)

# p1 = plt.bar(ind, old_stats)
# p2 = plt.bar(ind, new_stats, bottom=old_stats)

# plt.show()

#####################
#plotting pie
#####################

N = 3
ind = np.arange(N)    # the x locations for the groups

titles = ['ChEBI', 'B cell', 'T cell']
old_stats = [v2_chebi.shape[0], v2_b.shape[0], v2_t.shape[0]]
new_stats = [v1_chebi.shape[0], v1_b.shape[0], v1_t.shape[0]]

fig1, axes = plt.subplots(1,3)

axes[0].pie([old_stats[0], new_stats[0]], labels = [old_stats[0], new_stats[0]])
axes[0].title.set_text(titles[0])

axes[1].pie([old_stats[1], new_stats[1]], labels = [old_stats[1], new_stats[1]])
axes[1].title.set_text(titles[1])

axes[2].pie([old_stats[2], new_stats[2]], labels = [old_stats[2], new_stats[2]])
axes[2].title.set_text(titles[2])

plt.show()