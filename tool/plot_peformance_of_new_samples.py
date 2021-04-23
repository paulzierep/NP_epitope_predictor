import pandas as pd
import os
import numpy as np

import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.metrics import classification_report
from sklearn.metrics import precision_recall_curve

UPDATE_PATH1 = os.path.join('ML_data_updated', 'epitope_update_input', '05-11-2020') 
UPDATE_PATH2 = os.path.join('ML_data_updated', 'epitope_update_input', '10-29-2020') 

##################################
#Load and combine unseen data
##################################

df1 = pd.read_csv(os.path.join(UPDATE_PATH1, 'update_stats.csv'))
df2 = pd.read_csv(os.path.join(UPDATE_PATH2, 'update_stats.csv'))

df = pd.concat([df1,df2])

df = df.dropna()

# print((df.loc[:,'b_proba'] > 0.5).astype(int).value_counts())
# exit()

#################################
#Get AUC plots per cluster
################################

result_df = pd.DataFrame()

for cluster_id, cluster_df in df.groupby('cluster'):

	b_count = cluster_df.loc[(cluster_df['b_cell'] == 1),:].shape[0]
	t_count = cluster_df.loc[(cluster_df['t_cell'] == 1),:].shape[0]
	chebi_count = cluster_df.shape[0]

	if b_count == 0: 
		b_auc = np.NaN
	else:
		b_auc = roc_auc_score(cluster_df.loc[:,'b_cell'], cluster_df.loc[:,'b_proba'])

	if t_count == 0: 
		t_auc = np.NaN
	else:
		t_auc = roc_auc_score(cluster_df.loc[:,'t_cell'], cluster_df.loc[:,'t_proba'])

	result_df.loc[cluster_id, 'b'] = b_count
	result_df.loc[cluster_id, 't'] = t_count
	result_df.loc[cluster_id, 'b_auc'] = b_auc
	result_df.loc[cluster_id, 't_auc'] = t_auc
	result_df.loc[cluster_id, 'chebi'] = chebi_count

# print(result_df)

result_df.loc['all', 'b'] = df.loc[(df['b_cell'] == 1),:].shape[0]
result_df.loc['all', 't'] = df.loc[(df['t_cell'] == 1),:].shape[0]
result_df.loc['all', 'chebi'] = df.shape[0]

b_auc = roc_auc_score(df.loc[:,'b_cell'], df.loc[:,'b_proba'])
t_auc = roc_auc_score(df.loc[:,'t_cell'], df.loc[:,'t_proba'])

result_df.loc['all', 'b_auc'] = b_auc
result_df.loc['all', 't_auc'] = t_auc

result_df.to_csv('validation_dataset_auc.csv')

exit()

#################################
#Get AUC plots
################################

# print(df.loc[:,'b_cell'], df.loc[:,'b_proba'])

# b = roc_auc_score(df.loc[:,'b_cell'], df.loc[:,'b_proba'])
# t = roc_auc_score(df.loc[:,'t_cell'], df.loc[:,'t_proba'])
# print(classification_report(df.loc[:,'b_cell'], (df.loc[:,'b_proba'] > 0.5).astype(int)))
# print(classification_report(df.loc[:,'t_cell'], (df.loc[:,'t_proba'] > 0.5).astype(int)))

# fpr = dict()
# tpr = dict()

# fpr['b_cell'], tpr['b_cell'], _ = roc_curve(df.loc[:,'b_cell'], df.loc[:,'b_proba'])
# plt.plot(fpr['b_cell'], tpr['b_cell'], color='darkorange', label ='B cell (AUC: {0})'.format(round(b,3)))

# fpr['t_cell'], tpr['t_cell'], _ = roc_curve(df.loc[:,'t_cell'], df.loc[:,'t_proba'])
# plt.plot(fpr['t_cell'], tpr['t_cell'], color='green', label ='T cell (AUC: {0})'.format(round(t,3)))

# plt.legend()
# # plt.show()
# plt.savefig('update_AUC_10-29-2020_all.pdf')

# print(b)
# print(t)

##################################
#Get precision recall plot
################################

# precision, recall, _ = precision_recall_curve(df.loc[:,'b_cell'], df.loc[:,'b_proba'])
# plt.step(recall, precision, where='post', color='darkorange')
# # plt.plot(precision, recall, color='darkorange')
# precision, recall, _ = precision_recall_curve(df.loc[:,'t_cell'], df.loc[:,'t_proba'])
# plt.step(recall, precision, where='post', color='green')
# # plt.plot(precision, recall, color='green')

# plt.show()