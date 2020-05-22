import pandas as pd
from rdkit import Chem

df = pd.read_csv('features.csv', index_col = 0)
target_df = pd.read_csv('target.csv')
cluster_sel = target_df.loc[(target_df['clusters'] == 4),:]
# cluster_sel = cluster_sel.iloc[:10]

# cluster_sel = target_df

features = pd.DataFrame()
features['smarts'] = df['smart'].value_counts().index

X = pd.DataFrame()
for row_id, row in cluster_sel.iterrows():
	print(row_id)
	mol = Chem.MolFromSmiles(row['smiles'])

	for f_id, f_row in features.iterrows():
		smart = f_row['smarts']
		matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts(smart)))
		X.loc[row_id, f_id] = matches

X.to_csv('X.csv')