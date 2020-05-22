import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw

df = pd.read_csv('target.csv')

# print(df)
# exit()

#get only substructures for positive and cluster 4
cluster_sel = df.loc[(df['clusters'] == 4),:]
cluster_sel = cluster_sel.loc[(cluster_sel['b_cell'] == True),:]

runs = 1000

features_df = pd.DataFrame()

for i in range(runs):

	print(i)
	cluster_sample = cluster_sel.sample(frac=0.01)

	# print(cluster_sel['smiles'])
	mols = [Chem.MolFromSmiles(smi) for smi in cluster_sample['smiles']]

	MCS = rdFMCS.FindMCS(mols,timeout=1, threshold = 0.1, matchValences=True, completeRingsOnly = True)

	smart = MCS.smartsString

	features_df.loc[i,'smart'] = smart 

	mol = Chem.MolFromSmarts(smart)
	Draw.MolToFile(mol,'figures/{0}.png'.format(i))

features_df.to_csv('features.csv')