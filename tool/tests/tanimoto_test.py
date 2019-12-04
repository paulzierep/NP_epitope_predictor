from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

smiles1 = 'CCOP(=O)(O)OP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cc(C)c(=O)[nH]c2=O)C[C@@H]1O'
smiles2 = 'CCOP(=O)(O)OP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O'
smiles3 = 'CO[C@@H]1[C@@H](COP(C)(=O)O)O[C@@H](n2cnc3c(N(C)C)ncnc32)[C@@H]1O'

mol1 = Chem.MolFromSmiles(smiles1)
mol2 = Chem.MolFromSmiles(smiles2)
mol3 = Chem.MolFromSmiles(smiles3)

fp1 = AllChem.GetMorganFingerprint(
                 mol1,
                 3, 
                 #nBits = 1024,
                 useFeatures = False, 
                 useChirality = True,
                 )

fp2 = AllChem.GetMorganFingerprint(
                 mol2,
                 3, 
                 #nBits = 1024,
                 useFeatures = False, 
                 useChirality = True,
                 )

fp3 = AllChem.GetMorganFingerprint(
                 mol3,
                 3, 
                 #nBits = 1024,
                 useFeatures = False, 
                 useChirality = True,
                 )

print(DataStructs.TanimotoSimilarity(fp1, fp2))
print(DataStructs.TanimotoSimilarity(fp1, fp3))
print(DataStructs.TanimotoSimilarity(fp2, fp3))