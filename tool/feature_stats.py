import os
import pandas as pd

from statsmodels.stats.multitest import multipletests
from scipy import stats

import numpy as np

from rdkit_drawing import DrawBitFromSmiles

# import matplotlib.pyplot as plt

from sklearn.feature_selection import SelectKBest, chi2

#pd.reset_option('^display.', silent=True)
pd.set_option('display.precision',2)
pd.set_option('display.max_columns', None)
pd.set_option('display.float_format', lambda x: '%.5f' % x)

def fisher_exact_for_continuous_outcome(p, n):

    p = np.array(p)
    n = np.array(n)

    p_has_feature = np.count_nonzero(p)
    p_no_feature = np.count_nonzero(p == 0)

    n_has_feature = np.count_nonzero(n)
    n_no_feature = np.count_nonzero(n == 0)

    oddsratio, pvalue = stats.fisher_exact([[p_has_feature, n_has_feature], [p_no_feature, n_no_feature]])

    return(pvalue)

def chi2_undefined_categorical_observation(p, n):

    p = pd.Series(p, name = 'observed')
    n = pd.Series(n, name = 'expected')

    p_counts = p.value_counts()
    n_counts = n.value_counts()

    # print(p_counts)
    # print(n_counts)

    final_df = pd.concat([p_counts, n_counts], axis = 1).T

    final_df = final_df.fillna(0)

    # print(final_df)

    chi2, p_val, dof, expected = stats.chi2_contingency(final_df)

    return(p_val)

def Remove_Correlated_Features(X, threshold = 0.8):

    """
    Removes highly correlated features
    """
    
    # Create correlation matrix
    # https://stackoverflow.com/questions/48270953/pandas-corr-and-corrwith-very-slow
    # use np vectorization and avoid Nan lookup !
    # old: corr_matrix = X.corr().abs()
    arr = np.corrcoef(X.values, rowvar=False)
    corr_df = pd.DataFrame(arr, columns = X.columns, index = X.columns)
    corr_matrix = corr_df.abs()

    # Select upper triangle of correlation matrix
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))

    # Find index of feature columns with correlation greater than threshold
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]

    # Drop features 
    X = X.drop(X[to_drop], axis=1)

    return(X)

def Compute_Feature_Enrichment_DF(X, y):
    
    X = Remove_Correlated_Features(X)

    enrichment_df = pd.DataFrame()

    X_p = X.loc[y == 1,:]
    X_n = X.loc[y == 0,:]

    selector = SelectKBest(chi2, k = 8)

    selector.fit(X, y)    

    enrichment_df['scores'] = selector.scores_
    enrichment_df['p'] = selector.pvalues_

    enrichment_df.index = X.columns

    enrichment_df = enrichment_df.sort_values(by = 'p', ascending = True)

    # #corrected p
    _, enrichment_df["p corr."], _ , _ = multipletests(enrichment_df["p"],  method="bonferroni")

    # #how many samples do have the FP
    enrichment_df["Samples (P)"] = X_p.astype(bool).sum(axis=0)
    enrichment_df["Samples (N)"] = X_n.astype(bool).sum(axis=0)

    
    # #mean for those samples which do have a fingerprint
    enrichment_df["Mean (P)"] = X_p[X_p != 0].mean()
    enrichment_df["Mean (N)"] = X_n[X_n != 0].mean()

    enrichment_df["Samples (P) %"] = (enrichment_df["Samples (P)"] / X_p.shape[0]).round(2)
    enrichment_df["Samples (N) %"] = (enrichment_df["Samples (N)"] / X_n.shape[0]).round(2)

    enrichment_df["Fold"] = enrichment_df["Samples (P) %"] / enrichment_df["Samples (N) %"]
    enrichment_df["Mean diff."] = enrichment_df["Mean (P)"] - enrichment_df["Mean (N)"]
    
    enrichment_mod = enrichment_df.loc[:,['p corr.','Samples (P) %','Fold', "Mean diff."]]

    return(enrichment_mod)

###############################
#Test
################################

def get_figures_of_FPs(ids, iter_smiles, outputfolder):

    '''
    ids: list of fps ids to get figures for FPs

    iter_smiles : list of possible smiles which lead to that FP,
    '''

    # print(iter_smiles)

    id_svg_mapper = {}

    for FP_id in ids:

        for smiles in iter_smiles:
            # print(smiles)

            fp_svg = DrawBitFromSmiles(smiles, FP_id, weboutput = False)

            if fp_svg:
                id_svg_mapper[FP_id] = fp_svg
                break

    if outputfolder:

        os.makedirs(outputfolder, exist_ok = True)

        for fp_id, svg in id_svg_mapper.items():
            out_path = os.path.join(outputfolder, '{0}.svg'.format(fp_id))
            with open(out_path, 'w') as svg_img:
                svg_img.write(svg)

    else:
        return(id_svg_mapper)



# DATA_PATH = "ML_data"
# FP_PATH = os.path.join(DATA_PATH, "unfolded_fps")
# EXAMPLE_FP_PATH = os.path.join(FP_PATH, "4.csv")
# X = pd.read_csv(EXAMPLE_FP_PATH, index_col = "index")

# SMILES_PATH  = os.path.join(DATA_PATH, "target.csv")
# smiles_df = pd.read_csv(SMILES_PATH, index_col = "index")
# y = smiles_df.loc[X.index,"t_cell"]

# #####################
# #Feature importance
# #####################

# imp_feat = Compute_Feature_Enrichment_DF(X,y)

# imp_feat = imp_feat.iloc[:8,:]

# print(imp_feat)

# #get the figures for the important features
# get_figures_of_FPs(imp_feat.index, smiles_df.loc[X.index,'smiles'], 'fp_tests_ch2')

# exit()

# #####################
# #Hist plot
# #####################

# # exit()
# # print(X.loc[:,'2245384272'])

# p = X.loc[y == 1,'3994088662']
# n = X.loc[y == 0,'3994088662']

# fig, ax = plt.subplots(1,1,figsize = (5,5))

# bins = np.arange(max(set(list(p) + list(n))) + 2) - 0.5

# ax.hist(p, alpha=0.2,label = "P (Epitopes)", weights=np.ones(len(p)) / len(p), bins  = bins)
# ax.hist(n, alpha=0.2,label = "N (ChEBI)", weights=np.ones(len(n)) / len(n), bins  = bins)

# import matplotlib.ticker as ticker

# ax.yaxis.set_major_formatter(ticker.PercentFormatter(1))

# plt.legend()
# plt.show()



