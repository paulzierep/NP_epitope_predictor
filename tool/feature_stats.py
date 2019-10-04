import os
import pandas as pd

from statsmodels.stats.multitest import multipletests
from scipy import stats

import numpy as np

#pd.reset_option('^display.', silent=True)
pd.set_option('display.precision',2)
pd.set_option('display.max_columns', None)
pd.set_option('display.float_format', lambda x: '%.5f' % x)

def Remove_Correlated_Features(X, threshold = 0.8):

    """
    Removes highly correlated features
    """
    
    # Create correlation matrix
    corr_matrix = X.corr().abs()

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

    def ttest_apply(column):
        """
        Apply ttest, comparing the mean of each fp to the mean of the
        total population
        """
        

        #t, p = stats.ttest_1samp(column, X[column.name].mean())
        #t, p = stats.ttest_ind(column, X_n[column.name], equal_var = False)

        #t, p = stats.ranksums(column, X_n[column.name])

        #t, p = stats.ks_2samp(column, X_n[column.name])
        t, p = stats.mannwhitneyu(column, X[column.name])
        return(p)

    #p value of the stats test (currently ttest, positive vs negative)
    enrichment_df["p"] = X_p.apply(ttest_apply, axis = 0)

    #corrected p
    _, enrichment_df["p corr."], _ , _ = multipletests(enrichment_df["p"],  method="bonferroni")

    #how many samples do have the FP
    enrichment_df["Samples (P)"] = X_p.astype(bool).sum(axis=0)
    enrichment_df["Samples (N)"] = X_n.astype(bool).sum(axis=0)

    #mean of the FP count
    enrichment_df["MP"] = X_p.mean()
    enrichment_df["MN"] = X_n.mean()
    enrichment_df["M diff."] = enrichment_df["MP"] - enrichment_df["MN"]

    enrichment_df["Fold enrich."] = enrichment_df["MP"] / enrichment_df["MN"]

    enrichment_df["Samples (P) %"] = (enrichment_df["Samples (P)"] / X_p.shape[0]).round(2)
    enrichment_df["Samples (N) %"] = (enrichment_df["Samples (N)"] / X_n.shape[0]).round(2)

    #enrichment_df["Samples coverage %"] = (enrichment_df["Samples (P) %"] + enrichment_df["Samples (N) %"]) / 2

    #enrichment_df["p"] = enrichment_df["p"].apply(lambda x: '{:0.10e}'.format(x))
    #enrichment_df["p corr"] = enrichment_df["p corr"].round(2)

    enrichment_df = enrichment_df.loc[(enrichment_df["p corr."] < 0.05),:] 
    #enrichment_df = enrichment_df.loc[(enrichment_df["M diff."] > 0),:] 
    #enrichment_df = enrichment_df.loc[(enrichment_df["Samples (P)"] > 0.01),:] 
    #enrichment_df = enrichment_df.loc[(enrichment_df["Samples coverage %"] > 0.20),:] 

    #enrichment_df = enrichment_df.sort_values(by = ["s. p"], ascending = False)
    enrichment_df = enrichment_df.sort_values(by = ["p corr."], ascending = True)

    #enrichment_df = enrichment_df.sort_values(by = ["Fold enrich."], ascending = False)

    enrichment_df = enrichment_df.loc[:,["MP", "p corr.", "Samples (P) %", "Fold enrich.", "M diff."]]

    return(enrichment_df)

###############################
#Test
################################

# DATA_PATH = "ML_data"
# FP_PATH = os.path.join(DATA_PATH, "unfolded_fps")
# EXAMPLE_FP_PATH = os.path.join(FP_PATH, "0.csv")
# X = pd.read_csv(EXAMPLE_FP_PATH, index_col = "index")

# SMILES_PATH  = os.path.join(DATA_PATH, "target.csv")
# smiles_df = pd.read_csv(SMILES_PATH, index_col = "index")
# y = smiles_df.loc[X.index,"b_cell"]

# X = Compute_Feature_Enrichment_DF(X,y)
# print(X)
