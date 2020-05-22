import pandas as pd
import numpy as np
# print(pd.np.random.rand(5,30000))
# print(pd.np.corrcoef(pd.np.random.rand(5,30000)))

# df = pd.DataFrame(pd.np.random.rand(5,30000))
# print(df)
# print(df.corr())
# print(pd.np.corrcoef(pd.np.random.rand(5,30000)))

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

def Remove_Correlated_Features_np(X, threshold = 0.8):

    """
    Removes highly correlated features
    """
    
    # Create correlation matrix
    # https://stackoverflow.com/questions/48270953/pandas-corr-and-corrwith-very-slow
    # use np vectorization and avoid Nan lookup !

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

# np.random.seed(10)
df = pd.DataFrame(np.random.rand(20,10))
# df.columns = ['a','b','c']
print(df)

trim_df = Remove_Correlated_Features(df)
print(trim_df)
trim_df = Remove_Correlated_Features_np(df)
print(trim_df)