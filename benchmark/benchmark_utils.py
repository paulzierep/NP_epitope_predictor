import pandas as pd
import os
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import chi2
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score, cross_val_predict

import numpy as np

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

# def copy_predictor_for_classification():
#     pass

def get_X_y(predictor, cluster = 0, cell_type = 'b_cell'):
    '''
    Helper func. which gets the ML X, y data 
    from the epitope_predictor
    '''

    target_df = pd.read_csv(predictor.target_storage, index_col = 'index')

    for file in os.listdir(predictor.fp_unfolded_storage):
        if str(cluster) in file:

            #get the FPs
            fp_df = pd.read_csv(os.path.join(predictor.fp_unfolded_storage, file), index_col = 'index')
            fp_df.sort_index(inplace = True)
            fp_df.drop_duplicates(inplace = True) #remove doubs

            target_df = target_df.loc[fp_df.index,:] #do not take doublicates


    print('*'*50)
    print('Loading cluster {0} for cell type {1} with data size of {2}'.format(cluster, cell_type, target_df.shape[0]))
    print('Indices match: {0}'.format((fp_df.index == target_df.index).all()))

    X = fp_df
    y = target_df.loc[:,cell_type]

    print('Label distribution: \n', y.value_counts())

    return(X,y)

def evalutate_clf(X,y, clf):
    '''
    Evaluates a given clf, usinf feature correlation removal, and feature range.
    '''

    olf_f = X.shape[1]
    X = Remove_Correlated_Features(X)
    new_f = X.shape[1]

    print('Features before reduction: {0}; after reduction: {1}'.format(olf_f, new_f))

    cut_off = X.shape[1] #num of max features 
    f_range = np.around(np.geomspace(1, 2048, num=12)).astype(int) #log scale
    f_range = np.array(list(filter(lambda x: x < cut_off, f_range))) #log scale with cut-off

    f_vs_auc = pd.DataFrame()
    feature_vs_aucs = feature_vs_auc(X, y, f_range, clf = clf)

    mean = feature_vs_aucs.mean(axis = 0)
    yerr = feature_vs_aucs.sem(axis = 0)
    
    f_vs_auc["mean"] = mean
    f_vs_auc["err"] = yerr

    return(f_vs_auc)

def feature_vs_auc(X,y, f_range, clf = None, iterations = 3):
    '''
    Computes the AUC for a range of featues for a given clf
    '''

    feature_vs_aucs = pd.DataFrame()
    
    for k_best in f_range:
        
        print(k_best)
        
        f_selection = SelectKBest(chi2, k=k_best)
        multiple_Kfold_aucs = []
        
        for i in range(iterations):

            pipe_clf = Pipeline([
            #('sampling', RandomUnderSampler()),
            ("sel", f_selection),
            ('clf', clf),
            ])

            Kfold_aucs = cross_val_score(pipe_clf, 
                                         X, 
                                         y, 
                                         cv=StratifiedKFold(shuffle=True, n_splits = 5), 
                                         scoring='roc_auc', 
                                         n_jobs = -1)
            #print(Kfold_aucs)
            
            multiple_Kfold_aucs.extend(Kfold_aucs)
    
        feature_vs_aucs[k_best] = multiple_Kfold_aucs
        
    return(feature_vs_aucs)

def assign_random_y(y):
    '''
    Assign the same amount of epitopes but random, show that clf cannot learn arbitrary samples
    '''

    count = y[y == 1].shape[0]
    y.loc[:] = 0
    y.loc[y.sample(count).index] = 1

    return(y)


def evaluate_all_clusters(predictor, clf, storage_folder, mode = 'normal'):
    ''''
    Evaluates all clusters, and all cell types, stores the results
    '''
    os.makedirs(storage_folder, exist_ok = True)

    for cell_type in ['b_cell','t_cell']:

        OUTPUT_FOLDER = os.path.join(storage_folder, cell_type)
        os.makedirs(OUTPUT_FOLDER, exist_ok = True)

        for cluster in range(0,8):

            X,y = get_X_y(predictor,  cluster = cluster, cell_type = cell_type)

            if mode == 'random':
                y = assign_random_y(y)

            if y[y == 1].empty:
                print("Not enough targets for cluster {0}".format(cluster))

            else:
                f_vs_auc = evalutate_clf(X,y, clf)
                OUTPUT_PATH = os.path.join(OUTPUT_FOLDER, '{0}.csv'.format(cluster))
                f_vs_auc.to_csv(OUTPUT_PATH)