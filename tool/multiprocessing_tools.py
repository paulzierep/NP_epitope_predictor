from multiprocessing import  Pool
from functools import partial
import numpy as np
import pandas as pd

def df_split_apply_multip_combine(df, df_func, kwargs, num_of_processes = 2):
    '''
    splits a df into seqments, performes a func on each in parallel (num_of_processes), the 
    func can use the kwargs
    '''

    data_split = np.array_split(df, num_of_processes) #split df
    pool = Pool(num_of_processes) #start pool
    func = partial(df_func, kwargs = kwargs) #add kwargs to func
    data = pool.map(func, data_split) #apply func on each df-split
    pool.close()
    pool.join()

    data = pd.concat(data) #combine data

    return(data)

#############################
#Tests
#############################

# def df_func(df, kwargs):
#   print(kwargs)
#   df = df*2
#   return(df)

# df = pd.DataFrame([[0,1,2,3],[0,1,1,1]])
# kwargs = {'A':'1', 'B':'2'}

# data = df_split_apply_multip_combine(df, df_func, kwargs, num_of_processes = 10)

# print(data)