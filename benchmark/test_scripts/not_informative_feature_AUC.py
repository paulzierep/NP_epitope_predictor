import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
import pandas as pd
import numpy as np
np.random.seed(0)

import benchmark_utils

from sklearn.ensemble import RandomForestClassifier

X, y = load_iris(return_X_y = True)

X = pd.DataFrame(X)
y = pd.Series(y)

noise =  pd.DataFrame(np.random.random_sample((X.shape[0], 3000)))
# noise = pd.DataFrame(np.ones((X.shape[0], 100)))
X_noise = pd.concat([X, noise], axis = 1)
# print(X_noise)
# exit()

y = y[(y != 2)]
X = X_noise.loc[y.index,:]

clf = RandomForestClassifier(n_estimators = 100, random_state = 0)
f_vs_auc = benchmark_utils.evalutate_clf(X, y, clf)

f_vs_auc.plot(y = 'mean')
plt.show()