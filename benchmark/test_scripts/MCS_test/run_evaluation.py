import benchmark_utils
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier

target_df = pd.read_csv('target.csv')
X = pd.read_csv('X.csv', index_col = 0)
y = target_df.loc[X.index, 'b_cell']
# print(X)
# print(y)

clf = RandomForestClassifier(n_estimators = 100, random_state = 0)
f_vs_auc = benchmark_utils.evalutate_clf(X, y, clf)

f_vs_auc.plot(y = 'mean')
plt.show()