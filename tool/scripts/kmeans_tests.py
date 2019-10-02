from sklearn import datasets
from sklearn.cluster import KMeans
import pandas as pd

# import some data to play with
# X, y = datasets.make_classification()

# print(y)
# print(X)

X = pd.read_csv("../ML_data/fp_folded.csv", index_col = "index").iloc[:1000]

cluster_clf = KMeans(n_clusters=8, random_state=0, n_jobs = -1)
clusters = cluster_clf.fit_predict(X)

result = pd.DataFrame()
result["fit_pred"] = clusters
# result["fit"] = clusters2
#result["match"] = (result["fit_pred"] == result["fit"])

#print(result["match"].all())
print(result)