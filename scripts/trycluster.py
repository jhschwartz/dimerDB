#from sklearn.neighbors import kneighbors_graph
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import itertools

L = 25
X = np.random.rand(L,L)
X = (X + X.T)/2 # make sym
np.fill_diagonal(X, 0) # set diag 
clust = AgglomerativeClustering(distance_threshold=0.5, affinity='precomputed', linkage='complete', n_clusters=None)
result = clust.fit(X)
print(result.labels_)


for cluster in range(5):
    print(f'in cluster {cluster}...')
    cluster_indecies = np.where(result.labels_ == cluster)
    for x,y in itertools.combinations(cluster_indecies[0], 2):
        print(X[x,y], x, y)

# print(X)
