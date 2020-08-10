"""
NAM: Utilities.py (module)
DES: contains utilities & methods for fitting functions
"""
import numpy as np
from sklearn.decomposition import PCA

# --------- Utils for PCA ---------------------

def work_pca(data, n_components=2):
    data = data
    ndata = len(data)
    features = np.array(data).transpose()
    pca = PCA(n_components=n_components)
    pca.fit(features)
    pca_vectors = []
    results = []
    for i in range(n_components):
        pca_vectors.append(make_PCA_vector(data, pca.components_[i]))
    for i in range(ndata):  # loop over KIDs
        theData = data[i]
        res = []
#         res = [['#']]
        for j in range(len(data[0])):  # loop over TOD
            v = theData[j]
            res_raw = []
            pc_in = 0
            for k in range(n_components):  # loop over PCs
                dpc = pca_vectors[k][j] * pca.components_[k][i]
                # print i, j, k, dpc
                pc_in += dpc
            res.append([v - pc_in, pc_in, v])
        # ddsv.dump(res, resultFileKeyN.format(idMKIDs_use[i], n_components))
        results.append(res)

    return results

def make_PCA_vector(data, pca_comp):
    ndata = len(data)
    res = []
    for i in range(len(data[0])):
        v_pca = 0
        for j in range(ndata):
            v_pca += data[j][i] * pca_comp[j]
        res.append(v_pca)
    return res