import numpy as np
import pandas as pd
import scipy.stats
import collections
import operator
import sklearn.cluster
import sklearn.metrics
from scipy.cluster.vq import kmeans2
from os import listdir
from os.path import isfile, join
from functools import reduce
from random import sample


######## Utility functions for power iteration clustering ########

def delta_val(v,v2):
    return np.sum(np.fabs(v2-v))


def normalize(v):
    max=v.max()
    min=v.min()
    return (v-min)/(max-min)


def initVector(m):
    n = m.shape[0]
    ovec = np.transpose(np.ones(n))
    v = m @ ovec
    sinv = 1.0 / np.sum(v)
    return v * sinv


def pic(a, maxiter = 50000, eps = 1e-7):
    m = a
    d1 = np.array(np.matrix(np.diag(a.sum(0))).I)
    w = d1 @ m
    n = w.shape[0]
    v = initVector(m)
    for i in range(maxiter):
        v2 = w @ v
        ninv = 1.0/ np.linalg.norm(v2, 1)
        v2 = v2 * ninv
        delta = delta_val(v, v2)
        v = v2
        if (delta * n) < eps:
            break
    return normalize(v)

def get_spectral_coor(X, k=2):
    X = np.abs(X)
    N = X.shape[0]
    np.fill_diagonal(X, np.zeros((N)))
    d = X.sum(0)
    I = np.diag(np.ones((N)))
    L = I - np.diag(1/np.sqrt(d)) @ X @ np.diag(1/np.sqrt(d))
    eigval_, eigvec_ = scipy.linalg.eigh(L, eigvals = (1, k))
    return eigvec_

######## Gene-Fishing Main Function ########
# Note: meta_recorder must be defined in the global environment
# Seeking improvements for this...


def gene_fishing(bait_genes, pool_genes, alpha = 5, k = 2,
                 cluster = "Spectral", affinity = "Spearman"):
    num_of_subsample = round(len(pool_genes) / (len(bait_genes) * alpha))
    shuffled_idx = np.array_split(sample(pool_genes, len(pool_genes)), num_of_subsample)
    for i in range(num_of_subsample):
        bait_expr_mat = expr_mat[bait_genes]
        pool_expr_mat = expr_mat[shuffled_idx[i]]
        bait_pool_gene_expr_mat = bait_expr_mat.join(pool_expr_mat)
        if affinity == "Spearman":
            bait_pool_gene_corr_mat, _ = scipy.stats.spearmanr(bait_pool_gene_expr_mat)
        elif affinity == "Cosine":
            bait_pool_gene_corr_mat = \
            sklearn.metrics.pairwise.cosine_similarity(np.transpose(bait_pool_gene_expr_mat))
        else:
            print("Incorrect affinity matrix specification")
            return
        if cluster == "Spectral":
            cluster_result = np.array([kmeans2(get_spectral_coor(bait_pool_gene_corr_mat), 2)[1]])
        elif cluster == "K-means":
            cluster_result = np.array([kmeans2(bait_pool_gene_corr_mat, 2)[1]])
        elif cluster == "PIC":
            cluster_result = np.array([kmeans2(pic(bait_pool_gene_corr_mat), 2)[1]])
        else:
            print("Incorrect cluster algorithm specification")
            return
        cluster_df = pd.DataFrame((cluster_result), columns = bait_genes + list(shuffled_idx[i]))
        bait_majority = collections.Counter(cluster_df[bait_genes].values[0]).most_common(1)[0][0]
        check_pool = (cluster_df == bait_majority)[shuffled_idx[i]]
        fished_idx = list(check_pool.columns[check_pool.all()])

        if len(fished_idx) == 0:
            pass
        else:
            meta_recorder[fished_idx] = meta_recorder[fished_idx] + 1