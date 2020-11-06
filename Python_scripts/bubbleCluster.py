#######################
#     Bubble Cluster  #
#     @author zliu    #
#     11/4/2020       #
#######################

'''
    Bubble cluster take the relative path of absolute contact matrix as input
and ouput the umap reduction result, one has 2 dims for plot, the other has (cells-1)
dims for clustering.
    also create a ARIresult.txt with ARI for different dimentions.
'''

import os
import sys
import time
import numpy as np
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI
from multiprocessing import Pool

from multiprocessing import Pool
from scipy.sparse import csr_matrix
from scipy.stats import chi2_contingency
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI

import umap
import hdbscan

########config########
contactMatrixPath = sys.argv[1]
res=1000000
pad=1
rp=0.5
prct=20
ndim=20
ncpus=10 #I don't really know why but it seems that it will run out all of your cpus
######################

######Global Vars#####
status = 0
######################


mm9dim = [197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430]
hg19dim = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560]

#sometimes contactMatrix directory may have files other than contactMatrix
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f

def neighbor_ave_cpu(A, pad):
    if pad==0:
        return A
    ngene, _ = A.shape
    ll = pad * 2 + 1
    B, C, D, E = [np.zeros((ngene + ll, ngene + ll)) for i in range(4)]
    B[(pad + 1):(pad + ngene + 1), (pad + 1):(pad + ngene + 1)] = A[:]
    F = B.cumsum(axis = 0).cumsum(axis = 1)
    C[ll :, ll:] = F[:-ll, :-ll]
    D[ll:, :] = F[:-ll, :]
    E[:, ll:] = F[:, :-ll]
    return (np.around(F + C - D - E, decimals=8)[ll:, ll:] / float(ll * ll))

def random_walk_cpu(A, rp):
    ngene, _ = A.shape
    A = A - np.diag(np.diag(A))
    A = A + np.diag(np.sum(A, axis=0) == 0)
    P = np.divide(A, np.sum(A, axis = 0))
    Q = np.eye(ngene)
    I = np.eye(ngene)
    for i in range(30):
        Q_new = (1 - rp) * I + rp * np.dot(Q, P)
        delta = np.linalg.norm(Q - Q_new)
        Q = Q_new.copy()
        if delta < 1e-6:
            break
    return Q

def impute_cpu(args):
    global status
    cell, ngene, pad, rp = args
    D = np.loadtxt(cell)
    A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
    A = np.log2(A + A.T + 1)
    A = neighbor_ave_cpu(A, pad)
    if rp==-1:
        Q = A[:]
    else:
        Q = random_walk_cpu(A, rp)
    #output current status
    print("finish cell impute: " + str(status))
    status += 1
    return [cell, Q.reshape(ngene*ngene)]

def read_cell(filepath):
    '''
    input: filepath for contact matrix
    output: totalCellNumber, cellLabel in list, cell relative path
    '''
    totalCell = 0
    label = []
    cellRelativePath = []

    filenames = list(listdir_nohidden(filepath))
    totalCell = len(filenames)
    for filename in filenames:
        label.append(filename.split(sep='.')[1])
    
    for cellname in filenames:
        cellRelativePath.append(filepath + "/" + cellname)

    return totalCell,label,cellRelativePath

def umap_reduce_cluster(matrix,label):
    '''
    input: cell*n*n contact matrix
    output: reduce matrix to cell*2(save) and cell*cell, cluster result show by ARI
    '''
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(matrix)
    np.save("embedding_plot",embedding)

    reducer_cluster = umap.UMAP(n_components=50)
    embedding_cluster = reducer_cluster.fit_transform(matrix)
    
    predict = hdbscan.HDBSCAN().fit_predict(embedding_cluster)

    return [ARI(label, list(hdbscan.HDBSCAN().fit_predict(embedding_cluster[:, :ndim]))) for ndim in [2,5,10,20,50]]

def main():
    ari = []
    global res
    chrTotalLength = sum(hg19dim)
    ngene = int(chrTotalLength / res)+1

    start_time = time.time()

    totalCell,label,cellRelativePath = read_cell(contactMatrixPath)

    paras = [[cell, ngene, pad, rp] for cell in cellRelativePath]
    p = Pool(ncpus)
    result = p.map(impute_cpu, paras)
    p.close()

    index = {x[0]:j for j,x in enumerate(result)}

    Q_concat = np.array([result[index[x]][1] for x in cellRelativePath])

    if prct>-1:
        thres = np.percentile(Q_concat, 100 - prct, axis=1)
        Q_concat = (Q_concat > thres[:, None]) 

    ari = umap_reduce_cluster(Q_concat,label)
    
    print(ari)
    
    ARIresFile = open('ARIresult.txt', 'w')
    for i in ari:
        ARIresFile.write(str(i))
        ARIresFile.write('\n')
    ARIresFile.close()

    end_time = time.time()
    print('Load and impute all cells with', end_time - start_time, 'seconds')

if __name__ == "__main__":
    main()