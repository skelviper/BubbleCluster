#######################
#     Bubble Cluster  #
#     @author zliu    #
#     11/8/2020       #
#######################

'''
    Bubble cluster take the relative path of absolute-pos contact matrix as input
and ouput the reduction result for plot.
    If you are using benchmark mode, it will create a result.txt with ARI indicate the cluster performance.

    usage,check:
        python bubbleCluster.py -h
'''

import os
import sys
import time
import argparse
import ast

import numpy as np
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI
from multiprocessing import Pool
from scipy.sparse import csr_matrix
from scipy.stats import chi2_contingency
from sklearn.decomposition import PCA

import umap
import hdbscan
from sklearn.metrics import silhouette_score

########config########
parser = argparse.ArgumentParser(description="BubbleCluster, a software for clustering single cell Hi-C data")
parser.add_argument(
                    "-i","--input",
                    dest="contactMatrixPath",
                    type=str,
                    action="store",
                    help="input folder path of contactMatrix"
)

parser.add_argument(
                    "-r","--resolution",
                    dest="res",
                    type=int,
                    action="store",
                    default=1000000,
                    help="resolution for bubbleCluster analysis,default to 1 megabase"
)
parser.add_argument(
                    "-p","--pad",
                    dest="pad",
                    type=int,
                    default=1,
                    action="store",
                    help="pad for convolution,default to 1, equal to 1 megabase when you set resolution it's default"
)
parser.add_argument(
                    "--rp",
                    dest="rp",
                    type=int,
                    default=0.5,
                    action="store",
                    help="random walk restart probability,default is 0.5"
)
parser.add_argument(
                    "--prct",
                    dest="prct",
                    type=str,
                    default=20,
                    action="store",
                    help="prct of the most significant contacts you use for PCA"
)
parser.add_argument(
                    "-d","--ndim",
                    dest="ndim",
                    default=20,
                    type=int,
                    action="store",
                    help="max dimentions"
)
parser.add_argument(
                    "-t","--threads",
                    dest="ncpus",
                    type=int,
                    default=48,
                    action="store",
                    help="total threads of your machine, notice that bubbleCluster draws all threads you have, manage it with smk",
)
parser.add_argument(
    '-b','--benchmark',
    help='Benchmark mode, default to false. This would read cell type name from your pairs file and calc ARI. Input should be either "True" or "False".',
    type=ast.literal_eval,
    default=False,
    dest='benchmark',
)
parser.add_argument(
    '-l','--length',
    help = 'read chromosome length from file',
    type=str,
    dest='chrLengthPath',
)


args = parser.parse_args()

######################

########config########
contactMatrixPath = args.contactMatrixPath
res=args.res
pad=args.pad
rp=args.rp
prct=args.prct
ndim=args.ndim
ncpus=args.ncpus
benchmark=args.benchmark
chrLengthPath=args.chrLengthPath
######################

######################

#sometimes contactMatrix directory may have files other than contactMatrix,for example .ipycheckpoints created by jupyter.
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f

def read_chrlength(chrLengthPath):
    chrdim = []
    with open(chrLengthPath, 'r') as chrLength:
        lines = chrLength.readlines()
        for line in lines:
            chrdim.append(int(line.split(sep='\t')[1]))
    return chrdim

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
    cell, ngene, pad, rp = args
    D = np.loadtxt(cell)
    A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
    A = np.log2(A + A.T + 1)
    A = neighbor_ave_cpu(A, pad)
    if rp==-1:
        Q = A[:]
    else:
        Q = random_walk_cpu(A, rp)
    
    return [cell, Q.reshape(ngene*ngene)]

def read_cell(filepath):
    '''
    input: filepath for contact matrix
    output: totalCellNumber, cellLabel in list, cell relative path
    '''
    totalCell = 0
    cellRelativePath = []

    filenames = sorted(list(listdir_nohidden(filepath)))
    totalCell = len(filenames)

    for cellname in filenames:
        cellRelativePath.append(filepath + "/" + cellname)

    return totalCell,cellRelativePath

def read_label(filepath):
    '''
    read label is only use for testing performance, in this case, your contactMatrix(.conmat) file 
    should have cell type in it.
    '''
    label = []

    filenames = sorted(list(listdir_nohidden(filepath)))

    for filename in filenames:
        #label.append("".join(filename.split(sep='_')[1:2]).split(sep=".")[0]) # for Tan's dataset
        label.append(filename.split(sep='.')[1]) # for Ramani's dataset

    return label


def pca_reduce(matrix):
    '''
    input: ncell*n*n contact matrix 
    output: pca reduced matrix to ncell*ncell, a list contain percentage of each PCs.

    note: Come on, at least you have 30 cells right?
    '''
    reducer = PCA()
    matrix_reduce = reducer.fit_transform(matrix)

    return matrix_reduce,reducer.explained_variance_ratio_

def dicide_optimised_pcs(pcaMatrix):
    '''
    after PCA reduction, it's important to use only main PCs for clustering so that
    clustering results are not affected by random noise.
    input: PCA reduced matrix with dim ncell*necll
    output: min and max index for UMAP. 
    '''
    max_dim = 0
    min_dim = 0

    pca_dim_start = 10
    mostSuitableMax = []
    for i in range(pca_dim_start,round(len(pcaMatrix)*0.05)):
        reducer_cluster = umap.UMAP(random_state=42)
        embedding_cluster = reducer_cluster.fit_transform(pcaMatrix[:,0:i])
        mostSuitableMax.append(silhouette_score(embedding_cluster,list(hdbscan.HDBSCAN().fit_predict(embedding_cluster)),metric='euclidean'))

    max_dim = mostSuitableMax.index(max(mostSuitableMax)) + pca_dim_start
    #print(mostSuitableMax)


    mostSuitableMin = []
    for j in range(5):
        reducer_cluster = umap.UMAP(random_state=42)
        embedding_cluster = reducer_cluster.fit_transform(pcaMatrix[:,j:max_dim])
        mostSuitableMin.append(silhouette_score(embedding_cluster,list(hdbscan.HDBSCAN().fit_predict(embedding_cluster)),metric='euclidean'))

    min_dim = mostSuitableMin.index(max(mostSuitableMin))
    #print(mostSuitableMin)

    return min_dim,max_dim


def main():
    chrLength = read_chrlength(chrLengthPath)
    chrTotalLength = sum(chrLength)
    ngene = int(chrTotalLength / res)+1

    start_time = time.time()

    totalCell,cellRelativePath = read_cell(contactMatrixPath)

    paras = [[cell, ngene, pad, rp] for cell in cellRelativePath]

    result = []

    for i in range(totalCell//ncpus+1):
        with Pool(ncpus) as p:
            result += p.map(impute_cpu,paras[i*ncpus:i*ncpus+ncpus])

    index = {x[0]:j for j,x in enumerate(result)}
    Q_concat = np.array([result[index[x]][1] for x in cellRelativePath])
    del(result) #free memory

    # needed to optimise for memory consumption
    if prct>-1:
        thres = np.percentile(Q_concat, 100 - prct, axis=1)
        Q_concat = (Q_concat > thres[:, None]) 

    pcaMatrix,varianceRatio = pca_reduce(Q_concat)
    # end to modify
    
    # save for protential later use
    np.save("pcaMatrix",pcaMatrix)
    np.savetxt("varianceRatio",varianceRatio)
    
    min_dim,max_dim = dicide_optimised_pcs(pcaMatrix)

    #because i want an reproducible result,so set random_state to 42
    reducer_umap = umap.UMAP(random_state=42) 
    embedding_cluster = reducer_umap.fit_transform(pcaMatrix[:,min_dim:max_dim])

    np.save("umapMatrix",embedding_cluster)
    
    #benchmark
    if (benchmark==True):
        label = read_label(contactMatrixPath)
        ari = ARI(label, list(hdbscan.HDBSCAN().fit_predict(embedding_cluster)))

        with open('result.txt', 'w') as resFile:
            resFile.write("Benchmark mode, ARI is " + str(ari) + "\n")
            resFile.write("PCs selected are {} to {} \n".format(str(min_dim),str(max_dim)))
            end_time = time.time()
            resFile.write('Load and impute all cells with', end_time - start_time, 'seconds')

    
    

if __name__ == "__main__":
    main()