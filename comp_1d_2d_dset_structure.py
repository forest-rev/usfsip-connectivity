# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 16:25:46 2022

@author: pj276
"""

import networkit as nk
import time
import numpy as np
import sys, os
import h5py

# Get chunk ranges
# From https://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

# Create weighted graph
mg = nk.generators.MocnikGenerator(2, 50000, 2.6)
mgG = mg.generate()

nk.setNumberOfThreads(12)
print(f"There are {nk.getMaxNumberOfThreads()} threads available")

cLen = 500
os.chdir('C:/data/')

#%%
# 2D DATASET
f = h5py.File('cost_distances_2d7.h5', mode='a')
dset = f.create_dataset("dmat", (mgG.numberOfNodes(), mgG.numberOfNodes()),
                        compression="gzip", shuffle=True,
                        chunks=True, dtype='uint16', compression_opts=9)
#dset = f.create_dataset("dmat", (len(G8.nodes), len(G8.nodes)),
#                        compression="lzf", shuffle=True,
#                        chunks=True, dtype='f')
f.close()

# Timer
tic0 = time.perf_counter()
for i in chunker(list(range(mgG.numberOfNodes())),cLen):
    tic = time.perf_counter()
    spsp = nk.distance.SPSP(mgG, i)
    spsp.run()
    dists = spsp.getDistances()
    darr = np.array(dists)
    darr = np.tril(darr, i[0]-1) # Get array with elements above the k-th diagonal zeroed
    darr[darr==sys.float_info.max] = 0
    darr[darr > 500000] = 0
    #darr = np.rint(darr/10).astype('uint16')
    #print(darr)
    with h5py.File('cost_distances_2d7.h5', 'a', ) as f:
        f['dmat'][np.min(i):(np.max(i)+1),:] = darr
    toc = time.perf_counter()
    print(f"The subprocess took {toc - tic:0.4f} seconds")

toc0 = time.perf_counter()
print(f"The process took {toc0 - tic0:0.4f} seconds")

#%%
# 1D DATASET
f = h5py.File('cost_distances_1d7.h5', mode='a')
# 1D length (without zeros) is (nnodes^2)-nnodes
l1d = (mgG.numberOfNodes()**2)-mgG.numberOfNodes()
dset = f.create_dataset("dmat", (l1d,),
                        compression="gzip", shuffle=True,
                        chunks=True, dtype='uint16', compression_opts=9)
#dset = f.create_dataset("dmat", (len(G8.nodes), len(G8.nodes)),
#                        compression="lzf", shuffle=True,
#                        chunks=True, dtype='f')
f.close()

# Create vars to track start and end indices
sindex0 = [0]
eindex0 = [0]

# Timer
tic0 = time.perf_counter()
for i in chunker(list(range(mgG.numberOfNodes())),cLen):      
    tic = time.perf_counter()
    spsp = nk.distance.SPSP(mgG, i)
    spsp.run()
    dists = spsp.getDistances()
    for p,q in enumerate(i):
        dists[p] = dists[p][0:q]
    darr = np.array([item for sublist in dists for item in sublist])
    darr[darr==sys.float_info.max] = 0 # convert unreachable node dists to 0
    darr[darr > 500000] = 0 # convert nodes beyond max dist to 0    
    #darr = np.rint(darr/10).astype('uint16')
    #print(darr)

    # start and end indices
    if np.min(i) == 0:
        sindex = 0
        eindex = sindex + len(darr)
        sindex0.append(sindex)
        eindex0.append(eindex)
    else:
        sindex = eindex0[-1] # get last 
        eindex = sindex + len(darr)
        sindex0.append(sindex)
        eindex0.append(eindex)
    
    with h5py.File('cost_distances_1d7.h5', 'a', ) as f:
        f['dmat'][sindex:eindex,] = darr
    toc = time.perf_counter()
    print(f"The subprocess took {toc - tic:0.4f} seconds")

toc0 = time.perf_counter()
print(f"The process took {toc0 - tic0:0.4f} seconds")

#%%
# Check 2d dataset
f = h5py.File('cost_distances_2d.h5', 'r')
list(f.keys())
dset = f['dmat']
dset.shape # dimensions
dset[10,10] # diagonals are 0
dset[100,100] # diagonals are 0
f.close()

# Check 1d dataset
f = h5py.File('cost_distances_1d.h5', 'r')
list(f.keys())
dset = f['dmat']
dset.shape # dimensions

f.close()

#%%
# TESTING
X = np.array([[1,2,3],[4,5,6],[7,8,9]])
#array([[1, 2, 3],
#       [4, 5, 6],
#       [7, 8, 9]])

#get the upper triangular part of this matrix
v = X[np.tril_indices(X.shape[0], k = -1)]
print(v)
# [1 2 3 5 6 9]

# put it back into a 2D symmetric array
size_X = 3
X = np.zeros((size_X,size_X))
X[np.tril_indices(X.shape[0], k = 0)] = v
X = X + X.T - np.diag(np.diag(X))
#array([[1., 2., 3.],
#       [2., 5., 6.],
#       [3., 6., 9.]])

# Node indices
ni = np.array([0,1,2,3,4,5])
ni_1 = np.delete(ni-1,0)
sindex = np.cumsum(ni_1)
eindex = sindex + ni_1 


dists[0][0:500]
dists[1][0:501]

len(dists[0][0:500])
len(dists[1][0:501])

np.max(np.cumsum(range(500,1000)))
np.max(np.cumsum(range(0,5000)))



