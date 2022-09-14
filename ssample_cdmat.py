# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 15:00:41 2022

@author: pj276
"""
import numpy as np, math, time
from scipy import sparse
from scipy.sparse import diags
from pylab import imshow
import networkx as nx
import networkit as nk
import rasterio as rio
import h5py
import scipy.sparse as ss
import h5sparse
from osgeo import gdal
import pandas as pd
import random

#%%
nk.getMaxNumberOfThreads() # see maximum number of available threads

#%% 
# Sampling frequency, specify every Nth cell
cellSampling = 10

#%%
# Functions
# Minor modification from 
# https://stackoverflow.com/questions/30199070/how-to-create-a-4-or-8-connected-adjacency-matrix
# Diagonal weights set to sqrt(2). Orthogonal weights set to 1.
# Calculation for conductance is then conductance value/cellres*weight
def connected_adjacency(image, connect, patch_size=(1, 1)):
    """
    Creates an adjacency matrix from an image where nodes are considered adjacent 
    based on 4-connected or 8-connected pixel neighborhoods.
    :param image: 2 or 3 dim array
    :param connect: string, either '4' or '8'
    :param patch_size: tuple (n,m) used if the image will be decomposed into 
                   contiguous, non-overlapping patches of size n x m. The 
                   adjacency matrix will be formed from the smaller sized array
                   e.g. original image size = 256 x 256, patch_size=(8, 8), 
                   then the image under consideration is of size 32 x 32 and 
                   the adjacency matrix will be of size 
                   32**2 x 32**2 = 1024 x 1024
    :return: adjacency matrix as a sparse matrix (type=scipy.sparse.csr.csr_matrix)
    """
    r, c = image.shape[:2]
    r = int(r / patch_size[0])
    c = int(c / patch_size[1])
    if connect == '4':
        # constructed from 2 diagonals above the main diagonal
        d1 = np.tile(np.append(np.ones(c-1), [0]), r)[:-1]
        d2 = np.ones(c*(r-1))
        upper_diags = diags([d1, d2], [1, c])
        return upper_diags + upper_diags.T
    elif connect == '8':
        # constructed from 4 diagonals above the main diagonal
        d1 = np.tile(np.append(np.ones(c-1), [0]), r)[:-1]
        d2 = np.append([0], d1[:c*(r-1)])
        d3 = np.ones(c*(r-1))
        d4 = d2[1:-1]
        d4[d4==1] = 2.0**0.5
        upper_diags = diags([d1, d2, d3, d4], [1, c-1, c, c+1])
        return upper_diags + upper_diags.T
    else:
        raise ValueError('Invalid parameter \'connect\'={connect}, must be "4" or "8".'
                     .format(connect=repr(connect)))

# Get chunk ranges
# From https://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

#%%
# Read coordinates
# Set your working directory
dataDir = 'G:/My Drive/Projects/USFSIP_Connectivity/Data/GeoData/cleo_test/'
outDir = 'G:/My Drive/Projects/USFSIP_Connectivity/Data/GeoData/cleo_test/out'

# XYfile name
xyfilename = dataDir + '/xys/SP_base.xy'

# Xys
xys = pd.read_csv(xyfilename)

infile = dataDir + '/resistance/resist_base_roads.rsg'
# Use dataset index method to get cell numbers for each coordinate pair
with rio.open(infile) as ds:
    #cinds = [ds.index(i[1].X,i[1].Y) for i in xys.iterrows()]
    #sources = [(i[0]*ds.shape[1])+i[1] for i in cinds]
    A = ds.read(1)
    A[A==-9999]=1000
    tic = time.perf_counter()
    adj8 = connected_adjacency(ds, '8').astype("float32").tocsr()
    # Return the indices and values of the nonzero elements of a matrix
    # Output contains the row indices, column indices, and values of the nonzero matrix entries
    sinds = sparse.find(adj8)
    A = A.flatten()
    adj8[sinds[0],sinds[1]] = ((A[sinds[0]] + A[sinds[1]])/2)*sinds[2]*500
    G8 = nx.from_scipy_sparse_array(adj8)
    nx.set_node_attributes(G8, dict(zip(range(0,G8.number_of_nodes()), [{'value': i} for i in A])))
    
    node_values = nx.get_node_attributes(G8, 'value')
    G8.remove_nodes_from((n for n, w in node_values.items() if w == 1000))
    #edge_weights = nx.get_edge_attributes(G8,'weight')
    #G8.remove_edges_from((e for e, w in edge_weights.items() if w == 0))
    aa = nk.nxadapter.nx2nk(G8, weightAttr='weight')
    aa.indexEdges()
    toc = time.perf_counter()
    print(f"The process took {toc - tic:0.4f} seconds")

#%%
# Sample the resistance surface
with rio.open(infile) as ds:
    A = ds.read(1)
    # Create regular array, every Nth cell gets a value of 1
    x = np.zeros(A.shape,dtype=np.byte)
    x[0::cellSampling,::cellSampling] = 1
    # Get regularly sampled indices corresponding to the target patches
    cinds1 = np.argwhere((A!=-9999) & (x==1))
    cinds2 = np.argwhere((A!=-9999) & (x==1))
        
        
#%%
import sys, os
os.chdir('C:/data/')

#%%
# Write data in chunks to an array in the same dataset
# Takes about 45s per chunk (291 chunks total when running 1000 nodes at a time)
# Create a huge file because the array needs to be square

# Timer
tic0 = time.perf_counter()

# Number of nodes from which to calculate distances for each iteration
cLen = 500

# hdf file to hold distances
f = h5py.File('cost_distances2.h5', mode='a')
dset = f.create_dataset("dmat", (len(G8.nodes), len(G8.nodes)),
                        compression="gzip", shuffle=True,
                        chunks=True, dtype='uint16', compression_opts=9)
f.close()

xx = list(chunker(list(range(len(list(G8.nodes)))),cLen))[80:82]
for i in xx:
# i in this loop is simply a list of consecutive integers
#for i in chunker(list(range(len(list(G8.nodes)))),cLen):
    tic = time.perf_counter()
    spsp = nk.distance.SPSP(aa, i)
    spsp.run()
    dists = spsp.getDistances()
    darr = np.array(dists)
    darr = np.tril(darr, i[0]-1)
    darr[darr==sys.float_info.max] = 0
    darr[darr > 500000] = 0
    darr = np.rint(darr/10).astype('uint16')
    #print(darr)
    with h5py.File('cost_distances.h5', 'a', ) as f:
        f['dmat'][np.min(i):(np.max(i)+1),:] = darr
    toc = time.perf_counter()
    print(f"The subprocess took {toc - tic:0.4f} seconds")

toc0 = time.perf_counter()
print(f"The entire process took {toc0 - tic0:0.4f} seconds")

# Access
tic = time.perf_counter()
f = h5py.File('cost_distances.h5', mode='r')
f['dmat'][[0,5000,8000,10000,20000,30000,50000,100000,150000,200000],:]
#f['dmat'][[0,1],:]
f.close()
toc = time.perf_counter()
print(f"Accessing rows took {toc - tic:0.4f} seconds")

# Loop access
tic = time.perf_counter()
elist = []
with h5py.File('cost_distances.h5', 'r') as f:
    for p in [0,5000,8000,10000,20000,30000,50000,100000,150000,200000]:
        elist.append(f['dmat'][p,:])
toc = time.perf_counter()
print(f"Accessing rows took {toc - tic:0.4f} seconds")

# Loop access
tic = time.perf_counter()
elist = []
with h5py.File('cost_distances.h5', 'r') as f:
    for p in [(0,20),(5000,6000),(8000,200000),(10000,19000),
              (20000,40000),(30000,60000),(50000,120000),(100000,190000),
              (150000,200050),(200000,200060)]:
        elist.append(f['dmat'][p[1],p[0]])
toc = time.perf_counter()
print(f"Accessing rows took {toc - tic:0.4f} seconds")

# Verify upper triangle is zeros
tic = time.perf_counter()
with h5py.File('cost_distances.h5', 'r') as f:
    for p in [1000,91500]:
        print(np.max(f['dmat'][p,p:-1]))
        print(np.max(f['dmat'][p,:]))
toc = time.perf_counter()
print(f"Accessing rows took {toc - tic:0.4f} seconds")

