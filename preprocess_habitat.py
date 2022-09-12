# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:42:51 2022

@author: pj276
"""
import rasterio as rio
from pathlib import Path
from pylab import imshow

apath = Path("G:/My Drive/Projects/USFSIP_Local/data/borneo_pone/Resistance_2000.tif")

with rio.open(apath) as ds:
    A = ds.read(1)


A[A<=ds.nodatavals[0]] = 1000
imshow(A)

from numpy.random import default_rng
rng = default_rng()

s = rng.uniform(0,1.001,1000)


