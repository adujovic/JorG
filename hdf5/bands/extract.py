#!/usr/bin/python3

import h5py
import numpy as np
from sys import argv

if len(argv) < 2:
    print("Feed me >>band<<.hdf5")
    exit()


File = h5py.File(argv[1],'r')
#['distance', 'eigenvector', 'frequency', 'label', 'path']
#<HDF5 dataset "distance": shape (8, 101), type "<f8">
#<HDF5 dataset "eigenvector": shape (8, 101, 48, 48), type "<c16">
#<HDF5 dataset "frequency": shape (8, 101, 48), type "<f8">
#<HDF5 dataset "label": shape (9,), type "|S10">
#<HDF5 dataset "path": shape (8, 101, 3), type "<f8">

for path,bands in zip(File['path'],File['frequency']):
    for momentum,frequencies in zip(path,bands):
#        for k in momentum:
#            print(k,end=" ")
#        print("",end="\t")
        for omega in frequencies:
            print(omega, end=" ")
        print("")

