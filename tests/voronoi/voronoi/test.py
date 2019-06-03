#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,maxsize,path
path.insert(0,r'../../../')
import time
from voronoi.voronoi import Voronoi
import numpy as np
from JorG.loadsave import POSCARloader
from itertools import product

def main(**args):
    pass

if __name__ == '__main__':
    tracker  = -(time.time())

    loader   = POSCARloader('POSCARsmall','POSCAR','POSCARweird',spam=False)
    loader.parse()
    for i in range(len(loader)):
        voronoi = Voronoi(data=loader(i))
        voronoi.fill_points()
        voronoi.get_Voronoi_diagram()
        voronoi.show()
    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
