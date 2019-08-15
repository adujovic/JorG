#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../../')
import time
from geometry.voronoi import Voronoi
from POSCARloader import POSCARloader

def main():
    pass

if __name__ == '__main__':
    tracker  = -(time.time())

    loader   = POSCARloader('POSCARsmall','POSCAR','POSCARweird')
    loader.parse()
    for i in range(len(loader)):
        voronoi = Voronoi(data=loader(i))
        voronoi.get_voronoi_diagram(save=True,name="WSradia%03d.dat"%i)
        voronoi.show("output%03d.png"%i)
        voronoi.show()
    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
