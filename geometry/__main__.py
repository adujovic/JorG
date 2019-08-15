#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path,argv
path.insert(0,r'../')
from geometry.voronoi import Voronoi
from POSCARloader import POSCARloader
from aux.PeriodicTable import standardMass
import numpy as np

def main():
    pass

if __name__ == '__main__':
    loader   = POSCARloader(*argv[1:])
    loader.parse()
    volume = np.linalg.det(loader()['directions'])
    atoms  = [ atom[0] for atom in loader()['cell']]
    masses = [ standardMass[atom] for atom in atoms ]
    mass   = np.sum(masses)
    RHO    = mass/volume
    for atom,m in zip(loader()['cell'],masses):
        print(atom[0], "@", *atom[1], m, np.cbrt(3*m*volume/(mass*4*np.pi)))
    print(atoms,volume,mass,RHO)
    for i in range(len(loader)):
        voronoi = Voronoi(data=loader(i))
        voronoi.get_voronoi_diagram(save=True,name="%s_WS.txt"%argv[1+i])
        voronoi.show("%s.png"%argv[1+i])
        voronoi.show()
