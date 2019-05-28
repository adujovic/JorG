# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import numpy as np
from aux.periodic import maskFull,periodicTableElement
import spglib
from itertools import product

class Identity(dict):
    def __missing__(self, key):
            return key

class findFlips:
    wyckoffDict  = Identity() 
    logAccuracy  = 2
    crystal      = None 
    symmetry     = None
    atomTypeMask = maskFull
    Wyckoffs     = "abcdefghijklmnopqrstuvwz"
    distances    = set([])
    checker      = []

    def set_structure(self,cellSymmetry):
        self.symmetry = spglib.get_symmetry_dataset(cellSymmetry)
        self.crystal  = [(periodicTableElement[element-1],np.dot(cellSymmetry[0],atom))
                          for atom,element in zip(cellSymmetry[1],cellSymmetry[2])]

    def find_all_distances(self,referenceAtom,cutOff):
        self.distances = []
        for i,atom in enumerate(self.crystal):
           d = np.around(np.linalg.norm(atom[1]-referenceAtom[1]), decimals=self.logAccuracy)
           if ( self.wyckoffDict[self.symmetry['wyckoffs'][i]] in self.Wyckoffs and
                d > 0.0 and d <= cutOff and '$'+atom[0]+'$' in self.atomTypeMask):
                self.distances.append((d,atom[0]))
        
        self.distances = np.array(self.distances, dtype=[('distance', np.float), ('element', 'U3')])
        self.distances.sort(order='distance')

    def search_for_equivalent(self,distance,index,referenceAtom):
        for i in np.argwhere(self.symmetry['equivalent_atoms'] 
                          == self.symmetry['equivalent_atoms'][index]).flatten():
            if (distance == np.around(np.linalg.norm(self.crystal[i][1]-referenceAtom[1]),
                                     decimals=self.logAccuracy)
            and self.wyckoffDict[self.symmetry['wyckoffs'][i]] in self.Wyckoffs):
                self.checker[i] = True


    def unique(self,referenceAtom,cutOff):
        if self.symmetry is None:
            return None

        self.find_all_distances(referenceAtom,cutOff)
        
        self.checker = [ False for x in range(len(self.crystal))] 
        flipper = [] 
        case    = 1
        for ((distance,name),(i,atom)) in product(self.distances,enumerate(self.crystal)):
            if self.checker[i] or atom[0] != name:
                continue
            if distance == np.around(np.linalg.norm(atom[1]-referenceAtom[1]),
                                     decimals=self.logAccuracy):
                self.search_for_equivalent(distance,i,referenceAtom)
                flipper.append((i,atom,distance,self.wyckoffDict[self.symmetry['wyckoffs'][i]]))
                case += 1
        flipper = np.array(flipper, dtype=[('id', int), ('atom', list),
                                           ('distance', np.float), ('wyckoff', 'U1')])
        flipper.sort(order='distance')
        return flipper

    def all(self,referenceAtom,cutOff):
        flipper = []
        for i,atom in enumerate(self.crystal): 
            if "$"+atom[0]+"$" in self.atomTypeMask                         \
                and self.wyckoffDict[self.symmetry['wyckoffs'][i]] in self.Wyckoffs   \
                and np.linalg.norm(atom[1]-referenceAtom[1]) > 1e-3:
                flipper.append(i)

        return flipper


def find_all_distances(reference, 
                       crystal8,
                       cutOff, flipper,
                       logAccuracy = 2):
    acceptableFlips = np.append(flipper, reference)

    size = int(len(crystal8)/8)
    distances = []
    for i,flipA,j,flipB in product(range(8),
                                   acceptableFlips,
                                   range(8),
                                   acceptableFlips):
        d = np.round(np.linalg.norm(crystal8[flipA+size*i][1]
                                  - crystal8[flipB+size*j][1]),logAccuracy)
        if d not in distances and d <= cutOff:
            distances.append(d)
    distances = np.array(distances)
    distances.sort()
    return distances
            
