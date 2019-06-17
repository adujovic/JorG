# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import numpy as np
def get_number_of_pictures(directions,cutOff,referenceAtom=[0,np.zeros(3)]):
#    """
#        Finding the amount of copies
#        of the original SuperCell
#                                        """
    multipliers = [] # returned
    dDirs = np.tile(directions,(2,1))
    #
    #          .---------------------------.
    #         / \                           \
    #        /   \                           \
    #       /     \                           \
    #      /       ^---------------------------.
    #     /       /                           /
    #    /       /     reference point       /
    #   /       / C   .                     /
    #  ^       /                           /
    # B \     /                           /
    #    \   /                           /
    #     \ /                           /
    #      0--------------------------->
    #                   A
    # A,B,C -> Vectors
    # normalA = (BxC)/|BxC|
    # normalB = (CxA)/|CxA|
    # normalC = (AxB)/|AxB|

    for i in range(3):
        normal = np.cross(dDirs[i+1],dDirs[i+2])
        normal /= np.linalg.norm(normal)
        height = np.dot(dDirs[i],normal)
        relative = np.dot(referenceAtom[1],normal)
        if cutOff > relative:
            multipliers.append(int((cutOff-relative)/height)) # calculating multipliers
        else:
            multipliers.append(0)

    return np.array(multipliers)

#
#
#
#
#

from JorG.PeriodicTable import elementMagneticMoment
from itertools import product,chain

class CrystalGenerator:
#    multipliers
    def __init__(self, cell, directions,
                 atomNames, reference=0):
        self.cell       = cell
        self.directions = directions
        self.reference  = reference
        self.atomNames  = atomNames
        self.moments    = None

    def fix_moments(self):
        if len(self.moments) >= len(self.cell):
            self.moments = self.moments[:len(self.cell)]
        else:
            self.moments = [elementMagneticMoment[atom[0]] for atom in self.cell]

    def __call__(self,multipliers):
    #    """
    #        Generator of minimal required SuperCell
    #                                                """

        self.fix_moments()
        crystal = []
        for (name,mul,(atom,moment)) in product(self.atomNames,
                                              product(range(multipliers[0]+1),
                                                      range(multipliers[1]+1),
                                                      range(multipliers[2]+1)),
                                              zip(self.cell,self.moments)):
            if atom[0] == name:
                position = np.copy(atom[1])
                position += np.dot(mul,self.directions)
                crystal.append([atom[0],position,moment])

        newReference = None
        for i,atom in enumerate(crystal):
          if np.linalg.norm(atom[1] - self.cell[self.reference][1]) < 1e-2:
            newReference = i
            break

        return (crystal,newReference)

#
#
#
#
#

import spglib
def wyckoffs_dict(originalCell, neotericCell): 
#   """
#      originalCell= (originalDirections,
#                     originalCell, 
#                     originalAtoms)
#      neotericCell= (neotericDirections,
#                     neotericCell,
#                     neotericAtoms)
#                                                   """
    neotericSymmetry = spglib.get_symmetry_dataset(neotericCell)
    originalSymmetry = spglib.get_symmetry_dataset(originalCell)

    wyckoffPositionDict = {}
    for (i,neotericAtom),\
        (j,originalAtom) in product(enumerate(neotericCell[1]),
                                    enumerate(originalCell[1])):
        if np.linalg.norm(np.array(neotericAtom)\
                         -np.array(originalAtom)) < 1e-3:
            wyckoffPositionDict[\
                 neotericSymmetry['wyckoffs'][i]] =\
                         originalSymmetry['wyckoffs'][j]

    return wyckoffPositionDict,neotericSymmetry,originalSymmetry
#
#
#
#
#

from JorG.Masks         import maskFull
from JorG.PeriodicTable import periodicTableNumber
class NearestNeighborsGenerator:
    Wyckoffs='abcdefghijklmnopqrstuvwxyz'
    originalSymmetry    = None
    atomTypeMask        = maskFull
    moments             = None
    extraMultiplier     = np.zeros(3,dtype=int)
    wyckoffPositionDict = None
    distances           = []
    cutOff              = -1
    crystal             = None
    crystalSymmetry     = None
    newReference        = None
    multipliers         = None
    ISFOUND             = False

    def __init__(self,cell,referenceAtom,directions):
        self.cell            = cell
        self.referenceAtom   = referenceAtom
        self.directions      = directions
        self.nearestNeighbor = 1
        self.fix_names()
        originalSymmetryCell = NearestNeighborsGenerator.get_symmetry(cell,directions)
        originalSymmetry     = spglib.get_symmetry_dataset(originalSymmetryCell)

    def fix_names(self):
        # finding unique names with conserved order
        # Using: stackoverflow.com/questions/15637336/numpy-unique-with-order-preserved
        names  = np.array([atom[0] for atom in self.cell],dtype=str)
        _, idx = np.unique(names, return_index=True)
        self.atomNames = names[np.sort(idx)]

    def check_in_cell(self,cell,referenceAtom,directions):
        originalCell = self.get_symmetry(self.cell,self.directions)
        neotericCell = self.get_symmetry(     cell,     directions)
        (self.wyckoffPositionDict,
         symmetry,
         self.originalSymmetry) = wyckoffs_dict(originalCell,neotericCell)
        self.distances = []
        for i,(atom,wyck) in enumerate(zip(cell,symmetry['wyckoffs'])):
            distance = np.around(np.linalg.norm(atom[1]-referenceAtom[1]),2)
            if (distance                   not in self.distances  and
                self.wyckoffPositionDict[wyck] in self.Wyckoffs   and
                "$%s$"%atom[0]                 in self.atomTypeMask ):
                self.distances.append(distance)
        self.distances.sort()
        if np.ma.masked_greater(self.distances,self.cutOff).count() <= self.nearestNeighbor:
            return False
        return True

    def fix_moments(self):
        if self.moments is None:
            self.moments = []
            for atom in self.cell:
                self.moments.append(elementMagneticMoment[atom[0]])

    @staticmethod
    def get_symmetry(cell,directions):
        return (directions,
                [np.dot(row[1],np.linalg.inv(directions))
                    for row in cell],
                [periodicTableNumber[row[0]]
                    for row in cell])

    @staticmethod
    def get_extra_directions(multipliers,directions):
        return [(mul+1)*d for mul,d in
                     zip(multipliers,
                         directions)]

    def find_new_refernce(self):
        self.newReference = None
        self.newReferenceAtom = None
        for i,atom in enumerate(self.crystal):
            if np.linalg.norm(atom[1] - self.referenceAtom[1]) < 1e-3:
                self.newReference = i
                self.newReferenceAtom = atom

    @staticmethod
    def get_cutOffs(directions,nearestNeighbor):
        minDirection = 0.99*np.min([np.linalg.norm(d) for d in directions])
        return minDirection*np.sqrt(np.arange(1.0,np.power(nearestNeighbor+1,3),1.0))

    def __call__(self,nearestNeighbor):
        self.nearestNeighbor = nearestNeighbor
        self.fix_moments()
        self.crystal = []
        for self.cutOff in self.get_cutOffs(self.directions,self.nearestNeighbor):
            self.multipliers =\
                    get_number_of_pictures(self.directions,
                                           self.cutOff,
                                           self.referenceAtom) + self.extraMultiplier
      
            extraDirections  =\
                    NearestNeighborsGenerator.get_extra_directions(self.multipliers,self.directions)
            self.generate_crystal()
            self.find_new_refernce()
            self.ISFOUND = self.check_in_cell(self.crystal,
                           self.newReferenceAtom,extraDirections)
            if self.ISFOUND :
                self.crystalSymmetry = spglib.get_symmetry_dataset(
                        NearestNeighborsGenerator.get_symmetry(
                            self.crystal,extraDirections))
                self.cutOff = self.distances[self.nearestNeighbor]
                return (self.cutOff,       self.crystal,     self.crystalSymmetry,
                        self.newReference, self.multipliers, self.wyckoffPositionDict)
        return None

    def generate_crystal(self):
        self.crystal = []
        multipliers  = [range(m+1) for m in self.multipliers]
        for (name,mul,(atom,moment)) in product(
                                          self.atomNames,
                                          product(*multipliers),
                                          zip(self.cell,self.moments)):
            if atom[0] != name:
                continue
            self.crystal.append([atom[0],np.copy(atom[1]) + np.dot(mul,self.directions),moment])

#
#
#
#
#

def apply_mirrorsXYZ(dimensions,cell, cutOff=-1.0, reference=0):
    outputCell = []
    for p in product([0,-1],repeat=3):
        projection = np.array([p])
        translation = np.dot(projection,dimensions)[0]
        for i,atom in enumerate(cell):
#            if cutOff > 0:
#                if np.linalg.norm(atom[1]+translation - cell[reference][1]) > cutOff:
#                    continue
            outputCell.append([atom[0],atom[1]+translation,atom[2],i])
    return outputCell
