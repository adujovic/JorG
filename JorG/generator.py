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

    return multipliers

#
#
#
#
#

import numpy as np
from aux.periodic import elementMagneticMoment,periodicTableElement
from itertools import product
def generate_crystal(multipliers,cell,directions,atomNames,reference,moments=None):
#    """
#        Generator of minimal required SuperCell
#                                                """
    try:
        if len(moments) != len(cell):
            moments = None
    except TypeError:
        pass

    if moments is None:
        magneticMoments = moments
    elif len(moments) == len(cell):
        magneticMoments = moments
    else:
        magneticMoments = [periodicTableElement[atom[0]] for atom in cell]
    crystal = []
    for (name,mul,(atom,moment)) in product(atomNames,
                                          product(range(multipliers[0]+1),
                                                  range(multipliers[1]+1),
                                                  range(multipliers[2]+1)),
                                          zip(cell,magneticMoments)):
        if atom[0] == name:
            position = np.copy(atom[1])
            position += np.dot(mul,directions)
            crystal.append([atom[0],position,moment])

    newReference = None
    for i,atom in enumerate(crystal):
      if np.linalg.norm(atom[1] - cell[reference][1]) < 1e-2:
        newReference = i

    return (crystal,newReference)

#
#
#
#
#

import numpy as np
import spglib
from aux.periodic import periodicTableNumber
def wyckoffs_dict(originalCell,      cell,
                  originalDirections,directions,atomNames):
    symmetryCell = None
    try:
        symmetryCell = (directions,
                      [np.dot(row[1],np.linalg.inv(directions))
                          for row in cell],
                      [periodicTableNumber[atomNames[row[0]]]
                          for row in cell])
    except TypeError:
        symmetryCell = (directions,
                      [np.dot(row[1],np.linalg.inv(directions))
                          for row in cell],
                      [periodicTableNumber[row[0]]
                          for row in cell])

    symmetry = spglib.get_symmetry_dataset(symmetryCell)

    originalSymCell = None
    try:
        originalSymCell = (originalDirections,
                         [np.dot(row[1],np.linalg.inv(originalDirections))
                             for row in originalCell],
                         [periodicTableNumber[atomNames[row[0]]]
                             for row in originalCell])
    except TypeError:
        originalSymCell = (originalDirections,
                         [np.dot(row[1],np.linalg.inv(originalDirections))
                             for row in originalCell],
                         [periodicTableNumber[row[0]]
                             for row in originalCell])
    originalSym = spglib.get_symmetry_dataset(originalSymCell)

    wyckoffPositionDict = {}
    for i,atom in enumerate(cell):
        for j,originalAtom in enumerate(originalCell):
            if np.linalg.norm(atom[1]-originalAtom[1]) < 1e-3:
                wyckoffPositionDict[symmetry['wyckoffs'][i]] = originalSym['wyckoffs'][j]

    return wyckoffPositionDict,symmetry,originalSym
#
#
#
#
#

import numpy as np
import spglib
from aux.periodic import elementMagneticMoment,maskFull
from itertools import product

class generate_from_NN:
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


    def __init__(self,cell,
                 referenceAtom,
                 directions,
                 nearestNeighbor,
                 atomNames):
        self.cell            = cell
        self.referenceAtom   = referenceAtom
        self.directions      = directions
        self.nearestNeighbor = nearestNeighbor
        self.atomNames       = atomNames

        for atom in cell:
            try:
                atom[0] = atomNames[atom[0]]
            except Exception:
                pass

        try:
            originalSymmetryCell = generate_from_NN.get_symmetry(cell,directions)
            originalSymmetry = spglib.get_symmetry_dataset(originalSymmetryCell)
        except Exception:
            print("Failed to generate symmetry!")
            pass

        diagonal = np.sqrt(np.sum([np.dot(d-referenceAtom[1],d-referenceAtom[1]) for d in directions]))

    def check_in_cell(self,cell,referenceAtom,directions):
        (self.wyckoffPositionDict,
         symmetry,
         self.originalSymmetry) = wyckoffs_dict(self.cell,
                                           cell,
                                           self.directions,
                                           directions,
                                           self.atomNames)
        self.distances = []


        for i,(atom,wyck) in enumerate(
                             zip(cell,
                                 symmetry['wyckoffs'])):
            wyckoff = self.wyckoffPositionDict[wyck]
            atomType = atom[0]
            distance = np.around(np.linalg.norm(atom[1] - referenceAtom[1]),2)
            if ( distance     not in self.distances    and
                 wyckoff          in self.Wyckoffs     and
                 "$"+atomType+"$" in self.atomTypeMask    ):
                self.distances.append(distance)

        self.distances.sort()
        if np.ma.masked_greater(self.distances,self.cutOff).count() <= self.nearestNeighbor:
            return False
        return True

    def fix_moments(self):
        if self.moments is None:
            self.moments = []
            for atom in self.cell:
                self.moments.append(elementMagneticMoment[periodicTableElement[atom[0]]])

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

    def __call__(self):
        self.fix_moments()
        self.ISFOUND= self.check_in_cell(self.cell, self.referenceAtom, self.directions)
        if self.ISFOUND:
            symmetryCell = generate_from_NN.get_symmetry(self.cell,self.directions)
            self.crystalSymmetry = spglib.get_symmetry_dataset(symmetryCell)
            self.cutOff = self.distances[self.nearestNeighbor]
            return (self.cutOff,
                    self.cell,
                    self.crystalSymmetry,
                    self.newReference,
                    self.multipliers,
                    self.wyckoffPositionDict)

        self.crystal = []
        if len(self.distances) > 1:
            minDirection = self.distances[-1]
        else:
            minDirection = 0.99*np.min([np.linalg.norm(d) for d in self.directions])
        for cutOff in minDirection*np.sqrt(np.arange(1.0,np.power(self.nearestNeighbor+1,3),1.0)):
            self.cutOff = cutOff
            self.multipliers = get_number_of_pictures(self.directions,
                                                 self.cutOff,
                                                 self.referenceAtom)
            self.multipliers += self.extraMultiplier
            extraDirections = generate_from_NN.get_extra_directions(self.multipliers,self.directions)

            self.crystal = []
            for (name,mul,(atom,moment)) in product(
                                              self.atomNames,
                                              product(
                                                range(self.multipliers[0]+1),
                                                range(self.multipliers[1]+1),
                                                range(self.multipliers[2]+1)),
                                              zip(self.cell,self.moments)):
                if atom[0] == name:
                    position  = np.copy(atom[1])
                    position += np.dot(mul,self.directions)
                    self.crystal.append([atom[0],position,moment])


            self.find_new_refernce()

            self.ISFOUND = self.check_in_cell(
                           self.crystal,
                           self.newReferenceAtom,
                           extraDirections)
            if self.ISFOUND :
                symmetryCell = generate_from_NN.get_symmetry(self.crystal,extraDirections)
                self.crystalSymmetry = spglib.get_symmetry_dataset(symmetryCell)
                self.cutOff = self.distances[self.nearestNeighbor]
                return (self.cutOff,
                        self.crystal,
                        self.crystalSymmetry,
                        self.newReference,
                        self.multipliers,
                        self.wyckoffPositionDict)

        return None

#
#
#
#
#

from itertools import product
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
