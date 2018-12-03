import numpy as np
def get_number_of_pictures(directions,cutOff,referenceAtom=[0,np.zeros(3)]):
    """
        Finding the amount of copies
        of the original SuperCell
                                        """

    multipliers = [] # returned 
    dDirs = np.tile(directions,(2,1)) 

    for i in range(3):
        normal = np.cross(dDirs[i+1],dDirs[i+2])
        normal /= np.linalg.norm(normal)
        height = np.dot(dDirs[i],normal)
        relative = np.dot(dDirs[i]-referenceAtom[1],normal)
        multipliers.append(1 + int((cutOff-relative)/height)) # calculating multipliers
    
    return multipliers

#
#
#
#
#

import numpy as np
from aux.periodic import periodicTableNumber
from aux.periodic import elementMagneticMoment
def generate_crystal(multiplyers,cell,directions,atomNames,reference,moments=None):
    """
        Generator of minimal required SuperCell
                                                """
    crystal = []
    for name in atomNames:
        for x in range(multiplyers[0]+1):
            for y in range(multiplyers[1]+1):
                for z in range(multiplyers[2]+1):
                    if moments is None:
                        for atom in cell:
                            if atomNames[atom[0]] == name:
                                position = np.copy(atom[1])
                                for a,n in zip([x,y,z],directions):
                                    position += a*n
                                    
                                flag = "%s"%atomNames[int(atom[0])]
                                crystal.append([flag,position,elementMagneticMoment[periodicTableElement[atomNames[int(atom[0])]]]])    
                    else:
                        for atom,moment in zip(cell,moments):
                            if atomNames[atom[0]] == name:
                                position = np.copy(atom[1])
                                for a,n in zip([x,y,z],directions):
                                    position += a*n
                                    
                                flag = "%s"%atomNames[int(atom[0])]
                                crystal.append([flag,position,moment])

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
    except: 
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
    except:
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
from aux.periodic import periodicTableNumber
from aux.periodic import periodicTableElement
from aux.periodic import maskFull 
def check_in_cell(cell,referenceAtom,directions,nearestNeighbor,atomNames,
                  originalCell,originalDirections,Wyckoffs='abcdefghijklmnopqrstuvwxyz',
                  originalSymmetry=None,atomTypeMask=maskFull):

    (wyckoffPositionDict,
     symmetry,
     originalSymmetry) = wyckoffs_dict(originalCell,
                                       cell,
                                       originalDirections,
                                       directions,
                                       atomNames)
    distances = []

    for i,(atom,wyck) in enumerate(
                         zip(cell,
                             symmetry['wyckoffs'])):
        wyckoff = wyckoffPositionDict[wyck]
        atomType = atomNames[atom[0]]
        distance = np.around(np.linalg.norm(atom[1] - referenceAtom[1]),2)
        if ( distance     not in distances    and
             wyckoff          in Wyckoffs     and
             "$"+atomType+"$" in atomTypeMask    ):
            distances.append(distance)

    distances.sort()
    if len(distances) <= nearestNeighbor:
        return False,distances,None
    return True,distances,wyckoffPositionDict

from itertools import permutations as per
from aux.periodic import elementMagneticMoment
def generate_from_NN(cell,referenceAtom,directions,nearestNeighbor,atomNames,
                   Wyckoffs='abcdefghijklmnopqrstuvwxyz',
                    atomTypeMask=maskFull, moments=None):
    originalSymmetryCell = (directions,
                            [np.dot(row[1],np.linalg.inv(directions)) for row in cell],
                            [periodicTableNumber[atomNames[row[0]]] for row in cell])
    originalSymmetry = spglib.get_symmetry_dataset(originalSymmetryCell)

    newReference = None

    # Looking for the closest-to-carthesian representation
    bestScore = 0.0             # starting score
    bestPermutation = (0,1,2)   # starting permutation
    for p in per(range(3)):
        score = 1.0
        for i,d in enumerate(directions):
            score *= d[p[i]]
        if np.abs(score) > bestScore: # best score wins
            bestScore = np.abs(score)
            bestPermutation = p

    for i,p in enumerate(bestPermutation):
        almostCarthesian     = np.zeros(3) 
        almostCarthesian[i]  = directions[i][p]
    
    found = []
    index = [1,1,1]
    for iterator in range(nearestNeighbor**3):
        expandedDirections = np.array([m*np.linalg.norm(d) for m,d in zip(index,almostCarthesian)])
        switch = (expandedDirections == np.min(expandedDirections))
        for i,s in enumerate(switch):
            if s:
                index[i] += 1

        crystal = []
        extraDirections = [m*d for m,d in zip(index,directions)]
        for name in atomNames:
            for x in range(index[0]):
                for y in range(index[1]):
                    for z in range(index[2]):
                        if moments is None:
                            for atom in cell:
                                if atomNames[atom[0]] == name:
                                    position = np.copy(atom[1])
                                    for a,n in zip([x,y,z],directions):
                                        position += a*n
                                    if np.linalg.norm(position - referenceAtom[1]) < 1e-3:
                                        newReference = len(crystal)
                                    crystal.append([atom[0],position,elementMagneticMoment[periodicTableElement[atomNames[atom[0]]]]])    
                        else:
                            for atom,moment in zip(cell,moments):
                                if atomNames[atom[0]] == name:
                                    position = np.copy(atom[1])
                                    for a,n in zip([x,y,z],directions):
                                        position += a*n
                                    if np.linalg.norm(position - referenceAtom[1]) < 1e-3:
                                        newReference = len(crystal)
                                    crystal.append([atom[0],position,moment])    
            ISFOUND,distances,wyckoffPositionDict = check_in_cell(
                                                     crystal,
                                                     referenceAtom,
                                                     extraDirections,
                                                     nearestNeighbor,
                                                     atomNames,
                                                     cell,
                                                     directions,
                                                     Wyckoffs,
                                                     originalSymmetry,
                                                     atomTypeMask)

        if(ISFOUND):
            symmetryCell = (extraDirections,
                            [np.dot(row[1],np.linalg.inv(extraDirections))
                                for row in crystal],
                            [periodicTableNumber[atomNames[row[0]]]
                                for row in crystal])
            crystalSymmetry = spglib.get_symmetry_dataset(symmetryCell)
            cutOff = 0.01+distances[nearestNeighbor]
            for atom in crystal:
                atom[0] = atomNames[atom[0]]
            return (cutOff,
                    crystal,
                    crystalSymmetry,
                    newReference,
                    index,
                    wyckoffPositionDict)

    return None

from itertools import product
def apply_mirrorsXYZ(dimensions,cell):
    outputCell = []
    print(cell)
    for p in product([0,-1],repeat=3):
        projection = np.array([p])
        translation = np.dot(projection,dimensions)[0]
        for atom in cell: 
            outputCell.append([atom[0],atom[1]+translation])
    return outputCell
