import numpy as np
def get_number_of_pictures(directions,cutOff,referenceAtom=[0,np.zeros(3)]):
    """
        Finding the amount of copies
        of the original SuperCell
                                        """
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
        relative = np.dot(dDirs[i]-referenceAtom[1],normal)
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
from aux.periodic import periodicTableNumber
from aux.periodic import elementMagneticMoment
def generate_crystal(multiplyers,cell,directions,atomNames,reference,moments=None):
    """
        Generator of minimal required SuperCell
                                                """
    try:
        if len(moments) != len(cell):
            moments = None
    except:
        pass

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
                  originalCell,originalDirections, cutOff,
                  Wyckoffs='abcdefghijklmnopqrstuvwxyz',
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
    if np.ma.masked_greater(distances,cutOff).count() <= nearestNeighbor:
        return False,distances,None
    return True,distances,wyckoffPositionDict

#
#
#
#
#
#
#
from aux.periodic import elementMagneticMoment
def generate_from_NN(cell,
        referenceAtom,directions,nearestNeighbor,atomNames,
        Wyckoffs='abcdefghijklmnopqrstuvwxyz',
        atomTypeMask=maskFull, moments=None, extraMultiplier=np.zeros(3,dtype=int)):

    try:
        if len(moments) != len(cell):
            moments = None
    except:
       moments = None 

    try:
        originalSymmetryCell = (directions,
                            [np.dot(row[1],np.linalg.inv(directions)) for row in cell],
                            [periodicTableNumber[atomNames[row[0]]] for row in cell])
        originalSymmetry = spglib.get_symmetry_dataset(originalSymmetryCell)
    except:
        pass

    diagonal = np.sqrt(np.sum([np.dot(d-referenceAtom[1],d-referenceAtom[1]) for d in directions]))
    ISFOUND,distances,wyckoffPositionDict = \
      check_in_cell(cell,
                    referenceAtom,
                    directions,
                    nearestNeighbor,
                    atomNames,
                    cell,
                    directions,
                    diagonal,
                    Wyckoffs,
                    originalSymmetry,
                    atomTypeMask)

    crystal = []
    if len(distances) > 1:
        minDirection = distances[-1]
    else:
        minDirection = 0.99*np.min([np.linalg.norm(d) for d in directions])
    for cutOff in minDirection*np.sqrt(np.arange(1.0,np.power(nearestNeighbor+1,3),1.0)):
        multipliers = get_number_of_pictures(directions,
                                             cutOff,
                                             referenceAtom)
        multipliers += extraMultiplier
        extraDirections = [(mul+1)*d 
                           for mul,d in
                           zip(multipliers,
                               directions)]
        crystal = []                   
        for name in atomNames:
            for x in range(multipliers[0]+1):
                for y in range(multipliers[1]+1):
                    for z in range(multipliers[2]+1):
                        if moments is None:
                            for atom in cell:
                                if atomNames[atom[0]] == name:
                                    position = np.copy(atom[1])
                                    for a,n in zip([x,y,z],directions):
                                        position += a*n
                                    if isinstance(atom[0],int):
                                        crystal.append([atom[0],position,elementMagneticMoment[periodicTableElement[atomNames[atom[0]]]]])    
                                    else:
                                        crystal.append([atom[0],position,elementMagneticMoment[periodicTableElement[atom[0]]]])    
                        else:
                            for atom,moment in zip(cell,moments):
                                if atomNames[atom[0]] == name:
                                    position = np.copy(atom[1])
                                    for a,n in zip([x,y,z],directions):
                                        position += a*n
                                    crystal.append([atom[0],position,moment])    

        newReference = None
        newReferenceAtom = None
        for i,atom in enumerate(crystal):
            if np.linalg.norm(atom[1] - referenceAtom[1]) < 1e-3:
                newReference = i
                newReferenceAtom = atom

        ISFOUND,distances,wyckoffPositionDict = check_in_cell(
                                                     crystal,
                                                     newReferenceAtom,
                                                     extraDirections,
                                                     nearestNeighbor,
                                                     atomNames,
                                                     cell,
                                                     directions,
                                                     cutOff,
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
            cutOff = distances[nearestNeighbor]
            for atom in crystal:
                atom[0] = atomNames[atom[0]]
            return (cutOff,
                    crystal,
                    crystalSymmetry,
                    newReference,
                    multipliers,
                    wyckoffPositionDict)

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
