from itertools import permutations as per
def get_number_of_pictures(directions,cutOff):
    """
        Finding the amount of copies
        of the original SuperCell
                                        """
    multipliers = [] # returned 
    
    bestScore = 0.0             # starting score
    bestPermutation = (0,1,2)   # starting permutation
    for p in per(range(3)):
        score = 1.0
        for i,d in enumerate(directions):
            score *= d[p[i]]
        if np.abs(score) > bestScore: # best score wins
            bestScore = np.abs(score)
            bestPermutation = p
    
    for d,p in zip(directions,
                   bestPermutation):
        multipliers.append(int(cutOff/d[p])) # calculating multipliers
    
    return multipliers

#
#
#
#
#

import numpy as np
from aux.periodic import periodicTableNumber
def generate_crystal(multiplyers,cell,directions,atomNames):
    """
        Generator of minimal required SuperCell
                                                """
    crystal = []
    cellSymmetryFull = ([],[],[])

    invDirections = np.dot(np.diag(1.0/(1.0+np.array(multiplyers))),np.linalg.inv(directions))

    for m,d in zip(multiplyers,directions):
        cellSymmetryFull[0].append((m+1)*d)

    for name in atomNames:
        for x in range(multiplyers[0]+1):
            for y in range(multiplyers[1]+1):
                for z in range(multiplyers[2]+1):
                    for atom in cell:
                        if atomNames[atom[0]] == name:
                            position = np.copy(atom[1])
                            for a,n in zip([x,y,z],directions):
                                position += a*n
                            flag = "%s"%atomNames[int(atom[0])]
                            crystal.append([flag,position])    
    for atomName in atomNames:
        for atom in crystal:
            if atom[0]==atomName:
                cellSymmetryFull[1].append(tuple(np.dot(atom[1],invDirections)))
                cellSymmetryFull[2].append(periodicTableNumber[atom[0]])

    return crystal,cellSymmetryFull
