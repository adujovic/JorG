#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv
import numpy as np
import spglib
from aux.periodic import *
from aux.symmetry import *
from JorG.loadsave import * 
from JorG.generator import * 

if len(argv) == 1:
    workingFile = "POSCAR"
    cutOff = 0.0
elif len(argv) == 2:
    try:
        cutOff = float(argv[1])
    except:
        print("Broken input:\n\tCan't convert %s to float"%argv[1])
        raise
else:        
    workingFile = "POSCAR"
    try:
        cutOff = float(argv[1])
    except:
        print("Broken input:\n\tCan't convert %s to float"%argv[1])
        raise
    try:
        workingFile = argv[2]
        test = open(workingFile,'r+')
    except:
        print("Broken input:\n\tCan't open file %s"%workingFile)
        raise
    test.close()

""" Reading POSCAR file.
      TODO: bulletproofing """

readData = load("POSCAR")
comment       = readData['comment']
directions    = readData['directions']
cell          = readData['cell']
cellSymmetry  = readData['cellSymmetry']
cellVolume    = readData['cellVolume']
cellCenter    = readData['cellCenter']
cellAtoms     = readData['cellAtoms']
atomNames     = readData['atomNames']

symmetryCrude = spglib.get_symmetry_dataset( cellSymmetry )
refinedCell = (spglib.standardize_cell(cellSymmetry,
                                       to_primitive=1,
                                       no_idealize=0,
                                       symprec=1e-1))
symmetryRefined = spglib.get_symmetry_dataset( refinedCell )
write_report("(1) crude cell\n(2) refined cell",
             [symmetryCrude,symmetryRefined], cell,
             "output/input_report.txt");

copiesInEachDirection = get_number_of_pictures(directions,cutOff)

crystal, cellSymmetryFull = generate_crystal(copiesInEachDirection,cell,directions,atomNames)

save_xyz   ("output/crystal.xyz",crystal)
save_POSCAR("output/POSCAR",crystal,copiesInEachDirection,readData)

symmetryFull = spglib.get_symmetry_dataset( cellSymmetryFull)

write_report("generated cell",[symmetryFull],crystal,"output/output_report.txt");

"""
    Searching for unique atoms for calculations
                                                """

# Distances with supercell:
logAccuracy  = 2     # should be given by user, accuracy of distance is 10^(-logAccuracy)
reference    = 0     # should be given by user, number of atom in input cell ∈[0,N)
atomTypeMask = "MnO" # should be given by user, consists of elemental symbols

# several masks in aux.periodic; examples:
from aux.periodic import maskFull,maskD,maskF,mask2p,mask3d
atomTypeMask = maskFull  # all    elements 
atomTypeMask = maskD     # all  d elements 
atomTypeMask = maskF     # all  f elements 
atomTypeMask = mask2p    # all 2p elements 
atomTypeMask = mask3d    # all 3d elements 
referenceAtom = cell[reference]
distances = [ ]
for atom in crystal:
   d = np.around(np.linalg.norm(atom[1]-referenceAtom[1]),decimals=logAccuracy )
   if (d,atom[0]) not in distances and d > 0.0 and d <= cutOff:
       distances.append((d,atom[0]))

distances = np.array(distances, dtype=[('distance', np.float), ('element', 'U3')])
distances.sort(order='distance')

                                    #for d,n in distances:
                                    #    print(n,d)
                                    #for atom in crystal:
                                    #  print("%s %f %f %f"%(atom[0],atom[1][0],atom[1][1],atom[1][2]))

checker = [False for x in range(len(crystal))] 
case = 1
for (d,n) in distances:
    if n in atomTypeMask:
        for i,atom in enumerate(crystal):
            if not checker[i]:
                if d == np.around(np.linalg.norm(atom[1]-referenceAtom[1]),decimals=logAccuracy) and atom[0] == n:
                    for j in np.argwhere(symmetryFull['equivalent_atoms'] == symmetryFull['equivalent_atoms'][i]).flatten():
                        if d == np.around(np.linalg.norm(crystal[j][1]-referenceAtom[1]),decimals=logAccuracy):
                            print("Case no: %d, for %s with d= %f"%(case,crystal[j][0],d))
                            checker[j] = True
                    if checker[i]:
                        case += 1
print(checker)
