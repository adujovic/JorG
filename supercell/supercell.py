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
atomTypeMask = maskD     # all  d elements 
atomTypeMask = maskF     # all  f elements 
atomTypeMask = mask2p    # all 2p elements 
atomTypeMask = mask3d    # all 3d elements 
atomTypeMask = maskFull  # all    elements 
referenceAtom = cell[reference]

from JorG.equivalent import find_unique_flips
flipper = find_unique_flips(referenceAtom,crystal,symmetryFull,cutOff,atomTypeMask,logAccuracy)

output = ""
for i,rec in enumerate(flipper):
    if rec[0]:
        output += "%d\t"%(i+1) 
        output += str(rec[1][0])
        for x in rec[1][1]:
            output += ("  %.10f"%x)
        output += "\n"
print(output)
            
