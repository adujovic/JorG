#!/usr/bin/python3

from sys import argv
import numpy as np
import spglib
from aux.periodic import *
from aux.symmetry import *
from JorG.loadsave import * 

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
             [symmetryCrude,symmetryRefined],
             "output/input_report.txt");

multiplyers = get_number_of_pictures(directions,cutOff)

crystal, cellSymmetryFull = generate_crystal(multiplyers,cell,directions,atomNames)

write_xyz   ("output/crystal.xyz",crystal)
write_POSCAR("output/POSCAR",crystal,multiplyers,readData)

symmetryFull = spglib.get_symmetry_dataset( cellSymmetryFull)
write_report("generated cell",[symmetryFull],"output/output_report.txt");

