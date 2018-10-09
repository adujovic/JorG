#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv
import numpy as np
import spglib
from aux.periodic import *
from aux.symmetry import *
from aux.argv import options
from JorG.loadsave import * 
from JorG.generator import * 

def main(**args):
    pass

if __name__ == '__main__':
    
    currentOptions = options(*argv)
    workingFile    = currentOptions('input')
    cutOff         = currentOptions('cutOff')
    atomTypeMask   = currentOptions('mask')  
    reference      = currentOptions('reference')
    SYMMETRYRUN    = currentOptions('symmetry')
    wyckoffs       = currentOptions('Wyckoffs')
    
    """ Reading POSCAR file.
          TODO: bulletproofing """
    
    readData = load(workingFile)
    cell          = readData['cell']
    cellSymmetry  = readData['cellSymmetry']
    atomNames     = readData['atomNames']
    comment       = readData['comment']
    directions    = readData['directions']
    cellVolume    = readData['cellVolume']
    cellCenter    = readData['cellCenter']
    cellAtoms     = readData['cellAtoms']
    
    referenceAtom = cell[reference]
    
    if cutOff is None:
        cutOff = np.max(np.linalg.norm(directions))
    
    print("The reference was chosen to be atom No. %d"%reference)
    print("%s %f %f %f"%(atomNames[referenceAtom[0]],
                                   referenceAtom[1][0],
                                   referenceAtom[1][1],
                                   referenceAtom[1][2]))
    
    """ Checking the symmetry 
                    of the input """
    symmetryCrude   = spglib.get_symmetry_dataset(cellSymmetry)
    
    if(SYMMETRYRUN):
        refinedCell     = (spglib.standardize_cell(cellSymmetry,
                                               to_primitive=1,
                                               no_idealize=0,
                                               symprec=1e-1))
        symmetryRefined = spglib.get_symmetry_dataset(refinedCell)
        write_report("(1) crude cell\n(2) refined cell",
                 [symmetryCrude,symmetryRefined], cell, atomDict=atomNames)
    else:
        write_report("cell", [symmetryCrude], cell,
                     "output/input_report.txt", atomDict=atomNames);
    
    """ Generating output """
    
    copiesInEachDirection = get_number_of_pictures(directions,cutOff)
    
    crystal, cellSymmetryFull, newReference = generate_crystal(copiesInEachDirection,cell,directions,atomNames,reference=reference)
    
    """ Checking the symmetry 
                    of the output """
    symmetryFull = spglib.get_symmetry_dataset(cellSymmetryFull)
    
    if(SYMMETRYRUN):
        write_report("(3) crystal",
                 [symmetryFull], crystal)
        exit(1)
    
    write_report("generated cell",[symmetryFull],crystal,"output/output_report.txt");
    save_POSCAR("output/POSCAR",crystal,copiesInEachDirection,readData)
    
    """
        Searching for unique atoms for calculations
                                                    """
    logAccuracy  = 2     #  accuracy of distance is 10^(-logAccuracy)
    
    from JorG.equivalent import find_unique_flips
    flipper = find_unique_flips(crystal[newReference],crystal,symmetryFull,cutOff,atomTypeMask,Wyckoffs=wyckoffs,logAccuracy=logAccuracy)
    
    
    output = ""
    caseID = 1
    selected = [newReference]
    for i,rec in enumerate(flipper):
        if rec[0]:
            selected.append(i)
            output += "Case  %d:\t atom No. %d @ %s "%(caseID,i+1,rec[1]) 
            caseID += 1
            output += str(rec[2][0])
            for x in rec[2][1]:
                output += ("  %.10f"%x)
            output += "\n"
    print(output)
                
    save_xyz   ("output/crystal.xyz",crystal,selectedAtoms = selected)
