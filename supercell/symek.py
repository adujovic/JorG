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
    nearestNeighbor= currentOptions('neighbor')  
    atomTypeMask   = currentOptions('mask')  
    reference      = currentOptions('reference')
    SYMMETRYRUN    = currentOptions('symmetry')
    USEREFINED     = currentOptions('refined')
    wyckoffs       = currentOptions('Wyckoffs')
    
    """ Reading POSCAR file.
          TODO: bulletproofing """
#    
    readData = load(workingFile)
#    
    cell          = readData['cell']
    cellSymmetry  = readData['cellSymmetry']
    atomNames     = readData['atomNames']
    comment       = readData['comment']
    directions    = readData['directions']
    cellVolume    = readData['cellVolume']
    cellCenter    = readData['cellCenter']
    cellAtoms     = readData['cellAtoms']
    
    referenceAtom = cell[reference]
    
    print("The reference was chosen to be atom No. %d: %s %f %f %f"%(
          reference+1, 
          atomNames[referenceAtom[0]],
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
        write_report(["""Analysis of symmetry in:\n
            (1) the crude input cell""",
           "(2) the newly refined input cell"],
                [symmetryCrude,symmetryRefined],
                cell, atomDict=atomNames)
        exit(0)
    else:
        write_report(["Analysis of symmetry in the input cell"], [symmetryCrude], cell,
                     "output/input_report.txt", atomDict=atomNames);

    if USEREFINED: 
        refinedCell = (spglib.standardize_cell(cellSymmetry,
                                               to_primitive=1,
                                               no_idealize=0,
                                               symprec=1e-1))
        directions = np.array(refinedCell[0])
        for atom,refinedAtom in zip(cell,refinedCell[1]):
            newPosition = np.zeros(3)
            for x,d in zip(refinedAtom,refinedCell[0]):
                newPosition += x*np.array(d)
            atom[1] = newPosition

   

    """ Generating output """

    if cutOff is None:
        if nearestNeighbor is None:
            nearestNeighbor = 1
        (cutOff,
         crystal,
         symmetryFull, 
         newReference, 
         copiesInEachDirection,
         wyckoffDict) = generate_from_NN(cell,
                                         referenceAtom,
                                         directions,
                                         nearestNeighbor,
                                         atomNames,
                                         wyckoffs,
                                         atomTypeMask)
    else:
        copiesInEachDirection = get_number_of_pictures(directions,cutOff)
        extraDirections = [mul*d 
                           for mul,d in
                           zip(copiesInEachDirection,
                               directions)]
        crystal, newReference =\
                generate_crystal(copiesInEachDirection,
                                 cell,
                                 directions,
                                 atomNames,
                                 reference=reference)
        wyckoffDict, symmetryFull = wyckoffs_dict(cell, 
                                                      crystal,
                                                      directions,
                                                      extraDirections,
                                                      atomNames)        
#    
    """ Checking the symmetry 
                    of the output """
    write_report(["Analysis of symmetry in the generated cell"],
                 [symmetryFull],
                 crystal,
                 "output/output_report.txt");
    save_POSCAR("output/POSCAR",
                crystal,
                copiesInEachDirection,
                readData)
    
    """
        Searching for unique atoms for calculations
                                                    """
    logAccuracy  = 2     #  accuracy of distance is 10^(-logAccuracy)
    
    from JorG.equivalent import find_unique_flips
    flipper = find_unique_flips(crystal[newReference],crystal,
                                symmetryFull,cutOff,
                                atomTypeMask,Wyckoffs=wyckoffs,
                                wyckoffDict=wyckoffDict,
                                logAccuracy=logAccuracy)
    
    
    output = ""
    caseID = 1
    selected = [newReference]
    for i,rec in enumerate(flipper):
        if rec[0]:
            selected.append(i)
            output += "Case  %d:\t atom No. %d @ %s "%(caseID,i+1,rec[2]) 
            caseID += 1
            output += str(rec[1][0])
            for x in rec[1][1]:
                output += ("  %.10f"%x)
            output += "\n"
    print(output)
                
    save_xyz   ("output/crystal.xyz",crystal,selectedAtoms = selected)
