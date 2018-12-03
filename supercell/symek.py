#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv
from sys import maxsize
from os import system
import re
from datetime import datetime
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
    outDirName     = currentOptions('output')
   
    if outDirName == None:
      # if output directory is not given:  
      outDirName = "output/"+datetime.now().strftime("%Y%m%d%H%M%S")
    else:
      # remove multiple '/' and possible '/' at the end  
      outDirName = re.sub('/+','/',outDirName)
      outDirName = re.sub('/$','',outDirName)

    # cr4eating output path
    temporaryName = ""
    for partOfOutput in re.split('/',outDirName):
      temporaryName += partOfOutput
      system("mkdir -p %s"%temporaryName)
      temporaryName += "/"

    # clean output path ? SHOULD WE?
#    system("rm -r "+temporaryName+"*")
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

    referenceAtom = None
    if reference >= 0:
      referenceAtom = cell[reference]
    else:
      for i,atom in enumerate(cell):  
        if "$"+atomNames[atom[0]]+"$" in atomTypeMask:  
          referenceAtom = atom
          reference = i
          break
    if referenceAtom is None:
      print("Error: can not find any atoms (%s) in input file!"%re.sub('\$','',atomTypeMask))
      exit(-7)
    
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
                     outDirName+"/input_report.txt", atomDict=atomNames);

    if USEREFINED: 
        refinedCell = (spglib.standardize_cell(cellSymmetry,
                                               to_primitive=0,
                                               no_idealize=0,
                                               symprec=1e-1))
        directions = np.array(refinedCell[0])
        for atom,refinedAtom in zip(cell,refinedCell[1]):
            newPosition = np.zeros(3)
            for x,d in zip(refinedAtom,refinedCell[0]):
                newPosition += x*np.array(d)
            atom[1] = newPosition
   
    with open("INCAR","r") as INCARfile:
        incarData = INCARfile.read()
     
    oldMomentsText = re.search("\s*MAGMOM\s*=\s*(.*)\n",incarData)
    oldMoments = []
    print("magnetic moments read:")
    print(oldMomentsText.group(0))
    print("--------------------------------------------------------------------------------------------------------------------------------")
    if oldMomentsText is None:
        for atom in cell:
            oldMoments.append(elementMagneticMoment[atomNames[atom[0]]])
    else:
        for moment in oldMomentsText.group(1).split():
            oldMoments.append(np.float(moment))

    """ Generating output """
    if nearestNeighbor is None:
        nearestNeighbor = maxsize

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
                                         atomTypeMask,
                                         moments=oldMoments)
    else:
        copiesInEachDirection = get_number_of_pictures(directions,cutOff,referenceAtom)
        extraDirections = [(mul+1)*d 
                           for mul,d in
                           zip(copiesInEachDirection,
                               directions)]
        crystal, newReference =\
                generate_crystal(copiesInEachDirection,
                                 cell,
                                 directions,
                                 atomNames,
                                 reference=reference,
                                 moments=oldMoments)
        wyckoffDict, symmetryFull, symmetryOriginal = wyckoffs_dict(cell, 
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
                 outDirName+"/output_report.txt");
    save_POSCAR(outDirName+"/POSCAR",
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
    
    caseID = 1
    selected = [newReference]
    print("Reference atom in new system is %s at %f %f %f"%(crystal[newReference][0],*crystal[newReference][1]))
    for (i,atom,distance,wyck) in flipper:
        if caseID <= nearestNeighbor:
            selected.append(i)
            print("Case  %d:\t atom No. %d @ %s && %3.2f A | %s  %.5f %.5f %.5f"%(caseID,i+1,wyck,distance,atom[0],*atom[1]))
            caseID += 1
                
    if nearestNeighbor < len(flipper) :
        save_INCAR(outDirName,incarData,crystal,flipper[:nearestNeighbor])    
    else:
        save_INCAR(outDirName,incarData,crystal,flipper)    

    save_xyz   (outDirName+"/crystal.xyz",crystal,selectedAtoms = selected)
    save_xyz   (outDirName+"/crystalFull.xyz",apply_mirrorsXYZ(extraDirections,crystal),selectedAtoms = selected)
    system ("sed  -e 's/XXXXX/%f/g' -e 's/YYYYY/%f/g' -e 's/ZZZZZ/%f/g' -e 's/RRRRR/%f/g' script.template> %s"%(*crystal[newReference][1],cutOff,outDirName+"/script.jmol"))
