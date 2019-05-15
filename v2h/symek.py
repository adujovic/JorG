#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv
from sys import maxsize
from os import system,environ
import re
from datetime import datetime
import numpy as np
import spglib
from aux.periodic import *
from aux.symmetry import *
from aux.format import *
from aux.argv import options
from JorG.loadsave import * 
from JorG.generator import * 
from multiprocessing import Pool, TimeoutError


def main(**args):
    pass

if __name__ == '__main__':
    currentOptions = options(*argv)
    POSCARfile     = currentOptions('input')
    INCARfile      = currentOptions('incar')
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
    """ Reading POSCAR and INCAR files.
          TODO: bulletproofing """
#    
    readData             = load_POSCAR(POSCARfile)
    oldMoments,incarData = load_INCAR (readData['cell'],INCARfile)
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
    
    print("The reference was chosen to be atom No. %d:"%(reference+1))
    print_atom(referenceAtom)


    """ Checking the symmetry 
                    of the input """
    symmetryCrude   = spglib.get_symmetry_dataset(cellSymmetry)
    if(SYMMETRYRUN):
        standarizedCell  = (spglib.standardize_cell(cellSymmetry,
                                               to_primitive=1,
                                               no_idealize=0,
                                               symprec=1e-1))
        symmetryStandard = spglib.get_symmetry_dataset(standarizedCell)
        refinedCell      = (spglib.refine_cell(cellSymmetry,
                                               symprec=1e-1))
        symmetryRefined = spglib.get_symmetry_dataset(refinedCell)
        write_report(["""Analysis of symmetry in:\n
            (1) the crude input cell""",
           "(2) the standarized cell",
           "(3) the refined primitive cell"],
                [symmetryCrude,symmetryStandard,symmetryRefined],
                cell, atomDict=atomNames)
        niggli = spglib.niggli_reduce(cellSymmetry)
        print(niggli)
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
        for refinedAtom in refinedCell[1]:
            newPosition = np.zeros(3)
            for x,d in zip(refinedAtom,refinedCell[0]):
                newPosition += x*np.array(d)
            for atom in cell:
                if np.linalg.norm(atom[1] - newPosition) <= 1e-1:
                    atom[1] = newPosition
                    continue
    """

        Printing input data


    """

    print_label("INPUT")
    print_crystal(directions,cell,atomNames=atomNames)
    print_moments(oldMoments,cell=cell,atomNames=atomNames)

    """ 
    
        Generating output


    """
    if nearestNeighbor is None:
        nearestNeighbor = 2

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
        extraDirections = [(mul+1)*d 
                           for mul,d in
                           zip(copiesInEachDirection,
                               directions)]
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

    print_label("OUTPUT")
    print_crystal(extraDirections,crystal)

    
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
    print("Reference atom in the new system is No. %d:"%newReference)
    print_atom(crystal[newReference],vector=color.DARKCYAN)
    for (i,atom,distance,wyck) in flipper:
        if caseID <= nearestNeighbor:
            print_case(caseID,atom,i+1,wyck,distance)
            selected.append(i)
            caseID += 1
                
#    if nearestNeighbor < len(flipper) :
#        save_INCAR(outDirName,incarData,crystal,flipper[:nearestNeighbor])    
#    else:
#        save_INCAR(outDirName,incarData,crystal,flipper)    
#
    crystal8 = apply_mirrorsXYZ(extraDirections,crystal,
                                cutOff=cutOff,
                                reference=newReference)
    save_xyz   (outDirName+"/crystal.xyz",crystal,selectedAtoms = selected)
    save_xyz   (outDirName+"/crystalFull.xyz",crystal8,selectedAtoms = selected)
    system ("sed  -e 's/XXXXX/%f/g' -e 's/YYYYY/%f/g' -e 's/ZZZZZ/%f/g' -e 's/RRRRR/%f/g' script.template> %s"%(*crystal[newReference][1],cutOff,outDirName+"/script.jmol"))

    from JorG.equivalent import find_all_flips
    allFlippable = find_all_flips(crystal[newReference],crystal,
                                symmetryFull,cutOff,
                                atomTypeMask,Wyckoffs=wyckoffs,
                                wyckoffDict=wyckoffDict,
                                logAccuracy=logAccuracy)
 
    from itertools import product
    from JorG.configurations import *
    allOptions = []
    configurations = product([1,-1],repeat=len(allFlippable))
    numberOfConfigurations = int(2**len(allFlippable))
    print("Checking total number of configurations:",numberOfConfigurations)

    randomInteger = np.random.randint(10000000,99999999)

    with open(".input%d.dat"%randomInteger,"w+") as isingFile:
        for i in allFlippable:
            isingFile.write("%d %.8f %.8f %.8f %.2f\n"%(
            i,*crystal[i][1],crystal[i][2]))

    with open(".supercell%d.dat"%randomInteger,"w+") as supercellFile:
        for i,atom in enumerate(crystal):
            supercellFile.write("%d %.8f %.8f %.8f %.2f\n"%(
            i,*atom[1],atom[2]))

    with open(".directions%d.dat"%randomInteger,"w+") as dirFile:
        for d in extraDirections:
            dirFile.write("%.8f %.8f %.8f\n"%tuple(d))

    print("\n")
    system('cd asa/solver; make clean; make SITES=-D_SITESNUMBER=%d; cd ../../'%len(crystal))
    system('echo \"Running: ./asa/solver/start .directions%d.dat .supercell%d.dat .input%d.dat %d %d\"'%(randomInteger,randomInteger,randomInteger,newReference,nearestNeighbor))
    system('./asa/solver/start .directions%d.dat .supercell%d.dat .input%d.dat %d %d'%(randomInteger,randomInteger,randomInteger,newReference,nearestNeighbor))
    system('rm .*%d.dat'%randomInteger)

    flippingConfigurations = np.loadtxt('best.flips',bool)
    save_INCAR(outDirName,incarData,crystal,flippingConfigurations)
    exit()
    
#    with Pool(processes=4) as pool:
#        output = [pool.apply_async(
#                   get_configuration_penalty, (configuration, flipper,
#                                               crystal, crystal8,
#                                               allFlippable, newReference))
#              for configuration in configurations ]
#    for entry in output:
#        try:
#            print(entry.get(timeout=1))
#        except TimeoutError:
#            print("We lacked patience and got a multiprocessing.TimeoutError")
