#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,path
path.insert(0,r'../../')
from os import makedirs
from subprocess import call
import errno
import re
from datetime import datetime
import numpy as np
import spglib
import time
import shutil
from itertools import product

import JorG.symmetry
from JorG.argv import options
from JorG.format import color, print_case, print_vector
from JorG.format import print_crystal, print_moments, print_label
import JorG.loadsave as loadsave
import JorG.generator as generator

from JorG.pickup import EquationSolver,NaiveHeisenberg

from WiP.WiP import StreamHandler,JmolVisualization,TemporaryFiles

class errors:
    failed_to_generate = -201
    systemerror = -1

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
    extraDimentions= currentOptions('extra-dimentions')

    handler    = StreamHandler(outDirName)
    outDirName = handler()

#    """ Reading POSCAR and INCAR files
#                                       """

    load_POSCAR          = loadsave.POSCARloader(POSCARfile)
    load_POSCAR.parse()
    readData             = load_POSCAR(0)
    oldMoments,incarData = loadsave.load_INCAR (readData['cell'],INCARfile,atomNames=readData['atomNames'])

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
        if "$"+atom[0]+"$" in atomTypeMask:
          referenceAtom = atom
          reference = i
          break
    if referenceAtom is None:
      print("Error: can not find any atoms (%s) in input file!"%re.sub('\$','',atomTypeMask))
      exit(-7)

#    """ Checking the symmetry
#                    of the input """
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
        JorG.symmetry.write_report(["(1) the crude input cell",
                      "(2) the standarized cell",
                      "(3) the refined primitive cell"],
                [symmetryCrude,symmetryStandard,symmetryRefined],
                cell)
        exit(0)
    else:
        with open(outDirName+"/input_report.txt",'w+') as raport:
            JorG.symmetry.write_report(["Analysis of symmetry in the input cell"], [symmetryCrude], cell,stream=raport)
                     

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
#    """
#        Printing input data
#    """

    print_label("The reference was chosen to be atom No. %d:"%(reference+1),atoms=[referenceAtom],labelStyle=color.BF)

    print_label("INPUT",labelStyle=color.BF)
    print_crystal(directions,cell)
    print_moments(oldMoments,cell=cell)

#    """
#        Generating output
#    """
    if extraDimentions is not None:
        for separator in [' ', ',', ';', '.', '/']:
            extraMultiplier = np.fromstring(extraDimentions, dtype=int, sep=separator)
            if len(extraMultiplier) >= 3:
                extraMultiplier = extraMultiplier[:3]
                break
        if len(extraMultiplier) != 3:
            print("Error with the extraDimentions given. Omitting this functionality!")
            extraMultiplier = np.zeros(3,dtype=int)
    else:
        extraMultiplier = np.zeros(3,dtype=int)

    if nearestNeighbor is None:
        nearestNeighbor = -1

    if cutOff is None:
        if nearestNeighbor is -1:
            nearestNeighbor = 2

        generatorNN = generator.NearestNeighborsGenerator(cell,
                                     referenceAtom,
                                     directions)
        generatorNN.wyckoffs         = wyckoffs
        generatorNN.atomTypeMask     = atomTypeMask
        generatorNN.moments          = oldMoments
        generatorNN.extraMultiplier  = extraMultiplier

        try:
            (cutOff,
             crystal,
             symmetryFull,
             newReference,
             copiesInEachDirection,
             wyckoffDict           ) = generatorNN(nearestNeighbor)
        except Exception:
            print("Failed to generate crystal")
            exit(errors.failed_to_generate)

        extraDirections = [(mul+1)*d
                           for mul,d in
                           zip(copiesInEachDirection,
                               directions)]
    else:
        copiesInEachDirection = generator.get_number_of_pictures(directions,cutOff,referenceAtom)
        extraDirections = [(mul+extra+1)*d
                           for mul,extra,d in
                           zip(copiesInEachDirection,
                               extraMultiplier,
                               directions)]
        localGenerator = generator.CrystalGenerator(cell,directions,
                                                    atomNames,reference=reference)
        localGenerator.moments = oldMoments
        crystal, newReference = localGenerator(copiesInEachDirection)
        wyckoffDict, symmetryFull, symmetryOriginal = generator.wyckoffs_dict(cell,
                                                      crystal,
                                                      directions,
                                                      extraDirections)
#    """ Checking the symmetry
#                    of the output """
    with open(outDirName+"/output.txt",'w+') as raport:
        JorG.symmetry.write_report(["Analysis of symmetry in the generated cell"],
                     [symmetryFull],
                     crystal,
                     stream=raport)
    loadsave.save_POSCAR(outDirName+"/POSCAR",
                crystal,
                copiesInEachDirection+extraMultiplier,
                readData)
    realCopies = [copy+1 for copy in copiesInEachDirection]
    print_label("OUTPUT: %dx%dx%d"%(*realCopies,),labelStyle=color.BF)
    print_crystal(extraDirections,crystal)

#   """
#        Searching for unique atoms for calculations
#                                                    """
    logAccuracy  = 2     #  accuracy of distance is 10^(-logAccuracy)

    from JorG.equivalent import findFlips
    flipSearch = findFlips()
    flipSearch.symmetry     = symmetryFull
    flipSearch.crystal      = crystal
    flipSearch.atomTypeMask = atomTypeMask
    flipSearch.Wyckoffs     = wyckoffs
    flipSearch.wyckoffDict  = wyckoffDict
    flipSearch.logAccuracy  = logAccuracy
    flipper      = flipSearch.unique(crystal[newReference],cutOff)
    allFlippable = flipSearch.all(crystal[newReference],cutOff)

    caseID = 1
    selected = [newReference]

    print_label("Reference atom in the new system is No. %d:"%(newReference+1),atoms=[crystal[newReference]],vectorStyle=color.DARKCYAN,labelStyle=color.BF)

    for (i,atom,distance,wyck) in flipper:
        if caseID <= nearestNeighbor or nearestNeighbor < 0:
            print_case(caseID,atom,i+1,wyckoffPosition=wyck,distance=distance)
            selected.append(i)
            caseID += 1

    nearestNeighbor = caseID-1 

    crystal8 = generator.apply_mirrorsXYZ(extraDirections,crystal,
                                cutOff=cutOff,
                                reference=newReference)
    loadsave.save_xyz(outDirName+"/crystal.xyz",crystal,selectedAtoms = selected)
    loadsave.save_xyz(outDirName+"/crystalFull.xyz",crystal8,selectedAtoms = selected)
    JmolVisualization.create_script(outDirName,radius=cutOff,center=crystal[newReference][1])
    
    print("")
    numberOfConfigurations = int(2**len(allFlippable))
    print_label("Checking total number of configurations: %d"%numberOfConfigurations,labelStyle=color.BF+color.DARKRED)
    print_label("Preparing solver...",labelStyle=color.BF+color.BLUE)

    tmpFiles = TemporaryFiles()
    tmpFiles.write_input(allFlippable,crystal)
    tmpFiles.write_supercell(crystal)
    tmpFiles.write_directions(extraDirections)

    print("")
    call('cd asa/solver; make clean; make SITES=-D_SITESNUMBER=%d; cd ../../'%len(crystal), shell=True)
    print('Running: ./asa/solver/start %s %d %d'%(str(tmpFiles),newReference,4*nearestNeighbor+8))
    call('./asa/solver/start %s %d %d'%(str(tmpFiles),newReference,4*nearestNeighbor+8), shell=True)
    call('cd ./asa/solver; make clean; cd ../../', shell=True)
    del tmpFiles

    flippingConfigurations = np.loadtxt('best.flips',bool)
    try:
        np.shape(flippingConfigurations)[1]
    except IndexError:
        flippingConfigurations=[flippingConfigurations]


    gen = NaiveHeisenberg(flippingConfigurations,crystal,crystal8)
    systemOfEquations = gen.generate(atomTypeMask,[flip[2] for flip in flipper])

    eqs = EquationSolver(systemOfEquations,np.zeros(len(systemOfEquations)))
    systemOfEquations = eqs.remove_tautologies()
    remover = eqs.remove_linears()
    if remover:
        flippingConfigurations = np.delete(flippingConfigurations, tuple(remover), axis=0)
        systemOfEquations = eqs.equations
    if systemOfEquations.size == 0:
        print("ERROR! Not enough equations! Please rerun.")
        exit(-3)
    if not currentOptions('redundant'): # If the System of Equations is required to be consistent
        systemOfEquations,flippingConfigurations = eqs.remove_linear_combinations(flippingConfigurations)

    print_label("System of equations:",labelStyle=color.BF)
    for eq in systemOfEquations:
        print_vector(eq)
    if currentOptions('redundant'):
        print_label("Redundant system of equations.",labelStyle=color.BF)
        print_label("Least square method is to be used to obtain Heisenberg model.",labelStyle=color.BF)
        print_label("It may be better. But it may also mess everything.",labelStyle=color.BF)
    else:
        print_label("det SoE = %.1e"%np.linalg.det(systemOfEquations),labelStyle=color.BF)

    np.savetxt(outDirName+'/systemOfEquations.txt',systemOfEquations)
    saver = loadsave.INCARsaver(incarData,crystal)
    saver.save(outDirName,flippingConfigurations)
    exit()
