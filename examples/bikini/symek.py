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
from JorG.equivalent import findFlips

from JorG.pickup import EquationSolver,NaiveHeisenberg

from WiP.WiP import StreamHandler,JmolVisualization
from WiP.WiP import TemporaryFiles,errors
from WiP.WiP import VariableFixer

def main(**args):
    pass

if __name__ == '__main__':
    currentOptions = options(*argv)
    cutOff         = currentOptions('cutOff')
    nearestNeighbor= currentOptions('neighbor')
    reference      = currentOptions('reference')
    outDirName     = currentOptions('output')
    extraMultiplier = np.zeros(3,dtype=int)

    handler    = StreamHandler(outDirName)
    outDirName = handler()

#    """ Reading POSCAR and INCAR files
#                                       """
    readData,oldMoments,incarData =\
            handler.load_VASP(currentOptions('input'),currentOptions('incar'))
    cell          = readData['cell']
    cellSymmetry  = readData['cellSymmetry']
    atomNames     = readData['atomNames']
    comment       = readData['comment']
    directions    = readData['directions']
    referenceAtom,reference = VariableFixer.fix_reference(reference,cell,
                                                          currentOptions('mask'))
#    """ Checking the symmetry
#                    of the input """
    symmetryCrude   = spglib.get_symmetry_dataset(cellSymmetry)

    if(currentOptions('symmetry')):
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
                     

    if currentOptions('refined'):
        refinedCell = (spglib.standardize_cell(cellSymmetry, to_primitive=0, no_idealize=0, symprec=1e-1))
        cell,directions = VariableFixer.from_refined(refinedCell) 
#    """
#        Printing input data
#    """
    print_label("The reference was chosen to be atom No. %d:"%(reference+1),atoms=[referenceAtom],labelStyle=color.BF)
    print_label("INPUT",labelStyle=color.BF)
    print_crystal(directions,cell)
    print_moments(oldMoments,cell=cell)

    nearestNeighbor,cutOff = VariableFixer.fix_neighbor(nearestNeighbor,cutOff)

    if cutOff is None:
        generatorNN = generator.NearestNeighborsGenerator(cell,
                                     referenceAtom,
                                     directions)
        generatorNN.wyckoffs         = currentOptions('Wyckoffs')
        generatorNN.atomTypeMask     = currentOptions('mask')
        generatorNN.moments          = oldMoments
        generatorNN.extraMultiplier  = extraMultiplier

        try:
            (cutOff, crystal, symmetryFull, newReference,
                      copiesInEachDirection, wyckoffDict) = generatorNN(nearestNeighbor)
        except Exception:
            print("Failed to generate crystal")
            exit(errors.failed_to_generate)
        extraDirections = VariableFixer.fix_directions(copiesInEachDirection,directions)
    else:
        copiesInEachDirection = generator.get_number_of_pictures(directions,cutOff,referenceAtom)
        localGenerator = generator.CrystalGenerator(cell,directions, atomNames,reference=reference)
        localGenerator.moments = oldMoments
        crystal, newReference = localGenerator(copiesInEachDirection)
        extraDirections = VariableFixer.fix_directions(copiesInEachDirection,directions)
        wyckoffDict, symmetryFull, symmetryOriginal =\
                generator.wyckoffs_dict(cell, crystal, directions, extraDirections)

#    """ Checking the symmetry
#                    of the output """
    with open(outDirName+"/output.txt",'w+') as raport:
        JorG.symmetry.write_report(["Analysis of symmetry in the generated cell"],
                     [symmetryFull], crystal, stream=raport)
    loadsave.save_POSCAR(outDirName+"/POSCAR", crystal,
                copiesInEachDirection, readData)

    print_label("OUTPUT: %dx%dx%d"%(*VariableFixer.add_to_all(copiesInEachDirection),),labelStyle=color.BF)
    print_crystal(extraDirections,crystal)

#   """
#        Searching for unique atoms for calculations
#                                                    """
    logAccuracy  = 2     #  accuracy of distance is 10^(-logAccuracy)

    flipSearch              = findFlips()
    flipSearch.symmetry     = symmetryFull
    flipSearch.crystal      = crystal
    flipSearch.atomTypeMask = currentOptions('mask')
    flipSearch.Wyckoffs     = currentOptions('Wyckoffs')
    flipSearch.wyckoffDict  = wyckoffDict
    flipSearch.logAccuracy  = logAccuracy
    flipper                 = flipSearch.unique(crystal[newReference],cutOff)
    allFlippable            = flipSearch.all(crystal[newReference],cutOff)

    selected = [newReference]

    print_label("Reference atom in the new system is No. %d:"%(newReference+1),atoms=[crystal[newReference]],vectorStyle=color.DARKCYAN,labelStyle=color.BF)

    for caseID,(i,atom,distance,wyck) in enumerate(flipper):
        print_case(caseID+1,atom,i+1,wyckoffPosition=wyck,distance=distance)
        selected.append(i)

    crystal8 = generator.apply_mirrorsXYZ(extraDirections,crystal,
                                cutOff=cutOff, reference=newReference)
    loadsave.save_xyz(outDirName+"/crystal.xyz",crystal,selectedAtoms = selected)
    loadsave.save_xyz(outDirName+"/crystalFull.xyz",crystal8,selectedAtoms = selected)
    JmolVisualization.create_script(outDirName,radius=cutOff,center=crystal[newReference][1])
    
    print("")
    
    print_label("Checking total number of configurations: %d"%(int(2**len(allFlippable))),labelStyle=color.BF+color.DARKRED)
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
    systemOfEquations = gen.generate(currentOptions('mask'),[flip[2] for flip in flipper])

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
