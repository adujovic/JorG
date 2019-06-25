# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import re
import numpy as np
from os import makedirs,mkdir,rmdir
import errno
import shutil
import aux.PeriodicTable as periodic

class error:
    systemerror = -1
    unexcepted     = 10
    nonconvertable = 11
    vaspError      = 99


class INCARloader:
    settings = {'fileName'  : "INCAR"}

    def __init__(self,cell,**kwargs):
#    """
#        Loading INCAR
#                                        """
        self.settings.update(kwargs)
        self.INCARfile = open(self.settings['fileName'],"r")
        self.oldMoments = []
        self.incarData = self.INCARfile.read()
        self.cell      = cell

        self.oldMomentsText = re.search("\s*MAGMOM\s*=\s*(.*)\n",self.incarData)
        self.oldMoments = []
        if self.oldMomentsText is None:
            self.oldMoments = [periodic.elementMagneticMoment[atom[0]] for atom in self.cell]

    def __call__(self):
        if self.oldMoments:
            return self.oldMoments,self.incarData

        magmomLine = self.oldMomentsText.group(1)
        for record in re.findall("\s*([0-9]+)\s*\*\s*([\-0-9\.]+)",magmomLine):
            magmomLine = re.sub("{:s}\s*\*\s*{:s}".format(record[0],record[1]),(record[1]+" ")*int(record[0]),magmomLine)
        for moment in magmomLine.split():
            self.oldMoments.append(np.float(moment))
        return self.oldMoments,self.incarData

    def __del__(self):
        self.INCARfile.close()
#
#
#
#
#

from itertools import product
class save_xyz:
    settings = {'fileName'     : 'data4jmol.xyz',
                'numberOfAtoms': -1,
                'selectedAtoms': None}

    def __init__(self,crystal,**kwargs):
#    """
#        Saving data to xyz-style file
#                                        """
        self.settings.update(kwargs)
        if self.settings['numberOfAtoms'] < 0:
            self.settings['numberOfAtoms'] = len(crystal)

        self.xyzFile = open(self.settings['fileName'],"w+")
        self.xyzFile.write(str(self.settings['numberOfAtoms']))
        self.xyzFile.write("\n\n")
        for i,atom in enumerate(crystal):
            self.write_line(i,atom)
        self.xyzFile.write("\n")

    def write_line(self,i,atom):
        self.xyzFile.write("%s"%atom[0])
        self.xyzFile.write(" %.10f %.10f %.10f"%(*atom[1],))
        vector = 2*np.random.ranf(3)-1.0
        vector = 0.2*vector/np.linalg.norm(vector)
        try:
            if i == self.settings['selectedAtoms'][0]:
                vector = 2.0*vector
            elif i in self.settings['selectedAtoms']:
                self.xyzFile.write(" PatrialCharge(1.0) %f %f %f"%tuple(vector))
        except TypeError:
            print("Not selected!")
        self.xyzFile.write("\n")

    def __del__(self):
        self.xyzFile.close()
#
#
#
#
#
def save_vanilla_POSCAR(fileName,data):
#    """
#        Saving data to POSCAR file
#{'comment': 'Re/117 - (A3) - HCP [A2] A3 (117) (htqc', 'directions': [array([ 1.38476066, -2.39847581, -0.        ]), array([1.38476066, 2.39847581, 0.        ]), array([0.        , 0.        , 4.47184109])], 'cell': [[0, array([0., 0., 0.])], [0, array([ 1.38476066, -0.79949194,  2.23592055])]], 'cellSymmetry': ([(1.3847606573816698, -2.398475814907534, -0.0), (1.3847606573816698, 2.398475814907534, 0.0), (0.0, 0.0, 4.471841090100609)], [(0.0, -0.0, -0.0), (0.6666666666666714, 0.3333333333333428, 0.5)], [75, 75]), 'cellVolume': 29.704785298855388, 'cellCenter': array([1.38476066, 0.        , 2.23592055]), 'cellAtoms': array([2]), 'atomNames': ['Re']}
#                                    """
    with open(fileName,"w+") as vaspFile:
        vaspFile.write(data['comment'])
        vaspFile.write("\n1.0\n")
        for d in data['directions']:
            for field in d:
                vaspFile.write("  %.10f"%(field))
            vaspFile.write("\n")
        for atomName in data['atomNames']:
            vaspFile.write("%s "%atomName)
        vaspFile.write("\n")
        for atomNumber in data['cellAtoms']:
            vaspFile.write("%d "%atomNumber)
        vaspFile.write("\nDirect\n")
        for atom in data['cell']:
            for vasp in atom[1]:
                vaspFile.write(" %.10f "%vasp)
            try:
                vaspFile.write(" %s\n"%data['atomNames'][atom[0]])
            except TypeError:
                vaspFile.write(" %s\n"%atom[0])
        vaspFile.write("\n")

#
#
#
#
class save_POSCAR:
    settings = {'fileName'    : 'POSCAR',
                'multiplyers' : [1, 1, 1],
                'crystal'     : None}

    def __init__(self,data,**kwargs):
#    """
#        Saving data to POSCAR file
#                                    """
        self.settings.update(kwargs)
        if self.settings['crystal'] is None:
            self.settings['crystal'] = data['cell']

        self.multiplyers = np.array(self.settings['multiplyers']) + 1
        self.vaspFile = open(self.settings['fileName'],"w+")

        self.vaspFile.write(data['comment'])
        self.vaspFile.write("\n")

        self.write_directions(data['directions'])
        self.write_elements(data['atomNames'])
        self.write_populations(data['cellAtoms'])

        self.vaspFile.write("Carthesian\n")

        self.write_atoms(data['atomNames'])

    def write_directions(self,directions):
        self.vaspFile.write("1.0\n")
        for d in self.multiplyers*directions:
            self.vaspFile.write("  %.10f  %.10f  %.10f\n"%(*d,))

    def write_elements(self,atomNames):
        for atomName in atomNames:
            self.vaspFile.write("%s "%atomName)
        self.vaspFile.write("\n")

    def write_populations(self,cellAtoms):
        for atomNumber in cellAtoms:
            self.vaspFile.write("%d "%(np.prod(self.multiplyers)*atomNumber))
        self.vaspFile.write("\n")

    def write_atoms(self,atomNames):
        for atomName,atom in product(atomNames,
                       self.settings['crystal']):
            if atom[0] != atomName:
                continue
            self.vaspFile.write("  %.10f  %.10f  %.10f  %s\n"%(*atom[1],atom[0]))
        self.vaspFile.write("\n")

    def __del__(self):
        self.vaspFile.close()

#
#░▀█▀░█▀█░█▀▀░█▀█░█▀▄░░░█▀▀░█▀█░█░█░█▀▀░█▀▄
#░░█░░█░█░█░░░█▀█░█▀▄░░░▀▀█░█▀█░▀▄▀░█▀▀░█▀▄
#░▀▀▀░▀░▀░▀▀▀░▀░▀░▀░▀░░░▀▀▀░▀░▀░░▀░░▀▀▀░▀░▀
#

class INCARsaver:
    def __init__(self,oldINCAR,crystal):
        self.oldINCAR = oldINCAR
        self.crystal  = crystal

    @staticmethod
    def mkdir(fileName,flipName):
        try:
            mkdir(fileName+"/"+flipName)
        except OSError as err:
            if err.errno != errno.EEXIST:
                print("Creation of the directory %s/%s failed - does it exist?"%(fileName,flipName))
                exit(error.systemerror)

    @staticmethod
    def copy_POSCAR(fileName,flipName):
        try:
            shutil.copy2("%s/POSCAR"%fileName , "%s/%s/POSCAR"%(fileName,flipName))
        except OSError:
            print("Copying POSCAR to %s didn't work out!"%flipName)
            exit(error.systemerror)

    def write_INCAR(self,flipName,flip): 
        with open(self.fileName+"/"+flipName+"/INCAR","w+") as vaspFile:
            vaspFile.write(re.sub('\s*MAGMOM.*\n','\n',self.oldINCAR))
            vaspFile.write("MAGMOM = ")
            for bit,atom in zip(flip,self.crystal):
                if bit:
                    vaspFile.write("%f "%-atom[2])
                else:
                    vaspFile.write("%f "%atom[2])
            vaspFile.write("\n")

    def save(self,fileName,flips):
        self.fileName = fileName
        INCARsaver.mkdir(fileName,"noFlip")
        INCARsaver.copy_POSCAR(fileName,"noFlip")
        self.write_INCAR("noFlip",np.ones(len(self.crystal)))
        for i,flip in enumerate(flips):
            self.fileName = fileName
            INCARsaver.mkdir(fileName,"flip%05d"%i)
            INCARsaver.copy_POSCAR(fileName,"flip%05d"%i)
            self.write_INCAR("flip%05d"%i,flip)
