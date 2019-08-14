# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import re
import numpy as np
from os import mkdir
import errno
import shutil
import aux.PeriodicTable as periodic

class Error:
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
        self.incarFile = open(self.settings['fileName'],"r")
        self.oldMoments = []
        self.incarData = self.incarFile.read()
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
            magmomLine = re.sub("{:s}\s*\*\s*{:s}".format(record[0],record[1]),
                                (record[1]+" ")*int(record[0]),magmomLine)
        for moment in magmomLine.split():
            self.oldMoments.append(np.float(moment))
        return self.oldMoments,self.incarData

    def __del__(self):
        self.incarFile.close()
#
#
#
#
#

from itertools import product
class SaveXYZ:
    settings = {'fileName'     : 'data4jmol.xyz',
                'selectedAtoms': None}

    def __init__(self,crystal,**kwargs):
#    """
#        Saving data to xyz-style file
#                                        """
        self.settings.update(kwargs)

        self.xyzFile = open(self.settings['fileName'],"w+")
        self.xyzFile.write(str(len(crystal)))
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
        for direction in self.multiplyers*directions:
            self.vaspFile.write("  %.10f  %.10f  %.10f\n"%(*direction,))

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
                exit(Error.systemerror)

    @staticmethod
    def copy_POSCAR(filename,flipname):
        try:
            shutil.copy2("%s/POSCAR"%filename , "%s/%s/POSCAR"%(filename,flipname))
        except OSError:
            print("Copying POSCAR to %s didn't work out!"%flipname)
            exit(Error.systemerror)

    def write_incar(self,flipname,flip):
        with open(self.fileName+"/"+flipname+"/INCAR","w+") as vaspFile:
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
        self.write_incar("noFlip",np.zeros(len(self.crystal)))
        for i,flip in enumerate(flips):
            self.fileName = fileName
            INCARsaver.mkdir(fileName,"flip%05d"%i)
            INCARsaver.copy_POSCAR(fileName,"flip%05d"%i)
            self.write_incar("flip%05d"%i,flip)
