# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import re
import numpy as np
from os import makedirs,mkdir,rmdir
import errno
import shutil
import aux.periodic as periodic

class error:
    systemerror = -1
    unexcepted     = 10
    nonconvertable = 11
    vaspError      = 99


def load_INCAR(cell,INCARname="INCAR",atomNames=periodic.periodicTableElement):
    oldMoments = []
    with open(INCARname,"r") as INCARfile:
        incarData = INCARfile.read()

        oldMomentsText = re.search("\s*MAGMOM\s*=\s*(.*)\n",incarData)

        if oldMomentsText is None:
            for atom in cell:
                oldMoments.append(periodic.elementMagneticMoment[atom[0]])
        else:
            magmomLine = oldMomentsText.group(1)
            for record in re.findall("\s*([0-9]+)\s*\*\s*([\-0-9\.]+)",magmomLine):
                magmomLine = re.sub("{:s}\s*\*\s*{:s}".format(record[0],record[1]),(record[1]+" ")*int(record[0]),magmomLine)
            for moment in magmomLine.split():
                oldMoments.append(np.float(moment))

    return oldMoments,incarData
#
#
#
#
#

import numpy as np
def save_xyz(fileName,crystal,numberOfAtoms = -1, selectedAtoms = None):
    """
        Saving data to xyz-style file
                                        """
    if numberOfAtoms < 0:
        numberOfAtoms = len(crystal)


    with open(fileName,"w+") as xyzFile:
        xyzFile.write(str(numberOfAtoms))
        xyzFile.write("\n\n")
        for i,atom in enumerate(crystal):
            xyzFile.write("%s"%atom[0])
            for xyz in atom[1]:
                xyzFile.write(" %.10f"%xyz)
            if selectedAtoms is not None:
                vector = 2*np.random.ranf(3)-1.0
                vector = 0.2*vector/np.linalg.norm(vector)
                if i == selectedAtoms[0]:
                    vector = 2.0*vector
                    xyzFile.write(" PatrialCharge(1.0) %f %f %f"%tuple(vector))
                elif i in selectedAtoms[1:]:
                    xyzFile.write(" PatrialCharge(1.0) %f %f %f"%tuple(vector))
            xyzFile.write("\n")
        xyzFile.write("\n")

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
#
def save_POSCAR(fileName,crystal,multiplyers,data):
    """
        Saving data to POSCAR file
                                    """
    with open(fileName,"w+") as vaspFile:
        vaspFile.write(data['comment'])
        vaspFile.write("\n1.0\n")
        for m,d in zip(multiplyers,data['directions']):
            for field in d:
                vaspFile.write("  %.10f"%((m+1)*field))
            vaspFile.write("\n")
        for atomName in data['atomNames']:
            vaspFile.write("%s "%atomName)
        vaspFile.write("\n")
        for atomNumber in data['cellAtoms']:
            vaspFile.write("%d "%(np.prod(np.array(multiplyers)+1)*atomNumber))
        vaspFile.write("\nCarthesian\n")
        for atomName in data['atomNames']:
            for atom in crystal:
                if atom[0]==atomName:
                    for vasp in atom[1]:
                        vaspFile.write(" %.10f "%vasp)
                    vaspFile.write(" %s\n"%atom[0])
        vaspFile.write("\n")
#
#
#
#
#
import re
def save_INCAR(fileName,oldINCAR,crystal,flips):
    """
        Saving data to POSCAR file
                                    """
    try:
        mkdir("%s/noFlip"%fileName)
    except OSError as err:
        if err.errno != errno.EEXIST:
            print("Creation of the directory %s/noFlip failed - does it exist?"%fileName)
            exit(error.systemerror)

    try:
        shutil.copy2("%s/POSCAR"%fileName , "%s/noFlip/POSCAR"%fileName)
    except OSError:
        print("Copying POSCAR to noFlip didn't work out!")
        exit(error.systemerror)

    with open(fileName+"/noFlip/INCAR","w+") as vaspFile:
        vaspFile.write(re.sub('\s*MAGMOM.*\n','\n',oldINCAR))
        vaspFile.write("MAGMOM = ")
        for i,atom in enumerate(crystal):
            vaspFile.write("%f "%atom[2])
        vaspFile.write("\n")
    for i,flip in enumerate(flips):
        try:
            mkdir("%s/flip%d"%(fileName,i))
        except OSError as err:
            if err.errno != errno.EEXIST:
                print("Creation of the directory %s/flip%d failed - does it exist?"%(fileName,i))
            exit(error.systemerror)
        try:
            shutil.copy2("%s/POSCAR"%fileName , "%s/flip%d/POSCAR"%(fileName,i))
        except OSError:
            print("Copying POSCAR to flip%i didn't work out!"%i)
            exit(error.systemerror)
        with open(fileName+"/flip"+str(i)+"/INCAR","w+") as vaspFile:
            vaspFile.write(re.sub('\s*MAGMOM.*\n','\n',oldINCAR))
            vaspFile.write("MAGMOM = ")
            for bit,atom in zip(flip,crystal):
                if bit:
                    vaspFile.write("%f "%-atom[2])
                else:
                    vaspFile.write("%f "%atom[2])
            vaspFile.write("\n")
# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import re
import numpy as np
from os import system
import aux.periodic as periodic
import copy


class POSCARloader:
    data = []
    rawTxt = []

    def __init__(self,*args,**kwargs):
        for inputName in args:
            try:
                with open(inputName,"r+") as inFile:
                    # Remove white characters from the beggining
                    clearTxt = [re.sub("^\s+","",line) for line in inFile.readlines()]
                    # Remove white characters from the end
                    clearTxt = [re.sub("\s+$","",line) for line in clearTxt]
                    #remove double white characters
                    clearTxt = [re.sub("\s+"," ",line) for line in clearTxt]
                    self.rawTxt.append(clearTxt)
            except FileNotFoundError:
                print("File \"%s\" not found!"%inputName)
            except Exception:
                print("Unexcepted error!")
                exit(error.unexcepted)

    def __len__(self):
        return len(self.data)

    class fix_names(dict):
        def __init__(self):
            self.update(periodic.periodicTableElement)
        def __missing__(self,key):
            return key
    fix_names = fix_names()

    class fix_atomic_numbers(dict):
        def __init__(self):
            self.update(periodic.periodicTableNumber)
        def __missing__(self,key):
            return key
    fix_atomic_numbers = fix_atomic_numbers()

    @staticmethod
    def find_comment(text):
        return text[0]

    @staticmethod
    def find_directions(text):
        try:
            scale = np.float64(text[1])
        except TypeError:
            scale = 1.0
        directions = []
        for i in range(2,5):
            try:
                directions.append(scale*np.fromstring(text[i],sep=" ")) # crystal directions
            except ValueError: 
                print("Can't convert crystal directions in \"%s\""%text[i])
                exit(error.unconvertable)
            if(len(directions[-1]) != 3):
                print("Crystal directions in %s has %d != 3 dimensions!"%(text[i],len(directions[-1])))
                exit(error.unconvertable)
        return directions

    def check_type(self,character):
        self.ISDIRECT    = False
        self.ISSELECTIVE = False
        if character in "DdSs":
            self.ISDIRECT   = True
        if character in "Ss":
            self.ISSELECTIVE = True

    def find_single_atom(self,line,atomRead,atomType):
        if not self.ISSELECTIVE    \
           and atomRead == 0  \
           and self.atomNames[atomType] == atomType:
            self.atomNames[atomType] = POSCARloader.parse_atomName(line)

        atomCoordinates = POSCARloader.parse_atom(line)
        if self.ISDIRECT:
            self.cell.append((self.atomNames[atomType],
                         np.dot(self.directions,atomCoordinates)))
            self.center += np.dot(self.directions,atomCoordinates)
            self.cellSymmetry[1].append(tuple(atomCoordinates))
        else:
            self.cell.append((self.atomNames[atomType],atomCoordinates))
            self.cellSymmetry[1].append(tuple(
                    np.dot(np.linalg.inv(self.directions),atomCoordinates)))
            self.center += atomCoordinates

        self.cellSymmetry[2].append(self.fix_atomic_numbers[self.atomNames[atomType]])

    def find_atom_names(self,line):
        self.atomNames = line.split()

    def find_cell_size(self,line):
        self.cellAtoms = np.fromstring(line,sep=" ",dtype=np.int)
        self.cellSize = np.sum(self.cellAtoms)

    def find_cell(self,text):
        if re.match("\d",text[5]):
            text.insert(5,'')
        self.find_atom_names(text[5])
        self.find_cell_size (text[6])
        if not self.atomNames:
            self.atomNames = [i for i in range(len(self.cellAtoms))]

        self.directions = POSCARloader.find_directions(text)
        self.volume     = np.linalg.det(self.directions)

        self.center = np.zeros(3)
        self.cell = []
        self.cellSymmetry = ([tuple(d) for d in self.directions],
                             [],[]) # directions, direct units cell, atomic numbers

        cellInputType = text[7][0]
        self.check_type(cellInputType)

        atomRead = 0
        atomType = 0
        for i in range(8,8+self.cellSize):
            self.find_single_atom(text[i],atomRead,atomType)
            atomRead += 1
            if atomRead == self.cellAtoms[atomType]:
                atomType += 1
                atomRead = 0

        self.center /= self.cellSize
        self.atomNames = [POSCARloader.fix_names[atom] for atom in self.atomNames]

    @staticmethod
    def parse_atom(text):
        found = re.search("([\-\+]?\d+\.?\d*)\s([\-\+]?\d+\.?\d*)\s([\-\+]?\d+\.?\d*)",text)
        if found:
            return np.fromstring(found.group(0),sep=" ")
        return None

    @staticmethod
    def parse_atomName(text):
        found = re.search("[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s([a-zA-Z]+).*",text)
        if found:
            return(found.group(1))
        return None

    @staticmethod
    def parse_constrains(text):
        found = re.search("[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s([tTfF]).*\s([tTfF]).*\s([tTfF]).*",text)
        if not found:
            return None
        constrains = [False, False, False]
        for i in range(3):
            if 'T' in found.group(i+1) or \
               't' in found.group(i+1):
                constrains[i] = True
        return constrains

    def __call__(self,i=0):
        try:
            return self.data[i]
        except IndexError:
            print("Run parse first!")


    def parse_file(self,text):
        self.comment = POSCARloader.find_comment(text)
        self.find_cell(text)

    def parse(self):
        self.data = []
        for text in self.rawTxt:
            self.parse_file(text)
            self.data.append({})
            self.data[-1]['comment']       = self.comment
            self.data[-1]['directions']    = self.directions
            self.data[-1]['cell']          = self.cell
            self.data[-1]['cellSymmetry']  = self.cellSymmetry
            self.data[-1]['cellVolume']    = self.volume
            self.data[-1]['cellCenter']    = self.center
            self.data[-1]['cellAtoms']     = self.cellAtoms
            self.data[-1]['atomNames']     = self.atomNames
