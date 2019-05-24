# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import re
import numpy as np
from os import system
from aux.periodic import *

class error:
    unexcepted     = -1000
    nonconvertable = -1001
    vaspError      = -9999


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
            except:
                print("Unexcepted error!")
                exit(error.unexcepted)
    
    @staticmethod
    def find_comment(text):
        return text[0]

    @staticmethod
    def find_directions(text):
        try:
            scale = np.float64(text[1])
        except:
            scale = 1.0
        directions = []    
        for i in range(2,5):
            try:
                directions.append(scale*np.fromstring(text[i],sep=" ")) # crystal directions
            except:
                print("Can't convert crystal directions in \"%s\""%text[i])
                exit(error.unconvertable)
            if(len(directions[-1]) != 3):   
                print("Crystal directions in %s has %d != 3 dimensions!"%(text[i],len(directions[-1])))
                exit(error.unconvertable)
        return directions

    @staticmethod
    def find_cell(text):
        ISDIRECT    = False
        ISSELECTIVE = False

        if re.match("\d",text[5]):
            text.insert(5,'')
        atomNames = text[5].split(" ")
        cellAtoms = np.fromstring(text[6],sep=" ",dtype=np.int)
        cellSize = np.sum(cellAtoms)
        if text[5]=='':
            atomNames = [str(i) for i in range(len(cellAtoms))]

        return atomNames,cellAtoms,cellSize

        cellInputType = text[7][0]
        if cellInputType in "Dd":
            ISDIRECT = True
        elif cellInputType in "Ss":
            ISSELECTIVE = True

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
        if found:
            constrains = [False, False, False]
            for i in range(3):
                if 'T' in found.group(i+1) or \
                   't' in found.group(i+1):
                    constrains[i] = True
            return constrains
        return None


    def __call__(self,i=0):
        return self.rawTxt[i]



#    def find_cellSymmetry(self,text):
#    def find_cellCenter(self,text):
#    def find_cellAtomsCopy(self,text):
#    def find_atomNames(self,text):






def load_POSCAR(inputName,direct=False):
    """
        Reading POSCAR file. Extensive testing required.
                                                        """
    data = {}

    # Returned data:
    comment       = ""  # first line of POSCAR file
    directions    = []  # crystal directions
    cell          = []  # cell read from POSCAR file
    cellSymmetry  = ([],[],[])   # input for spglib symmetry refiner
    cellVolume    = 0.0          # volume of cell
    cellCenter    = np.zeros(3)  # center of volume
    cellAtomsCopy = np.array([],dtype=np.int) # number of atoms in cell
    atomNames     = []                        # name of atoms in cell

    """ Templates:
        according to vasp.wiki:
          direct coords are for 6th lines starting with 'D' or 'd'
      carthesian coords are for 6th lines starting with 'C', 'c', 'K' or 'k' """
    directTemplate     = "Dd"
    carthTemplate      = "CcKk"
    selectiveTemplate  = "Ss"

    # additional variables
    cellAtoms     = np.array([],dtype=np.int)
    cellSize = 1
    cellInputType = 'x'
    ISDIRECT = -1 
    invDirections = []  # inverted crystal directions
    offset = 0
    with open(inputName,"r+") as inFile:
        for i,raw in enumerate(inFile.readlines()):
            line = re.sub("^\s*","",raw)  # remove all blank characters from begining of the line
            line = re.sub("\s+"," ",line) # replace all blank characters in a row to a single space
            if i == 0:
                comment = line[:-1]
            elif i == 1:
                scale = np.float64(line[:-1]) # scaling factor
                if(scale < 0.0):
                    scale = 1.0
            elif i in range(2,5):
                try:
                    directions.append(scale*np.fromstring(line,sep=" ")) # crystal directions
                except:
                    print("Error reading file %s in line %d:\nCan't convert crystal directions."%(inputName,i))
                    exit(-2)
                if(len(directions[-1]) != 3):   
                    print("Error reading file %s in line %d:\nCrystal directions has %d != 3 dimensions!"%(inputName,i,len(directions[-1])))
                    exit(-3)
                cellSymmetry[0].append(tuple(directions[-1]))
                cellCenter += 0.5*directions[-1]
            elif i == 5:
                cellVolume = np.abs(np.linalg.det(np.array(directions)))
                try:
                    invDirections = np.linalg.inv(directions)
                except:
                    print("Error reading file %s in line %d:\nCrystal directions are not basis in 3D!"%(inputName,i))
                    exit(-4)
                if re.match("\D",line):
                    offset += 1
                    atomNames=line[:-1].split(" ")
                else:    
                    cellAtoms = np.fromstring(line,sep=" ",dtype=np.int)
                    cellAtomsCopy = np.copy(cellAtoms)
                    cellSize = np.sum(cellAtoms)
            if i == 6 and offset == 1:
                cellAtoms = np.fromstring(line,sep=" ",dtype=np.int)
                cellAtomsCopy = np.copy(cellAtoms)
                cellSize = np.sum(cellAtoms)
            elif i == 6 + offset:
                cellInputType = line[0]
                if cellInputType in carthTemplate:
                    ISDIRECT = False
                elif cellInputType in directTemplate:
                    ISDIRECT = True
                elif cellInputType in selectiveTemplate:
                    offset += 1
                else:
                    print("Error in POSCAR: unknown input in line %d: %s"%(i,line[:-1]))
                    exit(-1)
            elif ISDIRECT == -1 and i == 7:
                cellInputType = line[0]
                if cellInputType in carthTemplate:
                    ISDIRECT = False
                elif cellInputType in directTemplate:
                    ISDIRECT = True
                else:
                    print("Error in POSCAR: unknown input in line %d: %s"%(i,line[:-1]))
                    exit(-1)
            elif i in range(7+offset,7+offset+cellSize):
                found = re.search("([\-\+]?\d+\.?\d*)\s([\-\+]?\d+\.?\d*)\s([\-\+]?\d+\.?\d*)",line)
                if found:
                    atomType = np.flatnonzero(cellAtoms)
                    cellAtoms[atomType[0]] -= 1
                    atom  = [atomType[0],np.zeros(3)]
                    coords = np.fromstring(found.group(0),sep=" ")
                    if ISDIRECT:
                        cellSymmetry[1].append(tuple(coords))
                        if direct:
                            atom[1] = coords
                        else:
                            for coord,vector in zip(coords,directions):
                                atom[1] += coord*vector 
                    else:
                        cellSymmetry[1].append(tuple(np.dot(coords,invDirections)))
                        if direct:
                            atom[1] = np.dot(coords,invDirections)
                        else:
                            atom[1] = coords
                    cell.append(atom)    
                    cellSymmetry[2].append(atomType[0])
                found = re.search("[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s([a-zA-Z]+)",line)
                if found:
                    if found.group(1) not in atomNames:
                        atomNames.append(found.group(1))
                    
    for i,e in enumerate(cellSymmetry[2]):
        cellSymmetry[2][i] = periodicTableNumber[atomNames[e]]

    data['comment']       = comment
    data['directions']    = directions
    data['cell']          = cell
    data['cellSymmetry']  = cellSymmetry
    data['cellVolume']    = cellVolume
    data['cellCenter']    = cellCenter
    data['cellAtoms']     = cellAtomsCopy
    data['atomNames']     = atomNames
    return data 

#
#
#
#
#

import numpy as np
def load_INCAR(cell,INCARname="INCAR",atomNames=periodicTableElement):
    oldMoments = []
    with open(INCARname,"r") as INCARfile:
        incarData = INCARfile.read()
     
        oldMomentsText = re.search("\s*MAGMOM\s*=\s*(.*)\n",incarData)

        if oldMomentsText is None:
            for atom in cell:
                oldMoments.append(elementMagneticMoment[atomNames[atom[0]]])
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
    """
        Saving data to POSCAR file
{'comment': 'Re/117 - (A3) - HCP [A2] A3 (117) (htqc', 'directions': [array([ 1.38476066, -2.39847581, -0.        ]), array([1.38476066, 2.39847581, 0.        ]), array([0.        , 0.        , 4.47184109])], 'cell': [[0, array([0., 0., 0.])], [0, array([ 1.38476066, -0.79949194,  2.23592055])]], 'cellSymmetry': ([(1.3847606573816698, -2.398475814907534, -0.0), (1.3847606573816698, 2.398475814907534, 0.0), (0.0, 0.0, 4.471841090100609)], [(0.0, -0.0, -0.0), (0.6666666666666714, 0.3333333333333428, 0.5)], [75, 75]), 'cellVolume': 29.704785298855388, 'cellCenter': array([1.38476066, 0.        , 2.23592055]), 'cellAtoms': array([2]), 'atomNames': ['Re']}
                                    """
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
    system("mkdir -p %s/noFlip"%fileName)
    system("cp %s/POSCAR %s/noFlip/POSCAR"%(fileName,fileName))

    with open(fileName+"/noFlip/INCAR","w+") as vaspFile:
        vaspFile.write(re.sub('\s*MAGMOM.*\n','\n',oldINCAR))
        vaspFile.write("MAGMOM = ")
        for i,atom in enumerate(crystal):
            vaspFile.write("%f "%atom[2])
        vaspFile.write("\n")
    for i,flip in enumerate(flips):
        system("mkdir -p %s/flip%d"%(fileName,i))
        system("cp %s/POSCAR %s/flip%d/POSCAR"%(fileName,fileName,i))
        with open(fileName+"/flip"+str(i)+"/INCAR","w+") as vaspFile:
            vaspFile.write(re.sub('\s*MAGMOM.*\n','\n',oldINCAR))
            vaspFile.write("MAGMOM = ")
            for bit,atom in zip(flip,crystal):
                if bit:
                    vaspFile.write("%f "%-atom[2])
                else:
                    vaspFile.write("%f "%atom[2])
            vaspFile.write("\n")
