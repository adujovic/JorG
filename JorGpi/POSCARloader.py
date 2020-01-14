# -*- coding: utf-8 -*-
import re
import numpy as np
import JorGpi.aux.PeriodicTable as periodic

class Error:
    unexcepted     = 10
    nonconvertable = 11

class POSCARloader:
    def __init__(self,*args):
        self.data = []
        self.rawTxt = []
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
            except OSError:
                print("File \"%s\" not found!"%inputName)
            except Exception:
                print("Unexcepted error!")
                exit(Error.unexcepted)

    def __len__(self):
        return len(self.data)

    class FixNames(dict):
        def __init__(self):
            self.update(periodic.periodicTableElement)
        def __missing__(self,key):
            return key
    fix_names = FixNames()

    class FixAtomicNumbers(dict):
        def __init__(self):
            self.update(periodic.periodicTableNumber)
        def __missing__(self,key):
            return key
    fix_atomic_numbers = FixAtomicNumbers()

    @staticmethod
    def read_comment(text):
        return text[0]

    @staticmethod
    def read_directions(text):
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
                exit(Error.unconvertable)
        return directions

    @staticmethod
    def parse_atom(text):
        found = re.search("([\-\+]?\d+\.?\d*)\s([\-\+]?\d+\.?\d*)\s([\-\+]?\d+\.?\d*)",text)
        if found:
            return np.fromstring(found.group(0),sep=" ")
        return None

    @staticmethod
    def parse_atom_name(text):
        found = re.search(
"[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s([a-zA-Z]+).*",text)
        if found:
            return found.group(1)
        return None

    @staticmethod
    def parse_constrains(text):
        found = re.search(
"[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s([tTfF]).*\s([tTfF]).*\s([tTfF]).*",text)
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


    def parse(self):
        self.data = []
        for text in self.rawTxt:
            read = CellReader(text)
            self.data.append(read.read())
class CellReader:
    def __init__(self,text):
        self.text          = text

        self.read_atoms()
        self.directions    = POSCARloader.read_directions(text)
        self.volume        = np.linalg.det(self.directions)
        self.center        = np.zeros(3)
        self.cell          = []
        self.isdirect      = False
        self.isselective   = False
        self.cellSymmetry  = ([tuple(d) for d in self.directions],
                              [],[]) # directions, direct units cell, atomic numbers
        self.check_type(text[7][0])

    def check_type(self,character):
        if character in "DdSs":
            self.isdirect    = True
        if character in "Ss":
            self.isselective = True

    def read_single_atom(self,line):
        if not self.isselective    \
           and self.atomRead == 0  \
           and self.atomNames[self.atomType] == self.atomType:
            self.atomNames[self.atomType] = POSCARloader.parse_atom_name(line)

        atomCoordinates = POSCARloader.parse_atom(line)
        if self.isdirect:
            self.cell.append((self.atomNames[self.atomType],
                         np.dot(self.directions,atomCoordinates)))
            self.center += np.dot(self.directions,atomCoordinates)
            self.cellSymmetry[1].append(tuple(atomCoordinates))
        else:
            self.cell.append((self.atomNames[self.atomType],atomCoordinates))
            self.cellSymmetry[1].append(tuple(
                    np.dot(np.linalg.inv(self.directions),atomCoordinates)))
            self.center += atomCoordinates

        self.cellSymmetry[2].append(POSCARloader.fix_atomic_numbers[self.atomNames[self.atomType]])

    def read_atoms(self):
        if re.match("\d",self.text[5]):
            self.text.insert(5,'')
        self.atomNames = self.text[5].split()
        self.cellAtoms = np.fromstring(self.text[6],sep=" ",dtype=np.int)
        self.cellSize  = np.sum(self.cellAtoms)
        if not self.atomNames:
            self.atomNames = [i for i in range(len(self.cellAtoms))]

    def read(self):
        self.atomRead = 0
        self.atomType = 0
        for i in range(8,8+self.cellSize):
            self.read_single_atom(self.text[i])
            self.atomRead += 1
            if self.atomRead == self.cellAtoms[self.atomType]:
                self.atomType += 1
                self.atomRead  = 0
        self.center /= self.cellSize
        self.atomNames = [POSCARloader.fix_names[atom] for atom in self.atomNames]
        return {'comment'     : POSCARloader.read_comment(self.text),
                'directions'  : self.directions,
                'cell'        : self.cell,
                'cellSymmetry': self.cellSymmetry,
                'cellVolume'  : self.volume,
                'cellCenter'  : self.center,
                'cellAtoms'   : self.cellAtoms,
                'atomNames'   : self.atomNames}
