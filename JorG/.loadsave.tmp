# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import re
import numpy as np
from os import system
import aux.periodic as periodic
import copy

class error:
    unexcepted     = 10
    nonconvertable = 11
    vaspError      = 99


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

    def find_cell(self,text):
        self.center = np.zeros(3)
        if re.match("\d",text[5]):
            text.insert(5,'')
        self.atomNames = text[5].split(" ")
        cellAtoms = np.fromstring(text[6],sep=" ",dtype=np.int)
        cellSize = np.sum(cellAtoms)
        if text[5]=='':
            self.atomNames = [i for i in range(len(cellAtoms))]

        self.directions = POSCARloader.find_directions(text)
        self.volume     = np.linalg.det(self.directions)

        self.cell = []
        self.cellSymmetry = [[tuple(d) for d in self.directions],
                             [],[]] # directions, direct units cell, atomic numbers

        cellInputType = text[7][0]
        self.check_type(cellInputType)

        atomRead = 0
        atomType = 0
        for i in range(8,8+cellSize):
            self.find_single_atom(text[i],atomRead,atomType)
            atomRead += 1
            if atomRead == cellAtoms[atomType]:
                atomType += 1
                atomRead = 0

        self.center /= cellSize
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
        if found:
            constrains = [False, False, False]
            for i in range(3):
                if 'T' in found.group(i+1) or \
                   't' in found.group(i+1):
                    constrains[i] = True
            return constrains
        return None

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
            self.data[-1]['cellAtoms']     = copy.deepcopy(self.cell)
            self.data[-1]['atomNames']     = self.atomNames

