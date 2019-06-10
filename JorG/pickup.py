# -*- coding: utf-8 -*-

import re
from sys import argv,path
path.insert(0,r'../../')
from JorG.loadsave import POSCARloader
import numpy as np
from itertools import product

class error:
    unexcepted = 12

BRUTAL='Fe'
CUTOFF=2.46
NUMBER=1

class MAGMOMloader:
    rawTxt = []
    data   = []
    def __init__(self,*args,**kwargs):
        for inputName in args:
            try:
                with open(inputName,"r+") as inFile:
                    self.rawTxt.append(inFile.readlines())
            except FileNotFoundError:
                print("File \"%s\" not found!"%inputName)
            except Exception:
                print("Unexcepted error!")
                exit(error.unexcepted)

    def __len__(self):
        return len(self.data)

    @staticmethod
    def should_skip(line):
        for regex in ['^\s+$','^\s*#','^-+$',
                      '^\s*tot','magnetization \(x\)']:
            if re.search(regex,line):
                return True
        return False

    @staticmethod
    def should_stop__reading(line,ISTOBEREAD):
        if re.search('^\s*tot',line):
            return False
        return ISTOBEREAD

    @staticmethod
    def should_start_reading(line,ISTOBEREAD):
        if re.search(' \(x\)',line):
            return True
        return ISTOBEREAD

    @staticmethod
    def read_energy(line):
        if "energy  without entropy" in line:
            return float(line.split()[-1])

    @staticmethod
    def read_moments(text):
        moments={}
        energy = None
        ISTOBEREAD=False
        for line in text:
            ISTOBEREAD = MAGMOMloader.should_start_reading(line,ISTOBEREAD)
            ISTOBEREAD = MAGMOMloader.should_stop__reading(line,ISTOBEREAD)
            if MAGMOMloader.should_skip(line):
                continue
            elif ISTOBEREAD:
                moments[int(line.split()[0])] = float(line.split()[4])
            elif energy is None:
                energy = MAGMOMloader.read_energy(line)
        return {'moments': moments, 'energy': energy}

    def parse_text(self,text):
        self.data.append(MAGMOMloader.read_moments(text))

    def parse(self):
        self.data   = []
        for text in self.rawTxt:
            self.parse_text(text)
            print(self.data[-1])

#
#░█▀▀░█▄█░█▀█░█▀▄░▀█▀░░░█▀█░▀█▀░█▀▀░█░█░░░░░█░█░█▀█
#░▀▀█░█░█░█▀█░█▀▄░░█░░░░█▀▀░░█░░█░░░█▀▄░▄▄▄░█░█░█▀▀
#░▀▀▀░▀░▀░▀░▀░▀░▀░░▀░░░░▀░░░▀▀▀░▀▀▀░▀░▀░░░░░▀▀▀░▀░░
#
class SmartPickUp:
    POSCARs           = None
    MAGMOMs           = None
    systemOfEquations = []
    def __init__(self,   numberOfNeighbors,
                 cutoff, namesOfInteractingAtoms):
        self.numberOfNeighbors       = numberOfNeighbors
        self.namesOfInteractingAtoms = namesOfInteractingAtoms
        self.cutoff                  = cutoff

    def read_POSCARs(self,*args):
        POSCARs = [ "%s/POSCAR"%arg for arg in args ]
        self.POSCARs = POSCARloader(*POSCARs,spam=False)
        self.POSCARs.parse()

    def read_MAGMOMs(self,*args):
        OUTCARs = [ "%s/OUTCAR"%arg for arg in args ]
        self.MAGMOMs = MAGMOMloader(*OUTCARs,spam=False)
        self.MAGMOMs.parse()

    def read(self,*args):
        self.read_POSCARs(*args)
        self.read_MAGMOMs(*args)

    def make_crystal(self,idx=0):
        self.crystal = []
        for mul,(i,atom) in product(product(range(-1,1),repeat=3),
                                    enumerate(self.POSCARs(idx)['cell'])):
            boost = np.dot(mul,self.POSCARs(idx)['directions'])
            self.crystal.append([i+1,atom[0],atom[1]+boost])

    def map_distances(self,idx=0):
        self.distances = set([])
        self.make_crystal(idx)
        for atom in self.POSCARs(idx)['cell']:
            d = np.around(np.linalg.norm(atom[1]),decimals=2)
            if d<self.cutoff and atom[0] in self.namesOfInteractingAtoms:
                self.distances.add(d)

        self.distances = np.array([d for d in self.distances])
        self.distances = np.sort(self.distances)
        self.distances = self.distances[1:1+self.numberOfNeighbors]

    def get_system_of_equations(self):
        self.map_distances()

# return moments,directions,cell,energy
#data = []
#for name in argv[1:]:
#    data.append((read_flip(name)))
#
#recordA = data[0]
#crystal = get_crystal(*recordA,)
#
#distances= set([])
#
#
#
#
#
#
#
#
#systemOfEquations = []
#dE = []
#for i,recordB in enumerate(data[1:]):
#    dE.append(recordA[3]-recordB[3])
#    systemOfEquations.append(np.zeros(NUMBER))
#    flipped = []
#    for a,b in zip(recordA[0],recordB[0]):
#        scalar = recordA[0][a] * recordB[0][b]
#        if (abs(scalar) > 1e-5 and scalar < 0.0
#             and recordA[2][a-1][0] in BRUTAL
#             and recordB[2][b-1][0] in BRUTAL):
#            flipped.append(a)
#
#    for f,(j,atom) in product(flipped,enumerate(crystal)):
#        if atom[0] not in flipped and atom[1] in BRUTAL:
#            d = np.around(np.linalg.norm(atom[2]-recordB[2][f-1][1]),decimals=2)
#            if d in distances:
#                systemOfEquations[i][np.argwhere(distances==d)] += recordA[0][atom[0]]*recordB[0][f]
#
#if len(systemOfEquations) == 1:
#    print(500*dE[0]/systemOfEquations[0]/13.6056980659)
#    exit()
#systemOfEquations = np.array(systemOfEquations)
#dE                = np.array(dE)
#
#if systemOfEquations.shape[0] == systemOfEquations.shape[1]:
#    sol = np.linalg.solve(systemOfEquations,dE)
#else:
#    sol = np.linalg.lstsq(systemOfEquations,dE,rcond=None)[0]
#
#print("J in meV")
#print(*(500*sol))
