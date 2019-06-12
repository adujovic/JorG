# -*- coding: utf-8 -*-

import re
from sys import argv,path
path.insert(0,r'../../')
from JorG.loadsave import POSCARloader
import numpy as np
from itertools import product

class error:
    unexcepted = 12

class EnergyConverter:
    energyRatios = {'eV' :    1.0,
                    'meV':    1000.0,
                    'Ry' :    np.reciprocal(13.6056980659),
                    'mRy':    1000.0*np.reciprocal(13.6056980659),
                    'He' :    0.5*np.reciprocal(13.6056980659),
                    'mHe':    500.0*np.reciprocal(13.6056980659),
                    'K'  :    11604.51812}
    @staticmethod
    def convert(*args,**kwargs):
        try:
            return [EnergyConverter.energyRatios[kwargs['units']]*arg for arg in args]
        except KeyError:
            print("No unit defined! Values will be in eV")
            return args

class ReadMoments:
    @staticmethod
    def should_skip(line):
        for regex in ['^\s+$','^\s*#','^-+$',
                      '^\s*tot','magnetization \(x\)']:
            if re.search(regex,line):
                return True
        return False

    def should_stop__reading(self,line):
        if re.search('^\s*tot',line):
            self.ISTOBEREAD = False

    def should_start_reading(self,line):
        if re.search(' \(x\)',line):
            self.ISTOBEREAD = True

    @staticmethod
    def read_energy(line):
        if "energy  without entropy" in line:
            return float(line.split()[-1])

    def __init__(self):
        self.moments    = {}
        self.energy     = None
        self.ISTOBEREAD = False

    def __call__(self,text):
        for line in text:
            self.read_line(line)
        return {'moments': self.moments, 'energy': self.energy}

    def read_line(self,line):
        self.should_start_reading(line)
        self.should_stop__reading(line)
        if self.should_skip(line):
            return
        elif self.ISTOBEREAD:
            self.moments[int(line.split()[0])] = float(line.split()[4])
        elif self.energy is None:
            self.energy = self.read_energy(line)

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

    def parse_text(self,text):
        read_moments = ReadMoments()
        self.data.append(read_moments(text))

    def parse(self):
        self.data   = []
        for text in self.rawTxt:
            self.parse_text(text)

    def get_energy(self,idx=0):
        return self.data[idx]['energy']

    def get_moments(self,idx=0):
        return self.data[idx]['moments']

    def __call__(self,idx=0):
        return self.data[idx]

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
                 namesOfInteractingAtoms):
        self.numberOfNeighbors       = numberOfNeighbors
        self.namesOfInteractingAtoms = namesOfInteractingAtoms
        self.solution = None
        self.solver   = None

    def read_POSCARs(self,*args):
        self.solution = None
        POSCARs = [ "%s/POSCAR"%arg for arg in args ]
        self.POSCARs = POSCARloader(*POSCARs,spam=False)
        self.POSCARs.parse()

    def read_MAGMOMs(self,*args):
        self.solution = None
        OUTCARs = [ "%s/OUTCAR"%arg for arg in args ]
        self.MAGMOMs = MAGMOMloader(*OUTCARs,spam=False)
        self.MAGMOMs.parse()

    def read(self,*args):
        self.solution = None
        self.read_POSCARs(*args)
        self.read_MAGMOMs(*args)

    def make_crystal(self,idx=0):
        self.solution = None
        self.crystal = []
        for mul,(i,atom) in product(product(range(-1,1),repeat=3),
                                    enumerate(self.POSCARs(idx)['cell'])):
            boost = np.dot(mul,self.POSCARs(idx)['directions'])
            self.crystal.append([i,atom[0],atom[1]+boost])

    def map_distances(self,idx=0):
        self.solution = None
        self.distances = set([])
        self.make_crystal(idx)
        for atom in self.POSCARs(idx)['cell']:
            d = np.around(np.linalg.norm(atom[1]),decimals=2)
            if atom[0] in self.namesOfInteractingAtoms:
                self.distances.add(d)

        self.distances = np.sort(np.array([d for d in self.distances]))[1:1+self.numberOfNeighbors]

    # for sorted!
    def get_system_of_equations(self):
        self.solution = None
        self.map_distances()
        self.systemOfEquations = []
        self.dE = []
        for i in range(1,len(self.MAGMOMs)):
            try:
                self.dE.append(self.MAGMOMs(i)['energy']-self.MAGMOMs(0)['energy'])
            except TypeError:
                print("VASP hasn't finished this run (%d/%d)"%(i,len(self.MAGMOMs)-1))
                continue
            self.get_SOE_line(i)
        self.solver = EquationSolver(self.systemOfEquations,self.dE)

    def get_SOE_line(self,i):
        self.systemOfEquations.append(np.zeros(self.numberOfNeighbors))
        self.flipped = []
        for idx,atom in enumerate(self.POSCARs(0)['cell']):
            self.get_flips(i,idx,atom)
        self.get_SOE_elements(i)

    def get_flips(self,i,idx,atom):
        if atom[0] not in self.namesOfInteractingAtoms:
           return
        momentA = self.MAGMOMs(0)['moments'][idx+1]
        momentB = self.MAGMOMs(i)['moments'][idx+1]
        scalar = momentA * momentB
        if (abs(scalar) > 1e-5 and scalar < 0.0):
            self.flipped.append(idx)
      
    def get_SOE_elements(self,idx):
        for flip,atom in product(self.flipped,self.crystal):
            if atom[1] not in self.namesOfInteractingAtoms:
                continue
            elif atom[0] in self.flipped:
                continue
            self.get_SOE_element(idx,flip,atom)

    def get_SOE_element(self,idx,flip,atom):
        d = np.around(np.linalg.norm(atom[2]-self.POSCARs(0)['cell'][flip][1]),decimals=2)
        if d in self.distances:
            self.systemOfEquations[idx-1][np.argwhere(self.distances==d)] -= 2*self.MAGMOMs(0)['moments'][atom[0]+1]*self.MAGMOMs(idx)['moments'][flip+1]


    def solve(self,**kwargs):
        if self.solver is None:
            self.get_system_of_equations()
        return EnergyConverter.convert(*(self.solver.solve()),**kwargs)

class EquationSolver:
    def __init__(self,equations,vector):
        self.equations = np.array(equations)
        self.vector    = np.array(vector)
        self.solution  = None

    def solve(self,**kwargs):
        if self.solution is not None:
            return self.solution
        self.solve_system_of_equations(**kwargs)
        return self.solution

    def solve_system_of_equations(self, **kwargs):
        if len(self.equations) == 1:
            self.solution = self.vector[0]/self.equations[0]
        elif self.equations.shape[0] == self.equations.shape[1]:
            self.solution = np.linalg.solve(self.equations,self.vector,**kwargs)
        else:
            self.solution = np.linalg.lstsq(self.equations,self.vector,rcond=None,**kwargs)[0]


