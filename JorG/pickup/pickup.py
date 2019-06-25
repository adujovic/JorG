# -*- coding: utf-8 -*-

import re
from sys import argv,path
path.insert(0,r'../../')
from generator import apply_mirrorsXYZ
from POSCARloader import POSCARloader
from heisenberg import EquationSolver,NaiveHeisenberg
import numpy as np
from itertools import product
import argparse as ap
from copy import copy

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
    default = { 'groundMoment'  : 1.0,
                'excitedMoment' : 1.0}
    types = [ "       no moment",
              "  average moment",
              "geometric moment",
              " original moment",
              " neoteric moment" ]   
    @staticmethod
    def convert(*args,**kwargs):
        # Returns array of J values with different conventions (see types)
        settings = EnergyConverter.default
        settings.update(kwargs)
        try:
            notMomentSq = 1.0
            avgMomentSq = 0.25*(settings['groundMoment']+settings['excitedMoment'])**2
            geoMomentSq = settings['groundMoment']*settings['excitedMoment']
            orgMomentSq = settings['groundMoment']*settings['groundMoment']
            newMomentSq = settings['excitedMoment']*settings['excitedMoment'] 
            return [ [notMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args],
                     [avgMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args],
                     [geoMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args],
                     [orgMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args],
                     [newMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args]]

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
        if self.ISTOBEREAD:
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

    def get_average_magnitude(self,idx=0,indices=None):
        if indices is not None:
            return np.average(np.abs(np.take(list(self.data[idx]['moments'].values()),indices)))
        return np.average(np.abs(np.array(list(self.data[idx]['moments'].values()))))

    def __call__(self,idx=0):
        return self.data[idx]

#
#░█▀▀░█▄█░█▀█░█▀▄░▀█▀░░░█▀█░▀█▀░█▀▀░█░█░░░░░█░█░█▀█
#░▀▀█░█░█░█▀█░█▀▄░░█░░░░█▀▀░░█░░█░░░█▀▄░▄▄▄░█░█░█▀▀
#░▀▀▀░▀░▀░▀░▀░▀░▀░░▀░░░░▀░░░▀▀▀░▀▀▀░▀░▀░░░░░▀▀▀░▀░░
#
class SmartPickUp:
    def __init__(self,   numberOfNeighbors,
                 namesOfInteractingAtoms):
        self.numberOfNeighbors       = numberOfNeighbors
        self.namesOfInteractingAtoms = namesOfInteractingAtoms
        self.types                   = EnergyConverter.types

    def read_POSCARs(self,*args):
        lastBackSlashRemoved = [ re.sub('/$','',arg) for arg in args ]
        POSCARs = [ "%s/POSCAR"%arg for arg in lastBackSlashRemoved ]
        self.POSCARs = POSCARloader(*POSCARs,spam=False)
        self.POSCARs.parse()

    def read_MAGMOMs(self,*args):
        OUTCARs = [ "%s/OUTCAR"%arg for arg in args ]
        self.MAGMOMs = MAGMOMloader(*OUTCARs,spam=False)
        self.MAGMOMs.parse()

    def read(self,*args,**kwargs):
        self.read_MAGMOMs(*args)
        self.read_POSCARs(*args)
        if 'reference' in kwargs:
            self.reference = self.POSCARs(0)['cell'][kwargs['reference']][1]
            self.ref       = kwargs['reference']
        else:
            print("Warning: reference @ 0. Is that ok?")
            self.ref       = 0
            self.reference = self.POSCARs(0)['cell'][0][1]

    def make_crystal(self,idx=0):
        self.crystal  = self.POSCARs(idx)['cell']
        self.crystal  = [ [atom[0],atom[1],self.MAGMOMs.get_moments()[i+1]] for i,atom in enumerate(self.crystal) ]
        self.crystal8 = apply_mirrorsXYZ(self.POSCARs(0)['directions'],self.crystal,
                                         reference=self.ref)

    def map_distances(self,idx=0):
        self.distances = set([])
        self.make_crystal(idx)
        for atom in self.POSCARs(idx)['cell']:
            d = np.around(np.linalg.norm(atom[1]-self.reference),decimals=2)
            if atom[0] in self.namesOfInteractingAtoms:
                self.distances.add(d)
        self.distances = np.sort(np.array(list(self.distances)))[1:1+self.numberOfNeighbors]

    # for sorted!
    def get_system_of_equations(self):
        self.map_distances()
        self.systemOfEquations      = []
        self.dE                     = []
        self.flippingConfigurations = []
        for i in range(1,len(self.MAGMOMs)):
            try:
                self.dE.append(self.MAGMOMs(i)['energy']-self.MAGMOMs(0)['energy'])
            except TypeError:
                print("VASP hasn't finished this run (%d/%d)"%(i,len(self.MAGMOMs)-1))
                continue
            self.set_flipps(i)
        self.model             = NaiveHeisenberg(self.flippingConfigurations,self.crystal,self.crystal8)
        self.model.MAGMOMs     = self.MAGMOMs
        self.flipped           = np.unique(np.where(self.flippingConfigurations)[1])
        self.systemOfEquations = self.model.generate(self.namesOfInteractingAtoms,self.distances)
        self.solver            = EquationSolver(self.systemOfEquations,self.dE)
        self.solver.remove_tautologies()

    def set_flipps(self,i):
        self.flippingConfigurations.append([])
        for idx,atom in enumerate(self.POSCARs(0)['cell']):
            self.get_flip(i,idx,atom)

    def get_flip(self,i,idx,atom):
        if atom[0] not in self.namesOfInteractingAtoms:
           self.flippingConfigurations[-1].append(False)
           return
        momentA = self.MAGMOMs(0)['moments'][idx+1]
        momentB = self.MAGMOMs(i)['moments'][idx+1]
        scalar = momentA * momentB
        if (abs(scalar) > 1e-5 and scalar < 0.0):
            self.flippingConfigurations[-1].append(True)
            return
        self.flippingConfigurations[-1].append(False)

    def solve(self,**kwargs):
        try:
            self.solver
        except AttributeError:
            self.get_system_of_equations()
        averageMground = self.MAGMOMs.get_average_magnitude(0,self.flipped)
        averageMexcite = np.average([self.MAGMOMs.get_average_magnitude(i,self.flipped) for i in range(1,len(self.POSCARs))])
        return EnergyConverter.convert(*(self.solver.solve()),groundMoment=averageMground,excitedMoment=averageMexcite,**kwargs)

class Reference:
    def __init__(self,POSCAR):
        loader = POSCARloader(POSCAR)
        loader.parse()
        firstLine = loader()['comment']
        try:
            self.reference = int(re.search('NewRef:\s*([0-9]+),',firstLine).group(1))
        except AttributeError as err:
            print("Reference not found in POSCAR comment! Taking 0!")
            print(err)
            self.reference = 0
        except ValueError as err:
            print("Reference cannot be converted to int! Taking 0!")
            print(err)
            self.reference = 0
        if self.reference < 0:
            print("Wrong reference (%d < 0)! Taking 0!"%self.reference)
            self.reference = 0
    def __str__(self):
        return str(self.reference)
    def __call__(self):
        return self.reference

class CommandLineOptions:
    def __init__(self, *args):
        self.parser = ap.ArgumentParser(description='Finding Js')
        self.parser.add_argument('--number-of-interactions', '-J', type=int, default=1, metavar="#J",
                help='number of exchange-interaction magnitudes to-be-included in calculations')
        self.parser.add_argument('--reference', '--noFlip', '-R', metavar='dir', required=True,
                    help='reference directory (usually noFlip/)')
        self.parser.add_argument('--units', '-U', default='meV',
                    choices=['eV', 'meV', 'Ry', 'mRy', 'He', 'mHe', 'K'],
                    help='units of energy')
        self.parser.add_argument('--elements', '--atoms', '-E', metavar='symbol', 
                                 nargs="+", required=True,
                help='Symbol of elements taken into account in calculations')
        self.parser.add_argument('--directories', '-D', metavar='dir', nargs='+', required=True,
                    help='Directories containing flipped configurations (eg. flip00000/)')
        self.opt = self.parser.parse_args(args[1:])

    def __call__(self, key):
        try:
            return self.opt.__dict__[key]
        except KeyError:
            exit(-1)
