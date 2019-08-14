# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
from POSCARloader import POSCARloader
from heisenberg import EquationSolver,NaiveHeisenberg,apply_mirrorsXYZ
from pickup.read_vasprun import MAGMOMloaderXML
import numpy as np
import argparse as ap
import re

class error:
    unexcepted = 12

class one(dict):
    def __missing__(self,key):
        return 0.0

class EnergyConverter:
    energyRatios = {'eV' :    1.0,
                    'meV':    1000.0,
                    'Ry' :    np.reciprocal(13.6056980659),
                    'mRy':    1000.0*np.reciprocal(13.6056980659),
                    'He' :    0.5*np.reciprocal(13.6056980659),
                    'mHe':    500.0*np.reciprocal(13.6056980659),
                    'K'  :    11604.51812}
    default = { 'moments'  : one(),
                'units'         : 'meV' }
    types = [ "       no moment",
              "geometric moment" ]

    @staticmethod
    def multiply(arr,*args):
        return [[ arg*element for element in arr ] for arg in args ]

    @staticmethod
    def convert(*args,**kwargs):
        # Returns array of J values with different conventions (see types)
        settings = EnergyConverter.default
        settings.update(kwargs)
        data = np.array([np.copy(args),np.copy(args)])*EnergyConverter.energyRatios[kwargs['units']]
        for i,_ in enumerate(args):
            data[1][i] *= kwargs['moments'][i]
        return data

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
        vaspruns = [ "%s/vasprun.xml"%arg for arg in args ]
        self.MAGMOMs = MAGMOMloaderXML(*vaspruns,TRAPZ=True)
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
        try:
            self.crystal  = [ [atom[0],atom[1],self.MAGMOMs.get_moments()[i+1]] for i,atom in enumerate(self.crystal) ]
        except KeyError as err:
            print(self.MAGMOMs.get_moments())
            print(err)
            exit(-1)
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
        a = len(self.systemOfEquations)
        self.solver.remove_tautologies()
        if a>len(self.systemOfEquations):
            print("A tautology!")

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

        self.Js = np.array(EnergyConverter.convert(*(self.solver.solve()),
                           moments=self.model.avgMoments, **kwargs))
        return self.Js

    def __str__(self):
        try:
            strout  = '                       '
            strout += ''.join([ "%s | "%name for name in self.model.interactionNames ]) + '\n'
        except AttributeError:
            return 'None'
        try:
            strout += ''.join([ ("  %s:\t"+len(self.Js[0])*"  % 8.3f    "+"\n")\
                                %(typeName,*self.Js[i],) for i,typeName in enumerate(self.types) ])
            strout += '  mu1 * mu2:            '
            strout += ''.join([ "% .3f        "%mu for mu in self.model.avgMoments])
        except AttributeError:
            return 'None'
        return strout

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
        if key == 'elements':
            elements = ''.join(self.opt.__dict__[key])
            elements = re.search('([A-Z][a-z]?)(.*)$',elements)
            mask     = ''
            while elements is not None:
                mask    += elements.group(1)+"$"
                elements = re.search('([A-Z][a-z]?)(.*)$',elements.group(2))
            return mask
        try:
            return self.opt.__dict__[key]
        except KeyError:
            exit(-1)
