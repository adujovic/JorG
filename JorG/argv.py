# -*- coding: utf-8 -*-
import argparse as ap
import numpy as np
import re
from sys import path
path.insert(0,r'../')
import JorG.Masks as periodic

class error:
    access = 123

class options:
    keys = ["cutOff", "neighbor", "Wyckoffs", "reference",
            "input", "incar", "output", "mask", "symmetry",
            "refined", "redundant", "extra-dimentions", "spin-orbit"]

    def __init__(self, *args):
        self.parser = ap.ArgumentParser(description='Find minimal number of unique spin-flips')
        self.typeOfRange = self.parser.add_mutually_exclusive_group()
        self.add_arguments_io()
        self.add_arguments_geometric()
        self.add_arguments_elements()
        self.add_arguments_optional()

        self.opt = self.parser.parse_args(args[1:])
        self.fix_data()

    def fix_data(self):
        self.opt.__dict__['reference'] -= 1
        self.opt.__dict__['mask'] = self.generate_mask()

    def add_arguments_geometric(self):
        self.typeOfRange.add_argument('--cutOff', '-R', default=None, type=np.float,
                                 help='a cut-off distance (in Å) for calculations')
        self.typeOfRange.add_argument('--neighbor', '-N', default=None, type=int,
                                 help='a rank of the last Neighbor taken into account')
        self.parser.add_argument('--Wyckoffs', '-W', default='abcdefghijklmnopqrstuvwxyz',
                                 help='narrows down the atomic selection to the atoms in positions defined by string (eg. \'abc\')')
        self.parser.add_argument('--reference', '-r', default=-1, type=int,
                                 help='number of reference atom in inputFile')

    def add_arguments_io(self):
        self.parser.add_argument('--input', '-i', default='POSCAR',
                                 help='input POSCAR file')
        self.parser.add_argument('--incar', '--INCAR', '-I', default='INCAR',
                                 help='input INCAR file')
        self.parser.add_argument('--output', '-o', default=None,
                                 help='output directory')

    def add_arguments_elements(self):
        self.parser.add_argument('--elements','-E',
                                 help='string of all elements taken into account (eg. \'CuO\')')
        self.parser.add_argument('--group',choices=range(1,19), type=int, nargs='+',
                                 help='group number (eg. 1 <=> \'HLiNaKRbCsFr\')')
        self.parser.add_argument('--period',choices=periodic.periods, nargs='+',
                                 help='period name (eg. 3d <=> \'%s\')'%periodic.periods['3d'])
        self.parser.add_argument('--block',choices=periodic.blocks, nargs='+',
                                 help='block name (eg. P <=> \'%s\')'%periodic.blocks['P'])

    def add_arguments_optional(self):
        self.parser.add_argument('--symmetry', '-S', action='store_true',
                                 help='symmetry run only (default False)')
        self.parser.add_argument('--redundant', action='store_true',
                                 help='creates a redundant system of equations for final calculation of the Heisenberg exchange interaction (default False)')
        self.parser.add_argument('--spin-orbit', '--SOC', action='store_true',
                                 help='(work-in-progress) is sping-orbit coupling enabled (default False)')
        self.parser.add_argument('--refined',  action='store_true',
                                 help='should use refined supercell (default False)')
        self.parser.add_argument('--extra-dimentions','-X', default=None, action='store',dest='extra-dimentions',
                                 help='string \"X Y Z\" of extra cell copies in each directions (eg. \"0 0 1\")')

    def generate_separate(self):
        output=''
        for element in re.findall('[A-Z][a-z]?',self.opt.__dict__['elements']):
            output += '$'+element
        output += '$'
        return output

    @staticmethod
    def concatenate(keys,values):
        output = ''
        for key in keys:
            output += str(values[key])
        return output

    def generate_mask(self):
        output = ''
        if self.opt.__dict__['period'] is not None:
            output += self.concatenate(self.opt.__dict__['period'],
                                        periodic.periods)
        if self.opt.__dict__['group'] is not None:
             output += self.concatenate(self.opt.__dict__['group'],
                                        periodic.groups)
        if self.opt.__dict__['block'] is not None:
             output += self.concatenate(self.opt.__dict__['block'],
                                        periodic.blocks)
        if self.opt.__dict__['elements'] is not None:
            output += self.generate_separate()
        if output == '':
            return periodic.maskFull
        return output

    def __call__(self, key):
        try:
            return self.opt.__dict__[key]
        except KeyError:
            print("No key \"%s\" defined, please try: "%key)
            print("%s"%(str(self.keys)))
            exit(error.access)
