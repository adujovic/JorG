# -*- coding: utf-8 -*-
import argparse as ap
import numpy as np
import re
import aux.periodic as periodic

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
        self.parser.add_argument('--period',choices=periodic.periodNames, nargs='+',
                                 help='period name (eg. 2d <=> \'%s\')'%periodic.mask3d)
        self.parser.add_argument('--block',choices=periodic.blockNames, nargs='+',
                                 help='block name (eg. P <=> \'%s\')'%periodic.maskP)

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

    def generate_part(self,MendeleyevSet,extraChar=''):
        output=''
        for el in self.opt.__dict__[MendeleyevSet]:
            output += eval("periodic.mask%s%s"%(extraChar,str(el)))
        return output

    def generate_separate(self):
        output=''
        for element in re.findall('[A-Z][a-z]?',self.opt.__dict__['elements']):
            output += '$'+element
        output += '$'
        return output

    def generate_mask(self):
        output = ''
        if self.opt.__dict__['period'] is not None:
            output += self.generate_part('period')
        if self.opt.__dict__['group'] is not None:
            output += self.generate_part('group','G')
        if self.opt.__dict__['block'] is not None:
            output += self.generate_part('block')
        if self.opt.__dict__['elements'] is not None:
            output += self.generate_separate()
        if output == '':
            return periodic.maskFull
        return output

    def __call__(self, key):
        if key not in self.keys:
            print("No key \"%s\" defined, please try: "%key)
            print("%s"%(str(self.keys)))
            exit(-300)
        elif key == 'reference':
            return self.opt.__dict__['reference'] - 1;
        elif key in self.opt.__dict__:
            return self.opt.__dict__[key]
        elif key == 'mask':
            return self.generate_mask()
