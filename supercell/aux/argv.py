# -*- coding: utf-8 -*-
import argparse as ap
import numpy as np
import aux.periodic as periodic

class options:
    keys = ["cutOff", "neighbor", "Wyckoffs", "reference", "input", "output", "mask"]
    def __init__(self, *args):
        self.parser = ap.ArgumentParser(description='Find minimal number of the spin-flips')
        typeOfRange = self.parser.add_mutually_exclusive_group()
        typeOfRange.add_argument('--cutOff', '-R', default=-1.0, type=np.float,
                                 help='cut-off distance (in Å) for calculations')
        typeOfRange.add_argument('--neighbor', '-N', default=-1, type=int,
                                 help='rank of the last Neighbor taken into account')
        self.parser.add_argument('--Wyckoffs', '-W', default='abcdefghijklmnopqrstuvwxyz',
                                 help='narrow atoms selection to the ones occupied by positions given by string')
        self.parser.add_argument('--symmetry', '-S', action='store_true',
                                 help='symmetry run only')
        self.parser.add_argument('--reference', '-r', default=0, type=int,
                                 help='number of reference atom in inputFile (from 0)')
        self.parser.add_argument('--input', '-i', default='POSCAR',
                                 help='input POSCAR file')
        self.parser.add_argument('--output', '-o', default='output/POSCAR',
                                 help='output POSCAR file')
        typeOfMask = self.parser.add_mutually_exclusive_group()
        typeOfMask.add_argument('--elements','-E', default=periodic.maskFull,
                                help='string of all elements taken into account (eg. \'CuO\')')
        typeOfMask.add_argument('--group',choices=range(1,19), type=int,
                                help='group number (eg. 1 <=> \'HLiNaKRbCsFr\')')
        typeOfMask.add_argument('--period',choices=periodic.periodNames,
                                help='period name (eg. 2d <=> \'%s\')'%periodic.mask3d)
        typeOfMask.add_argument('--block',choices=periodic.blockNames,
                                help='block name (eg. P <=> \'%s\')'%periodic.maskP)

        self.opt = self.parser.parse_args(args[1:])

    def __str__(self):   
        from textwrap import wrap
        output = ""
        for opt in self.opt.__dict__:
            output += "\n".join(wrap("%s =  %s"%("{:<10}".format(str(opt)),str(self.opt.__dict__[opt])),70,subsequent_indent=14*" "))
            output += "\n"
        return output[:-1]
    
    def __call__(self, key):
        if key not in self.keys:
            print("No key %s defined, please try: "%key)
            output = ""
            for k in self.keys:
                output += k + ", "
            print(output[:-2])
            exit(-300)
        elif key == 'neighbor':
            neighbor = self.opt.__dict__['neighbor'] 
            if self.opt.__dict__['neighbor'] < 0:
                if self.opt.__dict__['cutOff'] < 0.0:
                    return 1
                else:
                    print("No neighbor defined: use cutOff")
                    return None                    
            return neighbor 
        elif key == 'cutOff':
            cutOff = self.opt.__dict__['cutOff'] 
            if cutOff < 0.0:
                print("No cutOff defined: use neighbor")
                return None                    
            return cutOff    
        elif key in self.opt.__dict__:
            return self.opt.__dict__[key]
        elif key == 'mask':
            if self.opt.__dict__['period'] is not None:
                return eval("periodic.mask"+self.opt.__dict__['period'])
            elif self.opt.__dict__['group'] is not None:
                return eval("periodic.maskG"+str(self.opt.__dict__['group'])) 
            elif self.opt.__dict__['block'] is not None:
                return eval("periodic.mask"+self.opt.__dict__['block'])
            else:
                return self.opt.__dict__['elements'] 
           
                
