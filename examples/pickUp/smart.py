#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,path
path.insert(0,r'../../')
import time
import argparse as ap
import numpy as np
import re
from JorG.pickup import SmartPickUp,Reference

def main(**args):
    pass

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


if __name__ == '__main__':
    tracker  = -(time.time())
    options = CommandLineOptions(*argv)

    elements = ''
    for e in options('elements'):
        elements += e
    ref = Reference(options('reference')+"/POSCAR")

    numOfJs = options('number_of_interactions')
    print("Running for NN=%d, \'%s\' from atom No %d:"%(numOfJs,elements,ref()))
    pickerUpper = SmartPickUp(numOfJs,elements)
    pickerUpper.read(options('reference'),*options('directories'),reference=ref())
    print("Exchange interaction magnitude(s) in %s:"%options('units'))
    Js = pickerUpper.solve(units=options('units'))
    print((len(Js)*"% 11.7f ")%(*Js,))
    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
