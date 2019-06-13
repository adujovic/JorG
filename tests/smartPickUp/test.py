#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,path
path.insert(0,r'../../')
import time
import numpy as np
from JorG.loadsave import POSCARloader
from JorG.pickup import SmartPickUp
from itertools import product

def main(**args):
    pass

if __name__ == '__main__':
    tracker  = -(time.time())

    POSCARs = [ "%s/POSCAR"%arg for arg in argv[1:] ]
    if not POSCARs:
        print("No files provided")
        exit(-1)

    loader   = POSCARloader(*POSCARs,spam=False)
    loader.parse()

    print("Running for NN=1, \'Fe\':")
    pickerUpper = SmartPickUp(1,'Fe')
    pickerUpper.read(*argv[1:])
    for unit in ['eV', 'meV', 'Ry', 'mRy', 'He', 'mHe', 'K']:
        print("Exchange interaction in %s:"%unit)
        Js = pickerUpper.solve(units=unit)
        print((len(Js)*"% 11.7f ")%(*Js,))

    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
