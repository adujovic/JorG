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

    loader   = POSCARloader(*POSCARs,spam=False)
    loader.parse()

    pickerUpper = SmartPickUp(1,2.46,'Fe')
    pickerUpper.read(*argv[1:])
    pickerUpper.get_system_of_equations()

    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
