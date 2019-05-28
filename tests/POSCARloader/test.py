#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,maxsize,path
path.insert(0,r'../../')
import re
from JorG._tmp_loadsave import * 
import time
from itertools import product

def main(**args):
    pass

if __name__ == '__main__':
    tracker  = -(time.time())

    print("Test: parsing atom data")
    print(POSCARloader.parse_atom      ("0.0000000000 0.0000000000 4.3870399043 Cu")) 
    print(POSCARloader.parse_atomName  ("0.0000000000 0.0000000000 4.3870399043 Cu")) 
    print(POSCARloader.parse_constrains("0.0000000000 0.0000000000 4.3870399043 Cu")) 

    print(POSCARloader.parse_atom      ("-0.00000000000000 0.00000000000000 0.36270017366243 True False False"))
    print(POSCARloader.parse_atomName  ("-0.00000000000000 0.00000000000000 0.36270017366243 True False False"))
    print(POSCARloader.parse_constrains("-0.00000000000000 0.00000000000000 0.36270017366243 True False False"))

    print("Test: parsing POSCAR files")
    loader   = POSCARloader('POSCAR_exp1','POSCAR_exp2','POSCAR_exp3','POSCAR_noExistent',spam=False)
    loader.parse()
    for data,i in product(('comment','directions','cell',
                           'cellSymmetry','cellVolume',
                           'cellCenter','cellAtoms','atomNames'),range(4)):
        print("Printing %s of %d file:"%(data,i))
        try:
            print(loader(i)[data])
        except TypeError:
            print("No data in %d"%i)
    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))

