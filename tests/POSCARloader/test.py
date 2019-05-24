#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,maxsize,path
path.insert(0,r'../../')
import re
from JorG.loadsave import * 
import time

def main(**args):
    pass

if __name__ == '__main__':
    tracker  = -(time.time())
    loader   = POSCARloader('POSCAR_exp1','POSCAR_exp2','POSCAR_exp3','POSCAR_noExistent',spam=False)
    print(loader.find_comment(loader()))
    print(loader.find_cell(loader()))
    print(loader.find_directions(loader()))

    print(loader.find_comment(loader(1)))
    print(loader.find_cell(loader(1)))
    print(loader.find_directions(loader(1)))

    print(POSCARloader.parse_atom      ("0.0000000000 0.0000000000 4.3870399043 Cu")) 
    print(POSCARloader.parse_atomName  ("0.0000000000 0.0000000000 4.3870399043 Cu")) 
    print(POSCARloader.parse_constrains("0.0000000000 0.0000000000 4.3870399043 Cu")) 

    print(POSCARloader.parse_atom      ("-0.00000000000000 0.00000000000000 0.36270017366243 True False False"))
    print(POSCARloader.parse_atomName  ("-0.00000000000000 0.00000000000000 0.36270017366243 True False False"))
    print(POSCARloader.parse_constrains("-0.00000000000000 0.00000000000000 0.36270017366243 True False False"))

    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))

