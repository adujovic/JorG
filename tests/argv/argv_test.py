#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,maxsize,path
path.insert(0,r'../../')
import numpy as np
from JorG.generator import generate_from_NN
import time
from datetime import datetime
import re
from aux.argv import options
import time


def main(**args):
    pass

if __name__ == '__main__':
    tracker  = -(time.time())

    currentOptions = options(*argv)

    print('Reading: cutOff', end='\t')
    print(currentOptions('cutOff'))
    print('Reading: neighbor', end='\t')
    print(currentOptions('neighbor'))
    print('Reading: Wyckoffs', end='\t')
    print(currentOptions('Wyckoffs'))
    print('Reading: reference', end='\t')
    print(currentOptions('reference'))
    print('Reading: input', end='\t')
    print(currentOptions('input'))
    print('Reading: incar', end='\t')
    print(currentOptions('incar'))
    print('Reading: output', end='\t')
    outDirName = currentOptions('output')
    if outDirName == None:
      # if output directory is not given:
      outDirName = "output/"+datetime.now().strftime("%Y%m%d%H%M%S")
    else:
      # remove multiple '/' and possible '/' at the end
      outDirName = re.sub('/+','/',outDirName)
      outDirName = re.sub('/$','',outDirName)
    print(outDirName)
    print('Reading: mask', end='\t')
    print(currentOptions('mask'))
    print('Reading: symmetry', end='\t')
    print(currentOptions('symmetry'))
    print('Reading: redundant', end='\t')
    print(currentOptions('redundant'))
    print('Reading: spin-orbit', end='\t')
    print(currentOptions('spin-orbit'))
    print('Reading: refined', end='\t')
    print(currentOptions('refined'))
    print('Reading: extra-dimentions', end='\t')
    print(currentOptions('extra-dimentions'))

    tracker += time.time()
    print("Runtime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
    exit(0)
