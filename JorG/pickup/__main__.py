#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv
import time
import numpy as np
from pickup.pickup import SmartPickUp,Reference,CommandLineOptions

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
    Js = np.array(Js)
    for i,typeName in enumerate(pickerUpper.types):
        print(("  %s:\t"+len(Js[0])*"% 11.7f ")%(typeName,*Js[i],))
    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
