#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv
import time
from pickup.pickup import SmartPickUp,Reference,CommandLineOptions

if __name__ == '__main__':
    tracker  = -(time.time())
    options = CommandLineOptions(*argv)

    elements = ''.join(options('elements'))
    ref = Reference(options('reference')+"/POSCAR")

    print("Running for NN=%d, \'%s\' from atom No %d:"%(options('number_of_interactions'),elements,ref()))
    pickerUpper = SmartPickUp(options('number_of_interactions'),elements)
    pickerUpper.read(options('reference'),*options('directories'),reference=ref())

    print("Exchange interaction magnitude(s) in %s:"%options('units'))
    pickerUpper.solve(units=options('units'))
    print(pickerUpper)

    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
