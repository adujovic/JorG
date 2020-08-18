#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv

from JorGpi.pickup.pickup import SmartPickUp,Reference,CommandLineOptions

if __name__ == '__main__':
    options = CommandLineOptions(*argv)

    elements = ''.join(options('elements'))
    ref = Reference(options('reference')+"/POSCAR")

    print("Running for NN=%d, \'%s\' from atom No %d:"%(options('number_of_interactions'),
                                                        elements,ref()))
    pickerUpper = SmartPickUp(options('number_of_interactions'),elements)
    pickerUpper.read(options('reference'),*options('directories'),reference=ref())

    print("Exchange interaction magnitude(s) in %s:"%options('units'))
    pickerUpper.solve(units=options('units'))
    print(pickerUpper)
