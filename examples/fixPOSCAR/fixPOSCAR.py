#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
from POSCARloader import POSCARloader
def save_vanilla_POSCAR(fileName,data):
#    """
#        Saving data to POSCAR file
#{'comment': 'Re/117 - (A3) - HCP [A2] A3 (117) (htqc', 'directions': [array([ 1.38476066, -2.39847581, -0.        ]), array([1.38476066, 2.39847581, 0.        ]), array([0.        , 0.        , 4.47184109])], 'cell': [[0, array([0., 0., 0.])], [0, array([ 1.38476066, -0.79949194,  2.23592055])]], 'cellSymmetry': ([(1.3847606573816698, -2.398475814907534, -0.0), (1.3847606573816698, 2.398475814907534, 0.0), (0.0, 0.0, 4.471841090100609)], [(0.0, -0.0, -0.0), (0.6666666666666714, 0.3333333333333428, 0.5)], [75, 75]), 'cellVolume': 29.704785298855388, 'cellCenter': array([1.38476066, 0.        , 2.23592055]), 'cellAtoms': array([2]), 'atomNames': ['Re']}
#                                    """
    with open(fileName,"w+") as vaspFile:
        vaspFile.write(data['comment'])
        vaspFile.write("\n1.0\n")
        for d in data['directions']:
            for field in d:
                vaspFile.write("  %.10f"%(field))
            vaspFile.write("\n")
        for atomName in data['atomNames']:
            vaspFile.write("%s "%atomName)
        vaspFile.write("\n")
        for atomNumber in data['cellAtoms']:
            vaspFile.write("%d "%atomNumber)
        vaspFile.write("\nDirect\n")
        for atom in data['cell']:
            for vasp in atom[1]:
                vaspFile.write(" %.10f "%vasp)
            try:
                vaspFile.write(" %s\n"%data['atomNames'][atom[0]])
            except TypeError:
                vaspFile.write(" %s\n"%atom[0])
        vaspFile.write("\n")

loader=POSCARloader('POSCAR')
loader.parse()
a = loader()
print(a)
save_vanilla_POSCAR('POSCAR_fixed',a)
