from sys import path
path.insert(0,r'../../')
from sys import argv
from POSCARloader import POSCARloader
import numpy as np


def save_vanilla_poscar(filename,data):
    with open(filename,"w+") as vaspFile:
        vaspFile.write(data['comment'])
        vaspFile.write("\n1.0\n")
        invDirections = np.linalg.inv(data['directions'])
        for direction in data['directions']:
            vaspFile.write(3*"  %.10f"%(*direction,)+"\n")
        for atomName in data['atomNames']:
            vaspFile.write("%s "%atomName)
        vaspFile.write("\n")
        for atomNumber in data['cellAtoms']:
            vaspFile.write("%d "%atomNumber)
        vaspFile.write("\nDirect\n")
        for atom in data['cell']:
            vaspFile.write(3*" %.10f "%(*np.dot(invDirections,atom[1]),))
            try:
                vaspFile.write(" %s\n"%data['atomNames'][atom[0]])
            except TypeError:
                vaspFile.write(" %s\n"%atom[0])
        vaspFile.write("\n")

if __name__ == '__main__':
    for arg in argv[1:]:
        loader = POSCARloader(arg)
        loader.parse()
        save_vanilla_poscar('%s.fxd'%arg,loader())
