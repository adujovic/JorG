from sys import path
path.insert(0,r'../../')
from sys import argv
from POSCARloader import POSCARloader
import numpy as np
import argparse as ap

class CommandLineOptions:
    def __init__(self, *args):
        self.parser = ap.ArgumentParser(description='Converting POSCAR files')
        self.parser.add_argument('--convert-to', '-T', default='Direct',
                    choices=['d','c','k','D','C','K','Direct','Carthesian'],
                    help='Convertion To')
        self.parser.add_argument('--files', '-F', metavar='POSCAR_file',
                                 nargs="+", required=True,
                help='POSCAR files to convert')
        self.opt = self.parser.parse_args(args[1:])

    def __call__(self, key):
        try:
            return self.opt.__dict__[key]
        except KeyError:
            exit(-1)

def save_vanilla_poscar(filename,data,option='D'):
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
        if option == 'D':
            vaspFile.write("\nDirect\n")
        else:
            vaspFile.write("\nCarthesian\n")
        for atom in data['cell']:
            if option == 'D':
                vaspFile.write(3*" %.10f "%(*np.dot(invDirections,atom[1]),))
            else:
                vaspFile.write(3*" %.10f "%(*atom[1],))
            try:
                vaspFile.write(" %s\n"%data['atomNames'][atom[0]])
            except TypeError:
                vaspFile.write(" %s\n"%atom[0])
        vaspFile.write("\n")

if __name__ == '__main__':
    options = CommandLineOptions(*argv)
    for poscar in options('files'):
        loader = POSCARloader(poscar)
        loader.parse()
        save_vanilla_poscar('%s.fxd'%poscar,loader(),options('convert_to').capitalize()[0])
