from .format import print_vector,print_atom,print_case,print_crystal,print_moments,print_label

def show_symmetry(symmetry):
    for i in range(symmetry['rotations'].shape[0]):
        print("  --------------- %4d ---------------" % (i + 1))
        rot = symmetry['rotations'][i]
        trans = symmetry['translations'][i]
        print("  rotation:")
        for x in rot:
            print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))
        print("  translation:")
        print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))

def show_lattice(lattice):
    print("Basis vectors:")
    for vec, axis in zip(lattice, ("a", "b", "c")):
        print("%s %10.5f %10.5f %10.5f" % (tuple(axis,) + tuple(vec)))

def show_cell(lattice, positions, numbers):
    show_lattice(lattice)
    print("Atomic points:")
    for p, s in zip(positions, numbers):
        print("%2d %10.5f %10.5f %10.5f" % ((s,) + tuple(p)))

from .PeriodicTable import periodicTableElement
from sys import stdout
import numpy as np

def print_line(raport,line,linewidth=88):
    raport.write(3*"*"+(line).center(linewidth-6)+3*"*"+"\n")

def get_equivalent_line(i,j,atom,wyck):
    output = "%s: "%(atom[0])
    output += " %d "%(i+1)+" -> "
    output += " %d "%(j+1)+" W: "
    output += "%s"%(wyck)
    return output

def write_single(comment, record, crystal,
                 raport, atomDict, linewidth=88):
    uniqueWyckoffs = set(record['wyckoffs'])
    wyckoffCount =  dict.fromkeys(uniqueWyckoffs,0)
    print_line(raport,comment,linewidth)
    raport.write(linewidth*"*"+"\n")

    raport.write(3*"*"+("Spacegroup: "
                      +"%s "%(record['international'])
                      +"(%d) "%(record['number'])
                      ).center(linewidth-6)
                      +3*"*"+'\n')
    print_line(raport,"Mapping to equivalent atoms with the Wyckoff positions:",linewidth)


    for i,(j,atom,wyck) in enumerate(zip(record['equivalent_atoms'],crystal,record['wyckoffs'])):
        output = get_equivalent_line(i,j,atom,wyck)
        print_line(raport,output,linewidth)
        wyckoffCount[wyck] += 1

    raport.write(linewidth*"*"+"\n")

    output = ""
    for wyck in uniqueWyckoffs:
        output += " #%s "%(wyck)
        output += " = "+" %d "%(wyckoffCount[wyck])
    raport.write(3*"*"+(output.center(linewidth-6)+3*"*"+"\n"))

    raport.write(linewidth*"*"+"\n")

def write_report(comments,data,crystal,fileName=None, atomDict=None, linewidth=88):
    try:
        raport = open(fileName,"w+")
    except TypeError:
        raport = stdout


    raport.write(linewidth*"*"+"\n")
    raport.write(3*"*"+("Symmetry analysis").center(linewidth-6)+3*"*"+"\n")
    raport.write(linewidth*"*"+"\n")

    for i,(comment,record) in enumerate(zip(comments,data)):
        write_single(comment,record,crystal,raport,atomDict,linewidth)
    if raport is not stdout:
        raport.close()
