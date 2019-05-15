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

from .periodic import periodicTableElement
from .format import *
from sys import stdout
import numpy as np
def write_report(comments,data,crystal,fileName=None, atomDict=None, linewidth=88):
    if fileName is None:
        raport = stdout
    else:
        raport = open(fileName,"w+")
   
    offset=len(color.BF+color.YELLOW+color.INV+color.END)
    raport.write(linewidth*"*"+"\n")
    raport.write(3*"*"+(color.BF+color.YELLOW+color.INV+"Symmetry analysis"+color.END).center(linewidth-6+offset)+3*"*"+"\n")
    raport.write(linewidth*"*"+"\n")

    for i,(comment,record) in enumerate(zip(comments,data)):
        uniqueWyckoffs = set(record['wyckoffs'])
        wyckoffCount =  dict.fromkeys(uniqueWyckoffs,0)
        raport.write(color.INV+3*"*"+(comment).center(linewidth-6)+3*"*"+color.END+"\n")
        raport.write(linewidth*"*"+"\n")
        
        raport.write(3*"*"+("Spacegroup: "+color.BF+color.DARKMAGENTA
                          +"%s "%(record['international'])
                          +color.DARKYELLOW+"(%d) "%(record['number'])
                          +color.END).center(linewidth-6+len(color.BF+color.DARKMAGENTA+color.DARKYELLOW+color.END))
                          +3*"*"+'\n')
        raport.write(3*"*"+("Mapping to equivalent atoms with the Wyckoff positions:".center(linewidth-6)+3*"*"+'\n'))

        offset=len(color.BF+color.DARKYELLOW+color.GRAY+color.DARKGRAY+color.GRAY+color.END+color.DARKBLUE+color.END) 
        if atomDict is None:
            for i, (x,atom,wyck) in enumerate(zip(record['equivalent_atoms'],crystal,record['wyckoffs'])):
                raport.write(3*"*"+(color.BF+color.DARKYELLOW+"%s: "%(atom[0])
                                   +color.GRAY+" %d "%(i+1)+color.DARKGRAY+" -> "
                                   +color.GRAY+" %d "%(x+1)+color.END+" W: "
                                   +color.DARKBLUE+"%s"%(wyck)+color.END).center(linewidth-6+offset)+3*"*"+"\n")
                wyckoffCount[wyck] += 1 
        else: 
            for i, (x,atom,wyck) in enumerate(zip(record['equivalent_atoms'],crystal,record['wyckoffs'])):
                raport.write(3*"*"+(color.BF+color.DARKYELLOW+"%s: "%(atomDict[atom[0]])
                                   +color.GRAY+" %d "%(i+1)+color.DARKGRAY+" -> "
                                   +color.GRAY+" %d "%(x+1)+color.END+" W: "
                                   +color.DARKBLUE+"%s"%(wyck)+color.END).center(linewidth-6+offset)+3*"*"+"\n")
                wyckoffCount[wyck] += 1 

        raport.write(linewidth*"*"+"\n")

        part=len(color.BF+color.DARKRED+color.END+color.BF+color.DARKGREEN+color.END)
        offset=0
        output = ""
        for wyck in uniqueWyckoffs:
            output += color.BF+color.DARKRED+" #%s "%(wyck)+color.END+" = "+color.BF+color.DARKGREEN+" %d "%(wyckoffCount[wyck])+color.END
            offset+=part
        raport.write(3*"*"+(output.center(linewidth-6+offset)+3*"*"+"\n"))

        raport.write(linewidth*"*"+"\n")
    if raport is not stdout:
        raport.close()
    
