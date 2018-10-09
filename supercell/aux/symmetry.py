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
from sys import stdout
def write_report(comment,data,crystal,fileName=None, atomDict=None):
    if fileName is None:
        raport = stdout
    else:
        raport = open(fileName,"w+")
    raport.write(str(comment))
    for i,record in enumerate(data):
        raport.write("\n\n*****************(%d)*********************\n\n"%(i+1))
        
        raport.write("Spacegroup is:\
                        \n   %s (%d)\n\n" % (record['international'],
                                             record['number']))
        raport.write("  Mapping to equivalent atoms with the Wyckoff positions:\n")
        if atomDict is None:
            for i, (x,atom,wyck) in enumerate(zip(record['equivalent_atoms'],crystal,record['wyckoffs'])):
                raport.write("%s:\t%d\t->\t%d\tw: %s\n" % (atom[0],i + 1, x + 1,wyck))
        else: 
            for i, (x,atom,wyck) in enumerate(zip(record['equivalent_atoms'],crystal,record['wyckoffs'])):
                raport.write("%s:\t%d\t->\t%d\tw: %s\n" % (atomDict[atom[0]],i + 1, x + 1,wyck))
    if raport is not stdout:
        raport.close()
    
