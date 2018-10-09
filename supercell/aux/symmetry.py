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
def write_report(comment,data,crystal,fileName):
    with open(fileName,"w+") as raportFile:
        raportFile.write(str(comment))
        for i,record in enumerate(data):
            raportFile.write("\n\n*****************(%d)*********************\n\n"%(i+1))
            
            raportFile.write("Spacegroup is:\
                            \n   %s (%d)\n\n" % (record['international'],
                                                 record['number']))
            raportFile.write("  Mapping to equivalent atoms with the Wyckoff positions:\n")
            for i, (x,atom,wyck) in enumerate(zip(record['equivalent_atoms'],crystal,record['wyckoffs'])):
                raportFile.write("%s:\t%d\t->\t%d\tw: %s\n" % (atom[0],i + 1, x + 1,wyck))
        
