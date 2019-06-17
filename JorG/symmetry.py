from JorG.format import print_label
from JorG.format import standard,line

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

from JorG.PeriodicTable import periodicTableElement
import numpy as np
def print_line(line,**kwargs):
    kwargs = standard.fix(**kwargs)
    kwargs['stream'].write("|"+(line).center(kwargs['linewidth'])+"|"+"\n")

def get_equivalent_line(i,j,atom,wyck):
    output = "%s: "%(atom)
    output += " %d "%(i+1)+" -> "
    output += " %d "%(j+1)+" W: "
    output += "%s"%(wyck)
    return output

def write_single(comment, record,
                 crystal, **kwargs):
    wyckoffCount   = dict.fromkeys(set(record['wyckoffs']),0)

    print_line(comment,**kwargs)
    line(**kwargs)
    print_line("Spacegroup: %s (%d) "%(record['international'],record['number']),**kwargs)
    print_line("Mapping to equivalent atoms with the Wyckoff positions:",**kwargs)

    for i,(j,atom,wyck) in enumerate(zip(record['equivalent_atoms'],crystal,record['wyckoffs'])):
        output = get_equivalent_line(i,j,atom,wyck)
        print_line(output,**kwargs)
        wyckoffCount[wyck] += 1

    line(**kwargs)
    output = ""
    for wyck in wyckoffCount:
        output += " #%s  =  %d "%(wyck,wyckoffCount[wyck])
    print_line(output)
    line(**kwargs)

def write_report(comments, data,
                 crystal, **kwargs):
    kwargs = standard.fix(**kwargs)
    line(**kwargs)
    print_label("Symmetry analysis",**kwargs)
    line(**kwargs)

    for i,(comment,record) in enumerate(zip(comments,data)):
        write_single(comment,record,
                     [periodicTableElement[
                         record['std_types'][
                             mapping]-1]
                         for mapping in record['mapping_to_primitive']],
                      **kwargs)
