from JorG.format import print_vector,print_atom,print_case,print_crystal,print_moments,print_label
from JorG.format import standard,empty_line,line

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
    kwargs['stream'].write(3*"*"+(line).center(kwargs['linewidth']-6)+3*"*"+"\n")

def get_equivalent_line(i,j,atom,wyck):
    output = "%s: "%(atom)
    output += " %d "%(i+1)+" -> "
    output += " %d "%(j+1)+" W: "
    output += "%s"%(wyck)
    return output

def write_single(comment, record,
                 crystal, **kwargs):
    uniqueWyckoffs = set(record['wyckoffs'])
    wyckoffCount   = dict.fromkeys(uniqueWyckoffs,0)
    print_line(comment,**kwargs)
    line(kwargs['linewidth'],kwargs['stream'])

    kwargs['stream'].write(3*"*"+("Spacegroup: "
                      +"%s "%(record['international'])
                      +"(%d) "%(record['number'])
                      ).center(kwargs['linewidth']-6)
                      +3*"*"+'\n')
    print_line("Mapping to equivalent atoms with the Wyckoff positions:",**kwargs)

    for i,(j,atom,wyck) in enumerate(zip(record['equivalent_atoms'],crystal,record['wyckoffs'])):
        output = get_equivalent_line(i,j,atom,wyck)
        print_line(output,**kwargs)
        wyckoffCount[wyck] += 1

    line(kwargs['linewidth'],kwargs['stream'])

    output = ""
    for wyck in uniqueWyckoffs:
        output += " #%s "%(wyck)
        output += " = "+" %d "%(wyckoffCount[wyck])
    kwargs['stream'].write(3*"*"+(output.center(kwargs['linewidth']-6)+3*"*"+"\n"))

    line(kwargs['linewidth'],kwargs['stream'])

def write_report(comments, data,
                 crystal, **kwargs):
    kwargs = standard.fix(**kwargs)
    line(kwargs['linewidth'],kwargs['stream'])
    print_label("Symmetry analysis",**kwargs)
    kwargs['stream'].write(3*"*"+("Symmetry analysis").center(kwargs['linewidth']-6)+3*"*"+"\n")
    line(kwargs['linewidth'],kwargs['stream'])

    for i,(comment,record) in enumerate(zip(comments,data)):
        write_single(comment,record,
                     [periodicTableElement[
                         record['std_types'][
                             mapping]-1]
                         for mapping in record['mapping_to_primitive']],
                      **kwargs)
