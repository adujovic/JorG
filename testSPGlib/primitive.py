#!/usr/bin/python3

import sys
import numpy as np
import spglib

def show_lattice(lattice):
    print("Basis vectors:")
    for vec, axis in zip(lattice, ("a", "b", "c")):
        print("%s %10.5f %10.5f %10.5f" % (tuple(axis,) + tuple(vec)))


def show_cell(lattice, positions, numbers):
    show_lattice(lattice)
    print("Atomic points:")
    for p, s in zip(positions, numbers):
        print("%2d %10.5f %10.5f %10.5f" % ((s,) + tuple(p)))


bcc_silicon_dist = ([(3.97, 0, 0),
                 (0, 4.03, 0),
                 (0, 0, 3.99)],
                [(0, 0, 0),
                 (0.501, 0.497, 0.5)],
                [14, 14]) # 14-> Si

print("[refine_cell]")
print(" Refine distorted bcc_silicon structure")
lattice, positions, numbers = spglib.refine_cell(bcc_silicon_dist,
                                                 symprec=1e-1)
show_cell(lattice, positions, numbers)
print('')

print("[find_primitive]")
print(" Find primitive distorted bcc_silicon structure")
lattice, positions, numbers = spglib.find_primitive(bcc_silicon_dist,
                                                    symprec=1e-15)
show_cell(lattice, positions, numbers)
print('')

print(" Find primitive from refined distorted bcc_silicon structure:")
bcc_silicon_std = spglib.refine_cell(bcc_silicon_dist,
                                     symprec=1e-1)
lattice, positions, numbers = spglib.find_primitive(bcc_silicon_std,
                                                    symprec=1e-15)
show_cell(lattice, positions, numbers)
print('')

print("[standardize_cell]")
print(" Standardize distorted bcc_silicon structure:")
print(" (to_primitive=0 and no_idealize=0)")
lattice, positions, numbers = spglib.standardize_cell(bcc_silicon_dist,
                                                      to_primitive=0,
                                                      no_idealize=0,
                                                      symprec=1e-1)
show_cell(lattice, positions, numbers)
print('')

print(" Standardize distorted bcc_silicon structure:")
print(" (to_primitive=0 and no_idealize=1)")
lattice, positions, numbers = spglib.standardize_cell(bcc_silicon_dist,
                                                      to_primitive=0,
                                                      no_idealize=1,
                                                      symprec=1e-1)
show_cell(lattice, positions, numbers)
print('')

print(" Standardize distorted bcc_silicon structure:")
print(" (to_primitive=1 and no_idealize=0)")
lattice, positions, numbers = spglib.standardize_cell(bcc_silicon_dist,
                                                      to_primitive=1,
                                                      no_idealize=0,
                                                      symprec=1e-1)
show_cell(lattice, positions, numbers)
print('')

print(" Standardize distorted bcc_silicon structure:")
print(" (to_primitive=1 and no_idealize=1)")
lattice, positions, numbers = spglib.standardize_cell(bcc_silicon_dist,
                                                      to_primitive=1,
                                                      no_idealize=1,
                                                      symprec=1e-1)
show_cell(lattice, positions, numbers)
print('')

