from utilities.KPOINTS.kpoints import KPOINTS
from sys import argv
import numpy as np

def print_set(string,data):
    display = 42
    print(string.center(display))
    KPOINTS.print_all(data,display)

if __name__ == '__main__':
    display = 42
    for arg in argv[1:]:
        points = KPOINTS(arg)
        print("KPOINT multiplier proposition for %s:"%arg.center(display))
        print_set("Good:"               , points(np.sqrt(1e-3)))
        print_set("Decent:"             , points(np.sqrt(1e-2)))
