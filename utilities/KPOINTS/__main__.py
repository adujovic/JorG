from utilities.KPOINTS.kpoints import KPOINTS
from sys import argv
import numpy as np

def print_set(string,data):
    DISPLAY = 42
    print(string.center(DISPLAY))
    KPOINTS.print_all(data,DISPLAY)

if __name__ == '__main__':
    display = 42
    for arg in argv[1:]:
        points = KPOINTS(arg)
        print("KPOINT multiplier proposition for %s:"%arg.center(display))
        print_set("Good:"               , points(np.sqrt(1e-3)))
        print_set("Decent:"             , points(np.sqrt(1e-2)))
        print_set("Plausible:"          , points(np.sqrt(1e-1)))
        print_set("Borderline terrible:", points(np.sqrt( 0.5)))
        print_set("Just don't:"         , points(np.sqrt( 1.0)))
