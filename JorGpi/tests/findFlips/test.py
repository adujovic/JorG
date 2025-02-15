#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import time

from JorGpi.generate.equivalent import FindFlips,find_all_distances

def main():
    pass

if __name__ == '__main__':
    crystal8      = [('Mn', np.array([ 0.000000, 0.000000, 0.000000])),
                     ('Mn', np.array([-2.791610, 0.000000, 0.000000])),
                     ('Mn', np.array([ 0.000000,-2.791610, 0.000000])),
                     ('Mn', np.array([ 0.000000, 0.000000,-2.791610])),
                     ('Mn', np.array([-2.791610,-2.791610, 0.000000])),
                     ('Mn', np.array([-2.791610, 0.000000,-2.791610])),
                     ('Mn', np.array([ 0.000000,-2.791610,-2.791610])),
                     ('Mn', np.array([-2.791610,-2.791610,-2.791610])),
                     ('Mn', np.array([ 1.395805, 1.395805, 1.395805])),
                     ('Mn', np.array([-1.395805, 1.395805, 1.395805])),
                     ('Mn', np.array([ 1.395805,-1.395805, 1.395805])),
                     ('Mn', np.array([ 1.395805, 1.395805,-1.395805])),
                     ('Mn', np.array([-1.395805,-1.395805, 1.395805])),
                     ('Mn', np.array([-1.395805, 1.395805,-1.395805])),
                     ('Mn', np.array([ 1.395805,-1.395805,-1.395805])),
                     ('Mn', np.array([-1.395805,-1.395805,-1.395805]))]
    symmetryCell  = ([(2.79161, 0.0, 0.0),
                      (0.0, 2.79161, 0.0),
                      (0.0, 0.0, 2.79161)],
                     [(0.0, 0.0, 0.0),
                      (0.5, 0.5, 0.5)],
                     [25, 25])
    reference     = 0
    referenceAtom = [0, np.array([0., 0., 0.])]
    cutOff        = 2.8
    print("Running test for:")
    for atom in crystal8:
        print("%s @ % .5f % .5f % .5f"%(atom[0],*atom[1]))

    tracker  = -(time.time())

    flipSearch = FindFlips()
    flipSearch.set_structure(symmetryCell)

    print("Unique flips:")
    flipper = flipSearch.unique(referenceAtom,cutOff)
    print(flipper)
    print("All flips")
    flipper = flipSearch.all(referenceAtom,cutOff)
    print(flipper)
    print("All possible distances until %f"%cutOff)
    print(find_all_distances(crystal8,cutOff,np.append(flipper,reference)))

    tracker += time.time()
    print("Runtime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
