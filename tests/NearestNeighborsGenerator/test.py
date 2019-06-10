#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,maxsize,path
path.insert(0,r'../../')
import numpy as np
from JorG.generator import NearestNeighborsGenerator
import time

def main(**args):
    pass

if __name__ == '__main__':
    totalTracker  = -(time.time())
    cell            = [['Mn', np.array([0., 0., 0.])], ['Mn', np.array([1.395805, 1.395805, 1.395805])]]
    referenceAtom   = [0, np.array([0., 0., 0.])]
    directions      = [np.array([2.79161, 0.     , 1.     ]), np.array([1.     , 2.79161, 0.     ]), np.array([0.     , 1.     , 2.79161])]
    nearestNeighbor = 3
    atomNames       = ['Mn']
    wyckoffs        = 'abcdefghijklmnopqrstuvwxyz'
    atomTypeMask    = '$Mn$'
    oldMoments      = [5.0, 5.0]
    extraMultiplier = [0, 0, 0]

    generator = NearestNeighborsGenerator(cell,
                                 referenceAtom,
                                    directions)
    generator.wyckoffs         = wyckoffs
    generator.atomTypeMask     = atomTypeMask
    generator.moments          = oldMoments
    generator.extraMultiplier  = extraMultiplier

    for i in range(1,5):
        print('Nearest Neighbor #%d\t'%i,end='')
        tracker  = -(time.time())
        (cutOff,
         crystal,
         symmetryFull,
         newReference,
         copiesInEachDirection,
         wyckoffDict           ) = generator(i)
#        try:
#            (cutOff,
#             crystal,
#             symmetryFull,
#             newReference,
#             copiesInEachDirection,
#             wyckoffDict           ) = generator(i)
#            print("Test succeeded")
#            print("\t\tCrystal size: %d atoms"%len(crystal))
#        except Exception as e:
#            print(e)
#            print("Test failed")
        tracker += time.time()
        print("\t\tRuntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))

    totalTracker += time.time()
    print("\nTotal runtime: %02d:%02d:%02d"%(int(totalTracker/3600),int(totalTracker/60),int(totalTracker)))
