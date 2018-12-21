#!/usr/bin/python3
# -*- coding: utf-8 -*-
import numpy as np
from itertools import product

def get_configuration_penalty(configuration, flipper,
                              crystal, crystal8,
                              allFlippable, newReference):
    singlePenalty = {distance: 0 for distance in flipper['distance']} 
    for i,flip1 in enumerate(configuration):
        if flip1 < 0:
            for (a,b) in product(range(0,8),repeat=2):
                d = np.round(
                      np.linalg.norm(
                        crystal8[newReference+len(crystal)*a][1]
                       -crystal8[allFlippable[i]+len(crystal)*b][1]
                      ),
                    logAccuracy)
                if d <= cutOff:
                    if d not in singlePenalty.keys():
                        singlePenalty[d] = 1
                    else:    
                        singlePenalty[d] += 1
        for j,flip2 in enumerate(configuration):
            if(j < i):
                if flip1*flip2 < 0:
                    for (a,b) in product(range(0,8),repeat=2):
                        d = np.round(
                              np.linalg.norm(
                                crystal8[allFlippable[i]+len(crystal)*a][1]
                               -crystal8[allFlippable[j]+len(crystal)*b][1]
                              ),
                            logAccuracy)
                        if d <= cutOff:
                            if d not in singlePenalty.keys():
                                singlePenalty[d] = 1
                            else:    
                                singlePenalty[d] += 1
    penalty = 0.0                                
    for distance,numberOfFlips in singlePenalty.items():
        penalty += numberOfFlips/np.power(distance,2)
    
    return [penalty,configuration]


