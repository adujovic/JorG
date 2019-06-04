#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,maxsize,path
path.insert(0,r'../../../')
import time
import numpy as np
from JorG.showcell import showCell

def main(**args):
    pass

if __name__ == '__main__':
    tracker  = -(time.time())

    points=np.array([[ 1, 1,-1],
                     [ 1, 1, 1],
                     [-1,-1, 1],
                     [-1,-1,-1]])

    plotter = showCell()
    plotter.add_sphere()
    plotter.add_polygon(points)
    plotter.show()

    for i in range(0,16,3):
        plotter = showCell(resolution=i)
        for j in range(-4,5):
            plotter.add_sphere(position=j*np.ones(3))
        plotter.show("plot_res%03d.png"%i)

    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
