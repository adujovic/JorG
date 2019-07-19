#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path,argv
path.insert(0,r'../')
from geometry.voronoi import Voronoi
from POSCARloader import POSCARloader

def main(**args):
    pass

if __name__ == '__main__':
    loader   = POSCARloader(*argv[1:],spam=False)
    loader.parse()
    for i in range(len(loader)):
        voronoi = Voronoi(data=loader(i))
        voronoi.get_Voronoi_diagram(save=True,name="%s_WS.txt"%argv[1+i])
        voronoi.show("%s.png"%argv[1+i])
        voronoi.show()
