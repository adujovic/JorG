#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
from vasprun import VaspRunXML

for d in ['Cr','CsF3Ni','Fe','Fe54','Fe54a','H','Nd']:
    print(d,'%s/vasprun.xml\t'%d,end='')
    parser = VaspRunXML('%s/vasprun.xml'%d)
    parser()
    print("E_f = %f"%parser.fermi_energy, end=' ')
    print("E_0 = %f"%parser.energy)
    for key in parser:
        print("Ion %d:\t s: % .7f p: % .7f d: % .7f f: % .7f m: % .7f"%(key,*parser[key],))
