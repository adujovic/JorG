#!/usr/bin/python3
# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../../')
from vasprun import VaspRunXML
import numpy as np

class MAGMOMloaderXML:
    def __init__(self,*args,**kwargs):
        self.vasprunFiles = args
        self.data         = [{} for arg in args]
        self.settings     = kwargs

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def parse(self):
        for i,vasprun in enumerate(self.vasprunFiles):
            parser = VaspRunXML(vasprun,**self.settings)
            parser()
            self.data[i]['fermi']   = parser.fermi_energy
            self.data[i]['energy']  = parser.energy
            self.data[i]['moments'] = parser.moments
            del parser

    def get_energy(self,idx=0):
        return self.data[idx]['energy']

    def get_moments(self,idx=0):
        return self.data[idx]['moments']

    def get_average_magnitude(self,idx=0,indices=None):
        if indices is not None:
            return np.average(np.abs(np.take(list(self.data[idx]['moments'].values()),indices)))
        return np.average(np.abs(np.array(list(self.data[idx]['moments'].values()))))

    def __call__(self,idx=0):
        return self.data[idx]
