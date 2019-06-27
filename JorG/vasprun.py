#!/usr/bin/python3
import re
import numpy as np
import xml.etree.ElementTree as ET

class VaspRunXML:
    def __init__(self,vasprun='vasprun.xml'):
        tree              = ET.parse(vasprun)
        self.root         = tree.getroot()
        self.fermi_energy = None
        self.ISOK2READ    = False
        self.ions         = {}

    def get_fermi_energy(self,field):
        if self.fermi_energy is not None:
            return
        try:
            if field.attrib['name'] == 'efermi':
                self.fermi_energy = float(field.text)
        except KeyError:
            return

    def get_energy(self,field): 
        if field.tag != 'energy':
            return
        for child in field.iter():
            self.get_energy_child(child)

    def get_energy_child(self,field):
        try:
            if field.attrib['name'] == 'alphaZ':
                self.ISOK2READ = True
                return
        except KeyError:
            return
        if not self.ISOK2READ:
            return
        try:
            if field.attrib['name'] == "e_0_energy":
                self.ISOK2READ = False
                self.energy    = float(field.text)
        except KeyError:
            return

    def get_partial_spin(self,field):
        if field.tag != 'partial':
            return
        for child in field.iter():
            self.get_ion_spin(child)

    def get_ion_spin(self,field):
        try:
            if field.tag == 'set' and 'ion' in field.attrib['comment']:
                self.index = int(re.search('\d+',field.attrib['comment']).group(0))
                self.ions[self.index] = { 1: [], 2: [] }
        except KeyError:
            return
        try:
            if field.tag == 'set' and 'spin' in field.attrib['comment']:
                self.spin = int(re.search('\d',field.attrib['comment']).group(0))
        except KeyError:
           return
        try:
            if field.tag == 'r':
                self.ions[self.index][self.spin].append(np.fromstring(field.text,sep=' '))
        except AttributeError:
           return

    def __call__(self):
        for child in self.root.iter():
            self.get_fermi_energy(child)
            self.get_energy(child)
            self.get_partial_spin(child)
        self.calculate_moments()

    def calculate_moments(self):
        self.moments = {}
        for ion in self.ions:
            self.calculate_moment(ion)

    def calculate_moment(self,ion):
        ups = np.delete(
                np.array(self.ions[ion][1]),
                         np.argwhere(
                             np.array(
                                 self.ions[ion][1])[:,0] > self.fermi_energy).flatten(),axis=0)
        dns = np.delete(
                np.array(self.ions[ion][2]),
                         np.argwhere(
                             np.array(
                                 self.ions[ion][2])[:,0] > self.fermi_energy).flatten(),axis=0)
        deltaEnergy = ups[1:,0]-ups[:-1,0]
        deltaEnergy = np.diag(
                        np.insert(deltaEnergy,0,deltaEnergy[0]))
        self.moments[ion]= ([np.sum(np.dot(deltaEnergy,ups[:,1   ] - dns[:,1   ])),   #moment @ s
                             np.sum(np.dot(deltaEnergy,ups[:,2:5 ] - dns[:,2:5 ])),   #moment @ p
                             np.sum(np.dot(deltaEnergy,ups[:,5:10] - dns[:,5:10])),   #moment @ d
                             np.sum(np.dot(deltaEnergy,ups[:,10: ] - dns[:,10: ])),   #moment @ f
                             np.sum(np.dot(deltaEnergy,ups[:,1:  ] - dns[:,1:  ]))])  #total moment

    def __getitem__(self, key):
        return self.moments[key]

    def __iter__(self):
        return iter(self.ions)
