#!/usr/bin/python3
import re
import numpy as np
import defusedxml.ElementTree as ET

class VaspRunXML:
    def __init__(self,vasprun='vasprun.xml'):
        self.root         = ET.parse(vasprun).getroot()
        self.fermi_energy = None
        self.ISOK2READ    = False
        self.partialDOS   = {}
        self.TRAPZ        = False

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
                self.partialDOS[self.index] = { 1: [], 2: [] }
        except KeyError:
            return
        try:
            if field.tag == 'set' and 'spin' in field.attrib['comment']:
                self.spin = int(re.search('\d',field.attrib['comment']).group(0))
        except KeyError:
           return
        try:
            if field.tag == 'r':
                self.partialDOS[self.index][self.spin].append(np.fromstring(field.text,sep=' '))
        except AttributeError:
           return

    def __len__(self):
        return len(self.partialDOS)

    def __call__(self):
        for child in self.root.iter():
            self.get_fermi_energy(child)
            self.get_energy(child)
            self.get_partial_spin(child)
        self.calculate_moments()

    def calculate_moments(self):
        self.moments = {}
        for ion in self.partialDOS:
            self.calculate_moment(ion)

    @staticmethod
    def DOS_below_ef(dos,ef):
        return np.delete(np.array(dos),
                np.argwhere(np.array(dos)[:,0] >= ef).flatten(),axis=0)

    def integrate(self,f,x):
        if self.TRAPZ:
            try:
                return np.trapz(np.sum(f,axis=1),x)
            except np.AxisError:
                return np.trapz(f,x)
        else:
            dx = x[1:]-x[:-1]
            dx = np.diag(np.insert(dx,0,dx[0]))
            return np.sum(np.dot(dx,f))

    def calculate_moment(self,ion):
        ups = VaspRunXML.DOS_below_ef(self.partialDOS[ion][1],self.fermi_energy)
        dns = VaspRunXML.DOS_below_ef(self.partialDOS[ion][2],self.fermi_energy)
        self.moments[ion]= ([self.integrate(ups[:,1   ] - dns[:,1   ],ups[:,0]),   #moment @ s
                             self.integrate(ups[:,2:5 ] - dns[:,2:5 ],ups[:,0]),   #moment @ p
                             self.integrate(ups[:,5:10] - dns[:,5:10],ups[:,0]),   #moment @ d
                             self.integrate(ups[:,10: ] - dns[:,10: ],ups[:,0]),   #moment @ f
                             self.integrate(ups[:,1:  ] - dns[:,1:  ],ups[:,0])])  #total moment

    def __getitem__(self, key):
        return self.moments[key]

    def __iter__(self):
        return iter(self.partialDOS)
