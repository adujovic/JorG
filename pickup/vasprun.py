# -*- coding: utf-8 -*-
import re
import numpy as np
import defusedxml.ElementTree as ET

class VaspRunXML:
    def __init__(self,vasprun='vasprun.xml',**kwargs):
        try:
            self.root     = ET.parse(vasprun).getroot()
        except ET.ParseError:
            print("File %s broken... Remove and try again!"%vasprun)
            exit(244)
        self.fermi_energy = None
        self.isOkToRead    = False
        self.partialDOS   = {}
        if 'trapez' in kwargs:
            self.trapez = kwargs['trapez']
        else:
            self.trapez        = False

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
                self.isOkToRead = True
                return
        except KeyError:
            return
        if not self.isOkToRead:
            return
        try:
            if field.attrib['name'] == "e_0_energy":
                self.isOkToRead = False
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

    def __call__(self,partial=False):
        for child in self.root.iter():
            self.get_fermi_energy(child)
            self.get_energy(child)
            self.get_partial_spin(child)
        self.calculate_moments(partial)

    def calculate_moments(self,partial):
        self.moments = {}
        for ion in self.partialDOS:
            self.calculate_moment(ion,partial)

    @staticmethod
    def DOS_below_ef(dos,ef):
        return np.delete(np.array(dos),
                np.argwhere(np.array(dos)[:,0] >= ef).flatten(),axis=0)

    def integrate(self,function,argument):
        if self.trapez:
            try:
                return np.trapz(np.sum(function,axis=1),argument)
            except np.AxisError:
                return np.trapz(function,argument)
        else:
            deltaArgument = argument[1:]-argument[:-1]
            deltaArgument = np.diag(np.insert(deltaArgument,0,deltaArgument[0]))
            return np.sum(np.dot(deltaArgument,function))

    def calculate_moment(self,ion,partial):
        ups = VaspRunXML.DOS_below_ef(self.partialDOS[ion][1],self.fermi_energy)
        dns = VaspRunXML.DOS_below_ef(self.partialDOS[ion][2],self.fermi_energy)
        if partial:
            self.moments[ion]= ([self.integrate(ups[:,1   ] - dns[:,1   ],ups[:,0]),   #moment @ s
                                 self.integrate(ups[:,2:5 ] - dns[:,2:5 ],ups[:,0]),   #moment @ p
                                 self.integrate(ups[:,5:10] - dns[:,5:10],ups[:,0]),   #moment @ d
                                 self.integrate(ups[:,10: ] - dns[:,10: ],ups[:,0]),   #moment @ f
                                 self.integrate(ups[:,1:  ] - dns[:,1:  ],ups[:,0])])  #total moment
        else:
            self.moments[ion]= self.integrate(ups[:,1:  ] - dns[:,1:  ],ups[:,0])

    def __getitem__(self, key):
        return self.moments[key]

    def __iter__(self):
        return iter(self.partialDOS)
