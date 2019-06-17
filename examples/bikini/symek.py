#!/usr/bin/python3
# -*- coding: utf-8 -*-
from sys import argv,path
path.insert(0,r'../../')
from JorGpi.run import JorGpi

if __name__ == '__main__':
    engine = JorGpi(*argv)
    engine.initialize_new_cell()
    engine.generate_possible_configurations()
    exit()
