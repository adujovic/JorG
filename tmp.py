#!/usr/bin/python3
import re

from setuptools import setup, Extension

requirements_jorgpi = []

if __name__ == '__main__':
    try:
        with open('requirements.txt','r') as reqs:
            for line in reqs.readlines():
                package = re.sub('[><=]+.*','',line)
                package = re.sub('\s','',package)
                if len(package)>0:
                    requirements_jorgpi.append(package)
    except FileNotFoundError:
        requirements_jorgpi = [ 'numpy', 'scipy', 'spglib', 'matplotlib', 'setuptools', 'defusedxml' ]
