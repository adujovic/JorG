import os
import re
import setuptools
import importlib
from abc import ABC, abstractmethod

class CrunSetup(ABC):
    compilerOptions = { 'include_dirs' : ['/usr/local/include'] }
    builddir = 'build'

    def __init__(self, *args, **kwargs):
        importlib.reload(setuptools)
        try:
            os.mkdir(self.builddir)
        except Exception as ext:
            print(ext,": ignoring")
        self.files = args
        self.compilerOptions.update(kwargs)

        if 'name' in kwargs:
            self.name = kwargs['name']
            del self.compilerOptions['name']
        else:
            self.name = re.sub('^.*/','',args[0])
            self.name = re.sub('\..*$','',self.name)

    @abstractmethod
    def __call__(self):
        pass

    def build(self, *args):
        self.module = setuptools.Extension(self.name, sources=list(args), **self.compilerOptions)
        setuptools.setup(
                script_args=[self.builddir,'--build-lib=%s-lib'%(self.name)],
                ext_modules=[self.module])

    def clean_module(self):
        setuptools.setup(name=self.name,script_args=['clean'],
                         ext_modules = [self.module])

    def __del__(self):
        setuptools.setup(script_args=['clean','--all','--build-lib=%s-lib'%(self.name)],
                         ext_modules = [self.module])
        importlib.reload(setuptools)
