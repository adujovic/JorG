# -*- coding: utf-8 -*-

from os import scandir
import re
from ctypes import CDLL,c_char_p
from JorGpi.generate.crunsetup import CrunSetup

class Crun(CrunSetup):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        super().build(*args)
        super().clean_module()
        
        self.search_for_library()

    def search_for_library(self):
        match = [re.match('^%s.*\.so'%self.name,entry.name)\
                 for entry in scandir('%s-lib'%(self.name))]
        try:
            while None in match:
                match.remove(None)
        except ValueError as err:
            if "x not in list" not in str(err):
                print('ValueError: ',str(err))
                exit(-1)
        try:
            self.library = CDLL('%s-lib/%s'%(self.name,match[0].string))
        except UnboundLocalError:
            print("I couldn't find directory")
            exit(-1)

    def __call__(self,function,*args):
        argStr = []
        for arg in args:
            if isinstance(arg,str):
                argStr.append(c_char_p(arg.encode('utf-8')))
            else:
                argStr.append(arg)
    
        f = eval("self.library.%s"%function)
        return f(*argStr)
