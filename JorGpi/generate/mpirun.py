from os import environ
from glob import glob
import subprocess as sp
from JorGpi.generate.crunsetup import CrunSetup

class Mpirun(CrunSetup):
    def __init__(self, *args, **kwargs):
        if 'nmpi' in kwargs:
            self.nmpi = kwargs['nmpi']
            del kwargs['nmpi']
        else:
            self.nmpi = 1

        if 'mpicompiler' in kwargs:
            self.compiler = kwargs['mpicompiler']
            del kwargs['mpicompiler']
        else:
            self.compiler = None

        super().__init__(*args, **kwargs)

        environ_backup = dict(environ)

        self.set_compiler()
        super().build(*args)

        environ.clear()
        environ.update(environ_backup)
        
        self.executable = 'start.exe'
        self.create_executable(**kwargs)

        super().clean_module()

    def set_compiler(self):
        mpicompilers = ['mpicxx', 'mpiCC', 'mpic++']
        if self.compiler not in mpicompilers:
            if 'CC' in environ and environ['CC'] in mpicompilers:
                self.compiler = environ['CC']
            elif 'CXX' in environ and environ['CXX'] in mpicompilers:
                self.compiler = environ['CXX']
            else:
                self.compiler = 'mpicxx'

        environ['CC'] = environ['CXX'] = self.compiler

    def create_executable(self, **kwargs):
        path = self.builddir + "/**/*.o"
        obj_files = glob(path, recursive=True)

        link_flags = kwargs['extra_link_args'] if 'extra_link_args' in kwargs else []
        cmd = [self.compiler] + obj_files + ['-o', self.name+'-lib/'+self.executable]\
               + link_flags

        self.run_process(cmd)

    def __call__(self, *args):
        cmd = ['mpirun', '-np', str(self.nmpi), 
               './'+self.name+'-lib/'+self.executable]\
               + [str(arg) for arg in args]

        self.run_process(cmd)

    def run_process(self, cmd):
        proc = sp.run(cmd, shell=False, capture_output=True, universal_newlines=True)
        print(proc.stdout)
