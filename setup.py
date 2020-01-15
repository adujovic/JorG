import re
from setuptools import setup, Extension
from os import system,makedirs,environ
from sys import argv,path

packages_jorgpi = ['JorGpi',
                   'JorGpi.generate',
                   'JorGpi.pickup',
                   'JorGpi.utilities',
                   'JorGpi.utilities.KPOINTS',
                   'JorGpi.aux',
                   'JorGpi.geometry']
executables_jorgpi = ['JorGpi/bin/JorGpi-kpoints',
                      'JorGpi/bin/JorGpi-pickup',
                      'JorGpi/bin/JorGpi-POSCAR',
                      'JorGpi/bin/JorGpi-startup']

requirements_jorgpi = []

if __name__ == '__main__':
    options=""
    try:
        if argv[1] != 'install':
            exit(-1)
    except IndexError:
        exit(-1)
    for arg in argv[2:]:
        options+=" %s"%arg

    options=re.sub('=',' ',options)
    options=re.sub('~',environ['HOME'],options)
    if '--install_reqs' in argv:
        print('Installing requirements')
        system('./install-requirements.sh'+options)
        argv.remove('--install_reqs')

    try:
        locDir = environ['PYTHONUSERBASE']
        try:
            makedirs(locDir,exist_ok=False)
            ISNEWPATH=True
        except FileExistsError:
            ISNEWPATH=False
        if ISNEWPATH:
            with open(environ['HOME']+'/.bashrc','a+') as bashrc:
                bashrc.write('\nexport PYTHONPATH=$PYTHONPATH:%s'%locDir)
    except KeyError:
        locDir = environ['HOME']+'/.local'
        print('Installing in %s'%(environ['HOME']))

    options+=' --prefix %s'%locDir
    try:
        with open('requirements.txt','r') as reqs:
            for line in reqs.readlines():
                package = re.sub('[><=]+.*','',line)
                package = re.sub('\s','',package)
                if len(package)>0:
                    requirements_jorgpi.append(package)
    except FileNotFoundError:
        requirements_jorgpi = [ 'numpy', 'scipy', 'spglib', 'matplotlib', 'setuptools', 'defusedxml' ]
    setup(name='JorGpi',
          version='0.1',
          description='JorGpi DFT-to-Heisenberg mapping module',
          long_description="""
          JorGpi DFT-to-Heisenberg mapping module
          """,
          author='Andrzej P. Kądzielawa',
          author_email='andrzej.piotr.kadzielawa@vsb.cz',
          url='https://github.com/Mellechowicz/JorG',
          packages=packages_jorgpi,
          install_requires=requirements_jorgpi,
          provides=['JorGpi'],
          scripts=executables_jorgpi,
          platforms=['POSIX'],
          license='',
          test_suite='nose.collector',
          tests_require=['nose'])

    system('./install-gsl.sh'+options)
    system('./install-asa.sh'+options)
