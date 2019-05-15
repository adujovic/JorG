#!/usr/bin/python3

from JorG.loadsave import load_POSCAR,save_vanilla_POSCAR

a = load_POSCAR('POSCAR',direct=True)
print(a)
save_vanilla_POSCAR('POSCAR_fixed',a)
