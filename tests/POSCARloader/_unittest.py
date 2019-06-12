import re
from JorG.loadsave import POSCARloader
import numpy as np
import unittest

class TestPOSCARloader(unittest.TestCase):
    atom_input = "0.0000000000 0.0000000000 4.3870399043 Cu"
    atom_proper = np.array([0.0000000000,0.0000000000,4.3870399043])
    atom_input_constr = "-0.00000000000000 0.00000000000000 0.36270017366243 True False False"
    atom_constr_proper = np.array([-0.00000000000000,0.00000000000000,0.36270017366243])

    def setUp(self):
        self.loader   = POSCARloader('tests/POSCARloader/POSCAR_exp1','tests/POSCARloader/POSCAR_exp2','tests/POSCARloader/POSCAR_exp3')

    def tearDown(self):
        del self.loader

    def test_parse_atom(self):
        atom = POSCARloader.parse_atom(self.atom_input)
        self.assertAlmostEqual(np.linalg.norm(atom-self.atom_proper),0.0)

    def test_parse_atom_constrains(self):
        atom = POSCARloader.parse_atom(self.atom_input_constr)
        self.assertAlmostEqual(np.linalg.norm(atom-self.atom_constr_proper),0.0)

    def test_parse_name(self):
        name = POSCARloader.parse_atomName(self.atom_input)
        self.assertEqual(name,'Cu')

    def test_parse_name_constrains(self):
        name = POSCARloader.parse_atomName(self.atom_input_constr)
        self.assertEqual(name,'True')

    def test_parse_constrains(self):
        constr = POSCARloader.parse_constrains(self.atom_input)
        self.assertIsNone(constr)

    def test_parse_constrains_nodata(self):
        constr = POSCARloader.parse_constrains(self.atom_input_constr)
        self.assertEqual(constr,[True, False, False])

    def test_parse_loader(self):
        self.loader.parse()

if __name__ == '__main__':
    unittest.main()

#Test: parsing POSCAR files
#File "POSCAR_noExistent" not found!

#    print("Test: parsing POSCAR files")
#    for data,i in product(('comment','directions','cell',
#                           'cellSymmetry','cellVolume',
#                           'cellCenter','cellAtoms','atomNames'),range(4)):
#        print("Printing %s of %d file:"%(data,i))
#        try:
#            print(loader(i)[data])
#        except TypeError:
#            print("No data in %d"%i)
#    tracker += time.time()
#    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
#Printing comment of 0 file:
#Ba2Cu3O6Y1 [TET,TET,tP12] (STD_PRIM do
#Printing comment of 1 file:
#Ba2Cu3O6Y1 [TET,TET,tP12] (STD_PRIM do
#Printing comment of 2 file:
#Ba2Cu3O6Y1 [TET,TET,tP12] (STD_PRIM do
#Printing comment of 3 file:
#Run parse first!
#No data in 3
#Printing directions of 0 file:
#[array([ 3.87753641, -0.        ,  0.        ]), array([-0.        ,  3.87753641,  0.        ]), array([-0.        ,  0.        , 12.09549987])]
#Printing directions of 1 file:
#[array([ 3.87753641, -0.        , -0.        ]), array([-0.        ,  3.87753641,  0.        ]), array([ 0.        ,  0.        , 12.09549987])]
#Printing directions of 2 file:
#[array([ 3.87753641, -0.        , -0.        ]), array([-0.        ,  3.87753641,  0.        ]), array([ 0.        ,  0.        , 12.09549987])]
#Printing directions of 3 file:
#Run parse first!
#No data in 3
#Printing cell of 0 file:
#[('Ba', array([1.94876821, 1.93876821, 2.34639831])), ('Ba', array([1.93856821, 1.93876821, 9.74910156])), ('Cu', array([0.005, 0.   , 0.   ])), ('Cu', array([0.       , 0.       , 4.3870399])), ('Cu', array([0.        , 0.        , 7.70845997])), ('O', array([0.006   , 0.      , 1.817221])), ('O', array([ 0.        ,  0.        , 10.27827887])), ('O', array([4.00000000e-03, 1.93876821e+00, 4.59432638e+00])), ('O', array([1.93876821, 0.04      , 4.59432638])), ('O', array([0.        , 1.93876821, 7.50117349])), ('O', array([1.93876821, 0.        , 7.50117349])), ('Y', array([1.93876821, 1.93876821, 6.04774994]))]
#Printing cell of 1 file:
#[('Ba', array([1.93876821, 1.93876821, 2.34639831])), ('Ba', array([1.93876821, 1.93876821, 9.74910156])), ('Cu', array([0., 0., 0.])), ('Cu', array([0.       , 0.       , 4.3870399])), ('Cu', array([0.        , 0.        , 7.70845997])), ('O', array([0.      , 0.      , 1.817221])), ('O', array([ 0.        ,  0.        , 10.27827887])), ('O', array([0.        , 1.93876821, 4.59432638])), ('O', array([1.93876821, 0.        , 4.59432638])), ('O', array([0.        , 1.93876821, 7.50117349])), ('O', array([1.93876821, 0.        , 7.50117349])), ('Y', array([1.93876821, 1.93876821, 6.04774994]))]
#Printing cell of 2 file:
#[(0, array([1.93876821, 1.93876821, 2.34639831])), (0, array([1.93876821, 1.93876821, 9.74910156])), (1, array([0., 0., 0.])), (1, array([0.       , 0.       , 4.3870399])), (1, array([0.        , 0.        , 7.70845997])), (2, array([0.      , 0.      , 1.817221])), (2, array([ 0.        ,  0.        , 10.27827887])), (2, array([0.        , 1.93876821, 4.59432638])), (2, array([1.93876821, 0.        , 4.59432638])), (2, array([0.        , 1.93876821, 7.50117349])), (2, array([1.93876821, 0.        , 7.50117349])), (3, array([1.93876821, 1.93876821, 6.04774994]))]
#Printing cell of 3 file:
#Run parse first!
#No data in 3
#Printing cellSymmetry of 0 file:
#([(3.8775364141, -0.0, 0.0), (-0.0, 3.8775364141, 0.0), (-0.0, 0.0, 12.0954998724)], [(0.5025789570959635, 0.4999999999871052, 0.19398936279219897), (0.49994842084544383, 0.4999999999871052, 0.8060106372160686), (0.0012894785415343496, 0.0, 0.0), (0.0, 0.0, 0.3627001736662843), (0.0, 0.0, 0.6372998263419832), (0.0015473742498412196, 0.0, 0.1502394297441653), (0.0, 0.0, 0.8497605702558346), (0.0010315828332274797, 0.4999999999871052, 0.3798376613176211), (0.4999999999871052, 0.010315828332274797, 0.3798376613176211), (0.0, 0.4999999999871052, 0.6201623386906465), (0.4999999999871052, 0.0, 0.6201623386906465), (0.4999999999871052, 0.4999999999871052, 0.5)], [56, 56, 29, 29, 29, 8, 8, 8, 8, 8, 8, 39])
#Printing cellSymmetry of 1 file:
#([(3.87753641405488, -0.0, -0.0), (-0.0, 3.87753641405488, 0.0), (0.0, 0.0, 12.09549987244948)], [(0.5, 0.5, 0.19398936278918), (0.5, 0.5, 0.80601063721082), (0.0, 0.0, 0.0), (-0.0, 0.0, 0.36270017366243), (0.0, -0.0, 0.63729982633757), (0.0, -0.0, 0.15023942974438), (-0.0, 0.0, 0.84976057025562), (0.0, 0.5, 0.37983766131592), (0.5, 0.0, 0.37983766131592), (-0.0, 0.5, 0.62016233868408), (0.5, -0.0, 0.62016233868408), (0.5, 0.5, 0.5)], [56, 56, 29, 29, 29, 8, 8, 8, 8, 8, 8, 39])
#Printing cellSymmetry of 2 file:
#([(3.87753641405488, -0.0, -0.0), (-0.0, 3.87753641405488, 0.0), (0.0, 0.0, 12.09549987244948)], [(0.5, 0.5, 0.19398936278918), (0.5, 0.5, 0.80601063721082), (0.0, 0.0, 0.0), (-0.0, 0.0, 0.36270017366243), (0.0, -0.0, 0.63729982633757), (0.0, -0.0, 0.15023942974438), (-0.0, 0.0, 0.84976057025562), (0.0, 0.5, 0.37983766131592), (0.5, 0.0, 0.37983766131592), (-0.0, 0.5, 0.62016233868408), (0.5, -0.0, 0.62016233868408), (0.5, 0.5, 0.5)], [0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3])
#Printing cellSymmetry of 3 file:
#Run parse first!
#No data in 3
#Printing cellVolume of 0 file:
#181.85933185893
#Printing cellVolume of 1 file:
#181.85933185544175
#Printing cellVolume of 2 file:
#181.85933185544175
#Printing cellVolume of 3 file:
#Run parse first!
#No data in 3
#Printing cellCenter of 0 file:
#[0.80988675 0.81115342 5.54377077]
#Printing cellCenter of 1 file:
#[0.80782009 0.80782009 5.54377077]
#Printing cellCenter of 2 file:
#[0.80782009 0.80782009 5.54377077]
#Printing cellCenter of 3 file:
#Run parse first!
#No data in 3
#Printing cellAtoms of 0 file:
#[2 3 6 1]
#Printing cellAtoms of 1 file:
#[2 3 6 1]
#Printing cellAtoms of 2 file:
#[2 3 6 1]
#Printing cellAtoms of 3 file:
#Run parse first!
#No data in 3
#Printing atomNames of 0 file:
#['Ba', 'Cu', 'O', 'Y']
#Printing atomNames of 1 file:
#['Ba', 'Cu', 'O', 'Y']
#Printing atomNames of 2 file:
#['H', 'He', 'Li', 'Be']
#Printing atomNames of 3 file:
#Run parse first!
#No data in 3
#
#
