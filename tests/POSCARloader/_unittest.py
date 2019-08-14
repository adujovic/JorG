from POSCARloader import POSCARloader
import numpy as np
import unittest

class TestPOSCARloader(unittest.TestCase):
    atom_input = "0.0000000000 0.0000000000 4.3870399043 Cu"
    atom_proper = np.array([0.0000000000,0.0000000000,4.3870399043])
    atom_input_constr = "-0.00000000000000 0.00000000000000 0.36270017366243 True False False"
    atom_constr_proper = np.array([-0.00000000000000,0.00000000000000,0.36270017366243])
    directions = [(3.87753641405488, -0.0, -0.0),
                 (-0.0, 3.87753641405488, 0.0),
                 (0.0, 0.0, 12.09549987244948)]
    cell       = [(0.5, 0.5, 0.19398936278918),
                  (0.5, 0.5, 0.80601063721082),
                  (0.0, 0.0, 0.0),
                  (-0.0, 0.0, 0.36270017366243),
                  (0.0, -0.0, 0.63729982633757),
                  (0.0, -0.0, 0.15023942974438),
                  (-0.0, 0.0, 0.84976057025562),
                  (0.0, 0.5, 0.37983766131592),
                  (0.5, 0.0, 0.37983766131592),
                  (-0.0, 0.5, 0.62016233868408),
                  (0.5, -0.0, 0.62016233868408),
                  (0.5, 0.5, 0.5)]
    atoms      = [56, 56, 29, 29, 29, 8, 8, 8, 8, 8, 8, 39]

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
        name = POSCARloader.parse_atom_name(self.atom_input)
        self.assertEqual(name,'Cu')

    def test_parse_name_constrains(self):
        name = POSCARloader.parse_atom_name(self.atom_input_constr)
        self.assertEqual(name,'True')

    def test_parse_constrains(self):
        constr = POSCARloader.parse_constrains(self.atom_input)
        self.assertIsNone(constr)

    def test_parse_constrains_nodata(self):
        constr = POSCARloader.parse_constrains(self.atom_input_constr)
        self.assertEqual(constr,[True, False, False])

    def test_loader_comment(self):
        self.loader.parse()
        self.assertEqual(self.loader(0)['comment'],"Ba2Cu3O6Y1 [TET,TET,tP12] (STD_PRIM do")
        self.assertEqual(self.loader(1)['comment'],"Ba2Cu3O6Y1 [TET,TET,tP12] (STD_PRIM do")
        self.assertEqual(self.loader(2)['comment'],"Ba2Cu3O6Y1 [TET,TET,tP12] (STD_PRIM do")

    def test_loader_directions(self):
        self.loader.parse()
        self.assertAlmostEqual(np.linalg.det(self.loader(0)['directions']),181.8593318589)
        self.assertAlmostEqual(np.linalg.det(self.loader(1)['directions']),181.85933185544175)
        self.assertAlmostEqual(np.linalg.det(self.loader(2)['directions']),181.85933185544175)

    def test_loader_cellVolume(self):
        self.loader.parse()
        self.assertAlmostEqual(self.loader(0)['cellVolume'],181.8593318589)
        self.assertAlmostEqual(self.loader(1)['cellVolume'],181.85933185544175)
        self.assertAlmostEqual(self.loader(2)['cellVolume'],181.85933185544175)

    def test_loader_cellCenter(self):
        centerOne = np.array([0.80988675,0.81115342,5.54377077])
        centerTwo = np.array([0.80782009,0.80782009,5.54377077])
        self.loader.parse()
        self.assertAlmostEqual(np.linalg.norm(self.loader(0)['cellCenter']-centerOne),0.0)
        self.assertAlmostEqual(np.linalg.norm(self.loader(1)['cellCenter']-centerTwo),0.0)
        self.assertAlmostEqual(np.linalg.norm(self.loader(2)['cellCenter']-centerTwo),0.0)

    def test_loader_cellAtoms(self):
        cellAtoms = [2,3,6,1]
        self.loader.parse()
        for i,number in enumerate(cellAtoms):
            self.assertEqual(self.loader(0)['cellAtoms'][i],number)
            self.assertEqual(self.loader(1)['cellAtoms'][i],number)
            self.assertEqual(self.loader(2)['cellAtoms'][i],number)

    def test_loader_atomNames(self):
        atomNamesOne = ['Ba', 'Cu', 'O', 'Y']
        atomNamesTwo = ['H', 'He', 'Li', 'Be']
        self.loader.parse()
        for i,name in enumerate(atomNamesOne):
            self.assertEqual(self.loader(0)['atomNames'][i],name)
            self.assertEqual(self.loader(1)['atomNames'][i],name)
        for i,name in enumerate(atomNamesTwo):
            self.assertEqual(self.loader(2)['atomNames'][i],name)

    def are_equal(self,cellA,cellB):
        for atomA,atomB in zip (cellA,cellB):
            self.assertEqual(atomA[0],atomB[0])
            self.assertAlmostEqual(np.linalg.norm(atomA[1]-atomB[1]),0.0,
                                   msg="%s != %s"%(str(atomA[1]),str(atomB[1])))

    def test_loader_cell_0(self):
        self.loader.parse()
        cell = [('Ba', np.array([1.94876821, 1.93876821,  2.34639831])),
                ('Ba', np.array([1.93856821, 1.93876821,  9.74910156])),
                ('Cu', np.array([0.00500000, 0.00000000,  0.00000000])),
                ('Cu', np.array([0.00000000, 0.00000000,  4.3870399 ])),
                ('Cu', np.array([0.00000000, 0.00000000,  7.70845997])),
                ('O' , np.array([0.00600000, 0.00000000,  1.817221  ])),
                ('O' , np.array([0.00000000, 0.00000000, 10.2782788690999993])),
                ('O' , np.array([0.00400000, 1.93876821,  4.59432638])),
                ('O' , np.array([1.93876821, 0.04000000,  4.59432638])),
                ('O' , np.array([0.00000000, 1.93876821,  7.50117349])),
                ('O' , np.array([1.93876821, 0.00000000,  7.50117349])),
                ('Y' , np.array([1.93876821, 1.93876821,  6.04774994]))]
        self.are_equal(cell,self.loader(0)[ 'cell'])

    def test_loader_cell_1(self):
        self.loader.parse()
        cell = [('Ba', np.array([1.93876821, 1.93876821,  2.34639831])),
                ('Ba', np.array([1.93876821, 1.93876821,  9.74910156])),
                ('Cu', np.array([0.00000000, 0.00000000,  0.00000000])),
                ('Cu', np.array([0.00000000, 0.00000000,  4.3870399])),
                ('Cu', np.array([0.00000000, 0.00000000,  7.70845997])),
                ('O' , np.array([0.00000000, 0.00000000,  1.817221])),
                ('O' , np.array([0.00000000, 0.00000000, 10.27827887])),
                ('O' , np.array([0.00000000, 1.93876821,  4.59432638])),
                ('O' , np.array([1.93876821, 0.00000000,  4.59432638])),
                ('O' , np.array([0.00000000, 1.93876821,  7.50117349])),
                ('O' , np.array([1.93876821, 0.00000000,  7.50117349])),
                ('Y' , np.array([1.93876821, 1.93876821,  6.04774994]))]
        self.are_equal(cell,self.loader(1)[ 'cell'])

    def test_loader_cell_2(self):
        self.loader.parse()
        cell = [(0, np.array([1.93876821, 1.93876821,  2.34639831])),
                (0, np.array([1.93876821, 1.93876821,  9.74910156])),
                (1, np.array([0.00000000, 0.00000000,  0.00000000])),
                (1, np.array([0.00000000, 0.00000000,  4.3870399])),
                (1, np.array([0.00000000, 0.00000000,  7.70845997])),
                (2, np.array([0.00000000, 0.00000000,  1.817221])),
                (2, np.array([0.00000000, 0.00000000, 10.27827887])),
                (2, np.array([0.00000000, 1.93876821,  4.59432638])),
                (2, np.array([1.93876821, 0.00000000,  4.59432638])),
                (2, np.array([0.00000000, 1.93876821,  7.50117349])),
                (2, np.array([1.93876821, 0.00000000,  7.50117349])),
                (3, np.array([1.93876821, 1.93876821,  6.04774994]))]
        self.are_equal(cell,self.loader(2)[ 'cell'])

    def test_loader_cellSymmetry_0(self):
        self.loader.parse()
        cell       = [(0.5025789570959635, 0.4999999999871052, 0.19398936279219897),
                      (0.49994842084544383, 0.4999999999871052, 0.8060106372160686),
                      (0.0012894785415343496, 0.0, 0.0),
                      (0.0, 0.0, 0.3627001736662843),
                      (0.0, 0.0, 0.6372998263419832),
                      (0.0015473742498412196, 0.0, 0.1502394297441653),
                      (0.0, 0.0, 0.8497605702558346),
                      (0.0010315828332274797, 0.4999999999871052, 0.3798376613176211),
                      (0.4999999999871052, 0.010315828332274797, 0.3798376613176211),
                      (0.0, 0.4999999999871052, 0.6201623386906465),
                      (0.4999999999871052, 0.0, 0.6201623386906465),
                      (0.4999999999871052, 0.4999999999871052, 0.5)]
        for directionSet,directionRead in zip(self.directions,
                                              self.loader(0)['cellSymmetry'][0]):
            self.assertAlmostEqual(np.linalg.norm(directionSet),
                                   np.linalg.norm(directionRead))
        for atomSet,atomRead in zip(cell,
                                    self.loader(0)['cellSymmetry'][1]):
            self.assertAlmostEqual(np.linalg.norm(atomSet),
                                   np.linalg.norm(atomRead))
        for nameSet,nameRead in zip(self.atoms,
                                    self.loader(0)['cellSymmetry'][2]):
            self.assertEqual(nameSet,nameRead)

    def test_loader_cellSymmetry_1(self):
        self.loader.parse()
        for directionSet,directionRead in zip(self.directions,
                                              self.loader(1)['cellSymmetry'][0]):
            self.assertAlmostEqual(np.linalg.norm(directionSet),
                                   np.linalg.norm(directionRead))
        for atomSet,atomRead in zip(self.cell,
                                    self.loader(1)['cellSymmetry'][1]):
            self.assertAlmostEqual(np.linalg.norm(atomSet),
                                   np.linalg.norm(atomRead))
        for nameSet,nameRead in zip(self.atoms,
                                    self.loader(1)['cellSymmetry'][2]):
            self.assertEqual(nameSet,nameRead)

    def test_loader_cellSymmetry_2(self):
        self.loader.parse()
        atoms      = [0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3]
        for directionSet,directionRead in zip(self.directions,
                                              self.loader(2)['cellSymmetry'][0]):
            self.assertAlmostEqual(np.linalg.norm(directionSet),
                                   np.linalg.norm(directionRead))
        for atomSet,atomRead in zip(self.cell,
                                    self.loader(2)['cellSymmetry'][1]):
            self.assertAlmostEqual(np.linalg.norm(atomSet),
                                   np.linalg.norm(atomRead))
        for nameSet,nameRead in zip(atoms,
                                    self.loader(2)['cellSymmetry'][2]):
            self.assertEqual(nameSet,nameRead)

if __name__ == '__main__':
    unittest.main()
