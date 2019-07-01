import unittest
import numpy as np
from JorGpi.generator import NearestNeighborsGenerator

class TestNearestNeighborsGenerator(unittest.TestCase):
    neighbour = 51
    size      = 686
    cutOff    = 17.04
    def setUp(self):
        self.cell            = [['Mn', np.array([0.0, 0.0, 0.0])],
                                ['Mn', np.array([1.3, 1.3, 1.3])]]
        self.referenceAtom   = [0, np.array([0., 0., 0.])]
        self.directions      = [np.array([2.6, 0.0, 1.0]),
                                np.array([1.0, 2.6, 0.0]),
                                np.array([0.0, 1.0, 2.6])]
        self.atomNames       = ['Mn']
        self.generator = NearestNeighborsGenerator(self.cell,
                                          self.referenceAtom,
                                             self.directions)
        self.generator.wyckoffs         = 'abcdefghijklmnopqrstuvwxyz'
        self.generator.atomTypeMask     = '$Mn$'
        self.generator.moments          = [5.0, 5.0]
        self.generator.extraMultiplier  = [0, 0, 0]

    def tearDown(self):
        del self.generator
        del self.cell
        del self.referenceAtom
        del self.directions
        del self.atomNames

    def test_generator_output(self):
        _,crystal,_,_,_,_ = self.generator(self.neighbour)
        self.assertEqual(len(crystal), self.size)

    def test_generator_cutOff(self):
        cutOff,_,_,_,_,_ = self.generator(self.neighbour)
        self.assertIsInstance(cutOff,np.float)
        self.assertAlmostEqual(cutOff,self.cutOff)

    def test_generator_symmetry(self):
        _,_,symmetry,_,_,_ = self.generator(self.neighbour)
        self.assertEqual(len(symmetry['equivalent_atoms']), self.size)

    def test_generator_new_reference(self):
        _,crystal,_,referenceNew,_,_ = self.generator(self.neighbour)
        self.assertEqual(np.linalg.norm(self.referenceAtom[1]-crystal[referenceNew][1]),0.0)

if __name__ == '__main__':
    unittest.main()
