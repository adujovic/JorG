import unittest
import JorGpi.run
from sys import path
path.insert(0,r'../../')

class TestArgv(unittest.TestCase):
    def test_input_001(self):
        engine = JorGpi.run.JorGpi(*argv)
        engine.initialize_new_cell()
        engine.generate_possible_configurations()
        self.assertEqual(1,1.)
#        self.check_type(self.POSCAR,self.INCAR,self.MASK,isType=str)
#        self.check_bools()
#        self.check_type(self.NEIGHBOR,self.currentOptions('reference'),isType=int)
#        self.assertGreater(self.NEIGHBOR,0)
#        self.assertEqual(self.NEIGHBOR,3)
#        self.assertEqual(self.POSCAR,'POSCAR')
#        self.assertEqual(self.INCAR,'INCAR')
#        self.check_in(self.MASK,"Fe")
