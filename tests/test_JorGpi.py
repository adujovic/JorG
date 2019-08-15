import unittest
import JorGpi.run
from sys import path
path.insert(0,r'../')

class TestArgv(unittest.TestCase):
    def test_input_001(self):
        engine = JorGpi.run.JorGpi()
        engine.initialize_new_cell()
        engine.possible_configurations()
        self.assertEqual(1,1.)
