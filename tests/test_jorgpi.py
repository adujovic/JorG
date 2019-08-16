import unittest
from sys import path
path.insert(0,r'../')
import JorGpi.run

class TestJorGpi(unittest.TestCase):
    def test_clear_run(self):
        engine = JorGpi.run.JorGpi()
        engine.initialize_new_cell()
        engine.possible_configurations(verbose='quiet')
        self.assertEqual(1,1.)

    def test_symmetry_run(self):
        engine = JorGpi.run.JorGpi(['test', '--symmetry'])
        engine.initialize_new_cell()
        engine.possible_configurations(verbose='quiet')
        self.assertEqual(1,1.)
