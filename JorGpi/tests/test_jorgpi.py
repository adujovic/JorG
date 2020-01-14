import unittest
from JorGpi.generate import run

class TestJorGpi(unittest.TestCase):
    def test_clear_run(self):
        engine = run.JorGpi()
        engine.initialize_new_cell()
        engine.possible_configurations()
        del engine
        self.assertEqual(1,1.)

    def test_symmetry_run(self):
        engine = run.JorGpi(['test', '--symmetry'])
        engine.initialize_new_cell()
        engine.possible_configurations()
        del engine
        self.assertEqual(1,1.)
