import unittest
from JorGpi.generate import run
from os import environ

class TestJorGpi(unittest.TestCase):
    def test_clear_run(self):
        engine = run.JorGpi()
        engine.initialize_new_cell()
        try:
            if environ['RUNCPP'] == '1':
                engine.possible_configurations()
        except KeyError:
            self.assertEqual(1,1.)
        del engine
        self.assertEqual(1,1.)

    def test_symmetry_run(self):
        engine = run.JorGpi(['test', '--symmetry'])
        engine.initialize_new_cell()
        try:
            if environ['RUNCPP'] == '1':
                engine.possible_configurations()
        except KeyError:
            self.assertEqual(1,1.)
        del engine
        self.assertEqual(1,1.)
