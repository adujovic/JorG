import unittest

from JorGpi.pickup.pickup import SmartPickUp,Reference,CommandLineOptions

class TestPickupIron(unittest.TestCase):
    @staticmethod
    def options(*args):
        return CommandLineOptions(*args)

    def test_iron_001(self):
        _input = "test -R _VASP/Fe/noFlip -D _VASP/Fe/flip00000 -E Fe -J1 -U mRy".split(" ")
        options = TestPickupIron.options(*_input)
        elements = ''.join(options('elements'))
        self.assertEqual(elements,'Fe$')

        ref = Reference(options('reference')+"/POSCAR")
        self.assertEqual(ref(),0)
        self.assertEqual(options('number_of_interactions'),1)

        pickerUpper = SmartPickUp(options('number_of_interactions'),elements)
        pickerUpper.read(options('reference'),*options('directories'),reference=ref())

        self.assertEqual(options('units'),'mRy')
        _J_ij = pickerUpper.solve(units=options('units')).flatten()
        self.assertEqual(_J_ij[0],1.1861042008301703)
        self.assertEqual(_J_ij[1],4.157645364906014)
