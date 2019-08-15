import unittest
from pickup.pickup import SmartPickUp,Reference,CommandLineOptions

class TestPickupIron(unittest.TestCase):
    def options(self,*args):
        return CommandLineOptions(*args)

    def test_iron_001(self):
        _input = "test -R _VASP/Fe/noFlip -D _VASP/Fe/flip00000 -E Fe -J1 -U mRy".split(" ") 
        options = self.options(*_input)
        elements = ''.join(options('elements'))
        self.assertEqual(elements,'Fe$')

        ref = Reference(options('reference')+"/POSCAR")
        self.assertEqual(ref(),0)
        self.assertEqual(options('number_of_interactions'),1)

        pickerUpper = SmartPickUp(options('number_of_interactions'),elements)
        pickerUpper.read(options('reference'),*options('directories'),reference=ref())

        self.assertEqual(options('units'),'mRy')
        Js = pickerUpper.solve(units=options('units')).flatten()
        self.assertEqual(Js[0],0.8724046188693274)
        self.assertEqual(Js[1],4.157645364906014)
