import unittest
from JorGpi.argv import Options

class TestArgv(unittest.TestCase):
    def set_up(self,_input):
        self.currentOptions = Options(*_input.split())
        self.mask     = str(self.currentOptions('mask'))
        self.output   = self.currentOptions('output')
        self.incar    = self.currentOptions('incar')
        self.poscar   = self.currentOptions('input')
        self.neighbor = self.currentOptions('neighbor')

    def check_type(self,*args,isType):
        for arg in args:
            self.assertIsInstance(arg,isType)

    def check_bools(self):
        self.check_type(self.currentOptions('symmetry'),
                        self.currentOptions('minimal_set'),
                        self.currentOptions('refined'),isType=bool)

    def check_in(self,where,*args):
        for arg in args:
            self.assertIn(arg,where)

    def test_input_001(self):
        _input="foo --incar INCAR --input POSCAR -N 3 -E Fe"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,3)
        self.assertEqual(self.poscar,'POSCAR')
        self.assertEqual(self.incar,'INCAR')
        self.check_in(self.mask,"Fe")

    def test_input_002(self):
        _input="foo --incar INCAR --input POSCAR -N 13 -E Ni"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,13)
        self.assertEqual(self.poscar,'POSCAR')
        self.assertEqual(self.incar,'INCAR')
        self.check_in(self.mask,"Ni")

    def test_input_003(self):
        _input="foo --incar INCAR --input POSCAR -N 3 -E Ni"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,3)
        self.assertEqual(self.poscar,'POSCAR')
        self.assertEqual(self.incar,'INCAR')
        self.check_in(self.mask,"Ni")
        self.assertGreater(self.neighbor,0)

    def test_input_004(self):
        _input="foo --incar INCAR_tst1 --input POSCAR_tst1 -N 3 -E Fe"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,3)
        self.assertEqual(self.poscar,'POSCAR_tst1')
        self.assertEqual(self.incar,'INCAR_tst1')
        self.check_in(self.mask,"Fe")
        self.assertGreater(self.neighbor,0)

    def test_input_005(self):
        _input="foo --symmetry -i POSCAR_Cs2F6Ni2 "
        self.set_up(_input)
        self.check_type(self.poscar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.currentOptions('reference'),isType=int)
        self.assertIsNone(self.neighbor)
        self.assertEqual(self.poscar,'POSCAR_Cs2F6Ni2')
        self.assertEqual(self.incar,None)
        self.assertEqual(self.currentOptions('symmetry'),True)

    def test_input_006(self):
        _input="foo -i POSCAR_tst1 "
        self.set_up(_input)
        self.check_type(self.poscar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.currentOptions('reference'),isType=int)
        self.assertIsNone(self.neighbor)
        self.assertEqual(self.poscar,'POSCAR_tst1')
        self.assertEqual(self.incar,None)
        self.check_in(self.mask,"Mn")

    def test_input_007(self):
        _input="foo --incar INCAR_tst1 --input POSCAR_tst1 -N 2 -E Fe1"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,2)
        self.assertEqual(self.poscar,'POSCAR_tst1')
        self.assertEqual(self.incar,'INCAR_tst1')
        self.check_in(self.mask,"Fe")

    def test_input_008(self):
        _input="foo --incar _INCARs/INCAR_tst1 --input _POSCARs/POSCAR_tst1 -N 1 -E Fe1"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,1)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_tst1')
        self.assertEqual(self.incar,'_INCARs/INCAR_tst1')
        self.check_in(self.mask,"Fe")

    def test_input_009(self):
        _input="foo --incar _INCARs/INCAR_CsNiF --input _POSCARs/POSCAR_CsNiF -N 1 -E Ni -o output/J1"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.output,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,1)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_CsNiF')
        self.assertEqual(self.incar,'_INCARs/INCAR_CsNiF')
        self.assertEqual(self.output,'output/J1')
        self.check_in(self.mask,"Ni")

    def test_input_010(self):
        _input="foo --symmetry -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF "
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.currentOptions('reference'),isType=int)
        self.assertIsNone(self.neighbor)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_CsNiF')
        self.assertEqual(self.incar,'_INCARs/INCAR_CsNiF')
        self.check_in(self.mask,"Fe")
        self.assertEqual(self.currentOptions('symmetry'),True)

    def test_input_011(self):
        _input="foo -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 6 -E Ni -o output/ASD --minimal-set"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.output,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,6)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_CsNiF')
        self.assertEqual(self.incar,'_INCARs/INCAR_CsNiF')
        self.assertEqual(self.output,'output/ASD')
        self.check_in(self.mask,"Ni")
        self.assertEqual(self.currentOptions('minimal_set'),True)

    def test_input_012(self):
        _input="foo -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -E Ni -o output/ASD --reference 21 -R 123.0"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.output,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.currentOptions('reference'),isType=int)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_CsNiF')
        self.assertEqual(self.incar,'_INCARs/INCAR_CsNiF')
        self.assertEqual(self.output,'output/ASD')
        self.check_in(self.mask,"Ni")
        self.assertEqual(self.currentOptions('reference'),21-1)
        self.check_type(self.currentOptions('cutOff'),isType=float)
        self.assertAlmostEqual(self.currentOptions('cutOff'),123.0)

    def test_input_013(self):
        _input="foo -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD --refined"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.output,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,2)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_CsNiF')
        self.assertEqual(self.incar,'_INCARs/INCAR_CsNiF')
        self.assertEqual(self.output,'output/ASD')
        self.check_in(self.mask,"Ni")
        self.assertEqual(self.currentOptions('refined'),True)

    def test_input_014(self):
        _input="foo -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD  --period 3d 5p"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.output,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,2)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_CsNiF')
        self.assertEqual(self.incar,'_INCARs/INCAR_CsNiF')
        self.assertEqual(self.output,'output/ASD')
        self.check_in(self.mask,"Ni","Ti","Mn","Te","In")

    def test_input_015(self):
        _input="foo -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD  --block D"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.output,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,2)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_CsNiF')
        self.assertEqual(self.incar,'_INCARs/INCAR_CsNiF')
        self.assertEqual(self.output,'output/ASD')
        self.check_in(self.mask,"Ni","Fe","Au","Hg","Cr","Hf")

    def test_input_016(self):
        _input="foo -i _POSCARs/POSCAR_CsNiF -I _INCARs/INCAR_CsNiF -N 2 -E Ni -o output/ASD  --group 2 5"
        self.set_up(_input)
        self.check_type(self.poscar,self.incar,self.output,self.mask,isType=str)
        self.check_bools()
        self.check_type(self.neighbor,self.currentOptions('reference'),isType=int)
        self.assertGreater(self.neighbor,0)
        self.assertEqual(self.neighbor,2)
        self.assertEqual(self.poscar,'_POSCARs/POSCAR_CsNiF')
        self.assertEqual(self.incar,'_INCARs/INCAR_CsNiF')
        self.assertEqual(self.output,'output/ASD')
        self.check_in(self.mask,"Ni","Be","Mg","Ca","Nb","Ta")

if __name__ == '__main__':
    unittest.main()
