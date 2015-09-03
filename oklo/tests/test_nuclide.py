import unittest

from oklo.core.ids import NuclideId, ReactionId
from oklo.core.nuclide import Nuclide

class TestNuclide(unittest.TestCase):

    def setUp(self):
        self.fixture = Nuclide(NuclideId(Z='Yttrium',A=96,M=0))

    def tearDown(self):
        del self.fixture
        
    def test_nuclide_Z(self):
        self.assertEqual(self.fixture.Z,39)

    def test_nuclide_A(self):
        self.assertEqual(self.fixture.A,96)
        
    def test_nuclide_M(self):
        self.assertEqual(self.fixture.M,0)

    def test_nuclide_N(self):
        self.assertEqual(self.fixture.N,96-39)

    def test_nuclide_element(self):
        self.assertEqual(self.fixture.element_name,'Yttrium')

    def test_nuclide_element_abbr(self):
        self.assertEqual(self.fixture.element_abbrev,'Y')

    def test_nuclide_name(self):
        self.assertEqual(self.fixture.name,'Yttrium_96')

if '__main__'==__name__:
    unittest.main()
