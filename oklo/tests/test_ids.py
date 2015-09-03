import unittest

from oklo.core.ids import NuclideId, ReactionId
from oklo.core.defs import ReactionType

class TestNuclideId(unittest.TestCase):

    def setUp(self):
        self.fixture = NuclideId(Z='Yttrium',A=96,M=0)

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

    def test_nuclide_compare1(self):
        self.assertTrue(self.fixture == NuclideId(Z=39,A=96,M=0))

    def test_nuclide_compare2(self):
        self.assertFalse(self.fixture == NuclideId(Z=39,A=97,M=0))
        
    def test_nuclide_compare3(self):
        self.assertTrue(self.fixture < NuclideId(Z=39,A=97,M=0))

    def test_nuclide_compare4(self):
        self.assertFalse(self.fixture == NuclideId(Z=39,A=96,M=1))

    def test_nuclide_endf(self):
        self.assertEqual(self.fixture.endf_name,'0390960')

    def test_nuclide_element(self):
        self.assertEqual(self.fixture.element_name,'Yttrium')

    def test_nuclide_element_abbr(self):
        self.assertEqual(self.fixture.element_abbrev,'Y')

    def test_nuclide_str(self):
        self.assertEqual(str(self.fixture),'Yttrium_96')

class TestReactionId(unittest.TestCase):

    def setUp(self):
        self.fixture = {'reac1': ReactionId(NuclideId('Yttrium_96'),
                                            ReactionType.BetaDecay),
                        'reac2': ReactionId(NuclideId('Strontium_96'),
                                            ReactionType.BetaDecay),
                        'reac3': ReactionId(NuclideId(Z='Yttrium',A=96),
                                            ReactionType.BetaDecay)}
                        

    def tearDown(self):
        del self.fixture
        
    def test_reac_id(self):
        self.assertEqual(self.fixture['reac1'],self.fixture['reac3'])

    def test_reac_final(self):
        self.assertEqual(self.fixture['reac1'].initial_nuclide_id,
                         self.fixture['reac2'].final_nuclide_id)
        
    def test_reac_str(self):
        self.assertEqual(str(self.fixture['reac1']),
                         'Yttrium_96_BetaDecay_to_Zirconium_96')

    def test_reac_type(self):
        self.assertEqual(self.fixture['reac1'].reaction_type,
                         ReactionType.BetaDecay)

    def test_reac_compare1(self):
        self.assertTrue(self.fixture['reac1'] == self.fixture['reac3'])

    def test_reac_compare2(self):
        self.assertFalse(self.fixture['reac1'] == self.fixture['reac2'])
        
    def test_reac_compare3(self):
        self.assertTrue(self.fixture['reac1'] > self.fixture['reac2'])
        
if '__main__'==__name__:
    unittest.main()

