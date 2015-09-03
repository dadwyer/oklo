import unittest

from oklo.core.ids import NuclideId, ReactionId
from oklo.core.defs import ReactionType
from oklo.core.reaction import Reaction

class TestReaction(unittest.TestCase):

    def setUp(self):
        self.fixture = Reaction(reac_id=ReactionId(NuclideId('Yttrium_96'),
                                                   ReactionType.BetaDecay))
                        
    def tearDown(self):
        del self.fixture
        
    def test_reac_id(self):
        self.assertEqual(self.fixture.id, ReactionId(NuclideId('Yttrium_96'),
                                                     ReactionType.BetaDecay))
        self.assertNotEqual(self.fixture.id, ReactionId(NuclideId('Yttrium_96'),
                                                    ReactionType.AlphaDecay))

if '__main__'==__name__:
    unittest.main()

