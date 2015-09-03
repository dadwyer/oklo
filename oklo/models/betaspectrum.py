from oklo.core.ids import NuclideId, ReactionId
from oklo.core.defs import ReactionType
from oklo.core.model import ReactionModel
from oklo.utils.parsers import parse_decays_ENDF
from oklo.utils.betadecay import BetaDecaySpectrum, BetaDecayBranch
##########################################################################

class BetaSpectrumENDF(ReactionModel):
    '''Add beta decay information and spectra to reactions'''
    def __init__(self, **kwargs):
        '''Constructor'''
        # Parse configuration and load data
        ReactionModel.__init__(self, **kwargs)
        self._decays_by_id = {}
        self._load_decays(kwargs['decay_data'])
        return

    def known_ids(self):
        '''Return the IDs of reactions known to this model'''
        return self._decays_by_id.keys()
    
    def process(self, reaction):
        '''Add beta decay information to this reaction'''
        if not self._decays_by_id.has_key(reaction.id):
            print '  No decay information for reaction %s' % reaction.id
            return
        reaction['beta_decay'] = self._decays_by_id[reaction.id]
        return
        
    def _load_decays(self,decay_data):
        '''Process the configuration for this model'''
        # Parse decay data
        decay_infos = parse_decays_ENDF(decay_data)
        # Reformat as a table by reaction id
        decays_by_id = {}
        for decay_info in decay_infos:
            # Step 1: Construct reaction ID
            reac_id = ReactionId(init_nucl_id=NuclideId(Z=decay_info['Z'],
                                                        A=decay_info['A'],
                                                        M=decay_info['M']),
                                 reac_type=ReactionType.BetaDecay)
            # Step 2: Create beta spectrum calculator for each branch
            branches = []
            for branch_info in decay_info['branch_infos']:
                branch = BetaDecayBranch(reaction_id=reac_id,
                                         e0=branch_info['e0'],
                                         sigma_e0=branch_info['sigma_e0'],
                                         fraction=branch_info['fraction'],
                                         sigma_fraction=branch_info['sigma_fraction'])
                branches.append(branch)
            decay = BetaDecaySpectrum(reaction_id=reac_id,
                                      q_value=decay_info['E0max'],
                                      half_life=decay_info['half_life'],
                                      branches=branches)
            # Step 3: Append to list of decays
            decays_by_id[reac_id] = decay
        # Save reformatted decay table
        self._decays_by_id = decays_by_id
        return

##########################################################################

