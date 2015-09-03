from oklo.core.data import DataContainer
##########################################################################

class Reaction(DataContainer):
    '''A single reaction, with unique ID and associated information
    '''
    def __init__(self, reac_id=None):
        DataContainer.__init__(self, reac_id)
        return

    @property
    def initial_nuclide_id(self):
        '''Return the initial state nuclide ID'''
        return self.id._init_nucl_id

    @property
    def final_nuclide_id(self):
        '''Return the final state nuclide ID'''
        return self.id._final_nucl_id

    @property
    def reaction_type(self):
        '''Return the final state nuclide ID'''
        return self.id._reac_type

    def is_reaction(self):
        '''Allow others to ask if this is a reaction'''
        return True

    
##########################################################################
