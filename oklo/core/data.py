##########################################################################

class DataContainer(dict):
    '''A dictionary-like container for data associated to a nuclide or reaction
    '''
    def __init__(self, id):
        '''Constructor'''
        dict.__init__(self)
        self._id = id

    @property
    def id(self):
        '''Return the unique ID for this data object'''
        return self._id

    def is_nuclide(self):
        '''Allow others to ask if this is nuclide data'''
        return False

    def is_reaction(self):
        '''Allow others to ask if this is reaction data'''
        return False
