##########################################################################

class BaseModel(object):
    '''Base class for all models'''
    def __init__(self, **kwargs):
        '''Constructor'''
        self._name = kwargs['name']
        return

    @property
    def name(self):
        '''Return the name of this model'''
        return self._name
    
    def known_ids(self):
        '''Return list of IDs known to this model'''
        return []
    
##########################################################################    

class NuclideModel(BaseModel):
    '''Base class for all nuclide models'''
    def __init__(self, **kwargs):
        '''Constructor'''
        BaseModel.__init__(self, **kwargs)
        return

    def process(self, nuclide):
        '''Add data to nuclide based on this model'''
        raise ValueError('This is an abstract method.  Must be redefined!')

##########################################################################    
    
class ReactionModel(BaseModel):
    '''Base class for all reaction models'''
    def __init__(self, **kwargs):
        '''Constructor'''
        BaseModel.__init__(self, **kwargs)
        return

    def process(self, reaction):
        '''Add data to reaction based on this model'''
        raise ValueError('This is an abstract method.  Must be redefined!')
    
