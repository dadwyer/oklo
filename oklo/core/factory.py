from oklo.core.ids import NuclideId, ReactionId
from oklo.core.nuclide import Nuclide
from oklo.core.reaction import Reaction
##########################################################################    

class Factory(object):
    '''Generic factory for building networks based on set of models'''
    def __init__(self, name='Unknown', model_list=[], extend_network=True):
        '''Constructor'''
        self._name = name
        self._default_model = None
        self._model_by_id = {}
        self._known_ids = None
        self._load_models(model_list)
        self._extend_network = extend_network
        return

    @property
    def name(self):
        return self._name
    
    def known_ids(self):
        '''Get the set of ids known to this factory'''
        return self._known_ids
    
    def process_element(self, element):
        '''Process a nuclide or reaction using the current configuration'''
        model = self.get_model(element.id)
        model.process(element)
        return

    def get_model(self, id):
        '''Get the model configured for handling the object with this id'''
        if self._model_by_id.has_key(id):
            return self._model_by_id[id]
        elif self._default_model:
            return self._default_model
        raise ValueError('No model for processing: %r' % id)

    def process(self, network):
        '''Process a reaction network with this factory'''
        # Extend network with new nuclides or reactions, if requested
        if self._extend_network:
            new_elements = []
            for id in self._known_ids:
                if id not in network.known_ids():
                    element = self._make_element(id)
                    new_elements.append(element)
            network.add(new_elements)
        # Process each element in network
        for element in self._get_elements(network):
            self.process_element(element)
        return
    
    def _load_models(self, model_list):
        '''Process model list.  Specifies mapping of ids to models'''
        # Build association between id and model
        known_ids = set()
        for model_info in model_list:
            # Lookup model from current active configuration
            model = model_info['model']
            # Set model scope
            scope = model_info['scope']
            if scope is 'default':
                self._default_model = model
                known_ids.update(model.known_ids())
            else:
                for name in scope:
                    id = _make_id(name)
                    self._model_by_id[id] = model
                    known_ids.add(id)
        self._known_ids = known_ids
        return

    def _make_element(self, id):
        '''Abstract function.  See NuclideFactory or ReactionFactory.'''
        raise ValueError('Abstract Function.  Should not be called.')

    def _get_elements(self, network):
        '''Abstract function.  See NuclideFactory or ReactionFactory.'''
        raise ValueError('Abstract Function.  Should not be called.')

    
##########################################################################

class NuclideFactory(Factory):
    '''A generic factory for generating nuclide data'''
    _make_id = NuclideId

    def _make_element(self, id):
        return Nuclide(nucl_id=id)

    def _get_elements(self, network):
        return network.nuclides
    
class ReactionFactory(Factory):
    '''A generic factory for generating reaction data'''
    _make_id = ReactionId

    def _make_element(self, id):
        return Reaction(reac_id=id)

    def _get_elements(self, network):
        return network.reactions

####################################################################
