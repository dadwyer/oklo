##########################################################################

class ReactionNetwork(object):
    '''Represent a reaction network as a set of nuclides and reactions
    between them.  (Essentially a 'graph data structure' in computer
    science lingo.)
    '''
    def __init__(self, name='Unknown', elements=[]):
        '''Constructor'''
        self._name = name
        self._elements = elements
        self._nuclides = []
        self._reactions = []
        self._reactions_from = {}
        self._reactions_to = {}
        self._id_map = {}
        self._update_index()
        return

    @property
    def name(self):
        return self._name

    @property
    def nuclides(self):
        return self._nuclides

    @property
    def reactions(self):
        return self._reactions

    def get(self, id):
        '''Retrieve nuclide or reaction by ID'''
        return self._id_map[id]

    def known_ids(self):
        '''Retrieve list of known nuclide and reaction ids'''
        return self._id_map.keys()
    
    def reactions_from(self, nucl_id):
        '''Return list of reactions leaving from this nuclide'''
        if not self._reactions_from.has_key(nucl_id):
            return None
        return self._reactions_from[nucl_id]

    def reactions_to(self, nucl_id):
        '''Return list of reactions going to this nuclide'''
        if not self._reactions_to.has_key(nucl_id):
            return None
        return self._reactions_to[nucl_id]

    def network_from(self, id):
        '''Return the sub-network extending from this nuclide or reaction'''
        # FIXME: Implement this function
        return None

    def network_to(self, id):
        '''Return the sub-network extending to this nuclide or reaction'''
        # FIXME: Implement this function
        return None

    def add(self, element_list):
        '''Add list of elements to this network'''
        if len(element_list) < 1: return
        for element in element_list:
            if element.id in self.known_ids():
                print 'Warning: Attempting to re-add %s to network.' % (
                    element.id)
                continue
            self._elements.append(element)
        self._update_index()
    
    def _update_index(self):
        '''Update indices of reactions and nuclides'''
        # Initialize indices
        self._nuclides = []
        self._reactions = []
        self._reactions_from = {}
        self._reactions_to = {}
        self._id_map = {}
        # Sort elements
        self._elements = sorted(self._elements, key=lambda elem: elem.id)
        # Add each element to indices
        for element in self._elements:
            if element.is_nuclide():
                self._nuclides.append(element)
            elif element.is_reaction():
                self._reactions.append(element)
            else:
                raise ValueError('Attempting to add invalid type "%s".' % (
                    element.__class__.__name__))
        for nuclide in self._nuclides:
            self._id_map[nuclide.id] = nuclide
        for reaction in self._reactions:
            self._id_map[reaction.id] = reaction
            init_nucl_id = reaction.initial_nuclide_id
            final_nucl_id = reaction.final_nuclide_id
            if not self._reactions_from.has_key(init_nucl_id):
                self._reactions_from[init_nucl_id] = []
            if not self._reactions_to.has_key(final_nucl_id):
                self._reactions_to[final_nucl_id] = []
            self._reactions_from[init_nucl_id].append(reaction)
            self._reactions_to[final_nucl_id].append(reaction)
        return

##########################################################################
