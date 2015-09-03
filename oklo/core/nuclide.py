from oklo.core.data import DataContainer
##########################################################################

class Nuclide(DataContainer):
    '''A single nuclide, with unique ID and associated information
    '''
    def __init__(self, nucl_id):
        '''Constructor'''
        DataContainer.__init__(self, nucl_id)
        return

    @property
    def name(self):
        '''Return the name for this nuclide'''
        return str(self.id)
        
    @property
    def Z(self):
        '''Return number of protons (Z) for this Nuclide'''
        return self.id.Z

    @property
    def A(self):
        '''Return number of nucleons (A) for this Nuclide'''
        return self.id.A

    @property
    def M(self):
        '''Return metastable nuclear isomer level (M) for this Nuclide'''
        return self.id.M

    @property
    def N(self):
        '''Return the neutron number (N) for this Nuclide'''
        return self.id.N

    @property
    def element_name(self):
        '''Return element name for this nuclide (e.g. Hydrogen)'''
        return self.id.element_name

    @property
    def element_abbrev(self):
        '''Return element abbreviation for this nuclide (e.g. H)'''
        return self.id.element_abbrev

    def is_nuclide(self):
        '''Allow others to ask if this is a nuclide'''
        return True
