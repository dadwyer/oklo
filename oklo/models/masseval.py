from oklo.core.ids import NuclideId
from oklo.core.model import NuclideModel
from oklo.utils.parsers import parse_isomers_table
####################################################################

class MassEvaluation(NuclideModel):
    '''Decorate nuclide with mass evaluation data'''
    def __init__(self, **kwargs):
        '''Constructor'''
        NuclideModel.__init__(self, **kwargs)
        # parse standard mass evaluation table
        self._mass_excess_by_id = {} 
        self._parse_table(kwargs['mass_data'])
        self._load_isomers(kwargs['isomer_data'])
        return

    def _parse_table(self, mass_data):
        '''Parse the mass evaluation data file'''
        from oklo.utils.parsers import parse_mass_eval_table
        nuclides_data = parse_mass_eval_table(mass_data)
        for nuclide_data in nuclides_data:
            nucl_id=nuclide_data['id']
            self._mass_excess_by_id[nucl_id] = nuclide_data['mass_excess']
        return

    def _load_isomers(self, isomer_pkg):
        '''Load the isomer data from python module'''
        isomer_data = parse_isomers_table(isomer_pkg)
        for isomer in isomer_data:
            # Estimate mass excess by adding isomer level to ground-state
            isomer_id = NuclideId(Z=isomer['Z'], A=isomer['A'], M=isomer['M'])
            ground_id = NuclideId(Z=isomer['Z'], A=isomer['A'], M=0)
            mass_excess = (self._mass_excess_by_id[ground_id]
                           + isomer['energy_level'])
            self._mass_excess_by_id[isomer_id] = mass_excess
        return
    
    def process(self, nuclide):
        '''Add mass information to this nuclide'''
        if not self._mass_excess_by_id.has_key(nuclide.id):
            print '  No mass data for nuclide %s' % nuclide.id
            return
        nuclide['mass_excess'] = self._mass_excess_by_id[nuclide.id]
        return

    def known_ids(self):
        '''Return the list of ids of nuclides known by this model'''
        return sorted(self._mass_excess_by_id.keys())
    
####################################################################

