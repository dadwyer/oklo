from oklo.core.model import NuclideModel
from oklo.utils.parsers import parse_yields_ENDFB
##########################################################################

class FissionYieldENDF(NuclideModel):
    '''Add cumulative fission yield data to nuclides'''
    def __init__(self, **kwargs):
        '''Constructor'''
        # Parse configuration and load data
        NuclideModel.__init__(self, **kwargs)
        self._yields_by_id = {}
        self._load_yields(kwargs['yield_data'])
        return

    def known_ids(self):
        '''Return the IDs of reactions known to this model'''
        return self._yields_by_id.keys()
    
    def process(self, nuclide):
        '''Add cumulative fission yields to this nuclide'''
        if not self._yields_by_id.has_key(nuclide.id):
            #print '  No known yield for nuclide %s' % nuclide.id
            return
        nuclide['cumulative_yield'] = self._yields_by_id[nuclide.id]
        return
        
    def _load_yields(self,yield_data):
        '''Process the configuration for this model'''
        # Parse fission yield data
        yields_by_parent = parse_yields_ENDFB(yield_data)
        # Reformat as a table for each fission daughter
        #  Step 1: Determine list of yielded daughter nuclides
        fiss_daughters = set()
        for yields in yields_by_parent.values():
            daughter_ids = yields.keys()
            fiss_daughters.update( set(daughter_ids) )
        #  Step 2: Prepare yield table for each daughter nuclide
        yields_by_id = {}
        for daught_id in fiss_daughters:
            daughter_yields = {}
            for (parent_id, yields) in yields_by_parent.iteritems():
                if yields.has_key(daught_id):
                    fiss_yield = yields[daught_id]['cumulative']
                    daughter_yields[parent_id] = fiss_yield
            yields_by_id[daught_id] = daughter_yields
        # Save reformatted yield table
        self._yields_by_id = yields_by_id
        return

##########################################################################

