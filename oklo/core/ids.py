from oklo.core.elements import element_name_table, element_Z_table
from oklo.core.defs import ReactionType
##########################################################################

class NuclideId(object):
    '''Represents a unique nuclide identifier, specified by:
       - Proton number, Z
       - Nucleon number, A
       - Isomer index, M (zero for ground state)
    '''
    def __init__(self, name=None, Z=None, A=None, M=None, id=None,
                 endf_id=None):
        '''Constructor.  Specify either name, (z,n,m), or ID'''
        if id:
            # Define by ID
            self._id = id
        elif None not in (Z, A):
            # Define by (Z,A,(M))
            if not M: M=0
            self._id = NuclideId._zam_to_id(Z,A,M)
        elif name:
            self._id = NuclideId._name_to_id(name)
        elif endf_id:
            self._id = NuclideId._endf_to_id(endf_id)
        else:
            raise ValueError('Invalid nuclide definition: Z=%r, A=%r, M=%r, id=%r' % (Z,A,M,id))
        self._element_name = element_name_table[self.Z]
        return
        
    @classmethod
    def _id_to_zam(cls, id):
        '''Convert Nuclide ID integer to (proton, nucleon, isomer) notation'''
        return (int(id / 10000), int((id % 10000)/10), id % 10)

    @classmethod
    def _zam_to_id(cls, Z, A, M):
        '''Convert (proton, nucleon, isomer) notation to Nuclide ID integer'''
        if isinstance(Z,str):
            # Convert element name or abbreviation to proton number Z
            Z = element_Z_table[Z]
        return int(Z*10000 + A*10 + M)

    @classmethod
    def _name_to_id(cls, name):
        '''Convert name 'Elem_A_mM' notation to Nuclide ID integer'''
        name_parts = name.split('_')
        Z = element_Z_table[name_parts[0]]
        A = int(name_parts[1])
        M = 0
        if len(name_parts)>2:
            M = int(name_parts[2].strip('m'))
        return cls._zam_to_id(Z,A,M)
    
    @classmethod
    def _endf_to_id(cls, endf):
        '''Convert ENDF (ZZZAAAM) notation to Nuclide ID integer'''
        return int(endf)

    #@property
    #def id(self):
    #    '''Return the Nuclide ID as integer'''
    #    return self._id

    @property
    def Z(self):
        '''Return number of protons (Z) for this Nuclide ID'''
        return NuclideId._id_to_zam(self._id)[0]

    @property
    def A(self):
        '''Return number of nucleons (A) for this Nuclide ID'''
        return NuclideId._id_to_zam(self._id)[1]

    @property
    def M(self):
        '''Return metastable nuclear isomer level (M) for this Nuclide ID'''
        return NuclideId._id_to_zam(self._id)[2]

    @property
    def N(self):
        '''Return the neutron number (N) for this Nuclide ID'''
        return self.A - self.Z

    def __hash__(self):
        '''Define hash comparison based on nuclide ID'''
        return self._id

    def __cmp__(self, other):
        '''Compare two Nuclide IDs'''
        return self._id.__cmp__(other._id)

    def __str__(self):
        '''Return string representation of nuclide ID'''
        name_str = self.element_name
        name_str += '_%d' % self.A
        if self.M > 0:
            name_str += '_m%d' % self.M
        return name_str

    @property
    def name(self):
        '''Return string representation (Element_AAA) of nuclide ID'''
        return str(self)
    
    @property
    def endf_name(self):
        '''Return ENDF string representation (ZZZAAAM) of nuclide ID'''
        return '%07d' % self._id

    @property
    def element_name(self):
        '''Return element name for this nuclide ID (e.g. Hydrogen)'''
        return self._element_name[0]

    @property
    def element_abbrev(self):
        '''Return element abbreviation for this nuclide ID (e.g. H)'''
        return self._element_name[1]

##########################################################################

class ReactionId(object):
    '''Unique key for identifying a specific transition between nuclides'''
    def __init__(self, init_nucl_id=None, reac_type=None, final_nucl_id=None):
        '''Constructor.  Specify reaction using:
        - initial nuclide name or ID, reaction_type, and final nuclide name/ID
        If final nuclide not provided, a reasonable guess will be made.
        (eg. 3H beta decay will result in 2H)
        '''
        # Convert nuclide name to ID, if needed
        if isinstance(init_nucl_id,str):
            init_nucl_id = NuclideId(init_nucl_id)
        # Convert reaction type name to ReactionType, if needed
        if isinstance(reac_type,str):
            reac_type = getattr(ReactionType,reac_type)
        # Convert nuclide name to ID, if needed
        if isinstance(final_nucl_id,str):
            final_nucl_id = NuclideId(final_nucl_id)
        if not final_nucl_id:
            # Automatically determine final nuclide based on reaction type
            final_nucl_id = ReactionId._determine_final_nuclide(init_nucl_id,
                                                                reac_type)
        self._init_nucl_id = init_nucl_id
        self._reac_type = reac_type
        self._final_nucl_id = final_nucl_id
        # Create a single, unique integer ID for this reaction
        self._id = None
        self._id = self.__hash__()
        return

    @classmethod
    def _determine_final_nuclide(cls, init_nucl_id, reac_type):
        '''Determine the most natural final nuclide for this reaction'''
        final_nucl_id = None
        [Zin, Ain] = [init_nucl_id.Z,init_nucl_id.A]
        if reac_type is ReactionType.AlphaDecay:
            final_nucl_id = NuclideId(Z=Zin-2,A=Ain-4,M=0)
        elif reac_type is ReactionType.BetaDecay:
            final_nucl_id = NuclideId(Z=Zin+1,A=Ain,M=0)
        elif reac_type is ReactionType.Beta+Decay:
            final_nucl_id = NuclideId(Z=Zin-1,A=Ain,M=0)
        elif reac_type is ReactionType.GammaDecay:
            final_nucl_id = NuclideId(Z=Zin,A=Ain,M=0)
        elif reac_type is ReactionType.NeutronCapture:
            final_nucl_id = NuclideId(Z=Zin,A=Ain+1,M=0)
        elif reac_type is ReactionType.p_n:
            final_nucl_id = NuclideId(Z=Zin+1,A=Ain,M=0)
        elif reac_type is ReactionType.n_p:
            final_nucl_id = NuclideId(Z=Zin-1,A=Ain,M=0)
        else:
            raise ValueError('Cannot determine final nuclide: init=%r, final=%r'
                             % (init_nucl_id, final_nucl_id))
        return final_nucl_id

    def __hash__(self):
        '''Define a unique hash integer for this reaction ID'''
        if not self._id:
            # Calculate once and store hash
            hash_id = ((self._init_nucl_id.__hash__() * 10000)
                       + self._reac_type * 100
                       + self._final_nucl_id.M)
            self._id = hash_id
        return self._id
                
    def __cmp__(self,other):
        '''Allow comparison of two reaction IDs'''
        return self._id.__cmp__(other._id)
    
    def __str__(self):
        '''Return string representation of this reaction ID'''
        return '%s_%s_to_%s' % (str(self._init_nucl_id),
                                ReactionType.to_name(self._reac_type),
                                str(self._final_nucl_id))

    @property
    def initial_nuclide_id(self):
        '''Return the initial state nuclide ID'''
        return self._init_nucl_id

    @property
    def final_nuclide_id(self):
        '''Return the final state nuclide ID'''
        return self._final_nucl_id

    @property
    def reaction_type(self):
        '''Return the final state nuclide ID'''
        return self._reac_type

    
