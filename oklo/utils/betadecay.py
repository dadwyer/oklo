from oklo.core.units import fm, hbarc, alphaFS, mp, me, MeV
from math import pi, log, sqrt, exp, atan, gamma
from numpy import zeros, array_equal
##########################################################################

def complexGammaSqApproxA(x,y):
    '''Approximation to |Gamma(x+iy)|^2, 
    Taken from Wilkinson, NIMA275, 378 (1989)'''
    return exp((x-1/2.)*log(x*x + y*y) 
               - 2*y*atan(y/x) 
               - 2*x 
               + log(2*pi) 
               + (1/6.)*x/(x*x + y*y))

def fermiG_Huber(Z, A, Te):
    '''Fermi function expression taken from P. Huber, arXiv:1106.0687 '''
    if Te < 0.001*MeV: Te=0.001 * MeV # Avoid divergence
    Ee = Te + me
    pe = sqrt(Ee*Ee - me*me)
    gam = sqrt( 1 - (alphaFS*Z)**2 )
    Rn = 0.0029*pow(A,1/3.) + 0.0063*pow(A,-1/3.) - 0.017/A
    gamNumSq = complexGammaSqApproxA(gam,alphaFS*Z*Ee/pe)
    gamDenomSq = (gamma(2*gam+1))**2
    fermiF = (2*(gam+1)*pow(2*(pe/me)*Rn,2*(gam-1))*exp(pi*alphaFS*Z*Ee/pe)
              *gamNumSq/gamDenomSq)
    return (pe/Ee) * fermiF

def alphaFermi(Z, Te):
    '''Alpha parameter for estimation of Fermi Function, taken from
    "A Simple Approximation of the Fermi Function in Nuclear Beta Decay"
    G.K.Schenter and P.Vogel, 1982'''
    if Te < 1.2*me:
        return -0.811 + 4.46e-2*Z + 1.08e-4*Z*Z
    return -8.46e-2 + 2.48e-2*Z + 2.37e-4*Z*Z

def betaFermi(Z, Te):
    '''Beta parameter for estimation of Fermi Function, taken from
    "A Simple Approximation of the Fermi Function in Nuclear Beta Decay"
    G.K.Schenter and P.Vogel, 1982'''
    if Te < 1.2*me:
        return 0.673 - 1.82e-2*Z + 6.38e-5*Z*Z
    return 1.15e-2 + 3.58e-4*Z - 6.17e-5*Z*Z

def fermiG_Vogel(Z, Te):
    ''' Approximation for the Fermi-related G-Function, taken from
    "A Simple Approximation of the Fermi Function in Nuclear Beta Decay"
    G.K.Schenter and P.Vogel, 1982'''
    import math
    Ee = Te + me
    return math.exp(alphaFermi(Z,Te) + betaFermi(Z,Te)*math.sqrt((Ee/me) - 1))

##########################################################################

class BetaDecayBranch:
    '''Calculate beta spectrum for one decay branch.  (Warning: currently
       only handles electrons, not positrons)
    '''
    def __init__(self, reaction_id, e0, sigma_e0, fraction, sigma_fraction,
                 decay_type=None):
        '''Constructor'''
        self._reaction_id = reaction_id
        self._e0 = e0
        self._sigma_e0 = sigma_e0
        self._fraction = fraction
        self._sigma_fraction = sigma_fraction
        # Default: Assume Allowed Gamow-Teller decay type.  Change if requested.
        self._decay_type = 'AllowedGT'
        if decay_type:
            self._decay_type = decay_type
        # Placeholders for cached spectrum evaluation
        self._cached_antinu_energies = None
        self._cached_antinu_spec = None
        # Pre-calculate some convenience variables
        self._A = self._reaction_id.initial_nuclide_id.A
        self._Zdaughter = self._reaction_id.final_nuclide_id.Z
        self._Tmax = self._e0
        self._deltaRadCoeff = (alphaFS/(2*pi))
        self._deltaRadBase = (3*log(mp / (2*(self._Tmax+me))) 
                              + (23/4.) 
                              - (4*pi*pi/3) )
        self._deltaRadBaseE = self._deltaRadBase - (23+3)/4. + (4-2)*pi*pi/3.
        rmoment = (36/35.)*(1.2*(self._A**(1/3.)))*fm
        self._deltaFSCoeff = -(3/2.)*self._Zdaughter*alphaFS*(rmoment / hbarc)
        mu_v = 4.7
        mn = mp
        gA = 1.27590 # arXiv:1007.3790
        self._deltaWMCoeff = (mu_v-(1/2.))/(mn*gA)
        # Set appropriate shape corrections
        if decay_type is "AllowedGT":
            self.weakMagnetismCorrectionNu = self.weakMagnetismCorrectionNu_AllowedGT
        elif decay_type is "NUForbGT_0m":
            self.shapeFactorNu = self.shapeFactorNu_NUForbGT_0m
        elif decay_type is "NUForbGT_1m":
            self.weakMagnetismCorrectionNu = self.weakMagnetismCorrectionNu_NUForbGT_1m
            self.shapeFactorNu = self.shapeFactorNu_NUForbGT_1m
        elif decay_type is "UForbGT_2m":
            self.weakMagnetismCorrectionNu = self.weakMagnetismCorrectionNu_UForbGT_2m
            self.shapeFactorNu = self.shapeFactorNu_UForbGT_2m
        elif decay_type is "NUForbF_1m":
            self.shapeFactorNu = self.shapeFactorNu_NUForbF_1m
        return

    @property
    def reaction_id(self):
        '''Return the reaction id for this beta decay branch'''
        return self._reaction_id
    
    @property
    def e0(self):
        '''Return maxmimum (endpoint) energy of this beta decay branch'''
        return self._e0
    
    @property
    def sigma_e0(self):
        '''Return the 1sig uncertainty in the maxmimum energy'''
        return self._sigma_e0

    @property
    def fraction(self):
        '''Return the branching fraction of this beta decay branch'''
        return self._fraction

    @property
    def sigma_fraction(self):
        '''Return the 1sig uncertainty in the branching fraction'''
        return self._sigma_fraction

    @property
    def decay_type(self):
        '''Return the allowed/forbidden decay type of this beta decay branch'''
        return self._decay_type
    
    def antineutrino_spectrum(self, energies):
        '''Return antineutrino spectrum evaluated at the given energies'''
        if not array_equal(self._cached_antinu_energies, energies):
            # Must recalculate spectrum
            self._cached_antinu_spec = None
        if self._cached_antinu_spec is not None:
            return self._cached_antinu_spec
        # Calculate and cache spectrum
        spectrum = zeros(len(energies))
        for idx in xrange(len(energies)):
            spectrum[idx] = self.dNdE_neutrino(energies[idx])
        # FIXME!!!: Normalization assumes equal spacing, and that
        # energies array spans complete beta decay spectrum.  Enforce
        # or improve!!!
        if energies[0] != 0:
            raise(ValueError, 'Beta decay calculation currently requires array spanning full spectral range.')
        if energies[-1] < self._e0:
            print 'Warning: Beta decay calculation currently assumes array spanning full spectral range.'
        norm = spectrum.sum() * (energies[1]-energies[0])
        if norm != 0:
            spectrum /= norm
        self._cached_antinu_spec = spectrum
        self._cached_antinu_energies = energies.copy()
        return self._cached_antinu_spec
    
    def dNdE_electron_base(self, Te):
        '''Simple allowed e- beta decay shape, with no corrections (Note
           change in phase space factor to match convention of fermiG)
        '''
        Ee = Te + me
        return ((self._Tmax-Te) * (self._Tmax-Te) * Ee * Ee
                * fermiG_Huber(self._Zdaughter,self._A,Te))

    def dNdE_neutrino_base(self, Tnu):
        '''Simple allowed antineutrino beta decay shape, with no corrections'''
        Te = self._Tmax - Tnu
        return self.dNdE_electron_base(Te)

    def dNdE_electron(self, Te):
        '''Complete electron spectrum, with corrections, unnormalized'''
        if (Te<0 or Te > self._Tmax): return 0
        Tnu = self._Tmax - Te
        deltaRad = self.radiativeCorrectionE(Te)
        deltaFS = self.finiteSizeCorrectionNu(Tnu)
        deltaWM = self.weakMagnetismCorrectionNu(Tnu)
        shapeFactor = self.shapeFactorNu(Tnu)
        totalSpec = (self.dNdE_electron_base(Te)
                     *shapeFactor
                     *(1+deltaRad+deltaFS+deltaWM))
        if totalSpec<0: totalSpec = 0
        return totalSpec

    def dNdE_neutrino(self, Tnu):
        '''Complete neutrino spectrum, with corrections, unnormalized'''
        if (Tnu<0 or Tnu >= self._Tmax): return 0
        deltaRad = self.radiativeCorrectionNu(Tnu)
        deltaFS = self.finiteSizeCorrectionNu(Tnu)
        deltaWM = self.weakMagnetismCorrectionNu(Tnu)
        shapeFactor = self.shapeFactorNu(Tnu)
        totalSpec = (self.dNdE_neutrino_base(Tnu)
                     *shapeFactor
                     *(1+deltaRad+deltaFS+deltaWM))
        if totalSpec<0: totalSpec = 0
        return totalSpec

    def radiativeCorrectionNu(self, Tnu):
        '''Calculate the radiative correction at this energy
           Taken from: A Sirlin, PRD84, 014021 (2011)
            Simple version: Eq. 11
           FIXME: confirm validity near antineutrino endpoint
        '''
        y = Tnu / (self._Tmax + me)
        deltaRad = self._deltaRadCoeff*(self._deltaRadBase - 3*log(1-y))
        return deltaRad

    def radiativeCorrectionE(self, Te):
        '''Calculate the radiative correction at this energy
           Taken from: A Sirlin, PRD84, 014021 (2011)
             Approximation taken from Eq. 13
          FIXME: confirm validity for low-energy electrons
        '''
        Eo = self._Tmax + me
        Ee = Te + me
        x = Ee / Eo
        lnx = log(x)
        xbar = (1-x)/x
        if xbar<0.001: xbar = 0.001 # Avoid singularity
        xbarSq = xbar*xbar
        lnxbar = log( xbar )
        lnEm = log(2*Eo/me)
        deltaRad = self._deltaRadCoeff*(self._deltaRadBaseE
                                        + (4*(lnx - 1)
                                           *((1-x)/(3*x) 
                                             - (3/2.) 
                                             + lnxbar))
                                        + lnx*xbarSq/6.
                                        + lnEm*(4*xbar/3.
                                                - 3
                                                + xbarSq/6.
                                                + 4*lnxbar))
        return deltaRad

    def finiteSizeCorrectionNu(self, Tnu):
        '''Calculate the final state correction at this energy
        Taken from A Hayes et. al., arXiv:1309.4146
        '''
        Te = self._Tmax - Tnu
        Ee = Te + me
        deltaFS = self._deltaFSCoeff * (Ee - (Tnu/27.) + ((me*me)/(3*Ee)))
        return deltaFS

    def weakMagnetismCorrectionNu_AllowedGT(self, Tnu):
        '''Weak magnetism correction for allowed GT transitions
        Taken from A Hayes et. al., arXiv:1309.4146
        '''
        Te = self._Tmax - Tnu
        Ee = Te + me
        pe = sqrt(Te*Te + 2*Te*me)
        betaE = pe/Ee
        deltaWM = (2/3.)*self._deltaWMCoeff*(Ee*betaE*betaE-Tnu)
        return deltaWM

    def weakMagnetismCorrectionNu_NUForbGT_1m(self, Tnu):
        '''Weak magnetism correction for non-unique forb. GT 1- transitions
        Taken from A Hayes et. al., arXiv:1309.4146
        '''
        Te = self._Tmax - Tnu
        Ee = Te + me
        EeSq = Ee**2
        peSq = EeSq - me**2
        betaESq = peSq/EeSq
        TnuSq = Tnu*Tnu
        kinFactor = ((((peSq+TnuSq)
                       *(betaESq*Ee-Tnu))
                      +2*betaESq*Ee*Tnu*(Tnu-Ee)/3.)
                     /(peSq+TnuSq-4*betaESq*Tnu*Ee/3.))
        deltaWM = self._deltaWMCoeff*kinFactor
        return deltaWM

    def weakMagnetismCorrectionNu_UForbGT_2m(self, Tnu):
        '''Weak magnetism correction for unique forb. GT 2- transitions
        Taken from A Hayes et. al., arXiv:1309.4146
        '''
        Te = self._Tmax - Tnu
        Ee = Te + me
        EeSq = Ee**2
        peSq = EeSq - me**2
        betaESq = peSq/EeSq
        TnuSq = Tnu*Tnu
        kinFactor = ((((peSq+TnuSq)
                       *(betaESq*Ee-Tnu))
                      +2*betaESq*Ee*Tnu*(Tnu-Ee)/3.)
                     /(peSq+TnuSq))
        deltaWM = (3/5.)*self._deltaWMCoeff*kinFactor
        return deltaWM

    def weakMagnetismCorrectionNu(self, Tnu):
        '''Calculate the weak magnetism correction at this energy'''
        return 0

    def shapeFactorNu_NUForbGT_0m(self, Tnu):
        '''Shape Factor for non-unique forbidden GT 0- transitions
        Taken from A Hayes et. al., arXiv:1309.4146
        '''
        Te = self._Tmax - Tnu
        Ee = Te + me
        EeSq = Ee**2
        peSq = EeSq - me**2
        betaESq = peSq/EeSq
        TnuSq = Tnu*Tnu
        shapeFactor = (peSq + TnuSq + 2*betaESq*Tnu*Ee)
        return shapeFactor

    def shapeFactorNu_NUForbGT_1m(self, Tnu):
        '''Shape Factor for non-unique forbidden GT 1- transitions
        Taken from A Hayes et. al., arXiv:1309.4146
        '''
        Te = self._Tmax - Tnu
        Ee = Te + me
        EeSq = Ee**2
        peSq = EeSq - me**2
        betaESq = peSq/EeSq
        TnuSq = Tnu*Tnu
        shapeFactor = (peSq + TnuSq - (4/3.)*betaESq*Tnu*Ee)
        return shapeFactor

    def shapeFactorNu_UForbGT_2m(self, Tnu):
        '''Shape Factor for unique forbidden GT 2- transitions
        Taken from A Hayes et. al., arXiv:1309.4146
        '''
        Te = self._Tmax - Tnu
        Ee = Te + me
        EeSq = Ee**2
        peSq = EeSq - me**2
        TnuSq = Tnu*Tnu
        shapeFactor = (peSq + TnuSq)
        return shapeFactor

    def shapeFactorNu_NUForbF_1m(self, Tnu):
        '''Shape Factor for non-unique forbidden Fermi 1- transitions
        Taken from A Hayes et. al., arXiv:1309.4146
        '''
        Te = self._Tmax - Tnu
        Ee = Te + me
        EeSq = Ee**2
        peSq = EeSq - me**2
        betaESq = peSq/EeSq
        TnuSq = Tnu*Tnu
        shapeFactor = (peSq + TnuSq + (2/3.)*betaESq*Tnu*Ee)
        return shapeFactor

    def shapeFactorNu(self, Tnu):
        '''Calculate the allowed/forbidden shape factor at this energy'''
        return 1

##########################################################################
    
class BetaDecaySpectrum:
    '''Calculate total beta spectrum for all branches.
    '''
    def __init__(self, reaction_id, q_value, half_life, branches):
        '''Constructor'''
        self._reaction_id = reaction_id
        self._q_value = q_value
        self._half_life = half_life
        self._branches = branches
        self._cached_energies = None
        self._cached_antinu_spec = None
        
    @property
    def reaction_id(self):
        '''Return the reaction id for this beta decay'''
        return self._reaction_id
    
    @property
    def q_value(self):
        '''Return maxmimum (endpoint) energy of this beta decay'''
        return self._q_value

    @property
    def half_life(self):
        '''Return the half life of this beta decay'''
        return self._half_life

    def branches(self):
        '''Return the branches for beta decay'''
        return self._branches
    
    def antineutrino_spectrum(self, energies):
        '''Return antineutrino spectrum evaluated at the given energies'''
        if not array_equal(self._cached_energies, energies):
            # Must recalculate spectrum
            self._cached_antinu_spec = None
        if self._cached_antinu_spec is not None:
            return self._cached_antinu_spec
        # Calculate and cache spectrum
        spectrum = zeros(len(energies))
        for branch in self._branches:
            spectrum += (branch.fraction *
                         branch.antineutrino_spectrum(energies))
        self._cached_antinu_spec = spectrum
        self._cached_energies = energies.copy()
        return spectrum    
