#!/usr/bin/env python

# Load core components
from oklo.core.factory import NuclideFactory, ReactionFactory
from oklo.core.network import ReactionNetwork
from oklo.core.ids import NuclideId, ReactionId
from oklo.core.defs import ReactionType
from oklo.core.units import MeV, seconds

# Load physical models and data
from oklo.models.masseval import MassEvaluation
from oklo.models.fissionyield import FissionYieldENDF
from oklo.models.betaspectrum import BetaSpectrumENDF

# Load other tools
from numpy import linspace, zeros, vectorize
from math import exp

#########################################################################
# Example: Setup a reaction network based on a set of models
#
# Prepare models
mass_model_ame = MassEvaluation(name='AtomicMassEval2012',
                                mass_data='data/mass.mas12',
                                isomer_data='data/isomers_nndc.py')

fission_model_endf = FissionYieldENDF(
    name='FissionYieldENDF_v7',
    yield_data=[
        'data/endfb_vii/nfpy_9228_92-U-235.dat',
        'data/endfb_vii/nfpy_9237_92-U-238.dat',
        'data/endfb_vii/nfpy_9437_94-Pu-239.dat',
        'data/endfb_vii/nfpy_9443_94-Pu-241.dat'])

decay_model_endf = BetaSpectrumENDF(
    name='BetaSpectrumENDF_v7',
    decay_data=['data/ensdf/beta_decays_ensdf6_ahayes.txt'])


# Prepare Factories
mass_factory = NuclideFactory(name='MassFactory',
                              model_list=[{'model':mass_model_ame,
                                           'scope':'default'},])

fission_factory = NuclideFactory(name='FissionYieldFactory',
                                 model_list=[{'model':fission_model_endf,
                                              'scope':'default'},])

decay_factory = ReactionFactory(name='BetaSpectrumFactory',
                                model_list=[{'model':decay_model_endf,
                                             'scope':'default'},])

factories = [mass_factory, fission_factory, decay_factory]

# Prepare Network
antinu_network = ReactionNetwork('AntinuSpectrumNetwork')

# Load data into reaction network
for factory in factories:
    factory.process(antinu_network)

#########################################################################
# Examples: Accessing Nuclide and Reaction data from the reaction network
#

# Define energy range and resolution for spectral calculation
energies = linspace(0*MeV, 15*MeV, 1501)

# Access data of one nuclide:
Y_96_id = NuclideId('Y_96')
Y_96 = antinu_network.get(id=Y_96_id)
print ' '
print 'Example nuclide:'
print '  nuclide:',Y_96_id
print '  mass_excess[MeV]:',Y_96['mass_excess'] / MeV
print '  cumulative fission yields:'
fiss_yields = Y_96['cumulative_yield']
fiss_parents = sorted(fiss_yields.keys())
for fiss_parent in fiss_parents:
    print '    %s: %f' % (fiss_parent, fiss_yields[fiss_parent])
    
# Access data of one reaction:
Y_96_beta_decay_id = ReactionId(init_nucl_id=Y_96_id,
                                reac_type=ReactionType.BetaDecay)
Y_96_beta_decay = antinu_network.get(id=Y_96_beta_decay_id)
print ' '
print 'Example reaction:'
print '  reaction:',Y_96_beta_decay_id
beta_decay_Y_96 = Y_96_beta_decay['beta_decay'] 
print '  q_value[MeV]:',beta_decay_Y_96.q_value / MeV
print '  half_life[s]:',beta_decay_Y_96.half_life / seconds
antinuspec_Y_96 = beta_decay_Y_96.antineutrino_spectrum(energies)
#print '  antinu_spec:',antinuspec_Y_96
print ''

#########################################################################
# Calculate total antineutrino spectrum per fission for a nominal reactor
#
#  Step 1: Provide relative fission rates of actinides in reactor
#          (Nominal fractions based on Daya Bay PWR reactors)
fission_fractions = {NuclideId('Uranium_235'):0.584,
                     NuclideId('Uranium_238'):0.076,
                     NuclideId('Plutonium_239'):0.290,
                     NuclideId('Plutonium_241'):0.050}
#
#  Step 2: Prepare total spectrum and other useful values
antinuspec_reactor = zeros(len(energies))
fission_daughters_total = 0
fission_daughters_included = 0
decay_rate_total = 0
decay_rate_included = 0
missing_decays = []
#
#  Step 3: Loop over known nuclides, summing antineutrino spectra
for nuclide in antinu_network.nuclides:
    # Skip nuclides which are not known fission daughters
    if not nuclide.has_key('cumulative_yield'): continue
    # Calculate decay rate for this nuclide in reactor, assuming equilibrium
    decay_rate = 0
    cumulative_yield = nuclide['cumulative_yield'] 
    for fiss_parent in fission_fractions.keys():
        # Skip fission parent if yield is not known
        if not cumulative_yield.has_key(fiss_parent): continue
        # Otherwise, use yield to estimate equilibrium decay rate
        fiss_yield = cumulative_yield[fiss_parent]
        decay_rate += fiss_yield*fission_fractions[fiss_parent]
    # Skip nuclide if decay rate is zero
    if decay_rate==0: continue
    # Prepare Reaction ID for beta decay of this nuclide
    beta_decay_id = ReactionId(init_nucl_id=nuclide.id,
                               reac_type=ReactionType.BetaDecay)
    # Check if beta decay above interaction threshold is allowed
    # according to energy conservation
    if beta_decay_id.final_nuclide_id in antinu_network.known_ids():
        # Retrieve final nuclide
        final_nuclide = antinu_network.get(beta_decay_id.final_nuclide_id)
        # Check that mass is known for initial and final states
        if (nuclide.has_key('mass_excess')
            and final_nuclide.has_key('mass_excess')):
            avail_energy = nuclide['mass_excess'] - final_nuclide['mass_excess']
            # Skip if energy below antineutrino detection threshold
            if avail_energy < 1.8 * MeV: continue
    # Count the number of relevant fission daughters
    fission_daughters_total += 1
    # Count the total fission daughter equilibrium decay rate
    decay_rate_total += decay_rate
    # Check if daughter has known beta decay spectrum
    if beta_decay_id not in antinu_network.known_ids():
        # No reaction information for this possible decay
        missing_decays.append( {'reac_id':beta_decay_id,
                                'decay_rate':decay_rate} )
        continue
    beta_reaction = antinu_network.get(beta_decay_id)
    if not beta_reaction.has_key('beta_decay'):
        # No decay spectrum information for this decay
        missing_decays.append( {'reac_id':beta_decay_id,
                                'decay_rate':decay_rate} )
        continue
    beta_decay = beta_reaction['beta_decay']
    # Add this decay to total reactor spectrum
    antinuspec_reactor += (decay_rate
                           * beta_decay.antineutrino_spectrum(energies))
    # Keep track of the fraction of decays with known spectrum
    fission_daughters_included += 1
    decay_rate_included += decay_rate
missing_decays = sorted(missing_decays, key=lambda elem: elem['decay_rate'],
                        reverse=True)
#
# Step 4: Print calculation results
print ''
print 'Reactor Antineutrino Spectrum Calculation:'
print '  Number of fission daughters in the calculation:'
print '   Total known, above interaction threshold:',fission_daughters_total
print '   With spectral data: %.3f (%.3f)' % (
    fission_daughters_included,
    fission_daughters_included/float(fission_daughters_total))
print '  Daughter nuclide decay rate, per initial actinide fission:'
print '   Total known, above interaction threshold:',decay_rate_total
print '   With spectral data:  %.3f (%.3f)' % (
    decay_rate_included,
    decay_rate_included/decay_rate_total)
print '  Rates for prominent missing decays:'
for missing_decay in missing_decays[:10]:
    print '    %0.3f [%0.2f%%]: %s' % (
        missing_decay['decay_rate'],
        100*(missing_decay['decay_rate']/decay_rate_total),
        missing_decay['reac_id'])
#print '  Energies:',energies
#print '  Spectrum:',antinuspec_reactor
print ''

# Generate plots, if matplotlib available
plot_avail = False
import imp
try:
    imp.find_module('matplotlib')
    plot_avail = True
except ImportError:
    plot_avail = False

def smooth_spec_func(energy):
    x = energy
    a = [4.73926e-01, 3.87703e-01, -3.61888e-01, 4.97195e-02, -2.99050e-03]
    return exp(a[0] + a[1]*x + a[2]*(x**2) + a[3]*(x**3) + a[4]*(x**4))

smooth_spec_vec = vectorize(smooth_spec_func)
smooth_spec = smooth_spec_vec(energies)

if plot_avail:
    import matplotlib.pyplot as plt
    from pylab import figure, axis, grid
    plt.ion()

    figure(1)
    plt.plot(energies,antinuspec_Y_96,'r-')
    axis([0*MeV, 10*MeV, 0, 0.25])
    grid(True)
    plt.xlabel('Antineutrino Energy [MeV]')
    plt.ylabel(r'dN/dE [MeV$^{-1}$]')
    plt.title('Yttrium-96 Spectrum')
    plt.show()

    figure(2)
    plt.semilogy(energies,antinuspec_reactor,'r-')
    plt.semilogy(energies,smooth_spec,'b--')
    axis([2*MeV, 8*MeV, 0.001, 1.5])
    grid(True)
    plt.xlabel('Antineutrino Energy [MeV]')
    plt.ylabel(r'dN/dE [MeV$^{-1}$ fission$^{-1}$]')
    plt.title('Nominal Reactor Spectrum')
    plt.show()

    figure(3)
    plt.plot(energies,antinuspec_reactor/smooth_spec,'r-')
    axis([2*MeV, 8*MeV, 0.85, 1.15])
    grid(True)
    plt.xlabel('Antineutrino Energy [MeV]')
    plt.ylabel('(dN/dE) / F(E)')
    plt.title('Nominal Reactor Spectrum')
    plt.show()
