OKLO: A toolkit for modeling nuclides and nuclear reactions
===========================================================

The OKLO package provides a set of convenient tools for modeling
nuclides and the reactions between them.  It provides methods to
uniquely identify the known nuclides, define unique reactions which
allow transitions from one nuclide to another, associate arbitrary
user-specified data to these objects, and treat collections of
nuclides and reactions as a single network.

While this package was initially designed for the purpose of
calculating the flux of antineutrinos emitted by nuclear reactors, it
is designed in a generic fashion for a broad range of applications
(e.g. solar physics, supernova physics, etc.).


Installation:
=============

Installion is most convenient via the pip utility:
::
 $ pip install oklo

Description:
============

Here follows a brief description of the core tools provided by the
OKLO package.

Nuclide:
--------

A nuclide is a unique nuclear state identified by the number of
protons and number of neutrons in the nucleus, as well as an optional
metastable isomeric energy level.  Each nuclide serves as a data
'whiteboard'.  Users can associate arbitrary data with a nuclide using
a standard (key, value) approach.

Example::
 >>> from oklo.core.ids import NuclideId
 >>> from oklo.core.nuclide import Nuclide 
 >>> C_12_id = NuclideId(Z=6,A=12,M=0)      # Create unique ID for Carbon-12
 >>> C_12 = Nuclide(nucl_id=C_12_id)        # Create nuclide object
 >>> C_12.Z
 6
 >>> C_12.A
 12
 >>> C_12.M
 0
 >>> C_12.N
 6
 >>> C_12.name
 'Carbon_12'
 >>> C_12.element_name
 'Carbon'
 >>> C_12.element_abbrev
 'C'
 >>> C_12['current_abundance'] = 0.8   # Associate user-defined data 
 >>> C_12['current_abundance']
 0.8


Reaction:
---------

A reaction is a unique type of nuclear transition between nuclides,
such as the beta decay of 12B to 12C.  Each reaction serves as a data
'whiteboard'.  Users can associate arbitrary data with a reaction
using a standard (key, value) approach.

Example::
 >>> from oklo.core.ids import NuclideId, ReactionId
 >>> from oklo.core.defs import ReactionType
 >>> from oklo.core.reaction import Reaction
 >>> B_12_beta_decay_id = ReactionId(init_nucl_id=NuclideId('Boron_12'), \
                                     reac_type=ReactionType.BetaDecay) 
 >>> B_12_beta_decay = Reaction(reac_id=B_12_beta_decay_id)
 >>> B_12_beta_decay.initial_nuclide_id.name
 Boron_12
 >>> B_12_beta_decay.final_nuclide_id.name
 Carbon_12
 >>> from oklo.core.units import Hz
 >>> B_12_beta_decay['current_rate'] = 1.0 * Hz
 >>> B_12_beta_decay['current_rate'] / Hz
 1.0


ReactionNetwork:
----------------

A ReactionNetwork is a collection of nuclides and reactions relating
these nuclides.  This class serves as a standard entry point for
calculation of quantities of interest.  For example, users can iterate
over the nuclides and reactions within the network to calculate total
quantities for the network.

This class is effectively a 'graph' data structure in computer science
lingo, where a set of nodes (nuclides) are connected by a set of links
(reactions).

Physical Models (NuclideModel and ReactionModel):
-------------------------------------------------

A physical model is a user-defined class which defines rules for
automating the addition of user-specified data to the whiteboards
Nuclides or Reactions.  For example, you could specify the relative
abundance of each nuclide present within the sun according to your
preferred physical model.

Users should define their own physical model class which inherits from
either NuclideModel or ReactionModel base classes, and implements a
'process(nuclide)' or 'process(reaction)' function which adds the
appropriate data to the whiteboard based on the unique nuclide or
reaction ID.

Factories (NuclideFactory and ReactionFactory):
-----------------------------------------------

Factories provide a convenient method to populate a ReactionNetwork
given one or a set of physical models.

Examples:
=========

For a more advance example, demonstrating the use of models and
factories, execute the following:

::
$ python -i oklo/examples/antineutrino_spectrum_endf.py

This builds a reaction network for modeling a nominal commercial PWR
reactor.  The network is populated with tabulated nuclear data on
cumulative fission yields and beta decay spectra.  The network is then
used to estimate the average antineutrino energy spectrum emitted per
fission in the reactor.

If matplotlib is installed, then associated figures will also be
generated.
