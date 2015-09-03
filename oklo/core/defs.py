##########################################################################
def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    enums['to_name'] = sequential.__getitem__
    return type('Enum', (), enums)

##########################################################################

ReactionType = enum('Unknown',
                    'AlphaDecay',
                    'BetaDecay',
                    'Beta+Decay',
                    'BetaNeutronDecay',
                    'GammaDecay',       # For metastable isomer decays only
                    'NeutronCapture',
                    'p_n',
                    'n_p')

