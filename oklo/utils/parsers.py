from oklo.core.ids import NuclideId
from oklo.core.units import seconds, keV, eV
##########################################################################

def parse_mass_eval_table(filename):
    '''Parse the Atomic Mass Evaluation table, and generate appropriate
    data for each nuclide.'''
    import pkg_resources
    if not pkg_resources.resource_exists('oklo',filename):
        raise ValueError('Mass evaluation file "%s" does not exist' % (
            filename))
    datafile = pkg_resources.resource_stream('oklo',filename)
    datalines = datafile.readlines()
    # Skip header lines
    n_header_lines = 39
    datalines = datalines[n_header_lines:]
    # Parse each line, and generate nuclide data
    nuclides_data = []
    for line in datalines:
        nucl_Z = int(line[9:14])
        nucl_A = int(line[14:19])
        mass_excess = float(line[27:41].strip().strip('#')) * keV
        #beta_decay_energy = None
        #bde_string = line[75:86].strip().strip('#').strip('*')
        #if len(bde_string)>0:
        #    beta_decay_energy = float(bde_string) * keV
        # Collate data
        nuclide_data = {'id' : NuclideId(Z=nucl_Z,A=nucl_A),
                        'mass_excess' : mass_excess}
        nuclides_data.append(nuclide_data)
    return nuclides_data

##########################################################################

def parse_isomers_table(filename):
    '''Parse the Isomers table, and generate appropriate data for each
    nuclide.'''
    import pkg_resources
    if not pkg_resources.resource_exists('oklo',filename):
        raise ValueError('Isomer data file "%s" does not exist' % (
            filename))
    full_filename = pkg_resources.resource_filename('oklo',filename)
    vars = {}
    execfile(full_filename, vars)
    return vars['isomer_data']

##########################################################################

def parse_yields_ENDFB(filenames):
    '''A function to parse the ENDF/B original fission yield files. It
    takes a list of ENDF filenames, loads their contents, and returns
    a dictionary of the cumulative fission yields by parent.'''
    import os.path
    yields_by_parent = {}
    for filename in filenames:
        import pkg_resources
        if not pkg_resources.resource_exists('oklo',filename):
            raise ValueError('Fission yield file "%s" does not exist' % (
                filename))
        datafile = pkg_resources.resource_stream('oklo',filename)
        fileraw = datafile.readlines()
        nuclLine = fileraw[5]
        nZ = int(nuclLine[1:3])
        nElem = (nuclLine[4:6]).strip()
        nA = int(nuclLine[7:10])
        fissParentId = NuclideId(Z=nZ, A=nA)
        # Parse cumulative yield data
        dataArray = []
        nLines = len(fileraw)
        for lineIdx in range(nLines):
            line = fileraw[lineIdx]
            line = line.strip()
            if line.endswith('8459    1'):
                # Found cumulative yield data block
                while True:
                    lineIdx += 1
                    line = fileraw[lineIdx]
                    line = line.strip()
                    if line.startswith('5.000000+5') or lineIdx==nLines:
                        #if line.startswith('2.530000-2') or lineIdx==nLines:
                        # Found block with reactor neutron energy distribution
                        lineIdx += 1
                        break
                if lineIdx >= nLines: break  # End of data
                # Should be at reactor section
                for lIdx in range(lineIdx, nLines):
                    line = fileraw[lIdx]
                    if line.startswith(' 1.400000+7'):
                        #if line.startswith(' 5.000000+5'):
                        # Found start of 14MeV data block.  Done.
                        break
                    data = line[1:66]
                    #print data
                    dataToks = data.split()
                    for dataTok in dataToks:
                        dataArray.append(dataTok)
            else:
                lineIdx += 1
        # Process data array
        dIdx = 0
        nData = len(dataArray)
        yields = {}
        while dIdx < nData-2:
            isotope = int(round( convertENDFField(dataArray[dIdx])*10 ))
            isomer = int(round( convertENDFField(dataArray[dIdx+1]) ))
            isotope += isomer
            if isotope != 0 and isotope >= 220660 and isotope <=721720:
                yld = convertENDFField( dataArray[dIdx+2] )
                sigYld = convertENDFField( dataArray[dIdx+3] )
                daught_id = NuclideId(endf_id=isotope)
                yields[daught_id] = {'cumulative': yld,
                                     'cumulative_unc': sigYld}
            dIdx += 4
        yields_by_parent[fissParentId] = yields
    # Return final data
    return yields_by_parent

def convertENDFField( data ):
    # Convert ENDF formatted number to actual value
    dataBase = 0.0
    dataExp = 0
    if data.find('+')>0:
        dataTok = data.split('+')
        dataBase = float(dataTok[0])
        dataExp = int(dataTok[1])
    elif data.find('-')>0:
        dataTok = data.split('-')
        dataBase = float(dataTok[0])
        dataExp = -1*int(dataTok[1])
    else:
        dataBase = float(data)
    return dataBase * pow(10,dataExp)

##########################################################################

def parse_decays_ENDF(filenames):
    '''A function to parse the ENDF beta decay data files. It
    takes a list of ENDF filenames, loads their contents, and returns
    a list of decay information blocks.'''
    import os.path
    data_lines = []
    for filename in filenames:
        import pkg_resources
        if not pkg_resources.resource_exists('oklo',filename):
            raise ValueError('Decay data file "%s" does not exist' % (
                filename))
        datafile = pkg_resources.resource_stream('oklo',filename)
        data_lines += datafile.readlines()
    decay_infos = []
    dataLen = len(data_lines)
    lineCounter = 0
    while True:
        try:
            # Process data block from one decay
            line = data_lines[lineCounter]
            # Skip comment lines
            if (line.strip())[0] == '#':
                lineCounter += 1
                continue
            numE = int(line.split()[4])
            data = data_lines[lineCounter:lineCounter+numE+1]
            nucInfo = data[0].split()
            decay_info = {
                'Z':int(nucInfo[0]),
                'A':int(nucInfo[1]),
                'half_life':float(nucInfo[2]) * seconds,
                'Q':float(nucInfo[3]) * eV,
                'E0max':float(nucInfo[5]) * eV,
                'M':int(nucInfo[6]),
                'branches':None,
            }
            # Process each beta decay branch
            branches = []
            for line in data[1:]:
                # Skip comment lines
                if (line.strip())[0] == '#': continue
                branch_data = line.split()
                branch = {
                    'e0':float(branch_data[0]) * eV,
                    'sigma_e0':float(branch_data[1]) * eV,
                    'fraction':float(branch_data[2]),
                    'sigma_fraction':float(branch_data[3]),
                    'forbiddeness':int(float(branch_data[4])),
                }
                branches.append(branch)
            decay_info['branch_infos'] = branches
            decay_infos.append(decay_info)
            lineCounter += 1 + numE
        except:
            break
    #print '  parse_decays_ENDF: Processed %d decays.' % (len(decay_infos))
    return decay_infos
