def get_combined_data(opt, period, inst, mods=None):
    if mods is None:
        mods = [opt]
    casenames = [m['casename'] for m in mods]

    analysis_files = []
    for mod in mods:
        _, filename = pkg_utils.get_process_path(mod, inst, period)
        analysis_files += [filename]

    # TODO: Harmonize this.  Waterlevel data is stored differently and this causes headaches.
    if inst in ['SST', 'MCTD']:
        data, stations = pkg_plot.get_analysis_data(opt, analysis_files, casenames)
        obs = None
    elif inst in ['TG', 'HTG']:
        # Load all pickles
        def loadpickle(mod):
            c4file = os.path.join(mod['dir_process'], inst, period, inst + '_class4_' + mod['casename'] + '.pickle')
            print('Loading ', c4file)
            with open(c4file, 'rb') as fid:
                return pickle.load(fid)

        data = [loadpickle(mod) for mod in mods]

        # Get station names common to all models under comparison (set intersection)
        def getstations(c4):
            return [x['obs']['station'] for x in c4]

        # stations = list(set.intersection(*[set(getstations(c4)) for c4 in data]))
        stationlists = [getstations(c4) for c4 in data]
        stations = list(pkg_utils.get_common_stations(opt, casenames, stationlists))

        # data is two level list of config, station.  Need to be dictionaries rather than lists.
        wl_data = {}

        def getsingle(d, key, stationlist):
            dd = {}
            for stn in d:
                code = stn[key]['station']
                if code in stationlist:
                    dd[code] = stn[key]
                else:
                    dd[code] = None
            return dd

        first = True
        for d, cfg in zip(data, mods):

            obs = getsingle(d, 'obs', stations)
            dd = getsingle(d, 'mod', stations)

            wl_data[cfg['casename']] = dd
            if first:
                wl_data['obs'] = obs
                first = False
            else:
                wl_data['obs'].update(obs)

        data = wl_data
        # TODO: this is repetitive, condense it in the routines higher up
        obs = wl_data['obs']

    elif inst in ['ADCP', 'CM', 'HADCP']:
        # data goes model, station, level, data.  data has both model and obs components.
        # Clunky, but need to harmonize how data is saved across instruments, and can consolidate then.
        data = {}
        for c in casenames:
            data[c] = {}
        stationlists = []
        obs = {}  # pull this out so there's a way to step through the levels relatively straightforwardly
        for c, f in zip(casenames, analysis_files):
            with open(f, 'rb') as fid:
                newdata = pickle.load(fid)
                for k in newdata.keys():  # instruments
                    # if k not in data[c].keys():
                    if newdata[k] is not None:
                        data[c][k] = newdata[k]  # newdata[k] is level, then data (mod and obs)
                        # print (data[c][k].keys())
                        # print (c, k)
                        if k not in obs.keys():
                            obs[k] = {}
                        for z in data[c][k].keys():
                            obs[k][z] = {'total': data[c][k][z]['total']['obs']}
                            if data[c][k][z]['residual'] is not None:
                                obs[k][z]['residual'] = data[c][k][z]['residual']['obs']
                            else:
                                obs[k][z]['residual'] = None

            stationlists.append(list(data[c].keys()))

        # Determine which stations are available in all sources (obs, reference model and candidate models)
        # Now with union!
        stations = pkg_utils.get_common_stations(opt, casenames, stationlists)

    return data, stations, obs
