from simulations import corr21cm

# TODO: MOVE ME ELSEWHERE -- THEORY LINE
#    zfilename = datapath_db.fetch('sim_15hr', intend_read=True,
#                                 pick='1')
#    zspace_cube = algebra.make_vect(algebra.load(zfilename))
#    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
#    bin_left, bin_center, bin_right = binning.bin_edges(bins, log=True)
#    pwrspec_input = simobj.get_pwrspec(bin_center)
#    if unitless:
#        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

#    agg_array = np.zeros((bin_center.size, len(filename[0])))
#    counter = 0
#    for spec_output in results:
#        (bin_left, bin_center, bin_right, counts_histo, binavg, pwrspec_3d) = spec_output
#        agg_array[:,counter] = binavg
#        counter += 1
#        #for specdata in zip(bin_left, bin_center,
#        #                    bin_right, counts_histo, binavg,
#        #                    pwrspec_input):
#        #    print ("%10.15g " * 6) % specdata
#
#    meanbin = np.mean(agg_array, axis=1)
#    stdbin = np.std(agg_array, axis=1)
#
#    for specdata in zip(bin_left, bin_center,
#                        bin_right, counts_histo, meanbin, stdbin,
#                        pwrspec_input):
#        print ("%10.15g " * 7) % specdata

def test_with_simulation(unitless=True):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()
    pfilename = datapath_db.fetch('simideal_15hr_physical', intend_read=True,
                                 pick='1')
    zfilename = datapath_db.fetch('simideal_15hr', intend_read=True,
                                 pick='1')
    print pfilename
    ofilename = "./physical_cube.npy"
    obfilename = "./physical_cube_beam.npy"
    pwindowfile = "./physical_window.npy"
    owindowfile = "./observed_window_physreg.npy"

    nbins=40
    bins = np.logspace(math.log10(0.00702349679605685),
                       math.log10(2.81187396154818),
                       num=(nbins + 1), endpoint=True)

    pwrspec_2d_p, pwrspec_2d_p, pwrspec_1d_p = \
                    pe.calculate_xspec_file(pfilename, pfilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)

    pwrspec_2d_o, pwrspec_2d_o, pwrspec_1d_o = \
                    pe.calculate_xspec_file(ofilename, ofilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)

    pwrspec_2d_ob, pwrspec_2d_ob, pwrspec_1d_ob = \
                    pe.calculate_xspec_file(obfilename, obfilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)

    zspace_cube = algebra.make_vect(algebra.load(zfilename))
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    # TODO: update this all to the new format
    #for specdata in zip(bin_left, bin_center,
    #                    bin_right, pcounts_histo, obinavg, obbinavg, pbinavg,
    #                    pwrspec_input):
    #    print ("%10.15g " * 8) % specdata


if __name__ == '__main__':
    #test_with_agg_simulation(unitless=True)
    #test_with_map_pair(unitless=True)
    #test_with_simulation(unitless=True)
