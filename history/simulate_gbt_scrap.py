# TODO: phase this out
def generate_delta_sim(input_file, output_file):
    r"""make the map with the temperature divided out (delta)"""
    print "reading %s -> %s (dividing by T_b(z))" % (input_file, output_file)

    simmap = algebra.make_vect(algebra.load(input_file))
    freq_axis = simmap.get_axis('freq') / 1.e6
    z_axis = units.nu21 / freq_axis - 1.0

    simobj = corr21cm.Corr21cm()
    T_b = simobj.T_b(z_axis)*1e-3

    simmap /= T_b[:, np.newaxis, np.newaxis]

    print "saving to" + output_file
    algebra.save(output_file, simmap)


# TODO: phase this out
def generate_proc_sim(input_file, weightfile, output_file,
                      meansub=False, degrade=False):
    r"""make the maps with various combinations of beam conv/meansub"""
    print "%s -> %s (beam, etc.)" % (input_file, output_file)
    simmap = algebra.make_vect(algebra.load(input_file))

    if degrade:
        print "performing common resolution convolution"
        beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
        freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
        freq_data *= 1.0e6
        beam_diff = sp.sqrt(max(1.1 * beam_data) ** 2 - (beam_data) ** 2)
        common_resolution = beam.GaussianBeam(beam_diff, freq_data)
        # Convolve to a common resolution.
        simmap = common_resolution.apply(simmap)

    if meansub:
        print "performing mean subtraction"
        noise_inv = algebra.make_vect(algebra.load(weightfile))
        means = sp.sum(sp.sum(noise_inv * simmap, -1), -1)
        means /= sp.sum(sp.sum(noise_inv, -1), -1)
        means.shape += (1, 1)
        simmap -= means
        # the weights will be zero in some places
        simmap[noise_inv < 1.e-20] = 0.

    # extra sanity?
    simmap[np.isinf(simmap)] = 0.
    simmap[np.isnan(simmap)] = 0.

    print "saving to" + output_file
    algebra.save(output_file, simmap)


