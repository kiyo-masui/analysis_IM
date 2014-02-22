def dot_mode_functions(path_key1, path_key2):
    (amplitudes1, mode_functions_l1, mode_functions_r1) = \
                                    process_mode_files(path_key1)
    (amplitudes2, mode_functions_l2, mode_functions_r2) = \
                                    process_mode_files(path_key2)

    mode_functions_l1_avg = np.mean(mode_functions_l1, axis=2)
    mode_functions_r1_avg = np.mean(mode_functions_r1, axis=2)
    mode_functions_l2_avg = np.mean(mode_functions_l2, axis=2)
    mode_functions_r2_avg = np.mean(mode_functions_r2, axis=2)
    (n_modes, n_freq) = mode_functions_l1_avg.shape

    filename = "plot_data/%s_dot_%s_modes.txt" % (path_key1, path_key2)
    fileobj = open(filename, "w")
    for index in range(0, n_modes):
        numer = np.dot(mode_functions_l1_avg[index, :],
                       mode_functions_l2_avg[index, :])

        denom = np.dot(mode_functions_l1_avg[index, :],
                       mode_functions_l1_avg[index, :])
        fileobj.write(("%10.5g" * 3 + "\n") % (numer, denom, numer/denom))
    fileobj.close()

    filename1 = "plot_data/%s_mode_functions.txt" % path_key1
    fileobj1 = open(filename1, "w")
    filename2 = "plot_data/%s_mode_functions.txt" % path_key2
    fileobj2 = open(filename2, "w")
    for index in range(0, n_freq):
        fileobj1.write(("%10.5g " * 20 + "\n") % tuple(mode_functions_l1_avg[0:20, index].tolist()))
        fileobj2.write(("%10.5g " * 20 + "\n") % tuple(mode_functions_l2_avg[0:20, index].tolist()))

    fileobj1.close()
    fileobj2.close()


def plot_six_svd_amplitudes(vals, fig_filename):
    n_vals = len(vals)

    plt.xlim((0.8,300))
    plt.loglog(abs(np.sort(-vals[0] / float(n_vals))), 'b.',
               abs(np.sort(-vals[1] / float(n_vals))), 'g.',
               abs(np.sort(-vals[2] / float(n_vals))), 'r.',
               abs(np.sort(-vals[3] / float(n_vals))), 'c.',
               abs(np.sort(-vals[4] / float(n_vals))), 'm.',
               abs(np.sort(-vals[5] / float(n_vals))), 'y.',
               linestyle='None')
    plt.title("Mode amplitudes")
    plt.xlabel("Mode number")
    plt.ylabel("Amplitude")
    plt.savefig(fig_filename, dpi=200)
    plt.clf()
    #print 'Mean noise: ', np.sum(vals) / n_vals
    #print 'Largest eigenvalues/vals: ',
    #print np.sort(vals / n_vals)[-10:]


def plot_six_svd_modes(vals, fig_filename):
    n_vals = len(vals)

    plt.xlim((0,231))
    plt.plot(vals[0], 'b-',
               vals[1], 'g-',
               vals[2], 'r-',
               vals[3], 'c-',
               vals[4], 'm-',
               vals[5], 'y-',
            )
    plt.title("Modes")
    plt.xlabel("Frequency index")
    plt.ylabel("Mode amplitude")
    plt.savefig(fig_filename, dpi=200)
    plt.clf()
    #print 'Mean noise: ', np.sum(vals) / n_vals
    #print 'Largest eigenvalues/vals: ',
    #print np.sort(vals / n_vals)[-10:]


def plot_svd_series(path_key):
    datapath_db = data_paths.DataPath()
    root = datapath_db.fetch(path_key)
    print root

    pairs = ["A_with_B", "A_with_C", "A_with_D",
             "B_with_C", "B_with_D", "C_with_D"]

    svdblob = [cPickle.load(open(root + "SVD_pair_%s.pkl" % pair, "r")) \
               for pair in pairs]
    (amp_ind, mode0_ind, mode1_ind) = (0,1,2)

    amp_set = [ svdblob[pair_ind][amp_ind] for pair_ind in range(0,6) ]
    plot_six_svd_amplitudes(amp_set, path_key + "_mode_amplitudes.png")

    moderange = range(0,15)
    filelist = [ "%s_svd_mode_%d.png" % (path_key, mode_index) for mode_index in moderange ]
    for (filename, index) in zip(filelist, moderange):
        print filename
        mode_set = [ svdblob[pair_ind][mode0_ind][index] for pair_ind in range(0,6) ]
        plot_six_svd_modes(mode_set, filename)

def plot_svd_series_comparison(path_key1, path_key2):
    datapath_db = data_paths.DataPath()
    root1 = datapath_db.fetch(path_key1)
    root2 = datapath_db.fetch(path_key2)

    pairs = ["A_with_B", "A_with_C", "A_with_D",
             "B_with_C", "B_with_D", "C_with_D"]

    svdblob1 = [cPickle.load(open(root1 + "SVD_pair_%s.pkl" % pair, "r")) \
               for pair in pairs]
    svdblob2 = [cPickle.load(open(root2 + "SVD_pair_%s.pkl" % pair, "r")) \
               for pair in pairs]
    n_modes = len(svdblob1[0][amp_ind])

    (amp_ind, mode0_ind, mode1_ind) = (0,1,2)

    #amp_set1 = [ np.sort(-(svdblob1[pair_ind][amp_ind])) for pair_ind in range(0,6) ]
    #amp_set2 = [ np.sort(-(svdblob2[pair_ind][amp_ind])) for pair_ind in range(0,6) ]

    moderange = range(0,100)
    filelist = [ "plot_data/%s_%s_svd_mode_%d.txt" % (path_key1, path_key2, mode_index) for mode_index in moderange ]
    for (filename, index) in zip(filelist, moderange):
        amplitudes1 = np.zeros((n_modes, 6))
        for pair_ind in range(0, 6):
            amplitudes1[:,pair_ind] = np.abs(np.sort(-svdblob1[pair_ind][mode_ind]))
            amplitudes2[:,pair_ind] = np.abs(np.sort(-svdblob2[pair_ind][mode_ind]))

        amp_out = np.mean(amplitudes, axis=1)
        #mode_set1 = [ svdblob1[pair_ind][mode0_ind][index] for pair_ind in range(0,6) ]
        #mode_set2 = [ svdblob2[pair_ind][mode0_ind][index] for pair_ind in range(0,6) ]
        #n_freq = len(mode_set1[0])
        #if index < 16:
        #    fileobj = open(filename, "w")
        #    for find in range(0, n_freq):
        #        fileobj.write(("%10.15f " * 2 + "\n") % (mode_set1[0][find], mode_set2[0][find]))
        #    fileobj.close()

        #print np.dot(mode_set1[0][index]




