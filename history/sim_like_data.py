makesimlikedata_init = {
        "shelvefile_in": "file",
        "shelvefile_out": "file",
        "treatments":["one", "two"],
        "outputdir": "dir",
    }

makesimlikedata_prefix = 'sld_'

class MakeSimLikeData(object):
    """convert the aggregated sims into a structure that looks just like the
    data. This is used to bin the 2d theory onto 1d identically to the real
    data.
    """

    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        np.seterr(invalid='raise')

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          makesimlikedata_init,
                                          prefix=makesimlikedata_prefix)

        print self.params["shelvefile_in"], "->", self.params["shelvefile_out"]
        self.stats_in = shelve.open(self.params["shelvefile_in"], "r")
        self.stats_out = shelve.open(self.params["shelvefile_out"], "w")
        self.treatments_in = self.stats_in["results"].keys()
        self.treatments_out = self.params["treatments"]

        print "MakeSimLikeData: input treatments: ", self.treatments_in
        print "MakeSimLikeData: output treatments: ", self.treatments_out

        if self.treatments_in[0] != "0modes":
            print "Transfer functions must be wrt only 0modes"
            return

    def execute(self, processes):
        """produce a list of files to combine and run"""
        transfer_compilation = {}
        for treatment in self.treatments_out:
            print treatment
            print self.stats_out["results"][treatment].keys()
            #stats_1d_in = self.stats_in["stats"]["0modes"]["pk_1d_stat"]["mean"]
            #stats_1d_out = self.stats_out["stats"][treatment]["pk_1d_stat"]["mean"]
            stats_2d_in = self.stats_in["stats"]["0modes"]["pk_2d_stat"]["mean"]
            stats_2d_out = self.stats_out["stats"][treatment]["pk_2d_stat"]["mean"]

            #counts_1d_in = self.stats_in["stats"]["0modes"]["pk_1d_counts"]["mean"]
            #counts_1d_out = self.stats_out["stats"][treatment]["pk_1d_counts"]["mean"]
            counts_2d_in = self.stats_in["stats"]["0modes"]["pk_2d_counts"]["mean"]
            counts_2d_out = self.stats_out["stats"][treatment]["pk_2d_counts"]["mean"]

            counts_prod = counts_2d_in * counts_2d_out
            transfer_2d = stats_2d_out / stats_2d_in
            transfer_2d[counts_prod == 0] = 0.

            #k_1d = self.stats_in["k_1d"]
            #k_1d_from_2d = self.stats_in["k_1d_from_2d"]
            kx_2d = self.stats_in["kx_2d"]
            ky_2d = self.stats_in["ky_2d"]

            transfer_2d_plot = copy.deepcopy(transfer_2d)
            transfer_2d_plot[transfer_2d_plot < 0.] = 0.
            transfer_2d_plot[transfer_2d_plot > 1.] = 1.

            outplot_file = "%s/transfer_2d_%s.png" % (self.params['outputdir'], treatment)
            logkx = np.log10(kx_2d['center'])
            logky = np.log10(ky_2d['center'])
            plot_slice.simpleplot_2D(outplot_file, transfer_2d_plot,
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "2D beam transfer", "T")

            transfer_compilation[treatment] = transfer_2d

        #transferfile = shelve.open(self.params["transferfile"], protocol=0)
        #transferfile["transfer_2d"] = transfer_compilation

        transferfile = h5py.File(self.params["transferfile"], "w")
        for treatment in self.treatments_out:
            transferfile[treatment] = transfer_compilation[treatment]

        transferfile.close()

        self.stats_in.close()
        self.stats_out.close()


