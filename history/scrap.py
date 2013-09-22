class CorrelateSingle():
    r"""Class to handle correlation of a single map pair"""
    def __init__(self, parameter_file_or_dict=None):
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      prefix=prefix)

        self.freq_list = sp.array(self.params['freq_list'], dtype=int)
        self.lags = self.params['lags']

        #self.output_root = self.datapath_db.fetch(self.params['output_root'],
        #                                          intend_write=True)
        self.output_root = self.params['output_root']
        self.output_filetag = self.params['output_filetag']
        self.inifile = self.output_root + self.output_filetag + ".ini"

        # Write parameter file.
        kiyopy.utils.mkparents(self.output_root)
        parse_ini.write_params(self.params, self.inifile,
                               prefix=prefix)

        # load the maps, N^-1
        map1 = algebra.make_vect(algebra.load(self.params['map1']))
        map2 = algebra.make_vect(algebra.load(self.params['map2']))
        noise_inv1 = algebra.make_vect(algebra.load(self.params['noise_inv1']))
        noise_inv2 = algebra.make_vect(algebra.load(self.params['noise_inv2']))
        self.pair = MapPair(map1, map2, noise_inv1, noise_inv2, self.freq_list)

    def execute(self):
        r"""calculate the correlation function of the given map pair"""
        print "stuff here"
        # SAVE TO SHELVE NOTPICKLE



