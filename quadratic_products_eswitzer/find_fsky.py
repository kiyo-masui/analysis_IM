import numpy as np
from utils import data_paths as dp
from core import algebra
from map import physical_gridding as pg
from utils import batch_handler as bh


def find_avg_fsky(map_key, tack_on=None, refinement=2, pad=5, order=1):
    """Take all the pairs that enter the autopower, open their weight files and
    find the fsky for each data treatment
    In a pipeline
        fsky = find_avg_fsky(self.params["map_key"],
                            tack_on=self.params["tack_on"],
                            refinement=self.params["refinement"],
                            pad=self.params["pad"],
                            order=self.params["order"])
    """
    fsky = {}
    datapath_db = dp.DataPath()

    map_cases = datapath_db.fileset_cases(map_key,
                                          "pair;type;treatment")

    # This code is essentially verbatim for the permutation in the real
    # autopower
    unique_pairs = dp.GBTauto_cross_pairs(map_cases['pair'],
                                          map_cases['pair'],
                                          cross_sym="_with_")

    treatment_list = map_cases['treatment']

    for treatment in treatment_list:
        for item in unique_pairs:
            dbkeydict = {}
            mapset0 = (map_key, item[0], treatment)
            mapset1 = (map_key, item[1], treatment)
            dbkeydict['noiseinv1_key'] = "%s:%s;noise_inv;%s" % mapset0
            dbkeydict['noiseinv2_key'] = "%s:%s;noise_inv;%s" % mapset1
            files = dp.convert_dbkeydict_to_filedict(dbkeydict,
                                                datapath_db=datapath_db,
                                                tack_on=tack_on)

            print files['noiseinv1_key'], files['noiseinv2_key']
            weight1 = algebra.make_vect(algebra.load(files['noiseinv1_key']))
            weight2 = algebra.make_vect(algebra.load(files['noiseinv2_key']))

            physweight1 = bh.repackage_kiyo(pg.physical_grid(
                                             weight1,
                                             refinement=refinement,
                                             pad=pad, order=order))

            physweight2 = bh.repackage_kiyo(pg.physical_grid(
                                             weight2,
                                             refinement=refinement,
                                             pad=pad, order=order))

            #fsky = np.sum(physweight1 * physweight2)**2
            #fsky /= np.sum(physweight1**2 * physweight2**2)
            #fsky /= float(physweight1.size)
            fsky = np.sum(weight1 * weight2)**2
            fsky /= np.sum(weight1**2 * weight2**2)
            fsky /= float(weight1.size)
            print "volume factor in noise weight: ", fsky

    return fsky



