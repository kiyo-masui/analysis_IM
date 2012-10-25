import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
from utils import data_paths as dp
import os


def process_mode_files(clean_path, outpath, n_modes=20, average=False):
    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    datapath_db = dp.DataPath()
    output_root = datapath_db.fetch(clean_path)
    filename = output_root + "/" + "SVD.hd5"
    print filename
    svd_data = h5py.File(filename, "r")

    compilation = {}
    for fieldname in svd_data:
        acc = []
        if average:
            for pair in svd_data[fieldname]:
                acc.append(svd_data[fieldname][pair].value)
        else:
            pair0 = svd_data[fieldname].keys()[0]
            acc.append(svd_data[fieldname][pair0].value)

        compilation[fieldname] = np.mean(np.dstack(acc), axis=-1)

        print fieldname, compilation[fieldname].shape

        # handle the eigenvalues
        if "val" in fieldname:
            filename = outpath + "/" + fieldname + ".dat"
            print "writing ", filename
            fileobj = open(filename, "w")
            #data = compilation[fieldname][0]
            data = np.abs(np.sort(-compilation[fieldname][0]))
            #data /= data[0]
            for index in range(data.shape[0]):
                fileobj.write("%10.15g\n" % data[index])

            fileobj.close()

        # handle the modes
        if "modes" in fieldname:
            filename = outpath + "/" + fieldname + ".dat"
            print "writing ", filename
            fileobj = open(filename, "w")
            data = compilation[fieldname][0:n_modes]

            for index in range(data.shape[1]):
                entry = tuple(data[:, index].tolist())
                fileobj.write(("%10.5g " * n_modes + "\n") % entry)

            fileobj.close()

    svd_data.close()


process_mode_files("GBT_15hr_map_oldcal_cleaned_path_Eric",
                   "./pwrspec_plots/GBT_15hr_map_oldcal/")

process_mode_files("GBT_15hr_map_oldcalpol_cleaned_path_Eric",
                   "./pwrspec_plots/GBT_15hr_map_oldcalpol/")
