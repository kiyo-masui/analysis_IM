from core import algebra
import numpy as np

rootdir = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal/"
filelist_in = ["secA_1hr_41-18_noise_weight_I_800_orig.npy",
               "secB_1hr_41-18_noise_weight_I_800_orig.npy",
               "secC_1hr_41-18_noise_weight_I_800_orig.npy",
               "secD_1hr_41-18_noise_weight_I_800_orig.npy"]

filelist_out = ["secA_1hr_41-18_noise_weight_I_800.npy",
                "secB_1hr_41-18_noise_weight_I_800.npy",
                "secC_1hr_41-18_noise_weight_I_800.npy",
                "secD_1hr_41-18_noise_weight_I_800.npy"]

for (file_in, file_out) in zip(filelist_in, filelist_out):
    filename_in = rootdir + file_in
    filename_out = rootdir + file_out
    print filename_in, filename_out

    weight_in = algebra.make_vect(algebra.load(filename_in))
    dec_axis = weight_in.get_axis("dec")
    goodregion = np.where(dec_axis < -0.2)
    weight_in[:,:,goodregion] = 0.
    algebra.save(filename_out, weight_in)
