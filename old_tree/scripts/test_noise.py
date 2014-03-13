from core import algebra
import numpy as np
from plotting import plot_cube as pc

if __name__ == '__main__':
    root = "/mnt/raid-project/gmrt/tcv/maps/"
    fullnoiseinv_file = root + "sec_A_15hr_41-90_noise_inv_I.npy"
    fullweight_file = "full_noise_inv_weight.npy"
    diagnoise_file = root + "sec_A_15hr_41-90_noise_diag_I.npy"
    diagweight_file = "diag_noise_inv_weight.npy"

    # load the diag(N^-1)
    #noise_inv = algebra.make_mat(algebra.open_memmap(fullnoiseinv_file, mode='r'))
    #noise_inv_diag = noise_inv.mat_diag()
    #algebra.save(fullweight_file, noise_inv_diag)

    # load the diag(N) and take the inverse, write to disk
    #diag_noise = algebra.make_vect(algebra.load(diagnoise_file))
    #diag_noise_inv = 1./diag_noise
    #diag_noise_inv[diag_noise < 1.e-20] = 0.
    #algebra.save(diagweight_file, 1./diag_noise)

    # make movies of both
    outputdir="/cita/d/www/home/eswitzer/movies/"
    pc.make_cube_movie(fullweight_file, "Full N^-1", pc.cube_frame_dir,
                        sigmarange=[0,1e7], outputdir=outputdir, multiplier=1000.,
                        transverse=False)

    pc.make_cube_movie(diagweight_file, "(diag N)^-1", pc.cube_frame_dir,
                        sigmarange=[0,1e7], outputdir=outputdir, multiplier=1000.,
                        transverse=False)

