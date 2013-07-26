This folder contains the scripts performing the tidal field reconstrction, as described in arXiv:1202.5804


General information
*******************
Please refer to analysis_IM/README.txt for general setup information


Run the estimator
*****************
To run the estimator, just run the init_file script. The command line should be :
    python kappa/init_file.py


Input and output files
**********************
The input file and output directory must be set directly in init_file.py
- To modify the input file, just modify the variable input_file, which contains the location of the input file
- In order to store information on physical size and axes of the box,
  the estimator requires a specific input type, slightly different from ndarray.
  The type used is 'vect' type, described in analysis_IM/algebra/core.py
- If the input map is just an ndarray, you have to set its information to match the required type.
  To do this, use the function set_info. For example, for an input map in a 500x500x500 Mpc/h box
  add the line :
    (...)
    delta = al.load(input_file)
 -> delta = set_info(delta,(500,500,500))
    delta = preprocess(delta)
    (...)

- The results of the estimator are saved in the repetory output_dir
- You can choose to same the results in binary or numpy format, via the variable out_type
- For example, with :
    output_dir = 'dir/'
    prefix = '2df_'
    out_type = 'binary'
  the content of the repertory dir/ after the run should look like this :
    2df_kappa.dat
    2df_clean.dat
    2df_weight.dat
    (...)


Parameters of the run
*********************
init_file.py contains every tunable parameters for the reconstruction.
Though most of them are not likely to be changed, here is a description of the different parameters,
all of them stored in the variable 'params'

variables :
- sigma_recon : sigma of the Gaussian window used to smooth the density field, before applying the logarithm
- sigma_smooth : sigma of the Gaussian window used to compare large scales in the input map 
                 and in the reconstructed kappa. Only used for visual comparison

list :
- save : list of maps to be saved at the end of the run.
         Every element of this list should be a valid attribute of the KappaEstimator class (see estimator.py)

flags :
- clean : True if you want to apply the weights in k space to the noisy reconstructed kappa.
          The clean map should only be used in 1D power spectra computations of for visual comparisons
- smooth : True if you want to apply the Gaussian smoothing with sigma = sigma_smooth
           The smooth maps shouldonly be used for visual comparisons
- interm_steps : True for the intermediary steps you want to keep. If these steps are not in the save list,
                 they will not be saved. Can be useful if some debug is needed
- fast : two versions of the quadratics estimators computation are implemented, however the fast version
         is likely to produce spurious results in certain cases. Since the faster version is only about
         1.5 to 2x faster than the other one,
         it is recommended to let this flag to False
- auto_rescale : because of the logarithmic smoothing, it is essential that the minimum value of the map
                 after smoothing with the Gaussian window is greater than -1.
                 If not, if this flag is turned to True, a normalization constant is automitically applied
                 so that this condition is respected. The common case is an approximation error during the
                 smoothing, so that the min value becomes -1 - epsilon. The rescaling avoids this problem.
                 Otherwise, if the flag is False, the algorithm raises an error and stops
- verbose : True if you want to be warned when a step of the reconstruction is complete.


I hope this helps.
