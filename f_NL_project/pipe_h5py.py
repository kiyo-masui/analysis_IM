import hdf5_to_algebra as ha

# Initialize JD settings
get_JD_algebra = True
JD_redshift = [0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500,
               1.600, 1.700, 1.800, 1.900, 2.000, 2.100, 2.200, 2.300,
               2.400, 2.500, 2.700, 3.000, 3.500, 4.000, 4.500, 5.000,
               6.000, 7.000, 8.000, 9.000, 10.000]
index = [0,1,2,3,4,5,6,7]
# Actually it is 512 but due to issues write 513
JD_bins = 513
JD_path = "/mnt/raid-project/gmrt/mufma/JD_h5py_catalogs/"
JD_file = "%.3fhalo_catalog000%d.hdf5"
JD_group_1 = 'T_b21'
JD_subgroup_1 = 'T_b21'
JD_group_2 = 'Galaxies_number'
JD_subgroup_2 = 'Galaxies_number'
JD_directory = "/mnt/raid-project/gmrt/mufma/JD_algebra_deltas_%d/"
JD_savefile_1 = 'T_b21_%dind%.3fred.npy'
JD_savefile_2 = 'Galaxies_%dind%.3f.npy'

# Initialize JOE settings
get_algebra = False
redshift = [0.042,0.086,0.130,0.220,0.267,0.316,0.416,0.468,0.523,0.636,
            0.696,0.758,0.889,0.958,1.030,1.185,1.267,1.353,1.538,1.638,
            1.742,1.969,2.092,2.222,2.505,2.660,2.825]
path = "/mnt/raid-project/gmrt/mufma/h5py_catalogs/"
folder_file = "simulation%d/%.3fhalo_catalog.hdf5"
num_bins = 161
catalogs = [range(1,66),range(67,77),range(78,101)]
group = 'Halo_Masses'
subgroup = 'Halo_Masses'
directory = "/mnt/raid-project/gmrt/mufma/algebra_deltas_%d/"
savefile  = 'dm_%dsim%.3fred.npy'

# Convert h5py catalogs to the algebra_npy files
if get_algebra == True:
    for list in catalogs:
        for num in list:
            for red in redshift:
                ha.converter(path + folder_file%(num,red), num_bins, group,
                             subgroup, directory, savefile%(num,red))
            print "finished processing catalog %d"%nums

if get_JD_algebra == True:
    for red in JD_redshift:
        print "Started to work with redshift", red, "\n"
        for ind in index:
            print "Started to work with index", ind
            ha.converter(JD_path + JD_file%(red,ind), JD_bins, JD_group_1,
                         JD_subgroup_1, JD_directory, JD_savefile_1%(ind,red))
            ha.converter(JD_path + JD_file%(red,ind), JD_bins, JD_group_2,
                         JD_subgroup_2, JD_directory, JD_savefile_2%(ind,red))
            print "Index ", ind, " done!"
        print "Done!", "\n"
