import h5py
import convert_binary_to_h5py as cb

# Initialize settings
JD = True
JOE = False

index = [0,1,2,3,4,5,6,7]
simulation_list = [range(1,66),range(67,77),range(78,101)]

JD_redshift = [0.700, 0.600, 0.500, 0.400, 0.300]
"""
JD_redshift = [0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700,
               1.800, 1.900, 2.000, 2.100, 2.200, 2.300, 2.400, 2.500, 2.600,
               2.700, 3.000, 3.500, 4.000, 4.500, 5.000, 6.000, 7.000, 8.000,
               9.000, 10.000]
"""
JOE_redshift = [0.042, 0.086, 0.130, 0.220, 0.267, 0.316, 0.416, 0.468, 0.523,
                0.636, 0.696, 0.758, 0.889, 0.958, 1.030, 1.185, 1.267, 1.353,
                1.538, 1.638, 1.742, 1.969, 2.092, 2.222, 2.505, 2.660, 2.825]

JD_path = "/mnt/scratch-3week/emberson/21cm/"
"""
JD_path = "/mnt/raid-project/gmrt/mufma/JD_halo_catalogs/"
"""
JD_h5py_path = '/mnt/raid-project/gmrt/mufma/JD_h5py_catalogs/'
JOE_path = "/mnt/raid-project/gmrt/mufma/halo_catalogs/"
JOE_h5py_path = '/mnt/raid-project/gmrt/mufma/h5py_catalogs/'

JD_constants = [1., 1.]
JOE_constants = [505. * (1./2048.), 7.9955, 1.3313e10]

JD_header = "i"
JOE_header = "f"

JD_body = "<" + "f"*4
JOE_body = "<" + "f"*28

JD_read_columns = [0, 1, 2, 3]
JOE_read_columns = [0, 1, 2, 6, 7, 8, 16]

JD_const_choice = [0, 0, 0, 1]
JOE_const_choice = [0, 0, 0, 1, 1, 1, 2]

JD_groups = ['Positions', 'Halo_Masses']
JOE_groups = ['Positions', 'Velocities', 'Halo_Masses']

JD_subgroups = [['x', 'y', 'z'],['Halo_Masses']]
JOE_subgroups = [['x', 'y', 'z'], ['x', 'y','z'],['Halo_Masses']]

# Converting JD binary files to the h5py catalogs
if JD == True:
    print "\n", "working on JD simulations", "\n"
    for red in JD_redshift:
        print "working on redshift:", red
        for ind in index:
            print "started index:", ind
            info = cb.binary_reader(JD_path + "%.3fhalos_GPC/"%red +
                                    "%.3fhalo000%d.bin"%(red,ind),
                                    JD_header, JD_body, JD_constants,
                                    JD_read_columns, JD_const_choice)
            """
            info = cb.binary_reader(JD_path + "%.3fhalos_GPC/"%red +
                                    "%.3fhalo000%d.bin"%(red,ind),
                                    JD_header, JD_body, JD_constants,
                                    JD_read_columns, JD_const_choice)
            """
            cb.convert_array_to_h5py(JD_h5py_path +
                                     "/%.3fhalo_catalog000%d.hdf5"%(red,ind),
                                     JD_groups, JD_subgroups, info)

# Converting JOE binary files to the h5py catalogs
if JOE == True:
    print "\n", "working on JOE simulations", "\n"
    for list in simulation_list:
        for num in list:
            print "working on simulation:", num
            for red in JOE_redshift:
                print "started redshift:", red
                info = cb.binary_reader(JOE_path + "simulation%d"%num +
                                        "/%.3fhalo0.dat"%red,
                                        JOE_header, JOE_body, JOE_constants,
                                        JOE_read_columns, JOE_const_choice)
                cb.convert_array_to_h5py(JOE_h5py_path + 'simulation%d'%num +
                                         "/%.3fhalo_catalog.hdf5"%red,
                                         JOE_groups, JOE_subgroups, info)

