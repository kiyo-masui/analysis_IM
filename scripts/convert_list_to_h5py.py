import h5py
import numpy as np

def convert_list_to_h5py(file_adress,directory_list,datasets_list,halo_catalog_adress):
    # part regarding HDF5 file and its structure
    file_name = file_adress
    create_file = h5py.File(file_name)
    directories = directory_list
    data_sets = datasets_list
    # part regarding data we are going to safe in HDF5 file
    halo_catalog_directory = halo_catalog_adress
    file_load = open(halo_catalog_directory)
    data_file = file_load.read()
    # part regarding the fact that type(data_file) is a "list"
    # task: create array out of it
    # -1 in summation since the last entry of data is empty ''
    data_file = data_file.split('\n')
    x_pos = np.zeros(len(data_file)-1)
    y_pos = np.zeros(len(data_file)-1)
    z_pos = np.zeros(len(data_file)-1)
    x_vel = np.zeros(len(data_file)-1)
    y_vel = np.zeros(len(data_file)-1)
    z_vel = np.zeros(len(data_file)-1)
    halo_masses = np.zeros(len(data_file)-1)
    for index in range(len(data_file)-1):
        data_file[index] = data_file[index].split(' ')
    for index in range(len(data)-1):
        x_pos[index] = data_file[index][0]
        y_pos[index] = data_file[index][1]
        z_pos[index] = data_file[index][2]
        x_vel[index] = data_file[index][3]
        y_vel[index] = data_file[index][4]
        z_vel[index] = data_file[index][5]
        halo_masses[index] = data_file[index][6]
    # creating structure of HDF5 file
    # and writing position and velocities arrays in data_sets
    for name in directories:
        create_file.create_group(name)
    for name in data_sets:
        create_file['Positions'].create_dataset(name, data = eval("%s_pos" % name))
        create_file['Velocities'].create_dataset(name, data = eval("%s_vel" % name))
    create_file['Halo_Masses'].create_dataset("Halo_Masses", data = halo_masses)

convert_list_to_h5py('Halo_Catalog.hdf5',['Positions','Velocities','Halo_Masses'],['x','y','z'],"/cita/h/home-2/mufma/code/analysis_IM/scripts/halo_data.dat")
