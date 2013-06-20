import struct
import numpy as np
import h5py
# description
"""
binary_reader(f,h,b,c) gets binary file located at f, which has structure:
header specified by h and body specified by b. c is just converting factors to
physical units.
It returns an array
convert_array_to_h5py() takes array from binary_reader and put it in h5py file
which has directories and datasets
"""

def binary_reader(file_directory, header_struct, body_struct, constants,
                  reading_param, const_param):
    # binary file we want to convert
    file_name = file_directory
    # set sequence in which binary file would be read
    record_header = struct.Struct(header_struct)
    halo_struct = struct.Struct(body_struct)
    # scaling factors to physical units
    # velocity - km/s mass - solar_mass
    # distance - Mpc/h

    # reading binary : opened_file.read(move_cursor.size)
    # it moves cursor along
    # the "opened_file" on "move_cursor.size" steps
    # the while loop works until the variable binary_info becomes False
    list = []
    for i in range(len(reading_param)):
        list.append([])
    with open(file_name, "rb") as halofile:
        header = halofile.read(record_header.size)
        binary_info = halofile.read(halo_struct.size)
        while binary_info:
            halo_info = halo_struct.unpack(binary_info)
            for index in range(len(reading_param)):
                list[index].append(constants[const_param[index]] *
                                   halo_info[reading_param[index]])
            binary_info = halofile.read(halo_struct.size)
        halofile.close()
    return list


def convert_array_to_h5py(file_name, directory_list, datasets_list, data):
    # part regarding HDF5 file and its structure
    create_file = h5py.File(file_name,'a')
    directories = directory_list
    data_sets = datasets_list
    data_save = data

    # creating structure of HDF5 file
    # and writing position and velocities arrays in data_sets
    reading_index = 0
    for name_index in range(len(directories)):
        create_file.create_group(directories[name_index])
        for data_index in range(len(data_sets[name_index])):
            create_file[directories[name_index]].create_dataset(
                data_sets[name_index][data_index],
                data=data_save[reading_index])
            reading_index = reading_index + 1
    create_file.close()


if __name__ == "__main__":
#####  JD'S CATALOGS  #########
    """
    redshift  = [0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500,
1.600, 1.700, 1.800, 1.900, 2.000, 2.100, 2.200, 2.300, 2.400, 2.500, 2.600,
2.700, 3.000, 3.500, 4.000, 4.500, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000]
    """
    redshift = [2.700, 3.000, 3.500, 4.000, 4.500, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000]
    index = [0,1,2,3,4,5,6,7]
    sim_dir = "/mnt/raid-project/gmrt/mufma/JD_halo_catalogs/"
    for red in redshift:
        print "we are at redshift:", red
        for ind in index:
            print ind
            info = binary_reader(sim_dir + "%.3fhalos_GPC/"%red +
                                "%.3fhalo000%d.bin"%(red,ind),
                                "i", "<" + "f"*4,
                                [1., 1.],
                                [0, 1, 2, 3], [0, 0, 0, 1])
            convert_array_to_h5py('/mnt/raid-project/gmrt/mufma/JD_h5py_catalogs/'
                                  + "/%.3fhalo_catalog000%d.hdf5"%(red,ind),
                                  ['Positions', 'Halo_Masses'],
                                  [['x', 'y', 'z'],['Halo_Masses']],
                                  info)


#####  JOE'S CATALOGS  #########
"""
    list = [0.042,0.086,0.130,0.220,0.267,0.316,0.416,0.468,0.523,0.636,
            0.696,0.758,0.889,0.958,1.030,1.185,1.267,1.353,1.538,1.638,
            1.742,1.969,2.092,2.222,2.505,2.660,2.825]
    sim_dir = "/mnt/raid-project/gmrt/mufma/halo_catalogs/"
    for num in [1]:
        print "Using catalog number : %d"%num
        for red in list:
            print "Reading binary file at redshift : %.3f"%red
            info = binary_reader(sim_dir + "simulation%d"%num +"/%.3fhalo0.dat"%red,
                                "f", "<" + "f"*28,
                                [505. * (1./2048.), 7.9955, 1.3313e10],
                                [0, 1, 2, 6, 7, 8, 16], [0, 0, 0, 1, 1, 1, 2])
            print "Converting binary to h5py at %.3f"%red
            convert_array_to_h5py('/mnt/raid-project/gmrt/mufma/h5py_catalogs/'
                                  + 'simulation%d'%num +"/%.3fhalo_catalog.hdf5"%red,
                                  ['Positions', 'Velocities', 'Halo_Masses'],
                                  [['x', 'y', 'z'], ['x', 'y','z'],['Halo_Masses']],
                                  info)
"""
