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
    data_list = []
    for i in range(len(reading_param)):
        data_list.append([])
    with open(file_name, "rb") as halofile:
        header = halofile.read(record_header.size)
        binary_info = halofile.read(halo_struct.size)
        while binary_info:
            halo_info = halo_struct.unpack(binary_info)
            for index in range(len(reading_param)):
                data_list[index].append(constants[const_param[index]] *
                                   halo_info[reading_param[index]])
            binary_info = halofile.read(halo_struct.size)
        halofile.close()
    return data_list


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
