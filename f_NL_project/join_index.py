import numpy as np
from core import algebra

def join_index(list_filename, savename):
    """
    Join all files from list_filename and save it to the savename
    """
    read_file = open(list_filename[0], "r")
    save_map = algebra.load(read_file)
    read_file.close()
    for num in range(1,len(list_filename)):
        read_file = open(list_filename[num], "r")
        add_map = algebra.load(read_file)
        save_map = save_map + add_map
        read_file.close()
    file = open(savename, "w")
    algebra.save(file, save_map)
    file.close()

