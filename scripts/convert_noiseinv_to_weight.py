import sys
from utils import data_paths as dp
from core import algebra

def convert_noiseinv_to_weight(mapkey):
    datapath_db = dp.DataPath()
    filedb = datapath_db.fetch(mapkey)[1]
    map_cases = datapath_db.fileset_cases(mapkey, "section;maptype")

    for section in map_cases['section']:
        noiseinv_file = filedb[section + ";noise_inv"]
        noiseweight_file = filedb[section + ";noise_weight"]
        print noiseinv_file, noiseweight_file

        noise_inv = algebra.make_mat(algebra.open_memmap(noiseinv_file, mode='r'))
        noise_inv_diag = noise_inv.mat_diag()
        algebra.save(noiseweight_file, noise_inv_diag)

if __name__ == '__main__':
    convert_noiseinv_to_weight(sys.argv[1])
