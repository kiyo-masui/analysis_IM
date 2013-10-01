import glob
import os

from core import fitsGBT

def replace(in_root, out_root):
    files = glob.glob(file_root + '/GBT*/*.fits')
    print files


if __name__ == '__main__':
    import sys
    in_root = sys.argv[0]
    out_root = sys.argv[1]
    replace(in_root, out_root)
