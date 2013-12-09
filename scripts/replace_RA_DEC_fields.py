import glob
import os

from core import fitsGBT
import kiyopy.utils

def replace(in_dir, out_dir):

    kiyopy.utils.mkdir_p(out_dir)
    files = glob.glob(in_dir + '/*.fits')
    data_dir = os.getenv('GBT_DATA')
    for file in files:
        out_fname = out_dir + file.split('/')[-1]
        data_fname = data_dir + '/'.join(file.split('/')[-2:])
        if not os.path.isfile(data_fname):
            msg = "Raw data file: %s missing" % data_fname
            raise ValueError(msg)
        raw_Reader = fitsGBT.Reader(data_fname)
        raw_data_list = raw_Reader.read()
        Reader = fitsGBT.Reader(file)
        data_list = Reader.read()
        if len(raw_data_list) != len(data_list):
            raise ValueError("Wrong number of scans in raw data.")
        for ii in range(len(data_list)):
            raw_data = raw_data_list[ii]
            ra = raw_data.field['RA']
            ra_format = raw_data.field_formats['RA']
            ra_axes = raw_data.field_axes['RA']
            data_list[ii].set_field('RA', ra, axis_names=ra_axes,
                    format=ra_format)
            
            dec = raw_data.field['DEC']
            dec_format = raw_data.field_formats['DEC']
            dec_axes = raw_data.field_axes['DEC']
            data_list[ii].set_field('DEC', dec, axis_names=dec_axes,
                    format=dec_format)
        Writer = fitsGBT.Writer(data_list)
        Writer.write(out_fname)
        



if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print ("Usage: python scripts/replace_RA_DEC_fields.py input_directory"
               " output_directory")
    else:
        in_dir = sys.argv[1]
        out_dir = sys.argv[2]
        replace(in_dir, out_dir)
