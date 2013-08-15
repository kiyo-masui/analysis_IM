from core import fitsGBT
import os
from numpy import *
import sys

# Need a .npy file and a .npy.meta file for the Map (combined sessions)
# There must also be a file fir_rebin_1hr_corrected_map_01-06_clean_map_I_775.npy.meta in this folder!
map_npy_file = '/cita/h/home-1/anat/analysis_IM/input/anat/npy_maps/fir_rebin_1hr_corrected_map_01-06_clean_map_I_775.npy'

# This is where the individual scan files are located.
in_path   = '/mnt/raid-project/gmrt/anat/maps/rebinned_corrected/GBT11B_055/'

# Where to put the corrected scan files.
out_path  = '/mnt/raid-project/gmrt/anat/maps/rebinned_corrected_twice_testing/GBT11B_055/'

def read_map_data(map_npy_file):

   """  Reads the median map.
        Need map file as .npy and also a .npy.meta file.
   """
   in_file = open(map_npy_file+'.meta','r')
   info = in_file.readline()
   in_file.close()
   info1 = safe_eval(info)
   data_map = load(map_npy_file)
   dec_map = info1['dec_delta']*(arange(data_map.shape[2]) - data_map.shape[2]/2) + info1['dec_centre']
   freq_map = info1['freq_delta']*(arange(data_map.shape[0]) - data_map.shape[0]/2) + info1['freq_centre']
   ra_map = info1['ra_delta']*(arange(data_map.shape[1]) - data_map.shape[1]/2) + info1['ra_centre']

   return freq_map,dec_map,ra_map,data_map


def get_val(ra_val,dec_val,ra_map,dec_map,data,fnum,scan_point):

   """    Returns Temp from the map closest to (RA,DEC) in the scan.
   """
   if (ra_val > ra_map[0]) or (ra_val < ra_map[len(ra_map)-1]):
      return scan_point   # Outside the map dimensions, so return the scan value rather than the map value.
   #
   #   Alert the user here that a few scan points are too close to the map boundary 
   #   to compute a derivative. 
   #
   if (dec_val < dec_map[0]) or (dec_val > dec_map[len(dec_map)-1]):
      return scan_point   # Outside the map dimensions.

   for i in range(0,len(ra_map)):
      if ra_map[i] < ra_val:
         break
   for j in range(0,len(dec_map)):
      if dec_map[j] > dec_val:
         break
  
   # Interpolate in RA
   avg_minus = data[fnum][i-1][j-1] + ( ((data[fnum][i][j-1]-data[fnum][i-1][j-1])/(ra_map[i]-ra_map[i-1]))*(ra_val-ra_map[i-1]) )
   avg_plus  = data[fnum][i-1][j] + ( ((data[fnum][i][j]-data[fnum][i-1][j])/(ra_map[i]-ra_map[i-1]))*(ra_val-ra_map[i-1]) )
   # Interpolate in DEC
   avg = avg_minus + ( ((avg_plus-avg_minus)/(dec_map[j]-dec_map[j-1])) * (dec_val-dec_map[j-1]) )


   return avg   

  
def get_pointing_error(fnum,scan_ra,scan_dec,map_ra,map_dec,map_data,scan_data):

   """ Computes the pointing errors, and returns the corrected scan.
       dRA and dDEC are arbitrary, but shouldn't be too small or too large.
       fnum can be 0,1,2,3 for 4 frequency bands 875,825,775,725 MHz.
   """
   dRA = 0.1
   dDEC = 0.1

   aver = mean(scan_data)
   scan_data = [val - aver for val in scan_data]
   temp_plus = [get_val(scan_ra[i]+dRA,scan_dec[i],map_ra,map_dec,map_data,fnum,scan_data[i]) for i in range(0,len(scan_ra))]
   temp_plus = [temp_plus[i] - mean(temp_plus) for i in range(0,len(temp_plus)) if temp_plus[i] != -1.] 
   temp = [get_val(scan_ra[i],scan_dec[i],map_ra,map_dec,map_data,fnum,scan_data[i]) for i in range(0,len(scan_ra))]
   temp = [temp[i] - mean(temp) for i in range(0,len(temp))]
   prime_RA = [(temp_plus[i]-temp[i])/dRA for i in range(0,len(scan_ra))]

   temp_plus_temp = [get_val(scan_ra[i],scan_dec[i]+dDEC,map_ra,map_dec,map_data,fnum,scan_data[i]) for i in range(0,len(scan_ra))]
   temp_plus = [val for val in temp_plus_temp if val != -1.]
   temp_plus = [temp_plus[i] - mean(temp_plus) for i in range(0,len(temp_plus))]
   prime_DEC = [(temp_plus[i]-temp[i])/dDEC for i in range(0,len(scan_ra))]

   deltaT = [scan_data[i]-temp[i] for i in range(0,len(scan_ra))]
   
   A1 = mean([deltaT[i]*prime_RA[i] for i in range(0,len(scan_ra))])
   A2 = mean([prime_RA[i]*prime_RA[i] for i in range(0,len(scan_ra))])
   A3 = mean([prime_RA[i]*prime_DEC[i] for i in range(0,len(scan_ra))]) 
   B1 = mean([deltaT[i]*prime_DEC[i] for i in range(0,len(scan_ra))])
   B2 = A3
   B3 = mean([prime_DEC[i]*prime_DEC[i] for i in range(0,len(scan_ra))])

   DeltaRA = (A1*B3 - B1*A3) / (A2*B3 - B2*A3)
   DeltaDEC = (A1*B2 - B1*A2) / (A3*B2 - B3*A2)
   new_scan = [aver + scan_data[i] - (prime_RA[i]*DeltaRA + prime_DEC[i]*DeltaDEC) for i in range(0,len(scan_ra))]
   return DeltaRA, DeltaDEC, new_scan
 

def main():

#   sessions = ['1', '2', '3', '4', '5', '6']
   sessions = ['1']
   for session in sessions:
      in_list_of_files  = [in_path  + name for name in os.listdir(in_path) if (name[0] == '0' and name[1] == session)]
      out_list_of_files = [out_path + name for name in os.listdir(in_path) if (name[0] == '0' and name[1] == session)]
 
      # Get info from the Map (combined sessions)
      freq_map,dec_map,ra_map,data_map = read_map_data(map_npy_file)
      for file_loop in range(0,len(in_list_of_files)):

         Reader = fitsGBT.Reader(in_list_of_files[file_loop])
         Writer = fitsGBT.Writer()

         for loop in range(0,10):
 
            Data = Reader.read(loop,0)   # Single scan.
            Data.calc_pointing()

            # Take care that the data doesn't contain NaN or -- etc, else you will get a warning.
            for fnum in range(0,4):
 
               scan_val = [Data.data[i][0][0][fnum] for i in range(0,len(Data.ra))]
               deltaRA,deltaDEC,new_scan = get_pointing_error(fnum,Data.ra,Data.dec,ra_map,dec_map,data_map,scan_val)
               # print deltaRA,deltaDE

               for i in range(0,len(Data.ra)):
                  Data.data[i][0][0][fnum] = new_scan[i]
            Writer.add_data(Data)
         Writer.write(out_list_of_files[file_loop])


if __name__ == "__main__":
   main()
