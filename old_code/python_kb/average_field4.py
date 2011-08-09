from numpy import *

fname1 = "/home/scratch/kbandura/GBT08B_037/np_data/field4_d"
fnameout = "/home/scratch/kbandura/GBT08B_037/np_data/field4_"
#dropped the incomplete field map on d3
#fname2 = ["1_","1_","2_","2_","2_","3_","3_","3_","3_","3_","5_","5_","5_","5_"]
fname2 = ["1_","1_","2_","2_","2_","3_","3_","3_","3_","5_","5_","5_","5_"]
window_names = ['695','725','755','785','815','845','875','905']
#fname3 = [63,69,55,61,67,9,15,21,27,33,30,36,42,48]
fname3 = [63,69,55,61,67,9,15,21,27,30,36,42,48]
window_indexs = arange(8)
names = ['dec','freq','mask_T_I','ra','T_I']
#nfiles = 14.0
nfiles = 13.0
strips = array([0,1,2,3,4,5])
window_index = 0
#window_name = window_names[0]

for window_name in window_names:
  for strip in strips:
    T0 = load((fname1 + fname2[0]  + window_name + '_' + str(fname3[0]+strip) + '_'+str(strip)+'_' + names[4]+'.npy'))
    mask0 = load((fname1 + fname2[0]  + window_name + '_' + str(fname3[0]+strip) + '_'+str(strip)+'_' + names[2]+'.npy'))
    count = (~mask0).astype(int)
    ra0 = load((fname1 + fname2[0]  + window_name + '_' + str(fname3[0]+strip) + '_'+str(strip)+'_' + names[3]+'.npy'))
    dec0 = load((fname1 + fname2[0]  + window_name + '_' + str(fname3[0]+strip) + '_'+str(strip)+'_' + names[0]+'.npy'))
    freq0 = load((fname1 + fname2[0]  + window_name + '_' + str(fname3[0]+strip) + '_'+str(strip)+'_' + names[1]+'.npy'))
    T_avg = T0
    mask = mask0
    ra = ra0
    dec = dec0
    freq = freq0
    for i in arange(1,int(nfiles)):
      print window_name, strip, i
      mask_temp = load((fname1 + fname2[i]  + window_name + '_' + str(fname3[i]+strip) + '_'+str(strip)+'_' + names[2]+'.npy'))
      T_temp = load((fname1 + fname2[i]  + window_name + '_' + str(fname3[i]+strip) + '_'+str(strip)+'_' + names[4]+'.npy'))
      ra_temp = load((fname1 + fname2[i]  + window_name + '_' + str(fname3[i]+strip) + '_'+str(strip)+'_' + names[3]+'.npy'))
      dec_temp = load((fname1 + fname2[i]  + window_name + '_' + str(fname3[i]+strip) + '_'+str(strip)+'_' + names[0]+'.npy'))
      freq_temp = load((fname1 + fname2[i]  + window_name + '_' + str(fname3[i]+strip) + '_'+str(strip)+'_' + names[1]+'.npy'))
      count = count + (~mask_temp).astype(int)
      #if ( (mask == False) & (mask_temp==False)):
      maskand = (mask == False) & (mask_temp==False)
      T_avg[maskand] = T_avg[maskand] + T_temp[maskand]
      #elif ( (mask==True) & (mask_temp == False) ):
      maskor = (mask==True) & (mask_temp == False)
      T_avg[maskor] = T_temp[maskor]
      mask = (mask & mask_temp)
      ra = ra+ra_temp
      dec = dec+dec_temp
      freq = freq+freq_temp

    T_avg = T_avg/count
    ra = ra/nfiles
    dec = dec/nfiles
    freq = freq/nfiles
    save((fnameout +window_name+'_'+str(strip)+'_'+names[4]+'.npy'), T_avg)
    save((fnameout +window_name+'_'+str(strip)+'_'+names[2]+'.npy'), mask)
    save((fnameout +window_name+'_'+str(strip)+'_'+names[3]+'.npy'), ra)
    save((fnameout +window_name+'_'+str(strip)+'_'+names[0]+'.npy'), dec)
    save((fnameout +window_name+'_'+str(strip)+'_'+names[1]+'.npy'), freq)
