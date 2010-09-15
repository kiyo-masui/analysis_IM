from numpy import *

#fname1 = "/home/scratch/kbandura/GBT08B_037/np_data/field4_"
fname1 = "field4_"
window_names = ['695','725','755','785','815','845','875','905']
names = ['dec','freq','mask_T_I','ra','T_I']
strips = array([0,1,2,3,4,5])

strip = 0
dec1 = load((fname1+'all_windows_' + str(strip) +'_'+names[0]+ '.npy'))
freq1 = load((fname1+'all_windows_' + str(strip) +'_'+names[1]+ '.npy'))
mask1=load((fname1+'all_windows_' + str(strip) +'_'+names[2]+ '.npy'))
ra1=load((fname1+'all_windows_' + str(strip) +'_'+names[3]+ '.npy'))
T1=load((fname1+'all_windows_' + str(strip) +'_'+names[4]+ '.npy')) 

dec = ones((strips.shape[0],dec1.shape[0]))
ra = ones((strips.shape[0],ra1.shape[0]))
freq = ones((strips.shape[0],freq1.shape[0]))
T =  ones((strips.shape[0],T1.shape[0],T1.shape[1]))
mask = T > 2

for strip in strips:
  dec[strip] = load((fname1+'all_windows_' + str(strip) +'_'+names[0]+ '.npy'))
  freq[strip] = load((fname1+'all_windows_' + str(strip) +'_'+names[1]+ '.npy'))
  mask[strip]=load((fname1+'all_windows_' + str(strip) +'_'+names[2]+ '.npy'))
  ra[strip]=load((fname1+'all_windows_' + str(strip) +'_'+names[3]+ '.npy'))
  T[strip]=load((fname1+'all_windows_' + str(strip) +'_'+names[4]+ '.npy'))

