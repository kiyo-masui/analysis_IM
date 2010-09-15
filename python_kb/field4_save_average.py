from numpy import *
import pylab

#fname1 = "/home/scratch/kbandura/GBT08B_037/np_data/field4_"
fname1 = "field4_"
window_names = ['695','725','755','785','815','845','875','905']
names = ['dec','freq','mask_T_I','ra','T_I']
strips = array([0,1,2,3,4,5])

#window_name = window_names[3]
#strip = strips[0]
T_plot = ones((299,8*1230))
mask_plot = T_plot > 2
freq_plot = ones(8*1230)
for strip in strips:
 T_plot = ones((299,8*1230))
 mask_plot = T_plot > 2
 freq_plot = ones(8*1230)
 for i,window_name in enumerate(window_names): 
  T_avg = load((fname1 +window_name+'_'+str(strip)+'_'+names[4]+'.npy'))
  mask = load((fname1 +window_name+'_'+str(strip)+'_'+names[2]+'.npy'))
  ra = load((fname1 +window_name+'_'+str(strip)+'_'+names[3]+'.npy'))
  dec = load((fname1 +window_name+'_'+str(strip)+'_'+names[0]+'.npy'))
  freq = load((fname1 +window_name+'_'+str(strip)+'_'+names[1]+'.npy'))
  cut = (T_avg > 2) | (T_avg < -2)
  mask = mask | cut
  T_avg[mask] = T_avg[~mask].min()
  print "T_avg",T_avg[~mask].min()
  print "<-3",(T_avg[~mask] < -3.0).sum()
  print "> 3",(T_avg[~mask] > 3.0).sum()
  T_plot[:,i*1230:(i*1230+1230)] = T_avg[:,409:-409][:,::-1]
  mask_plot[:,i*1230:(i*1230+1230)] = mask[:,409:-409][:,::-1]
  freq_plot[i*1230:i*1230+1230] = freq[409:-409][::-1]

 print T_plot[~mask_plot].std()
 save((fname1+'all_windows_' + str(strip) +'_'+names[0]+ '.npy'), dec)
 save((fname1+'all_windows_' + str(strip) +'_'+names[1]+ '.npy'), freq_plot)
 save((fname1+'all_windows_' + str(strip) +'_'+names[2]+ '.npy'), mask_plot)
 save((fname1+'all_windows_' + str(strip) +'_'+names[3]+ '.npy'), ra)
 save((fname1+'all_windows_' + str(strip) +'_'+names[4]+ '.npy'), T_plot) 
 #pylab.imshow(T_plot, cmap=pylab.cm.hot, extent=(freq_all.min(),freq_all.max(), ra.max(),ra.min()), aspect=50000000)
 #pylab.xlabel("Frequency")
 #pylab.ylabel("Right Ascension")
 #pylab.show()
 #pylab.savefig("T_plot_"+str(strip)+".png", dpi=300)
