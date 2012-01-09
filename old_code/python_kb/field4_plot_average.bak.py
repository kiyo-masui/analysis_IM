from numpy import *
import pylab

fname1 = "/home/scratch/kbandura/GBT08B_037/np_data/field4_"
window_names = ['695','725','755','785','815','845','875','905']
names = ['dec','freq','mask_T_I','ra','T_I']
strips = array([0,1,2,3,4,5])

#window_name = window_names[3]
strip = strips[0]
T_plot = ones((299,8*1230))
freq_all = ones(8*1230)

for i,window_name in enumerate(window_names): 
  T_avg = load((fname1 +window_name+'_'+str(strip)+'_'+names[4]+'.npy'))
  mask = load((fname1 +window_name+'_'+str(strip)+'_'+names[2]+'.npy'))
  ra = load((fname1 +window_name+'_'+str(strip)+'_'+names[3]+'.npy'))
  dec = load((fname1 +window_name+'_'+str(strip)+'_'+names[0]+'.npy'))
  freq = load((fname1 +window_name+'_'+str(strip)+'_'+names[1]+'.npy'))
  cut = (T_avg > 1) | (T_avg < -1)
  mask = mask | cut
  T_avg[mask] = T_avg[~mask].min()
  print "T_avg",T_avg[~mask].min()
  print "<-3",(T_avg[~mask] < -3.0).sum()
  print "> 3",(T_avg[~mask] > 3.0).sum()
  T_plot[:,i*1230:(i*1230+1230)] = T_avg[:,409:-409]
  freq_all[i*1230:i*1230+1230] = freq[409:-409]
  #pylab.plot(freq,i*ones(freq.shape[0]))

#pylab.show()

#pylab.plot(ra,T_avg[:,950:1050].mean(axis=1))
#pylab.show()

pylab.imshow(T_plot, cmap=pylab.cm.hot, extent=(freq_all.min(),freq_all.max(), ra.max(),ra.min()), aspect=100000000)
pylab.show()
