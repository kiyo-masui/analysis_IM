
# An example.

from numpy import *
import rfi   # The RFI clean module.

len = 2048


f = open('f2048')
freq = [double(P) for P in f.readlines()]
f.close()

f = open('XX2048')
XX = [double(P) for P in f.readlines()]
f.close()

f = open('YY2048')
YY = [double(P) for P in f.readlines()]
f.close()

f = open('XY2048')
XY = [double(P) for P in f.readlines()]
f.close()

f = open('YX2048')
YX = [double(P) for P in f.readlines()]
f.close()

cross = zeros(len, double)
for i in range(0,len):
   cross[i] = math.sqrt( (XY[i]**2) + (YX[i]**2) ) / math.sqrt( XX[i]*YY[i] )


fit   = zeros(len, double) # holds the piecewise-linear fit to the cross waveform.
mask  = zeros(len, int)    # The boolean(int) mask containing flagged frequencies.

# define some SWIG arrays.
freq_array = rfi.new_doublearray(len)
XX_array = rfi.new_doublearray(len)
cross_array = rfi.new_doublearray(len)
fit_array = rfi.new_doublearray(len)
mask_array = rfi.new_intarray(len)

# fill the SWIG arrays with your data.
# use rfi.doublearray_setitem to insert double.
# use rfi.intarray_setitem to insert int.

for i in range(0,len):
   rfi.doublearray_setitem(freq_array,i,freq[i]);
   rfi.doublearray_setitem(XX_array,i,XX[i]);
   rfi.doublearray_setitem(cross_array,i,cross[i]);

# Make a rough fit to the cross correlation.
# You could fit to the mean of the cross correlation instead of fitting each array.
rfi.get_fit(len,cross_array,freq_array,fit_array)


# Some parameters used by clean().

sig = 3.0        # Threshold is 3 sigma.
tol = 3          # Number of adjacent bins to discard.
flat = 1         # 1=yes, you want the cross to be flat.
spike = 1        # 1=yes, throw away missed noise by computing dT/df.
dTdf_limit = 10  # throw away the 10 largest dT/df. Ignored if spike = 0.
dTdf_tol = 2     # Number of adjacent bins to discard. Ignored if spike = 0. 



# Call the clean() function. This returns "mask_array" which is a True/False array.
# 1 means noisy, 0 means clean.
rfi.clean(len,sig,tol,flat,spike,dTdf_limit,dTdf_tol,fit_array,cross_array,XX_array,freq_array,mask_array)

# use rfi.doublearray_getitem to retrieve double
# use rfi.intarray_getitem to retrieve int
# Filling in a regular array. 
for i in range(0,len):
   mask[i] = rfi.intarray_getitem(mask_array,i)

# Print the noisy data.
for i in range(0,len):
   print freq[i],XX[i]

# Print the cleaned data. 
#for i in range(0,len):
#   if mask[i] == 0:
#      print freq[i],XX[i]
