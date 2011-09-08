
import rfi4096   # The RFI clean module.

len = 4096   # length of the array, i.e. # of frequencies per IF band.

# cross_array = 1D, 4096 array with the cross correlation r = sqrt(XY^2 + YX^2) / sqrt(XX * YY)
# freq_array = 1D, 4096 array of frequencies.
# fit_array = 1D, 4096 array, holds the fit to the cross correlation.

rfi4096.get_fit4096(cross_array,freq_array,fit_array)


# Some parameters used by clean().

sig = 3.0        # Threshold is 3 sigma.
tol = 3          # Number of adjacent bins to discard.
flat = 1         # 1=yes, you want the cross to be flat.
spike = 1        # 1=yes, throw away missed noise by computing dT/df.
dTdf_limit = 10  # throw away the 10 largest dT/df. Ignored if spike = 0.
dTdf_tol = 2     # Number of adjacent bins to discard. Ignored if spike = 0. 



# Call the clean() function. This returns "mask_array" which is a True/False array.
# 1 means noisy, 0 means clean.
rfi4096.clean4096(sig,tol,flat,spike,dTdf_limit,dTdf_tol,fit_array,cross_array,XX_array,freq_array,mask_array)

# Print the cleaned data. 
for i in range(0,len):
   if mask[i] == 0:
      print freq[i],XX[i]
