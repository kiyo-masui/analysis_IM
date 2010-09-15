import pyfits
import pylab
#import azel2radec

file = pyfits.open('/home/scratch/kbandura/GBT08A_081/friday/field3_ralong_165-170_onescan.raw.acs.fits')
gbt_data = file[1].data
ra_gbt = gbt_data.field('crval2')
dec_gbt = gbt_data.field('crval3')

file2 = pyfits.open('zcat.dr3.v1_0.uniq.fits')
DEEP2_data = file2[1].data
radeep = DEEP2_data.field('RA')
decdeep = DEEP2_data.field('DEC')
maskra = (radeep > 350 ) & ( radeep < 360 )  # field 3
#maskra = (radeep > 35 ) & ( radeep < 40 )   # field 4

#pylab.plot(radeep[maskra] , decdeep[maskra], 'go')
#pylab.plot(ra_gbt , dec_gbt, 'ro', hold=True)
#pylab.savefig('field3_ralong.png')
myfile = open('observing_positions.txt','w')
for i in range(len(ra_gbt)):
   myfile.write(str(ra_gbt[i]) + ' ' + str(dec_gbt[i]) + '\n')

myfile.close 
