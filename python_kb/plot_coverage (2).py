import pyfits
import pylab
import azel2radec

file2 = pyfits.open('l:/zcat.dr3.v1_0.uniq.fits')
DEEP2_data = file2[1].data
radeep = DEEP2_data.field('RA')
decdeep = DEEP2_data.field('DEC')
#maskra = (radeep > 350 ) & ( radeep < 360 )  # field 3
maskra = (radeep > 35 ) & ( radeep < 40 )   # field 4

pylab.plot(radeep[maskra] , decdeep[maskra], 'go')
pylab.plot(azel2radec.rag_d , azel2radec.decg_d, 'ro', hold=True)
pylab.plot(azel2radec.ra_elaz_d, azel2radec.dec_elaz_d, 'b.', hold=True)
pylab.xlabel('RA (Deg)')
pylab.ylabel('DEC (deg)')
pylab.title('Field 4 Coverage plot')
pylab.show()