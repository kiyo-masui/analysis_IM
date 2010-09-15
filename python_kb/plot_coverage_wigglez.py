import pyfits
import pylab
import azel2radec

file2 = pyfits.open('wigglez_dr1.fits')
wigglez_data = file2[1].data
rawigglez = wigglez_data.field('RA')
decwigglez = wigglez_data.field('DEC')
#maskra = (rawigglez > 320 ) & ( rawigglez < 330 )  # 22hr wigglez
maskra = (rawigglez > 10 ) & ( rawigglez < 20 )   # 1hr wigglez

pylab.plot(rawigglez[maskra] , decwigglez[maskra], 'go')
pylab.plot(azel2radec.rag_d , azel2radec.decg_d, 'r-', hold=True)
pylab.savefig('wigglez1hr_last_azel.pdf')
pylab.show()
