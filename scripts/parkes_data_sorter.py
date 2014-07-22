#import fnmatch
#import os
import numpy as np
import pyfits
import codecs
from operator import itemgetter

list = open('datalist_2012_yc.txt', 'r').read().splitlines()

newlist = []
for item in list:
    hdu_list = pyfits.open(item)
    data = hdu_list[1]
    newlist.append((item, data.data['date-obs'][0], np.average(data.data['time'])))
a = sorted(newlist, key=itemgetter(2))
b = sorted(a, key=itemgetter(1))


thefile = open('sorted_datalist_2012_yc.txt', 'w')
#print list
#print newlist
thefile.write('\n'.join(item[0] for item in b))
thefile.close()
