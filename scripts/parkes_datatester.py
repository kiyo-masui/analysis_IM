import numpy as np
import pyfits

matches  = open('datalist_2012.txt', 'r').read().splitlines()

test = pyfits.open('/mnt/raid-project/gmrt/kiyo/parkes_data/hipsr_2012-10-27_1825-P641_west1_1315_P641.sdfits')
tester = test[1]

differences = []
for file in matches:
    hdu_list = pyfits.open(file)
    data = hdu_list[1]
    differences.append(np.max(tester.data['data'] -data.data['data']))
print np.max(differences)
