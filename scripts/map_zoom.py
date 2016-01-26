#Script for zooming in on requested map region
# Also including the Abell cluster locations that overlap with Parkes fields here, since I see no better place for it.

'''
|name      |ra     |dec   |_count|bmtype|redshift|rich|dist|vmag |
|ABELL1142S|23 41.3|-30 13|    20|II-III|   0.070|0   |4   | 15.7|
            -18.7
|ABELL1238|11 23.0|+01 05|    63|III    |   0.072|1   |4   | 16.0|
|ABELL1650|12 58.8|-01 45|   114|I-II   |   0.084|2   |5   | 17.0|
|ABELL1750|13 30.9|-01 50|    40|II-III:|   0.086|0   |4   | 15.9|
|ABELL912 |10 01.2|-00 06|    36|       |   0.089|0   |4   | 15.9|
|ABELL1399|11 51.2|-03 05|    82|III    |   0.091|2   |4   | 16.0|
'''

import core.algebra as al
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

f0 = 1420405751.786

def make_abell_list():
    # Returns freq, ra, dec of Abell clusters seen in HI
    clusters = []
    clusters.append([z_to_f(0.070), ra_to_deg(0, -18.7), to_deg(-30, 13) ])
    clusters.append([z_to_f(0.072), ra_to_deg(11, 23.0), to_deg(1, 5) ])
    clusters.append([z_to_f(0.084), ra_to_deg(12, 58.8), to_deg(-1, 45) ])
    clusters.append([z_to_f(0.086), ra_to_deg(13, 30.9), to_deg(-1, 50) ])
    clusters.append([z_to_f(0.089), ra_to_deg(10, 1.2), to_deg(0, -6) ])
    clusters.append([z_to_f(0.091), ra_to_deg(11, 51.2), to_deg(-3, 5) ])
    return clusters

def ra_to_deg(hour, min):
    return hour*15 + 15*min/60.

def to_deg(deg, arcmin):
    return deg + math.copysign(1, deg)*arcmin/60.

def z_to_f(z):
    return f0/(z+1)

def f_to_z(f):
    return f0/f - 1

def zoom_cube(map, center, size):
    #Returns new map centered at given center, with given size
    #map: info_array object, dimensions (freq, ra, dec)
    #center: list with freq, ra, dec coordinates
    #size: list with cube size in freq, ra, dec
    info = map.info
    delta = []
    delta.append(info['freq_delta'])
    delta.append(info['ra_delta'])
    delta.append(info['dec_delta'])
    pix_center = []
    pix_center.append(round((center[0] - info['freq_centre'])/delta[0]) + map.shape[0]//2)
    pix_center.append(round((center[1] - info['ra_centre'])/delta[1]) + map.shape[1]//2)
    pix_center.append(round((center[2] - info['dec_centre'])/delta[2]) + map.shape[2]//2)
    range = [round(pair[0]/(abs(pair[1])*2)) for pair in zip(size, delta)]
    print pix_center
    print range
    orig_shape = map.shape
    map = map[pix_center[0]-range[0]: pix_center[0]+range[0]+1, pix_center[1]-range[1]: pix_center[1]+range[1]+1, pix_center[2]-range[2]: pix_center[2]+range[2]+1 ]
    map.info['freq_centre'] = info['freq_centre'] + (pix_center[0]-orig_shape[0]//2)*delta[0]
    map.info['ra_centre'] = info['ra_centre'] + (pix_center[1]-orig_shape[1]//2)*delta[1]
    map.info['dec_centre'] = info['dec_centre'] + (pix_center[2]-orig_shape[2]//2)*delta[2]
    return map 
   
def map_pdf(map, z_clust, name):
    f_center = map.info['freq_centre']
    df = map.info['freq_delta']
    pdf = PdfPages(name + '.pdf')
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 2}

    matplotlib.rc('font', **font)
    for f in xrange(map.shape[0]):
        data = np.transpose(map[f,:,:])
        ra_step = map.info['ra_delta']
        dec_step = map.info['dec_delta']
        ra_range = [map.info['ra_centre'] + ra_step*data.shape[1]/2, map.info['ra_centre'] - ra_step*data.shape[1]/2]
        dec_range = [map.info['dec_centre'] + dec_step*data.shape[0]/2, map.info['dec_centre'] - dec_step*data.shape[0]/2]
        data = np.fliplr(data)
        data = np.flipud(data)
        plt.figure(figsize=(2, 2))
        plt.imshow(data, extent = [ra_range[0], ra_range[1], dec_range[0], dec_range[1]])
        freq = f_center + (f - map.shape[0]//2)*df
        z = f_to_z(freq)
        dz = z - z_clust
        dv = dz*300000
        deg = (np.pi/180)*hubble_d(z)
        plt.title('z ' + trunc(z, 5) + ', ' + trunc(dv, 2) + ' km/s,' + ' 1 deg is ' + trunc(deg,2) + ' h-1Mpc')
        pdf.savefig()
        plt.close()
    pdf.close()

def trunc(num, digits):
   sp = str(num).split('.')
   return '.'.join([sp[0], sp[1][:digits]])

def hubble_d(z):
   # in h-1Mpc
   return 300000*z/100
