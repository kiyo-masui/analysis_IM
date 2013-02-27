import fnmatch
import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt

def find_pattern(pattern,root_dir):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

def saveplot(hist, extent, j, path, name):
    currentplot=plt.imshow(hist, vmax=2, vmin=-2, extent=extent, interpolation='nearest')
    plt.colorbar()
    plt.savefig(path + name + '{0:03}'.format(j), bbox_inches=0)
    plt.close()


#matches = find_pattern("*east*.sdfits","/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/rawdata/")
matches = open('final_datalist.txt', 'r').read().splitlines()
matches2 = open('sorted_datalist_2012_ycf.txt', 'r').read().splitlines()
#a=np.arange(len(matches))
#b=str(a.tolist()).strip()
#print b
#print "\n".join(matches)

def fix(list):
    wh_wrap = (list>180)
    list[wh_wrap] = list[wh_wrap] - 360
    return list

def ra_max(list):
    ra_max_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        ra_vec = fix(hdu_data.data['CRVAL3'])
        ra_max_file = np.max(ra_vec)
        ra_max_list.append(ra_max_file)
    max_ra = max(ra_max_list)
    return max_ra

def ra_min(list):
    ra_min_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        ra_vec = fix(hdu_data.data['CRVAL3'])
        ra_min_file = np.min(ra_vec)
        ra_min_list.append(ra_min_file)
    min_ra = min(ra_min_list)
    return min_ra

def dec_max(list):
    dec_max_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        dec_vec = hdu_data.data['CRVAL4']
        dec_max_file = np.max(dec_vec)
        dec_max_list.append(dec_max_file)
    max_dec = max(dec_max_list)
    return max_dec

def dec_min(list):
    dec_min_list = []
    for file in list:
        hdulist = pyfits.open(file)
        hdu_data = hdulist[1]
        dec_vec = hdu_data.data['CRVAL4']
        dec_min_file = np.min(dec_vec)
        dec_min_list.append(dec_min_file)
    min_dec = min(dec_min_list)
    return min_dec

def mapper(x, y, matches, index):
    hitmap = np.zeros((y, x))
    weight_map = hitmap
    ran=ra_min(matches)
    rax=ra_max(matches)
    decn=dec_min(matches)
    decx=dec_max(matches)
    #i=1
    for file in matches:
        for beam in range(1,14):
    #for beam in range(1,14):
        #hitmap = np.zeros((90,630))
        #weight_map = hitmap
        #for file in matches:
            hdulist = pyfits.open(file)
            hdu_data = hdulist[1]
            print hdu_data.data['BEAM']
            beam_mask = hdu_data.data['BEAM'] == beam
            print np.sum(beam_mask)
            ra_vec = fix(hdu_data.data['CRVAL3'][beam_mask])
            dec_vec = hdu_data.data['CRVAL4'][beam_mask]
            #map_data = np.mean(hdu_data.data['data'], axis=1)
            map_data = hdu_data.data['data'][beam_mask,index]
            map_data -= np.mean(map_data)
            #print map_data.shape
            (skycov, dec_edge, ra_edge) = np.histogram2d(dec_vec, ra_vec,                                                                                                                                range=[[decn,decx],[ran,rax]],
                                                     bins=hitmap.shape, weights=None)
            (weighted, dec_edge, ra_edge) = np.histogram2d(dec_vec, ra_vec,
                                    range=[[decn,decx],[ran,rax]],
                                    bins=weight_map.shape,
                                    normed=False,                                                                                                                                   weights=map_data)
            hitmap = hitmap + skycov
            weight_map = weight_map + weighted
        #hitrap = (hitmap < 1)
        #hitmap[hitrap] = hitmap[hitrap]+1
        #map=weight_map/hitmap
        #extent=[ra_edge[0], ra_edge[-1], dec_edge[0], dec_edge[-1]]
        #extent=[25, 35, dec_edge[0], dec_edge[-1]]
        #saveplot(hitmap, extent , beam)
    hitrap = (hitmap < 1)
    hitmap[hitrap] = hitmap[hitrap]+1
    map=weight_map/hitmap
    extent=[ra_edge[0], ra_edge[-1], dec_edge[0], dec_edge[-1]]
    extent2=[40, 55, dec_edge[0], dec_edge[-1]]
    return {'extent':extent, 'map':map, 'extent2':extent2} 
    extent2=[40, 55, dec_edge[0], dec_edge[-1]]
mapper(630, 90, matches, 900)
saveplot(mapper(630, 90, matches, 900)['map'], mapper(630, 90, matches, 900)['extent2'], 1,                                                                       '/cita/h/home-2/anderson/anderson/parkes_analysis_IM/parkes_roughmaps/',                                                                        'partialsky_ueli' )
saveplot(mapper(630, 90, matches, 900)['map'], mapper(630, 90, matches, 900)['extent'], 1,                                                                      '/cita/h/home-2/anderson/anderson/parkes_analysis_IM/parkes_roughmaps/',                                                                        'fullsky_ueli' )
saveplot(mapper(630, 90, matches2, 900)['map'], mapper(630, 90, matches2, 900)['extent'], 1,
                '/cita/h/home-2/anderson/anderson/parkes_analysis_IM/parkes_roughmaps/',
                'fullsky_2012' )
saveplot(mapper(630, 90, matches2, 900)['map'], mapper(630, 90, matches2, 900)['extent2'], 1,
                '/cita/h/home-2/anderson/anderson/parkes_analysis_IM/parkes_roughmaps/',
                'partialsky_2012' )
#saveplot(map,extent, 2)
#saveplot(map, extent , 4)
#start_ind = int(map.shape[1] * (40 - ra_edge[0]) / (ra_edge[-1] - ra_edge[0]))
#end_ind = int(map.shape[1] * (55 - ra_edge[0]) / (ra_edge[-1] - ra_edge[0]))
#plt.imshow(map,  vmax=2, vmin=-2, extent=extent, interpolation='nearest')
#plt.imshow(map[:,start_ind:end_ind],  vmax=2, vmin=-2, extent=extent, interpolation='nearest')
#print weight_map
#print hitmap
#print weight_map.shape
#plt.colorbar()
#plt.show()
