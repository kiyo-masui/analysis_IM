import numpy as np
import scipy as sp
import scipy.optimize
import matplotlib
matplotlib.use('Agg')
import pylab
from matplotlib.backends.backend_pdf import PdfPages
from scripts import map_zoom as mz
import matplotlib.pyplot as plt
import core.algebra as al
import numpy.ma as ma
#from scripts import parkes_corr_maps as cm
from map import beam

# list of [ra, dec (degrees), [flux in mJy, wavelength in cm]]

sources = {'ra30decn31': [[30.05125, -30.88917, [1280.00, 6.], [2900.00, 20.], [2280.000, 11.]],
    [30.05, -30.89083, [1440.00, 6],[3600.00, 20],[2370.00,11.1]]],
 'ra28decn30': [[27.64167, -29.53333, [810.00, 6.],[2900.00, 20.],[1510.000, 11.1]]],
 'ran10decn27': [[-10.02125, -27.46694,[1300.00,6.],[3230.00,20.],[1930.00,11.]],
    [-10.02125, -27.46667,[1260.00,6.],[2900.00, 20.],[1910.000,11.1]]],
 'ra187dec2': [[187.27708, 2.04778, [44940.0, 6.],[50100., 20.],[65243., 82.]], [187.2775, 2.0525, [40442.50, 6.],[42212.50, 20.],[38970., 11.1]],
[187.27792, 2.0525, [36700., 6.],[42000., 20.],[40900., 11.1]]],
 #'ra230dec5':[[230.30958, 4.50583, [1030., 6.],[4500., 20.],[2300., 11.1]],
#[230.30625, 4.50722, [1246., 6.],[4384., 20.],[3474., 82.]],
#[230.30958, 4.50583, [1060., 6.],[4345., 20.],[2365., 82.]]],
 'ra190decn5': [[190.58083, -4.77167, [975., 6.],[3533.33, 20.],[1860., 11.]]],
 'ra165dec2': [[164.62292, 1.56639, [3216., 6.],[3465., 20.],[2850., 11.]]],
 'ra78decn0': [[177.6825, -0.39833, [1872.5, 6.],[3197.5, 20.],[2495., 11.]],
               [177.68292, -0.39833, [1900., 6.],[3100., 20.],[2400., 11.1]]],
 'ra209dec1': [[209.25625, 1.07806,[710.,6.],[2300.,20.],[1350.,11.1]]],
 'ra228dec1': [[228.1075, 1.34333,[858.,6.],[2275.,20.],[6733.,82.]],
               [228.10583, 1.3525,[680.,6.],[2220.,20.],[1280.,11.1]]],
 'ra231dec3': [[231.45417, 3.14055,[984.,6.],[2178.,20.],[4098.,82.]]],
 #'ra181decn4': [[181.00833, -4.38167,[970.,6.],[2000.,20.],[1320.,11.1]]],
 'ran32decn28': [[-31.98417, -28.47444,[1320.,6.],[2800.,20.],[2000.,11.1]],
                 [-31.98458, -28.475,[1360.,6.],[2500.,20.],[2000.,11]]]
}

def fit_pl_lsq(fluxes, freqs):
    #Will fit to power law in linear units, rather than to line in log units
    #Variables fluxes and freqs should be numpy arrays.
    def model(pars):
        #pars[0] is the amplitude at 1400 MHz.
        #pars[1] is the power law index
        model = pars[0]*(freqs/1400.)**pars[1]
        return model
    init_index = np.log(fluxes[0]/fluxes[-1])/np.log(freqs[0]/freqs[-1])
    init_amp = fluxes[0]*(1400./freqs[0])**init_index
    #print init_amp
    #print init_index
    init_pars = [init_amp, init_index]
    #print model(init_pars)
    fit_pars, err = scipy.optimize.leastsq(lambda q: model(q) - fluxes, init_pars, ftol=1e-5, xtol=1e-5)
    return [fit_pars, err]

def fit_index_given(fluxes, freqs, index):
    #fluxes in K, freqs in MHz
    def model(amp):
        model = amp*(freqs/1400.)**index
        return model
    init_amp = fluxes[0]*(1400./freqs[0])**index
    fit_pars, err = scipy.optimize.leastsq(lambda q: model(q) - fluxes, init_amp, ftol=1e-5, xtol=1e-5)
    return [fit_pars, index, err]

def plot_fluxes(fluxes, freq, labels, fname):
    for pair in zip(fluxes, labels):
        pylab.plot(freq, pair[0], label = pair[1])
    pylab.rc('font', family='serif', size=10)
    pylab.xlabel('Frequency')
    pylab.ylabel('Temp in K')
    pylab.legend()
    pylab.savefig(fname + '.png')
    pylab.clf()

def plot_page(fluxes, freqs, labels, title, pdf, units = 'Jy'):
    for triple in zip(fluxes, freqs, labels):
        pylab.plot(triple[1], triple[0], label = triple[2])
    pylab.rc('font', family='serif', size=10)
    pylab.xlabel('Frequency')
    pylab.ylabel('Flux in ' + units)
    pylab.legend()
    pylab.title(title)
    pdf.savefig()
    pylab.clf()

def plot_map(map, title, pdf, f_avg = True, white_nohits= True, convolve = False,f_range = [0,64]):
    if convolve:
        parkes_res = degrade_parkes(1.)
        map = parkes_res.apply(al.make_vect(map))
    map = map[f_range[0]:f_range[1],:,:]
    if white_nohits:
        map[map == 0.0] = np.nan
    if f_avg:
        data = np.transpose(np.mean(map, axis=0))
        ra_step = map.info['ra_delta']
        dec_step = map.info['dec_delta']
        ra_range = [map.info['ra_centre'] + ra_step*data.shape[1]/2, map.info['ra_centre'] - ra_step*data.shape[1]/2]
        dec_range = [map.info['dec_centre'] + dec_step*data.shape[0]/2, map.info['dec_centre'] - dec_step*data.shape[0]/2]
        data = np.fliplr(data)
        data = np.flipud(data)
        im = pylab.imshow(data, extent = [ra_range[0], ra_range[1], dec_range[0], dec_range[1]], interpolation='nearest')
        pylab.colorbar(im)
        pylab.title(title)
        pdf.savefig()
        pylab.clf()

def degrade_parkes(factor):
    freq_data = sp.array([1283.5, 1300, 1325, 1350, 1430], dtype=float)
    beam_data = sp.array([14.4, 14.4, 14.4, 14.4, 14.4])/60.
    beam_data = beam_data*1420./freq_data
    freq_data *= 1.0e6
    beam_diff = sp.sqrt(max(factor*beam_data)**2-(beam_data)**2)
    common_resolution = beam.GaussianBeam(beam_diff, freq_data)
    return common_resolution

def make_submap_dict(sources, input_maps, sub_size):
    sub_maps = {}
    for key in sources.keys():
        catalogs = sources[key]
        #Each catalog entry stored by the same key should have nearly identical
        #right ascension and declination, so just use the first entry.
        pos = [catalogs[0][0], catalogs[0][1]]
        for map in input_maps:
            if source_in_map(pos, map):
                print pos
                zoom = mz.zoom_cube(map, [1315500000.,pos[0],pos[1]], sub_size)
                sub_maps[key + '_from_ra_centre_' + str(map.info['ra_centre'])]=zoom
    return sub_maps

def source_in_map(position, map):
    #returns True if source is in map, False if not
    #position[0] should be ra, postion[1] should be dec
    info = map.info
    ra_delta = info['ra_delta']
    dec_delta = info['dec_delta']
    ra_centre = info['ra_centre']
    dec_centre = info['dec_centre']
    ra_range = np.sort(np.array([ra_centre - ra_delta*map.shape[1]//2, ra_centre + ra_delta*map.shape[1]//2]))
    dec_range = np.sort(np.array([dec_centre - dec_delta*map.shape[2]//2, dec_centre + dec_delta*map.shape[2]//2]))
    print ra_range
    print dec_range
    ra_bool = position[0] >= ra_range[0] and position[0] <= ra_range[1]
    dec_bool = position[1] >= dec_range[0] and position[1] <= dec_range[1] 
    bool = ra_bool and dec_bool
    return bool

def plot_zooms(zoom_dict, filename, f_avg=True, convolve = False, f_range = [0,64]):
    #Will average over freqs if f_avg=True.
    with PdfPages(filename + '.pdf') as pdf:
        #for key in zoom_dict.keys():
        for key in [ el[0] for el in sorted(zoom_dict.iteritems(), key = lambda x: x[1].info['ra_centre'])]:
            plot_map(zoom_dict[key], key, pdf, f_avg = f_avg, convolve = convolve, f_range = f_range)

def plot_sources(sources, filename, fit=True, type = 'catalog', units = 'Jy', convert_to_temp = False, fit_only = False, freq_edges = [0,64], fit_cat = None):
    #Will plot the fluxes of every source in sources dictionary,
    #and then will plot the power law fit to that flux.
    #Returns dictionary with same key as source, containing power law fits.
    #If convert_to_temp is true, the index gains an additional -2.0
    output = {}
    with PdfPages(filename + '.pdf') as pdf:
        #for key in sources.keys():
        for key in [ el[0] for el in sorted(sources.iteritems(), key = lambda x: x[1])]:
            output[key]=[]
            version = sources[key]
            fluxes = []
            freqs = []
            labels = []
            for el in version:
                ra = el[0]
                dec = el[1]
                if type == 'catalog':
                    print version
                    flux = np.array([el[2][0],el[3][0],el[4][0]])/1000.
                    freq = 1420.*21.*np.array([el[2][1]**-1,el[3][1]**-1,el[4][1]**-1])
                    order = np.argsort(freq)
                    flux = flux[order]
                    freq = freq[order]
                    if not fit_only:
                        fluxes.append(flux)
                        freqs.append(freq)
                        labels.append('ra ' + str(ra) + ', dec ' + str(dec))
                else:
                    #Must be real freq. spectrum from a map zoomed
                    #on a point source.
                    flux = np.array(el[2])
                    freq = 1315.5 - 32.0 + np.array(range(64))
                    freq = freq[freq_edges[0]:freq_edges[1]]
                    if not fit_only:
                        fluxes.append(flux)
                        freqs.append(freq)
                        labels.append('ra ' + str(ra) + ', dec ' + str(dec))
                if fit:
                    the_fits = fit_pl_lsq(flux, freq)
                    amp = the_fits[0][0]
                    index = the_fits[0][1]
                    freq_list = np.linspace(freq[0], freq[-1], num=20)
                    if convert_to_temp:
                        amp = amp*1315.5**2
                        index -= 2.0
                        freq_list = 1315.5 - 32.0 + np.array(range(64))
                        freq_list = freq_list[freq_edges[0]:freq_edges[1]]
                    freqs.append(freq_list)
                    fluxes.append(model([amp,index],freq_list))
                    labels.append('ra ' + str(ra) + ', dec ' + str(dec) + ', index ' + str(round(index, 3)))
                    output[key].append([fluxes, freqs, labels]) 
                if type == 'map' and fit_cat != None:
                    #Plot fit using catalog index
                    #Will use index of first entry.  Add -2.0 to go to K.
                    dat = fit_cat[key[0:key.find('_')]][0]
                    flux_cat = np.array([dat[2][0],dat[3][0],dat[4][0]])/1000.
                    freq_cat = 1420.*21.*np.array([dat[2][1]**-1,dat[3][1]**-1,dat[4][1]**-1])
                    index = fit_pl_lsq(flux_cat, freq_cat)[0][1]
                    freq = 1315.5 - 32.0 + np.array(range(64))
                    freq = freq[freq_edges[0]:freq_edges[1]]
                    the_fit = fit_index_given(el[2], freq, index-2.0)
                    fluxes.append(model([the_fit[0][0],index-2.0],freq))
                    freqs.append(freq)
                    labels.append('index ' + str(index-2.0) + ' fit from catalog')
            plot_page(fluxes, freqs, labels, key, pdf, units = units)
    return output

def model(pars, freqs):
    #pars[0] is the amplitude at 1400 MHz.
    #pars[1] is the power law index
    model = pars[0]*(freqs/1400.)**pars[1]
    return model

def fit_gauss_const(map, width):
    #Fits each freq. slice of map to a Gaussian at the center plus a constant.
    return map

if __name__ == '__main__':
    #First, plot the spectrum and power law fit for each source.
    #freq_edges = [10,59]
    cat_fits = plot_sources(sources, 'catalog_sources_temp', convert_to_temp=True, fit_only = True)
    #Note, calibration is a factor of 2 too high in this folder.
    map_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/cal_nobp_jansky/combined/I/'
    
    map27_dir ='/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_bp_divide/hitconv_sync27_mbcal/I/'

    map20_dir = '/scratch2/p/pen/nluciw/parkes/maps/bandpass/I/'

    map_nobp_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_to_share/'

    map_meanbp_dir = '/scratch2/p/pen/nluciw/parkes/maps/meanbandpass/I/'

    svd_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_to_share/only_good_beams/conv_fixed/'

    svd_subs = ['ran18', 'ra33', 'ra165', 'ra182', 'ra199', 'ra216']

    sub_dirs = ['ran18', 'ra33', 'ra165','ra199','ra216']
    maps = []
    maps27 = []
    maps20 = []
    maps_meanbp = []
    maps_nobp = {'XX':{'beam1':[],
                 'beam6':[]},
                 'YY':{'beam1':[],
                 'beam6':[]}}
    labels = [['ran18decn30','488by106'],['ra33decn30','568by106'],
              ['ra216dec0','440by136'],['ra199dec0','440by136'],
              ['ra165dec0','440by136']]
    beginning = 'parkes_parallel_thread_'
    middle = '_p08_'
    end = '_beam_all__clean_map_bp_div_I_1316.npy'
    #nobp_mid = '_beam'
    nobp_end_xx = '_clean_map_XX_1316.npy'
    nobp_end_yy = '_clean_map_YY_1316.npy'
    #Get common beam resolution
    parkes_res = degrade_parkes(1.01)
    submaps_nobp = {'XX': {}, 'YY':{}}

    mode0_maps = {}
    for sub in svd_subs:
        mode0_maps[sub+'_XX'] = al.load(svd_dir + 'XX/' + sub + '/' + 'sec_beam1_modes_clean_map_I_with_beam2_15modes.npy')
    plot_zooms(mode0_maps, 'mode1_XX_beam1', f_range = [1,2])

    print 'plotting SVD mode 0.'

    for el in labels:
        #Multiply by 0.5 to fix calibration.
        #maps.append(0.5*al.load(map_dir + beginning + el[0] + middle + el[1] + end))
        #temp = al.load(map_nobp_dir + beginning + el[0] + middle + el[1] + '_' + 'beam6' + nobp_end_xx)
        #temp = al.load(map_nobp_dir + beginning + el[0] + middle + el[1] + '_' + 'beam1' + nobp_end_xx)
        for beam in range(1,14):
            for end in [nobp_end_xx, nobp_end_yy]: 
                temp = al.load(map_nobp_dir + beginning + el[0] + middle + el[1] + '_' + 'beam'+ str(beam) + end)
                #Go to common beam resolution 
                temp = parkes_res.apply(al.make_vect(temp))
                key1 = end[11:13]
                key2 = beam
                if not key2 in maps_nobp[key1]:
                    maps_nobp[key1][key2] = [] 
                #maps_nobp['XX']['beam6'].append(temp)
                #maps_nobp['XX']['beam'+str(beam)].append(temp)
                maps_nobp[key1][key2].append(temp)
                submaps_nobp[end[11:13]][beam] = make_submap_dict(sources,maps_nobp[end[11:13]][beam],[64000000.,0.5,0.5])
    for dir in sub_dirs:
        maps27.append(al.load(map27_dir + dir + '/' + 'combined_clean_map_0modes.npy'))
        maps20.append(al.load(map20_dir + dir + '/' + 'combined_clean_map_0modes.npy'))
        maps_meanbp.append(al.load(map_meanbp_dir + dir + '/' + 'combined_clean_map_0modes.npy'))
        #temp = al.load(map_nobp_dir + beginning + el[0] + middle + el[1] + '_' + 'beam1' + nobp_end_xx)
        #Go to common beam resolution 
        
        #maps_nobp['XX']['beam1'].append(al.load(map_nobp_dir + beginning + el[0] + middle + el[1] + '_' + 'beam1' + nobp_end_xx))
    #submaps = make_submap_dict(sources,maps27,[64000000.,1.2,1.2])
    #submaps_nobp = make_submap_dict(sources,maps_nobp['XX']['beam6'],[64000000.,2.2,2.2])
    #submaps = make_submap_dict(sources,maps20,[64000000.,0.5,0.5])
    #submaps_nobp = make_submap_dict(sources,maps_nobp['XX']['beam1'],[64000000.,0.5,0.5])
    #submaps_nobp = make_submap_dict(sources,maps_nobp['XX']['beam1'],[64000000.,1.2,1.2])
    #plot_zooms(submaps, 'test_submaps_point5zoom')
    #plot_zooms(submaps_nobp, 'beam6_XX_submaps_commonbeam', convolve = False)
    #plot_zooms(submaps_nobp, 'beam1_XX_submaps_commonbeam', convolve = False, f_range = [10,59])
    #submaps = make_submap_dict(sources,maps27,[64000000.,0.5,0.5])
    #submaps = make_submap_dict(sources,maps_meanbp,[64000000.,0.5,0.5])
    #submaps = make_submap_dict(sources,maps20,[64000000.,0.5,0.5])
    the_submaps = submaps_nobp

    #Now, make spectra and fits for each source in the index -2.7 maps.
    source_spectra = {}
    #for key in submaps.keys():

    #freq_edges = [0, 63]
    freq_edges = [10, 59]

    bp_est_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/bp_estimates/pt_source_catalog_estimates/'

    for key1 in the_submaps.keys():
        for key2 in the_submaps[key1].keys():
            submaps = the_submaps[key1][key2]
            for key in [ el[0] for el in sorted(submaps.iteritems(), key = lambda x: x[1].info['ra_centre'])]:
                map = submaps[key]
                #map = map[freq_edges[0]:freq_edges[1]]
                hit_bool = map == 0.0
                hit_bool = np.logical_or(hit_bool, np.isnan(map))
                ma_map = ma.MaskedArray(map, mask = hit_bool)
                #Optionally, subtract minimum from each frequency slice
                #temp_min = ma.min(map, axis = 1)
                #t_map = ma.transpose(ma_map) - ma.min(ma.MaskedArray(temp_min, mask = np.isnan(temp_min)), axis = 1)
                t_map = ma.transpose(ma_map)
                flux = ma.mean(ma.mean(t_map, axis =0), axis=0)
                cat_pl = cat_fits[key[0:key.find('_')]]
                source_spectra[key] = [[map.info['ra_centre'], map.info['dec_centre'],flux]]
                source_spectra[key][0][2] /= cat_pl[0][0][0]
                source_spectra[key][0][2] /= ma.mean(source_spectra[key][0][2])
            plot_sources(source_spectra, bp_est_dir + str(key1) + 'map' + 'beam' + str(key2) + '_bpestimate_from_source_spectra_width_0point5', fit = False, type = 'map', units = 'Normalized')
      

    for key in [ el[0] for el in sorted(submaps.iteritems(), key = lambda x: x[1].info['ra_centre'])]:
        print key
        map = submaps[key]
        map = map[freq_edges[0]:freq_edges[1]]
        print map.shape
        #print map[0,:,:]
        #print map[20,:,:]
        #parkes_res = degrade_parkes(1.)
        #map = parkes_res.apply(al.make_vect(map))
        #Mask any values equal to exactly zero or nan, which means no hits.
        hit_bool = map == 0.0
        hit_bool = np.logical_or(hit_bool, np.isnan(map))
        #print hit_bool[0,:,:]
        #print hit_bool[20,:,:]
        ma_map = ma.MaskedArray(map, mask = hit_bool)
        
        #Optionally, subtract minimum from each frequency slice
        #temp_min = ma.min(map, axis = 1)
        #t_map = ma.transpose(ma_map) - ma.min(ma.MaskedArray(temp_min, mask = np.isnan(temp_min)), axis = 1)

        t_map = ma.transpose(ma_map)
        #print t_map
        #print t_map[:,:,0]
        #print t_map.shape
        flux = ma.mean(ma.mean(t_map, axis =0), axis=0)
        print flux.shape
        #print flux
        #print flux
        #cat_pl = cat_fits[key[0:key.find('_')]]
        #print cat_pl
        source_spectra[key] = [[map.info['ra_centre'], map.info['dec_centre'],flux]]
        print source_spectra[key][0][2].shape
        #print cat_pl[0][0]
        #source_spectra[key][0][2] /= cat_pl[0][0][0]
        #source_spectra[key][0][2] /= ma.mean(source_spectra[key][0][2])

        #print source_spectra[key][0][2]

        #print flux.shape
    #plot_sources(source_spectra, 'maps_meanbp_source_spectra_width_0point5_fcut', type = 'map', units = 'K', freq_edges = freq_edges)
    plot_sources(source_spectra, 'maps_meanbp_source_spectra_width_0point5_fcut_cat', type = 'map', units = 'K', freq_edges = freq_edges, fit_cat = sources)
    #plot_sources(source_spectra, 'map20_source_spectra', type = 'map', units = 'K')
    #plot_sources(source_spectra, 'map_meanbp_source_spectra', type = 'map', units = 'K')
    #plot_sources(source_spectra, 'mapXXbeam6_bpestimate_from_source_spectra', fit = False, type = 'map', units = 'K')
    #plot_sources(source_spectra, 'mapXXbeam1_bpestimate_from_source_spectra_width_0point5', fit = False, type = 'map', units = 'K')
