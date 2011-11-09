import scipy as sp
from scipy.interpolate import interp1d
import core.algebra as algebra

def axis_param(axis):
    r"""return min, max and delta assuming unif. spacing"""
    (min_axis, max_axis) = (min(axis), max(axis))
    extent = max_axis-min_axis
    return (min_axis, max_axis, abs(axis[1]-axis[0]), extent, len(axis))


def find_map_region(min_ra, max_ra, min_dec, max_dec, target_sample=0.25,
                    multiplier=16, search_start=0):
    r"""target 0.25 pixels/FWHM"""
    # do multiples of 16 so that 256 freq * 16N is a multiple of 4096

    ra_sampling = target_sample*2.  # initial value that will not fail
    n_ra = search_start
    while ra_sampling > target_sample:
        n_ra += multiplier
        ra_sampling = find_map_dimensions(define_map_region(min_ra, max_ra,
                                                     min_dec, max_dec, n_ra, 32), silent=True)
        ra_sampling = ra_sampling[0]

    dec_sampling = target_sample*2.
    n_dec = search_start
    while dec_sampling > target_sample:
        n_dec += multiplier
        dec_sampling = find_map_dimensions(define_map_region(min_ra, max_ra,
                                           min_dec, max_dec, n_ra, n_dec), silent=True)
        dec_sampling = dec_sampling[1]

    ra_sample_ratio = target_sample/ra_sampling
    dec_sample_ratio = target_sample/dec_sampling

    print "n_ra=%d, samp=%g, ratio=%g" % (n_ra, ra_sampling,
                                          target_sample/ra_sampling)

    print "n_dec=%d, samp=%g, ratio=%g" % (n_dec, dec_sampling,
                                          target_sample/dec_sampling)

    print "original ra=(%g,%g), dec=(%g,%g)" % (min_ra, max_ra,
                                                min_dec, max_dec)

    template_map = define_map_region(min_ra, max_ra, min_dec, max_dec, n_ra, n_dec)
    # now expand the map a bit so that pixels are exactly 0.25 deg
    info = template_map.info

    blank = sp.zeros((256, n_ra, n_dec))
    info['ra_delta'] *= ra_sample_ratio
    info['dec_delta'] *= dec_sample_ratio
    map_prod = algebra.make_vect(blank, axis_names=('freq', 'ra', 'dec'))
    map_prod.info = info

    return map_prod


def define_map_region(min_ra, max_ra, min_dec, max_dec, n_ra, n_dec):
    blank = sp.zeros((256, n_ra, n_dec))

    #axis = (delta*(sp.arange(len) - len//2) + centre)
    ra_delta = (max_ra - min_ra)/float(n_ra-1)
    ra_centre = min_ra + ra_delta*float((n_ra)//2)

    dec_delta = (max_dec - min_dec)/float(n_dec-1)
    dec_centre = min_dec + dec_delta*float((n_dec)//2)

    info = {'ra_delta': ra_delta,
            'dec_delta': dec_delta,
            'dec_centre': dec_centre,
            'axes': ('freq', 'ra', 'dec'),
            'ra_centre': ra_centre,
            'freq_centre': 799609375.0,
            'freq_delta': -781250.0,
            'type': 'vect'}

    map_prod = algebra.make_vect(blank, axis_names=('freq', 'ra', 'dec'))
    map_prod.info = info
    return map_prod


def find_map_dimensions(map_in, silent=False):
    r"""print various aspects of the maps"""
    if isinstance(map_in, str):
        map_in = algebra.make_vect(algebra.load(map_in))

    ra_axis = map_in.get_axis('ra')
    mmd_ra = axis_param(ra_axis)

    dec_axis = map_in.get_axis('dec')
    mmd_dec = axis_param(dec_axis)

    freq_axis = map_in.get_axis('freq')/1.e6
    mmd_freq = axis_param(freq_axis)

    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    beam_FWHM = interp1d(freq_data, beam_data)

    beam_minfreq = beam_FWHM(mmd_freq[0])
    beam_maxfreq = beam_FWHM(mmd_freq[1])

    ra_fact = sp.cos(sp.pi * map_in.info['dec_centre'] / 180.0)
    dracosdec = mmd_ra[2]*ra_fact

    pixelfrac_ramin = dracosdec/beam_minfreq
    pixelfrac_decmin = mmd_dec[2]/beam_minfreq
    pixelfrac_ramax = dracosdec/beam_maxfreq
    pixelfrac_decmax = mmd_dec[2]/beam_maxfreq
    if not silent:
        infolist = ['ra_centre', 'ra_delta',
                    'dec_centre', 'dec_delta',
                    'freq_centre', 'freq_delta']
        for key in infolist:
            print "%s: %s" % (key, map_in.info[key])

        shp = map_in.shape
        print shp
        # 2^30 or 10^9...
        covsize = shp[0] * shp[1] ** 2. * shp[2] ** 2. * 4. / 2.**30.
        print "covariance size in GB: %g" % covsize
        print "RA min=%g deg, max=%g deg, delta=%g deg, extent=%g deg, N=%d" % mmd_ra
        print "Dec min=%g deg, max=%g deg, delta=%g deg, extent=%g deg, N=%d" % mmd_dec
        print "Freq min=%g MHz, max=%g MHz, delta=%g MHz, extent=%g MHz, N=%d" % mmd_freq

        print "Note that these fractions are in delta RA * cos(Dec)"
        print "beam FWHM at %g MHz=%g deg; pixel frac RA, Dec=(%g, %g)" % \
              (mmd_freq[0], beam_minfreq, pixelfrac_ramin, pixelfrac_decmin)

        print "beam FWHM at %g MHz=%g deg; pixel frac RA, Dec=(%g, %g)" % \
              (mmd_freq[1], beam_maxfreq, pixelfrac_ramax, pixelfrac_decmax)
        print "-"*80

    # return the sampling at the finest resolution
    return (pixelfrac_ramax, pixelfrac_decmax)


def print_map_summary():
    #print "former 15hr field dimensions"
    #find_map_dimensions('/mnt/raid-project/gmrt/kiyo/wiggleZ/maps/sec_A_15hr_41-73_clean_map_I.npy')
    print "15hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/tcv/maps/sec_A_15hr_41-90_clean_map_I.npy')
    print "22hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/tcv/maps/sec_A_22hr_41-90_clean_map_I.npy')
    print "1hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/tcv/maps/sec_A_1hr_none_clean_map_I.npy')


def new_map_templates(target_sample=0.25, multiplier=16, search_start=16):

    #print "proposed 15hr field dimensions"
    #template_15hr = find_map_region(220.3, 215.5, 0.7, 3.3,
    #                target_sample=target_sample, multiplier=multiplier,
    #                search_start=search_start)
    #find_map_dimensions(template_15hr)

    #print "proposed 22hr field dimensions"
    #template_22hr = find_map_region(327.9, 323., -1.5, 1.5,
    #                target_sample=target_sample, multiplier=multiplier,
    #                search_start=search_start)
    #find_map_dimensions(template_22hr)

    print "proposed 1hr field dimensions"
    #template_1hr = find_map_region(18., 8., -1.6, 4.4,
    #                target_sample=target_sample, multiplier=multiplier,
    #                search_start=search_start)
    #template_1hr = find_map_region(18., 8., -0.7, 4.4,
    #                target_sample=target_sample, multiplier=multiplier,
    #                search_start=search_start)
    #find_map_dimensions(template_1hr)
    # 100 GB
    template_1hr = find_map_region(17., 9., -0.5, 4.3,
                    target_sample=target_sample, multiplier=multiplier,
                    search_start=search_start)
    find_map_dimensions(template_1hr)

    #template_1hr = find_map_region(17.3, 8.6, -0.6, 4.4,
    #                target_sample=target_sample, multiplier=multiplier,
    #                search_start=search_start)
    #find_map_dimensions(template_1hr)


if __name__ == '__main__':
    new_map_templates(target_sample=0.25, multiplier=1, search_start=2)
    #print "previous map dimensions"
    #print "="*80
    #print_map_summary()
    #new_map_templates(target_sample=0.25, multiplier=16, search_start=16)
