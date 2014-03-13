import scipy as sp
from scipy.interpolate import interp1d
import core.algebra as algebra
from utils.cosmology import Cosmology
from utils import units
from utils import misc
from core import constants as cc
# TODO: clean up the variable names, shorten code

def axis_param(axis):
    r"""return min, max and delta assuming unif. spacing"""
    (min_axis, max_axis) = (min(axis), max(axis))
    extent = max_axis-min_axis
    return (min_axis, max_axis, abs(axis[1]-axis[0]), extent, len(axis))


def find_map_region(min_ra, max_ra, min_dec, max_dec, target_sample=0.25,
                    multiplier=16, search_start=0,
                    max_freq=700.391, min_freq=899.609, n_freq=256,
                    exact_freq=True):
    r"""target 0.25 pixels/FWHM"""
    # do multiples of 16 so that 256 freq * 16N is a multiple of 4096

    ra_sampling = target_sample*2.  # initial value that will not fail
    n_ra = search_start
    while ra_sampling > target_sample:
        n_ra += multiplier
        ra_sampling = find_map_dimensions(define_map_region(min_freq, max_freq, min_ra, max_ra,
                                                     min_dec, max_dec, n_freq,
                                                     n_ra, 32, exact_freq=True), silent=True)
        ra_sampling = ra_sampling[0]

    dec_sampling = target_sample*2.
    n_dec = search_start
    while dec_sampling > target_sample:
        n_dec += multiplier
        dec_sampling = find_map_dimensions(define_map_region(min_freq, max_freq, min_ra, max_ra,
                                           min_dec, max_dec, n_freq,
                                           n_ra, n_dec, exact_freq=True), silent=True)
        dec_sampling = dec_sampling[1]

    ra_sample_ratio = target_sample/ra_sampling
    dec_sample_ratio = target_sample/dec_sampling

    print "n_ra=%d, samp=%g, ratio=%g" % (n_ra, ra_sampling,
                                          target_sample/ra_sampling)

    print "n_dec=%d, samp=%g, ratio=%g" % (n_dec, dec_sampling,
                                          target_sample/dec_sampling)

    print "original ra=(%g,%g), dec=(%g,%g)" % (min_ra, max_ra,
                                                min_dec, max_dec)

    template_map = define_map_region(min_freq, max_freq, min_ra, max_ra,
                                     min_dec, max_dec, n_freq, n_ra, n_dec)
    # now expand the map a bit so that pixels are exactly 0.25 deg
    info = template_map.info

    blank = sp.zeros((n_freq, n_ra, n_dec))
    info['ra_delta'] *= ra_sample_ratio
    info['dec_delta'] *= dec_sample_ratio
    map_prod = algebra.make_vect(blank, axis_names=('freq', 'ra', 'dec'))
    map_prod.info = info

    return map_prod


def define_map_region(min_freq, max_freq, min_ra, max_ra, min_dec, max_dec,
                      n_freq, n_ra, n_dec, exact_freq=False):
    """target: n_freq=256, 'freq_centre': 799609375.0, 'freq_delta': -781250.0
    Freq min=700.391 MHz, max=899.609 MHz, delta=0.78125 MHz, extent=199.219 MHz, N=256
    new: Freq min=700.391 MHz, max=899.609 MHz, delta=0.781247 MHz, extent=199.218 MHz
    new: freq_centre: 800390623.529, freq_delta: 781247.058824
    """
    min_freq *= 1.e6
    max_freq *= 1.e6

    #axis = (delta*(sp.arange(len) - len//2) + centre)
    freq_delta = (max_freq - min_freq)/float(n_freq-1)
    freq_centre = min_freq + freq_delta*float((n_freq)//2)

    ra_delta = (max_ra - min_ra)/float(n_ra-1)
    ra_centre = min_ra + ra_delta*float((n_ra)//2)

    dec_delta = (max_dec - min_dec)/float(n_dec-1)
    dec_centre = min_dec + dec_delta*float((n_dec)//2)

    if exact_freq:
        freq_delta = -781250.0
        freq_centre = 799609375.0
        n_freq = 256

    info = {'ra_delta': ra_delta,
            'dec_delta': dec_delta,
            'dec_centre': dec_centre,
            'axes': ('freq', 'ra', 'dec'),
            'ra_centre': ra_centre,
            'freq_centre': freq_centre,
            'freq_delta': freq_delta,
            'type': 'vect'}

    blank = sp.zeros((n_freq, n_ra, n_dec))
    map_prod = algebra.make_vect(blank, axis_names=('freq', 'ra', 'dec'))
    map_prod.info = info
    return map_prod


def find_map_dimensions(map_in, silent=False):
    r"""print various aspects of the maps"""
    if isinstance(map_in, str):
        map_in = algebra.make_vect(algebra.load(map_in))
    cosmology = Cosmology()

    ra_axis = map_in.get_axis('ra')
    #print ra_axis
    mmd_ra = axis_param(ra_axis)

    dec_axis = map_in.get_axis('dec')
    #print dec_axis
    mmd_dec = axis_param(dec_axis)

    freq_axis = map_in.get_axis('freq')/1.e6
    mmd_freq = axis_param(freq_axis)

    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    beam_FWHM = interp1d(freq_data, beam_data)

    if (mmd_freq[0] < freq_data.min()):
        beam_minfreq = beam_FWHM(freq_data.min())
        #print "warning: no beam data for this freq, using %s" % beam_minfreq
    else:
        beam_minfreq = beam_FWHM(mmd_freq[0])

    if (mmd_freq[1] > freq_data.max()):
        beam_maxfreq = beam_FWHM(freq_data.max())
        #print "warning: no beam data for this freq, using %s" % beam_maxfreq
    else:
        beam_maxfreq = beam_FWHM(mmd_freq[1])

    ra_fact = sp.cos(sp.pi * map_in.info['dec_centre'] / 180.0)
    dracosdec = mmd_ra[2]*ra_fact

    pixelfrac_ramin = dracosdec/beam_minfreq
    pixelfrac_decmin = mmd_dec[2]/beam_minfreq
    pixelfrac_ramax = dracosdec/beam_maxfreq
    pixelfrac_decmax = mmd_dec[2]/beam_maxfreq

    # cosmological sizes
    dfreq = mmd_freq[2]
    z1 = cc.freq_21cm_MHz/mmd_freq[0] - 1.  # z at min freq
    z2 = cc.freq_21cm_MHz/mmd_freq[1] - 1.  # z at max freq
    dz1 = cc.freq_21cm_MHz/(mmd_freq[0]+dfreq) - 1.
    dz2 = cc.freq_21cm_MHz/(mmd_freq[1]-dfreq) - 1.
    d1 = cosmology.proper_distance(z1)
    d2 = cosmology.proper_distance(z2)
    c1 = cosmology.comoving_distance(z1)
    c2 = cosmology.comoving_distance(z2)
    dc1 = cosmology.comoving_distance(dz1)
    dc2 = cosmology.comoving_distance(dz2)

    ra_extent_z1 = mmd_ra[3] * d1 * units.degree * ra_fact
    ra_pixel_extent_z1 = mmd_ra[2] * d1 * units.degree * ra_fact
    ra_extent_z2 =  mmd_ra[3] * d2 * units.degree * ra_fact
    ra_pixel_extent_z2 = mmd_ra[2] * d2 * units.degree * ra_fact
    dec_extent_z1 = mmd_dec[3] * d1 * units.degree
    dec_pixel_extent_z1 = mmd_dec[2] * d1 * units.degree
    dec_extent_z2 =  mmd_dec[3] * d2 * units.degree
    dec_pixel_extent_z2 = mmd_dec[2] * d2 * units.degree

    redshift_extent = c1 - c2
    redshift_pixel_extent_z1 = c1 - dc1
    redshift_pixel_extent_z2 = dc2 - c2

    if not silent:
        infolist = ['ra_centre', 'ra_delta',
                    'dec_centre', 'dec_delta',
                    'freq_centre', 'freq_delta']
        for key in infolist:
            print "%s: %s" % (key, map_in.info[key])

        sexagesimal = misc.radec_to_sexagesimal(map_in.info["ra_centre"],
                                                 map_in.info["dec_centre"])
        (ralong, declong, decsign) = sexagesimal
        print "RA: %dH:%dm:%10.5fs, dec %s %dd:%dm:%10.5fs" % \
              (ralong[0], ralong[1], ralong[2], decsign, declong[0], declong[1],
               declong[2])

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


        print "redshift range: z=%g-%g" % (z2, z1)
        print "proper distance range (h^-1 cMpc): d=%g-%g" % (d2, d1)
        print "comoving distance range (h^-1 cMpc): d=%g-%g, extent = %g" % \
              (c2, c1, redshift_extent)
        print "at %g MHz, width in RA=%g Dec=%g h^-1 cMpc" % \
              (mmd_freq[0], ra_extent_z1, dec_extent_z1)
        print "at %g MHz, pixel width in RA=%g Dec=%g h^-1 cMpc" % \
              (mmd_freq[0], ra_pixel_extent_z1, dec_pixel_extent_z1)
        print "the freq bin starting at %g MHz has width %g" % \
                (mmd_freq[0], redshift_pixel_extent_z1)

        print "at %g MHz, width in RA=%g Dec=%g h^-1 cMpc" % \
              (mmd_freq[1], ra_extent_z2, dec_extent_z2)
        print "at %g MHz, pixel width in RA=%g Dec=%g h^-1 cMpc" % \
              (mmd_freq[1], ra_pixel_extent_z2, dec_pixel_extent_z2)
        print "the freq bin starting at %g MHz has width %g" % \
                (mmd_freq[1], redshift_pixel_extent_z2)

        print "-"*80

    # return the sampling at the finest resolution
    return (pixelfrac_ramax, pixelfrac_decmax)


def print_map_summary():
    print "15hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/sec_A_15hr_41-90_clean_map_I.npy')
    #print "22hr field dimensions"
    #find_map_dimensions('/mnt/raid-project/gmrt/tcv/maps/sec_A_22hr_41-90_clean_map_I.npy')
    print "1hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal/secA_1hr_41-18_clean_map_I_800.npy')


def new_map_templates(target_sample=0.25, multiplier=16, search_start=16):

    #print "proposed 15hr field dimensions"
    template_15hr = find_map_region(220.3, 215.5, 0.7, 3.3,
                    target_sample=target_sample, multiplier=multiplier,
                    search_start=search_start)
    find_map_dimensions(template_15hr)

    #print "proposed 22hr field dimensions"
    template_22hr = find_map_region(327.9, 323., -1.5, 1.5,
                    target_sample=target_sample, multiplier=multiplier,
                    search_start=search_start)
    find_map_dimensions(template_22hr)

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


def complete_wigglez_regions(target_sample=0.25, multiplier=16,
                             search_start=16):
    """15hr:
    ([214.00001499999999, 222.99998500000001], [3.0000000000000001e-06,
    3.9999989999999999], [676383691.31904757, 946937167.84666669])

    22hr:
    ([322.00003099999998, 329.99996900000002], [-1.9999990000000001,
    1.9999979999999999], [676383691.31904757, 946937167.84666669])

    1hr
    ([9.0000020000000003, 15.999999000000001], [-1.9999979999999999, 1.999997],
    [676383691.31904757, 946937167.84666669])
    """

    #print "proposed 15hr field dimensions"
    template_15hr = find_map_region(223., 214., 0., 4,
                    target_sample=target_sample, multiplier=multiplier,
                    search_start=search_start, exact_freq=False,
                    max_freq=676383691.31904757/1.e6,
                    min_freq=946937167.84666669/1.e6, n_freq=512)
    find_map_dimensions(template_15hr)
    algebra.save("templates/wigglez_15hr_complete.npy", template_15hr)

    #print "proposed 22hr field dimensions"
    template_22hr = find_map_region(330., 322., -2., 2.,
                    target_sample=target_sample, multiplier=multiplier,
                    search_start=search_start, exact_freq=False,
                    max_freq=676383691.31904757/1.e6,
                    min_freq=946937167.84666669/1.e6, n_freq=512)
    find_map_dimensions(template_22hr)
    algebra.save("templates/wigglez_22hr_complete.npy", template_22hr)

    #print "proposed 1hr field dimensions"
    template_1hr = find_map_region(16., 9., -2., 2.,
                    target_sample=target_sample, multiplier=multiplier,
                    search_start=search_start, exact_freq=False,
                    max_freq=676383691.31904757/1.e6,
                    min_freq=946937167.84666669/1.e6, n_freq=512)
    find_map_dimensions(template_1hr)
    algebra.save("templates/wigglez_1hr_complete.npy", template_1hr)


if __name__ == '__main__':
    #new_map_templates(target_sample=0.25, multiplier=1, search_start=2)
    #complete_wigglez_regions(target_sample=0.25, multiplier=1, search_start=2)
    #print "previous map dimensions"
    #print "="*80
    print_map_summary()
    #new_map_templates(target_sample=0.25, multiplier=16, search_start=16)
