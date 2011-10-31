import scipy as sp
from scipy.interpolate import interp1d
import core.algebra as algebra

def axis_param(axis):
    r"""return min, max and delta assuming unif. spacing"""
    (min_axis, max_axis) = (min(axis), max(axis))
    extent = max_axis-min_axis
    return (min_axis, max_axis, abs(axis[1]-axis[0]), extent, len(axis))


def test_15hr_map():
    blank = sp.zeros((256, 64, 32))

    info = {'ra_delta': -0.075045715822411624,
            'dec_delta': 0.075,
            'dec_centre': 2.0,
            'axes': ('freq', 'ra', 'dec'),
            'ra_centre': 217.9,
            'freq_centre': 799609375.0,
            'freq_delta': -781250.0,
            'type': 'vect'}

    map_prod = algebra.make_vect(blank, axis_names=('freq', 'ra', 'dec'))
    map_prod.info = info
    return map_prod


def find_map_dimensions(map_in, print_infofield=True):
    r"""print various aspects of the maps"""
    if isinstance(map_in, str):
        map_in = algebra.make_vect(algebra.load(map_in))

    print map_in.shape

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

    if print_infofield:
        for key in map_in.info:
            print "%s: %s" % (key, map_in.info[key])

    print "RA min=%g deg, max=%g deg, delta=%g deg, extent=%g deg, N=%d" % mmd_ra
    print "Dec min=%g deg, max=%g deg, delta=%g deg, extent=%g deg, N=%d" % mmd_dec
    print "Freq min=%g MHz, max=%g MHz, delta=%g MHz, extent=%g MHz, N=%d" % mmd_freq

    beam_minfreq = beam_FWHM(mmd_freq[0])
    beam_maxfreq = beam_FWHM(mmd_freq[1])
    pixelfrac_ramin = mmd_ra[2]/beam_minfreq
    pixelfrac_decmin = mmd_dec[2]/beam_minfreq
    pixelfrac_ramax = mmd_ra[2]/beam_maxfreq
    pixelfrac_decmax = mmd_dec[2]/beam_maxfreq
    print "beam FWHM at %g MHz=%g deg; pixel frac RA, Dec=(%g, %g)" % \
          (mmd_freq[0], beam_minfreq, pixelfrac_ramin, pixelfrac_decmin)

    print "beam FWHM at %g MHz=%g deg; pixel frac RA, Dec=(%g, %g)" % \
          (mmd_freq[1], beam_maxfreq, pixelfrac_ramax, pixelfrac_decmax)
    print "-"*80

if __name__ == '__main__':
    print "former 15hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/kiyo/wiggleZ/maps/sec_A_15hr_41-73_clean_map_I.npy')
    print "15hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/tcv/maps/sec_A_15hr_41-90_clean_map_I.npy')
    print "22hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/tcv/maps/sec_A_22hr_41-90_clean_map_I.npy')
    print "1hr field dimensions"
    find_map_dimensions('/mnt/raid-project/gmrt/tcv/maps/sec_A_1hr_none_clean_map_I.npy')

    print "proposed 15hr field dimensions"
    find_map_dimensions(test_15hr_map())
