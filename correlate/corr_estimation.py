import numpy as np
from core import algebra


def corr_est(map1, map2, noise1, noise2, freq1, freq2,
             lags=(), speedup=False, verbose=False):
    r"""Calculate the cross correlation function of the maps.

    The cross correlation function is a function of f1, f2 and angular lag.
    The angular lag bins are passed, all pairs of frequencies are
    calculated.

    Parameters
    ----------
    lags: array like
        Angular lags bins (upper side bin edges).
    speedup: boolean
        Speeds up the correlation. This works fine, yes? Should be the
        normal way if so.

    Returns
    -------
    corr: array
        The correlation between 2 maps.
    counts: array
        The weighting of the correlation based on the maps' weights.

    """
    map1_ra = map1.get_axis('ra')
    map2_ra = map2.get_axis('ra')
    map1_dec = map1.get_axis('dec')
    map2_dec = map2.get_axis('dec')

    input_map1 = map1[freq1, :, :]
    input_map2 = map2[freq2, :, :]
    input_noise1 = noise1[freq1, :, :]
    input_noise2 = noise2[freq2, :, :]

    # Noise weight
    input_map1 *= input_noise1
    input_map2 *= input_noise2

    nlags = len(lags)
    nfreq = len(freq1)
    corr = sp.zeros((nfreq, nfreq, nlags), dtype=float)
    counts = sp.zeros(corr.shape, dtype=float)
    # Noting that if DEC != 0, then a degree of RA is less than a degree.
    ra_fact = sp.cos(sp.pi * map1.info['dec_centre'] / 180.0)

    # Calculate the pairwise lags.
    dra = (map1_ra[:, None] - map2_ra[None, :]) * ra_fact
    ddec = map1_dec[:, None] - map2_dec[None, :]
    lag = dra[:, None, :, None] ** 2 + ddec[None, :, None, :] ** 2
    lag = sp.sqrt(lag)
    # Bin this up.
    lag_inds = sp.digitize(lag.flatten(), lags)

    if speedup:
        print "Starting Correlation (sparse version)"

        (nr1, nd1) = (len(map1_ra), len(map1_dec))
        (nr2, nd2) = (len(map2_ra), len(map2_dec))
        (r1ind, d1ind) = (sp.arange(nr1), sp.arange(nd1))
        (r2ind, d2ind) = (sp.arange(nr2), sp.arange(nd2))
        ra1_pairind = r1ind.repeat(nr2 * nd1 * nd2)
        ra2_pairind = sp.tile(r2ind.repeat(nd2), (1, nr1 * nd1)).flatten()
        dec1_pairind = sp.tile(d1ind.repeat(nr2 * nd2), (1, nr1)).flatten()
        dec2_pairind = sp.tile(d2ind, (1, nr1 * nr2 * nd1)).flatten()

        # precalculate the pair indices for a given lag
        # could also imagine calculating the map slices here
        posmaskdict = {}
        for klag in range(nlags):
            mask = (lag_inds == klag)
            posmaskdict[repr(klag)] = (ra1_pairind[mask],
                                       ra2_pairind[mask],
                                       dec1_pairind[mask],
                                       dec2_pairind[mask])

        for if1 in range(len(freq1)):
            for jf2 in range(len(freq2)):
                start = time.time()

                data1 = input_map1[if1, :, :]
                data2 = input_map2[jf2, :, :]
                weights1 = input_noise1[if1, :, :]
                weights2 = input_noise2[jf2, :, :]

                for klag in range(nlags):
                    (r1m, r2m, d1m, d2m) = posmaskdict[repr(klag)]
                    dprod = data1[r1m, d1m] * data2[r2m, d2m]
                    wprod = weights1[r1m, d1m] * weights2[r2m, d2m]
                    corr[if1, jf2, klag] += sp.sum(dprod)
                    counts[if1, jf2, klag] += sp.sum(wprod)

                if verbose:
                    print if1, jf2, (time.time() - start)
                    print counts[if1, jf2, :]
    else:
        print "Starting Correlation (full version)"
        for if1 in range(len(freq1)):
            for jf2 in range(len(freq2)):
                start = time.time()
                # Calculate the pairwise products.
                data1 = input_map1[if1, :, :]
                data2 = input_map2[jf2, :, :]
                weights1 = input_noise1[if1, :, :]
                weights2 = input_noise2[jf2, :, :]
                dprod = data1[..., None, None] * data2[None, None, ...]
                wprod = weights1[..., None, None] * \
                        weights2[None, None, ...]
                for klag in range(nlags):
                    mask = (lag_inds == klag)
                    corr[if1, jf2, klag] += sp.sum(dprod.flatten()[mask])
                    counts[if1, jf2, klag] += sp.sum(wprod.flatten()[mask])

                if verbose:
                    print if1, jf2, (time.time() - start)
                    print counts[if1, jf2, :]

    mask = (counts < 1e-20)
    counts[mask] = 1
    corr /= counts
    corr[mask] = 0
    counts[mask] = 0

    return corr, counts

