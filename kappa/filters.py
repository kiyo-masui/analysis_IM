import numpy as np


def spherical_separation_filter(map, k_axes, hi_k, lo_k=None, high=True, low=True, isReal=True):
    r"""Function performing a cutoff in k space and returning two maps :
    one with all k smaller than the cutoff, the other with all k higher than the cutoff"""
    if low & (lo_k is None): lo_k = hi_k
    l, m, n = map.shape
    if isReal:
        f_map = np.fft.ifftshift(map)
        f_map = np.fft.fftn(f_map)
    else:
        f_map = map

    kz, kx, ky = k_axes
    rad_3d = (kz*kz)[:, None, None] * np.ones((1, m, n))
    rad_3d += (kx*kx)[None, :, None] * np.ones((l, 1, n))
    rad_3d += (ky*ky)[None, None, :] * np.ones((l, m, 1))

    res = []

    if high:
        hi_map = rad_3d > (hi_k*hi_k)
        hi_map = np.fft.ifftshift(hi_map)
        hi_map = f_map*hi_map
        hi_map = np.fft.ifftn(hi_map)
        hi_map = np.fft.fftshift(hi_map.real)
        res.append(hi_map)

    if low:
        lo_map = rad_3d <= (lo_k*lo_k)
        lo_map = np.fft.ifftshift(lo_map)
        lo_map = f_map*lo_map
        lo_map = np.fft.ifftn(lo_map)
        lo_map = np.fft.fftshift(lo_map.real)
        res.append(lo_map)

    if len(res) == 1: return res[0]
    else: return res


def scale_separation(map_s, n=20):
    r"""Find the correlation coefficient between large scales in the input map
    and kappa reconstructed from small scales, f different values of the cutoff""" 
    r = np.zeros(n)
    r_ref = np.zeros(n)
    w = map_s.gaussian_window(8.)
    for j in range(n):
        k = j*0.48/30.
        print "step", j, ": k=", k
        map_s.modify_filters(k, k, reconstruct=True, kz_filt=True)
        map_s.kappa_hi = fft_convolve(map_s.kappa_hi, w)
        if map_s.mask is not None:
            map_s.kappa_hi = map_s.mask*map_s.kappa_hi
        m1m1 = np.sum(map_s.kappa_hi * map_s.kappa_hi)
        m2m2 = np.sum(map_s.lo_map * map_s.lo_map)
        m1m2 = np.sum(map_s.kappa_hi * map_s.lo_map)
        r[j] = m1m2/np.sqrt(m1m1*m2m2)
        print r[j]
        m1m1b = np.sum(map_s.hi_map * map_s.hi_map)
        m1m2b = np.sum(map_s.hi_map * map_s.lo_map)
        r_ref[j] = m1m2b/np.sqrt(m1m1b*m2m2)
    return r, r_ref


def kappa_filter(kappa, hi_k, axes):
    return spherical_separation_filter(kappa, axes, hi_k, hi_k, high=False)


def kappa_autocorr(map_s, dk, mask=None, kappa=None):
    if kappa is None: kappa = map_s.kappa
    kz = map_s.get_k_axis('z')
    kx = map_s.get_k_axis('x')
    ky = map_s.get_k_axis('y')
    r = np.zeros(20)
    for j in range(20):
        k = j*dk
        print "step", j, ": k=", k
        map_s.modify_filters(k, k)
        kap_tmp = kappa_filter(kappa, k, [kz, kx, ky])
        if mask is not None:
            kap_tmp = mask*kap_tmp
        m1m1 = np.sum(kap_tmp * kap_tmp)
        m2m2 = np.sum(map_s.lo_map * map_s.lo_map)
        m1m2 = np.sum(kap_tmp * map_s.lo_map)
        r[j] = m1m2/np.sqrt(m1m1*m2m2)
        print r[j]
    return r
