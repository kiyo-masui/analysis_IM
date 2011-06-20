import numpy as np

import radialprofile

def ps_azimuth(img, width = None, bwperp = 3, bwpar = 3, kmodes = True):

    fq = np.abs(np.fft.fftshift(np.fft.fftn(img)))**2

    hpar = fq[0].shape[0] / 2
    fqh = fq[hpar:]

    binsperp = (fq.shape[1] / 2) / bwpar
    binspar = fqh.shape[0] / bwpar

    ta = np.zeros((binspar, fq.shape[1], fq.shape[2]))
    tp = np.zeros((binspar, binsperp))

    for i in range(binspar):
        sl = np.mean(fqh[bwpar*i:bwpar*(i+1),:,:], axis=0)
        ta[i,:,:] = sl

        kperp, tp[i,:] = radialprofile.azimuthalAverage(sl, bw = bwperp)

    kpar = (np.arange(binspar)) * bwpar * 2*np.pi / width[0]
    kperp *= 2*np.pi / width[0]
    tp /= np.array(fq.shape).prod()**2
    if width != None:
        tp *= np.array(width).prod()

    #return ta, tp, kpar, kperp
    if kmodes:
        return tp, kpar, kperp
    else:
        return tp
    

    

    
