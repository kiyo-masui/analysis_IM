import core.algebra as al
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
from itertools import combinations
import pylab
from matplotlib.colors import LogNorm
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages

def plot_parkes_fwhm(fwhm_file, output, beam, compare_model = True):
    fwhms = np.load(fwhm_file)
    fwhm = fwhms[beam,:]
    freqs = np.linspace(1315.5-32,1315.5+32, fwhm.shape[0])
    pylab.plot(freqs, fwhm)
    if compare_model:
        model = (14.4/60.)*(1420./freqs)
    pylab.plot(freqs,model)
    pylab.ylim([0.2,0.4])
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('FWHM (degrees)')
    pylab.savefig(output + '.png')
    pylab.clf() 

def plot_parkes_fwhm_pdf(fwhm_file, output):
    fwhms = np.load(fwhm_file)
    freqs = np.linspace(1315.5-32,1315.5+32, fwhms.shape[-1])
    model = (14.4/60.)*(1420./freqs)
    pdf = PdfPages(output + '.pdf')
    #font = {'family' : 'normal',
    #    'weight' : 'normal',
    #    'size'   : 2}
    #matplotlib.rc('font', **font)
    matplotlib.rc('xtick', labelsize=12) 
    matplotlib.rc('ytick', labelsize=12) 
    for b in xrange(13):
        pylab.clf()
        pylab.plot(freqs,fwhms[b,1,:], label = 'YY')
        pylab.plot(freqs,fwhms[b,0,:], label = 'XX')
        pylab.plot(freqs,model, label = 'Model')
        pylab.ylim([0.2,0.4])
        pylab.xlabel('Frequency (MHz)', fontsize = 12)
        pylab.ylabel('FWHM (degrees)', fontsize = 12)
        pylab.legend(loc='upper right',prop={'size':8})
        pylab.title(r'Beam ' + str(b+1), fontsize = 12)
        pdf.savefig()
        pylab.close()
    pdf.close()
