import os
from subprocess import Popen
import ConfigParser
import numpy as np
from scipy.interpolate import interp1d
import shlex


# Open params_orig.ini to read comments
# Since this code remove them
def modify_params(cosmo_params):
    """
    Open params.ini file and set parameters to given in
    a dictionary cosmo_params
    """
    home = os.environ['HOME']
    cambway = '/code/CAMB/'
    config = ConfigParser.ConfigParser()
    print home + cambway + 'params.ini'
    config.read(home + cambway + 'params.ini')
    part = 'Cosmology'
    for key in cosmo_params.keys():
        config.set(part, key, cosmo_params[key])
    with open(home + cambway + 'params.ini', 'w') as configfile:
        config.write(configfile)

def camb_call():
    """
    Call CAMB code to output .dat files, with a params.ini file.
    """
    home = os.environ['HOME']
    cambway = '/code/CAMB/'
    CAMB = home + cambway + 'camb'
    PARAMS = home + cambway + 'params.ini'
    process = Popen([CAMB, PARAMS])
    process.communicate()

def interpolate(filename, column1, column2):
    """
    Calculate the interpolated log-log function. To call it after
    use 10**f(10**k) where k represents wave vector
    """
    data = np.genfromtxt(filename)
    dat_1 = data[:,column1]
    dat_2 = data[:,column2]
    f = interp1d(np.log10(dat_1), np.log10(dat_2))
    return f

def camb_pipe(params):
    """
    Modify params at .ini file. Call camb code with updated .ini file.
    Return interpolated powerspectrum and transfer function.
    """
    modify_params(params)
    camb_call()
    P = interpolate('test_matterpower.dat', 0, 1)
    T = interpolate('test_transfer_out.dat', 0, 6)
    return P, T


def set_initial_params():
    """
    Set parameters at the .ini file to the initial ones.
    """
    params_dict = {'hubble': 70, 'transfer_redshift(1)': 0, 'omk': 0,
                   'omega_baryon': 0.046, 'omega_cdm': 0.224,
                   'omega_lambda': 0.73, 'ombh2': 0.022, 'omnuh2': 0.,
                   'massive_neutrinos': 0., 'scalar_amp(1)': 2.3e-9}
    modify_params(params_dict)
