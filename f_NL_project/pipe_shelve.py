import subprocess
import matplotlib.pyplot as plt
import shelve_outplot as so
import numpy as np

# Initialize settings
shelve_calc = False
pwrspec = False
r_val = True
manager_path = '/cita/h/home-2/mufma/code/analysis_IM/pipeline/manager.py'
ini_path = 'input/mm/autopower/kappa_dm_sim.ini'

# Produce shelve file from the algebra.npy file
if shelve_calc == True:
    subprocess.Popen(['python', manager_path, ini_path])

# Set up maps for reading
path = "/mnt/raid-project/gmrt/mufma/JD_DM_Sim/gal_21_pwrspec_512/"
cross_map = path + "gal_21_cross_0.500red.shelve"
gal_map = path + "gal_auto0.500red.shelve"
T_b21_map = path + "auto_21_0.500red.shelve"

# Plot the r value
if pwrspec == True:
    from matplotlib import rc
    rc('text', usetex=True)

    p1 = so.plot_pwrspec(cross_map)
    p2 = so.plot_pwrspec(gal_map)
    p4 = so.plot_pwrspec(T_b21_map)
    p3 = so.plot_camb()
    so.plot_legend([p1,p2,p4,p3],[r'$P_{gal\times21}$',
                   r'$P_{gal\times gal}$',
                   r'$P_{21\times 21}$', r'$P_{CAMB}$'])
    plt.show()
    """
    p1 = so.plot_bias(gal_map)
    p2 = so.plot_bias(T_b21_map)
    so.plot_legend([p1,p2],[r'$b_{gal}$', r'$b_{21}*T^{2}_{b}$'])
    plt.show()
    """
if r_val == True:
    so.plot_r_val(cross_map, gal_map, T_b21_map)
    plt.show()
