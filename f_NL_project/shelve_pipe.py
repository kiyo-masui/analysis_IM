import subprocess
import shelve_outplot as so

# Initialize settings
shelve_calc = False
plotter = True
manager_path = '/cita/h/home-2/mufma/code/analysis_IM/pipeline/manager.py'
ini_path = 'input/mm/autopower/kappa_dm_sim.ini'

# Produce shelve file from the algebra.npy file
if shelve_calc == True:
    subprocess.Popen(['python', manager_path, ini_path])

# Set up maps for reading
cross_map = '/mnt/raid-project/gmrt/mufma/goblot_pwrspec/cross_pwr_kap_del.shelve'
kappa_map = '/mnt/raid-project/gmrt/mufma/goblot_pwrspec/pwr_delta.shelve'
delta_map = '/mnt/raid-project/gmrt/mufma/goblot_pwrspec/pwr_kappa.shelve'

# Plot the r value
if plotter == True:
    so.plot_r_val(cross_map, delta_map, kappa_map)
