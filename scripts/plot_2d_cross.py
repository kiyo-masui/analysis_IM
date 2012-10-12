from quadratic_products import power_spectrum

root = "/mnt/raid-project/gmrt/eswitzer/GBT/pwrspec/GBT_15hr_map_oldcal_x_WiggleZ_15hr_blackman_order1/"
datafile = root + "GBT_15hr_map_oldcal_x_WiggleZ_15hr_data.shelve"
mockfile = root + "GBT_15hr_map_oldcal_x_WiggleZ_15hr_mock.shelve"

dataobj = power_spectrum.PowerSpectrum(datafile)
mockobj = power_spectrum.PowerSpectrum(mockfile)

dataobj.summarize_2d_agg_pwrspec("20modes", "xspec_data.dat")
mockobj.summarize_2d_agg_pwrspec("20modes", "xspec_mock.dat")

