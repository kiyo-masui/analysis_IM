import abundance_matching as am
import galaxy_calc as gc

# Initialize settings
HI = True
galaxy = True
redshift = [0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600,
            1.700, 1.800, 1.900, 2.000, 2.100, 2.200, 2.300, 2.400,
            2.500, 2.700, 3.000, 3.500, 4.000, 4.500, 5.000, 6.000,
            7.000, 8.000, 9.000, 10.000]
index = [0,1,2,3,4,5,6,7]
path = "/mnt/raid-project/gmrt/mufma/JD_h5py_catalogs/"
print "Redshifts are:"
print redshift


# Calculate number of galaxies
if galaxy == True:
    print "Processing galaxies"
    for red in redshift:
        print "starting processing redshift", red
        for ind in index:
            gc.galaxy_calc(path + '%.3fhalo_catalog000%d.hdf5'%(red,ind))
            print "finished processing index", ind

# Separation between tasks
print 20*"#", 7*" ", 20*"#", 7*" ", 20*"#"

# Calculate 21 temperature brightness
if HI == True:
    print "Processing 21 T_b"
    f = am.interpol_f()
    for red in redshift:
        print "starting processing redshift", red
        for ind in index:
            am.M_halo_to_M_HI(path + '%.3fhalo_catalog000%d.hdf5'%(red, ind),
                              f, 'r+')
            am.M_HI_to_21(path + '%.3fhalo_catalog000%d.hdf5'%(red, ind),
                          red,'r+')
            print "finished processing index", ind
