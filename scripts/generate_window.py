from utils import data_paths
from core import algebra
from utils import fftutil
from map import physical_gridding

def generate_windows(window="blackman"):
    datapath_db = data_paths.DataPath()
    # first generate a window for the full physical volume
    filename = datapath_db.fetch('simideal_15hr_physical', intend_read=True,
                                 pick='1')
    print filename
    pcube = algebra.make_vect(algebra.load(filename))
    pwindow = algebra.make_vect(fftutil.window_nd(pcube.shape, name=window),
                                axis_names=('freq', 'ra', 'dec'))
    pwindow.copy_axis_info(pcube)
    print pwindow.shape
    algebra.save("physical_window.npy", pwindow)

    # now generate one for the observed region and project onto the physical
    # volume.
    filename = datapath_db.fetch('simideal_15hr_beam', intend_read=True,
                                 pick='1')
    print filename
    ocube = algebra.make_vect(algebra.load(filename))
    owindow = algebra.make_vect(fftutil.window_nd(ocube.shape, name=window),
                                axis_names=('freq', 'ra', 'dec'))
    owindow.copy_axis_info(ocube)
    print owindow.shape
    print owindow.axes
    algebra.save("observed_window.npy", owindow)
    pwindow = physical_gridding.physical_grid(owindow, refinement=2)
    print pwindow.shape
    algebra.save("observed_window_physreg.npy", pwindow)

if __name__ == '__main__':
    generate_windows()

