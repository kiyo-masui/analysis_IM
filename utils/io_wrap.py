"""wrap IO through shelves and pickle files"""
import os
import cPickle

def save_pickle(pickle_data, filename):
    r"""cPickle the `pickle_data` to file with name `filename`.
    """

    pickle_out_root = os.path.dirname(filename)
    if not os.path.isdir(pickle_out_root):
        os.mkdir(pickle_out_root)

    pickle_handle = open(filename, 'w')
    cPickle.dump(pickle_data, pickle_handle)
    pickle_handle.close()


def load_pickle(filename):
    r"""Return `pickle_data` saved in file with name `filename`.
    """
    pickle_handle = open(filename, 'r')
    pickle_data = cPickle.load(pickle_handle)
    pickle_handle.close()
    return pickle_data



