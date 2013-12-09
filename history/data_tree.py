class DataTree(object):
    r"""intermediate between shelve/dict/h5py
    Due to what appears to be a bug in hd5py/hdf5, functional recursion does
    not seem to work (it is implemented below)
    """

    def __init__(self, filename):
        self.filename = filename
        self.dataset = {}
        self.root = h5py.File(filename, "w")

    def close(self):
        self.root.flush()
        self.root.close()
        time.sleep(1)
        print "done writing hd5file %s" % self.filename

    def absorb_data_dict(self, data_dict):
        self.traverse_data_dict(data_dict)

    def traverse_data_dict(self, data_dict, path=""):
        for data_key in data_dict:
            if isinstance(data_dict[data_key], np.ndarray):
                current_path = path + "/" + data_key + ".npy"
                self.dataset[current_path] = \
                        self.root.create_dataset(current_path,
                                                 data=data_dict[data_key])
                print "add numpy data:", current_path, data_dict[data_key].shape
                self.root.flush()
            elif not isinstance(data_dict[data_key], dict):
                pass
            else:
                current_path = path + "/oo" + data_key
                print "adding path ", current_path
                self.traverse_data_dict(data_dict[data_key], current_path)


class DataTree_old(object):
    r"""intermediate between shelve/dict/h5py
    Due to what appears to be a bug in hd5py/hdf5, functional recursion does
    not seem to work (it is implemented below)
    """

    def __init__(self, filename):
        self.filename = filename
        self.branches = {}
        self.branches["root"] = h5py.File(filename, "w")

    def close(self):
        print self.filename, " size ", self.branches["root"].userblock_size
        self.branches["root"].flush()
        self.branches["root"].close()
        print "done writing hd5file %s" % self.filename

    def absorb_data_dict(self, data_dict):
        self.traverse_data_dict(data_dict)

    def traverse_data_dict(self, data_dict, path="root"):
        for data_key in data_dict:
            if isinstance(data_dict[data_key], np.ndarray):
                current_path = path + "/" + data_key + ".npy"
                self.branches[current_path] = \
                    self.branches[path].create_dataset(data_key,
                                                data=data_dict[data_key])
                #self.branches[path] = data_dict[data_key]
                print "add numpy data:", current_path, data_dict[data_key].shape
                self.branches["root"].flush()
            elif not isinstance(data_dict[data_key], dict):
                pass
            else:
                current_path = path + "/" + data_key
                self.branches[current_path] = \
                        self.branches[path].create_group(data_key)

                print "make hd5 dir: ", self.branches[current_path].name

                self.traverse_data_dict(data_dict[data_key], current_path)

def class_convert_numpytree_hdf5(data_dict, filename):
    print "writing hd5file %s from dict tree" % filename
    datatree = DataTree(filename)
    datatree.absorb_data_dict(data_dict)
    datatree.close()

