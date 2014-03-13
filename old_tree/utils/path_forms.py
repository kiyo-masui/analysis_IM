r"""the PathForms class constructs the database and provides several
functions to build groups of common files for convenience.
"""
import sys


class PathForms(object):
    r"""provide some macros for handling the path database
    desc: [required] one-sentence description of file/code used to produce it
    status: [required]; either 'does not exist', 'in development' or 'vetted'
    notes: notes about how the file was generated; markdown format (optional)
    Then, the path is specified depending on the data type
    type:
        * 'file' single file
        * 'path' directory for outputs
        * 'filelist' group of files with the same desc, etc.
        the file list is a dict where the keys specify the entry in that file
        series: e.g. '44' is the 44th of 100 simulations. The dictionary is
        un-ordered, so these objects also have a 'listindex' which gives the
        keys in a particular order.

    """
    def __init__(self, verbose=False):
        self.groups = {}
        self.pdb = {}
        self.grouplist = []
        self.verbose = verbose

        self.register_empty_group('Paths')  # always a group to hold paths

    def register_empty_group(self, group_name):
        r"""establish a new groups of files"""
        self.groups[group_name] = []
        self.grouplist.append(group_name)
        if self.verbose:
            print "registered group %s as []" % group_name

    def register_list_empty_groups(self, group_name_list):
        r"""establish a list of new groups of files"""
        for groupitem in group_name_list:
            self.register_empty_group(groupitem)

    def register_path(self, key, pathname, desc, notes=None, status=None):
        r"""register a path
        environment variables can be enclosed in [] in the pathname and will
        be interpreted later."""
        self.groups['Paths'].append(key)
        val = {'desc': desc}
        val['path'] = pathname

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        if self.verbose:
            print "registered path: " + repr(val)

        self._set_key(key, val)

    def register_envpath(self, env_var, tag, notes=None, status=None):
        r"""register a path from an environment variable"""
        desc = "env path to %s (%s) directory" % (env_var, tag)
        pathname = "[%s]/" % env_var
        self.register_path(env_var, pathname, desc, notes=notes, status=status)

    def fetch_path(self, key):
        r"""extract the directory associated with some path key in the db"""
        if key not in self.groups['Paths']:
            print "ERROR: %s not an identified path" % key
            sys.exit()

        pathentry = self.pdb[key]
        if self.verbose:
            print "fetched path: " + pathentry['path']

        return pathentry['path']

    def _set_key(self, key, value):
        r"""set a key in the database"""
        if key in self.pdb:
            print "ERROR: key %s is already set" % key
            sys.exit()

        self.pdb[key] = value

    def register_file(self, key, group_key, parent, filename, desc,
                      notes=None, status=None):
        r"""register a single file with the database"""
        val = {'desc': desc}
        self.groups[group_key].append(key)
        val['file'] = self.fetch_path(parent) + filename

        val['group_key'] = group_key
        val['parent'] = parent

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        if self.verbose:
            print "registered file: " + repr(val)

        self._set_key(key, val)

    def register_file_set(self, key, group_key, parent, prefix, indices, desc,
                 notes=None, status=None, suffix=".npy"):
        r""" macro to assign big file series
        `parent` is the DB tag to the parent directory
        `prefix` is the file prefix before the series index
        `indices` is the list of indices for the file list
        `desc` is the description
        """
        val = {'desc': desc}
        self.groups[group_key].append(key)
        fileroot = self.fetch_path(parent)

        val['group_key'] = group_key
        val['parent'] = parent

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        minmax = repr(min(indices)) + ", " + repr(max(indices))
        val['notes'] = '`%s%s[%s]%s`' % (fileroot, prefix, minmax, suffix)
        val['status'] = status

        filelist = {}
        listindex = []
        for index in indices:
            filelist[repr(index)] = '%s%s%03d%s' % \
                                    (fileroot, prefix, index, suffix)
            listindex.append(repr(index))

        val['filelist'] = filelist
        val['listindex'] = listindex
        if self.verbose:
            print "registered fileset: " + repr(val)

        self._set_key(key, val)

    def register_maprun(self, key, group_key, parent, field_tag, desc,
                        notes=None, status=None, sectag="sec",
                        skip_firstpass=False):
        r"""register all the files produced in a mapping run
        """
        val = {'desc': desc}
        self.groups[group_key].append(key)
        fileroot = self.fetch_path(parent)

        val['group_key'] = group_key
        val['parent'] = parent

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        if skip_firstpass:
            sections = ['A', 'B', 'C', 'D']
        else:
            sections = ['A', 'B', 'C', 'D', 'firstpass']
        #suffixes = ['history.hist', 'noise_diag_I.npy', 'clean_map_I.npy',
        #            'params.ini', 'noise_inv_I.npy', 'dirty_map_I.npy']
        #suffixdesc = ['log', 'noise_diag', 'clean_map', 'params', 'noise_inv',
        #          'dirty_map']
        suffixes = ['noise_diag_I.npy', 'clean_map_I.npy',
                    'params.ini', 'noise_inv_I.npy', 'dirty_map_I.npy']
        suffixdesc = ['noise_diag', 'clean_map', 'params', 'noise_inv',
                  'dirty_map']

        filelist = {}
        listindex = []

        for sec in sections:
            for (suffix, sdesc) in zip(suffixes, suffixdesc):
                filekey = "%s;%s" % (sec, sdesc)
                if (sec != 'firstpass'):
                    filelist[filekey] = "%s%s_%s_%s_%s" % \
                                        (fileroot, sectag, sec, field_tag, suffix)
                    listindex.append(filekey)
                else:
                    if suffix != "noise_diag_I.npy":
                        filelist[filekey] = "%s%s_%s" % \
                                            (fileroot, field_tag, suffix)
                        listindex.append(filekey)

        val['listindex'] = listindex
        val['filelist'] = filelist
        if self.verbose:
            print "registered map run set: " + repr(val)

        self._set_key(key, val)

    # secC_15hr_41-90_noise_inv_I_799.npy
    # vs old: sec_B_15hr_41-90_dirty_map_I.npy
    # /mnt/raid-project/gmrt/kiyo/gbt_out/maps/dec29.2011/secA_15hr_41-90_chol_I_737.npy
    # /mnt/raid-project/gmrt/kiyo/gbt_out/maps/dec29.2011/secA_15hr_41-90_clean_map_I_737.npy
    # /mnt/raid-project/gmrt/kiyo/gbt_out/maps/dec29.2011/secA_15hr_41-90_dirty_map_I_737.npy
    # /mnt/raid-project/gmrt/kiyo/gbt_out/maps/dec29.2011/secA_15hr_41-90_noise_diag_I_737.npy
    # /mnt/raid-project/gmrt/kiyo/gbt_out/maps/dec29.2011/secA_15hr_41-90_noise_inv_I_737.npy
    def register_optimalmap_section_run(self, key, group_key, parent, field_tag, band, desc,
                                        notes=None, status=None,
                                        no_weight=False):
        r"""register all the files produced in a mapping run
        """
        val = {'desc': desc}
        self.groups[group_key].append(key)
        fileroot = self.fetch_path(parent)

        val['group_key'] = group_key
        val['parent'] = parent

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        sections = ['A', 'B', 'C', 'D']
        if no_weight:
            suffixes = ['clean_map_I', 'dirty_map_I',
                        'noise_diag_I', 'noise_inv_I']
            suffixdesc = ['clean_map', 'dirty_map',
                          'noise_diag', 'noise_inv']
        else:
            suffixes = ['clean_map_I', 'dirty_map_I',
                        'noise_diag_I', 'noise_inv_I',
                        'noise_weight_I']
            suffixdesc = ['clean_map', 'dirty_map',
                          'noise_diag', 'noise_inv',
                          'noise_weight']

        filelist = {}
        listindex = []

        for sec in sections:
            for (suffix, sdesc) in zip(suffixes, suffixdesc):
                filekey = "%s;%s" % (sec, sdesc)
                filelist[filekey] = "%ssec%s_%s_%s_%s.npy" % \
                                    (fileroot, sec, field_tag, suffix, band)
                listindex.append(filekey)

        val['listindex'] = listindex
        val['filelist'] = filelist
        if self.verbose:
            print "registered map run set: " + repr(val)

        self._set_key(key, val)


    #/mnt/raid-project/gmrt/kiyo/gbt_out/maps/jan16.2012/secA_15hr_41-90_clean_map_I_all.npy
    #/mnt/raid-project/gmrt/kiyo/gbt_out/maps/jan16.2012/secA_15hr_41-90_noise_diag_I_all.npy
    #/mnt/raid-project/gmrt/kiyo/gbt_out/maps/jan16.2012/secA_15hr_41-90_noise_inv_diag_I_all.npy
    def register_optimalmap_glued_run(self, key, group_key, parent, field_tag, desc,
                                        notes=None, status=None):
        r"""register all the files produced in a mapping run
        """
        val = {'desc': desc}
        self.groups[group_key].append(key)
        fileroot = self.fetch_path(parent)

        val['group_key'] = group_key
        val['parent'] = parent

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        sections = ['A', 'B', 'C', 'D']
        suffixes = ['clean_map_I', 'noise_diag_I', 'noise_inv_diag_I']
        suffixdesc = ['clean_map', 'noise_diag', 'noise_inv']

        filelist = {}
        listindex = []

        for sec in sections:
            for (suffix, sdesc) in zip(suffixes, suffixdesc):
                filekey = "%s;%s" % (sec, sdesc)
                filelist[filekey] = "%ssec%s_%s_%s_all.npy" % \
                                    (fileroot, sec, field_tag, suffix)
                listindex.append(filekey)

        val['listindex'] = listindex
        val['filelist'] = filelist
        if self.verbose:
            print "registered glued optimal map run set: " + repr(val)

        self._set_key(key, val)


    def register_fourway_list(self, key, group_key, parent, desc, modelist,
                 notes=None, status=None, paramfile="params.ini", tag="",
                 register_modes=True, register_pickles=False,
                 register_corrsvd=True, register_separate_weights=False):
        r"""make a database set for map pair cleaning runs
        a typical run and key pairs might be:
        X_with_Y;map;#modes       sec_X_cleaned_clean_map_I_with_Y_#modes.npy
        X_with_Y;noise_inv;#modes sec_X_cleaned_noise_inv_I_with_Y_#modes.npy
        X_with_Y;weight;#modes    sec_X_cleaned_weight_I_with_Y_#modes.npy
        X_with_Y;modes;#modes     sec_X_modes_clean_map_I_with_Y_#modes.npy
        X_with_Y;fore_corr foreground_corr_pair_X_with_Y.pkl
        X_with_Y;SVD       SVD_pair_X_with_Y.pkl
        param              params.ini
        The files sec_X_15hr_41-90_noise_inv_I_diag.npy which are cached by the
        cleaning algorithm are not tracked.
        """
        val = {'desc': desc}
        self.groups[group_key].append(key)
        fileroot = self.fetch_path(parent)

        val['group_key'] = group_key
        val['parent'] = parent

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        filelist = {}
        listindex = []

        pairs = [('A', 'B'), ('A', 'C'), ('A', 'D'),
                 ('B', 'A'), ('B', 'C'), ('B', 'D'),
                 ('C', 'A'), ('C', 'B'), ('C', 'D'),
                 ('D', 'A'), ('D', 'B'), ('D', 'C')]
        corrpairs = [('A', 'B'), ('A', 'C'), ('A', 'D'),
                     ('B', 'C'), ('B', 'D'), ('C', 'D')]

        for mode_num in modelist:
            modetag = "%dmodes" % mode_num
            for (left, right) in pairs:
                pairname = "%s_with_%s" % (left, right)
                pairmap = pairname + ";map;%s" % modetag
                pairnoise = pairname + ";noise_inv;%s" % modetag
                listindex.extend([pairmap, pairnoise])
                prefix = "%ssec_%s_%s" % (fileroot, left, tag)
                suffix = "_with_%s_%s.npy" % (right, modetag)
                filelist[pairmap] = "%scleaned_clean_map_I%s" % (prefix, suffix)
                filelist[pairnoise] = "%scleaned_noise_inv_I%s" % (prefix, suffix)

                if register_separate_weights:
                    pairweight = pairname + ";weight;%s" % modetag
                    listindex.extend([pairweight])
                    filelist[pairweight] = "%scleaned_weight_I%s" % (prefix, suffix)

                if register_modes:
                    pairmodes = pairname + ";modes;%s" % modetag
                    listindex.extend([pairmodes])
                    filelist[pairmodes] = "%smodes_clean_map_I%s" % \
                                          (prefix, suffix)

        if register_corrsvd:
            for (left, right) in corrpairs:
                pairname = "%s_with_%s" % (left, right)
                paircorr = pairname + ";fore_corr"
                pair_svd = pairname + ";SVD"
                listindex.extend([paircorr, pair_svd])
                filelist[paircorr] = \
                        "%sforeground_corr_pair_%s%s_with_%s.pkl" % \
                        (fileroot, tag, left, right)

                filelist[pair_svd] = "%sSVD_pair_%s%s_with_%s.pkl" % \
                                      (fileroot, tag, left, right)

        if paramfile is not None:
            listindex.append('param')
            filelist['param'] = fileroot + paramfile

        val['listindex'] = listindex
        val['filelist'] = filelist
        if self.verbose:
            print "registered fourway list: " + repr(val)

        self._set_key(key, val)

    def register_combined_maprun(self, key, group_key, parent, desc,
                                 modelist, notes=None, status=None):
        r"""register all the files produced in a mapping run
        """
        val = {'desc': desc}
        self.groups[group_key].append(key)
        fileroot = self.fetch_path(parent)

        val['group_key'] = group_key
        val['parent'] = parent

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        filelist = {}
        listindex = []

        for mode_num in modelist:
            modetag = "%dmodes" % mode_num

            filekey = 'map;%s' % modetag
            filename = '%scombined_clean_map_%s.npy' % (fileroot, modetag)
            filelist[filekey] = filename
            listindex.append(filekey)

            filekey = 'product;%s' % modetag
            filename = '%scombined_clean_product_%s.npy' % (fileroot, modetag)
            filelist[filekey] = filename
            listindex.append(filekey)

            filekey = 'weight;%s' % modetag
            filename = '%scombined_clean_weight_%s.npy' % (fileroot, modetag)
            filelist[filekey] = filename
            listindex.append(filekey)

            filekey = 'ones;%s' % modetag
            filename = '%scombined_clean_ones_%s.npy' % (fileroot, modetag)
            filelist[filekey] = filename
            listindex.append(filekey)

        val['listindex'] = listindex
        val['filelist'] = filelist
        if self.verbose:
            print "registered combined map run set: " + repr(val)

        self._set_key(key, val)

    def register_combined_batchsimrun(self, key, group_key, parent, desc,
                                      modelist, simlist, notes=None, status=None):
        r"""register all the files produced in a mapping run
        """
        val = {'desc': desc}
        self.groups[group_key].append(key)
        fileroot = self.fetch_path(parent)

        val['group_key'] = group_key
        val['parent'] = parent

        if notes is not None:
            val['notes'] = notes

        if status is not None:
            val['status'] = status

        filelist = {}
        listindex = []

        for sim_num in simlist:
            for mode_num in modelist:
                modetag = "%dmodes" % mode_num
                simtag = "sim%d" % sim_num

                filekey = 'map;%s;%s' % (modetag, simtag)
                filename = '%scombined_%s_clean_map_%s.npy' % \
                            (fileroot, simtag, modetag)

                filelist[filekey] = filename
                listindex.append(filekey)

                filekey = 'product;%s;%s' % (modetag, simtag)
                filename = '%scombined_%s_clean_product_%s.npy' % \
                            (fileroot, simtag, modetag)

                filelist[filekey] = filename
                listindex.append(filekey)

                filekey = 'weight;%s;%s' % (modetag, simtag)
                filename = '%scombined_%s_clean_weight_%s.npy' % \
                            (fileroot, simtag, modetag)

                filelist[filekey] = filename
                listindex.append(filekey)

                filekey = 'ones;%s;%s' % (modetag, simtag)
                filename = '%scombined_%s_clean_ones_%s.npy' % \
                            (fileroot, simtag, modetag)

                filelist[filekey] = filename
                listindex.append(filekey)

        val['listindex'] = listindex
        val['filelist'] = filelist
        if self.verbose:
            print "registered combined map run set: " + repr(val)

        self._set_key(key, val)

