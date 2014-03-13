"""History class for keeping track of what has been done to a set of files."""

import kiyopy.custom_exceptions as ce

class History(dict) :
    """Class that contains the history of a piece of data."""

    def add(self, history_entry, details=()) :

        local_details = details
        # Input checks.
        #if len(history_entry) > 70 :
        #    raise ValueError('History entries limited to 70 characters.')
        if type(details) is str :
        #    if len(details) > 70 :
        #        raise ValueError('History details limited to 70 characters.')
            local_details = (details, )
        for detail in details :
            if not type(detail) is str :
                raise TypeError('History details must be a squence of strings'
                                ' or a single string.')
            #if len(detail) > 70 :
            #    raise ValueError('History details limited to 70 characters.')

        n_entries = len(self)
        # '+' operator performs input type check.
        hist_str = ('%03d: ' % n_entries) + history_entry
        self[hist_str] = tuple(local_details)

    def display(self) :
        """Prints the data history in human readable format."""
    
        history_keys = self.keys()
        history_keys.sort()
        for history in history_keys :
            details = self[history]
            print history
            for detail in details :
                print '    ' + detail

    def merge(self, *args) :
        """Merge this History object with ones passed to this function."""

        for obj in args :
            if hasattr(obj, 'history') :
                thishistory = obj.history
            else :
                thishistory = obj
            for entry, details in thishistory.iteritems() :
                for detail in details :
                    try :
                        if not detail in self[entry] :
                            self[entry] = self[entry] + (detail, )
                    except KeyError :
                        raise ce.DataError("Histories to be merged must have"
                                               " identical keys.")

    def write(self, fname) :
        """Write this history to disk"""

        f = open(fname, 'w')
        
        try :
            f.write('{\n')
            keys = self.keys()
            keys.sort()
            for history in keys :
                f.write(repr(history))
                f.write(' : (\n    ')
                for detail in self[history] :
                    f.write(repr(detail))
                    f.write(',\n    ')
                f.write('),\n')
            f.write('}\n')
        finally :
            f.close()

def read(fname) :
    """Read a History object from file."""

    f = open(fname, 'r')
    try :
        filestring = f.read()
    finally :
        f.close()
    return History(eval(filestring))

def merge_histories(*args) :
    """Merges DataBlock histories.

    This function accepts an arbitray number of History objects (or classes
    containing a history object in a 'history' attribute), and returns a 
    history dictionary that is a merger of them.  History keys must match; 
    details are added."""
    
    if hasattr(args[0], 'history') :
        history = History(args[0].history)
    else :
        history = History(args[0])
    history.merge(*args[1:])
        
    return history

