"""Code that is shared between fitsGBT.py and fits_map.py.

Since Reading a map and a table for a fits file is so different, the fits
handling code for maps and tables cannot inherit form a common class (without
be doing an insane amount of work).  Here I simply provide some functions that
overlap to reduce code redundancy as best I can.
"""

import scipy as sp
import numpy.ma as ma
import pyfits

import kiyopy.custom_exceptions as ce
import base_data as bd

# These globals are the cards for (our custom) history entries in a fits header
card_hist = 'DB-HIST'
card_detail = 'DB-DET'

def get_history_header(prihdr) :
    """Gets the history from a pyfits primary header.
    
    This function accepts the primary header of a pyfits hdulist and reads
    the data 'history' from it.  This is the history that is tracked by this
    code, with cards DB-HIST and DB-DET (not the normal fits HISTORY cards).
    """
    
    # Initialize a blank history object
    history = bd.History()
    # If there is no history, return.
    try :
        ii = prihdr.ascardlist().index_of(card_hist)
    except KeyError :
        return history
    n_cards = len(prihdr.ascardlist().keys())
    while ii < n_cards :
        if prihdr.ascardlist().keys()[ii] == card_hist :
            hist_entry = prihdr[ii]
            details = []
        elif prihdr.ascardlist().keys()[ii] == card_detail :
            details.append(prihdr[ii])
        ii = ii + 1
        if ii == n_cards or prihdr.ascardlist().keys()[ii] == card_hist :
            history.add(hist_entry, details)

    return history

def write_history_header(prihdr, history) :
    """Puts a puts a data history into a pyfits header.

    history is a bd.History object, that is stored at the end of the pyfits
    header using the DB-HIST and DB-DET cards.
    """

    history_keys  = history.keys()
    history_keys.sort()
    for hist in history_keys :
        details = history[hist]
        # Chop off the number, since they are already sorted.
        hcard = pyfits.Card(card_hist, hist[5:])
        prihdr.ascardlist().append(hcard)
        for detail in details :
            dcard = pyfits.Card(card_detail, detail)
            prihdr.ascardlist().append(dcard)



