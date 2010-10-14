"""Data container for maps."""

import base_data

class DataMap(base_data.BaseData) :
    """This class holds data from a map."""

    axes = ('long', 'lat', 'freq')


# Clone some extra functions:

merge_histories = base_data.merge_histories


print_history = base_data.print_history

