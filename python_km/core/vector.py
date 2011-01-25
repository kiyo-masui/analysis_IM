"""This module provides vector and matrix interfaces tylored to the needs of
intensity mapping analysis.

The considerations for intensity mapping analysis can be summarized as follows:
    - Elements in data vectors often have a significant amount of assotiated
      meta data.  For instance in map space, a data point is a pixel and it has
      pointing and frequency information that needs to be stored with it.
    - Matricies can generally refer to data vectors for meta data.
    - Data vectors can normally be held in the memory of a single machine, but
      matricies cannot.
    - Matricies are occationally Block diagonal and we would like to take
      advantage of this.
    - We will evenetually need to do operations massivly in parralelle with
      SCALLAPACK.

To satisfy these needs, this module provides:
    - A data structure for vectors and infrastructure to read and write them to
      disk.  This is implemented using numpu record arrays. Memory mapping is
      fully supported.
    - Infrastructure for using numpy array objects as matricies, with memory
      mapping a being a consideration throughout.
    - To take advantage of block diagonal matricies, all matricies are 3 
      dimentional arrays, with the 0th dimention indicating the block.  The 0th
      dimension may be of lenth 1 (for a dense array).
    - Infrastructure for maniputlating these formats and performing standard
      operations.
"""

import scipy as sp
import numpy.lib.format as format

import kiyopy.custom_exceptions as ce


