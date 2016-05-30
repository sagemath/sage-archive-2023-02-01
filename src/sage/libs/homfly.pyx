r"""
Cython wrapper for libhomfly library


AUTHORS:

- Miguel Marco (2015-03-24): initial version.


This is used to call the libhomfly library directly from python. Knots
and Links are passed following the convention in libhomfly. It is basically
the oriented Gauss code, represented as a string of integers separated
by spaces as follows:

- how many strings,

    - for each string, how many crossings, then

        - for each crossing, the cross name, then 1 if over, -1 if under

- for each crossing, the name of the crossing and 1 if right, -1 if left.

If there are n crossings, they must be named 0..n-1.
"""

#*****************************************************************************
#       Copyright (C) 2015 Miguel Marco  <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#clib homfly
#clib gc

include 'cysignals/signals.pxi'

cdef extern from "homfly.h":
    char* homfly(char *argv)

def homfly_polynomial(link):
    """
    Return the HOMFLY polynomial of a link.

    INPUT:

    - ``link`` -- a string of space-separated integers representing the link.

    OUTPUT:

    A string with the HOMFLY polynomial in the variables `M` and `L`

    EXAMPLES::

        sage: from sage.libs.homfly import homfly_polynomial
        sage: trefoil = '1 6 0 1  1 -1  2 1  0 -1  1 1  2 -1 0 1 1 1 2 1'
        sage: homfly_polynomial(trefoil) # optional - libhomfly
        ' - L^-4 - 2L^-2 + M^2L^-2'

    """
    cdef char* c_string = link
    sig_on()
    cdef char* c_output = homfly(c_string)
    sig_off()
    output = <bytes> c_output
    return output
