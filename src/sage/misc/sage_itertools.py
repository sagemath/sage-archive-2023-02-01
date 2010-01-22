"""
Miscellaneous functions which should eventually be moved upstream into
Python's standard itertools module.
"""
#*****************************************************************************
#  Copyright (C) 2010     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

import itertools, heapq

def unique_merge(*lists):
    """
    INPUT:

     - ``lists``: sorted lists (or iterables)

    Return an iterator over the elements of each list in ``lists``, in
    sorted order, with duplicates removed.

        sage: from sage.misc.sage_itertools import unique_merge
        sage: list(unique_merge([1,2,2,3,4,7,9], [0,2,4], [2,5]))
        [0, 1, 2, 3, 4, 5, 7, 9]

    Inspired from: http://rosettacode.org/wiki/Create_a_Sequence_of_unique_elements#Python
    """
    return (k for k,g in itertools.groupby(heapq.merge(*lists)))

