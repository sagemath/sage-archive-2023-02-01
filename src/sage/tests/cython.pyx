"""
This file collects tests requiring Cython.
"""
#*****************************************************************************
#       Copyright (C) 2012 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#       Copyright (C) 2012 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_singleton cimport FastHashable_class
cdef class ClassWithLargeHash(FastHashable_class):
    """
    This class tests against a bug with :class:`FastHashable_class`
    (in an earlier version of the patch at Trac #11900) that occurred
    on systems where ``sys.maxint`` does not fit into a C int.

    TESTS::

        sage: import sage.tests.cython
        sage: C = sage.tests.cython.ClassWithLargeHash(); C
        A successfully created object with a very large hash
        sage: hash(C) == sys.maxint
        True
    """
    def __init__(self):
        import sys
        self._hash = sys.maxint
    def __repr__(self):
        return 'A successfully created object with a very large hash'
