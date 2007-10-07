#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../ext/stdsage.pxi"
include "../ext/python_sequence.pxi"
include "../ext/python_list.pxi"
include "../ext/python_tuple.pxi"

def prod(x, z=None):
    """
    Return the product of the elements in the list x.  If optional
    argument z is not given, start the product with the first element
    of the list, otherwise use z.  The empty product is the int 1 if z
    is not specified, and is z if given.

    This assumes that your multiplication is associative; we don't promise
    which end of the list we start at.

    EXAMPLES:
        sage: prod([1,2,34])
        68
        sage: prod([2,3], 5)
        30
        sage: prod((1,2,3), 5)
        30
        sage: F = factor(-2006); F
        -1 * 2 * 17 * 59
        sage: prod(F)
        -2006

    AUTHORS:
        Joel B. Mohler (2007-10-03 -- Reimplemented in Cython and optimized)
    """
    if not PyList_CheckExact(x) and not PyTuple_CheckExact(x):
        try:
            return x.prod()
        except AttributeError:
            try:
                return x.mul()
            except AttributeError:
                pass

        x = list(x)

    cdef Py_ssize_t j
    cdef Py_ssize_t i = 0

    if z is None:
        if len(x) == 0:
            import sage.rings.integer
            return sage.rings.integer.Integer(1)
        z = x[0]
        i = 1

    # TODO: Change this to use a balanced tree in some cases, e.g.,
    # if input is a list.

    # DEFINITELY:  Change it to a balanced tree.  It is vastly faster,
    # e.g., when multiplying a list of integers to compute n!, doing
    # a balanced tree makes it possible to use asymptotic multiplication,
    # which is way way faster.  -- William Stein
    #

    for j from i <= j < len(x):
        z *= x[j]
    return z
