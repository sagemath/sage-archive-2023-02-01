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

cdef extern from *:
    bint PyGen_Check(x)


def running_total(L, start=None):
    """
    Returns a list where the i-th entry is the sum of all entries up to (and incling) i.

    INPUT:
        L     -- the list
        start -- (optional) a default start value

    EXAMPLES:
        sage: running_total(range(5))
        [0, 1, 3, 6, 10]
        sage: running_total("abcdef")
        ['a', 'ab', 'abc', 'abcd', 'abcde', 'abcdef']
        sage: running_total([1..10], start=100)
        [101, 103, 106, 110, 115, 121, 128, 136, 145, 155]
    """
    cdef bint first = 1
    for x in L:
        if first:
            total = L[0] if start is None else L[0]+start
            running = [total]
            first = 0
            continue
        total += x
        PyList_Append(running, total)
    return running


def prod(x, z=None, Py_ssize_t recursion_cutoff = 5):
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
        Robert Bradshaw (2007-10-26) -- Balanced product tree, other optimizations, (lazy) generator support
    """
    if not PyList_CheckExact(x) and not PyTuple_CheckExact(x):

        if PyGen_Check(x):
            # lazy list, do lazy product
            try:
                prod = x.next() if z is None else z * x.next()
                for a in x:
                    prod *= a
                return prod
            except StopIteration:
                x = []

        else:

            try:
                return x.prod()
            except AttributeError:
                pass

            try:
                return x.mul()
            except AttributeError:
                pass

            x = list(x)

    cdef Py_ssize_t n = len(x)

    if n == 0:
        if z is None:
            import sage.rings.integer
            return sage.rings.integer.Integer(1)
        else:
            return z

    prod = balanced_list_prod(x, 0, n, recursion_cutoff)

    if z is not None:
        prod = z*prod

    return prod


cdef balanced_list_prod(L, Py_ssize_t offset, Py_ssize_t count, Py_ssize_t cutoff):
    """
    INPUT:
        L      -- the terms (MUST be a tuple or list)
        off    -- offset in the list from which to start
        count  -- how many terms in the product
        cutoff -- the minimum count to recurse on

    OUTPUT:
        L[offset] * L[offset+1] * ... * L[offset+count-1]

    NOTE: The parameter cutoff must be at least 1, and there is no reason to
          ever make it less than 3. However, there are at least two advantages
          to setting it higher (and consequently not recursing all the way
          down the tree). First, one avoids the overhead of the function
          calls at the base of the tree (which is the majority of them) and
          second, it allows one to save on object creation if inplace
          operations are used. The asymptotic gains should usually be at the
          top of the tree anyway.
    """
    cdef Py_ssize_t k
    if count <= cutoff:
        prod = <object>PySequence_Fast_GET_ITEM(L, offset)
        for k from offset < k < offset+count:
            prod *= <object>PySequence_Fast_GET_ITEM(L, k)
        return prod
    else:
        k = (1+count) >> 1
        return balanced_list_prod(L, offset, k, cutoff) * balanced_list_prod(L, offset+k, count-k, cutoff)


class NonAssociative:
    """
    This class is to test the balance nature of prod.

    EXAMPLES:
        sage: from sage.misc.misc_c import NonAssociative
        sage: L = [NonAssociative(label) for label in 'abcdef']
        sage: prod(L)
        (((a*b)*c)*((d*e)*f))
        sage: L = [NonAssociative(label) for label in range(20)]
        sage: prod(L, recursion_cutoff=5)
        ((((((0*1)*2)*3)*4)*((((5*6)*7)*8)*9))*(((((10*11)*12)*13)*14)*((((15*16)*17)*18)*19)))
        sage: prod(L, recursion_cutoff=1)
        (((((0*1)*2)*(3*4))*(((5*6)*7)*(8*9)))*((((10*11)*12)*(13*14))*(((15*16)*17)*(18*19))))
        sage: L = [NonAssociative(label) for label in range(14)]
        sage: prod(L, recursion_cutoff=1)
        ((((0*1)*(2*3))*((4*5)*6))*(((7*8)*(9*10))*((11*12)*13)))
    """
    def __init__(self, left, right=None):
        """
        EXAMPLES:
            sage: from sage.misc.misc_c import NonAssociative
            sage: NonAssociative('a')
            a
            sage: NonAssociative('a','b')
            (a*b)
        """
        self.left = left
        self.right = right

    def __repr__(self):
        """
        EXAMPLES:
            sage: from sage.misc.misc_c import NonAssociative
            sage: NonAssociative(1)
            1
            sage: NonAssociative(2,3)
            (2*3)
        """
        if self.right is None:
            return str(self.left)
        else:
            return "(%s*%s)" % (self.left, self.right)

    def __mul__(self, other):
        """
        EXAMPLES:
            sage: from sage.misc.misc_c import NonAssociative
            sage: a, b, c = [NonAssociative(label) for label in 'abc']
            sage: (a*b)*c
            ((a*b)*c)
            sage: a*(b*c)
            (a*(b*c))
        """
        return NonAssociative(self, other)
