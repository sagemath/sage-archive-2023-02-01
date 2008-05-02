"""
Matrix Constructor.
"""


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

import sage.rings.all as rings
import sage.matrix.matrix_space as matrix_space
from sage.structure.element import is_Vector
from sage.structure.sequence import Sequence
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.integer_ring import ZZ
from sage.misc.misc_c import running_total
from matrix import is_Matrix

import sage.categories.pushout

def matrix(*args, **kwds):
    """
    Create a matrix.

    INPUT:

        The matrix command takes the entries of a matrix, optionally
        preceded by a ring and the dimensions of the matrix, and returns
        a matrix.

        The entries of a matrix can be specified as a flat list of
        elements, a list of lists (i.e., a list of rows), a list of
        Sage vectors, or a dictionary having positions as keys and
        matrix entries as values (see the examples).  You can create a
        matrix of zeros by passing an empty list or the integer zero
        for the entries.  To construct a multiple of the identity
        ($cI$), you can specify square dimensions and pass in $c$.
        Calling matrix() with a Sage object may return something that
        makes sense.  Calling matrix() with a numpy array will convert
        the array to a matrix.

        The ring, number of rows, and number of columns of the matrix
        can be specified by setting the ring, nrows, or ncols
        parameters or by passing them as the first arguments to the
        function in the order ring, nrows, ncols.  The ring defaults
        to ZZ if it is not specified or cannot be determined from the
        entries.  If the numbers of rows and columns are not specified
        and cannot be determined, then an empty 0x0 matrix is
        returned.

        ring -- the base ring for the entries of the matrix.

        nrows -- the number of rows in the matrix.

        ncols -- the number of columns in the matrix.

        sparse -- create a sparse matrix.  This defaults to True when
        the entries are given as a dictionary, otherwise defaults to False.


    OUTPUT:

        a matrix

    EXAMPLES:
        sage: m=matrix(2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

        sage: m=matrix(2,3); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

        sage: m=matrix(QQ,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        sage: v1=vector((1,2,3))
        sage: v2=vector((4,5,6))
        sage: m=matrix([v1,v2]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

        sage: m=matrix(QQ,2,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        sage: m=matrix(QQ,2,3,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        sage: m=matrix({(0,1): 2, (1,1):2/5}); m; m.parent()
        [  0   2]
        [  0 2/5]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field

        sage: m=matrix(QQ,2,3,{(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field

        sage: import numpy
        sage: n=numpy.array([[1,2],[3,4]],float)
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0]
        [3.0 4.0]
        Full MatrixSpace of 2 by 2 dense matrices over Real Double Field

        sage: v = vector(ZZ, [1, 10, 100])
        sage: m=matrix(v); m; m.parent()
        [  1  10 100]
        Full MatrixSpace of 1 by 3 dense matrices over Integer Ring
        sage: m=matrix(GF(7), v); m; m.parent()
        [1 3 2]
        Full MatrixSpace of 1 by 3 dense matrices over Finite Field of size 7

        sage: g = graphs.PetersenGraph()
        sage: m = matrix(g); m; m.parent()
        [0 1 0 0 1 1 0 0 0 0]
        [1 0 1 0 0 0 1 0 0 0]
        [0 1 0 1 0 0 0 1 0 0]
        [0 0 1 0 1 0 0 0 1 0]
        [1 0 0 1 0 0 0 0 0 1]
        [1 0 0 0 0 0 0 1 1 0]
        [0 1 0 0 0 0 0 0 1 1]
        [0 0 1 0 0 1 0 0 0 1]
        [0 0 0 1 0 1 1 0 0 0]
        [0 0 0 0 1 0 1 1 0 0]
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring

        sage: matrix(ZZ, 10, 10, range(100), sparse=True).parent()
        Full MatrixSpace of 10 by 10 sparse matrices over Integer Ring

        sage: R = PolynomialRing(QQ, 9, 'x')
        sage: A = matrix(R, 3, 3, R.gens()); A
        [x0 x1 x2]
        [x3 x4 x5]
        [x6 x7 x8]
        sage: det(A)
        -x2*x4*x6 + x1*x5*x6 + x2*x3*x7 - x0*x5*x7 - x1*x3*x8 + x0*x4*x8


    TESTS:
        sage: m=matrix(); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m=matrix(QQ); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m=matrix(QQ,2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: m=matrix(QQ,2,3); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m=matrix([]); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m=matrix(QQ,[]); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m=matrix(2,2,1); m; m.parent()
        [1 0]
        [0 1]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: m=matrix(QQ,2,2,1); m; m.parent()
        [1 0]
        [0 1]
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: m=matrix(2,3,0); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m=matrix(QQ,2,3,0); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m=matrix([[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m=matrix(QQ,2,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m=matrix(QQ,3,[[1,2,3],[4,5,6]]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: Number of rows does not match up with specified number.
        sage: m=matrix(QQ,2,3,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m=matrix(QQ,2,4,[[1,2,3],[4,5,6]]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: Number of columns does not match up with specified number.
        sage: m=matrix([(1,2,3),(4,5,6)]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m=matrix([1,2,3,4,5,6]); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Integer Ring
        sage: m=matrix((1,2,3,4,5,6)); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Integer Ring
        sage: m=matrix(QQ,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Rational Field
        sage: m=matrix(QQ,3,2,[1,2,3,4,5,6]); m; m.parent()
        [1 2]
        [3 4]
        [5 6]
        Full MatrixSpace of 3 by 2 dense matrices over Rational Field
        sage: m=matrix(QQ,2,4,[1,2,3,4,5,6]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: entries has the wrong length
        sage: m=matrix(QQ,5,[1,2,3,4,5,6]); m; m.parent()
        Traceback (most recent call last):
        ...
        TypeError: entries has the wrong length
        sage: m=matrix({(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        sage: m=matrix(QQ,{(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m=matrix(QQ,3,{(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        [0 0 0]
        Full MatrixSpace of 3 by 3 sparse matrices over Rational Field
        sage: m=matrix(QQ,3,4,{(1,1): 2}); m; m.parent()
        [0 0 0 0]
        [0 2 0 0]
        [0 0 0 0]
        Full MatrixSpace of 3 by 4 sparse matrices over Rational Field
        sage: m=matrix(QQ,2,{(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m=matrix(QQ,1,{(1,1): 2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list
        sage: m=matrix({}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m=matrix(QQ,{}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Rational Field
        sage: m=matrix(QQ,2,{}); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m=matrix(QQ,2,3,{}); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field
        sage: m=matrix(2,{}); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        sage: m=matrix(2,3,{}); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Integer Ring
        sage: m=matrix(0); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m=matrix(0,2); m; m.parent()
        []
        Full MatrixSpace of 0 by 2 dense matrices over Integer Ring
        sage: m=matrix(2,0); m; m.parent()
        []
        Full MatrixSpace of 2 by 0 dense matrices over Integer Ring
        sage: m=matrix(0,[1]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: entries has the wrong length
        sage: m=matrix(1,0,[]); m; m.parent()
        []
        Full MatrixSpace of 1 by 0 dense matrices over Integer Ring
        sage: m=matrix(0,1,[]); m; m.parent()
        []
        Full MatrixSpace of 0 by 1 dense matrices over Integer Ring
        sage: m=matrix(0,[]); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m=matrix(0,{}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m=matrix(0,{(1,1):2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list
        sage: m=matrix(2,0,{(1,1):2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list
        sage: import numpy
        sage: n=numpy.array([[numpy.complex(0,1),numpy.complex(0,2)],[3,4]],complex)
        sage: m=matrix(n); m; m.parent()
        [1.0*I 2.0*I]
        [  3.0   4.0]
        Full MatrixSpace of 2 by 2 dense matrices over Complex Double Field
        sage: n=numpy.array([[1,2],[3,4]],'int32')
        sage: m=matrix(n); m; m.parent()
        [1 2]
        [3 4]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: n = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'float32')
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Real Double Field
        sage: n=numpy.array([[1,2,3],[4,5,6],[7,8,9]],'float64')
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Real Double Field
        sage: n=numpy.array([[1,2,3],[4,5,6],[7,8,9]],'complex64')
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Complex Double Field
        sage: n=numpy.array([[1,2,3],[4,5,6],[7,8,9]],'complex128')
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Complex Double Field
        sage: a = matrix([[1,2],[3,4]])
        sage: b = matrix(a.numpy()); b
        [1 2]
        [3 4]
        sage: a == b
        True
        sage: c = matrix(a.numpy('float32')); c
        [1.0 2.0]
        [3.0 4.0]
        sage: v = vector(ZZ, [1, 10, 100])
        sage: m=matrix(ZZ['x'], v); m; m.parent()
        [  1  10 100]
        Full MatrixSpace of 1 by 3 dense matrices over Univariate Polynomial Ring in x over Integer Ring
        sage: matrix(ZZ, 10, 10, range(100)).parent()
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring
        sage: m = matrix(GF(7), [[1/3,2/3,1/2], [3/4,4/5,7]]); m; m.parent()
        [5 3 4]
        [6 5 0]
        Full MatrixSpace of 2 by 3 dense matrices over Finite Field of size 7
        sage: m = matrix([[1,2,3], [RDF(2), CDF(1,2), 3]]); m; m.parent()
        [        1.0         2.0         3.0]
        [        2.0 1.0 + 2.0*I         3.0]
        Full MatrixSpace of 2 by 3 dense matrices over Complex Double Field
        sage: m=matrix(3,3,1/2); m; m.parent()
        [1/2   0   0]
        [  0 1/2   0]
        [  0   0 1/2]
        Full MatrixSpace of 3 by 3 dense matrices over Rational Field
        sage: matrix([[1],[2,3]])
        Traceback (most recent call last):
        ...
        ValueError: List of rows is not valid (rows are wrong types or lengths)
        sage: matrix([[1],2])
        Traceback (most recent call last):
        ...
        ValueError: List of rows is not valid (rows are wrong types or lengths)
        sage: matrix(vector(RR,[1,2,3])).parent()
        Full MatrixSpace of 1 by 3 dense matrices over Real Field with 53 bits of precision

    AUTHORS:
        -- ??: Initial implementation
        -- Jason Grout: almost a complete rewrite, with bits and pieces
           from the original implementation (Mar 2008)
    """
    args = list(args)
    sparse = kwds.get('sparse',False)
    # if the first object already knows how to make itself into a
    # matrix, try that, defaulting to a matrix over the integers.
    if len(args) == 1 and hasattr(args[0], '_matrix_'):
        try:
                return args[0]._matrix_(sparse=sparse)
        except TypeError:
                return args[0]._matrix_()
    elif len(args) == 2:
        if hasattr(args[0], '_matrix_'):
            try:
                return args[0]._matrix_(args[1], sparse=sparse)
            except TypeError:
                return args[0]._matrix_(args[1])
        elif hasattr(args[1], '_matrix_'):
            try:
                return args[1]._matrix_(args[0], sparse=sparse)
            except TypeError:
                return args[1]._matrix_(args[0])

    if len(args) == 0:
        # if nothing was passed return the empty matrix over the
        # integer ring.
        return matrix_space.MatrixSpace(rings.ZZ, 0, 0, sparse=sparse)([])

    if len(args) >= 1 and rings.is_Ring(args[0]):
        # A ring is specified
        if kwds.get('ring', args[0]) != args[0]:
            raise ValueError, "Specified rings are not the same"
        else:
            ring = args[0]
            args.pop(0)
    else:
        ring = kwds.get('ring', None)

    if len(args) >= 1:
        # check to see if the number of rows is specified
        try:
            nrows = int(args[0])
            args.pop(0)
            if kwds.get('nrows', nrows) != nrows:
                raise ValueError, "Number of rows specified in two places and they are not the same"
        except TypeError:
            nrows = kwds.get('nrows', None)
    else:
        nrows = kwds.get('nrows', None)

    if len(args) >= 1:
        # check to see if additionally, the number of columns is specified
        try:
            ncols = int(args[0])
            args.pop(0)
            if kwds.get('ncols', ncols) != ncols:
                raise ValueError, "Number of columns specified in two places and they are not the same"
        except TypeError:
            ncols = kwds.get('ncols', None)
    else:
        ncols = kwds.get('ncols', None)


    # Now we've taken care of initial ring, nrows, and ncols arguments.
    # We've also taken care of the Sage object case.

    # Now the rest of the arguments are a list of
    # rows, a flat list of entries, a dict, a numpy array, or a single
    # value.
    if len(args) == 0:
        # If no entries are specified, pass back a zero matrix
        entries = 0
        entry_ring = rings.ZZ
    elif len(args) == 1:
        if isinstance(args[0], (list, tuple)):
            if len(args[0]) == 0:
                # no entries are specified, pass back the zero matrix
                entries = 0
                entry_ring = rings.ZZ
            elif isinstance(args[0][0], (list, tuple)) or is_Vector(args[0][0]):
		# Ensure we have a list of lists, each inner list having the same number of elements
		first_len = len(args[0][0])
		if not all( (isinstance(v, (list, tuple)) or is_Vector(v)) and len(v) == first_len for v in args[0]):
		    raise ValueError, "List of rows is not valid (rows are wrong types or lengths)"
                # We have a list of rows or vectors
                if nrows is None:
                    nrows = len(args[0])
                elif nrows != len(args[0]):
                    raise ValueError, "Number of rows does not match up with specified number."
                if ncols is None:
		    ncols = len(args[0][0])
                elif ncols != len(args[0][0]):
                    raise ValueError, "Number of columns does not match up with specified number."

                entries = sum([list(v) for v in args[0]], [])

            else:
                # We have a flat list; figure out nrows and ncols
                if nrows is None:
                    nrows = 1

                if nrows > 0:
                    if ncols is None:
                        ncols = len(args[0]) // nrows
                    elif ncols != len(args[0]) // nrows:
                        raise ValueError, "entries has the wrong length"
                elif len(args[0]) > 0:
                    raise ValueError, "entries has the wrong length"

                entries = args[0]

            if nrows > 0 and ncols > 0 and ring is None:
                entries, ring = prepare(entries)

        elif isinstance(args[0], dict):
            # We have a dictionary
            # default to sparse
            sparse = kwds.get('sparse', True)
            if len(args[0]) == 0:
                # no entries are specified, pass back the zero matrix
                entries = 0
            else:
                entries, entry_ring = prepare_dict(args[0])
                if nrows is None:
                    nrows = nrows_from_dict(entries)
                    ncols = ncols_from_dict(entries)
                # note that ncols can still be None if nrows is set --
                # it will be assigned nrows down below.

            # See the construction after the numpy case below.
        else:
            import numpy
            if isinstance(args[0], numpy.ndarray):
                num_array = args[0]
                str_dtype = str(num_array.dtype)

                if not( num_array.flags.c_contiguous is True or num_array.flags.f_contiguous is True):
                    raise TypeError('numpy matrix must be either c_contiguous or f_contiguous')
                if str_dtype.count('float32')==1:
                    m=matrix(RDF,num_array.shape[0],num_array.shape[1],0)
                    m._replace_self_with_numpy32(num_array)

                elif str_dtype.count('float64')==1:
                    m=matrix(RDF,num_array.shape[0],num_array.shape[1],0)
                    m._replace_self_with_numpy(num_array)

                elif str_dtype.count('complex64')==1:
                    m=matrix(CDF,num_array.shape[0],num_array.shape[1],0)
                    m._replace_self_with_numpy32(num_array)

                elif str_dtype.count('complex128')==1:
                    m=matrix(CDF,num_array.shape[0],num_array.shape[1],0)
                    m._replace_self_with_numpy(num_array)

                elif str_dtype.count('int') == 1:
                    m = matrix(ZZ, [list(row) for row in list(num_array)])

                elif str_dtype.count('object') == 1:
                    #Get the raw nested list from the numpy array
                    #and feed it back into matrix
                    try:
                        return matrix( [list(row) for row in list(num_array)])
                    except TypeError:
                        raise TypeError("cannot convert numpy matrix to SAGE matrix")
                else:
                    raise TypeError("cannot convert numpy matrix to SAGE matrix")

                if num_array.flags.c_contiguous:
                    return m
                else:
                    return m.transpose()
            elif nrows is not None and ncols is not None:
                # assume that we should just pass the thing into the
                # MatrixSpace constructor and hope for the best
                # This is not documented, but it is doctested
                if ring is None:
                    entry_ring = args[0].parent()
                entries = args[0]
            else:
                raise ValueError, "Invalid matrix constructor.  Type matrix? for help"
    else:
        raise ValueError, "Invalid matrix constructor.  Type matrix? for help"

    if nrows is None:
        nrows = 0
    if ncols is None:
        ncols = nrows


    if ring is None:
        try:
            ring = entry_ring
        except NameError:
            ring = rings.ZZ

    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse=sparse)(entries)



def prepare(w):
    if 0 == len(w):
        return Sequence([], rings.ZZ), rings.ZZ
    entries = Sequence(w)
    ring = entries.universe()
    if ring is int or ring is long:
        ring = rings.ZZ
    elif ring is float:
        ring = rings.RDF
    elif ring is complex:
        ring = rings.CDF
    elif not rings.is_Ring(ring):
        raise TypeError, "unable to find a common ring for all elements"
    return entries, ring

def prepare_dict(w):
    Z = w.items()
    X = [x for _, x in Z]
    entries, ring = prepare(X)
    return dict([(Z[i][0],entries[i]) for i in xrange(len(entries))]), ring

def nrows_from_dict(d):
    if 0 == len(d):
        return 0
    return max([0] + [ij[0] for ij in d.keys()]) + 1

def ncols_from_dict(d):
    if 0 == len(d):
        return 0
    return max([0] + [ij[1] for ij in d.keys()]) + 1

Matrix = matrix


def random_matrix(R, nrows, ncols=None, sparse=False, density=1, *args, **kwds):
    """
    Return a random matrix with entries in the ring R.

    INPUT:
        R -- a ring
        nrows -- integer; number of rows
        ncols -- (default: None); number of columns; if None defaults to nrows
        sparse -- (default; False); whether or not matrix is sparse.
        density -- integer (default: 1)
        *args, **kwds -- passed on to randomize function

    EXAMPLES:
        sage: A = random_matrix(ZZ,50,x=2^16)    # entries are up to 2^16 i size
        sage: A
        50 x 50 dense matrix over Integer Ring
    """
    if ncols is None:
        ncols = nrows
    A = matrix_space.MatrixSpace(R, nrows, ncols, sparse=sparse).zero_matrix()
    A.randomize(density=density, *args, **kwds)
    return A

def diagonal_matrix(arg0=None, arg1=None, arg2=None, sparse=None):
    """
    INPUT:
    Supported formats
        1. matrix(diagonal_entries, [sparse=True]):
               matrix with each row constructed from the list_of_rows
        2. matrix(nrows, diagonal_entries, [sparse=True]):
               matrix with each row constructed from the list_of_rows
        3. matrix(ring, diagonal_entries, [sparse=True]):
               matrix with each row constructed from the list_of_rows
        4. matrix(ring, nrows, diagonal_entries, [sparse=True]):
               matrix with given number of rows and flat list of entries
    The sparse option is optional, must be explicitly named (i.e.,
    sparse=True), and may be either True or False.

    EXAMPLES:
    Input format 1.
        sage: diagonal_matrix([1,2,3])
        [1 0 0]
        [0 2 0]
        [0 0 3]

    Input format 2.
        sage: diagonal_matrix(GF(3), [1,2,3])
        [1 0 0]
        [0 2 0]
        [0 0 0]

    Input format 3:
        sage: diagonal_matrix(3, [1,2])
        [1 0 0]
        [0 2 0]
        [0 0 0]

    Input format 4:
        sage: diagonal_matrix(GF(3), 3, [8,2])
        [2 0 0]
        [0 2 0]
        [0 0 0]
    """
    ring = None
    if isinstance(arg0, (list, tuple)):
        # Format 1
        v = arg0
        nrows = len(v)
    elif isinstance(arg0, (int, long, rings.Integer)):
        nrows = arg0
        v = arg1
    elif rings.is_Ring(arg0):
        ring = arg0
        if isinstance(arg1, (list, tuple)):
            v = arg1
            nrows = len(v)
        else:
            nrows = arg1
            v = arg2

    if isinstance(v, list):
        w = {}
        for i in range(len(v)):
            w[(i,i)] = v[i]
    else:
        w = v

    if ring is None:
        return matrix(nrows, nrows, w, sparse=sparse)
    else:
        return matrix(ring, nrows, nrows, w, sparse=sparse)


def identity_matrix(ring, n=0, sparse=False):
    r"""
    Return the $n \times n$ identity matrix over the given ring.

    The default ring is the integers.

    EXAMPLES:
        sage: M = identity_matrix(QQ, 2); M
        [1 0]
        [0 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: M = identity_matrix(2); M
        [1 0]
        [0 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: M = identity_matrix(3, sparse=True); M
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: M.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring
    """
    if isinstance(ring, (int, long, rings.Integer)):
        n = ring
        ring = rings.ZZ
    return matrix_space.MatrixSpace(ring, n, n, sparse).identity_matrix()


def zero_matrix(ring, nrows, ncols=None, sparse=False):
    """
    Return the $nrows \times ncols$ zero matrix over the given ring.

    The default ring is the integers.

    EXAMPLES:
        sage: M = zero_matrix(QQ, 2); M
        [0 0]
        [0 0]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: M = zero_matrix(2, 3); M
        [0 0 0]
        [0 0 0]
        sage: M.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: M = zero_matrix(3, 1, sparse=True); M
        [0]
        [0]
        [0]
        sage: M.parent()
        Full MatrixSpace of 3 by 1 sparse matrices over Integer Ring
    """
    if isinstance(ring, (int, long, rings.Integer)):
        nrows, ncols = (ring, nrows)
        ring = rings.ZZ
    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse).zero_matrix()


def block_matrix(sub_matrices, nrows=None, ncols=None, subdivide=True):
    """
    Returns a larger matrix made by concatinating the sub_matrices
    (rows first, then columns). For example, the matrix

        [ A B ]
        [ C D ]

    is made up of submatrices A, B, C, and D.

    INPUT:
        sub_matrices -- matrices (must be of the correct size, or constants)
        nrows         -- (optional) the number of block rows
        ncols         -- (optional) the number of block cols
        subdivide    -- boolean, whether or not to add
                        subdivision information to the matrix

    EXAMPLES:
        sage: A = matrix(QQ, 2, 2, [3,9,6,10])
        sage: block_matrix([A, -A, ~A, 100*A])
        [    3     9|   -3    -9]
        [    6    10|   -6   -10]
        [-----------+-----------]
        [-5/12   3/8|  300   900]
        [  1/4  -1/8|  600  1000]

    One can use constant entries:
        sage: block_matrix([1, A, 0, 1])
        [ 1  0| 3  9]
        [ 0  1| 6 10]
        [-----+-----]
        [ 0  0| 1  0]
        [ 0  0| 0  1]

    One can specify the number of rows or columns (optional for square number of matrices):
        sage: block_matrix([A, -A, ~A, 100*A], ncols=4)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

        sage: block_matrix([A, -A, ~A, 100*A], nrows=1)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

    It handle baserings nicely too:
        sage: R.<x> = ZZ['x']
        sage: block_matrix([1/2, A, 0, x-1])
        [  1/2     0|    3     9]
        [    0   1/2|    6    10]
        [-----------+-----------]
        [    0     0|x - 1     0]
        [    0     0|    0 x - 1]
        sage: block_matrix([1/2, A, 0, x-1]).parent()
        Full MatrixSpace of 4 by 4 dense matrices over Univariate Polynomial Ring in x over Rational Field

    Subdivisions are optional:
        sage: B = matrix(QQ, 2, 3, range(6))
        sage: block_matrix([~A, B, B, ~A], subdivide=False)
        [-5/12   3/8     0     1     2]
        [  1/4  -1/8     3     4     5]
        [    0     1     2 -5/12   3/8]
        [    3     4     5   1/4  -1/8]
    """
    # determine the block dimensions
    n = ZZ(len(sub_matrices))
    if nrows is None:
        if ncols is None:
            if n.is_square():
                nrows = ncols = n.sqrt()
            else:
                raise ValueError, "Must specify rows or cols for non-square block matrix."
        else:
            nrows = int(n/ncols)
    elif ncols is None:
        ncols = int(n/nrows)
    if nrows * ncols != n:
        raise ValueError, "Given number of rows (%s), columns (%s) incompatable with number of submatrices (%s)" % (nrows, ncols, n)

    # empty matrix
    if n == 0:
      ans = Matrix(ZZ,0,0,[])
      if subdivide:
        ans.subdivide([0]*(nrows-1),[0]*(ncols-1))
      return ans

    # determine the sub-block dimensions
    row_heights = [None] * nrows
    col_widths = [None] * ncols
    for i in range(nrows):
        for j in range(0, ncols):
            M = sub_matrices[i*ncols+j]
            if is_Matrix(M):
                if row_heights[i] is None:
                    row_heights[i] = M.nrows()
                if col_widths[j] is None:
                    col_widths[j] = M.ncols()

    if None in row_heights or None in col_widths:
        for i in range(nrows):
            for j in range(0, ncols):
                x = sub_matrices[i*ncols+j]
                if not is_Matrix(x) and x: # must be square matrix
                    if row_heights[i] is None:
                        row_heights[i] = col_widths[j]
                    if col_widths[j] is None:
                        col_widths[j] = row_heights[i]

        if None in row_heights or None in col_widths:
            raise ValueError, "Insufficient information to determine dimensions."

    # determine the base ring
    base = ZZ
    for M in sub_matrices:
        R = M.base_ring() if is_Matrix(M) else M.parent()
        if R is not ZZ:
            base = sage.categories.pushout.pushout(base, R)

    # finally concatinate
    for i in range(nrows):
        for j in range(ncols):
            # coerce
            M = sub_matrices[i*ncols+j]
            if is_Matrix(M):
                if M.base_ring() is not base:
                    M = M.change_ring(base)
            else:
                M = matrix(base, row_heights[i], col_widths[j], M)
            # append
            if j == 0:
                row = M
            else:
                row = row.augment(M)
        if i == 0:
            big = row
        else:
            big = big.stack(row)

    if subdivide:
        big.subdivide(running_total(row_heights[:-1]),
                      running_total(col_widths[:-1]))

    return big

def block_diagonal_matrix(*sub_matrices, **kwds):
    """
    Create a block matrix whose diagonal block entries are given by sub_matrices,
    with zero elsewhere.

    See also \code{block_matrix}.

    EXAMPLES:
        sage: A = matrix(ZZ, 2, [1,2,3,4])
        sage: block_diagonal_matrix(A, A)
        [1 2|0 0]
        [3 4|0 0]
        [---+---]
        [0 0|1 2]
        [0 0|3 4]

    The sub-matrices need not be square:
        sage: B = matrix(QQ, 2, 3, range(6))
        sage: block_diagonal_matrix(~A, B)
        [  -2    1|   0    0    0]
        [ 3/2 -1/2|   0    0    0]
        [---------+--------------]
        [   0    0|   0    1    2]
        [   0    0|   3    4    5]
    """
    if isinstance(sub_matrices, (list, tuple)) and len(sub_matrices) == 1:
        sub_matrices = sub_matrices[0]
    n = len(sub_matrices)
    entries = [ZZ.zero_element()] * n**2
    for i in range(n):
        entries[n*i+i] = sub_matrices[i]
    return block_matrix(entries, **kwds)

def jordan_block(eigenvalue, size, sparse=False):
    r"""
    Form the Jordan block with the specified size associated with the
    eigenvalue.

    INPUT:
        eigenvalue -- eigenvalue for the diagonal entries of the block
        size -- size of the Jordan block
        sparse -- (default False) if True, return a sparse matrix

    EXAMPLE:
        sage: jordan_block(5, 3)
        [5 1 0]
        [0 5 1]
        [0 0 5]

    """
    block = diagonal_matrix([eigenvalue]*size, sparse=sparse)
    for i in xrange(size-1):
        block[i,i+1]=1
    return block

