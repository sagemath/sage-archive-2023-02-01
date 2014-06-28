"""
Matrix Constructor
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
import types
import sage.rings.all as rings
from sage.rings.ring import is_Ring
import sage.matrix.matrix_space as matrix_space
from sage.modules.free_module_element import vector
from sage.structure.element import is_Vector
from sage.structure.sequence import Sequence
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.all import ZZ, QQ
from sage.misc.misc_c import running_total
from matrix import is_Matrix
from copy import copy

import sage.categories.pushout

def matrix_method(func=None, name=None):
    """
    Allows a function to be tab-completed on the global matrix
    constructor object.

    INPUT:

    - ``*function`` -- a single argument. The function that is being
      decorated.

    - ``**kwds`` -- a single optional keyword argument
      ``name=<string>``. The name of the corresponding method in the
      global matrix constructor object. If not given, it is derived
      from the function name.

    EXAMPLES::

        sage: from sage.matrix.constructor import matrix_method
        sage: def foo_matrix(n): return matrix.diagonal(range(n))
        sage: matrix_method(foo_matrix)
        <function foo_matrix at ...>
        sage: matrix.foo(5)
        [0 0 0 0 0]
        [0 1 0 0 0]
        [0 0 2 0 0]
        [0 0 0 3 0]
        [0 0 0 0 4]
        sage: matrix_method(foo_matrix, name='bar')
        <function foo_matrix at ...>
        sage: matrix.bar(3)
        [0 0 0]
        [0 1 0]
        [0 0 2]
    """
    if func is not None:
        if name is None:
            name = func.__name__.replace('matrix', '').strip('_')
        prefix = "    This function is available as %s(...) and matrix.%s(...)." % (
            func.__name__, name)
        func.__doc__ = "%s\n\n%s" % (prefix, func.__doc__)
        setattr(matrix, name, func)
        return func
    else:
        return lambda func: matrix_method(func, name=name)

def _matrix_constructor(*args, **kwds):
    """
    Create a matrix.

    This implements the ``matrix`` constructor::

        sage: matrix([[1,2],[3,4]])
        [1 2]
        [3 4]

    It also contains methods to create special types of matrices, see
    ``matrix.[tab]`` for more options. For example::

        sage: matrix.identity(2)
        [1 0]
        [0 1]

    INPUT:

    The matrix command takes the entries of a matrix, optionally
    preceded by a ring and the dimensions of the matrix, and returns a
    matrix.

    The entries of a matrix can be specified as a flat list of
    elements, a list of lists (i.e., a list of rows), a list of Sage
    vectors, a callable object, or a dictionary having positions as
    keys and matrix entries as values (see the examples). If you pass
    in a callable object, then you must specify the number of rows and
    columns. You can create a matrix of zeros by passing an empty list
    or the integer zero for the entries.  To construct a multiple of
    the identity (`cI`), you can specify square dimensions and pass in
    `c`. Calling matrix() with a Sage object may return something that
    makes sense. Calling matrix() with a NumPy array will convert the
    array to a matrix.

    The ring, number of rows, and number of columns of the matrix can
    be specified by setting the ring, nrows, or ncols parameters or by
    passing them as the first arguments to the function in the order
    ring, nrows, ncols. The ring defaults to ZZ if it is not specified
    or cannot be determined from the entries. If the numbers of rows
    and columns are not specified and cannot be determined, then an
    empty 0x0 matrix is returned.


    -  ``ring`` - the base ring for the entries of the
       matrix.

    -  ``nrows`` - the number of rows in the matrix.

    -  ``ncols`` - the number of columns in the matrix.

    -  ``sparse`` - create a sparse matrix. This defaults
       to True when the entries are given as a dictionary, otherwise
       defaults to False.


    OUTPUT:

    a matrix

    EXAMPLES::

        sage: m=matrix(2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

    ::

        sage: m=matrix(2,3); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    ::

        sage: m=matrix(QQ,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m = matrix(QQ, 3, 3, lambda i, j: i+j); m
        [0 1 2]
        [1 2 3]
        [2 3 4]
        sage: m = matrix(3, lambda i,j: i-j); m
        [ 0 -1 -2]
        [ 1  0 -1]
        [ 2  1  0]

    ::

        sage: matrix(QQ,2,3,lambda x, y: x+y)
        [0 1 2]
        [1 2 3]
        sage: matrix(QQ,3,2,lambda x, y: x+y)
        [0 1]
        [1 2]
        [2 3]

    ::

        sage: v1=vector((1,2,3))
        sage: v2=vector((4,5,6))
        sage: m=matrix([v1,v2]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    ::

        sage: m=matrix(QQ,2,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m=matrix(QQ,2,3,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m=matrix({(0,1): 2, (1,1):2/5}); m; m.parent()
        [  0   2]
        [  0 2/5]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field

    ::

        sage: m=matrix(QQ,2,3,{(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field

    ::

        sage: import numpy
        sage: n=numpy.array([[1,2],[3,4]],float)
        sage: m=matrix(n); m; m.parent()
        [1.0 2.0]
        [3.0 4.0]
        Full MatrixSpace of 2 by 2 dense matrices over Real Double Field

    ::

        sage: v = vector(ZZ, [1, 10, 100])
        sage: m=matrix(v); m; m.parent()
        [  1  10 100]
        Full MatrixSpace of 1 by 3 dense matrices over Integer Ring
        sage: m=matrix(GF(7), v); m; m.parent()
        [1 3 2]
        Full MatrixSpace of 1 by 3 dense matrices over Finite Field of size 7

    ::

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

    ::

        sage: matrix(ZZ, 10, 10, range(100), sparse=True).parent()
        Full MatrixSpace of 10 by 10 sparse matrices over Integer Ring

    ::

        sage: R = PolynomialRing(QQ, 9, 'x')
        sage: A = matrix(R, 3, 3, R.gens()); A
        [x0 x1 x2]
        [x3 x4 x5]
        [x6 x7 x8]
        sage: det(A)
        -x2*x4*x6 + x1*x5*x6 + x2*x3*x7 - x0*x5*x7 - x1*x3*x8 + x0*x4*x8

    TESTS::

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
        TypeError: cannot construct an element of
        Full MatrixSpace of 5 by 1 dense matrices over Rational Field
        from [1, 2, 3, 4, 5, 6]!
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
        sage: matrix(numpy.array([[5]]))
        [5]
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
        sage: matrix(ZZ, [[0] for i in range(10^5)]).is_zero() # see #10158
        True

    AUTHORS:

    - ??: Initial implementation

    - Jason Grout (2008-03): almost a complete rewrite, with bits and
      pieces from the original implementation
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

    if len(args) >= 1 and is_Ring(args[0]):
        # A ring is specified
        if kwds.get('ring', args[0]) != args[0]:
            raise ValueError("Specified rings are not the same")
        else:
            ring = args[0]
            args.pop(0)
    else:
        ring = kwds.get('ring', None)

    if len(args) >= 1:
        # check to see if the number of rows is specified
        try:
            import numpy
            if isinstance(args[0], numpy.ndarray):
                raise TypeError
            nrows = int(args[0])
            args.pop(0)
            if kwds.get('nrows', nrows) != nrows:
                raise ValueError("Number of rows specified in two places and they are not the same")
        except TypeError:
            nrows = kwds.get('nrows', None)
    else:
        nrows = kwds.get('nrows', None)

    if len(args) >= 1:
        # check to see if additionally, the number of columns is specified
        try:
            import numpy
            if isinstance(args[0], numpy.ndarray):
                raise TypeError
            ncols = int(args[0])
            args.pop(0)
            if kwds.get('ncols', ncols) != ncols:
                raise ValueError("Number of columns specified in two places and they are not the same")
        except TypeError:
            ncols = kwds.get('ncols', None)
    else:
        ncols = kwds.get('ncols', None)


    # Now we've taken care of initial ring, nrows, and ncols arguments.
    # We've also taken care of the Sage object case.

    # Now the rest of the arguments are a list of rows, a flat list of
    # entries, a callable, a dict, a numpy array, or a single value.
    if len(args) == 0:
        # If no entries are specified, pass back a zero matrix
        entries = 0
        entry_ring = rings.ZZ
    elif len(args) == 1:
        if isinstance(args[0], (types.FunctionType, types.LambdaType, types.MethodType)):
            if ncols is None and nrows is None:
                raise ValueError("When passing in a callable, the dimensions of the matrix must be specified")
            if ncols is None:
                ncols = nrows
            elif nrows is None:
                nrows = ncols

            f = args[0]
            args[0] = [[f(i,j) for j in range(ncols)] for i in range(nrows)]

        if isinstance(args[0], (list, tuple)):
            if len(args[0]) == 0:
                # no entries are specified, pass back the zero matrix
                entries = 0
                entry_ring = rings.ZZ
            elif isinstance(args[0][0], (list, tuple)) or is_Vector(args[0][0]):
                # Ensure we have a list of lists, each inner list having the same number of elements
                first_len = len(args[0][0])
                if not all( (isinstance(v, (list, tuple)) or is_Vector(v)) and len(v) == first_len for v in args[0]):
                    raise ValueError("List of rows is not valid (rows are wrong types or lengths)")
                # We have a list of rows or vectors
                if nrows is None:
                    nrows = len(args[0])
                elif nrows != len(args[0]):
                    raise ValueError("Number of rows does not match up with specified number.")
                if ncols is None:
                    ncols = len(args[0][0])
                elif ncols != len(args[0][0]):
                    raise ValueError("Number of columns does not match up with specified number.")

                entries = []
                for v in args[0]:
                    entries.extend(v)

            else:
                # We have a flat list; figure out nrows and ncols
                if nrows is None:
                    nrows = 1

                if nrows > 0:
                    if ncols is None:
                        ncols = len(args[0]) // nrows
                    elif ncols != len(args[0]) // nrows:
                        raise ValueError("entries has the wrong length")
                elif len(args[0]) > 0:
                    raise ValueError("entries has the wrong length")

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
                        raise TypeError("cannot convert NumPy matrix to Sage matrix")
                else:
                    raise TypeError("cannot convert NumPy matrix to Sage matrix")

                return m
            elif nrows is not None and ncols is not None:
                # assume that we should just pass the thing into the
                # MatrixSpace constructor and hope for the best
                # This is not documented, but it is doctested
                if ring is None:
                    entry_ring = args[0].parent()
                entries = args[0]
            else:
                raise ValueError("Invalid matrix constructor.  Type matrix? for help")
    else:
        raise ValueError("Invalid matrix constructor.  Type matrix? for help")

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


class MatrixFactory(object):
    """
    The class of the ``matrix`` object.

    See :func:`_matrix_constructor` for the implementation of the call
    interface.  Implemented as a class for nicer tab completion.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getdoc, sage_getsource
        sage: sage_getdoc(matrix)       # used in output of matrix?
        '   Create a matrix.\n\n   This implements the "matrix" ...'
        sage: sage_getsource(matrix)    # used in output of matrix??
        'def _matrix_constructor(*args, **kwds):...'
    """

    __doc__ = _matrix_constructor.__doc__
    __call__ = staticmethod(_matrix_constructor)

    def _sage_src_(self):
        """
        For introspection of the global matrix object.

        TESTS::

            sage: from sage.misc.sageinspect import sage_getsource
            sage: sage_getsource(matrix) == sage_getsource(sage.matrix.constructor._matrix_constructor)
            True
        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(_matrix_constructor)

Matrix = matrix = MatrixFactory()


def prepare(w):
    """
    Given a list w of numbers, find a common ring that they all
    canonically map to, and return the list of images of the elements
    of w in that ring along with the ring.

    This is for internal use by the matrix function.

    INPUT:

    - ``w`` - list

    OUTPUT:

    list, ring

    EXAMPLES::

        sage: sage.matrix.constructor.prepare([-2, Mod(1,7)])
        ([5, 1], Ring of integers modulo 7)

    Notice that the elements must all canonically coerce to a common
    ring (since Sequence is called)::

        sage: sage.matrix.constructor.prepare([2/1, Mod(1,7)])
        Traceback (most recent call last):
        ...
        TypeError: unable to find a common ring for all elements
    """
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
    elif not is_Ring(ring):
        raise TypeError("unable to find a common ring for all elements")
    return entries, ring

def prepare_dict(w):
    """
    Given a dictionary w of numbers, find a common ring that they all
    canonically map to, and return the dictionary of images of the
    elements of w in that ring along with the ring.

    This is for internal use by the matrix function.

    INPUT:

    - ``w`` - dict

    OUTPUT:

    dict, ring

    EXAMPLES::

        sage: sage.matrix.constructor.prepare_dict({(0,1):2, (4,10):Mod(1,7)})
        ({(0, 1): 2, (4, 10): 1}, Ring of integers modulo 7)
    """
    Z = w.items()
    X = [x for _, x in Z]
    entries, ring = prepare(X)
    return dict([(Z[i][0],entries[i]) for i in xrange(len(entries))]), ring

def nrows_from_dict(d):
    """
    Given a dictionary that defines a sparse matrix, return the number
    of rows that matrix should have.

    This is for internal use by the matrix function.

    INPUT:

    - ``d`` - dict

    OUTPUT:

    integer

    EXAMPLES::

        sage: sage.matrix.constructor.nrows_from_dict({})
        0

    Here the answer is 301 not 300, since there is a 0-th row.

    ::
        sage: sage.matrix.constructor.nrows_from_dict({(300,4):10})
        301
    """
    if 0 == len(d):
        return 0
    return max([0] + [ij[0] for ij in d.keys()]) + 1

def ncols_from_dict(d):
    """
    Given a dictionary that defines a sparse matrix, return the number
    of columns that matrix should have.

    This is for internal use by the matrix function.

    INPUT:

    - ``d`` - dict

    OUTPUT:

    integer

    EXAMPLES::

        sage: sage.matrix.constructor.ncols_from_dict({})
        0

    Here the answer is 301 not 300, since there is a 0-th row.

    ::

        sage: sage.matrix.constructor.ncols_from_dict({(4,300):10})
        301
    """
    if 0 == len(d):
        return 0
    return max([0] + [ij[1] for ij in d.keys()]) + 1


@matrix_method
def column_matrix(*args, **kwds):
    r"""
    Constructs a matrix, and then swaps rows for columns and columns for rows.

    .. note::

        Linear algebra in Sage favors rows over columns.  So,
        generally, when creating a matrix, input vectors and lists are
        treated as rows.  This function is a convenience that turns
        around this convention when creating a matrix.  If you are not
        familiar with the usual :class:`matrix <MatrixFactory>`
        constructor, you might want to consider it first.

    INPUT:

    Inputs are almost exactly the same as for the :class:`matrix
    <MatrixFactory>` constructor, which are documented there.  But see
    examples below for how dimensions are handled.

    OUTPUT:

    Output is exactly the transpose of what the :class:`matrix
    <MatrixFactory>` constructor would return.  In other words, the
    ``matrix`` constructor builds a matrix and then this function
    exchanges rows for columns, and columns for rows.

    EXAMPLES:

    The most compelling use of this function is when you have a
    collection of lists or vectors that you would like to become the
    columns of a matrix. In almost any other situation, the
    :class:`matrix <MatrixFactory>` constructor can probably do the
    job just as easily, or easier. ::

        sage: col_1 = [1,2,3]
        sage: col_2 = [4,5,6]
        sage: column_matrix([col_1, col_2])
        [1 4]
        [2 5]
        [3 6]

        sage: v1 = vector(QQ, [10, 20])
        sage: v2 = vector(QQ, [30, 40])
        sage: column_matrix(QQ, [v1, v2])
        [10 30]
        [20 40]

    If you only specify one dimension along with a flat list of entries,
    then it will be the number of columns in the result (which is different
    from the behavior of the ``matrix`` constructor).  ::

        sage: column_matrix(ZZ, 8, range(24))
        [ 0  3  6  9 12 15 18 21]
        [ 1  4  7 10 13 16 19 22]
        [ 2  5  8 11 14 17 20 23]

    And when you specify two dimensions, then they should be number of
    columns first, then the number of rows, which is the reverse of how
    they would be specified for the ``matrix`` constructor. ::

        sage: column_matrix(QQ, 5, 3, range(15))
        [ 0  3  6  9 12]
        [ 1  4  7 10 13]
        [ 2  5  8 11 14]

    And a few unproductive, but illustrative, examples. ::

        sage: A = matrix(ZZ, 3, 4, range(12))
        sage: B = column_matrix(ZZ, 3, 4, range(12))
        sage: A == B.transpose()
        True

        sage: A = matrix(QQ, 7, 12, range(84))
        sage: A == column_matrix(A.columns())
        True

        sage: A=column_matrix(QQ, matrix(ZZ, 3, 2, range(6)) )
        sage: A
        [0 2 4]
        [1 3 5]
        sage: A.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
    """
    return matrix(*args, **kwds).transpose()


@matrix_method
def random_matrix(ring, nrows, ncols=None, algorithm='randomize', *args, **kwds):
    r"""
    Return a random matrix with entries in a specified ring, and possibly with additional properties.

    INPUT:

    -  ``ring`` - base ring for entries of the matrix

    -  ``nrows`` - Integer; number of rows

    -  ``ncols`` - (default: ``None``); number of columns; if ``None``
       defaults to ``nrows``

    -  ``algorithm`` - (default: ``randomize``); determines what properties
       the matrix will have.  See examples below for possible additional
       arguments.

       -  ``randomize`` - randomize the elements of the matrix, possibly
          controlling the density of non-zero entries.

       -  ``echelon_form`` - creates a matrix in echelon form

       -  ``echelonizable`` - creates a matrix that has a predictable
          echelon form

       - ``subspaces`` - creates a matrix whose four subspaces, when
         explored, have reasonably sized, integral valued, entries.

       - ``unimodular`` - creates a matrix of determinant 1.

       - ``diagonalizable`` - creates a diagonalizable matrix whose
         eigenvectors, if computed by hand, will have only integer
         entries.

    -  ``*args, **kwds`` - arguments and keywords to describe additional properties.
       See more detailed documentation below.

    .. warning::

        An upper bound on the absolute value of the entries may be set
        when the ``algorithm`` is ``echelonizable`` or ``unimodular``.
        In these cases it is possible for this constructor to fail with
        a ``ValueError``.   If you *must* have this routine return
        successfully, do not set ``upper_bound``.  This behavior can
        be partially controlled by a ``max_tries`` keyword.

    .. note::

        When constructing matrices with random entries and no additional properties
        (i.e. when ``algorithm='randomize'``), most of the randomness
        is controlled by the ``random_element`` method for elements of the
        base ring of the matrix, so the documentation of that method may be
        relevant or useful.  Also, the default is to not create zero entries,
        unless the ``density`` keyword is set to something strictly less
        than one.

    EXAMPLES:

    Random integer matrices.  With no arguments, the majority of the entries
    are -1 and 1, never zero, and rarely "large." ::

        sage: random_matrix(ZZ, 5, 5)
        [ -8   2   0   0   1]
        [ -1   2   1 -95  -1]
        [ -2 -12   0   0   1]
        [ -1   1  -1  -2  -1]
        [  4  -4  -6   5   0]

    The ``distribution`` keyword  set to ``uniform`` will limit values
    between -2 and 2, and never zero. ::

        sage: random_matrix(ZZ, 5, 5, distribution='uniform')
        [ 1  0 -2  1  1]
        [ 1  0  0  0  2]
        [-1 -2  0  2 -2]
        [-1 -1  1  1  2]
        [ 0 -2 -1  0  0]

    The ``x`` and ``y`` keywords can be used to distribute entries uniformly.
    When both are used ``x`` is the minimum and ``y`` is one greater than the the maximum.
    But still entries are never zero, even if the range contains zero. ::

        sage: random_matrix(ZZ, 4, 8, x=70, y=100)
        [81 82 70 81 78 71 79 94]
        [80 98 89 87 91 94 94 77]
        [86 89 85 92 95 94 72 89]
        [78 80 89 82 94 72 90 92]

        sage: random_matrix(ZZ, 3, 7, x=-5, y=5)
        [-3  3  1 -5  3  1  2]
        [ 3  3  0  3 -5 -2  1]
        [ 0 -2 -2  2 -3 -4 -2]

    If only ``x`` is given, then it is used as the upper bound of a range starting at 0. ::

        sage: random_matrix(ZZ, 5, 5, x=25)
        [20 16  8  3  8]
        [ 8  2  2 14  5]
        [18 18 10 20 11]
        [19 16 17 15  7]
        [ 0 24  3 17 24]

    To allow, and control, zero entries use the ``density`` keyword at a value
    strictly below the default of 1.0, even if distributing entries across an
    interval that does not contain zero already.  Note that for a square matrix it
    is only necessary to set a single dimension. ::

        sage: random_matrix(ZZ, 5, x=-10, y=10, density=0.75)
        [-6  1  0  0  0]
        [ 9  0  0  4  1]
        [-6  0  0 -8  0]
        [ 0  4  0  6  0]
        [ 1 -9  0  0 -8]

        sage: random_matrix(ZZ, 5, x=20, y=30, density=0.75)
        [ 0 28  0 27  0]
        [25 28 20  0  0]
        [ 0 21  0 21  0]
        [ 0 28 22  0  0]
        [ 0  0  0 26 24]

    It is possible to construct sparse matrices, where it may now be advisable
    (but not required) to control the density of nonzero entries. ::

        sage: A=random_matrix(ZZ, 5, 5)
        sage: A.is_sparse()
        False
        sage: A=random_matrix(ZZ, 5, 5, sparse=True)
        sage: A.is_sparse()
        True

        sage: random_matrix(ZZ, 5, 5, density=0.3, sparse=True)
        [ 4  0  0  0 -1]
        [ 0  0  0  0 -7]
        [ 0  0  2  0  0]
        [ 0  0  1  0 -4]
        [ 0  0  0  0  0]

    For algorithm testing you might want to control the number of bits,
    say 10,000 entries, each limited to 16 bits.  ::

        sage: A = random_matrix(ZZ, 100, 100, x=2^16); A
        100 x 100 dense matrix over Integer Ring (use the '.str()' method to see the entries)

    Random rational matrices.  Now ``num_bound`` and ``den_bound`` control the
    generation of random elements, by specifying limits on the absolute value of
    numerators and denominators (respectively).  Entries will be positive and
    negative (map the absolute value function through the entries to get all
    positive values), and zeros are avoided unless the density is set.  If either
    the numerator or denominator bound (or both) is not used, then the values
    default to the distribution for `ZZ` described above that is most frequently
    positive or negative one. ::

        sage: random_matrix(QQ, 2, 8, num_bound=20, den_bound=4)
        [ -1/2     6    13   -12  -2/3  -1/4     5     5]
        [ -9/2   5/3    19  15/2  19/2  20/3 -13/4     0]

        sage: random_matrix(QQ, 4, density = 0.5, sparse=True)
        [    0    71     0  -1/2]
        [    0     0     0     0]
        [31/85     0 -31/2     0]
        [    1  -1/4     0     0]

        sage: A = random_matrix(QQ, 3, 10, num_bound = 99, den_bound = 99)
        sage: positives = map(abs, A.list())
        sage: matrix(QQ, 3, 10, positives)
        [61/18 47/41  1/22   1/2 75/68   6/7     1   1/2 72/41   7/3]
        [33/13   9/2 40/21 45/46 17/22     1 70/79 97/71  7/24  12/5]
        [ 13/8  8/25   1/3 61/14 92/45  4/85  3/38 95/16 82/71   1/5]

        sage: random_matrix(QQ, 4, 10, den_bound = 10)
        [  -1    0  1/8  1/6  2/9 -1/6  1/5 -1/8  1/5 -1/5]
        [ 1/9  1/5   -1  2/9  1/4 -1/7  1/8 -1/9    0    2]
        [ 2/3    2  1/8   -2    0    0   -2    2    0 -1/2]
        [   0    2    1 -2/3    0    0  1/6    0 -1/3 -2/9]

    Random matrices over other rings.  Several classes of matrices have specialized
    ``randomize()`` methods.  You can locate these with the Sage command::

        search_def('randomize')

    The default implementation of :meth:`~sage.matrix.matrix2.randomize` relies
    on the ``random_element()`` method for the base ring.  The ``density`` and
    ``sparse`` keywords behave as described above. ::

        sage: K.<a>=FiniteField(3^2)
        sage: random_matrix(K, 2, 5)
        [      1       a       1 2*a + 1       2]
        [    2*a   a + 2       0       2       1]

        sage: random_matrix(RR, 3, 4, density=0.66)
        [ 0.000000000000000 -0.806696574554030 -0.693915509972359  0.000000000000000]
        [ 0.629781664418083  0.000000000000000 -0.833709843116637  0.000000000000000]
        [ 0.922346867410064  0.000000000000000  0.000000000000000 -0.940316454178921]

        sage: A = random_matrix(ComplexField(32), 3, density=0.8, sparse=True); A
        [                 0.000000000  0.399739209 + 0.909948633*I                  0.000000000]
        [-0.361911424 - 0.455087671*I -0.687810605 + 0.460619713*I  0.625520058 - 0.360952012*I]
        [                 0.000000000                  0.000000000 -0.162196416 - 0.193242896*I]
        sage: A.is_sparse()
        True

    Random matrices in echelon form.  The ``algorithm='echelon_form'`` keyword,
    along with a requested number of non-zero rows (``num_pivots``) will return
    a random matrix in echelon form.  When the base ring is ``QQ`` the result has integer
    entries.  Other exact rings may be also specified. ::

        sage: A=random_matrix(QQ, 4, 8, algorithm='echelon_form', num_pivots=3); A # random
        [ 1 -5  0 -2  0  1  1 -2]
        [ 0  0  1 -5  0 -3 -1  0]
        [ 0  0  0  0  1  2 -2  1]
        [ 0  0  0  0  0  0  0  0]
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (4, 8)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True
        sage: A.rank()
        3
        sage: A==A.rref()
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_rref_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_rref_matrix

    Random matrices with predictable echelon forms.  The ``algorithm='echelonizable'``
    keyword, along with a requested rank (``rank``) and optional size control
    (``upper_bound``) will return a random matrix in echelon form.  When the
    base ring is ``ZZ`` or ``QQ`` the result has integer entries, whose magnitudes
    can be limited by the value of ``upper_bound``, and the echelon form of the
    matrix also has integer entries.  Other exact rings may be also
    specified, but there is no notion of controlling the size.  Square matrices
    of full rank generated by this function always have determinant one, and
    can be constructed with the ``unimodular`` keyword. ::

        sage: A=random_matrix(QQ, 4, 8, algorithm='echelonizable', rank=3, upper_bound=60); A # random
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (4, 8)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True
        sage: A.rank()
        3
        sage: all([abs(x)<60 for x in A.list()])
        True
        sage: A.rref() in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_echelonizable_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_echelonizable_matrix

    Random diagonalizable matrices.  The ``algorithm='diagonalizable'`` keyword,
    along with a requested matrix size (``size``) and optional lists of
    eigenvalues (``eigenvalues``) and the corresponding eigenspace
    dimensions (``dimensions``) will return a random diagonalizable matrix.
    When the eigenvalues and dimensions are not specified the result will have
    randomly generated values for both that fit with the designated size. ::

        sage: A=random_matrix(QQ, 5, algorithm='diagonalizable', eigenvalues=[2,3,-1], dimensions=[1,2,2]); A # random
        sage: all([x in ZZ for x in (A-(2*identity_matrix(5))).rref().list()])
        True
        sage: all([x in ZZ for x in (A-(3*identity_matrix(5))).rref().list()])
        True
        sage: all([x in ZZ for x in (A-(-1*identity_matrix(5))).rref().list()])
        True
        sage: A.jordan_form()
        [ 2| 0| 0| 0| 0]
        [--+--+--+--+--]
        [ 0| 3| 0| 0| 0]
        [--+--+--+--+--]
        [ 0| 0| 3| 0| 0]
        [--+--+--+--+--]
        [ 0| 0| 0|-1| 0]
        [--+--+--+--+--]
        [ 0| 0| 0| 0|-1]

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_diagonalizable_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_diagonalizable_matrix

    Random matrices with predictable subspaces.  The ``algorithm='subspaces'``
    keyword, along with an optional rank (``rank``) will return
    a matrix whose natural basis vectors for its four fundamental subspaces, if computed as
    described in the documentation of the :func:`~sage.matrix.constructor.random_subspaces_matrix`
    contain only integer entries.  If ``rank``, is not set, the
    rank of the matrix will be generated randomly. ::

        sage: B=random_matrix(QQ, 5, 6, algorithm='subspaces', rank=3); B #random
        sage: B_expanded=B.augment(identity_matrix(5)).rref()
        sage: (B.nrows(), B.ncols())
        (5, 6)
        sage: all([x in ZZ for x in B_expanded.list()])
        True
        sage: C=B_expanded.submatrix(0,0,B.nrows()-B.nullity(),B.ncols())
        sage: L=B_expanded.submatrix(B.nrows()-B.nullity(),B.ncols())
        sage: B.right_kernel()==C.right_kernel()
        True
        sage: B.row_space()==C.row_space()
        True
        sage: B.column_space()==L.right_kernel()
        True
        sage: B.left_kernel()==L.row_space()
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_subspaces_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_subspaces_matrix

    Random unimodular matrices.  The ``algorithm='unimodular'``
    keyword, along with an optional entry size control (``upper_bound``)
    will return a matrix of determinant 1. When the base ring is ``ZZ``
    or ``QQ`` the result has integer entries, whose magnitudes
    can be limited by the value of ``upper_bound``. ::

        sage: C=random_matrix(QQ, 5, algorithm='unimodular', upper_bound=70); C # random
        sage: det(C)
        1
        sage: C.base_ring()
        Rational Field
        sage: (C.nrows(), C.ncols())
        (5, 5)
        sage: all([abs(x)<70 for x in C.list()])
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_unimodular_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_unimodular_matrix

    TESTS:

    We return an error for a bogus value of ``algorithm``::

        sage: random_matrix(ZZ, 5, algorithm = 'bogus')
        Traceback (most recent call last):
        ...
        ValueError: random matrix algorithm "bogus" is not recognized

    AUTHOR:

    - William Stein (2007-02-06)

    - Rob Beezer (2010-08-25) Documentation, code to allow additional types of output
    """
    if ncols is None:
        ncols = nrows
    sparse = kwds.pop('sparse', False)
    # Construct the parent of the desired matrix
    parent = matrix_space.MatrixSpace(ring, nrows, ncols, sparse=sparse)
    if algorithm == 'randomize':
        density = kwds.pop('density', None)
        # zero matrix is immutable, copy is mutable
        A = copy(parent.zero_matrix())
        if density is None:
            A.randomize(density=float(1), nonzero=False, *args, **kwds)
        else:
            A.randomize(density=density, nonzero=True, *args, **kwds)
        return A
    elif algorithm == 'echelon_form':
        return random_rref_matrix(parent, *args, **kwds)
    elif algorithm == 'echelonizable':
        return random_echelonizable_matrix(parent, *args, **kwds)
    elif algorithm == 'diagonalizable':
        return random_diagonalizable_matrix(parent, *args, **kwds)
    elif algorithm == 'subspaces':
        return random_subspaces_matrix(parent, *args, **kwds)
    elif algorithm == 'unimodular':
        return random_unimodular_matrix(parent, *args, **kwds)
    else:
        raise ValueError('random matrix algorithm "%s" is not recognized' % algorithm)


@matrix_method
def diagonal_matrix(arg0=None, arg1=None, arg2=None, sparse=True):
    r"""
    Return a square matrix with specified diagonal entries, and zeros elsewhere.

    FORMATS:

      1. diagonal_matrix(entries)

      2. diagonal_matrix(nrows, entries)

      3. diagonal_matrix(ring, entries)

      4. diagonal_matrix(ring, nrows, entries)

    INPUT:

    - ``entries`` - the values to place along the diagonal
      of the returned matrix.  This may be a flat list, a
      flat tuple, a vector or free module element, or
      a one-dimensional NumPy array.

    - ``nrows`` - the size of the returned matrix, which
      will have an equal number of columns

    - ``ring`` - the ring containing the entries of the
      diagonal entries.  This may not be specified in
      combination with a NumPy array.

    - ``sparse`` - default: ``True`` - whether or not
      the result has a sparse implementation.

    OUTPUT:

    A square matrix over the given ``ring`` with a size
    given by ``nrows``.  If the ring is not given it
    is inferred from the given entries.  The values on
    the diagonal of the returned matrix come from ``entries``.
    If the number of entries is not enough to fill the whole
    diagonal, it is padded with zeros.

    EXAMPLES:

    We first demonstrate each of the input formats with various
    different ways to specify the entries.

    Format 1: a flat list of entries.  ::

        sage: A = diagonal_matrix([2, 1.3, 5]); A
        [ 2.00000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000  1.30000000000000 0.000000000000000]
        [0.000000000000000 0.000000000000000  5.00000000000000]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Real Field with 53 bits of precision

    Format 2: size specified, a tuple with initial entries. Note that a short list of entries
    is effectively padded with zeros.  ::

        sage: A = diagonal_matrix(3, (4, 5)); A
        [4 0 0]
        [0 5 0]
        [0 0 0]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring

    Format 3: ring specified, a vector of entries. ::

        sage: A = diagonal_matrix(QQ, vector(ZZ, [1,2,3])); A
        [1 0 0]
        [0 2 0]
        [0 0 3]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Rational Field

    Format 4: ring, size and list of entries. ::

        sage: A = diagonal_matrix(FiniteField(3), 3, [2, 16]); A
        [2 0 0]
        [0 1 0]
        [0 0 0]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Finite Field of size 3

    NumPy arrays may be used as input. ::

        sage: import numpy
        sage: entries = numpy.array([1.2, 5.6]); entries
        array([ 1.2,  5.6])
        sage: A = diagonal_matrix(3, entries); A
        [1.2 0.0 0.0]
        [0.0 5.6 0.0]
        [0.0 0.0 0.0]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Real Double Field

        sage: j = numpy.complex(0,1)
        sage: entries = numpy.array([2.0+j, 8.1, 3.4+2.6*j]); entries
        array([ 2.0+1.j ,  8.1+0.j ,  3.4+2.6j])
        sage: A = diagonal_matrix(entries); A
        [2.0 + 1.0*I         0.0         0.0]
        [        0.0         8.1         0.0]
        [        0.0         0.0 3.4 + 2.6*I]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Complex Double Field

        sage: entries = numpy.array([4, 5, 6])
        sage: A = diagonal_matrix(entries); A
        [4 0 0]
        [0 5 0]
        [0 0 6]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring

        sage: entries = numpy.array([4.1, 5.2, 6.3])
        sage: A = diagonal_matrix(ZZ, entries); A
        Traceback (most recent call last):
        ...
        TypeError: Cannot convert non-integral float to integer

    By default returned matrices have a sparse implementation.  This can be changed
    when using any of the formats.  ::

        sage: A = diagonal_matrix([1,2,3], sparse=False)
        sage: A.parent()
        Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

    An empty list and no ring specified defaults to the integers. ::

        sage: A = diagonal_matrix([])
        sage: A.parent()
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring

    Giving the entries improperly may first complain about not having a length.  ::

        sage: diagonal_matrix(QQ, 5, 10)
        Traceback (most recent call last):
        ...
        TypeError: unable to determine number of entries for diagonal matrix construction

    Giving too many entries will raise an error. ::

        sage: diagonal_matrix(QQ, 3, [1,2,3,4])
        Traceback (most recent call last):
        ...
        ValueError: number of diagonal matrix entries (4) exceeds the requested matrix size (3)

    A negative size sometimes causes the error that there are too many elements. ::

        sage: diagonal_matrix(-2, [2])
        Traceback (most recent call last):
        ...
        ValueError: number of diagonal matrix entries (1) exceeds the requested matrix size (-2)

    Types for the entries are limited, even though they may have a length.  ::

        sage: diagonal_matrix(x^2)
        Traceback (most recent call last):
        ...
        TypeError: diagonal matrix entries are not a supported type (list, tuple, vector, or NumPy array)

    AUTHOR:

        - Rob Beezer (2011-01-11): total rewrite
    """
    # Roll arguments leftward
    #
    # Leads with a ring?
    # Formats 3, 4, else remains None
    ring = None
    if is_Ring(arg0):
        ring = arg0
        arg0 = arg1
        arg1 = arg2
    # Size of matrix specified?
    # Formats 2, 4
    nrows = None
    if isinstance(arg0, (int, long, rings.Integer)):
        nrows = arg0
        arg0 = arg1
    # Object holding entries
    # Formats 1, 2, 3, 4
    entries = arg0

    # Reconcile matrix size and number of entries
    try:
        nentries = len(entries)
    except TypeError:
        raise TypeError('unable to determine number of entries for diagonal matrix construction')
    # sometimes catches a negative size
    if not nrows is None and nentries > nrows:
        raise ValueError('number of diagonal matrix entries (%s) exceeds the requested matrix size (%s)' % (nentries, nrows))
    if nrows is None:
        nrows = nentries

    # provide a default ring for an empty list
    if len(entries) == 0 and ring is None:
      ring = rings.ZZ

    # Sanity check on entries (partially, e.g. a list of lists will survive this check)
    from numpy import ndarray
    if not any([isinstance(entries, (list, tuple)), isinstance(entries, ndarray), is_Vector(entries)]):
        raise TypeError('diagonal matrix entries are not a supported type (list, tuple, vector, or NumPy array)')

    # Convert entries to a list v over a common ring
    from sage.modules.free_module_element import prepare
    v, ring = prepare(entries, ring)

    # Create a "diagonal" dictionary for matrix constructor
    # If nentries < nrows, diagonal is effectively padded with zeros at end
    w = {}
    for i in range(len(v)):
        w[(i,i)] = v[i]

    # Ship ring, matrix size, dictionary to matrix constructor
    if ring is None:
        return matrix(nrows, nrows, w, sparse=sparse)
    else:
        return matrix(ring, nrows, nrows, w, sparse=sparse)


@matrix_method
def identity_matrix(ring, n=0, sparse=False):
    r"""
    Return the `n \times n` identity matrix over the given
    ring.

    The default ring is the integers.

    EXAMPLES::

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
        sage: M.is_mutable()
        True
        sage: M = identity_matrix(3, sparse=True); M
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: M.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring
        sage: M.is_mutable()
        True
    """
    if isinstance(ring, (int, long, rings.Integer)):
        n = ring
        ring = rings.ZZ
    return matrix_space.MatrixSpace(ring, n, n, sparse)(1)


@matrix_method
def zero_matrix(ring, nrows, ncols=None, sparse=False):
    r"""
    Return the `nrows \times ncols` zero matrix over the given
    ring.

    The default ring is the integers.

    EXAMPLES::

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
        sage: M.is_mutable()
        True
        sage: M = zero_matrix(3, 1, sparse=True); M
        [0]
        [0]
        [0]
        sage: M.parent()
        Full MatrixSpace of 3 by 1 sparse matrices over Integer Ring
        sage: M.is_mutable()
        True
    """
    if isinstance(ring, (int, long, rings.Integer)):
        nrows, ncols = (ring, nrows)
        ring = rings.ZZ
    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse)(0)


@matrix_method
def ones_matrix(ring, nrows=None, ncols=None, sparse=False):
    r"""
    Return a matrix with all entries equal to 1.

    CALL FORMATS:

    In each case, the optional keyword ``sparse`` can be used.

      1. ones_matrix(ring, nrows, ncols)
      2. ones_matrix(ring, nrows)
      3. ones_matrix(nrows, ncols)
      4. ones_matrix(nrows)

    INPUT:

    - ``ring`` - default: ``ZZ`` - base ring for the matrix.
    - ``nrows`` - number of rows in the matrix.
    - ``ncols`` - number of columns in the matrix.
      If omitted, defaults to the number of rows, producing a square matrix.
    - ``sparse`` - default: ``False`` - if ``True`` creates a sparse representation.

    OUTPUT:

    A matrix of size ``nrows`` by ``ncols`` over the ``ring`` with every
    entry equal to 1.  While the result is far from sparse, you may wish
    to choose a sparse representation when mixing this matrix with
    other sparse matrices.

    EXAMPLES:

    A call specifying the ring and the size.  ::

        sage: M= ones_matrix(QQ, 2, 5); M
        [1 1 1 1 1]
        [1 1 1 1 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 5 dense matrices over Rational Field

    Without specifying the number of columns, the result is square. ::

        sage: M = ones_matrix(RR, 2); M
        [1.00000000000000 1.00000000000000]
        [1.00000000000000 1.00000000000000]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Real Field with 53 bits of precision

    The ring defaults to the integers if not given. ::

        sage: M = ones_matrix(2, 3); M
        [1 1 1]
        [1 1 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    A lone integer input produces a square matrix over the integers. ::

        sage: M = ones_matrix(3); M
        [1 1 1]
        [1 1 1]
        [1 1 1]
        sage: M.parent()
        Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

    The result can have a sparse implementation. ::

        sage: M = ones_matrix(3, 1, sparse=True); M
        [1]
        [1]
        [1]
        sage: M.parent()
        Full MatrixSpace of 3 by 1 sparse matrices over Integer Ring

    Giving just a ring will yield an error. ::

        sage: ones_matrix(CC)
        Traceback (most recent call last):
        ...
        ValueError: constructing an all ones matrix requires at least one dimension
    """
    if isinstance(ring, (int, long, rings.Integer)):
        nrows, ncols = (ring, nrows)
        ring = rings.ZZ
    if nrows is None:
        raise ValueError("constructing an all ones matrix requires at least one dimension")
    if ncols is None:
        nents = nrows**2
    else:
        nents = nrows*ncols
    one = ring(1)
    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse).matrix([one]*nents)


@matrix_method
def elementary_matrix(arg0, arg1=None, **kwds):
    r"""
    Creates a square matrix that corresponds to a row operation or a column operation.

    FORMATS:

    In each case, ``R`` is the base ring, and is optional. ``n`` is the size
    of the square matrix created.  Any call may include the ``sparse`` keyword
    to determine the representation used.  The default is ``False`` which
    leads to a dense representation.  We describe the matrices by their
    associated row operation, see the output description for more.

    - ``elementary_matrix(R, n, row1=i, row2=j)``

      The matrix which swaps rows ``i`` and ``j``.

    - ``elementary_matrix(R, n, row1=i, scale=s)``

      The matrix which multiplies row ``i`` by ``s``.

    - ``elementary_matrix(R, n, row1=i, row2=j, scale=s)``

      The matrix which multiplies row ``j`` by ``s``
      and adds it to row ``i``.

    Elementary matrices representing column operations are created
    in an entirely analogous way, replacing ``row1`` by ``col1`` and
    replacing ``row2`` by ``col2``.

    Specifying the ring for entries of the matrix is optional.  If it
    is not given, and a scale parameter is provided, then a ring containing
    the value of ``scale`` will be used.  Otherwise, the ring defaults
    to the integers.

    OUTPUT:

    An elementary matrix is a square matrix that is very close to being
    an identity matrix.  If ``E`` is an elementary matrix and ``A`` is any
    matrix with the same number of rows, then ``E*A`` is the result of
    applying a row operation to ``A``.  This is how the three types
    created by this function are described.  Similarly, an elementary matrix
    can be associated with a column operation, so if ``E`` has the same number
    of columns as ``A`` then ``A*E`` is the result of performing a column
    operation on ``A``.

    An elementary matrix representing a row operation is created if ``row1``
    is specified, while an elementary matrix representing a column operation
    is created if ``col1`` is specified.

    EXAMPLES:

    Over the integers, creating row operations. Recall that row
    and column numbering begins at zero.  ::

        sage: A = matrix(ZZ, 4, 10, range(40)); A
        [ 0  1  2  3  4  5  6  7  8  9]
        [10 11 12 13 14 15 16 17 18 19]
        [20 21 22 23 24 25 26 27 28 29]
        [30 31 32 33 34 35 36 37 38 39]

        sage: E = elementary_matrix(4, row1=1, row2=3); E
        [1 0 0 0]
        [0 0 0 1]
        [0 0 1 0]
        [0 1 0 0]
        sage: E*A
        [ 0  1  2  3  4  5  6  7  8  9]
        [30 31 32 33 34 35 36 37 38 39]
        [20 21 22 23 24 25 26 27 28 29]
        [10 11 12 13 14 15 16 17 18 19]

        sage: E = elementary_matrix(4, row1=2, scale=10); E
        [ 1  0  0  0]
        [ 0  1  0  0]
        [ 0  0 10  0]
        [ 0  0  0  1]
        sage: E*A
        [  0   1   2   3   4   5   6   7   8   9]
        [ 10  11  12  13  14  15  16  17  18  19]
        [200 210 220 230 240 250 260 270 280 290]
        [ 30  31  32  33  34  35  36  37  38  39]

        sage: E = elementary_matrix(4, row1=2, row2=1, scale=10); E
        [ 1  0  0  0]
        [ 0  1  0  0]
        [ 0 10  1  0]
        [ 0  0  0  1]
        sage: E*A
        [  0   1   2   3   4   5   6   7   8   9]
        [ 10  11  12  13  14  15  16  17  18  19]
        [120 131 142 153 164 175 186 197 208 219]
        [ 30  31  32  33  34  35  36  37  38  39]

    Over the rationals, now as column operations. Recall that row
    and column numbering begins at zero.  Checks now have the
    elementary matrix on the right.  ::

        sage: A = matrix(QQ, 5, 4, range(20)); A
        [ 0  1  2  3]
        [ 4  5  6  7]
        [ 8  9 10 11]
        [12 13 14 15]
        [16 17 18 19]

        sage: E = elementary_matrix(QQ, 4, col1=1, col2=3); E
        [1 0 0 0]
        [0 0 0 1]
        [0 0 1 0]
        [0 1 0 0]
        sage: A*E
        [ 0  3  2  1]
        [ 4  7  6  5]
        [ 8 11 10  9]
        [12 15 14 13]
        [16 19 18 17]

        sage: E = elementary_matrix(QQ, 4, col1=2, scale=1/2); E
        [  1   0   0   0]
        [  0   1   0   0]
        [  0   0 1/2   0]
        [  0   0   0   1]
        sage: A*E
        [ 0  1  1  3]
        [ 4  5  3  7]
        [ 8  9  5 11]
        [12 13  7 15]
        [16 17  9 19]

        sage: E = elementary_matrix(QQ, 4, col1=2, col2=1, scale=10); E
        [ 1  0  0  0]
        [ 0  1 10  0]
        [ 0  0  1  0]
        [ 0  0  0  1]
        sage: A*E
        [  0   1  12   3]
        [  4   5  56   7]
        [  8   9 100  11]
        [ 12  13 144  15]
        [ 16  17 188  19]

    An elementary matrix is always nonsingular.  Then repeated row
    operations can be represented by products of elementary matrices,
    and this product is again nonsingular.  If row operations are to
    preserve fundamental properties of a matrix (like rank), we do not
    allow scaling a row by zero.  Similarly, the corresponding elementary
    matrix is not constructed.  Also, we do not allow adding a multiple
    of a row to itself, since this could also lead to a new zero row.  ::

        sage: A = matrix(QQ, 4, 10, range(40)); A
        [ 0  1  2  3  4  5  6  7  8  9]
        [10 11 12 13 14 15 16 17 18 19]
        [20 21 22 23 24 25 26 27 28 29]
        [30 31 32 33 34 35 36 37 38 39]

        sage: E1 = elementary_matrix(QQ, 4, row1=0, row2=1)
        sage: E2 = elementary_matrix(QQ, 4, row1=3, row2=0, scale=100)
        sage: E = E2*E1
        sage: E.is_singular()
        False
        sage: E*A
        [  10   11   12   13   14   15   16   17   18   19]
        [   0    1    2    3    4    5    6    7    8    9]
        [  20   21   22   23   24   25   26   27   28   29]
        [1030 1131 1232 1333 1434 1535 1636 1737 1838 1939]

        sage: E3 = elementary_matrix(QQ, 4, row1=3, scale=0)
        Traceback (most recent call last):
        ...
        ValueError: scale parameter of row of elementary matrix must be non-zero

        sage: E4 = elementary_matrix(QQ, 4, row1=3, row2=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: cannot add a multiple of a row to itself

    If the ring is not specified, and a scale parameter is given, the
    base ring for the matrix is chosen to contain the scale parameter.
    Otherwise, if no ring is given, the default is the integers. ::

        sage: E = elementary_matrix(4, row1=1, row2=3)
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Integer Ring

        sage: E = elementary_matrix(4, row1=1, scale=4/3)
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Rational Field

        sage: E = elementary_matrix(4, row1=1, scale=I)
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Symbolic Ring

        sage: E = elementary_matrix(4, row1=1, scale=CDF(I))
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Complex Double Field

        sage: E = elementary_matrix(4, row1=1, scale=QQbar(I))
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Algebraic Field

    Returned matrices have a dense implementation by default,
    but a sparse implementation may be requested.  ::

        sage: E = elementary_matrix(4, row1=0, row2=1)
        sage: E.is_dense()
        True

        sage: E = elementary_matrix(4, row1=0, row2=1, sparse=True)
        sage: E.is_sparse()
        True

    And the ridiculously small cases.  The zero-row matrix cannot be built
    since then there are no rows to manipulate. ::

        sage: elementary_matrix(QQ, 1, row1=0, row2=0)
        [1]
        sage: elementary_matrix(QQ, 0, row1=0, row2=0)
        Traceback (most recent call last):
        ...
        ValueError: size of elementary matrix must be 1 or greater, not 0

    TESTS::

        sage: E = elementary_matrix('junk', 5, row1=3, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        TypeError: optional first parameter must be a ring, not junk

        sage: E = elementary_matrix(5, row1=3, scale='junk')
        Traceback (most recent call last):
        ...
        TypeError: scale must be an element of some ring, not junk

        sage: E = elementary_matrix(ZZ, 5, row1=3, col2=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: received an unexpected keyword: col2=3

        sage: E = elementary_matrix(QQ, row1=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: size of elementary matrix must be given

        sage: E = elementary_matrix(ZZ, 4/3, row1=3, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        TypeError: size of elementary matrix must be an integer, not 4/3

        sage: E = elementary_matrix(ZZ, -3, row1=3, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: size of elementary matrix must be 1 or greater, not -3

        sage: E = elementary_matrix(ZZ, 5, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: row1 or col1 must be specified

        sage: E = elementary_matrix(ZZ, 5, row1=3, col1=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: cannot specify both row1 and col1

        sage: E = elementary_matrix(ZZ, 5, row1=4/3, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        TypeError: row of elementary matrix must be an integer, not 4/3

        sage: E = elementary_matrix(ZZ, 5, col1=5, col2=1, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: column of elementary matrix must be positive and smaller than 5, not 5

        sage: E = elementary_matrix(ZZ, 5, col1=3, col2=4/3, scale=12)
        Traceback (most recent call last):
        ...
        TypeError: column of elementary matrix must be an integer, not 4/3

        sage: E = elementary_matrix(ZZ, 5, row1=3, row2=-1, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: row of elementary matrix must be positive and smaller than 5, not -1

        sage: E = elementary_matrix(ZZ, 5, row1=3, row2=1, scale=4/3)
        Traceback (most recent call last):
        ...
        TypeError: scale parameter of elementary matrix must an element of Integer Ring, not 4/3

        sage: E = elementary_matrix(ZZ, 5, row1=3)
        Traceback (most recent call last):
        ...
        ValueError: insufficient parameters provided to construct elementary matrix

        sage: E = elementary_matrix(ZZ, 5, row1=3, row2=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: cannot add a multiple of a row to itself

        sage: E = elementary_matrix(ZZ, 5, col1=3, scale=0)
        Traceback (most recent call last):
        ...
        ValueError: scale parameter of column of elementary matrix must be non-zero

    AUTHOR:

    - Rob Beezer (2011-03-04)
    """
    import sage.structure.element
    # determine ring and matrix size
    if not arg1 is None and not is_Ring(arg0):
        raise TypeError('optional first parameter must be a ring, not {0}'.format(arg0))
    scale = kwds.pop('scale', None)
    if is_Ring(arg0):
        R = arg0
        arg0 = arg1
    elif scale is not None:
        if not sage.structure.element.is_RingElement(scale):
            raise TypeError('scale must be an element of some ring, not {0}'.format(scale))
        R = scale.parent()
    else:
        R = rings.ZZ
    if arg0 is None:
        raise ValueError('size of elementary matrix must be given')
    try:
        n = rings.Integer(arg0)
    except TypeError:
        raise TypeError('size of elementary matrix must be an integer, not {0}'.format(arg0))
    if n <= 0:
        raise ValueError('size of elementary matrix must be 1 or greater, not {0}'.format(n))
    # row operations or column operations?
    # column operation matrix will be transpose of a row operation matrix
    row1 = kwds.pop('row1', None)
    col1 = kwds.pop('col1', None)
    if row1 is None and col1 is None:
        raise ValueError('row1 or col1 must be specified')
    if not row1 is None and not col1 is None:
        raise ValueError('cannot specify both row1 and col1')
    rowop = not row1 is None
    if rowop:
        opstring = "row"
        row2 = kwds.pop('row2', None)
    else:
        opstring = "column"
        row1 = col1
        row2 = kwds.pop('col2', None)
    sparse = kwds.pop('sparse', False)
    if kwds:
        extra = kwds.popitem()
        raise ValueError('received an unexpected keyword: {0}={1}'.format(extra[0], extra[1]))

    # analyze parameters to determine matrix type
    try:
        row1 = rings.Integer(row1)
    except TypeError:
        raise TypeError('{0} of elementary matrix must be an integer, not {1}'.format(opstring, row1))
    if row1 < 0 or row1 >= n :
        raise ValueError('{0} of elementary matrix must be positive and smaller than {1}, not {2}'.format(opstring, n, row1))
    if not row2 is None:
        try:
            row2 = rings.Integer(row2)
        except TypeError:
            raise TypeError('{0} of elementary matrix must be an integer, not {1}'.format(opstring, row2))
        if row2 < 0 or row2 >= n :
            raise ValueError('{0} of elementary matrix must be positive and smaller than {1}, not {2}'.format(opstring, n, row2))
    if not scale is None:
        try:
            scale = R(scale)
        except Exception:
            raise TypeError('scale parameter of elementary matrix must an element of {0}, not {1}'.format(R, scale))

    # determine type of matrix and adjust an identity matrix
    # return row operation matrix or the transpose as a column operation matrix
    elem = identity_matrix(R, n, sparse=sparse)
    if row2 is None and scale is None:
        raise ValueError('insufficient parameters provided to construct elementary matrix')
    elif not row2 is None and not scale is None:
        if row1 == row2:
            raise ValueError('cannot add a multiple of a {0} to itself'.format(opstring))
        elem[row1, row2] = scale
    elif not row2 is None and scale is None:
        elem[row1, row1] = 0
        elem[row2, row2] = 0
        elem[row1, row2] = 1
        elem[row2, row1] = 1
    elif row2 is None and not scale is None:
        if scale == 0:
            raise ValueError('scale parameter of {0} of elementary matrix must be non-zero'.format(opstring))
        elem[row1, row1] = scale
    if rowop:
        return elem
    else:
        return elem.transpose()

def _determine_block_matrix_grid(sub_matrices):
    r"""
    For internal use. This tries to determine the dimensions
    of rows/columns when assembling the matrices in sub_matrices in a
    rectangular grid. It returns a pair of lists containing
    respectively the sizes of rows and columns.

    sub_matrices must be a list of lists of matrices. All sublists
    are expected to be the same size.

    Non-zero scalars are considered to be square matrices of any size,
    and zeroes are considered to be zero matrices of any size.

    A ValueError is raised if there is insufficient or
    conflicting information.

    TESTS::

        sage: from sage.matrix.constructor import _determine_block_matrix_grid
        sage: A = matrix(QQ, 2, 2, [3,9,6,10])
        sage: _determine_block_matrix_grid([[A, A], [A, A]])
        ([2, 2], [2, 2])

        sage: B = matrix(QQ, 1, 1, [ 1 ] )
        sage: C = matrix(QQ, 2, 2, [ 2, 3, 4, 5 ] )
        sage: _determine_block_matrix_grid([[B, 0], [0, C]])
        ([1, 2], [1, 2])
    """

    nrows = len(sub_matrices)
    if nrows == 0:
        return ([], [])
    ncols = len(sub_matrices[0])

    if ncols == 0:
        return ([0] * nrows, [])

    row_heights = [None] * nrows
    col_widths = [None] * ncols

    changing = True
    while changing:
        changing = False
        for i in range(nrows):
            for j in range(ncols):
                M = sub_matrices[i][j]
                sub_width = None
                sub_height = None
                if is_Matrix(M):
                    sub_width = M.ncols()
                    sub_height = M.nrows()
                elif M: # non-zero scalar is interpreted as a square matrix
                    if row_heights[i] is None:
                        sub_width = col_widths[j]
                    else:
                        sub_width = row_heights[i]
                    sub_height = sub_width
                if sub_width is not None:
                    if col_widths[j] is None:
                        changing = True
                        col_widths[j] = sub_width
                    elif col_widths[j] != sub_width:
                        raise ValueError("incompatible submatrix widths")
                if sub_height is not None:
                    if row_heights[i] is None:
                        changing = True
                        row_heights[i] = sub_height
                    elif row_heights[i] != sub_height:
                        raise ValueError("incompatible submatrix heights")

    if None in row_heights or None in col_widths:
        if None in row_heights or None in col_widths:
            raise ValueError("insufficient information to determine dimensions.")

    return (row_heights, col_widths)

def _determine_block_matrix_rows(sub_matrices):
    """
    For internal use. This tests if the matrices in sub_matrices
    fit in a rectangular matrix when assembled a row at a time.

    sub_matrices must be a list of lists of matrices.

    It returns a pair (row_heights, zero_widths, width) where
    row_heights is the list of row heights, zero_widths is the
    total width filled up by zero matrices per row, and width
    is the total width of the resulting matrix.

    Non-zero scalars are considered to be square matrices of any size,
    and zeroes are considered to be zero matrices of any size.

    A ValueError is raised if there is insufficient or
    conflicting information.

    TESTS::

        sage: from sage.matrix.constructor import _determine_block_matrix_rows
        sage: A = Matrix(ZZ, 1, 4, [1, 2, 3, 4])
        sage: _determine_block_matrix_rows([ [1, 1], [ A ] ])
        ([2, 1], [0, 0], 4)

        sage: B = Matrix(ZZ, 2, 2, [1, 2, 3, 4])
        sage: _determine_block_matrix_rows([ [B, B], [B, 1] ])
        ([2, 2], [0, 0], 4)
    """

    total_height = 0
    total_width = None

    row_heights = [ None ] * len(sub_matrices)
    zero_widths = [ 0 ] * len(sub_matrices)

    # We first do a pass to see if we can determine the width
    unknowns = False
    for i in range(len(sub_matrices)):
        R = sub_matrices[i]
        height = None
        # We first do a pass to see if we can determine the height
        # of this row
        found_zeroes = False
        for M in R:
            if is_Matrix(M):
                if height is None:
                    height = M.nrows()
                elif height != M.nrows():
                    raise ValueError("incompatible submatrix heights")
            elif not M:
                found_zeroes = True
        if len(R) == 0:
            height = 0

        # If we have a height, then we know the dimensions of any
        # non-zero scalars, and can maybe compute the width
        if height is not None and not found_zeroes:
            width = 0
            for M in R:
                if is_Matrix(M):
                    width += M.ncols()
                else:
                    # non-zero scalar
                    width += height
            if total_width is None:
                total_width = width
            elif total_width != width:
                raise ValueError("incompatible submatrix widths")
            row_heights[i] = height
        else:
            # We don't set height here even if we know it,
            # to signal this row hasn't been fit yet.
            unknowns = True

    if total_width is None:
        raise ValueError("insufficient information to determine submatrix widths")

    if unknowns:
        # Do a second pass and see if the remaining rows can be
        # determined now that we know the width of the matrix.
        for i in range(len(sub_matrices)):
            if row_heights[i] is not None:
                continue
            R = sub_matrices[i]
            zero_state = 0
            # 0: no zeroes found
            # 1: consecutive zeroes found
            # 2: consecutive zeroes followed by non-zero found
            # 3: non-consecutive zeroes found
            scalars = 0
            width = 0
            height = None
            for j in range(len(R)):
                M = R[j]
                if is_Matrix(M):
                    height = M.nrows()
                    width += M.ncols()
                    if zero_state == 1:
                        zero_state = 2
                elif not M:
                    if zero_state == 0:
                        zero_state = 1
                    elif zero_state == 2:
                        zero_state = 3
                else:
                     scalars += 1

            remaining_width = total_width - width
            # This remaining width has to be split over the
            # zeroes and (non-zero) scalars

            if height is not None:
                remaining_width -= scalars * height
                if remaining_width < 0:
                    raise ValueError("incompatible submatrix widths")
                if remaining_width > 0 and zero_state == 3:
                    raise ValueError("insufficient information to determine submatrix widths")
                if remaining_width > 0 and zero_state == 0:
                    raise ValueError("incompatible submatrix widths")
                # otherwise, things fit
                row_heights[i] = height
                zero_widths[i] = remaining_width
            elif zero_state != 0:
                # if we don't know the height, and there are zeroes,
                # we can't determine the height
                raise ValueError("insufficient information to determine submatrix heights")
            elif total_width % len(R) != 0:
                raise ValueError("incompatible submatrix widths")
            else:
                height = int(total_width / len(R))
                row_heights[i] = height

    # If we got this far, then everything fits
    return (row_heights, zero_widths, total_width)


@matrix_method
def block_matrix(*args, **kwds):
    r"""
    Returns a larger matrix made by concatenating submatrices
    (rows first, then columns). For example, the matrix

    ::

        [ A B ]
        [ C D ]

    is made up of submatrices A, B, C, and D.

    INPUT:

    The block_matrix command takes a list of submatrices to add
    as blocks, optionally preceded by a ring and the number of block rows
    and block columns, and returns a matrix.

    The submatrices can be specified as a list of matrices (using
    ``nrows`` and ``ncols`` to determine their layout), or a list
    of lists of matrices, where each list forms a row.

    -  ``ring`` - the base ring

    -  ``nrows`` - the number of block rows

    -  ``ncols`` - the number of block cols

    -  ``sub_matrices`` - matrices (see below for syntax)

    -  ``subdivide`` - boolean, whether or not to add
       subdivision information to the matrix

    -  ``sparse`` - boolean, whether to make the resulting matrix sparse


    EXAMPLES::

        sage: A = matrix(QQ, 2, 2, [3,9,6,10])
        sage: block_matrix([ [A, -A], [~A, 100*A] ])
        [    3     9|   -3    -9]
        [    6    10|   -6   -10]
        [-----------+-----------]
        [-5/12   3/8|  300   900]
        [  1/4  -1/8|  600  1000]

    If the number of submatrices in each row is the same,
    you can specify the submatrices as a single list too::

        sage: block_matrix(2, 2, [ A, A, A, A ])
        [ 3  9| 3  9]
        [ 6 10| 6 10]
        [-----+-----]
        [ 3  9| 3  9]
        [ 6 10| 6 10]


    One can use constant entries::

        sage: block_matrix([ [1, A], [0, 1] ])
        [ 1  0| 3  9]
        [ 0  1| 6 10]
        [-----+-----]
        [ 0  0| 1  0]
        [ 0  0| 0  1]

    A zero entry may represent any square or non-square zero matrix::

        sage: B = matrix(QQ, 1, 1, [ 1 ] )
        sage: C = matrix(QQ, 2, 2, [ 2, 3, 4, 5 ] )
        sage: block_matrix([ [B, 0], [0, C] ])
        [1|0 0]
        [-+---]
        [0|2 3]
        [0|4 5]

    One can specify the number of rows or columns as keywords too::

        sage: block_matrix([A, -A, ~A, 100*A], ncols=4)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

        sage: block_matrix([A, -A, ~A, 100*A], nrows=1)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

    It handles base rings nicely too::

        sage: R.<x> = ZZ['x']
        sage: block_matrix(2, 2, [1/2, A, 0, x-1])
        [  1/2     0|    3     9]
        [    0   1/2|    6    10]
        [-----------+-----------]
        [    0     0|x - 1     0]
        [    0     0|    0 x - 1]

        sage: block_matrix(2, 2, [1/2, A, 0, x-1]).parent()
        Full MatrixSpace of 4 by 4 dense matrices over Univariate Polynomial Ring in x over Rational Field

    Subdivisions are optional. If they are disabled, the columns need not line up::

        sage: B = matrix(QQ, 2, 3, range(6))
        sage: block_matrix([ [~A, B], [B, ~A] ], subdivide=False)
        [-5/12   3/8     0     1     2]
        [  1/4  -1/8     3     4     5]
        [    0     1     2 -5/12   3/8]
        [    3     4     5   1/4  -1/8]

    Without subdivisions it also deduces dimensions for scalars if possible::

        sage: C = matrix(ZZ, 1, 2, range(2))
        sage: block_matrix([ [ C, 0 ], [ 3, 4 ], [ 5, 6, C ] ], subdivide=False )
        [0 1 0 0]
        [3 0 4 0]
        [0 3 0 4]
        [5 6 0 1]

    If all submatrices are sparse (unless there are none at all), the result
    will be a sparse matrix. Otherwise it will be dense by default. The
    ``sparse`` keyword can be used to override this::

        sage: A = Matrix(ZZ, 2, 2, [0, 1, 0, 0], sparse=True)
        sage: block_matrix([ [ A ], [ A ] ]).parent()
        Full MatrixSpace of 4 by 2 sparse matrices over Integer Ring
        sage: block_matrix([ [ A ], [ A ] ], sparse=False).parent()
        Full MatrixSpace of 4 by 2 dense matrices over Integer Ring

    Consecutive zero submatrices are consolidated.  ::

        sage: B = matrix(2, range(4))
        sage: C = matrix(2, 8, range(16))
        sage: block_matrix(2, [[B,0,0,B],[C]], subdivide=False)
        [ 0  1  0  0  0  0  0  1]
        [ 2  3  0  0  0  0  2  3]
        [ 0  1  2  3  4  5  6  7]
        [ 8  9 10 11 12 13 14 15]

    Ambiguity is not tolerated.  ::

        sage: B = matrix(2, range(4))
        sage: C = matrix(2, 8, range(16))
        sage: block_matrix(2, [[B,0,B,0],[C]], subdivide=False)
        Traceback (most recent call last):
        ...
        ValueError: insufficient information to determine submatrix widths

    Historically, giving only a flat list of submatrices, whose number
    was a perfect square, would create a new matrix by laying out the submatrices
    in a square grid.  This behavior is now deprecated.  ::

        sage: A = matrix(2, 3, range(6))
        sage: B = matrix(3, 3, range(9))
        sage: block_matrix([A, A, B, B])
        doctest:...: DeprecationWarning: invocation of block_matrix with just a list whose length is a perfect square is deprecated. See the documentation for details.
        [0 1 2|0 1 2]
        [3 4 5|3 4 5]
        [-----+-----]
        [0 1 2|0 1 2]
        [3 4 5|3 4 5]
        [6 7 8|6 7 8]

    Historically, a flat list of matrices whose number is not a perfect square,
    with no specification of the number of rows or columns, would raise an error.
    This behavior continues, but could be removed when the deprecation above is
    completed.  ::

        sage: A = matrix(2, 3, range(6))
        sage: B = matrix(3, 3, range(9))
        sage: block_matrix([A, A, A, B, B, B])
        Traceback (most recent call last):
        ...
        ValueError: must specify nrows or ncols for non-square block matrix.
    """
    args = list(args)
    sparse = kwds.get('sparse',None)

    if len(args) == 0:
        if sparse is not None:
            return matrix_space.MatrixSpace(rings.ZZ, 0, 0, sparse=sparse)([])
        else:
            return matrix_space.MatrixSpace(rings.ZZ, 0, 0)([])

    if len(args) >= 1 and is_Ring(args[0]):
        # A ring is specified
        if kwds.get('ring', args[0]) != args[0]:
            raise ValueError("base ring specified twice and they are different")
        ring = args[0]
        args.pop(0)
    else:
        ring = kwds.get('ring', None)

    if len(args) >= 1:
        try:
            nrows = int(args[0])
            args.pop(0)
            if kwds.get('nrows', nrows) != nrows:
                raise ValueError("number of rows specified twice and they are different")
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
                raise ValueError("number of columns specified twice and they are different")
        except TypeError:
            ncols = kwds.get('ncols', None)
    else:
        ncols = kwds.get('ncols', None)

    # Now we've taken care of initial ring, nrows, and ncols arguments.

    # Now the rest of the arguments are a list of rows, a flat list of
    # matrices, or a single value.

    if len(args) == 0:
        args = [[]]
    if len(args) > 1:
        print args
        raise TypeError("invalid block_matrix invocation")

    sub_matrices = args[0]

    if is_Matrix(sub_matrices):
        # a single matrix (check nrows/ncols/ring)
        if (nrows is not None and nrows != 1) or \
           (ncols is not None and ncols != 1):
            raise ValueError("invalid nrows/ncols passed to block_matrix")
        if ring is not None:
            M = M.change_ring(ring)
        if sparse is not None and M.is_sparse() != sparse:
            M = M.sparse_matrix() if sparse else M.dense_matrix()
        return M

    if not isinstance(sub_matrices, (list, tuple)):
        raise TypeError("invalid block_matrix invocation")

    subdivide = kwds.get('subdivide', True)

    # Will we try to place the matrices in a rectangular grid?
    try_grid = True


    if len(sub_matrices) == 0:
        if (nrows is not None and nrows != 0) or \
           (ncols is not None and ncols != 0):
            raise ValueError("invalid nrows/ncols passed to block_matrix")
    elif isinstance(sub_matrices[0], (list, tuple)):
        # A list of lists: verify all elements are lists, and if
        # ncols is set, the lengths match.
        if nrows is not None and len(sub_matrices) != nrows:
            raise ValueError("invalid nrows passed to block_matrix")
        first_len = len(sub_matrices[0])
        if ncols is not None and first_len != ncols:
            raise ValueError("invalid ncols passed to block_matrix")
        same_length = all(isinstance(v, (list, tuple)) and len(v) == first_len for v in sub_matrices)
        if subdivide and not same_length:
            raise ValueError("list of rows is not valid (rows are wrong types or lengths)")
        try_grid = same_length
    else:
        # A flat list
        # determine the block dimensions
        n = ZZ(len(sub_matrices))
        if nrows is None:
            if ncols is None:
                if n.is_square():
                    import warnings
                    warnings.resetwarnings()
                    warnings.warn("invocation of block_matrix with just a list whose length is a perfect square is deprecated. See the documentation for details.", DeprecationWarning, stacklevel=2)
                    nrows = ncols = n.sqrt()
                else:
                    # this form (ie just a flat list) could be allowed once deprecated invocation (above) goes away
                    raise ValueError("must specify nrows or ncols for non-square block matrix.")
            else:
                nrows = int(n/ncols)
        elif ncols is None:
            ncols = int(n/nrows)
        if nrows * ncols != n:
            raise ValueError("given number of rows (%s), columns (%s) incompatible with number of submatrices (%s)" % (nrows, ncols, n))
        # Now create a list of lists from this
        sub_matrices = [ sub_matrices[i*ncols : (i+1)*ncols] for i in range(nrows) ]

    # At this point sub_matrices is a list of lists

    # determine the base ring and sparsity
    if ring is None:
        ring = ZZ
        for row in sub_matrices:
            for M in row:
                R = M.base_ring() if is_Matrix(M) else M.parent()
                if R is not ZZ:
                    ring = sage.categories.pushout.pushout(ring, R)

    if sparse is None:
        sparse = True
        for row in sub_matrices:
            for M in row:
                if sparse and is_Matrix(M) and not M.is_sparse():
                    sparse = False

    row_heights = None
    col_widths = None
    zero_widths = None
    total_width = None

    # We first try to place the matrices in a rectangular grid
    if try_grid:
        try:
            (row_heights, col_widths) = _determine_block_matrix_grid(sub_matrices)
        except ValueError as e:
            if subdivide:
                raise ValueError(e)


    if col_widths is None:
        # Try placing the matrices in rows instead
        # (Only if subdivide is False)
        (row_heights, zero_widths, total_width) = _determine_block_matrix_rows(sub_matrices)


    # Success, so assemble the final matrix

    big = None
    for i in range(len(sub_matrices)):
        R = sub_matrices[i]
        row = None
        for j in range(len(R)):
            M = R[j]

            if is_Matrix(M):
                if M.base_ring() is not ring:
                    M = M.change_ring(ring)
                if M.is_sparse() != sparse:
                    M = M.sparse_matrix() if sparse else M.dense_matrix()
            elif not M and zero_widths is not None:
                if zero_widths[i] > 0:
                    M = matrix(ring, row_heights[i], zero_widths[i], 0, sparse=sparse)
                    zero_widths[i] = 0
                else:
                    continue # zero-width matrix
            else:
                if zero_widths is not None:
                    M = matrix(ring, row_heights[i], row_heights[i], M, sparse=sparse)
                else:
                    M = matrix(ring, row_heights[i], col_widths[j], M, sparse=sparse)

            # append M to this row
            if row is None:
                row = M
            else:
                row = row.augment(M)

        # append row to final matrix
        if big is None:
            big = row
        else:
            big = big.stack(row)

    if big is None:
        if ring is None:
            ring = ZZ
        big = Matrix(ring, 0, 0)

    if subdivide:
        big.subdivide(running_total(row_heights[:-1]),
                      running_total(col_widths[:-1]))

    return big


@matrix_method
def block_diagonal_matrix(*sub_matrices, **kwds):
    """
    Create a block matrix whose diagonal block entries are given by
    sub_matrices, with zero elsewhere.

    See also :meth:`block_matrix`.

    EXAMPLES::

        sage: A = matrix(ZZ, 2, [1,2,3,4])
        sage: block_diagonal_matrix(A, A)
        [1 2|0 0]
        [3 4|0 0]
        [---+---]
        [0 0|1 2]
        [0 0|3 4]

    The sub-matrices need not be square::

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
    return block_matrix(n, n, entries, **kwds)


@matrix_method
def jordan_block(eigenvalue, size, sparse=False):
    r"""
    Returns the Jordan block for the given eigenvalue with given size.

    INPUT:

    -  ``eigenvalue`` - eigenvalue for the diagonal entries of the block
    -  ``size`` - size of the square matrix
    -  ``sparse`` - (default: False) - if True, return a sparse matrix


    EXAMPLE::

        sage: jordan_block(5, 3)
        [5 1 0]
        [0 5 1]
        [0 0 5]

    TESTS::

        sage: jordan_block(6.2, 'junk')
        Traceback (most recent call last):
        ...
        TypeError: size of Jordan block needs to be an integer, not junk
        sage: jordan_block(6.2, -1)
        Traceback (most recent call last):
        ...
        ValueError: size of Jordan block must be non-negative, not -1
    """
    try:
        size = ZZ(size)
    except TypeError:
        msg = "size of Jordan block needs to be an integer, not {0}"
        raise TypeError(msg.format(size))
    if size < 0:
        msg = "size of Jordan block must be non-negative, not {0}"
        raise ValueError(msg.format(size))
    block = diagonal_matrix([eigenvalue]*size, sparse=sparse)
    for i in xrange(size-1):
        block[i,i+1]=1
    return block


@matrix_method
def companion_matrix(poly, format='right'):
    r"""
    Create a companion matrix from a monic polynomial.

    INPUT:

    - ``poly`` - a univariate polynomial, or an iterable containing
      the coefficients of a polynomial, with low-degree coefficients first.
      The polynomial (or the polynomial implied by the coefficients) must
      be monic.  In other words, the leading coefficient must be one.
      A symbolic expression that might also be a polynomial is not
      proper input, see examples below.

    - ``format`` - default: 'right' - specifies one of four
      variations of a companion matrix.  Allowable values are
      'right', 'left', 'top' and 'bottom', which indicates which
      border of the matrix contains the negatives of the coefficients.

    OUTPUT:

    A square matrix with a size equal to the degree of the polynomial.
    The returned matrix has ones above, or below the diagonal, and the
    negatives of the coefficients along the indicated border of the
    matrix (excepting the leading one coefficient).
    See the first examples below for precise illustrations.

    EXAMPLES:

    Each of the four possibilities.  Notice that the coefficients are
    specified and their negatives become the entries of the matrix.  The
    leading one must be given, but is not used.  The permutation matrix
    ``P`` is the identity matrix, with the columns reversed.  The last three
    statements test the general relationships between the four variants.  ::

        sage: poly = [-2, -3, -4, -5, -6, 1]
        sage: R = companion_matrix(poly, format='right'); R
        [0 0 0 0 2]
        [1 0 0 0 3]
        [0 1 0 0 4]
        [0 0 1 0 5]
        [0 0 0 1 6]
        sage: L = companion_matrix(poly, format='left'); L
        [6 1 0 0 0]
        [5 0 1 0 0]
        [4 0 0 1 0]
        [3 0 0 0 1]
        [2 0 0 0 0]
        sage: B = companion_matrix(poly, format='bottom'); B
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        [2 3 4 5 6]
        sage: T = companion_matrix(poly, format='top'); T
        [6 5 4 3 2]
        [1 0 0 0 0]
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]

        sage: perm = Permutation([5, 4, 3, 2, 1])
        sage: P = perm.to_matrix()
        sage: L == P*R*P
        True
        sage: B == R.transpose()
        True
        sage: T == P*R.transpose()*P
        True

    A polynomial may be used as input, however a symbolic expression,
    even if it looks like a polynomial, is not regarded as such when used
    as input to this routine.  Obtaining the list of coefficients from a
    symbolic polynomial is one route to the companion matrix. ::

        sage: x = polygen(QQ, 'x')
        sage: p = x^3 - 4*x^2 + 8*x - 12
        sage: companion_matrix(p)
        [ 0  0 12]
        [ 1  0 -8]
        [ 0  1  4]

        sage: y = var('y')
        sage: q = y^3 -2*y + 1
        sage: companion_matrix(q)
        Traceback (most recent call last):
        ...
        TypeError: input must be a polynomial (not a symbolic expression, see docstring), or other iterable, not y^3 - 2*y + 1

        sage: coeff_list = [q(y=0)] + [q.coeff(y^k) for k in range(1, q.degree(y)+1)]
        sage: coeff_list
        [1, -2, 0, 1]
        sage: companion_matrix(coeff_list)
        [ 0  0 -1]
        [ 1  0  2]
        [ 0  1  0]

    The minimal polynomial of a companion matrix is equal to the
    polynomial used to create it.  Used in a block diagonal
    construction, they can be used to create matrices with
    any desired minimal polynomial, or characteristic polynomial.  ::

        sage: t = polygen(QQ, 't')
        sage: p = t^12 - 7*t^4 + 28*t^2 - 456
        sage: C = companion_matrix(p, format='top')
        sage: q = C.minpoly(var=t); q
        t^12 - 7*t^4 + 28*t^2 - 456
        sage: p == q
        True

        sage: p = t^3 + 3*t - 8
        sage: q = t^5 + t - 17
        sage: A = block_diagonal_matrix( companion_matrix(p),
        ...                              companion_matrix(p^2),
        ...                              companion_matrix(q),
        ...                              companion_matrix(q) )
        sage: A.charpoly(var=t).factor()
        (t^3 + 3*t - 8)^3 * (t^5 + t - 17)^2
        sage: A.minpoly(var=t).factor()
        (t^3 + 3*t - 8)^2 * (t^5 + t - 17)

    TESTS::

        sage: companion_matrix([4, 5, 1], format='junk')
        Traceback (most recent call last):
        ...
        ValueError: format must be 'right', 'left', 'top' or 'bottom', not junk

        sage: companion_matrix(sin(x))
        Traceback (most recent call last):
        ...
        TypeError: input must be a polynomial (not a symbolic expression, see docstring), or other iterable, not sin(x)

        sage: companion_matrix([2, 3, 896])
        Traceback (most recent call last):
        ...
        ValueError: polynomial (or the polynomial implied by coefficients) must be monic, not a leading coefficient of 896

        sage: F.<a> = GF(2^2)
        sage: companion_matrix([4/3, a+1, 1])
        Traceback (most recent call last):
        ...
        TypeError: unable to find common ring for coefficients from polynomial

        sage: A = companion_matrix([1])
        sage: A.nrows(); A.ncols()
        0
        0

        sage: A = companion_matrix([])
        Traceback (most recent call last):
        ...
        ValueError: polynomial cannot be specified by an empty list

    AUTHOR:

    - Rob Beezer (2011-05-19)
    """
    import sage.matrix.constructor
    if format not in ['right', 'left', 'top', 'bottom']:
        raise ValueError("format must be 'right', 'left', 'top' or 'bottom', not {0}".format(format))
    try:
        poly = list(poly)
    except TypeError:
        raise TypeError('input must be a polynomial (not a symbolic expression, see docstring), or other iterable, not {0}'.format(poly))
    n = len(poly) - 1
    if n == -1:
        raise ValueError('polynomial cannot be specified by an empty list')
    if not poly[n] == 1:
        raise ValueError('polynomial (or the polynomial implied by coefficients) must be monic, not a leading coefficient of {0}'.format(poly[n]))
    entries = [0]*(n*n)
    # 1's below diagonal, or above diagonal
    if format in ['right', 'top']:
        for i in range(n-1):
            entries[(i+1)*n + i] = 1
    else:
        for i in range(n-1):
            entries[i*n + i+1] = 1
    # right side, left side (reversed), bottom edge, top edge (reversed)
    if format == 'right':
        for i in range(n):
            entries[i*n + n-1] = -poly[i]
    elif format == 'left':
        for i in range(n):
            entries[(n-1-i)*n + 0] = -poly[i]
    elif format == 'bottom':
        for i in range(n):
            entries[(n-1)*n + i] = -poly[i]
    elif format == 'top':
        for i in range(n):
            entries[0*n + n-1-i] = -poly[i]
    try:
        M = sage.matrix.constructor.matrix(n, n, entries)
    except TypeError:
        raise TypeError("unable to find common ring for coefficients from polynomial")
    return M

@matrix_method
def random_rref_matrix(parent, num_pivots):
    r"""
    Generate a matrix in reduced row-echelon form with a specified number of non-zero rows.

    INPUT:

    - ``parent`` - A matrix space specifying the base ring, dimensions and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``num_pivots`` - The number of non-zero rows in the result, i.e. the rank.

    OUTPUT:

    A matrix in reduced row echelon form with ``num_pivots`` non-zero rows. If the
    base ring is `ZZ` or `QQ` then the entries are all integers.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='echelon_form'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    Matrices generated are in reduced row-echelon form with specified rank. If the
    base ring is `QQ` the result has only integer entries.  ::

        sage: from sage.matrix.constructor import random_rref_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5, 6)
        sage: A=random_rref_matrix(matrix_space, num_pivots=4); A # random
        [ 1  0  0 -6  0 -3]
        [ 0  1  0  2  0  3]
        [ 0  0  1 -4  0 -2]
        [ 0  0  0  0  1  3]
        [ 0  0  0  0  0  0]
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (5, 6)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 5, 6)
        True
        sage: A.rank()
        4
        sage: A==A.rref()
        True

    Matrices can be generated over other exact rings. ::

        sage: B=random_matrix(FiniteField(7), 4, 4, algorithm='echelon_form', num_pivots=3); B # random
        [1 0 0 0]
        [0 1 0 6]
        [0 0 1 4]
        [0 0 0 0]
        sage: B.rank() == 3
        True
        sage: B.base_ring()
        Finite Field of size 7
        sage: B==B.rref()
        True

    TESTS:

    Rank of a matrix must be an integer. ::

        sage: random_matrix(QQ, 120, 56, algorithm='echelon_form', num_pivots=61/2)
        Traceback (most recent call last):
        ...
        TypeError: the number of pivots must be an integer.

    Matrices must be generated over exact fields. ::

        sage: random_matrix(RR, 40, 88, algorithm='echelon_form', num_pivots=39)
        Traceback (most recent call last):
        ...
        TypeError: the base ring must be exact.

    Matrices must have the number of pivot columns be less than or equal to the number of rows. ::

        sage: C=random_matrix(ZZ, 6,4, algorithm='echelon_form', num_pivots=7); C
        Traceback (most recent call last):
        ...
        ValueError: number of pivots cannot exceed the number of rows or columns.

    Matrices must have the number of pivot columns be less than or equal to the number of columns. ::

        sage: D=random_matrix(QQ, 1,3, algorithm='echelon_form', num_pivots=5); D
        Traceback (most recent call last):
        ...
        ValueError: number of pivots cannot exceed the number of rows or columns.

    Matrices must have the number of pivot columns be greater than zero. ::

        sage: random_matrix(QQ, 5, 4, algorithm='echelon_form', num_pivots=-1)
        Traceback (most recent call last):
        ...
        ValueError: the number of pivots must be zero or greater.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    import sage.gsl.probability_distribution as pd
    from sage.misc.prandom import randint

    try:
        num_pivots=ZZ(num_pivots)
    except TypeError:
        raise TypeError("the number of pivots must be an integer.")
    if num_pivots<0:
        raise ValueError("the number of pivots must be zero or greater.")
    ring = parent.base_ring()
    if not ring.is_exact():
        raise TypeError("the base ring must be exact.")
    num_row = parent.nrows()
    num_col = parent.ncols()
    if num_pivots>num_row or num_pivots>num_col:
        raise ValueError("number of pivots cannot exceed the number of rows or columns.")
    else:
        one=ring.one()
        # Create a matrix of the desired size to be modified and then returned.
        return_matrix=copy(parent.zero_matrix())
        pivots=[0] #Force first column to be a pivot. No harm if no pivots at all.
        # Probability distribution for the placement of leading one's.
        pivot_generator=pd.RealDistribution("beta",[1.6,4.3])
        while len(pivots)<num_pivots:
            pivot_column=int(pivot_generator.get_random_element()*num_col)
            if pivot_column not in pivots:
                pivots.append(pivot_column)
        pivots.sort()
        pivot_row=0
        # Use the list of pivot columns to set the pivot entries of the return_matrix to leading ones.
        while pivot_row<num_pivots:
            return_matrix[pivot_row,pivots[pivot_row]]=one
            pivot_row+=1
        if ring==QQ or ring==ZZ:
            # Keep track of the non-pivot columns by using the pivot_index, start at the first column to
            # the right of the initial pivot column, go until the first column to the left of the next
            # pivot column.
            for pivot_index in range(num_pivots-1):
                for non_pivot_column_index in range(pivots[pivot_index]+1,pivots[pivot_index+1]):
                    entry_generator1=pd.RealDistribution("beta",[6,4])
                    # Experimental distribution used to generate the values.
                    for non_pivot_column_entry in range(pivot_index+1):
                        sign1=(2*randint(0,1)-1)
                        return_matrix[non_pivot_column_entry,non_pivot_column_index]=sign1*int(entry_generator1.get_random_element()*((1-non_pivot_column_entry/return_matrix.ncols())*7))
            # Use index to fill entries of the columns to the right of the last pivot column.
            for rest_non_pivot_column in range(pivots[num_pivots-1]+1,num_col):
                entry_generator2=pd.RealDistribution("beta",[2.6,4])
                # experimental distribution to generate small values.
                for rest_entries in range(num_pivots):
                    sign2=(2*randint(0,1)-1)
                    return_matrix[rest_entries,rest_non_pivot_column]=sign2*int(entry_generator2.get_random_element()*5)
        else:
            for pivot_index in range(num_pivots-1):
                for non_pivot_column_index in range(pivots[pivot_index]+1,pivots[pivot_index+1]):
                    for non_pivot_column_entry in range(pivot_index+1):
                            return_matrix[non_pivot_column_entry,non_pivot_column_index]=ring.random_element()
            for rest_non_pivot_column in range(pivots[num_pivots-1]+1,num_col):
                for rest_entries in range(num_pivots):
                    return_matrix[rest_entries,rest_non_pivot_column]=ring.random_element()
    return return_matrix

@matrix_method
def random_echelonizable_matrix(parent, rank, upper_bound=None, max_tries=100):
    r"""
    Generate a matrix of a desired size and rank, over a desired ring, whose reduced
    row-echelon form has only integral values.

    INPUT:

    - ``parent`` - A matrix space specifying the base ring, dimensions and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``rank`` - Rank of result, i.e the number of non-zero rows in the
      reduced row echelon form.

    - ``upper_bound`` - If designated, size control of the matrix entries is desired.
      Set ``upper_bound`` to 1 more than the maximum value entries can achieve.
      If None, no size control occurs. But see the warning below.  (default: None)

    - ``max_tries`` - If designated, number of tries used to generate each new random row;
      only matters when upper_bound!=None. Used to prevent endless looping. (default: 100)

    OUTPUT:

    A matrix not in reduced row-echelon form with the desired dimensions and properties.

    .. warning::

        When ``upper_bound`` is set, it is possible for this constructor to
        fail with a ``ValueError``.  This may happen when the ``upper_bound``,
        ``rank`` and/or matrix dimensions are all so small that it becomes
        infeasible or unlikely to create the requested matrix.  If you *must*
        have this routine return successfully, do not set ``upper_bound``.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='echelonizable'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    Generated matrices have the desired dimensions, rank and entry size. The
    matrix in reduced row-echelon form has only integer entries. ::

        sage: from sage.matrix.constructor import random_echelonizable_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5, 6)
        sage: A=random_echelonizable_matrix(matrix_space, rank=4, upper_bound=40); A # random
        [  1  -1   1  -3  -4   6]
        [  5  -4   0   8   4  19]
        [ -3   3  -2   4   7 -16]
        [ -4   5  -7  26  31 -31]
        [  2  -3   4 -11 -14  17]
        sage: A.rank()
        4
        sage: max(map(abs,A.list()))<40
        True
        sage: A.rref()==A.rref().change_ring(ZZ)
        True

    An example with default settings (i.e. no entry size control). ::

        sage: C=random_matrix(QQ, 6, 7, algorithm='echelonizable', rank=5); C # random
        [  1   0   5  -2 -26 -16   0]
        [ -3   1 -19   6  97  61   1]
        [  0   4 -15  -1  71  50   3]
        [  2   4  -9   0  39  25   8]
        [  2   2   3  -3 -18  -9   3]
        [ -3  -4  -2  14  14  -6   4]
        sage: C.rank()
        5
        sage: C.rref()==C.rref().change_ring(ZZ)
        True

    A matrix without size control may have very large entry sizes. ::

        sage: D=random_matrix(ZZ, 7, 8, algorithm='echelonizable', rank=6); D # random
        [    9   -53  -255    45 -1519  4043  9819  3324]
        [    3   -14   -64     8  -369   972  2350   810]
        [    2   -14   -65     9  -377  1000  2420   829]
        [    4   -24  -116    21  -693  1846  4485  1516]
        [   -3    14    68   -16   426 -1134 -2767  -919]
        [   -5    21    92   -13   548 -1432 -3466 -1183]
        [    1    -9   -42     7  -254   670  1624   547]

    Matrices can be generated over any exact ring. ::

        sage: F.<a>=GF(2^3)
        sage: B=random_matrix(F, 4, 5, algorithm='echelonizable', rank=4, upper_bound=None); B # random
        [      a + 1 a^2 + a + 1         a^2           0           1]
        [          1           a       a + 1         a^2     a^2 + a]
        [          a         a^2 a^2 + a + 1       a + 1           1]
        [          1           0 a^2 + a + 1           0     a^2 + 1]
        sage: B.rank()
        4

    Square matrices over ZZ or QQ with full rank are always unimodular. ::

        sage: E=random_matrix(QQ, 7, 7, algorithm='echelonizable', rank=7); E # random
        [  1  -1   5  12 -24 -41  47]
        [  0   1  -1   3   0 -11  40]
        [  1  -1   6   6 -19 -20 -11]
        [ -2   1 -10 -12  35  44  -4]
        [  3  -1   9   7 -35 -40 -18]
        [  0   0   0  -4   4  13 -32]
        [  3  -3  11   6 -33 -31 -35]
        sage: det(E)
        1

    TESTS:

    Matrices must have a rank zero or greater, and less than
    both the number of rows and the number of columns. ::

        sage: random_matrix(QQ, 3, 4, algorithm='echelonizable', rank=-1)
        Traceback (most recent call last):
        ...
        ValueError: matrices must have rank zero or greater.
        sage: random_matrix(QQ, 3, 8, algorithm='echelonizable', rank=4)
        Traceback (most recent call last):
        ...
        ValueError: matrices cannot have rank greater than min(ncols,nrows).
        sage: random_matrix(QQ, 8, 3, algorithm='echelonizable', rank=4)
        Traceback (most recent call last):
        ...
        ValueError: matrices cannot have rank greater than min(ncols,nrows).

    The base ring must be exact. ::

        sage: random_matrix(RR, 3, 3, algorithm='echelonizable', rank=2)
        Traceback (most recent call last):
        ...
        TypeError: the base ring must be exact.

    Works for rank==1, too. ::

        sage: random_matrix( QQ, 3, 3, algorithm='echelonizable', rank=1).ncols()
        3


    AUTHOR:

    Billy Wonderly (2010-07)
    """

    from sage.misc.prandom import randint

    ring = parent.base_ring()
    rows = parent.nrows()
    if rank<0:
        raise ValueError("matrices must have rank zero or greater.")
    if rank>min(rows,parent.ncols()):
        raise ValueError("matrices cannot have rank greater than min(ncols,nrows).")
    matrix = random_rref_matrix(parent, rank)

    # Entries of matrices over the ZZ or QQ can get large, entry size is regulated by finding the largest
    # entry of the resultant matrix after addition of scalar multiple of a row.
    if ring==QQ or ring==ZZ:
        # If upper_bound is not set, don't control entry size.
        if upper_bound is None:
        # If size control is not desired, the routine will run slightly faster, particularly with large matrices.
            for pivots in range(rank-1,-1,-1):
                row_index=0
                while row_index<rows:
                    if pivots==row_index:
                        row_index+=1
                    if pivots!=row_index and row_index!=rows:
                        matrix.add_multiple_of_row(row_index,matrix.pivot_rows()[pivots],randint(-5,5))
                        row_index+=1
            if rows>1:
                matrix.add_multiple_of_row(0,randint(1,rows-1),randint(-3,3))
        else:
            if rank==1:  # would be better just to have a special generator...
               tries = 0
               while max(map(abs,matrix.list()))>=upper_bound:
                  matrix = random_rref_matrix(parent, rank)
                  tries += 1
                  if tries > max_tries: # to prevent endless attempts
                     raise ValueError("tried "+str(max_tries)+" times to get a rank 1 random matrix. Try bigger upper_bound?")
               matrix_copy = matrix

            rrr = range(len(matrix.pivots())-1,-1,-1)
            for pivots in rrr:
            # keep track of the pivot column positions from the pivot column with the largest index to
            # the one with the smallest.
                row_index=0
                tries = 0
                while row_index<rows:
                    # To each row in a pivot column add a scalar multiple of the pivot row.
                    # for full rank, square matrices, using only this row operation preserves the determinant of 1.
                    if pivots!=row_index:
                    # To ensure a leading one is not removed by the addition of the pivot row by its
                    # additive inverse.
                        matrix_copy=matrix.with_added_multiple_of_row(row_index,matrix.pivot_rows()[pivots],randint(-5,5))
                        tries += 1
                        # Range for scalar multiples determined experimentally.
                    if max(map(abs,matrix_copy.list()))<upper_bound:
                    # Continue if the the largest entry after a row operation is within the bound.
                        matrix=matrix_copy
                        row_index+=1
                        tries = 0
                    if tries > max_tries: # to prevent endless unsuccessful row adding
                        raise ValueError("tried "+str(max_tries)+" times to get row number "+str(row_index)+". Try bigger upper_bound?")
            # The leading one in row one has not been altered, so add a scalar multiple of a random row
            # to row one.
            row1=0
            if rows>1:
                while row1<1:
                    matrix_copy=matrix.with_added_multiple_of_row(0,randint(1,rows-1),randint(-3,3))
                    if max(map(abs,matrix_copy.list()))<upper_bound:
                        matrix=matrix_copy
                        row1+=1
    # If the matrix generated over a different ring, random elements from the designated ring are used as and
    # the routine is run similarly to the size unchecked version for rationals and integers.
    else:
        for pivots in range(rank-1,-1,-1):
            row_index=0
            while row_index<rows:
                if pivots==row_index:
                    row_index+=1
                if pivots!=row_index and row_index!=rows:
                    matrix.add_multiple_of_row(row_index,matrix.pivot_rows()[pivots],ring.random_element())
                    row_index+=1
        if rows>1:
            matrix.add_multiple_of_row(0,randint(1,rows-1),ring.random_element())
    return matrix

@matrix_method
def random_subspaces_matrix(parent, rank=None):
    r"""
    Create a matrix of the designated size and rank whose right and
    left null spaces, column space, and row space have desirable
    properties that simplify the subspaces.

    INPUT:

    - ``parent`` - A matrix space specifying the base ring, dimensions, and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``rank`` - The desired rank of the return matrix (default: None).

    OUTPUT:

    A matrix whose natrual basis vectors for its four subspaces, when
    computed, have reasonably sized, integral valued, entries.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='subspaces'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A 6x8 matrix with designated rank of 3.  The four subspaces are
    determined using one simple routine in which we augment the
    original matrix with the equal row dimension identity matrix.  The
    resulting matrix is then put in reduced row-echelon form and the
    subspaces can then be determined by analyzing subdivisions of this
    matrix. See the four subspaces routine in [BEEZER]_ for more. ::

        sage: from sage.matrix.constructor import random_subspaces_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 6, 8)
        sage: B=random_subspaces_matrix(matrix_space, rank=3); B # random
        [ 113  339  -46  218 -243  641 -269 -306]
        [ -33  -99   13  -63   69 -185   77   90]
        [  35  105  -14   67  -74  197  -82  -95]
        [ -18  -54    7  -34   37 -100   41   49]
        [ -26  -78   10  -49   53 -144   59   71]
        [  -8  -24    3  -15   16  -44   18   22]
        sage: B.rank()
        3
        sage: B.nullity()
        3
        sage: (B.nrows(), B.ncols())
        (6, 8)
        sage: all([x in ZZ for x in B.list()])
        True
        sage: B_expanded=B.augment(identity_matrix(6)).rref()
        sage: all([x in ZZ for x in B_expanded.list()])
        True
        sage: B_expanded # random
        [  1   3   0   0   1   1   3  -2   0   0  -3   0  -9  16]
        [  0   0   1   0   3  -2  -1  -3   0   0   2   0  11 -27]
        [  0   0   0   1  -1   2  -3  -1   0   0   2   0   7 -14]
        [  0   0   0   0   0   0   0   0   1   0  -5   0  -3   2]
        [  0   0   0   0   0   0   0   0   0   1   1   0   1  -3]
        [  0   0   0   0   0   0   0   0   0   0   0   1  -1   1]

    Check that we fixed Trac #10543 (echelon forms should be immutable)::

        sage: B_expanded.is_immutable()
        True

    We want to modify B_expanded, so replace it with a copy::

        sage: B_expanded = copy(B_expanded)
        sage: B_expanded.subdivide(B.nrows()-B.nullity(),B.ncols());B_expanded # random
        [  1   3   0   0   1   1   3  -2|  0   0  -3   0  -9  16]
        [  0   0   1   0   3  -2  -1  -3|  0   0   2   0  11 -27]
        [  0   0   0   1  -1   2  -3  -1|  0   0   2   0   7 -14]
        [-------------------------------+-----------------------]
        [  0   0   0   0   0   0   0   0|  1   0  -5   0  -3   2]
        [  0   0   0   0   0   0   0   0|  0   1   1   0   1  -3]
        [  0   0   0   0   0   0   0   0|  0   0   0   1  -1   1]
        sage: C=B_expanded.subdivision(0,0)
        sage: C # random
        [ 1  3  0  0  1  1  3 -2]
        [ 0  0  1  0  3 -2 -1 -3]
        [ 0  0  0  1 -1  2 -3 -1]
        sage: L=B_expanded.subdivision(1,1)
        sage: L # random
        [ 1  0 -5  0 -3  2]
        [ 0  1  1  0  1 -3]
        [ 0  0  0  1 -1  1]
        sage: B.right_kernel()==C.right_kernel()
        True
        sage: B.row_space()==C.row_space()
        True
        sage: B.column_space()==L.right_kernel()
        True
        sage: B.left_kernel()==L.row_space()
        True

    A matrix to show that the null space of the L matrix is the column space of the starting matrix. ::

        sage: A=random_matrix(QQ, 5, 7, algorithm='subspaces', rank=None); A # random
        [-31  12  -9 -27  21   2 -15]
        [105 -24   6 103 -30 -34  79]
        [ 29  -9   5  26 -14  -5  17]
        [233 -55  16 228 -71 -73 173]
        [-42  10  -3 -41  13  13 -31]
        sage: (A.nrows(), A.ncols())
        (5, 7)
        sage: all([x in ZZ for x in A.list()])
        True
        sage: A.nullity() # random
        1
        sage: A_expanded=A.augment(identity_matrix(5)).rref()
        sage: A_expanded # random
        [  1   0   0   0   0   1   0   0   3   7  25 151]
        [  0   1   0   0   1   2   1   0   5  21  84 493]
        [  0   0   1   0  -1   2   0   0   2  13  53 308]
        [  0   0   0   1   0  -1   1   0  -2  -3  -9 -57]
        [  0   0   0   0   0   0   0   1  -3   1   1  -2]
        sage: all([x in ZZ for x in A_expanded.list()])
        True
        sage: C=A_expanded.submatrix(0,0,A.nrows()-A.nullity(),A.ncols())
        sage: L=A_expanded.submatrix(A.nrows()-A.nullity(),A.ncols())
        sage: A.right_kernel()==C.right_kernel()
        True
        sage: A.row_space()==C.row_space()
        True
        sage: A.column_space()==L.right_kernel()
        True
        sage: A.left_kernel()==L.row_space()
        True

    TESTS:

    The designated rank of the L matrix cannot be greater than the
    number of desired rows, nor can the rank be negative. ::

        sage: random_matrix(QQ, 19, 20, algorithm='subspaces', rank=21)
        Traceback (most recent call last):
        ...
        ValueError: rank cannot exceed the number of rows or columns.
        sage: random_matrix(QQ, 19, 20, algorithm='subspaces', rank=-1)
        Traceback (most recent call last):
        ...
        ValueError: matrices must have rank zero or greater.

    REFERENCES:

        .. [BEEZER] `A First Course in Linear Algebra <http://linear.ups.edu/>`_.
           Robert A. Beezer, accessed 15 July 2010.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    import sage.gsl.probability_distribution as pd

    ring = parent.base_ring()
    rows = parent.nrows()
    columns = parent.ncols()

    # If rank is not designated, generate using probability distribution
    # skewing to smaller numbers, always at least 1.
    if rank is None:
        left_nullity_generator = pd.RealDistribution("beta", [1.4, 5.5])
        nullity = int(left_nullity_generator.get_random_element()*(rows-1) + 1)
        rank = rows - nullity
    if rank<0:
        raise ValueError("matrices must have rank zero or greater.")
    if rank > rows or rank > columns:
        raise ValueError("rank cannot exceed the number of rows or columns.")
    nullity = rows - rank
    B = random_matrix(ring, rows, columns, algorithm='echelon_form',
            num_pivots=rank)

    # Create a nonsingular matrix whose columns will be used to stack a matrix
    # over the L matrix, forming a nonsingular matrix.
    K_nonzero_columns = random_matrix(ring, rank, rank,
            algorithm='echelonizable', rank=rank)
    K = matrix(QQ, rank, rows)
    L = random_matrix(ring, nullity, rows, algorithm='echelon_form',
            num_pivots=nullity)
    for column in range(len(L.nonpivots())):
        for entry in range(rank):
            K[entry, L.nonpivots()[column]] = K_nonzero_columns[entry, column]
    J = K.stack(L)

    # By multiplying the B matrix by J.inverse() we hide the B matrix of the
    # solution using row operations required to change the solution K matrix to
    # the identity matrix.
    return J.inverse()*B

@matrix_method
def random_unimodular_matrix(parent, upper_bound=None, max_tries=100):
    r"""
    Generate a random unimodular (determinant 1) matrix of a desired size over a desired ring.

    INPUT:

    - ``parent`` - A matrix space specifying the base ring, dimensions
      and representation (dense/sparse) for the result.  The base ring
      must be exact.

    - ``upper_bound`` - For large matrices over QQ or ZZ,
      ``upper_bound`` is the largest value matrix entries can achieve.  But
      see the warning below.

    - ``max_tries`` - If designated, number of tries used to generate each new random row;
      only matters when upper_bound!=None. Used to prevent endless looping. (default: 100)

    A matrix not in reduced row-echelon form with the desired dimensions and properties.

    OUTPUT:

    An invertible matrix with the desired properties and determinant 1.

    .. warning::

        When ``upper_bound`` is set, it is possible for this constructor to
        fail with a ``ValueError``.  This may happen when the ``upper_bound``,
        ``rank`` and/or matrix dimensions are all so small that it becomes
        infeasible or unlikely to create the requested matrix.  If you *must*
        have this routine return successfully, do not set ``upper_bound``.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='unimodular'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A matrix size 5 over QQ. ::

        sage: from sage.matrix.constructor import random_unimodular_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5)
        sage: A=random_unimodular_matrix(matrix_space); A # random
        [  -8   31   85  148 -419]
        [   2   -9  -25  -45  127]
        [  -1   10   30   65 -176]
        [  -3   12   33   58 -164]
        [   5  -21  -59 -109  304]
        sage: det(A)
        1

    A matrix size 6 with entries no larger than 50. ::

        sage: B=random_matrix(ZZ, 7, algorithm='unimodular', upper_bound=50);B # random
        [ -1   0   3   1  -2   2   9]
        [  1   2  -5   0  14  19 -49]
        [ -3  -2  12   5  -6  -4  24]
        [  1   2  -9  -3   3   4  -7]
        [ -2  -1   7   2  -8  -5  31]
        [  2   2  -6  -3   8  16 -32]
        [  1   2  -9  -2   5   6 -12]
        sage: det(B)
        1

    A matrix over the number Field in `y` with defining polynomial `y^2-2y-2`. ::

        sage: y = var('y')
        sage: K=NumberField(y^2-2*y-2,'y')
        sage: C=random_matrix(K, 3, algorithm='unimodular');C # random
        [   2*y - 33 681*y - 787   31*y - 37]
        [      y + 6 -155*y + 83    -7*y + 4]
        [         -y   24*y + 51       y + 3]
        sage: det(C)
        1

    TESTS:

    Unimodular matrices are square. ::

        sage: random_matrix(QQ, 5, 6, algorithm='unimodular')
        Traceback (most recent call last):
        ...
        TypeError: a unimodular matrix must be square.

    Only matrices over ZZ and QQ can have size control. ::

        sage: F.<a>=GF(5^7)
        sage: random_matrix(F, 5, algorithm='unimodular', upper_bound=20)
        Traceback (most recent call last):
        ...
        TypeError: only matrices over ZZ or QQ can have size control.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    ring=parent.base_ring()
    size=parent.nrows()
    if parent.nrows()!=parent.ncols():
        raise TypeError("a unimodular matrix must be square.")
    if upper_bound is not None and (ring!=ZZ and ring!=QQ):
        raise TypeError("only matrices over ZZ or QQ can have size control.")
    if upper_bound is None:
        # random_echelonizable_matrix() always returns a determinant one matrix if given full rank.
        return random_matrix(ring, size, algorithm='echelonizable', rank=size)
    elif upper_bound is not None and (ring==ZZ or ring==QQ):
        return random_matrix(ring, size,algorithm='echelonizable',rank=size, upper_bound=upper_bound, max_tries=max_tries)


@matrix_method
def random_diagonalizable_matrix(parent,eigenvalues=None,dimensions=None):
    """
    Create a random matrix that diagonalizes nicely.

    To be used as a teaching tool.  Return matrices have only real
    eigenvalues.

    INPUT:

    If eigenvalues and dimensions are not specified in a list,
    they will be assigned randomly.

    - ``parent`` - the desired size of the square matrix.

    - ``eigenvalues`` - the list of desired eigenvalues (default=None).

    - ``dimensions`` - the list of dimensions corresponding to each
      eigenspace (default=None).

    OUTPUT:

    A square, diagonalizable, matrix with only integer entries. The
    eigenspaces of this matrix, if computed by hand, give basis
    vectors with only integer entries.

    .. note::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='diagonalizable'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A diagonalizable matrix, size 5. ::

        sage: from sage.matrix.constructor import random_diagonalizable_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5)
        sage: A=random_diagonalizable_matrix(matrix_space); A # random
        [ 10  18   8   4 -18]
        [ 20  10   8   4 -16]
        [-60 -54 -22 -12  18]
        [-60 -54 -24  -6   6]
        [-20 -18  -8  -4   8]
        sage: A.eigenvalues() # random
        [10,6,2,-8,-10]
        sage: S=A.right_eigenmatrix()[1]; S # random
        [ 1  1  1  1  0]
        [ 1  1  1  0  1]
        [-3 -3 -4 -3 -3]
        [-3 -4 -3 -3 -3]
        [-1 -1 -1 -1 -1]
        sage: S_inverse=S.inverse(); S_inverse # random
        [ 1  1  1  1 -5]
        [ 0  0  0 -1  3]
        [ 0  0 -1  0  3]
        [ 0 -1  0  0 -1]
        [-1  0  0  0 -1]
        sage: S_inverse*A*S # random
        [ 10   0   0   0   0]
        [  0   6   0   0   0]
        [  0   0   2   0   0]
        [  0   0   0  -8   0]
        [  0   0   0   0 -10]

    A diagonalizable matrix with eigenvalues and dimensions designated,
    with a check that if eigenvectors were calculated by hand
    entries would all be integers. ::

        sage: B=random_matrix(QQ, 6, algorithm='diagonalizable', eigenvalues=[-12,4,6],dimensions=[2,3,1]); B # random
        [ -52   32  240 -464  -96 -520]
        [   6    4  -48   72   36   90]
        [  46  -32 -108  296  -12  274]
        [  24  -16  -64  164    0  152]
        [  18  -16    0   72  -48   30]
        [   2    0  -16   24   12   34]
        sage: all([x in ZZ for x in (B-(-12*identity_matrix(6))).rref().list()])
        True
        sage: all([x in ZZ for x in (B-(4*identity_matrix(6))).rref().list()])
        True
        sage: all([x in ZZ for x in (B-(6*identity_matrix(6))).rref().list()])
        True
        sage: S=B.right_eigenmatrix()[1]; S_inverse=S.inverse(); S_inverse*B*S # random
        [  6   0   0   0   0   0]
        [  0 -12   0   0   0   0]
        [  0   0 -12   0   0   0]
        [  0   0   0   4   0   0]
        [  0   0   0   0   4   0]
        [  0   0   0   0   0   4]

    TESTS:

    Eigenvalues must all be integers. ::

        sage: random_matrix(QQ,3,algorithm='diagonalizable', eigenvalues=[2+I,2-I,2],dimensions=[1,1,1])
        Traceback (most recent call last):
        ...
        TypeError: eigenvalues must be integers.

    Diagonal matrices must be square. ::

        sage: random_matrix(QQ, 5, 7, algorithm='diagonalizable', eigenvalues=[-5,2,-3], dimensions=[1,1,3])
        Traceback (most recent call last):
        ...
        TypeError: a diagonalizable matrix must be square.

    A list of eigenvalues must be accompanied with a list of dimensions. ::

        sage: random_matrix(QQ,10,algorithm='diagonalizable',eigenvalues=[4,8])
        Traceback (most recent call last):
        ...
        ValueError: the list of eigenvalues must have a list of dimensions corresponding to each eigenvalue.

    A list of dimensions must be accompanied with a list of eigenvalues. ::

        sage: random_matrix(QQ, 10,algorithm='diagonalizable',dimensions=[2,2,4,2])
        Traceback (most recent call last):
        ...
        ValueError: the list of dimensions must have a list of corresponding eigenvalues.

    The sum of the eigenvalue dimensions must equal the size of the matrix. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6,-1],dimensions=[2,3,5,1])
        Traceback (most recent call last):
        ...
        ValueError: the size of the matrix must equal the sum of the dimensions.

    Each eigenspace dimension must be at least 1. ::

        sage: random_matrix(QQ,9,algorithm='diagonalizable',eigenvalues=[-15,22,8,-4,90,12],dimensions=[4,2,2,4,-3,0])
        Traceback (most recent call last):
        ...
        ValueError: eigenspaces must have a dimension of at least 1.

    Each eigenvalue must have a corresponding eigenspace dimension. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6,-1],dimensions=[4,3,5])
        Traceback (most recent call last):
        ...
        ValueError: each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.

    Each dimension must have an eigenvalue paired to it. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6],dimensions=[2,3,5,2])
        Traceback (most recent call last):
        ...
        ValueError: each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.

    TODO:

    Modify the routine to allow for complex eigenvalues.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    from sage.misc.prandom import randint

    size=parent.nrows()
    if parent.nrows()!=parent.ncols():
        raise TypeError("a diagonalizable matrix must be square.")
    if eigenvalues is not None and dimensions is None:
        raise ValueError("the list of eigenvalues must have a list of dimensions corresponding to each eigenvalue.")
    if eigenvalues is None and dimensions is not None:
        raise ValueError("the list of dimensions must have a list of corresponding eigenvalues.")
    if eigenvalues is None and dimensions is None:
        values=[]
        #create a list with "size" number of entries
        for eigen_index in range(size):
            eigenvalue=randint(-10,10)
            values.append(eigenvalue)
        values.sort()
        dimensions=[]
        eigenvalues=[]
        #create a list with no duplicate values to be the eigenvalues
        for eigenvalue in range(size):
            if values[eigenvalue] not in eigenvalues:
                eigenvalues.append(values[eigenvalue])
        for dimension in range(len(eigenvalues)):
            #dimension is equal to how many times an eigenvalue was generated in the 'values' list
            dimensions.append(values.count(eigenvalues[dimension]))
    size_check=0
    for check in range(len(dimensions)):
        size_check=size_check+dimensions[check]
    if not map(lambda x: x in ZZ, eigenvalues)==[True]*len(eigenvalues):
        raise TypeError("eigenvalues must be integers.")
    if size!=size_check:
        raise ValueError("the size of the matrix must equal the sum of the dimensions.")
    if min(dimensions)<1:
        raise ValueError("eigenspaces must have a dimension of at least 1.")
    if len(eigenvalues)!=len(dimensions):
        raise ValueError("each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.")
    #sort the dimensions in order of increasing size, and sort the eigenvalues list in an identical fashion, to maintain corresponding values.
    dimensions_sort=sorted(zip(dimensions,eigenvalues))
    dimensions=[x[0] for x in dimensions_sort]
    eigenvalues=[x[1] for x in dimensions_sort]
    #Create the matrix of eigenvalues on the diagonal.  Use a lower limit and upper limit determined by the eigenvalue dimensions.
    diagonal_matrix=matrix(QQ,size)
    up_bound=0
    low_bound=0
    for row_index in range(len(dimensions)):
        up_bound=up_bound+dimensions[row_index]
        for entry in range(low_bound,up_bound):
            diagonal_matrix[entry,entry]=eigenvalues[row_index]
        low_bound=low_bound+dimensions[row_index]
    # Create a matrix to hold each of the eigenvectors as its columns, begin with an identity matrix so that after row and column
    # operations the resulting matrix will be unimodular.
    eigenvector_matrix=matrix(QQ,size,size,1)
    upper_limit=0
    lower_limit=0
    #run the routine over the necessary number of columns corresponding eigenvalue dimension.
    for dimension_index in range(len(dimensions)-1):
        upper_limit=upper_limit+dimensions[dimension_index]
        lowest_index_row_with_one=size-dimensions[dimension_index]
        #assign a one to the row that is the eigenvalue dimension rows up from the bottom row then assign ones diagonally down to the right.
        for eigen_ones in range(lower_limit,upper_limit):
            eigenvector_matrix[lowest_index_row_with_one,eigen_ones]=1
            lowest_index_row_with_one+=1
        lower_limit=lower_limit+dimensions[dimension_index]
    #Create a list to give the eigenvalue dimension corresponding to each column.
    dimension_check=[]
    for i in range(len(dimensions)):
        for k in range(dimensions[i]):
            dimension_check.append(dimensions[i])
    #run routine over the rows that are in the range of the protected ones.  Use addition of column multiples to fill entries.
    for dimension_multiplicity in range(max(dimensions),min(dimensions),-1):
        highest_one_row=size-dimension_multiplicity
        highest_one_column=0
        #find the column with the protected one in the lowest indexed row.
        while eigenvector_matrix[highest_one_row,highest_one_column]==0:
            highest_one_column+=1
        #dimension_check determines if column has a low enough eigenvalue dimension to take a column multiple.
        for bottom_entry_filler in range(len(dimension_check)):
            if dimension_check[bottom_entry_filler]<dimension_multiplicity and eigenvector_matrix[highest_one_row,bottom_entry_filler]==0:
                # randint range determined experimentally to keep entries manageable.
                eigenvector_matrix.add_multiple_of_column(bottom_entry_filler,highest_one_column,randint(-4,4))
    #Fill remaining rows using scalar row addition.
    for row in range(size-max(dimensions),size):
        for upper_row in range(size-max(dimensions)):
            # range of multiplier determined experimentally so that entries stay manageable for small matrices
            eigenvector_matrix.add_multiple_of_row(upper_row,row,randint(-4,4))
    return eigenvector_matrix*diagonal_matrix*(eigenvector_matrix.inverse())


@matrix_method
def vector_on_axis_rotation_matrix(v, i, ring=None):
    r"""
    Return a rotation matrix `M` such that `det(M)=1` sending the vector
    `v` on the i-th axis so that all other coordinates of `Mv` are zero.

    .. NOTE::

        Such a matrix is not uniquely determined. This function returns one
        such matrix.

    INPUT:

    - ``v``` - vector
    - ``i`` - integer
    - ``ring`` - ring (optional, default: None) of the resulting matrix

    OUTPUT:

    A matrix

    EXAMPLES::

        sage: from sage.matrix.constructor import vector_on_axis_rotation_matrix
        sage: v = vector((1,2,3))
        sage: vector_on_axis_rotation_matrix(v, 2) * v
        (0, 0, sqrt(14))
        sage: vector_on_axis_rotation_matrix(v, 1) * v
        (0, sqrt(14), 0)
        sage: vector_on_axis_rotation_matrix(v, 0) * v
        (sqrt(14), 0, 0)

    ::

        sage: x,y = var('x,y')
        sage: v = vector((x,y))
        sage: vector_on_axis_rotation_matrix(v, 1)
        [ y/sqrt(x^2 + y^2) -x/sqrt(x^2 + y^2)]
        [ x/sqrt(x^2 + y^2)  y/sqrt(x^2 + y^2)]
        sage: vector_on_axis_rotation_matrix(v, 0)
        [ x/sqrt(x^2 + y^2)  y/sqrt(x^2 + y^2)]
        [-y/sqrt(x^2 + y^2)  x/sqrt(x^2 + y^2)]
        sage: vector_on_axis_rotation_matrix(v, 0) * v
        (x^2/sqrt(x^2 + y^2) + y^2/sqrt(x^2 + y^2), 0)
        sage: vector_on_axis_rotation_matrix(v, 1) * v
        (0, x^2/sqrt(x^2 + y^2) + y^2/sqrt(x^2 + y^2))

    ::

        sage: v = vector((1,2,3,4))
        sage: vector_on_axis_rotation_matrix(v, 0) * v
        (sqrt(30), 0, 0, 0)
        sage: vector_on_axis_rotation_matrix(v, 0, ring=RealField(10))
        [ 0.18  0.37  0.55  0.73]
        [-0.98 0.068  0.10  0.14]
        [ 0.00 -0.93  0.22  0.30]
        [ 0.00  0.00 -0.80  0.60]
        sage: vector_on_axis_rotation_matrix(v, 0, ring=RealField(10)) * v
        (5.5, 0.00098, 0.00098, 0.00)

    AUTHORS:

        Sebastien Labbe (April 2010)
    """
    dim = len(v)
    v = vector(v)
    m = identity_matrix(dim, sparse=True)
    L = range(i-1, -1, -1) + range(dim-1,i,-1)
    for i in L:
        rot = ith_to_zero_rotation_matrix(v, i, ring=ring)
        v = rot * v
        m = rot * m
    return m


@matrix_method
def ith_to_zero_rotation_matrix(v, i, ring=None):
    r"""
    Return a rotation matrix that sends the i-th coordinates of the
    vector v to zero by doing a rotation with the (i-1)-th coordinate.

    INPUT:

    - ``v``` - vector
    - ``i`` - integer
    - ``ring`` - ring (optional, default: None) of the resulting matrix

    OUTPUT:

    A matrix

    EXAMPLES::

        sage: from sage.matrix.constructor import ith_to_zero_rotation_matrix
        sage: v = vector((1,2,3))
        sage: ith_to_zero_rotation_matrix(v, 2)
        [             1              0              0]
        [             0  2/13*sqrt(13)  3/13*sqrt(13)]
        [             0 -3/13*sqrt(13)  2/13*sqrt(13)]
        sage: ith_to_zero_rotation_matrix(v, 2) * v
        (1, sqrt(13), 0)

    ::

        sage: ith_to_zero_rotation_matrix(v, 0)
        [ 3/10*sqrt(10)              0 -1/10*sqrt(10)]
        [             0              1              0]
        [ 1/10*sqrt(10)              0  3/10*sqrt(10)]
        sage: ith_to_zero_rotation_matrix(v, 1)
        [ 1/5*sqrt(5)  2/5*sqrt(5)            0]
        [-2/5*sqrt(5)  1/5*sqrt(5)            0]
        [           0            0            1]
        sage: ith_to_zero_rotation_matrix(v, 2)
        [             1              0              0]
        [             0  2/13*sqrt(13)  3/13*sqrt(13)]
        [             0 -3/13*sqrt(13)  2/13*sqrt(13)]

    ::

        sage: ith_to_zero_rotation_matrix(v, 0) * v
        (0, 2, sqrt(10))
        sage: ith_to_zero_rotation_matrix(v, 1) * v
        (sqrt(5), 0, 3)
        sage: ith_to_zero_rotation_matrix(v, 2) * v
        (1, sqrt(13), 0)

    Other ring::

        sage: ith_to_zero_rotation_matrix(v, 2, ring=RR)
        [  1.00000000000000  0.000000000000000  0.000000000000000]
        [ 0.000000000000000  0.554700196225229  0.832050294337844]
        [ 0.000000000000000 -0.832050294337844  0.554700196225229]
        sage: ith_to_zero_rotation_matrix(v, 2, ring=RDF)
        [            1.0             0.0             0.0]
        [            0.0  0.554700196225  0.832050294338]
        [            0.0 -0.832050294338  0.554700196225]

    On the symbolic ring::

        sage: x,y,z = var('x,y,z')
        sage: v = vector((x,y,z))
        sage: ith_to_zero_rotation_matrix(v, 2)
        [                 1                  0                  0]
        [                 0  y/sqrt(y^2 + z^2)  z/sqrt(y^2 + z^2)]
        [                 0 -z/sqrt(y^2 + z^2)  y/sqrt(y^2 + z^2)]
        sage: ith_to_zero_rotation_matrix(v, 2) * v
        (x, y^2/sqrt(y^2 + z^2) + z^2/sqrt(y^2 + z^2), 0)

    TESTS::

        sage: ith_to_zero_rotation_matrix((1,0,0), 0)
        [ 0  0 -1]
        [ 0  1  0]
        [ 1  0  0]
        sage: ith_to_zero_rotation_matrix((1,0,0), 1)
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: ith_to_zero_rotation_matrix((1,0,0), 2)
        [1 0 0]
        [0 1 0]
        [0 0 1]

    AUTHORS:

        Sebastien Labbe (April 2010)
    """
    if not ring is None:
        # coerce the vector so that computations
        # are done in that ring
        v = vector(ring, v)
    dim = len(v)
    i = i % dim
    j = (i-1) % dim
    a, b = v[j], v[i]
    if b == 0:
        return identity_matrix(dim, sparse=True)
    from sage.functions.all import sqrt
    norm = sqrt(a*a + b*b)
    aa = a / norm
    bb = b / norm
    entries = {}
    for k in range(dim):
        entries[(k, k)] = 1
    entries.update({(j,j):aa, (j,i):bb, (i,j):-bb, (i,i):aa})
    return matrix(entries, nrows=dim, ring=ring)



