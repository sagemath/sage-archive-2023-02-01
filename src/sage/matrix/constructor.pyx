"""
General matrix Constructor
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import types
import sage.rings.all as rings
from sage.rings.ring import is_Ring
import sage.matrix.matrix_space as matrix_space
from sage.structure.coerce import py_scalar_parent
from sage.structure.element import is_Vector
from sage.structure.sequence import Sequence
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.all import ZZ


class MatrixFactory(object):
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
    def __call__(self, *args, **kwds):
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

    TESTS:

    Check that :trac:`19920` is fixed::

        sage: import numpy
        sage: matrix([[numpy.int8(1)]])
        [1]
    """
    if not w:
        return Sequence([], rings.ZZ), rings.ZZ
    entries = Sequence(w)
    ring = entries.universe()
    if isinstance(ring,type):
        ring = py_scalar_parent(ring)
    if not is_Ring(ring):
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


from .special import *
