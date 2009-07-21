r"""
Matrices over an arbitrary ring

AUTHORS:

- William Stein

- Martin Albrecht: conversion to Pyrex

- Jaap Spies: various functions

- Gary Zablackis: fixed a sign bug in generic determinant.

- William Stein and Robert Bradshaw - complete restructuring.

- Rob Beezer - refactor kernel functions.

Elements of matrix spaces are of class ``Matrix`` (or a
class derived from Matrix). They can be either sparse or dense, and
can be defined over any base ring.

EXAMPLES:

We create the `2\times 3` matrix

.. math::

        \left(\begin{matrix} 1&2&3\\4&5&6 \end{matrix}\right)


as an element of a matrix space over `\QQ`::

    sage: M = MatrixSpace(QQ,2,3)
    sage: A = M([1,2,3, 4,5,6]); A
    [1 2 3]
    [4 5 6]
    sage: A.parent()
    Full MatrixSpace of 2 by 3 dense matrices over Rational Field

Alternatively, we could create A more directly as follows (which
would completely avoid having to create the matrix space)::

    sage: A = matrix(QQ, 2, [1,2,3, 4,5,6]); A
    [1 2 3]
    [4 5 6]

We next change the top-right entry of `A`. Note that matrix
indexing is `0`-based in Sage, so the top right entry is
`(0,2)`, which should be thought of as "row number
`0`, column number 2".

::

    sage: A[0,2] = 389
    sage: A
    [  1   2 389]
    [  4   5   6]

Also notice how matrices print. All columns have the same width and
entries in a given column are right justified. Next we compute the
reduced row echelon form of `A`.

::

    sage: A.echelon_form()
    [      1       0 -1933/3]
    [      0       1  1550/3]

We save and load a matrix::

    sage: A = matrix(Integers(8),3,range(9))
    sage: loads(dumps(A)) == A
    True

MUTABILITY: Matrices are either immutable or not. When initially
created, matrices are typically mutable, so one can change their
entries. Once a matrix `A` is made immutable using
``A.set_immutable()`` the entries of `A`
cannot be changed, and `A` can never be made mutable again.
However, properties of `A` such as its rank, characteristic
polynomial, etc., are all cached so computations involving
`A` may be more efficient. Once `A` is made
immutable it cannot be changed back. However, one can obtain a
mutable copy of `A` using ``A.copy()``.

EXAMPLES::

    sage: A = matrix(RR,2,[1,10,3.5,2])
    sage: A.set_immutable()
    sage: A.copy() is A
    False

The echelon form method always returns immutable matrices with
known rank.

EXAMPLES::

    sage: A = matrix(Integers(8),3,range(9))
    sage: A.determinant()
    0
    sage: A[0,0] = 5
    sage: A.determinant()
    1
    sage: A.set_immutable()
    sage: A[0,0] = 5
    Traceback (most recent call last):
    ...
    ValueError: matrix is immutable; please change a copy instead (use self.copy()).

Implementation and Design
-------------------------

Class Diagram (an x means that class is currently supported)::

    x Matrix
    x   Matrix_sparse
    x     Matrix_generic_sparse
    x     Matrix_integer_sparse
    x     Matrix_rational_sparse
          Matrix_cyclo_sparse
    x     Matrix_modn_sparse
          Matrix_RR_sparse
          Matrix_CC_sparse
          Matrix_RDF_sparse
          Matrix_CDF_sparse

    x  Matrix_dense
    x     Matrix_generic_dense
    x     Matrix_integer_dense
          Matrix_integer_2x2_dense
    x     Matrix_rational_dense
          Matrix_cyclo_dense    -- idea: restrict scalars to QQ, compute charpoly there, then factor
    x     Matrix_modn_dense
          Matrix_RR_dense
          Matrix_CC_dense
    x     Matrix_real_double_dense
    x     Matrix_complex_double_dense

The corresponding files in the sage/matrix library code directory
are named

::

              [matrix] [base ring] [dense or sparse].

See the files ``matrix_template.pxd`` and
``matrix_template.pyx``.

::

    New matrices types can only be implemented in Cython.

    *********** LEVEL 1  **********
    NON-OPTIONAL
    For each base field it is *absolutely* essential to completely
    implement the following functionality for that base ring:

       * __cinit__     -- should use sage_malloc from ext/stdsage.pxi (only
                          needed if allocate memory)
       * __init__      -- this signature: 'def __init__(self, parent, entries, copy, coerce)'
       * __dealloc__   -- use sage_free (only needed if allocate memory)
       * set_unsafe(self, size_t i, size_t j, x) -- doesn't do bounds or any other checks; assumes x is in self._base_ring
       * get_unsafe(self, size_t i, size_t j) -- doesn't do checks
       * __richcmp__    -- always the same (I don't know why its needed -- bug in PYREX).
    *********** LEVEL 2  **********

    IMPORTANT (and *highly* recommended):

    After getting the special class with all level 1 functionality to
    work, implement all of the following (they should not change
    functionality, except speed (always faster!) in any way):

       * def _pickle(self):
              return data, version
       * def _unpickle(self, data, int version)
              reconstruct matrix from given data and version; may assume _parent, _nrows, and _ncols are set.
              Use version numbers >= 0 so if you change the pickle strategy then
              old objects still unpickle.
       * cdef _list -- list of underlying elements (need not be a copy)
       * cdef _dict -- sparse dictionary of underlying elements
       * cdef _add_ -- add two matrices with identical parents
       * _matrix_times_matrix_c_impl -- multiply two matrices with compatible dimensions and
                                        identical base rings (both sparse or both dense)
       * cdef _cmp_c_impl -- compare two matrices with identical parents
       * cdef _lmul_c_impl -- multiply this matrix on the right by a scalar, i.e., self * scalar
       * cdef _rmul_c_impl -- multiply this matrix on the left by a scalar, i.e., scalar * self
       * __copy__
       * __neg__

    The list and dict returned by _list and _dict will *not* be changed
    by any internal algorithms and are not accessible to the user.

    *********** LEVEL 3  **********
    OPTIONAL:

       * cdef _sub_
       * __invert__
       * _multiply_classical
       * __deepcopy__

    Further special support:
       * Matrix windows -- to support Strassen multiplication for a given base ring.
       * Other functions, e.g., transpose, for which knowing the
         specific representation can be helpful.

    .. note::

       - For caching, use self.fetch and self.cache.

       - Any method that can change the matrix should call
         ``check_mutability()`` first.  There are also many fast cdef'd bounds checking methods.

       - Kernels of matrices
         Implement only a left_kernel() or right_kernel() method, whichever requires
         the least overhead (usually meaning little or no transposing).  Let the
         methods in the matrix2 class handle left, right, generic kernel distinctions.
"""
