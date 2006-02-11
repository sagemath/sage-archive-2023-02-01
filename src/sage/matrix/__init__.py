"""
Classes for matrix algebra.

The most important classes are matrix and matrix_space.

   matrix_space -- top level definition of MatrixSpace class
   matrix -- definition of Matrix class

The other modules define sometimes strange implementations of
various aspects of matrix arithmetic.  The matrix and matrix_space
modules define a CLEAN UNIFIED interface to all this other
functionality.

The other modules are as follows:
   dense_matrix_pyx -- dense matrices over Z/pZ and Q in pyrex
   sparse_matrix -- sparse matrices over Z/pZ and Q
   sparse_matrix_pyx -- sparse matrices over Z/pZ and Q in pyrex
"""
