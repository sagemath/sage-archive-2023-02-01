"""
Python package for convex optimization

CVXOPT is a free software package for convex optimization based on the
Python programming language. It can be used with the interactive Python
interpreter, on the command line by executing Python scripts, or
integrated in other software via Python extension modules. Its main
purpose is to make the development of software for convex optimization
applications straightforward by building on Python's extensive standard
library and on the strengths of Python as a high-level programming
language.
"""

# Copyright 2004-2009 J. Dahl and L. Vandenberghe.
#
# This file is part of CVXOPT version 1.1.2
#
# CVXOPT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# CVXOPT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License

import base

def normal(nrows, ncols=1, mean=0.0, std=1.0):
    '''
    Randomly generates a matrix with normally distributed entries.

    normal(nrows, ncols=1, mean=0, std=1)

    PURPOSE
    Returns a matrix with typecode 'd' and size nrows by ncols, with
    its entries randomly generated from a normal distribution with mean
    m and standard deviation std.

    ARGUMENTS
    nrows     number of rows

    ncols     number of columns

    mean      approximate mean of the distribution
    std       standard deviation of the distribution
    '''

    from cvxopt import gsl
    return gsl.normal(nrows, ncols, mean, std)

def uniform(nrows, ncols=1, a=0, b=1):
    '''
    Randomly generates a matrix with uniformly distributed entries.

    uniform(nrows, ncols=1, a=0, b=1)

    PURPOSE
    Returns a matrix with typecode 'd' and size nrows by ncols, with
    its entries randomly generated from a uniform distribution on the
    interval (a,b).

    ARGUMENTS
    nrows     number of rows

    ncols     number of columns

    a         lower bound

    b         upper bound
    '''

    from cvxopt import gsl
    return gsl.uniform(nrows, ncols, a, b)

def setseed(val = 0):
    '''
    Sets the seed value for the random number generator.

    setseed(val = 0)

    ARGUMENTS
    value     integer seed.  If the value is 0, the current system time
              is used.
    '''

    from cvxopt import gsl
    gsl.setseed(val)


def getseed():
    '''
    Returns the seed value for the random number generator.

    getseed()
    '''

    from cvxopt import gsl
    return gsl.getseed()


import __builtin__
omax = __builtin__.max
omin = __builtin__.min

def max(*args):
    '''
    Elementwise max for matrices.

    PURPOSE
    max(a1, ..., an) computes the elementwise max for matrices.  The arguments
    must be matrices of compatible dimensions,  and scalars.  The elementwise
    max of a matrix and a scalar is a new matrix where each element is
    defined as max of a matrix element and the scalar.

    max(iterable)  where the iterator generates matrices and scalars computes
    the elementwise max between the objects in the iterator,  using the
    same conventions as max(a1, ..., an).
    '''

    if len(args) == 1 and type(args[0]).__name__ in \
            ['list', 'tuple', 'xrange', 'generator']:
        return +reduce(base.emax, *args)
    elif len(args) == 1 and type(args[0]) is base.matrix:
        return omax(args[0])
    elif len(args) == 1 and type(args[0]) is base.spmatrix:
        if len(args[0]) == mul(args[0].size):
            return omax(args[0])
        else:
            return omax(omax(args[0]), 0.0)
    else:
        return +reduce(base.emax, args)


def min(*args):
    '''
    Elementwise min for matrices.

    PURPOSE
    min(a1, ..., an) computes the elementwise min for matrices.  The arguments
    must be matrices of compatible dimensions,  and scalars.  The elementwise
    min of a matrix and a scalar is a new matrix where each element is
    defined as min of a matrix element and the scalar.

    min(iterable)  where the iterator generates matrices and scalars computes
    the elementwise min between the objects in the iterator,  using the
    same conventions as min(a1, ..., an).
    '''

    if len(args) == 1 and type(args[0]).__name__ in \
            ['list', 'tuple', 'xrange', 'generator']:
        return +reduce(base.emin, *args)
    elif len(args) == 1 and type(args[0]) is base.matrix:
        return omin(args[0])
    elif len(args) == 1 and type(args[0]) is base.spmatrix:
        if len(args[0]) == mul(args[0].size):
            return omin(args[0])
        else:
            return omin(omin(args[0]), 0.0)
    else:
        return +reduce(base.emin, args)

def mul(*args):
    '''
    Elementwise multiplication for matrices.

    PURPOSE
    mul(a1, ..., an) computes the elementwise multiplication for matrices.
    The arguments must be matrices of compatible dimensions,  and scalars.
    The elementwise multiplication of a matrix and a scalar is a new matrix
    where each element is
    defined as the multiplication of a matrix element and the scalar.

    mul(iterable)  where the iterator generates matrices and scalars computes
    the elementwise multiplication between the objects in the iterator,
    using the same conventions as mul(a1, ..., an).
    '''

    if len(args) == 1 and type(args[0]).__name__ in \
            ['list', 'tuple', 'xrange', 'generator']:
        return +reduce(base.emul, *args)
    else:
        return +reduce(base.emul, args)

def div(*args):
    '''
    Elementwise division for matrices.

    PURPOSE
    div(a1, ..., an) computes the elementwise division for matrices.
    The arguments must be matrices of compatible dimensions,  and scalars.
    The elementwise division of a matrix and a scalar is a new matrix
    where each element is defined as the division between a matrix element
    and the scalar.

    div(iterable)  where the iterator generates matrices and scalars computes
    the elementwise division between the objects in the iterator,
    using the same conventions as div(a1, ..., an).
    '''

    if len(args) == 1 and type(args[0]).__name__ in \
            ['list', 'tuple', 'xrange', 'generator']:
        return +reduce(base.ediv, *args)
    else:
        return +reduce(base.ediv, args)

base.normal, base.uniform = normal, uniform
base.setseed, base.getseed = setseed, getseed
base.mul, base.div = mul, div

from cvxopt import printing
matrix_str    = printing.matrix_str_default
matrix_repr   = printing.matrix_repr_default
spmatrix_str  = printing.spmatrix_str_default
spmatrix_repr = printing.spmatrix_repr_default

from cvxopt.base import matrix, spmatrix, sparse, spdiag, sqrt, sin, cos, \
    exp, log

__all__ = [ 'blas', 'lapack', 'amd', 'umfpack', 'cholmod', 'solvers',
    'modeling', 'printing', 'info', 'matrix', 'spmatrix', 'sparse',
    'spdiag', 'sqrt', 'sin', 'cos', 'exp', 'log', 'min', 'max', 'mul',
    'div', 'normal', 'uniform', 'setseed', 'getseed' ]
