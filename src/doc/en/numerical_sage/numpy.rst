NumPy
=====

NumPy is not imported into sage initially so first at your sage
prompt (or in the notebook), import numpy.

::

    sage: import numpy

The basic object of computation in NumPy is arrays. Do the
following

.. link

::

    sage: l=numpy.array([1,2,3])
    sage: l
    array([1, 2, 3], dtype=object)

Note the dtype (data type) of object. NumPy arrays can store any
type of python object. But this is almost never what you want. If
you do the following

.. link

::

    sage: l=numpy.array([1,2,3],dtype=float)
    sage: l
    array([ 1.,  2.,  3.])

Now

.. link

::

    sage: l.dtype
    dtype('float64')

should report a data-type of float. If the data type is float then
the array is stored as an array of machine floats which takes up
much less space and which can be operated on much faster. Usually
you want your NumPy arrays to be arrays of floats. Another subtlety
to watch out for is that if you do

.. link

::

    sage: l=numpy.array([1.0,2.0,3.0])
    sage: l
    array([1.00000000000000, 2.00000000000000, 3.00000000000000], dtype=object)

then l.dtype will still report a data type of object. If no data
type is explicitly shown when NumPy prints the object it is either
float or int, meaning machine floats or ints. If you want floats
you need to specify dtype=float when creating your array.
Alternatively do

.. link

::

    sage: RealNumber=float

Now

.. link

::

    sage: l=numpy.array([1.0,2.0,3.0])

will create an array of floats like you wanted without explicitly
requiring you to specify the data type.

You can access elements of a NumPy array just like any list, as
well as take slices

.. link

::

    sage: l=numpy.array(range(10),dtype=float)
    sage: l[3]
    3.0
    sage: l[3:6]
    array([ 3.,  4.,  5.])

You can do basic arithmetic operations

.. link

::

    sage: l+l
    array([  0.,   2.,   4.,   6.,   8.,  10.,  12.,  14.,  16.,  18.])
    sage: 2.5*l
    array([  0. ,   2.5,   5. ,   7.5,  10. ,  12.5,  15. ,  17.5,  20. ,  22.5])

I am assuming you have set RealNumber=float, otherwise in the
second line you would need float(2.5)\*l. Note that l\*l will
multiply the elements of l componentwise. To get a dot product do
numpy.dot

.. link

::

    sage: l*l
    array([  0.,   1.,   4.,   9.,  16.,  25.,  36.,  49.,  64.,  81.])
    sage: numpy.dot(l,l)
    285.0

We can also create two dimensional arrays

.. link

::

    sage: m = numpy.array([[1,2],[3,4]],dtype=float)
    sage: m
    array([[ 1.,  2.],
           [ 3.,  4.]])
    sage: m[1,1]
    4.0

This is basically equivalent to the following

.. link

::

    sage: m=numpy.matrix([[1,2],[3,4]],dtype=float)
    sage: m
    matrix([[ 1.,  2.],
            [ 3.,  4.]])
    sage: m[0,1]
    2.0

The difference is that with numpy.array, m is treated as just an
array of data. In particular m\*m will multiply componentwise,
however with numpy.matrix m\*m will do matrix multiplication. We
can also do matrix vector multiplication, and matrix addition

.. link

::

    sage: n = numpy.matrix([[1,2],[3,4]],dtype=float)
    sage: v = numpy.array([[1],[2]],dtype=float)
    sage: n*v
    matrix([[  5.],
            [ 11.]])
    sage: n+n
    matrix([[ 2.,  4.],
            [ 6.,  8.]])

If n were created with numpy.array, then to do matrix vector
multiplication you do numpy.dot(n,v).

All NumPy arrays have a shape attribute. This is a useful attribute
to manipulate

.. link

::

    sage: n = numpy.array(range(25),dtype=float)
    sage: n
    array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
            11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.,  20.,  21.,
            22.,  23.,  24.])
    sage: n.shape=(5,5)
    sage: n
    array([[  0.,   1.,   2.,   3.,   4.],
           [  5.,   6.,   7.,   8.,   9.],
           [ 10.,  11.,  12.,  13.,  14.],
           [ 15.,  16.,  17.,  18.,  19.],
           [ 20.,  21.,  22.,  23.,  24.]])

This changes the 1 dimensional array, into a 5x5 array.

NumPy arrays can be sliced as well

.. link

::

    sage: n=numpy.array(range(25),dtype=float)
    sage: n.shape=(5,5)
    sage: n[2:4,1:3]
    array([[ 11.,  12.],
           [ 16.,  17.]])

It is important to note that the sliced matrices are references to
the original

.. link

::

    sage: m=n[2:4,1:3]
    sage: m[0,0]=100
    sage: n
    array([[   0.,    1.,    2.,    3.,    4.],
           [   5.,    6.,    7.,    8.,    9.],
           [  10.,  100.,   12.,   13.,   14.],
           [  15.,   16.,   17.,   18.,   19.],
           [  20.,   21.,   22.,   23.,   24.]])

You will note that the original matrix changed. This may or may not
be what you want. If you want to change the sliced matrix without
changing the original you should make a copy

.. link

::

    m=n[2:4,1:3].copy()

Some particularly useful commands are

.. link

::

    sage: x=numpy.arange(0,2,.1,dtype=float)
    sage: x
    array([ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1. ,
            1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9])

You can see that arange creates an array of floats increasing by .1
from 0 to 2. There is a useful command r\_ that is best explained by example

.. link

::

    sage: from numpy import r_
    sage: j=numpy.complex(0,1)
    sage: RealNumber=float
    sage: Integer=int
    sage: n=r_[0.0:5.0]
    sage: n
    array([ 0.,  1.,  2.,  3.,  4.])
    sage: n=r_[0.0:5.0, [0.0]*5]
    sage: n
    array([ 0.,  1.,  2.,  3.,  4.,  0.,  0.,  0.,  0.,  0.])

r\_ provides a shorthand for constructing NumPy arrays efficiently.
Note in the above 0.0:5.0 was shorthand for 0.0,1.0,2.0,3.0,4.0.
Suppose we want to divide the interval from 0 to 5 into 10
intervals. We can do this as follows

.. link

::

    sage: r_[0.0:5.0:11*j]
    array([ 0. ,  0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ])

The notation ``0.0:5.0:11\*j`` expands to a list of 11 equally space
points between 0 and 5 including both endpoints. Note that j is the
NumPy imaginary number, but it has this special syntax for creating
arrays. We can combine all of these techniques

.. link

::

    sage: n=r_[0.0:5.0:11*j,int(5)*[0.0],-5.0:0.0]
    sage: n
    array([ 0. ,  0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ,
            0. ,  0. ,  0. ,  0. ,  0. , -5. , -4. , -3. , -2. , -1. ])

Another useful command is meshgrid, it produces meshed grids. As an
example suppose you want to evaluate :math:`f(x,y)=x^2+y^2` on a
an equally spaced grid with :math:`\Delta x = \Delta y = .25` for
:math:`0\le x,y\le 1`. You can do that as follows

::

    sage: import numpy
    sage: RealNumber=float
    sage: j=numpy.complex(0,1)
    sage: def f(x,y):
    ...       return x**2+y**2
    sage: from numpy import meshgrid
    sage: x=numpy.r_[0.0:1.0:5*j]
    sage: y=numpy.r_[0.0:1.0:5*j]
    sage: xx,yy= meshgrid(x,y)
    sage: xx
    array([[ 0.  ,  0.25,  0.5 ,  0.75,  1.  ],
           [ 0.  ,  0.25,  0.5 ,  0.75,  1.  ],
           [ 0.  ,  0.25,  0.5 ,  0.75,  1.  ],
           [ 0.  ,  0.25,  0.5 ,  0.75,  1.  ],
           [ 0.  ,  0.25,  0.5 ,  0.75,  1.  ]])
    sage: yy
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.25],
           [ 0.5 ,  0.5 ,  0.5 ,  0.5 ,  0.5 ],
           [ 0.75,  0.75,  0.75,  0.75,  0.75],
           [ 1.  ,  1.  ,  1.  ,  1.  ,  1.  ]])
    sage: f(xx,yy)
    array([[ 0.    ,  0.0625,  0.25  ,  0.5625,  1.    ],
           [ 0.0625,  0.125 ,  0.3125,  0.625 ,  1.0625],
           [ 0.25  ,  0.3125,  0.5   ,  0.8125,  1.25  ],
           [ 0.5625,  0.625 ,  0.8125,  1.125 ,  1.5625],
           [ 1.    ,  1.0625,  1.25  ,  1.5625,  2.    ]])

You can see that meshgrid produces a pair of matrices, here denoted
:math:`xx` and :math:`yy`, such that
:math:`(xx[i,j],yy[i,j])` has coordinates :math:`(x[i],y[j])`.
This is useful because to evaluate :math:`f` over a grid, we only
need to evaluate it on each pair of entries in :math:`xx`,
:math:`yy`. Since NumPy automatically performs arithmetic
operations on arrays componentwise, it is very easy to evaluate
functions over a grid with very little code.

A useful module is the linalg module. If you want to solve an
equation ax=b do

::

    sage: import numpy
    sage: from numpy import linalg
    sage: A=numpy.random.randn(5,5)
    sage: b=numpy.array(range(1,6),dtype=float)
    sage: x=linalg.solve(A,b)
    sage: numpy.dot(A,x)
    array([ 1.,  2.,  3.,  4., 5.])

This creates a random 5x5 matrix A, and solves Ax=b where
b=[0.0,1.0,2.0,3.0,4.0]. There are many other routines in the
linalg module that are mostly self explanatory. For example there
are routines for doing qr, and lu decompositions named qr, and lu.
There is also a command eigs for computing eigenvalues of a matrix.
You can always do :math:`\langle` function
name :math:`\rangle`? to get the documentation which is quite
good for these routines.

Hopefully this gives you a sense of what NumPy is like. You should
explore the package as there is quite a bit more functionality.
