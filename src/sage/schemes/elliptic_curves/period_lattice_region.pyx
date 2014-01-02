r"""
This module is used to represent sub-regions of a fundamental parallelogram
of the period lattice of an elliptic curve, used in computing minimum height
bounds.
"""

##############################################################################
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
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
##############################################################################


import numpy as np
cimport numpy as np

from sage.rings.all import CIF
from sage.rings.arith import gcd

cdef class PeriodicRegion:

    # The generators of the lattice, a complex numbers.
    cdef public w1
    cdef public w2

    # The bitmap of this region.
    cdef public data

    # Whether bitmap represents the whole fundamental parallelogram,
    # or only half (in which case it is symmetric about the center).
    cdef public bint full

    def __init__(self, w1, w2, data, full=True):
        if data.dtype is not np.int8:
            data = data.astype(np.int8)
        full = int(full)
        self.w1 = w1
        self.w2 = w2
        self.data = data
        self.full = full

    def is_empty(self):
        return self.data.sum() == 0

    def _ensure_full(self):
        if not self.full:
            rows, cols = self.data.shape
            new_data = np.ndarray((rows, 2*cols), self.data.dtype)
            new_data[:,0:cols] = self.data
            new_data[::-1,-1:cols-1:-1] = self.data
            self.data = new_data
            self.full = True
        return self

    def ds(self):
        m, n = self.data.shape
        return self.w1/m, self.w2/(n * (2-self.full))

    def verify(self, condition):
        for line in self.border(False):
            if not condition(line):
                return False
        return True

    def refine(self, condition=None, int times=1):
        if times <= 0:
            return self
        cdef int i, j, m, n
        cdef np.ndarray[np.npy_int8, ndim=2] less, fuzz, new_data
        m, n = self.data.shape
        if condition is None:
            new_data = np.ndarray((2*m, 2*n), self.data.dtype)
            new_data[0::2,0::2] = self.data
            new_data[1::2,0::2] = self.data
            new_data[0::2,1::2] = self.data
            new_data[1::2,1::2] = self.data
        else:
            more = self.expand().data
            less = self.contract().data
            fuzz = less ^ more
            new_data = np.zeros((2*m, 2*n), self.data.dtype)
            dw1, dw2 = self.ds()
            dw1 /= 2
            dw2 /= 2
            for i in range(2*m):
                for j in range(2*n):
                    if less[i/2, j/2]:
                        new_data[i,j] = True
                    elif fuzz[i/2, j/2]:
                        new_data[i,j] = condition(dw1*(i+.5) + dw2*(j+.5))
        return PeriodicRegion(self.w1, self.w2, new_data, self.full).refine(condition, times-1)

    def expand(self, bint corners=True):
        cdef int i, j, m, n
        m, n = self.data.shape
        cdef np.ndarray[np.npy_int8, ndim=2] framed, new_data
        framed = frame_data(self.data, self.full)
        new_data = np.zeros((m+2, n+2), self.data.dtype)
        for i in range(m):
            for j in range(n):
                if framed[i,j]:
                    new_data[i  , j  ] = True
                    new_data[i-1, j  ] = True
                    new_data[i+1, j  ] = True
                    new_data[i  , j-1] = True
                    new_data[i  , j+1] = True
                    if corners:
                        new_data[i-1, j-1] = True
                        new_data[i+1, j-1] = True
                        new_data[i+1, j+1] = True
                        new_data[i-1, j+1] = True
        return PeriodicRegion(self.w1, self.w2, unframe_data(new_data, self.full), self.full)

    def contract(self, corners=True):
        return ~(~self).expand(corners)

    def __div__(self, int n):
        cdef int i, j, a, b, rows, cols, g
        cdef double d, e
        if n == 1:
            return self
        self._ensure_full()
        cdef np.ndarray[np.npy_int8, ndim=2] data, new_data
        data = self.data
        rows, cols = self.data.shape
        g = gcd([n, rows, cols])
        if g != 1:
            new_data = np.zeros(self.data.shape, self.data.dtype)
            for a in range(g):
                for b in range(g):
                    new_data[(rows/g)*a:(rows/g)*(a+1), (cols/g)*b:(cols/g)*(b+1)] |= data[a::g, b::g]
            n /= g
            data = new_data
        if n != 1:
            new_data = np.zeros(self.data.shape, self.data.dtype)
            d = 1.0/n
            e = 1-d/2
            for i in range(rows):
                for j in range(cols):
                    if data[i,j]:
                        for a in range(n):
                            for b in range(n):
                                new_data[<int>((a*rows+i  )*d), <int>((b*cols+j  )*d)] = True
                                new_data[<int>((a*rows+i+e)*d), <int>((b*cols+j  )*d)] = True
                                new_data[<int>((a*rows+i+e)*d), <int>((b*cols+j+e)*d)] = True
                                new_data[<int>((a*rows+i  )*d), <int>((b*cols+j+e)*d)] = True
            data = new_data
        return PeriodicRegion(self.w1, self.w2, data)

    def __invert__(self):
        return PeriodicRegion(self.w1, self.w2, 1-self.data, self.full)

    def __and__(left, right):
        assert isinstance(left, PeriodicRegion) and isinstance(right, PeriodicRegion)
        assert left.w1 == right.w1 and left.w2 == right.w2
        if left.full ^ right.full:
            left._ensure_full()
            right._ensure_full()
        return PeriodicRegion(left.w1, left.w2, left.data & right.data, left.full)

    def __xor__(left, right):
        assert isinstance(left, PeriodicRegion) and isinstance(right, PeriodicRegion)
        assert left.w1 == right.w1 and left.w2 == right.w2
        if left.full ^ right.full:
            left._ensure_full()
            right._ensure_full()
        return PeriodicRegion(left.w1, left.w2, left.data ^ right.data, left.full)

    def border(self, raw=True):
        cdef np.ndarray[np.npy_int8, ndim=2] framed = frame_data(self.data, self.full)
        cdef int m, n
        m, n = self.data.shape
        cdef int i, j
        L = []
        for i in range(m):
            for j in range(n):
                if framed[i,j]:
                    if not framed[i-1, j]:
                        L.append((i, j, 0))
                    if not framed[i+1, j]:
                        L.append((i+1, j, 0))
                    if not framed[i, j-1]:
                        L.append((i, j, 1))
                    if not framed[i, j+1]:
                        L.append((i, j+1, 1))
        if not raw:
            dw1, dw2 = self.ds()
            for ix, (i, j, dir) in enumerate(L):
                L[ix] = CIF(i*dw1 + j*dw2).union(CIF((i+dir)*dw1 + (j+(1-dir))*dw2))
        return L

    def innermost_point(self):
        if self.is_empty():
            raise ValueError, "empty"
        inside = self
        while not inside.is_empty():
            self, inside = inside, inside.contract()
        i, j = tuple(np.transpose(self.data.nonzero())[0])
        dw1, dw2 = self.ds()
        return float(i + 0.5) * dw1 + float(j + 0.5) * dw2

    def plot(self, **kwds):
        from sage.all import line
        dw1, dw2 = self.ds()
        L = []
        F = line([(0,0), tuple(self.w1), tuple(self.w1+self.w2), tuple(self.w2), (0,0)])
        if not self.full:
            F += line([tuple(self.w2/2), tuple(self.w1+self.w2/2)])
        if 'color' not in kwds:
            kwds['color'] = 'red'
        for i, j, dir in self.border():
            ii, jj = i+dir, j+(1-dir)
            L.append(line([tuple(i*dw1 + j*dw2), tuple(ii*dw1 + jj*dw2 )], **kwds))
        return sum(L, F)


def frame_data(data, full=True):
    m, n = data.shape
    framed = np.empty((m+2, n+2), data.dtype)
    # center
    framed[:-2,:-2] = data
    # top and bottom
    if full:
        framed[:-2,-1] = data[:,-1]
        framed[:-2,-2] = data[:, 0]
    else:
        framed[:-2,-1] = data[::-1, 0]
        framed[:-2,-2] = data[::-1,-1]
    # left and right
    framed[-2,:] = framed[ 0,:]
    framed[-1,:] = framed[-3,:]
    return framed

def unframe_data(framed, full=True):
    framed = framed.copy()
    framed[ 0,:] |= framed[-2,:]
    framed[-3,:] |= framed[-1,:]
    if full:
        framed[:-2,-3] |= framed[:-2,-1]
        framed[:-2, 0] |= framed[:-2,-2]
    else:
        framed[-3::-1, 0] |= framed[:-2,-1]
        framed[-3::-1,-3] |= framed[:-2,-2]
    return framed[:-2,:-2]
