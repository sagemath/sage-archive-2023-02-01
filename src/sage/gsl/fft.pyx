"""
Fast Fourier Transforms
"""

import sage.plot.all
import sage.libs.pari.all

def FastFourierTransform(size, base_ring=None):
    """
    EXAMPLES:
        sage: a = FastFourierTransform(128)
        sage: for i in range(1, 11):
        ...    a[i] = 1
        ...    a[128-i] = 1
        sage.: show(plot(a), ymin=0)
        sage: a.forward_transform()
        sage.: show(plot(a))
    """
    return FastFourierTransform_radix2_complex(size)

FFT = FastFourierTransform

cdef class FastFourierTransform_base:
    pass

cdef class FastFourierTransform_radix2_complex(FastFourierTransform_base):

    def __init__(self, size_t n, size_t stride=1):
        self.n = n
        self.stride = stride
        self.data = <double*> PyMem_Malloc(sizeof(double)*(2*n))
        cdef int i
        for i from 0 <= i < 2*n:
            self.data[i] = 0

    def  __dealloc__(self):
        PyMem_Free(self.data)

    def __len__(self):
        return self.n

    def __setitem__(self, size_t i, xy):
        # just set real for now
        if i < 0 or i >= self.n:
            raise IndexError
        if isinstance(xy, tuple):
            self.data[2*i] = xy[0]
            self.data[2*i+1] = xy[1]
        else:
            self.data[2*i] = xy

    def __getitem__(self, size_t i):
        if i < 0 or i >= self.n:
            raise IndexError
        cdef int j
        j = 2*i
        return self.data[2*i], self.data[2*i+1]

    def __getslice__(self, i, j):
        # Todo -- make this actually fast.
        return list(self)[i:j]

    def __repr__(self):
        return str(list(self))

    def _plot_polar(self, xmin, xmax, **args):
        cdef int i
        v = []

        pari = sage.libs.pari.all.pari
        point = sage.plot.all.point
        pi = pari('Pi')
        s = 1/(3*pi)   # so arg gets scaled between -1/3 and 1/3.

        for i from xmin <= i < xmax:
            z = pari('%s + I*%s'%(self.data[2*i], self.data[2*i+1]))
            mag = z.abs()
            if mag > 0:
                arg = z.arg()*s
            else:
                arg = 0
            if i > 0:
                v.append(point([(i-1, prev_mag), (i,mag)], hue=arg, **args))
            prev_mag = mag
        return sum(v)

    def _plot_rect(self, xmin, xmax, **args):
        cdef int i
        cdef double pr_x, x, h
        v = []

        point = sage.plot.all.point

        for i from xmin <= i < xmax:
            x = self.data[2*i]
            h = self.data[2*i+1]
            if i > 0:
                v.append(point([(i-1, pr_x), (i,x)], hue=h, **args))
            pr_x = x
        return sum(v)

    def plot(self, style='rect', xmin=None, xmax=None, **args):
        """
        INPUT:
            style -- 'rect': height represents real part, color represents imaginary part
                  -- 'polar': height represents absolute value
            **args -- passed on to the line plotting function.
        """
        if xmin is None:
            xmin = 0
        else:
            xmin = int(xmin)
        if xmax is None:
            xmax = self.n
        else:
            xmax = int(xmax)
        if style == 'rect':
            return self._plot_rect(xmin, xmax, **args)
        elif style == 'polar':
            return self._plot_polar(xmin, xmax, **args)
        else:
            raise ValueError, "unknown style '%s'"%style

    def forward_transform(self):
        """
        Compute the in-place forward Fourier transform of this data
        using the Cooley-Tukey algorithm.
        """
        gsl_fft_complex_radix2_forward(self.data, self.stride, self.n)






cdef class FourierTransform_complex:
    pass

cdef class FourierTransform_radix2_real:
    pass

cdef class FourierTransform_real:
    pass
