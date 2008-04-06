"""
Fast Fourier Transforms using GSL.


AUTHORS:
    William Stein (2006-9) - initial file (radix2)
    D. Joyner (2006-10) - Minor modifications (from radix2 to general case
                          and some documentation).

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
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

include '../ext/stdsage.pxi'

import sage.plot.all
import sage.libs.pari.all
from sage.rings.integer import Integer
from sage.rings.complex_number import ComplexNumber

def FastFourierTransform(size, base_ring=None):
    """
    EXAMPLES:
        sage: a = FastFourierTransform(128)
        sage: for i in range(1, 11):
        ...    a[i] = 1
        ...    a[128-i] = 1
        sage: a.plot().show(ymin=0)
        sage: a.forward_transform()
        sage: a.plot().show()
    """
    return FastFourierTransform_complex(int(size))

FFT = FastFourierTransform

cdef class FastFourierTransform_base:
    pass

cdef class FastFourierTransform_complex(FastFourierTransform_base):

    def __init__(self, size_t n, size_t stride=1):
        self.n = n
        self.stride = stride
        self.data = <double*> sage_malloc(sizeof(double)*(2*n))
        cdef int i
        for i from 0 <= i < 2*n:
            self.data[i] = 0

    def  __dealloc__(self):
        sage_free(self.data)

    def __len__(self):
        return self.n

    def __setitem__(self, size_t i, xy):
        # just set real for now
        if i < 0 or i >= self.n:
            raise IndexError
        if isinstance(xy, (tuple, ComplexNumber)):
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

    def __getslice__(self, Py_ssize_t i, Py_ssize_t j):
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
        using the Cooley-Tukey algorithm. If the number of sample
        points in the input is a power of 2 then the function
        gsl_fft_complex_radix2_forward is automatically called.
        Otherwise, gsl_fft_complex_forward is called.
        """
        cdef gsl_fft_complex_wavetable * wt
        cdef gsl_fft_complex_workspace * mem
        N = Integer(self.n)
        e = N.exact_log(2)
        if N==2**e:
            gsl_fft_complex_radix2_forward(self.data, self.stride, self.n)
        if N!=2**e:
            #cdef gsl_fft_complex_wavetable * wt
            #cdef gsl_fft_complex_workspace * mem
            mem = gsl_fft_complex_workspace_alloc(self.n)
            wt = gsl_fft_complex_wavetable_alloc(self.n)
            gsl_fft_complex_forward(self.data, self.stride, self.n, wt, mem)
            gsl_fft_complex_workspace_free(mem)
            gsl_fft_complex_wavetable_free(wt)

    def inverse_transform(self):
        """
        Compute the in-place forward Fourier transform of this data
        using the Cooley-Tukey algorithm. If the number of sample
        points in the input is a power of 2 then the function
        gsl_fft_complex_radix2_inverse is automatically called.
        Otherwise, gsl_fft_complex_inverse is called.

        EXAMPLES:
            sage: a = FastFourierTransform(125)
            sage: b = FastFourierTransform(125)
            sage: for i in range(1, 60): a[i]=1
            sage: for i in range(1, 60): b[i]=1
            sage: a.forward_transform()
            sage: a.inverse_transform()
            sage: (a.plot()+b.plot()).save('a.png', ymin=0)
        """
        cdef gsl_fft_complex_wavetable * wt
        cdef gsl_fft_complex_workspace * mem
        N = Integer(self.n)
        e = N.exact_log(2)
        if N==2**e:
            gsl_fft_complex_inverse(self.data, self.stride, self.n, wt, mem)
        if N!=2**e:
            mem = gsl_fft_complex_workspace_alloc(self.n)
            wt = gsl_fft_complex_wavetable_alloc(self.n)
            gsl_fft_complex_inverse(self.data, self.stride, self.n, wt, mem)
            gsl_fft_complex_workspace_free(mem)
            gsl_fft_complex_wavetable_free(wt)

    def backward_transform(self):
        """
        Compute the in-place backwards Fourier transform of this data
        using the Cooley-Tukey algorithm. This is the same as "inverse"
        but lacks normalization so that backwards*forwards(f) = n*f.

        EXAMPLES:
            sage: a = FastFourierTransform(125)
            sage: b = FastFourierTransform(125)
            sage: for i in range(1, 60): a[i]=1
            sage: for i in range(1, 60): b[i]=1
            sage: a.forward_transform()
            sage: a.backward_transform()
            sage: (a.plot()+b.plot()).save('a.png', ymin=0)

        """
        cdef gsl_fft_complex_wavetable * wt
        cdef gsl_fft_complex_workspace * mem
        N = Integer(self.n)
        e = N.exact_log(2)
        if N==2**e:
            gsl_fft_complex_backward(self.data, self.stride, self.n, wt, mem)
        if N!=2**e:
            mem = gsl_fft_complex_workspace_alloc(self.n)
            wt = gsl_fft_complex_wavetable_alloc(self.n)
            gsl_fft_complex_backward(self.data, self.stride, self.n, wt, mem)
            gsl_fft_complex_workspace_free(mem)
            gsl_fft_complex_wavetable_free(wt)

cdef class FourierTransform_complex:
    pass

cdef class FourierTransform_real:
    pass


