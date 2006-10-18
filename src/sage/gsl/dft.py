""" 10-15-2006
This file contains functions useful for computing discrete Fourier
transforms and probability distribution functions for discrete random
variables for sequences of elements from QQ or CC, indexed by
a range(..) or ZZ/NZZ or an AbelianGroup or the conjugacy classes
of a permutation group or the conjugacy classes of a matrix group.

This file implements:

*  __eq__
*  __rmul__ (for right multiplication by a scalar)
*  plotting (plot, plot_histogram)
*  dft  -  computes the discrete Fourier transform for the
           following cases:
           * a sequence (over QQ or CyclotomicField) indexed by range(N) or ZZ/NZZ
           * a sequence (as above) indexed by a finite AbelianGroup
           * a sequnce (as above) indexed by a complete set of representatives of
             the conjugacy classes of a finite permutation group
           * a sequnce (as above) indexed by a complete set of representatives of
             the conjugacy classes of a finite matrix group
*  dct, dst  (for discrete Fourier/Cosine/Sine transform)
*  convolution (in convolution and convolution_periodic)
*  fft, ifft - (fast fourier transforms) wrapping GSL's gsl_fft_complex_forward, gsl_fft_complex_inverse,
         using William Stein's FastFourierTransform class
*  dwt, idwt - (fast wavelet transforms) wrapping GSL's gsl_dwt_forward, gsl_dwt_backward
         using Joshua Kantor's WaveletTransform class. Allows for waveltes of
         type: "haar", "daubechies", "daubechies_centered",
          "haar_centered", "bspline", "bspline_centered".

TODO:
 - "filtered" DFTs.

AUTHOR: David Joyner (10-2006)

"""

##########################################################################
#  Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL):
#
#                  http://www.gnu.org/licenses/
##########################################################################

from sage.rings.number_field.number_field import CyclotomicField
from sage.plot.plot import PolygonFactory
from sage.groups.abelian_gps.dual_abelian_group import DualAbelianGroup
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.matrix_gps.matrix_group import MatrixGroup
from sage.rings.integer_ring import IntegerRing
ZZ = IntegerRing()

class Collection:
    def __init__(self, list, indexset):
        r"""
        \code{indexset} must be a SAGE object with an _iter_ method
        containing the same number of elements as self, which is a
        list of elements taken from a field.

        EXAMPLES:
            sage: s = Collection(A,I)
            sage: s.dict()

            {0: 1/10,
             1: 1/10,
             2: 1/10,
             3: 1/10,
             4: 1/10,
             5: 1/10,
             6: 1/10,
             7: 1/10,
             8: 1/10,
             9: 1/10}
            sage: s.list()
            [1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10]
            sage: s.index_set()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: s.base_ring()
            Rational Field
        """
        self._base_ring = list[0].parent()
        self._index_set = indexset
        self._list = list
        dict = {}
        for i in range(len(indexset)):
            dict[indexset[i]] = list[i]
        self._dict = dict

    def dict(self):
        return self._dict

    def list(self):
        return self._list

    def base_ring(self):
        """
        This just returns the common parent R of the N list
        elements. In some applications (say, when computing the
        discrete Fourier transform, dft), it is more accurate to think
        of the base_ring as the group ring QQ(zeta_N)[R].
        """
        return self._base_ring

    def index_set(self):
        return self._index_set

    def __repr__(self):
        return "Collection of elements "+str(self.list())+"\n indexed by "+str(self.index_set())

    def __str__(self):
        """
        Implements print method.

        EXAMPLES:
            sage: A = [ZZ(i) for i in range(3)]
            sage: I = range(3)
            sage: s = Collection(A,I)
            sage: s
            Collection of elements [0, 1, 2]
             indexed by [0, 1, 2]
            sage: print s
            Collection of elements [0, 1, 2] indexed by [0, 1, 2]
            sage: I = GF(3)
            sage: A = [i^2 for i in I]
            sage: s = Collection(A,I)
            sage: s
            Collection of elements [0, 1, 1]
             indexed by Finite Field of size 3
        """
        return "Collection of elements "+str(self.list())+" indexed by "+str(self.index_set())

    def plot_histogram(self):
        """
        Plots the histogram plot of the sequence, which is assumed to be real
        or from a finite field, with a real indexing set I coercible into RR.

        EXAMPLES:
            sage: J = GF(3)
            sage: A = [ZZ(i^2)+1 for i in J]
            sage: s = Collection(A,J)
            sage: P = s.plot_histogram()

        Now type show(P) to view this in a browser.
        """
        F = self.base_ring()   ## elements must be coercible into RR
        I = self.index_set()
        N = len(I)
        S = self.list()
        P = [polygon([[RR(I[i]),0],[RR(I[i]),RR(S[i])],[RR(I[i+1]),RR(S[i])],[RR(I[i+1]),0],[RR(I[i]),0]], rgbcolor=(1,0,0)) for i in range(N-1)]
        return sum(P)

    def plot(self):
        """
        Plots the points the sequence, which is assumed to be real
        or from a finite field, with a real indexing set I = range(len(self)).

        EXAMPLES:
            sage: I = GF(3).list()
            sage: A = [ZZ(i^2)+1 for i in I]
            sage: s = Collection(A,I)
            sage: P = s.plot()

        Now type show(P) to view this in a browser.
        """
        F = self.base_ring()   ## elements must be coercible into RR
        I = self.index_set()
        N = len(I)
        S = self.list()
        P = line([[RR(I[i]),RR(S[i])] for i in range(N-1)])
        return P

    def dft(self, chi = lambda x: x):
        """
        Implements a discrete Fourier transform "over QQ" using exact
        N-th roots of unity.

        EXAMPLES:
            sage: J = range(5)
            sage: A = [ZZ(1) for i in J]
            sage: s = Collection(A,J)
            sage: s.dft()
            Collection of elements [5, 0, 0, 0, 0]
              indexed by [0, 1, 2, 3, 4]
            sage: G = SymmetricGroup(3)
            sage: J = G.conjugacy_classes_representatives()
            sage: J
            [(), (1,2), (1,2,3)]
            sage: s = Collection([1,2,3],J)
            sage: s.dft()
            Collection of elements [8, 2, 2]
             indexed by [(), (1,2), (1,2,3)]
            sage: F = AbelianGroup(2,[2,3],names='ab')
            sage: s = Collection([1,2,3,4,5,6],F)
            sage: s.dft()
            Collection of elements [21.000000000000000,-3.0000000000000027 - 1.7320508075688741*I,-2.9999999999999964 + 1.7320508075688821*I,-9.0000000000000000 + 0.0000000000000020751629580000000*I,0.0000000000000013322676295501878 - 0.0000000000000026645352591003757*I,-0.00000000000000088817841970012523 - 0.0000000000000026645352591003757*I]
             indexed by Multiplicative Abelian Group isomorphic to C2 x C3
            sage: p = 7; J = range(p); A = [kronecker_symbol(j,p) for j in J]
            sage: s = Collection(A,J); s
            Collection of elements [0, 1, 1, -1, 1, -1, -1]
             indexed by [0, 1, 2, 3, 4, 5, 6]
            sage: Fs = s.dft().list()
            sage: s.list()
            [0, 1, 1, -1, 1, -1, -1]

        The DFT of the values of the quadratic residue symbol is itself, up to
        a constant factor.

        TODO: Read the parent of the elements of S; if QQ or CC leave as
        is; if AbelianGroup, use abelian_group_dual; if some other
        implemented Group (permutation, matrix), call .characters()
        and test if the index list is the set of conjugacy classes.
        """
        J = self.index_set()   ## index set of length N
        N = len(J)
        S = self.list()
        F = self.base_ring()   ## elements must be coercible into QQ(zeta_N)
        if J[0] in ZZ and F.is_field() and F.base_ring()==QQ:
            ## assumes J is range(N)
            zeta = CyclotomicField(N).gen()
            FT = [sum([S[i]*chi(zeta**(i*j)) for i in J]) for j in J]
        elif (J[0].parent()).is_abelian() and F == ZZ or (F.is_field() and F.base_ring()==QQ):
            ## assumes J is AbelianGroup(...)
            G = (J[0].parent()).dual_group()
            FT = [sum([S[i]*chi(J[i]) for i in range(N)]) for chi in G]
        elif (J[0].parent()).is_finite() and F == ZZ or (F.is_field() and F.base_ring()==QQ):
            ## assumes J is the list of conj class represetatives of a
            ## PermuationGroup(...) or Matrixgroup(...)
            chi = (J[0].parent()).character_table()
            FT = [sum([S[i]*chi[i][j] for i in range(N)]) for j in range(N)]
        else:
            raise ValueError,"list elements must be in QQ(zeta_"+str(N)+")"
        return Collection(FT,J)

    def idft(self):
        """
        Implements a discrete Fourier transform

        EXAMPLES:
            sage: J = range(5)
            sage: A = [ZZ(1) for i in J]
            sage: s = Collection(A,J)
            sage: fs = s.dft()
            sage: s == fs.idft()
            1
        """
        F = self.base_ring()   ## elements must be coercible into QQ(zeta_N)
        J = self.index_set()   ## must be = range(N)
        N = len(J)
        S = self.list()
        zeta = CyclotomicField(N).gen()
        iFT = [sum([S[i]*zeta^(-i*j) for i in J]) for j in J]
        return Collection(iFT,J)*(1/N)

    def dct(self):
        """
        Implements a discrete Cosine transform

        EXAMPLES:
            sage: J = range(5)
            sage: A = [exp(-2*pi*i*I/5) for i in J]
            sage: s = Collection(A,J)
            sage: s.dct()
            [2.5000000000000004 - 0.00000000000000011102230246251565*I,
             2.5000000000000004 - 0.00000000000000011102230246251565*I,
             2.5000000000000004 - 0.00000000000000011102230246251565*I,
             2.5000000000000004 - 0.00000000000000011102230246251565*I,
             2.5000000000000004 - 0.00000000000000011102230246251565*I]
        """
        F = self.base_ring()   ## elements must be coercible into RR
        J = self.index_set()   ## must be = range(N)
        N = len(J)
        S = self.list()
        FT = [sum([S[i]*cos(2*pi*i/N) for i in J]) for j in J]
        return Collection(FT,J)

    def dst(self):
        """
        Implements a discrete Sine transform

        EXAMPLES:
            sage: J = range(5)
            sage: A = [exp(-2*pi*i*I/5) for i in J]
            sage: s = Collection(A,J)
            sage: s.dst()
            [0.00000000000000016653345369377348 - 2.5000000000000000*I,
             0.00000000000000016653345369377348 - 2.5000000000000000*I,
             0.00000000000000016653345369377348 - 2.5000000000000000*I,
             0.00000000000000016653345369377348 - 2.5000000000000000*I,
             0.00000000000000016653345369377348 - 2.5000000000000000*I]
        """
        F = self.base_ring()   ## elements must be coercible into RR
        J = self.index_set()   ## must be = range(N)
        N = len(J)
        S = self.list()
        FT = [sum([S[i]*sin(2*pi*i/N) for i in J]) for j in J]
        return Collection(FT,J)

    def convolution(self, other):
        """
        Convolves two sequences of the same length (automatically expands
        the shortest one by extending it by 0 if they have different lengths).
        If {a_n} and {b_n} are sequences of length N (n=0,1,...,N-1), extended
        by zero for all n in ZZ, then the convolution is

                 c_j = \sum_{i=0}^{N-1} a_ib_{j-i}.

        INPUT:
            self, other   --  a collection of elements of a ring with
                              index set a finite abelian group (under +)

        OUTPUT:
            self*other -- the Dirichlet convolution

        EXAMPLES:
            sage: I = range(5)
            sage: A = [ZZ(1) for i in I]
            sage: B = [ZZ(1) for i in I]
            sage: s = Collection(A,I)
            sage: t = Collection(B,I)
            sage: s.convolution(t)
            [1, 2, 3, 4, 5, 4, 3, 2, 1]

        AUTHOR: David Joyner (9-2006)
        """
        S = self.list()
        T = other.list()
        I0 = self.index_set()
        J0 = other.index_set()
        F = self.base_ring()
        E = other.base_ring()
        if F!=E:
            raise TypeError,"Collections must have same base ring"
        if I0!=J0:
            raise TypeError,"Collections must have same index set"
        M = len(S)
        N = len(T)
        if M<N:                    ## first, extend by 0 if necessary
            a = [S[i] for i in range(M)]+[F(0) for i in range(2*N)]
            b = T+[E(0) for i in range(2*M)]
        if M>N:                  ## python trick - a[-j] is really j from the *right*
            b = [T[i] for i in range(N)]+[E(0) for i in range(2*M)]
            a = S+[F(0) for i in range(2*M)]
        if M==N:                ## so need only extend by 0 to the *right*
            a = S+[F(0) for i in range(2*M)]
            b = T+[E(0) for i in range(2*M)]
        N = max(M,N)
        c = [sum([a[i]*b[j-i] for i in range(N)]) for j in range(2*N-1)]
        #print [[b[j-i] for i in range(N)] for j in range(N)]
        return c

    def convolution_periodic(self, other):
        """
        Convolves two collections indexed by a range(...) of the same length (automatically
        expands the shortest one by extending it by 0 if they have different lengths).
        If {a_n} and {b_n} are sequences of length N (n=0,1,...,N-1), extended
        periodically for all n in ZZ, then the convolution is

                 c_j = \sum_{i=0}^{N-1} a_ib_{j-i}.

        INPUT:
            self, other   --  a sequence of elements of CC, RR or GF(q)

        OUTPUT:
            self*other -- the Dirichlet convolution

        EXAMPLES:
            sage: I = range(5)
            sage: A = [ZZ(1) for i in I]
            sage: B = [ZZ(1) for i in I]
            sage: s = Collection(A,I)
            sage: t = Collection(B,I)
            sage: s.convolution_periodic(t)
            [5, 5, 5, 5, 5]

        AUTHOR: David Joyner (9-2006)
        """
        S = self.list()
        T = other.list()
        I = self.index_set()
        J = other.index_set()
        F = self.base_ring()
        E = other.base_ring()
        if F!=E:
           raise TypeError,"Collections must have same parent"
        if I!=J:
            raise TypeError,"Collections must have same index set"
        M = len(S)
        N = len(T)
        if M<N:                    ## first, extend by 0 if necessary
            a = [S[i] for i in range(M)]+[F(0) for i in range(N-M)]
            b = other
        if M>N:
            b = [T[i] for i in range(N)]+[E(0) for i in range(M-N)]
            a = self
        if M==N:
            a = S
            b = T
        N = max(M,N)
        c = [sum([a[i]*b[(j-i)%N] for i in range(N)]) for j in range(2*N-1)]
        return c

    def __rmul__(self,other):
        """
        Implements scalar multiplication (on the right).

        EXAMPLES:
            sage: J = range(5)
            sage: A = [ZZ(1) for i in J]
            sage: s = Collection(A,J)
            sage: t = s*(1/3); t; t.base_ring()
            Collection of elements [1/3, 1/3, 1/3, 1/3, 1/3]
             indexed by [0, 1, 2, 3, 4]
             Rational Field

        """
        S = self.list()
        J = self.index_set()
        F = self.base_ring()
        #if not(other in F):
        #    raise TypeError,"The base rings must be consistent"
        S1 = [S[i]*other for i in J]
        return Collection(S1,J)

    def __eq__(self,other):
        """
        Implements boolean ==.

        EXAMPLES:
            sage: J = range(5)
            sage: A = [ZZ(1) for i in J]
            sage: s = Collection(A,J)
            sage: t = s*(1/3)
            sage: t*3==s
            1
        """
        S = self.list()
        T = other.list()
        I = self.index_set()
        J = other.index_set()
        F = self.base_ring()
        E = other.base_ring()
        if I!=J:
            return 0
        for i in I:
            if abs(S[i]-T[i])> 10^(-8):  ## tests if they differ as reals
                return 0
        #if F!=E:               ## omitted this test since it
        #    return 0           ## doesn't take into account coercions
        return 1

    def fft(self):
        """
        Wraps the gsl FastFourierTransform.forward in fft.pyx (written
        by William Stein). If the length is a power of 2 then this automatically
        uses the radix2 method. If the number of sample
        points in the input is a power of 2 then the wrapper for the
        GSL function gsl_fft_complex_radix2_forward is automatically called.
        Otherwise, gsl_fft_complex_forward is used.

        EXAMPLES:
            sage: J = range(5)
            sage: A = [RR(1) for i in J]
            sage: s = Collection(A,J)
            sage: t = s.fft(); t
            Collection of elements [5.0000000000000000, 0, 0, 0, 0]
             indexed by [0, 1, 2, 3, 4]

        """
        F = self.base_ring()   ## elements must be coercible into RR
        J = self.index_set()   ## must be = range(N)
        N = len(J)
        S = self.list()
        a = FastFourierTransform(N)
        for i in range(N):
            a[i] = S[i]
        a.forward_transform()
        return Collection([a[j][0]+I*a[j][1] for j in J],J)

    def ifft(self):
        """
        Implements the gsl FastFourierTransform.inverse in fft.pyx.
        If the number of sample points in the input is a power of 2
        then the wrapper for the GSL function
        gsl_fft_complex_radix2_inverse is automatically called.
        Otherwise, gsl_fft_complex_inverse is used.

        EXAMPLES:
            sage: J = range(5)
            sage: A = [RR(1) for i in J]
            sage: s = Collection(A,J)
            sage: t = s.fft(); t
            Collection of elements [5.0000000000000000, 0, 0, 0, 0]
             indexed by [0, 1, 2, 3, 4]
            sage: t.ifft()
            [(1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0)]
            sage: t.ifft() == s
            1

        """
        F = self.base_ring()   ## elements must be coercible into RR
        J = self.index_set()   ## must be = range(N)
        N = len(J)
        S = self.list()
        a = FastFourierTransform(N)
        for i in range(N):
            a[i] = S[i]
        a.inverse_transform()
        return Collection([a[j][0]+I*a[j][1] for j in J],J)

    def dwt(self,other="haar",wavelet_k=2):
        """
        Wraps the gsl WaveletTransform.forward in dwt.pyx (written
        by Johua Kantor). Assumes the length of the sample is a power of 2.
        Uses the GSL function gsl_wavelet_transform_forward.

        other -- the wavelet_type:   the name of the type of wavelet,
                                     valid choices are:
                                     'daubechies','daubechies_centered',
                                     'haar' (default),'haar_centered',
                                     'bspline', and 'bspline_centered'.

        wavelet_k -- For daubechies wavelets, wavelet_k specifies a
                     daubechie wavelet with k/2 vanishing moments.
                     k = 4,6,...,20 for k even are the only ones implemented.
                     For Haar wavelets, wavelet_k must be 2.
                     For bspline wavelets,
                     wavelet_k = 103,105,202,204,206,208,301,305, 307,309
                     will give biorthogonal B-spline wavelets of order (i,j) where
                     wavelet_k=100*i+j.

        The wavelet transform uses J=log_2(n) levels.

        EXAMPLES:
            sage: J = range(7)
            sage: A = [RR(1) for i in J]
            sage: s = Collection(A,J)
            sage: t = s.dwt(); t
            Collection of elements [2.8284271247500001,-0.000000000000000091940344226800001,0.000000000000000085326710974600004, 0.000000000000000085326710974600004, 0.00000000000000000, 0.00000000000000000, 0.00000000000000000, 0.00000000000000000]
            indexed by [0, 1, 2, 3, 4, 5, 6, 7]

        """
        F = self.base_ring()   ## elements must be coercible into RR
        J = self.index_set()   ## must be = range(N)
        N = len(J)             ## must be 1 minus a power of 2
        S = self.list()
        if other=="haar" or other=="haar_centered":
            if k in [2]:
                a = WaveletTransform(N,other,wavelet_k)
            else:
                raise ValueError,"wavelet_k must be = 2"
        if other=="debauchies" or other=="debauchies_centered":
            if k in [4,6,8,10,12,14,16,18,20]:
                a = WaveletTransform(N,other,wavelet_k)
            else:
                raise ValueError,"wavelet_k must be in {4,6,8,10,12,14,16,18,20}"
        if other=="bspline" or other=="bspline_centered":
            if k in [103,105,202,204,206,208,301,305,307,309]:
                a = WaveletTransform(N,other,103)
            else:
                raise ValueError,"wavelet_k must be in {103,105,202,204,206,208,301,305,307,309}"
        for i in range(N):
            a[i] = S[i]
        a.forward_transform()
        return Collection([RR(a[j]) for j in J],J)

    def idwt(self,other="haar",wavelet_k=2):
        """
        Implements the gsl WaveletTransform.backward in dwt.pyx.
        other must be an element of
         {"haar", "daubechies", "daubechies_centered",
          "haar_centered", "bspline", "bspline_centered"}.
        Assumes the length of the sample is a power of 2. Uses the
        GSL function gsl_wavelet_transform_backward.
        See the above docstring for dwt for further details.

        EXAMPLES:
            sage: J = range(8)
            sage: A = [RR(1) for i in J]
            sage: s = Collection(A,J)
            sage: t = s.dwt(); t
            Collection of elements [2.8284271247500001,
             -0.000000000000000091940344226800001,
             0.000000000000000085326710974600004,
             0.000000000000000085326710974600004,
             0.00000000000000000, 0.00000000000000000,
             0.00000000000000000, 0.00000000000000000]
             indexed by [0, 1, 2, 3, 4, 5, 6, 7]
            sage: t.idwt()
            Collection of elements [1.0000000000000000, 1.0000000000000000,
             1.0000000000000000, 1.0000000000000000, 1.0000000000000000,
             1.0000000000000000, 1.0000000000000000, 1.0000000000000000]
             indexed by [0, 1, 2, 3, 4, 5, 6, 7]
            sage: t.idwt() == s
             1
            sage: J = range(16)
            sage: A = [RR(1) for i in J]
            sage: s = Collection(A,J)
            sage: t = s.dwt("bspline"); t
             Collection of elements [4.0000000000000000,
             -0.00000000000000014333152720300001,
             -0.000000000000000091940344226800001,
             -0.000000000000000091940344226800001,
             0.000000000000000085326710974600004,
             0.000000000000000085326710974600004,
             0.000000000000000085326710974600004,
             0.000000000000000085326710974600004, 0.00000000000000000,
             0.00000000000000000, 0.00000000000000000,
             0.00000000000000000, 0.00000000000000000, 0.00000000000000000,
             0.00000000000000000, 0.00000000000000000]
             indexed by [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            sage: t.idwt("bspline") == s
             1

        """
        F = self.base_ring()   ## elements must be coercible into RR
        J = self.index_set()   ## must be = range(N)
        N = len(J)             ## must be 1 minus a power of 2
        S = self.list()
        if other=="haar" or other=="haar_centered":
            if k in [2]:
                a = WaveletTransform(N,other,wavelet_k)
            else:
                raise ValueError,"wavelet_k must be = 2"
        if other=="debauchies" or other=="debauchies_centered":
            if k in [4,6,8,10,12,14,16,18,20]:
                a = WaveletTransform(N,other,wavelet_k)
            else:
                raise ValueError,"wavelet_k must be in {4,6,8,10,12,14,16,18,20}"
        if other=="bspline" or other=="bspline_centered":
            if k in [103,105,202,204,206,208,301,305,307,309]:
                a = WaveletTransform(N,other,103)
            else:
                raise ValueError,"wavelet_k must be in {103,105,202,204,206,208,301,305,307,309}"
        for i in range(N):
            a[i] = S[i]
        a.backward_transform()
        return Collection([RR(a[j]) for j in J],J)

