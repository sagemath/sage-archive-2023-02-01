"""
Counting Primes

AUTHORS:
    * R. Andrew Ohana
    * William Stein

TESTS::

    sage: z = sage.functions.prime_pi.PrimePi()
    sage: loads(dumps(z))
    Function that counts the number of primes up to x
    sage: loads(dumps(z)) == z
    True
"""

#*****************************************************************************
#       Copyright (C) 2009 R. Andrew Ohana <ohanar@u.washington.edu>
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../ext/interrupt.pxi'

from sage.rings.integer import Integer
from sage.rings.fast_arith import prime_range
import sage.plot.all

cdef extern from "math.h":
    double sqrt(double)
    double floor(double)
    float sqrtf(float)
    float floorf(float)

cdef inline long sqrt_longlong(long long N) except -1:
    """
    Return the truncated integer part of the square root of N, i.e., the
    floor of sqrt(N).
    """
    # code from Bill Hart.
    cdef long long m = <long> floor(sqrt(<double> N))
    cdef long long res = m + ((m+1)*(m+1) <= N)  # add 1 if too small.
    res = res - (res*res > N)  # subtract 1 if too big.
    if res*res <= N and (res+1)*(res+1) > N:
        return res
    else:
        raise RuntimeError, "there is a bug in prime_pi.pyx's sqrt_longlong (res=%s, N=%s). Please report!"%(res,N)

cdef inline long sqrt_long(long N) except -1:
    """
    Return the truncated integer part of the square root of N, i.e., the
    floor of sqrt(N).
    """
    if sizeof(long) == 8:
        return sqrt_longlong(N)
    # code from Bill Hart.
    cdef long m = <long> floorf(sqrtf(<float> N))
    cdef long res = m + ((m+1)*(m+1) <= N)  # add 1 if too small.
    res = res - (res*res > N)  # subtract 1 if too big.
    if res*res <= N and (res+1)*(res+1) > N:
        return res
    else:
        raise RuntimeError, "there is a bug in prime_pi.pyx's sqrt_long (res=%s, N=%s). Please report!"%(res,N)


cdef class PrimePi:
    r"""
    Return the number of primes $\leq x$.

    EXAMPLES::

        sage: prime_pi(7)
        4
        sage: prime_pi(100)
        25
        sage: prime_pi(1000)
        168
        sage: prime_pi(100000)
        9592
        sage: prime_pi(0.5)
        0
        sage: prime_pi(-10)
        0
        sage: prime_pi(500509)
        41581

    The prime_pi function allows for use of additional memory::

        sage: prime_pi(500509, 8)
        41581

    The prime_pi function also has a special plotting method, so it plots
    quickly and perfectly as a step function::

        sage: P = plot(prime_pi, 50,100)
    """
    def __repr__(self):
        """
        EXAMPLES::

            sage: prime_pi.__repr__()
            'Function that counts the number of primes up to x'
        """
        return 'Function that counts the number of primes up to x'

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: P = sage.functions.prime_pi.PrimePi()
            sage: P == prime_pi
            True
            sage: P != prime_pi
            False
        """
        if isinstance(other, PrimePi):
            return 0
        return cmp(type(self), type(other))

    cdef long* primes
    cdef long primeCount, maxPrime

    def __call__(self, long long x, long mem_mult = 1):
        r"""
        EXAMPLES::

            sage: prime_pi.__call__(7)
            4
            sage: prime_pi.__call__(100)
            25
            sage: prime_pi.__call__(1000)
            168
            sage: prime_pi.__call__(100000)
            9592
            sage: prime_pi.__call__(0.5)
            0
            sage: prime_pi.__call__(-10)
            0
            sage: prime_pi.__call__(500509)
            41581
            sage: prime_pi.__call__(500509, 8)
            41581

        TESTS:

        Make sure we actually compute correct results::

            sage: for n in (30..40):
            ...       prime_pi(2**n)
            54400028
            105097565
            203280221
            393615806
            762939111
            1480206279
            2874398515
            5586502348
            10866266172
            21151907950
            41203088796

        We know this implementation is broken at least on some 32-bit
        systems for `2^{46}`, so we are capping the maximum allowed value::

            sage: prime_pi(2^40+1)
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of prime_pi() greater 2**40 not implemented

        """
        if mem_mult < 1:
            raise ValueError, "mem_mult must be positive"
        if x < 2:
            return 0
        if x > 1099511627776L:
            raise NotImplementedError, "computation of prime_pi() greater 2**40 not implemented"
        x += x & 1
        # m_max is the current sieving value, for prime counting - this value is sqrt(x)
        cdef long m_max = sqrt_longlong(x) + 1
        cdef long prime
        cdef long i = 0
        # mem_mult is the multiplier that determines how large a list of primes built
        tempList = prime_range(mem_mult*m_max)
        self.primeCount = len(tempList)
        if self.primeCount > 1:
            self.maxPrime = tempList[self.primeCount - 1]
        self.primes = <long*> sage_malloc(self.primeCount*sizeof(long))
        for prime in tempList:
            self.primes[i] = prime
            i += 1
        _sig_on
        s = Integer(self.prime_phi_large(x, m_max))
        _sig_off
        sage_free(self.primes)
        return s

    cdef long long prime_phi_large(self, long long N, long m_max = 0):
        """
        Uses Legendre's Formula for computing pi(x). Both this and the prime_phi_small
        methods are specialized versions of Legendre's Formula that takes advantage of
        its recursive properties. Hans Riesel's Prime "Numbers and Computer Methods for
        Factorization" discusses Legnedre's Formula, however there are many places
        that discuss it as nearly every prime counting algorithm uses it.
        """
        try:
            return self.prime_phi_small(long(N), m_max)
        except OverflowError:
            pass
        # m_max is the current sieving value
        if m_max == 0:
            m_max = sqrt_longlong(N) + 1
        cdef long long sum = ( (N >> 1) - (N-4)/6 - (N-16)/10 + (N-16)/30 - (N-8)/14 \
                               + (N-22)/42 + (N-106)/70 - (N-106)/210 + 2 )
        cdef long prime = 11
        cdef bint preTransition = True
        cdef long long y
        cdef long i
        for i in range(4, self.primeCount):
            if prime >= m_max: break
            y = (N+prime-1)/(2*prime) << 1
            if preTransition:
                if prime*prime <= y:
                    sum -= self.prime_phi_large(y, prime) - i
                else:
                    sum -= self.prime_phi_large(y) - i
                    preTransition = False
            else:
                sum -= self.prime_phi_large(y) - i
            prime = self.primes[i + 1]
        return sum

    cdef long prime_phi_small(self, long N, long m_max = 0):
        if N < 2:
            return 0
        N += N & 1
        if N < 9:
            return N >> 1
        if N < 25:
            return (N >> 1) - (N-4)/6
        if N < 49:
            return (N >> 1) - (N-4)/6 - (N-16)/10 + (N-16)/30
        if N < 121:
            return (N >> 1) - (N-4)/6 - (N-16)/10 + (N-16)/30 - (N-36)/14 + (N-22)/42
        # m_max is the current sieving value
        if m_max == 0:
            if N < self.maxPrime:
                return self.bisect(N) + 1
            m_max = sqrt_long(N) + 1
        cdef long sum = ( (N >> 1) - (N-4)/6 - (N-16)/10 + (N-16)/30 - (N-8)/14 \
                          + (N-22)/42 + (N-106)/70 - (N-106)/210 + 2 )
        cdef long prime = 11
        cdef bint preTransition = True
        cdef long x, i
        for i in range(4, self.primeCount):
            if prime >= m_max: break
            x = N / prime
            if preTransition:
                if prime*prime <= x:
                    sum -= self.prime_phi_small(x, prime) - i
                else:
                    sum -= self.prime_phi_small(x) - i
                    preTransition = False
            else:
                sum -= self.prime_phi_small(x) - i
            prime = self.primes[i + 1]
        return sum

    cdef long bisect(self, long N):
        cdef long size = self.primeCount >> 1
        if N >> 2 < size:
            size = N >> 2
        cdef long position = size
        cdef long prime
        while size > 0:
            prime = self.primes[position]
            size >>= 1
            if prime < N:
                position += size
            elif prime > N:
                position -= size
            else: break
        if position < self.primeCount - 1:
            prime = self.primes[position + 1]
            while prime < N:
                position += 1
                prime = self.primes[position + 1]
        prime = self.primes[position]
        while prime > N:
            position -= 1
            prime = self.primes[position]
        return position

    def plot(self, xmin=0, xmax=100, vertical_lines=True, **kwds):
        """
        Draw a plot of the prime counting function from xmin to xmax.
        All additional arguments are passed on to the line command.

        WARNING: we draw the plot of prime_pi as a stairstep function
        with explicitly drawn vertical lines where the function
        jumps. Technically there should not be any vertical lines, but
        they make the graph look much better, so we include them.
        Use the option ``vertical_lines=False`` to turn these off.

        EXAMPLES::

            sage: plot(prime_pi, 1, 100)
            sage: prime_pi.plot(-2,50,thickness=2, vertical_lines=False)
        """
        y = self(xmin)
        v = [(xmin, y)]
        for p in prime_range(xmin+1, xmax+1):
            y += 1
            v.append((p,y))
        v.append((xmax,y))
        from sage.plot.step import plot_step_function
        return plot_step_function(v, vertical_lines=vertical_lines, **kwds)

#############
prime_pi = PrimePi()

