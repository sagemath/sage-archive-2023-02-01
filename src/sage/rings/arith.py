"""
Miscellaneous arithmetic functions
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import math

import sage.misc.misc as misc
import sage.misc.search
from sage.libs.pari.all import pari
import sage.rings.rational_field
import sage.rings.integer_ring
import sage.rings.integer
import sage.rings.rational
import sage.rings.real_field
import sage.rings.complex_field
import sage.rings.complex_number
import sage.rings.real_mpfr
import sage.structure.factorization as factorization
from sage.structure.element import RingElement, canonical_coercion, bin_op
from sage.interfaces.all import gp

##################################################################
# Elementary Arithmetic
##################################################################

def algdep(z, n):
    """
    Returns a polynomial of degree at most $n$ which is approximately
    satisfied by the number $z$.  Note that the returned polynomial
    need not be irreducible, and indeed usually won't be if $z$ is a good
    approximation to an algebraic number of degree less than $n$.

    ALGORITHM: Uses the PARI C-library algdep command.

    INPUT:
        z -- real, complex, or $p$-adic number
        n -- an integer

    EXAMPLES:
        sage: algdep(1.888888888888888, 1)
        9*x - 17
        sage: algdep(0.12121212121212,1)
        33*x - 4
        sage: algdep(sqrt(2),2)
        x^2 - 2

    This example involves a complex number.
        sage: z = (1/2)*(1 + sqrt(3) *CC.0); z
        0.500000000000000 + 0.866025403784438*I
        sage: p = algdep(z, 6); p
        x^6 + 2*x^3 + 1                      # 32-bit
        x^5 - x^4 + x^3 + x^2 - x + 1        # 64-bit
        sage: p.factor()
        (x + 1)^2 * (x^2 - x + 1)^2          # 32-bit
        (x + 1) * (x^2 - x + 1)^2            # 64-bit
        sage: z^2 - z + 1
        0.000000000000000111022302462515

    This example involves a $p$-adic number.
        sage: K = pAdicField(3)
        sage: a = K(7/19); a
        1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
        sage: algdep(a, 1)
        19*x - 7
    """

    # TODO -- change to use PARI C library???
    import sage.rings.polynomial_ring
    x = sage.rings.polynomial_ring.PolynomialRing(
        sage.rings.integer_ring.IntegerRing(), 'x').gen()

    if isinstance(z, (int, long, sage.rings.integer.Integer)):
        return x - sage.rings.integer.Integer(z)

    n = sage.rings.integer.Integer(n)

    if isinstance(z, (sage.rings.rational.Rational)):
        return z.denominator()*x   -   z.numerator()

    if isinstance(z, float):
        z = sage.rings.real_field.RealField()(z)
    elif isinstance(z, complex):
        z = sage.rings.complex_field.ComplexField()(z)

    if misc.is_64_bit and isinstance(z, (sage.rings.real_mpfr.RealNumber,
                                         sage.rings.complex_number.ComplexNumber)):
        bits = int(float(z.prec()/3))
        if bits == 0:
            bits = 1
        f = pari(z).algdep(n, bits)
    else:
        y = pari(z)
        f = y.algdep(n)

    return x.parent()(list(reversed(eval(str(f.Vec())))))


algebraic_dependency = algdep

def bernoulli(n, algorithm='pari'):
    r"""
    Return the n-th Bernoulli number, as a rational number.

    INPUT:
        n -- an integer
        algorithm:
            'pari' -- (default) use the PARI C library, which is
                      by *far* the fastest.
            'gap'  -- use GAP
            'gp'   -- use PARI/GP interpreter
            'magma' -- use MAGMA (optional)
            'python' -- use pure Python implementation

    EXAMPLES:
        sage: bernoulli(12)
        -691/2730
        sage: bernoulli(50)
        495057205241079648212477525/66

    We use of each of the alternative algorithms:
        sage: bernoulli(12, algorithm='gap')
        -691/2730
        sage: bernoulli(12, algorithm='gp')
        -691/2730
        sage: bernoulli(12, algorithm='magma')           # optional
        -691/2730
        sage: bernoulli(12, algorithm='pari')
        -691/2730
        sage: bernoulli(12, algorithm='python')
        -691/2730

    \note{If $n>50000$ then algorithm = 'gp' is used instead of
    algorithm = 'pari', since the C-library interface to PARI
    is limited in memory for individual operations.}

    AUTHOR: David Joyner and William Stein
    """
    from sage.rings.all import Integer, Rational
    n = Integer(n)
    if n > 50000 and algorithm == 'pari':
        algorithm = 'gp'
    if algorithm == 'pari':
        x = pari(n).bernfrac()         # Use the PARI C library
        return Rational(x)
    elif algorithm == 'gap':
        import sage.interfaces.gap
        x = sage.interfaces.gap.gap('Bernoulli(%s)'%n)
        return Rational(x)
    elif algorithm == 'magma':
        import sage.interfaces.magma
        x = sage.interfaces.magma.magma('Bernoulli(%s)'%n)
        return Rational(x)
    elif algorithm == 'gp':
        import sage.interfaces.gp
        x = sage.interfaces.gp.gp('bernfrac(%s)'%n)
        return Rational(x)
    elif algorithm == 'python':
        import sage.rings.bernoulli
        return sage.rings.bernoulli.bernoulli_python(n)
    else:
        raise ValueError, "invalid choice of algorithm"

def Li(x):
    r"""
    Return value of the function Li(x), which is by definition
    $$
       \int_2^{x} dt / \log(t).
    $$

    The function Li(x) is an approximation for the number
    of primes up to $x$.  In fact, the famous Riemann
    Hypothesis is equivalent to the statement that for
    $x \geq 2.01$ we have
    $$
        |\pi(x) - Li(x)| \leq \sqrt{x} \log(x).
    $$
    For ``small'' $x$, $Li(x)$ is always slightly bigger than
    $\pi(x)$.  However it is a theorem that there are (very large,
    e.g., around $10^{316}$) values of $x$ so that $\pi(x) > Li(x)$.
    See ``A new bound for the smallest x with $\pi(x) > li(x)$'',
    Bays and Hudson, Mathematics of Computation, 69 (2000) 1285--1296.

    ALGORITHM: Computed numerically using PARI.

    INPUT:
        x -- a real number >= 2.

    OUTPUT:
        x -- a real double

    EXAMPLES:
        sage: pari.init_primes(10^6)   # needed to compute prime_pi(10^6)
        sage: for n in range(1,7):
        ...    print '%-10s%-10s%-20s'%(10^n, prime_pi(10^n), Li(10^n))
        10        4         5.12043572467
        100       25        29.080977804
        1000      168       176.56449421
        10000     1229      1245.09205212
        100000    9592      9628.76383727
        1000000   78498     78626.5039957
    """
    from real_double import RDF
    x = RDF(x)
    return RDF(gp('intnum(t=2,%s,1/log(t))'%x))

def prime_pi(x):
    """
    Return the number of primes $\leq x$.

    EXAMPLES:
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
    """
    if x < 2:
        return sage.rings.integer.Integer(0)
    return sage.rings.integer.Integer(pari(x).primepi())

number_of_primes = prime_pi


def factorial(n, algorithm='gmp'):
    r"""
    Compute the factorial of $n$, which is the product
    $1\cdot 2\cdot 3 \cdots (n-1) n$.

    INPUT:
        n -- an integer
        algorithm -- string (default: 'gmp')
             'gmp' -- use the GMP C-library factorial function
             'pari' -- use PARI's factorial function

    OUTPUT:
        an integer

    EXAMPLES:
        sage: factorial(0)
        1
        sage: factorial(4)
        24
        sage: factorial(10)
        3628800
        sage: factorial(1) == factorial(0)
        True
        sage: factorial(6) == 6*5*4*3*2
        True
        sage: factorial(1) == factorial(0)
        True
        sage: factorial(71) == 71* factorial(70)
        True
        sage: factorial(-32)
        Traceback (most recent call last):
        ...
        ValueError: factorial -- must be nonnegative

    PERFORMANCE:
    This discussion is valid as of April 2006.  All timings
    below are on a Pentium Core Duo 2Ghz MacBook Pro running Linux
    with a 2.6.16.1 kernel.

    \begin{itemize}
       \item It takes less than a minute to compute the factorial of
          $10^7$ using the GMP algorithm, and the factorial of $10^6$
          takes less than 4 seconds.

       \item The GMP algorithm is faster and more memory efficient
          than the PARI algorithm.  E.g., PARI computes $10^7$
          factorial in 100 seconds on the core duo 2Ghz.

       \item For comparison, computation in Magma $\leq$ 2.12-10 of
             $n!$ is best done using \code{&*[1..n]}.
             It takes 113 seconds to compute the factorial of $10^7$
             and 6 seconds to compute the factorial of $10^6$.
             Mathematica V5.2 compute the factorial of $10^7$ in
             136 seconds and the factorial of $10^6$ in 7 seconds.
             (Mathematica is notably very efficient at memory usage
             when doing factorial calculations.)
    \end{itemize}

    """
    if n < 0:
        raise ValueError, "factorial -- must be nonnegative"
    Z = sage.rings.integer.Integer
    if algorithm == 'gmp':
        return Z(n).factorial()
    elif algorithm == 'pari':
        return Z(pari('%s!'%Z(n)))
    else:
        raise ValueError, 'unknown algorithm'

def is_prime(n, flag=0):
    r"""
    Returns True if $x$ is prime, and False otherwise.  The result
    is proven correct -- {\em this is NOT a pseudo-primality test!}.

    INPUT:
        flag -- int
                0 (default): use a combination of algorithms.
                1: certify primality using the Pocklington-Lehmer Test.
                2: certify primality using the APRCL test.
    OUTPUT:
        bool -- True or False

    \note{We do not consider negatives of prime numbers as prime.}

    EXAMPLES::
        sage: is_prime(389)
        True
        sage: is_prime(2000)
        False
        sage: is_prime(2)
        True
        sage: is_prime(-1)
        False
        sage: factor(-6)
        -1 * 2 * 3
        sage: is_prime(1)
        False
        sage: is_prime(-2)
        False

    IMPLEMENTATION: Calls the PARI isprime function.
    """
    n = sage.rings.integer.Integer(n)
    return pari(n).isprime()

def is_pseudoprime(n, flag=0):
    r"""
    Returns True if $x$ is a pseudo-prime, and False otherwise.  The result
    is \em{NOT} proven correct -- {\em this is a pseudo-primality test!}.

    INPUT:
        flag -- int
                0 (default): checks whether x is a Baillie-Pomerance-Selfridge-Wagstaff pseudo prime (strong Rabin-Miller pseudo prime for base 2, followed by strong Lucas test for the sequence (P,-1), P smallest positive integer such that P^2 - 4 is not a square mod x).
                > 0: checks whether x is a strong Miller-Rabin pseudo prime for flag randomly chosen bases (with end-matching to catch square roots of -1).

    OUTPUT:
        bool -- True or False

    \note{We do not consider negatives of prime numbers as prime.}

    EXAMPLES::
        sage: is_pseudoprime(389)
        True
        sage: is_pseudoprime(2000)
        False
        sage: is_pseudoprime(2)
        True
        sage: is_pseudoprime(-1)
        False
        sage: factor(-6)
        -1 * 2 * 3
        sage: is_pseudoprime(1)
        False
        sage: is_pseudoprime(-2)
        False

    IMPLEMENTATION: Calls the PARI ispseudoprime function.
    """
    n = sage.rings.integer.Integer(n)
    return pari(n).ispseudoprime()

def is_prime_power(n, flag=0):
    r"""
    Returns True if $x$ is a prime power, and False otherwise.  The result
    is proven correct -- {\em this is NOT a pseudo-primality test!}.

    INPUT:
        n -- an integer
        flag (for primality testing) -- int
                0 (default): use a combination of algorithms.
                1: certify primality using the Pocklington-Lehmer Test.
                2: certify primality using the APRCL test.

    OUTPUT:
        bool -- True or False

    IMPLEMENTATION: Calls the PARI isprime and ispower functions.

    EXAMPLES::
        sage: is_prime_power(389)
        True
        sage: is_prime_power(2000)
        False
        sage: is_prime_power(2)
        True
        sage: is_prime_power(1024)
        True
        sage: is_prime_power(-1)
        False
        sage: is_prime_power(1)
        True
        sage: is_prime_power(997^100)
        True
    """
    Z = sage.rings.integer.Integer
    n = Z(n)
    if n == 1:
        return True
    if n < 1:
        return False
    if is_prime(n, flag):
        return True
    k, g = pari(n).ispower()
    if not k:
        return False
    return g.isprime(flag)

def valuation(m, p):
    """
    The exact power of p>0 that divides the integer m.
    We do not require that p be prime, and if m is 0,
    then this function returns rings.infinity.

    EXAMPLES::

        sage: valuation(512,2)
        9
        sage: valuation(1,2)
        0

    Valuation of 0 is defined, but valuation with respect to 0 is not::

        sage: valuation(0,7)
        Infinity
        sage: valuation(3,0)
        Traceback (most recent call last):
        ...
        ValueError: valuation at 0 not defined

    Here are some other example::

        sage: valuation(100,10)
        2
        sage: valuation(200,10)
        2
        sage: valuation(243,3)
        5
        sage: valuation(243*10007,3)
        5
        sage: valuation(243*10007,10007)
        1
    """
    if p <= 0:
        raise ValueError, "valuation at 0 not defined"
    if m == 0:
        import sage.rings.all
        return sage.rings.all.infinity
    r=0
    power=p
    while m%power==0:
        r += 1
        power *= p
    return r




def prime_range(start, stop=None, leave_pari=False):
    r"""
    List of all primes between start and stop-1, inclusive.  If the
    second argument is omitted, returns the primes up to the first
    argument.

    Use this function when both start and stop are not too large,
    since in all cases this function makes a table of primes up to
    stop.  If both are large, use the primes iterator function
    instead.

    INPUT:
        start -- lower bound
        stop -- upper bound
        leave_pari -- bool (default: False) if True the returned list
                    is a PARI list; this is *vastly* faster since the
                    time of prime_range is dominated by conversion
                    from PARI to SAGE integers.   However, PARI integers
                    are much different than SAGE integers.
                    If you use this option the lower bound must be 2.

    You can also call this function with \code{prime_range(bound)} to
    get all primes up to bound.

    EXAMPLES:
        sage: prime_range(10)
        [2, 3, 5, 7]
        sage: prime_range(7)
        [2, 3, 5]
        sage: prime_range(2000,2020)
        [2003, 2011, 2017]
        sage: prime_range(2,2)
        []
        sage: prime_range(2,3)
        [2]
        sage: prime_range(10)
        [2, 3, 5, 7]
    """
    if stop is None:
        start, stop = 2, start
    v = pari.primes_up_to_n(stop-1)
    Z = sage.rings.integer.Integer
    if leave_pari:
        if start != 2:
            raise ValueError, "lower bound must be 2 if leave_pari is True"
        return v
    if start <= 2:
        return [Z(p) for p in v]     # this dominates runtime!
    start = pari(start)
    return [Z(p) for p in v if p >= start]     # this dominates runtime!

def primes_first_n(n, leave_pari=False):
    r"""
    Return the first $n$ primes.

    INPUT:
        leave_pari -- bool (default: False) if True the returned list
                    is a PARI list; this is *vastly* (10 times!)
                    faster since the time of prime_range is dominated
                    by conversion from PARI to SAGE integers.
                    However, PARI integers are much different than
                    SAGE integers.  If you use this option the lower
                    bound must be 2.
    OUTPUT:
        a list of the first $n$ prime numbers.

    EXAMPLES:
        sage: primes_first_n(10)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        sage: len(primes_first_n(1000))
        1000

    This is very fast, because we leave the output as a PARI object:
        sage: v = primes_first_n(10^6, leave_pari=True)
        sage: len(v)
        1000000
    """
    v = pari.prime_list(n)
    Z = sage.rings.integer.Integer
    if leave_pari:
        return v
    return [Z(p) for p in v]     # this dominates runtime!


#
# This is from
#    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/366178
# It's impressively fast given that it's in Pure Python.
#
def eratosthenes(n):
    r"""
    Return a list of the primes $\leq n$.

    This is extremely slow and is for educational purposes only.
    """
    n = int(n)
    if n == 2:
        return [2]
    elif n<2:
        return []
    s = range(3,n+3,2)
    mroot = n ** 0.5
    half = (n+1)/2
    i = 0
    m = 3
    while m <= mroot:
        if s[i]:
            j = (m*m-3)/2
            s[j] = 0
            while j < half:
                s[j] = 0
                j += m
        i = i+1
        m = 2*i+3
    return [sage.rings.integer.Integer(2)] + [sage.rings.integer.Integer(x) for x in s if x and x <= n]

# My old versions; not as fast as the above.
## def eratosthenes(n):
##     """
##     Returns a list of the primes up to n, computed
##     using the Sieve of Eratosthenes.
##     Input:
##         n -- a positive integer
##     Output:
##         list -- a list of the primes up to n
##     Examples:
##     sage: eratosthenes(7)
##     [2, 3, 5, 7]
##     sage: eratosthenes(45)
##     [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
##     """
##     if n <= 1: return []
##     X = [i for i in range(3,n+1) if i%2 != 0]
##     P = [2]
##     sqrt_n = sqrt(n)
##     while len(X) > 0 and X[0] <= sqrt_n:
##         p = X[0]
##         P.append(p)
##         X = [a for a in X if a%p != 0]
##     return P + X

prange = prime_range

def primes(start, stop=None):
    r"""
    Returns an iterator over all primes between start and stop-1,
    inclusive.  This is much slower than \code{prime_range}, but
    potentially uses less memory.

    This command is like the xrange command, except it only iterates
    over primes.  In some cases it is better to use primes than
    prime_range, because primes does not build a list of all primes in
    the range in memory all at once.  However it is potentially much
    slower since it simply calls the \code{next_prime} function
    repeatedly, and \code{next_prime} is slow, partly because it
    proves correctness.

    EXAMPLES:
        sage.: for p in primes(5,10):
        ...     print p
        ...
        5
        7
        sage: list(primes(10000000000, 10000000100))
        [10000000019, 10000000033, 10000000061, 10000000069, 10000000097]
    """

    start = sage.rings.integer.Integer(start)
    if stop == None:
        stop = start
        start = sage.rings.integer.Integer(2)
    else:
        stop = sage.rings.integer.Integer(stop)
    n = start - 1
    while True:
        n = next_prime(n)
        if n <= stop:
            yield n
        else:
            return

def next_prime_power(n):
    """
    The next prime power greater than the integer n.  If n is a prime
    power, then this function does not return n, but the next prime
    power after n.

    EXAMPLES:
        sage: next_prime_power(-10)
        1
        sage: is_prime_power(1)
        True
        sage: next_prime_power(0)
        1
        sage: next_prime_power(1)
        2
        sage: next_prime_power(2)
        3
        sage: next_prime_power(10)
        11
        sage: next_prime_power(7)
        8
        sage: next_prime_power(99)
        101
    """
    if n < 0:   # negatives are not prime.
        return sage.rings.integer.Integer(1)
    if n == 2:
        return sage.rings.integer.Integer(3)
    n = sage.rings.integer.Integer(n) + 1
    while not is_prime_power(n):  # pari isprime is provably correct
        n += 1
    return n

def next_prime(n, proof=True):
    """
    The next prime greater than the integer n.  If n is prime, then
    this function does not return n, but the next prime after n.  If
    the optional argument proof is False (the default), this function
    only returns a pseudo-prime, as defined by the PARI nextprime
    function.

    INPUT:
        n -- integer
        proof -- bool (default: True)

    EXAMPLES:
        sage: next_prime(-100)
        2

    Notice that the next_prime(5) is not 5 but 7.
        sage: next_prime(5)
        7
        sage: next_prime(2004)
        2011
    """
    if n < 2:   # negatives are not prime.
        return sage.rings.integer.Integer(2)
    if n == 2:
        return sage.rings.integer.Integer(3)
    if not proof:  # pari nextprime is probabilistic (according to their docs)
        return sage.rings.integer.Integer((eval(str(pari(n+1).nextprime()))))

    if n % 2 == 0:
        n += 1
    else:
        n += 2
    while not is_prime(n):  # pari isprime is provably correct
        n += 2
    return sage.rings.integer.Integer(n)

def previous_prime(n):
    """
    The largest prime < n.  The result is provably
    correct.   If n <= 2, this function raises a ValueError.

    EXAMPLES:
        sage: previous_prime(10)
        7
        sage: previous_prime(7)
        5
        sage: previous_prime(8)
        7
        sage: previous_prime(7)
        5
        sage: previous_prime(5)
        3
        sage: previous_prime(3)
        2
        sage: previous_prime(2)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime
        sage: previous_prime(1)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime
        sage: previous_prime(-20)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime
    """
    n = sage.rings.integer.Integer(n)-1
    if n <= 1:
        raise ValueError, "no previous prime"
    if n <= 3:
        return sage.rings.integer.Integer(n)
    if n%2 == 0:
        n -= 1
    while not is_prime(n):
        n -= 2
    return sage.rings.integer.Integer(n)

def previous_prime_power(n):
    r"""
    The largest prime power $< n$.  The result is provably
    correct. If $n \leq 2$, this function returns $-x$,
    where $x$ is prime power and $-x < n$ and no larger negative
    of a prime power has this property.

    EXAMPLES:
        sage: previous_prime_power(2)
        1
        sage: previous_prime_power(10)
        9
        sage: previous_prime_power(7)
        5
        sage: previous_prime_power(127)
        125

        sage: previous_prime_power(0)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime power
        sage: previous_prime_power(1)
        Traceback (most recent call last):
        ...
        ValueError: no previous prime power

        sage: n = previous_prime_power(2^16 - 1)
        sage: while is_prime(n):
        ...    n = previous_prime_power(n)
        sage: factor(n)
        251^2
    """
    n = sage.rings.integer.Integer(n)-1
    if n <= 0:
        raise ValueError, "no previous prime power"
    while not is_prime_power(n):
        n -= 1
    return n

def random_prime(n):
    """
    Returns a random prime p satisfying between 2 and p (i.e. 2 <= p <= n).

    INPUT:
        n -- an integer >= 2.

    EXAMPLES:
        sage: random_prime(100000)
        54601
        sage: random_prime(2)
        2


    AUTHOR:
        -- Jon Hanke: 2006-08-08  (with standard Stein cleanup)
    """
    import random    # since we don't want random to get
                     # pulled when you say "from sage.arith import *".
    n = sage.rings.integer.Integer(n)
    if n < 2:
        raise ValueError, "n must be >= 2."
    elif n == 2:
        return sage.rings.integer.Integer(n)
    else:
        previous_prime(random.randint(3,n))

    return previous_prime(random.randint(3,n))



def divisors(n):
    """
    Returns a list of all positive integer divisors
    of the nonzero integer n.

    EXAMPLES:
        sage: divisors(-3)
        [1, 3]
        sage: divisors(6)
        [1, 2, 3, 6]
        sage: divisors(28)
        [1, 2, 4, 7, 14, 28]
        sage: divisors(2^5)
        [1, 2, 4, 8, 16, 32]
        sage: divisors(100)
        [1, 2, 4, 5, 10, 20, 25, 50, 100]
        sage: divisors(1)
        [1]
        sage: divisors(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be nonzero
        sage: divisors(2^3 * 3^2 * 17)
        [1, 2, 3, 4, 6, 8, 9, 12, 17, 18, 24, 34, 36, 51, 68, 72, 102, 136, 153, 204, 306, 408, 612, 1224]
    """
    if n == 0:
        raise ValueError, "n must be nonzero"
    if n < 0:  # make positive, or, since (-1,1) is a factor, will get wrong answer.
        n*=-1
    F = factor(n)
    r = [0 for i in range(len(F))]
    e = [m for _, m in F]     # max exponents
    ans = []
    x = 1
    while r != e:
        ans.append(sage.rings.integer.Integer(x))
        r[0] += 1
        if r[0] <= e[0]:
            x *= F[0][0]
        else:  # carry
            i = 0
            while i < len(F) and r[i] > e[i]:
                x /= F[i][0]**F[i][1]
                r[i] = 0
                i += 1
                if i < len(F):
                    r[i] += 1
                    if r[i] <= e[i]:
                        x *= F[i][0]
        #endif
    #endwhile
    ans.append(sage.rings.integer.Integer(n))
    ans.sort()
    return ans

def sigma(n, k=1):
    """
    Return the sum of the k-th powers of the divisors of n.

    INPUT:
        n -- integer
        k -- integer (default: 1)

    OUTPUT:
        integer

    EXAMPLES:
        sage: sigma(5)
        6
        sage: sigma(5,2)
        26
    """
    n = sage.rings.integer.Integer(n)
    k = sage.rings.integer.Integer(k)
    return sum([d**k for d in divisors(n)])

def gcd(a, b=0, integer=False):
    """
    The greatest commond divisor of a and b.

    INPUT:
        a -- number
        b -- number (optional)
        integer -- (default: False); if True, do an integer GCD
    or
        v -- vector
        integer -- (default: False); if True, do an integer GCD
            NOTE -- this is *vastly* faster than doing the generic GCD

    EXAMPLES:
        sage: GCD(97,100)
        1
        sage: GCD(97*10^15, 19^20*97^2)
        97
        sage: GCD(2/3, 4/3)
        2/3
        sage: GCD([2,4,6,8])
        2
        sage: GCD(srange(0,10000,10))  # fast  !!
        10
    """
    if integer:
        if isinstance(a,list):
            return sage.rings.integer.GCD_list(a)
        else:
            return sage.rings.integer.Integer(a).gcd(\
                sage.rings.integer.Integer(b))
    if isinstance(a,list):
        return __GCD_list(a)
    if not isinstance(a, RingElement):
        a = sage.rings.integer.Integer(a)
    return a.gcd(b)

GCD = gcd

def lcm(a, b=None, integer=False):
    """
    The least common multiple of a and b, or if a is a list and b is
    omitted the least common multiple of all elements of v.

    NOTE: Use integer=True to make this vastly faster if you are
    working with lists of integers.

    INPUT:
        a -- number
        b -- number (optional)
        integer -- (default: False); if True, do an integer LCM
    or
        v -- vector
        integer -- (default: False); if True, do an integer LCM
            NOTE -- this is *vastly* faster than doing the generic LCM

    EXAMPLES:
        sage: LCM(97,100)
        9700
        sage: LCM(0,2)
        0
        sage: LCM(-3,-5)
        15
        sage: LCM([1,2,3,4,5/3])
        60
        sage: v = LCM(range(1,10000),integer=True)   # *very* fast!
        sage: len(str(v))
        4349
    """
    if integer:
        if isinstance(a,list):
            return sage.rings.integer.LCM_list(a)
        else:
            return sage.rings.integer.Integer(a).lcm(\
                sage.rings.integer.Integer(b))
    if isinstance(a, list):
        return __LCM_list(a)
    if not isinstance(a, RingElement):
        a = sage.rings.integer.Integer(a)
    return a.lcm(b)

LCM = lcm

def __LCM_list(v):
    if len(v) == 0:
        return 1
    x = v[0]
    for i in range(1,len(v)):
        x = LCM(x, v[i])
    return x

## def GCD_python(a, b=0):
##     """This function should behave exactly the same as GCD,
##     but is implemented in pure python."""
##     if isinstance(a,list):
##         return __GCD_list(a)
##     if a == 0:
##         return abs(b)
##     if b == 0:
##         return abs(a)
##     if a < 0:
##         a = -a
##     if b < 0:
##         b = -b
##     while b != 0:
##         c = a % b
##         a = b; b = c
##     return a

def __GCD_list(v):
    if len(v) == 0:
        return 1
    if len(v) == 1:
        return v[0]
    g = v[0]
    for i in range(1,len(v)):
        g = GCD(g, v[i])
    return g


def xgcd(a, b):
    """
    Returns triple of integers (g,s,t) such that g = s*a+t*b =
    gcd(a,b).

        sage: xgcd(56, 44)
        (4, 4, -5)
        sage: 4*56 + (-5)*44
        4
    """
    if not isinstance(a, RingElement):
        a = sage.rings.integer.Integer(a)
    return a.xgcd(b)

XGCD = xgcd

## def XGCD_python(a, b):
##     """
##     Returns triple (g,p,q) such that g = p*a+b*q = GCD(a,b).
##     This function should behave exactly the same as XGCD,
##     but is implemented in pure python.
##     """
##     if a == 0 and b == 0:
##         return (0,0,1)
##     if a == 0:
##         return (abs(b), 0, b/abs(b))
##     if b == 0:
##         return (abs(a), a/abs(a), 0)
##     psign = 1
##     qsign = 1
##     if a < 0:
##         a = -a
##         psign = -1
##     if b < 0:
##         b = -b
##         qsign = -1
##     p = 1; q = 0; r = 0; s = 1
##     while b != 0:
##         c = a % b
##         quot = a/b
##         a = b; b = c
##         new_r = p - quot*r
##         new_s = q - quot*s
##         p = r; q = s
##         r = new_r; s = new_s
##     return (a, p*psign, q*qsign)

def inverse_mod(a, m):
    """
    The inverse of the integer a modulo the integer m.
    sage: inverse_mod(7,1)
    0
    sage: inverse_mod(5,14)
    3
    sage: inverse_mod(3,-5)
    2
    """
    if m<0:
        m *= -1
    if m==1:
        return 0
    return sage.rings.integer.Integer((~(pari(a).Mod(m))).lift())

# def sqrt_mod(a, m):
#     """A square root of a modulo m."""

# def xxx_inverse_mod(a, m):
#     """The inverse of a modulo m."""
#     g,s,t = XGCD(a,m)
#     if g != 1:
#         raise "inverse_mod(a=%s,m=%s), error since GCD=%s"%(a,m,g)
#     return s

def power_mod(a,m,n):
    """The m-th power of a modulo the integer n.
    sage: power_mod(0,0,5)
    1
    sage: power_mod(2,390,391)
    285
    sage: power_mod(2,-1,7)
    4
    """
    if n==0:
        raise ZeroDivisionError, "Modulus must be nonzero."
    if n==1:
        return 0
    if m < 0:
        ainv = inverse_mod(a,n)
        return power_mod(ainv, -m, n)
    if m==0:
        return 1
    power = 1
    i = 0
    apow2 = a
    while ((m>>i) > 0):
        if((m>>i) & 1):
            power = (power * apow2) % n
        apow2 = (apow2 * apow2) % n
        i += 1
    return power

def generic_power(a, m, one=1):
    """
    The m-th power of a, where m is a non-negative
    integer and a is a Python object on which
    multiplication is defined.  The exponentiation
    is done using the standard binary powering algorithm.

    EXAMPLES:
        sage: generic_power(2,5)
        32
        sage: generic_power(RealField()('2.5'),4)
        39.0625000000000
        sage: generic_power(0,0)
        1
        sage: generic_power(2,-3)
        1/8
    """
    if a == one:
        return a
    if m < 0:
        a = ~a
        m = -m
    if m == 0:
        return one
    power = None
    i = 0
    apow2 = a
    while (m>>i) > 0:
        if (m>>i) & 1:
            if power is None:
                power = apow2
            else:
                power *= apow2
        apow2 *= apow2
        i += 1
    if power is None:
        return one
    return power


def rational_reconstruction(a, m, algorithm='fast'):
    """
    This function tries to compute x/y, where x/y is rational number
    is lowest terms such that reduction of x/y modulo m is equal to a
    and the absolute values of x and y are both <= sqrt(m/2).  If such
    x/y exists, that pair is unique and this function returns it.  If
    no such pair exists, this function raises ZeroDivisionError.

    An efficient algorithm for computing rational reconstruction is
    very similar to the extended Euclidean algorithm.  For more
    details, see Knuth, Vol 2, 3rd ed, pages 656-657.

     Input:  a -- an integer
             m -- a modulus
             algorithm -- (default: 'fast')
                  fast -- a fast compiled implementation
                  python -- a slow pure python implementation

     Output: Numerator and denominator n, d of the unique rational
             number r=n/d, if it exists, with
             |n| and |d| <= sqrt(N/2).
             Return (0,0) if no such number exists.

    The algorithm for rational reconstruction is described (with a
    complete nontrivial proof) on pages 656-657 of Knuth, Vol 2, 3rd
    ed. as the solution to exercise 51 on page 379.  See in particular
    the conclusion paragraph right in the middle of page 657, which
    describes the algorithm thus:

    This discussion proves that the problem can be solved efficiently by
    applying Algorithm 4.5.2X with u=m and v=a, but with the following
    replacement for step X2: If v3<=sqrt(m/2), the algorithm terminates.
    The pair (x,y)=(|v2|,v3*sign(v2)) is then the unique solution,
    provided that x and y are coprime and x<=sqrt(m/2); otherwise there is
    no solution.   (Alg 4.5.2X is the extended Euclidean algorithm.)

    Knuth remarks that this algorithm is due to Wang, Kornerup, and
    Gregory from around 1983.


    EXAMPLES::
        sage: m = 100000
        sage: (119*inverse_mod(53,m))%m
        11323
        sage: rational_reconstruction(11323,m)
        119/53

        sage: rational_reconstruction(400,1000)
        Traceback (most recent call last):
        ...
        ValueError: Rational reconstruction of 400 (mod 1000) does not exist.

        sage: rational_reconstruction(3,292393, algorithm='python')
        3
        sage: a = Integers(292393)(45/97); a
        204977
        sage: rational_reconstruction(a,292393, algorithm='python')
        45/97
        sage: a = Integers(292393)(45/97); a
        204977
        sage: rational_reconstruction(a,292393, algorithm='fast')
        45/97
        sage: rational_reconstruction(293048,292393, algorithm='fast')
        Traceback (most recent call last):
        ...
        ValueError: Rational reconstruction of 655 (mod 292393) does not exist.
        sage: rational_reconstruction(293048,292393, algorithm='python')
        Traceback (most recent call last):
        ...
        ValueError: Rational reconstruction of 655 (mod 292393) does not exist.
    """
    if algorithm == 'fast':
        return sage.rings.integer.Integer(a).rational_reconstruction(m)
    elif algorithm == 'python':
        return _rational_reconstruction_python(a,m)
    else:
        raise ValueError, "unknown algorithm"

def _rational_reconstruction_python(a,m):
    a = int(a); m = int(m)
    a %= m
    if a == 0 or m==0:
        return sage.rings.integer.Integer(0)/sage.rings.integer.Integer(1)
    if m < 0:
        m = -m
    if a < 0:
        a = m-a
    if a == 1:
        return sage.rings.integer.Integer(1)/sage.rings.integer.Integer(1)
    u = m
    v = a
    bnd = math.sqrt(m/2)
    U = (1,0,u)
    V = (0,1,v)
    while abs(V[2]) > bnd:
        q = U[2]/V[2]  # floor is implicit
        T = (U[0]-q*V[0], U[1]-q*V[1], U[2]-q*V[2])
        U = V
        V = T
    x = abs(V[1])
    y = V[2]
    if V[1] < 0:
        y *= -1
    if x <= bnd and GCD(x,y) == 1:
        return sage.rings.integer.Integer(y) / sage.rings.integer.Integer(x)
    raise ValueError, "Rational reconstruction of %s (mod %s) does not exist."%(a,m)

def mqrr_rational_reconstruction(u, m, T):
    """
    Maximal Quotient Rational Reconstruction.

    FOR research purposes only -- this is pure Python, so slow.

    Input:
        u, m, and T are integers and
        m > u>=0, T>0.
    Output:
        Either integer n,d such that d>0, gcd(n,d)=1, n/d=u (mod m),
        and T*abs(n)*d < m, or None.

    Reference: Monagan, Maximal Quotient Rational Reconstruction: An
               Almost Optimal Algorithm for Rational Reconstruction (page 11)

    This algorithm is probabilistic.
    """
    if u == 0:
        if m > T:
            return (0,1)
        else:
            return None
    n, d = 0, 0
    t0, r0 = 0, m
    t1, r1 = 1, u
    while r1 != 0 and r0 > T:
        q = r0/r1   # C division implicit floor
        if q > T:
            n, d, T = r1, t1, q
        r0, r1 = r1, r0 - q*r1
        t0, t1 = t1, t0 - q*t1
    if d != 0 and GCD(n,d) == 1:
        return (n,d)
    return None


######################


def trial_division(n, bound=None):
    """
    Return the smallest prime divisor <= bound of the positive integer
    n, or n if there is no such prime.  If the optional argument bound
    is omitted, then bound=n.

    INPUT:
        n -- a positive integer
        bound - (optional) a positive integer

    OUTPUT:
        int -- a prime p<=bound that divides n, or n if
               there is no such prime.

    EXAMPLES:
        sage: trial_division(15)
        3
        sage: trial_division(91)
        7
        sage: trial_division(11)
        11
        sage: trial_division(387833, 300)
        387833
        sage: # 300 is not big enough to split off a
        sage: # factor, but 400 is.
        sage: trial_division(387833, 400)
        389
    """
    if n == 1: return 1
    for p in [2, 3, 5]:
        if n%p == 0: return p
    if bound == None: bound = n
    dif = [6, 4, 2, 4, 2, 4, 6, 2]
    m = 7; i = 1
    while m <= bound and m*m <= n:
        if n%m == 0:
            return m
        m += dif[i%8]
        i += 1
    return n

def __factor_using_trial_division(n):
    """
    Returns the factorization of the integer n as
    a sorted list of tuples (p,e).

    INPUT:
        n -- an integer
    OUTPUT:
        list -- factorization of n
    """
    if n in [-1, 0, 1]: return []
    if n < 0: n = -n
    F = []
    while n != 1:
        p = trial_division(n)
        e = 1
        n /= p
        while n%p == 0:
            e += 1; n /= p
        F.append((p,e))
    F.sort()
    return F

def __factor_using_pari(n, int_=False, debug_level=0):
    if int_:
        Z = int
    else:
        import sage.rings.integer_ring
        Z = sage.rings.integer_ring.IntegerRing()
    prev = pari.get_debug_level()
    pari.set_debug_level(debug_level)
    F = pari(n).factor()
    B = F[0]
    e = F[1]
    v = [(Z(B[i]),Z(e[i])) for i in xrange(len(B))]
    if debug_level > 0:
        pari.set_debug_level(prev)
    return v


#todo: add a limit option to factor, so it will only split off
# primes at most a given limit.

def factor(n, proof=True, int_=False, algorithm='pari', verbose=0):
    """
    Returns the factorization of the integer n as a sorted list of
    tuples (p,e).

    INPUT:
        n -- an nonzero integer
        proof -- bool (default: True)
        int_ -- bool (default: False) whether to return answers as Python ints
        algorithm -- string
                 * 'pari' -- (default)  use the PARI c library
                 * 'kash' -- use KASH computer algebra system (requires
                             the optional kash package be installed)
        verbose -- integer (default 0); pari's debug variable is set to this;
                   e.g., set to 4 or 8 to see lots of output during factorization.
    OUTPUT:
        factorization of n

    NOTES:
        The qsieve and ecm commands give access to highly optimized
        implementations of algorithms for doing certain integer
        factorization problems.  These implementation are not used by
        the generic factor command, which currently just calls PARI
        (note that PARI also implements sieve and ecm algorithms, but
        they aren't as optimized).  Thus you might consider using them
        instead for certain numbers.

    EXAMPLES:
        sage: factor(500)
        2^2 * 5^3
        sage: factor(-20)
        -1 * 2^2 * 5

        sage: factor(500, algorithm='kash')     # requires optional kash package
        2^2 * 5^3

        sage: factor(0)
        Traceback (most recent call last):
        ...
        ArithmeticError: Prime factorization of 0 not defined.
        sage: factor(1)
        1
        sage: factor(-1)
        -1
        sage: factor(2004)
        2^2 * 3 * 167

        sage: factor(2^197 + 1)       # takes a long time
        3 * 197002597249 * 1348959352853811313 * 251951573867253012259144010843
    """
    Z = sage.rings.integer.Integer
    if not isinstance(n, (int,long, Z)):
        try:
            return n.factor()
        except AttributeError:
            raise TypeError, "unable to factor n"
    #n = abs(n)
    n = Z(n)
    if n < 0:
        unit = Z(-1)
        n = -n
    else:
        unit = Z(1)

    if n == 0:
        raise ArithmeticError, "Prime factorization of 0 not defined."
    if n == 1:
        return factorization.Factorization([], unit)
    #if n < 10000000000: return __factor_using_trial_division(n)
    if algorithm == 'pari':
        return factorization.Factorization(__factor_using_pari(n,
                                   int_=int_, debug_level=verbose), unit)
    elif algorithm == 'kash':
        from sage.interfaces.all import kash
        F = kash.eval('Factorization(%s)'%n)
        i = F.rfind(']') + 1
        F = F[:i]
        F = F.replace("<","(").replace(">",")")
        F = eval(F)
        if not int_:
            F = [(Z(a), Z(b)) for a,b in F]
        return factorization.Factorization(F, unit)
    else:
        raise ValueError, "Algorithm is not known"


def prime_divisors(n):
    """
    The prime divisors of the integer n, sorted in increasing order.  If n
    is negative, we do *not* include -1 among the prime divisors, since -1 is
    not a prime number.

    sage: prime_divisors(1)
    []
    sage: prime_divisors(100)
    [2, 5]
    sage: prime_divisors(-100)
    [2, 5]
    sage: prime_divisors(2004)
    [2, 3, 167]
    """
    v = [p for p,_ in factor(n) if p != -1]
    v.sort()
    return v

prime_factors = prime_divisors

def odd_part(n):
    """
    The odd part of the integer $n$.  This is $n / 2^v$,
    where $v =$ \code{valuation(n,2)}.
    """
    n = sage.rings.integer.Integer(n)
    v = valuation(n,2)
    return n / (2**v)


def prime_to_m_part(n,m):
    """
    Returns the prime-to-m part of n, i.e., the largest divisor
    of n that is coprime to m.

    INPUT:
        n -- Integer (nonzero)
        m -- Integer
    OUTPUT:
        Integer
    """
    if n == 0:
        raise ValueError, "n must be nonzero."
    if m == 0:
        return sage.rings.integer.Integer(1)
    n = sage.rings.integer.Integer(n); m = sage.rings.integer.Integer(m)
    while True:
        g = gcd(n,m)
        if g == 1:
            return n
        n = n // g


def is_square(n, root=False):
    """
    Returns whether or not n is square, and if n is a square
    also returns the square root.  If n is not square, also
    returns None.
    INPUT:
        n -- an integer
        root -- whether or not to also return a square root (default: False)
    OUTPUT:
        bool -- whether or not a square
        object --
    """
    t, x = pari(n).issquare(find_root=True)
    if root:
        if t:
            if hasattr(n, 'parent'):
                x = n.parent()(str(x))
            else:
                x = x.python()
        return t, x
    return t


def is_squarefree(n):
    """
    Returns True if and only if n is not divisible by the square of an integer > 1.
    """
    if n==0:
        return False
    for p, r in factor(n):
        if r>1:
            return False
    return True


#################################################################
# Euler phi function
#################################################################
def euler_phi(n):
    """
    Return the value of the Euler phi function on the integer n.  We
    defined this to be the number of positive integers <= n that are
    relatively prime to n.  Thus if n<=0 then \code{euler_phi(n)} is
    defined and equals 0.

    INPUT:
        n -- an integer

    EXAMPLES:

        sage: euler_phi(1)
        1
        sage: euler_phi(2)
        1
        sage: euler_phi(3)
        2
        sage: euler_phi(12)
        4
        sage: euler_phi(37)
        36

    Notice that euler_phi is defined to be 0 on negative numbers and 0.

        sage: euler_phi(-1)
        0
        sage: euler_phi(0)
        0

    We verify directly that the phi function is correct for 21.

       sage: euler_phi(21)
       12
       sage: [i for i in range(21) if gcd(21,i) == 1]
       [1, 2, 4, 5, 8, 10, 11, 13, 16, 17, 19, 20]

    The length of the list of integers 'i' in range(n) such that
    the gcd(i,n) == 1 equals euler_phi(n).

       sage: len([i for i in range(21) if gcd(21,i) == 1]) == euler_phi(21)
       True

    AUTHORS:
        - William Stein
        - Alex Clemesha (2006-01-10): some examples
    """
    if n<=0:
        return 0
    if n<=2:
        return 1
    return sage.rings.integer.Integer(pari(n).phi())
    #return misc.mul([(p-1)*p**(r-1) for p, r in factor(n)])

def crt(a,b=0,m=1,n=1):
    """
    Use the Chinese Remainder Theorem to find some integer x such
    that x=a (mod m) and x=b (mod n).   Note that x is only well-defined
    modulo m*n.

    sage: crt(2, 1, 3, 5)
    -4
    sage: crt(13,20,100,301)
    -2087
    """
    if isinstance(a,list):
        return CRT_list(a,b)
    g, alpha, beta = XGCD(m,n)
    if g != 1:
        raise ValueError, "arguments a and b must be coprime"
    return a+(b-a)*alpha*m

CRT = crt

def CRT_list(v, moduli):
    if len(v) == 0:
        return 0
    x = v[0]
    m = moduli[0]
    for i in range(1,len(v)):
        x = CRT(x,v[i],m,moduli[i])
        m *= moduli[i]
    return x%m

def CRT_basis(moduli):
    """
    Return a list of integers a[i] such that CRT to the given moduli
    of numbers x[0],...,x[n-1] is a[0]*x0 + ... + a[n-1]*x[n-1].

    INPUT:
        list -- list of integers
    """
    n = len(moduli)
    if n == 0:
        return []
    v = [0 for _ in range(n)]
    a = list(v)  # copy
    for i in range(n):
        v[i] = 1
        a[i] = CRT_list(v, moduli)
        v[i] = 0
    return a

def CRT_vectors(X, moduli):
    """
    INPUT:
        X -- list of lists of the same length
        moduli -- list of len(X) moduli
    OUTPUT:
        list -- application of CRT componentwise.
    """
    # First find the CRT basis:
    if len(X) == 0 or len(X[0]) == 0:
        return []
    n = len(X)
    if n != len(moduli):
        raise ValueError, "number of moduli must equal length of X"
    a = CRT_basis(moduli)
    modulus = misc.prod(moduli)
    return [sum([a[i]*X[i][j] for i in range(n)]) % modulus for j in range(len(X[0]))]

def binomial(x,m):
    r"""
    Return the binomial coefficient
    $$
       x (x-1) \cdots (x-m+1) / m!
    $$
    which is defined for $m \in \Z$ and any $x$.
    If $m<0$ return $0$.

    INPUT::
        x -- number
        m -- integer

    OUTPUT::
        number

    EXAMPLES::
        sage: binomial(5,2)
        10
        sage: binomial(2,0)
        1
        sage: binomial(1/2, 0)
        1
        sage: binomial(3,-1)
        0
        sage: binomial(20,10)
        184756
        sage: binomial(RealField()('2.5'), 2)
        1.87500000000000
    """
    if not isinstance(m, (int, long, sage.rings.integer.Integer)):
        raise TypeError, 'm must be an integer'
    if isinstance(x, (int, long, sage.rings.integer.Integer)):
        return sage.rings.integer.Integer(pari(x).binomial(m))
    try:
        P = x.parent()
    except AttributeError:
        P = type(x)
    if m < 0:
        return P(0)
    return misc.prod([x-i for i in xrange(m)]) / P(factorial(m))

def gaussian_binomial(n,k,q):
    r"""
    Return the gaussian binomial
    $$
       \binom{n}{k}_q = \frac{(1-q^m)(1-q^{m-1})\cdots (1-q^{m-r+1})}
                             {(1-q)(1-q^2)\cdots (1-q^r)}.
    $$

    EXAMPLES:
        sage: gaussian_binomial(5,1,2)
        31

    AUTHOR: David Joyner and William Stein
    """
    n = sage.rings.integer.Integer(misc.prod([1 - q**i for i in range((n-k+1),n+1)]))
    d = sage.rings.integer.Integer(misc.prod([1 - q**i for i in range(1,k+1)]))
    return n / d

def kronecker_symbol(x,y):
    """
    The Kronecker symbol (x|y).

    INPUT:
        x -- integer
        y -- integer

    EXAMPLES:
        sage: kronecker(3,5)
        -1
        sage: kronecker(3,15)
        0
        sage: kronecker(2,15)
        1
        sage: kronecker(-2,15)
        -1

    IMPLEMENTATION: Using Pari.
    """
    return sage.rings.integer.Integer(pari(x).kronecker(y).python())

def kronecker(x,y):
    r"""
    Synonym for \code{kronecker_symbol}.
    """
    return kronecker_symbol(x,y)

def primitive_root(n):
    """
    Return a generator for the multiplicative group of integers
    modulo $n$, if one exists.

    EXAMPLES:
        sage: primitive_root(23)
        5
        sage: print [primitive_root(p) for p in primes(100)]
        [1, 2, 2, 3, 2, 2, 3, 2, 5, 2, 3, 2, 6, 3, 5, 2, 2, 2, 2, 7, 5, 3, 2, 3, 5]
    """
    Z = sage.rings.integer.Integer
    try:
        return Z(pari(Z(n)).znprimroot())
    except RuntimeError:
        raise ArithmeticError, "There is no primitive root modulo n"

def discrete_log_generic(b, a, ord=None):
    """
    Return an integer $n$ such that $b^n = a$, assuming that ord is a
    multiple of the multiplicative order of $a$.  If ord is not
    specified an attempt is made to compute it.

    The elements a and b must support exponentiation to a negative
    power.

    If no such $x$ exits, this function raises a ValueError exception.

    ALGORITHM: Baby step giant step.

    EXAMPLES:
        sage: b = Mod(2,37);  a = b^20
        sage: discrete_log_generic(b, a)
        20
        sage: b = Mod(2,997);  a = b^20
        sage: discrete_log_generic(b, a)
        20

        sage: K = GF(3^6,'b')
        sage: b = K.gen()
        sage: a = b^210
        sage: discrete_log_generic(b, a, K.order()-1)
        210

        sage: b = Mod(1,37);  a = Mod(2,37)
        sage: discrete_log_generic(b, a)
        Traceback (most recent call last):
        ...
        ValueError: Log of a to the base b does not exist.
        sage: b = Mod(1,997);  a = Mod(2,997)
        sage: discrete_log_generic(b, a)
        Traceback (most recent call last):
        ...
        ValueError: Log of a to the base b does not exist.

    AUTHOR: William Stein and David Joyner (2005-01-05)
    """
    Z = sage.rings.integer.Integer

    if b == 0:
        if a == 0:
            return Integer(1)
        else:
            raise ValueError, "Log of a to the base b does not exist."
    elif a == 0:
        if b == 0:
            return Integer(1)
        else:
            raise ValueError, "Log of a to the base b does not exist."

    if ord is None:
        ord = b.multiplicative_order()
    ord = Z(ord)
    if ord < 100:
        c = 1
        for i in range(ord):
            if c == a:        # is b^i
                return Z(i)
            c *= b
        raise ValueError, "Log of a to the base b does not exist."

    m = ord.isqrt()
    g = [a]
    c = b**(-m)
    S2 = [1]
    for i in range(m):
        g.append(g[i]*c)
        if i < m-1:
            S2.append(S2[i]*b)
    for y in g:
        if y in S2:
            x = S2.index(y)
            return Z(m*(g.index(y)) + x)

    raise ValueError, "Log does not exist."



def quadratic_residues(n):
    r"""
    Return a sorted list of all squares modulo the integer $n$ in the
    range $0\leq x < |n|$.

    EXAMPLES:
        sage: quadratic_residues(11)
        [0, 1, 3, 4, 5, 9]
        sage: quadratic_residues(1)
        [0]
        sage: quadratic_residues(2)
        [0, 1]
        sage: quadratic_residues(8)
        [0, 1, 4]
        sage: quadratic_residues(-10)
        [0, 1, 4, 5, 6, 9]
        sage: v = quadratic_residues(1000); len(v);
        159
    """
    n = abs(int(n))
    Z = sage.rings.integer.Integer
    X = list(set([Z((a*a)%n) for a in range(n/2+1)]))
    X.sort()
    return X

## This much slower than above, for obvious reasons.
## def quadratic_residues2(p):
##     return [x for x in range(p-1) if kronecker_symbol(x,p)==1]

def Max(x):
    """
    The maximum of a list of objects, on which a binary max operation
    is defined.
    """
    assert isinstance(x,list), "Argument must be a list."
    if len(x)==0:
        return 0
    m=x[0]
    for i in range(1,len(x)):
        m=max(m,x[i])
    return m

def Min(x):
    """
    The minimum of a list of objects, on which a binary min operation
    is defined.
    """
    assert isinstance(x,list), "Argument must be a list."
    if len(x)==0:
        return 0
    m=x[0]
    for i in range(1,len(x)):
        m=min(m,x[i])
    return m

def moebius(n):
    r"""
    Returns the value of the Moebius function of abs(n), where n is an integer.

    DEFINITION:
        $\mu(n)$ is 0 if $n$ is not square free, and otherwise equals $(-1)^r$,
        where $n$ has $r$ distinct prime factors.

    INPUT:
        n -- an integer
    OUTPUT:
        0, 1, or -1
    EXAMPLES:
        sage: moebius(-5)
        -1
        sage: moebius(9)
        0
        sage: moebius(12)
        0
        sage: moebius(-35)
        1
        sage: moebius(-1)
        1
        sage: moebius(7)
        -1
    """
    if n < 0:
        n = -n
    F = factor(n)
    for _, e in F:
        if e >= 2:
            return 0
    return (-1)**len(F)

def farey(v, lim):
    """
    Return the Farey sequence associated to the floating point
    number v.

    INPUT:
       v -- float (automatically converted to a float)
       lim --  maximum denominator.
    OUTPUT:
       Results are (numerator, denominator); (1, 0) is"infinity".

    AUTHOR: Scott David Daniels, Python Cookbook, 2nd Ed., Recipe 18.13
    """
    v = float(v)
    if v < 0:
        n, d = farey(-v, lim)
        return -n, d
    z = lim - lim    # Get a "0 of the right type" for denominator
    lower, upper = (z, z+1), (z+1, z)
    while True:
        mediant = (lower[0] + upper[0]), (lower[1] + upper[1])
        if v * mediant[1] > mediant[0]:
            if lim < mediant[1]:
                return upper
            lower = mediant
        elif v * mediant[1] == mediant[0]:
            if lim >= mediant[1]:
                return mediant
            if lower[1] < upper[1]:
                return lower
            return upper
        else:
            if lim < mediant[1]:
                return lower
            upper = mediant

def number_of_partitions(n):
    """
    Return the number of partitions of the integer $n$.

    To compute all the partitions of $n$ use \code{partitions(n)}.

    EXAMPLES:
        sage: number_of_partitions(3)
        3
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(40)
        37338
        sage: number_of_partitions(100)
        190569292
        sage: number_of_partitions(-5)
        0
        sage: number_of_partitions(0)
        1
    """
    ZZ = sage.rings.integer.Integer
    return ZZ(pari(ZZ(n)).numbpart())

def partitions(n):
    """
    Generator of all the partitions of the integer $n$.

    To compute the number of partitions of $n$ use
    \code{number_of_partitions(n)}.

    INPUT:
        n -- int

    EXAMPLES:
        >> partitions(3)
        <generator object at 0xab3b3eac>
        sage: list(partitions(3))
        [(1, 1, 1), (1, 2), (3,)]

    AUTHOR: David Eppstein, Jan Van lent, George Yoshida; Python Cookbook 2, Recipe 19.16.
    """
    n == sage.rings.integer.Integer(n)
    # base case of the recursion: zero is the sum of the empty tuple
    if n == 0:
        yield ( )
        return
    # modify the partitions of n-1 to form the partitions of n
    for p in partitions(n-1):
        yield (1,) + p
        if p and (len(p) < 2 or p[1] > p[0]):
            yield (p[0] + 1,) + p[1:]

def continued_fraction(x, partial_convergents=False):
    r"""
    Returns the continued fraction of x.

    \begin{note}
    This may be slow since it's implemented in pure
    Python for real input.  For rational number input the PARI C
    library is used.
    \end{note}

    EXAMPLES:
        sage: continued_fraction(45/17)
        [2, 1, 1, 1, 5]
        sage: continued_fraction(sqrt(2))
        [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1]
        sage: continued_fraction(RR(pi), partial_convergents=True)
        ([3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3],
         [(3, 1),
          (22, 7),
          (333, 106),
          (355, 113),
          (103993, 33102),
          (104348, 33215),
          (208341, 66317),
          (312689, 99532),
          (833719, 265381),
          (1146408, 364913),
          (4272943, 1360120),
          (5419351, 1725033),
          (80143857, 25510582),
          (245850922, 78256779)])
        sage: continued_fraction(e)
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 11]
        sage: continued_fraction(RR(e))
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 11]
        sage: print continued_fraction(RealField(200)(e))
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 14, 1, 1, 16, 1, 1, 18, 1, 1, 20, 1, 1, 22, 1, 1, 24, 1, 1, 26, 1, 1, 28, 1, 1, 30, 1, 1, 32, 1, 1, 34, 1, 1, 36, 1, 1, 38, 1, 1]
    """
    if isinstance(x, (sage.rings.integer.Integer, sage.rings.rational.Rational,
                      int, long)):
        return pari(x).contfrac().python()
    x_in = x
    v = []
    w = [(0,1), (1,0)] # keep track of convergents
    start = x
    i = 0
    try:
        while True:
            i += 1
            a = sage.rings.integer.Integer(int(x.floor()))
            v.append(a)
            n = len(v)-1
            pn = v[n]*w[n+1][0] + w[n][0]
            qn = v[n]*w[n+1][1] + w[n][1]
            w.append((pn, qn))
            x -= a
            if abs(start - pn/qn) == 0:
                del w[0]; del w[0]
                if partial_convergents:
                    return v, w
                else:
                    return v
            x = 1/x
    except (AttributeError, NotImplementedError, TypeError), msg:
        raise NotImplementedError, "%s\ncomputation of continued fraction of x not implemented; try computing continued fraction of RR(x) instead."%msg

def convergent(v, n):
    """
    Return the n-th continued fraction convergent of the continued
    fraction defined by the sequence of integers v.  We assume
    $n \geq 0$.

    INPUT:
        v -- list of integers
        n -- integer
    OUTPUT:
        a rational number

    If the continued fraction integers are
    $$
            v = [a_0, a_1, a_2, \ldots, a_k]
    $$
    then \code{convergent(v,2)} is the rational number
    $$
              a_0 + 1/a_1
    $$
    and \code{convergent(v,k)} is the rational number
    $$
            a1 + 1/(a2+1/(...) ... )
    $$
    represented by the continued fraction.

    EXAMPLES:
        sage: convergent([2, 1, 2, 1, 1, 4, 1, 1], 7)
        193/71
    """
    Q = sage.rings.rational_field.RationalField()
    i = int(n)
    x = Q(v[i])
    i -= 1
    while i >= 0:
        x = Q(v[i]) + 1/x
        i -= 1
    return x



## def convergents_pnqn(x):
##     """
##     Return the pairs (pn,qn) that are the numerators and denominators
##     of the partial convergents of the continued fraction of x.  We
##     include (0,1) and (1,0) at the beginning of the list (these are
##     the -2 and -1 th convergents).
##     """
##     v = pari(x).contfrac()
##     w = [(0,1), (1,0)]
##     for n in range(len(v)):
##         pn = w[n+1][0]*v[n] + w[n][0]
##         qn = w[n+1][1]*v[n] + w[n][1]
##         w.append(int(pn), int(qn))
##     return w


def convergents(v):
    """
    Return all the partial convergents of a continued fraction
    defined by the sequence of integers v.

    If v is not a list, compute the continued fraction of v and return
    its convergents (this is potentially much faster than calling
    continued_fraction first, since continued fractions are
    implemented using PARI and there is overhead moving the answer
    back from PARI).

    INPUT:
        v -- list of integers or a rational number
    OUTPUT:
        list -- of partial convergents, as rational numbers

    EXAMPLES:
        sage: convergents([2, 1, 2, 1, 1, 4, 1, 1])
        [2, 3, 8/3, 11/4, 19/7, 87/32, 106/39, 193/71]
    """
    Q = sage.rings.rational_field.RationalField()
    if not isinstance(v, list):
        v = pari(v).contfrac()
    w = [(0,1), (1,0)]
    for n in range(len(v)):
        pn = w[n+1][0]*v[n] + w[n][0]
        qn = w[n+1][1]*v[n] + w[n][1]
        w.append((pn, qn))
    return [Q(x) for x in w[2:]]


def number_of_divisors(n):
    """
    Return the number of divisors of the integer n.
    """
    m = sage.rings.integer.Integer(n)
    if m.is_zero():
        raise ValueError, "input must be nonzero"
    return sage.rings.integer.Integer(pari(m).numdiv())



def hilbert_symbol(a, b, p, algorithm="pari"):
    """
    Returns 1 if $ax^2 + by^2$ $p$-adically represents a nonzero
    square, otherwise returns $-1$.  If either a or b is 0, returns 0.

    INPUT:
        a, b -- integers
        p -- integer; either prime or -1 (which represents the archimedean place)
        algorithm -- string
                   'pari' -- (default) use the PARI C library
                   'direct' -- use a Python implementation
                   'all' -- use both PARI and direct and check that
                            the results agree, then return the common answer
    OUTPUT:
        integer (0, -1, or 1)

    EXAMPLES:
        sage: hilbert_symbol (-1, -1, -1, algorithm='all')
        -1
        sage: hilbert_symbol (2,3, 5, algorithm='all')
        1
        sage: hilbert_symbol (4, 3, 5, algorithm='all')
        1
        sage: hilbert_symbol (0, 3, 5, algorithm='all')
        0
        sage: hilbert_symbol (-1, -1, 2, algorithm='all')
        -1
        sage: hilbert_symbol (1, -1, 2, algorithm='all')
        1
        sage: hilbert_symbol (3, -1, 2, algorithm='all')
        -1

    AUTHORS:
       -- William Stein and David Kohel (2006-01-05)
    """
    Integer = sage.rings.integer.Integer

    p = Integer(p)
    if p != -1 and not p.is_prime():
        raise ValueError, "p must be prime or -1"
    a = Integer(a)
    b = Integer(b)

    if algorithm == "pari":

        from sage.libs.all import pari
        return Integer(pari(a).hilbert(b,p))

    elif algorithm == 'direct':
        if a == 0 or b == 0:
            return Integer(0)

        p = Integer(p)
        one = Integer(1)

        if p != -1:
            p_sqr = p**2
            while a%p_sqr == 0: a //= p_sqr
            while b%p_sqr == 0: b //= p_sqr

        if p != 2 and True in ( kronecker(x,p) == 1 for x in (a,b,a+b) ):
            return one
        if a%p == 0:
            if b%p == 0:
                return hilbert_symbol(p,-(b//p),p)*hilbert_symbol(a//p,b,p)
            elif p == 2 and (b%4) == 3:
                if kronecker(a+b,p) == -1:
                    return -one
            elif kronecker(b,p) == -1:
                return -one
        elif b%p == 0:
            if p == 2 and (a%4) == 3:
                if kronecker(a+b,p) == -1:
                    return -one
            elif kronecker(a,p) == -1:
                return -one
        elif p == 2 and (a%4) == 3 and (b%4) == 3:
            return -one
        return one
    elif algorithm == 'all':
        ans_pari = hilbert_symbol(a,b,p,algorithm='pari')
        ans_direct = hilbert_symbol(a,b,p,algorithm='direct')
        if ans_pari != ans_direct:
            raise RuntimeError, "There is a bug in hilbert_symbol; two ways of computing the Hilbert symbol (%s,%s)_%s disagree"%(a,b,p)
        return ans_pari
    else:
        raise ValueError, "Algorithm %s not defined"%algorithm



##############################################################################
##  falling and rising factorials
##  By Jaap Spies
##
##       Copyright (C) 2006 Jaap Spies <j.spies@hccnet.nl>
##      Copyright (C) 2006 William Stein <wstein@gmail.com>
##
## Distributed under the terms of the GNU General Public License (GPL)
##                  http://www.gnu.org/licenses/
##############################################################################


def falling_factorial(x, a):
    r"""
    Returns the falling factorial $(x)_a$.

    The notation in the literature is a mess: often $(x)_a$, but there
    are many other notations: GKP: Concrete Mathematics uses
    $x^{\underline{a}}$.

    Definition: for integer $a \ge 0$ we have $x(x-1) \cdots (x-a+1)$.
    In all other cases we use the GAMMA-function:
    $\frac {\Gamma(x+1)} {\Gamma(x-a+1)}$.

    INPUT:
        x -- element of a ring
        a -- a non-negative integer
      or
        x and a -- any numbers

    OUTPUT:
        the falling factorial

    EXAMPLES:
        sage: falling_factorial(10, 3)
        720
        sage: falling_factorial(10, RR('3.0'))
        720.000000000000
        sage: falling_factorial(10, RR('3.3'))
        1310.11633396600
        sage: falling_factorial(10, 10)
        3628800
        sage: factorial(10)
        3628800
        sage: falling_factorial(1+I, I)
        0.652965496420166 + 0.343065839816545*I
        sage: falling_factorial(1+I, 4)
        2.00000000000000 + 4.00000000000000*I
        sage: falling_factorial(I, 4)
        -10.0000000000000

        sage: M = MatrixSpace(ZZ, 4, 4)
        sage: A = M([1,0,1,0,1,0,1,0,1,0,10,10,1,0,1,1])
        sage: falling_factorial(A, 2) # A(A - I)
        [  1   0  10  10]
        [  1   0  10  10]
        [ 20   0 101 100]
        [  2   0  11  10]

        sage: x = ZZ['x'].0
        sage: falling_factorial(x, 4)
        x^4 - 6*x^3 + 11*x^2 - 6*x

    AUTHOR:
        -- Jaap Spies (2006-03-05)
    """
    if isinstance(a, (sage.rings.integer.Integer, int, long)) and a >= 0:
        return misc.prod([(x - i) for i in range(a)])
    from sage.functions.transcendental import gamma
    return gamma(x+1) / gamma(x-a+1)

def rising_factorial(x, a):
    r"""
    Returns the rising factorial $(x)^a$.

    The notation in the literature is a mess: often $(x)^a$, but there
    are many other notations: GKP: Concrete Mathematics uses
    $x^{\overline{a}}$.

    The rising factorial is also known as the Pochhammer symbol, see
    Maple and Mathematica.

    Definition: for integer $a \ge 0$ we have $x(x+1) \cdots (x+a-1)$.
    In all other cases we use the GAMMA-function:
    $\frac {\Gamma(x+a)} {\Gamma(x)}$.

    INPUT:
        x -- element of a ring
        a -- a non-negative integer
      or
        x and a -- any numbers

    OUTPUT:
        the rising factorial

    EXAMPLES:
        sage: rising_factorial(10,3)
        1320

        sage: rising_factorial(10,RR('3.0'))
        1320.00000000000

        sage: rising_factorial(10,RR('3.3'))
        2826.38895824964

        sage: rising_factorial(1+I, I)
        0.266816390637832 + 0.122783354006371*I

        sage: rising_factorial(I, 4)
        -10.0000000000000

    See falling_factorial(I, 4).

        sage: R = ZZ['x']
        sage: rising_factorial(x, 4)
        x^4 + 6*x^3 + 11*x^2 + 6*x

    AUTHOR:
        -- Jaap Spies (2006-03-05)
    """
    if isinstance(a, (sage.rings.integer.Integer, int, long)) and a >= 0:
        return misc.prod([(x + i) for i in range(a)])
    from sage.functions.transcendental import gamma
    return gamma(x+a) / gamma(x)




def ceil(x):
    """
    Return the ceiling of x.
    """
    try:
        return sage.rings.all.Integer(x.ceil())
    except AttributeError:
        try:
            return sage.rings.all.Integer(int(math.ceil(float(x))))
        except TypeError:
            pass
    raise NotImplementedError, "computation of floor of %s not implemented"%x

ceiling = ceil

def floor(x):
    r"""
    Return the largest integer $\leq x$.

    INPUT:
        x -- an object that has a floor method or is coercible to int

    OUTPUT:
        an Integer

    EXAMPLES:
        sage: floor(5.4)
        5
        sage: floor(float(5.4))
        5
        sage: floor(-5/2)
        -3
        sage: floor(RDF(-5/2))
        -3
    """
    try:
        return sage.rings.all.Integer(x.floor())
    except AttributeError:
        try:
            return sage.rings.all.Integer(int(math.floor(float(x))))
        except TypeError:
            pass
    raise NotImplementedError, "computation of floor of %s not implemented"%x



def two_squares(n, algorithm='gap'):
    """
    Write the integer n as a sum of two integer squares if
    possible; otherwise raise a ValueError.

    EXAMPLES:
        sage: two_squares(389)
        (10, 17)
        sage: two_squares(7)
        Traceback (most recent call last):
        ...
        ValueError: 7 is not a sum of two squares
        sage: a,b = two_squares(2009); a,b
        (28, 35)
        sage: a^2 + b^2
        2009

    TODO: Create an implementation using PARI's continued
    fraction implementation.
    """
    from sage.rings.all import Integer
    n = Integer(n)

    if algorithm == 'gap':
        import sage.interfaces.gap as gap
        a = gap.gap.eval('TwoSquares(%s)'%n)
        if a == 'fail':
            raise ValueError, "%s is not a sum of two squares"%n
        x, y = eval(a)
        return Integer(x), Integer(y)
    else:
        raise RuntimeError, "unknown algorithm '%s'"%algorithm
