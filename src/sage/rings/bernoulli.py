"""
Computation of Bernoulli numbers
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


def bernoulli_python(m):
    r"""
    Returns the Bernoulli number $B_m$ computed using
    a Python implementation.
    """
    import sage.rings.all as rings
    import sage.misc.all as misc

    m = int(m)
    if m == 1:
        return rings.Rational('-1/2')
    if m % 2:
        return 0

    tm = misc.verbose('bernoulli_python setup')
    RR = rings.RealField(4*4*m)
    tm = misc.verbose('computing pi...', tm)
    pi = RR.pi()
    tm = misc.verbose('computing factorial...', tm)
    m_factorial = rings.factorial(m)
    tm = misc.verbose('computing pi pow...', tm)
    K = 2*m_factorial/((2*pi)**m)
    tm = misc.verbose('computing P...', tm)
    P = rings.prime_range(m+2)

    # IDIOTIC -- should compute by factoring m!
    d = misc.prod([p for p in P if m % (p-1) == 0])

    tm = misc.verbose('computing N...', tm)
    #N = long((K*d)**(1.0/(m-1.0)) + 1)
    dK = d*K
    n = str(dK).find('.') + 1
    N = 10**(float(n)/float(m-1))
    tm = misc.verbose('N = %s; computing product...'%N, tm)
    assert N < m
    z = 1
    z_prev = z
    Rl = rings.RealField()
    for p in P:
        if p > N:
            break
        z /= 1 - 1/(RR(p)**RR(m))
        diff = abs(z - z_prev).log()
        print p, Rl(diff), -16*m
        if diff < -16*m:
            break
        z_prev =z
    a = long((dK*z).ceil())
    if m % 4 == 0:
        a = -a
    return rings.Rational(a)/rings.Rational(d)
    #return K



def bernoulli_cf(m):
    r"""
    Returns the Bernoulli number $B_m$.
    """
    import sage.rings.all as rings
    import sage.misc.all as misc

    m = int(m)
    if m == 1:
        return rings.Rational('-1/2')
    if m % 2:
        return 0

    tm = misc.verbose('bernoulli_python setup')
    RR = rings.RealField(4*4*m)
    tm = misc.verbose('computing pi...', tm)
    pi = RR.pi()
    tm = misc.verbose('computing factorial...', tm)
    m_factorial = rings.factorial(m)
    tm = misc.verbose('computing pi pow...', tm)
    pi_pow = (2*pi)**m
    tm = misc.verbose('computing quo...', tm)
    K = 2*m_factorial/pi_pow
    tm = misc.verbose('computing zeta times K...', tm)
    s = RR(m).zeta() * K
    tm = misc.verbose('done', tm)
    return s
