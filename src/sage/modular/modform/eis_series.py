# -*- coding: utf-8 -*-
"""
Eisenstein Series
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import sage.misc.all as misc

import sage.modular.dirichlet as dirichlet

from sage.rings.all import ComplexField, RealField, Integer

from sage.rings.all import (bernoulli, CyclotomicField, prime_range,
                            is_FiniteField, ZZ, QQ, Integer, divisors,
                            LCM, is_squarefree)
from sage.rings.power_series_ring import PowerSeriesRing
from eis_series_cython import eisenstein_series_poly

def eisenstein_series_qexp(k, prec = 10, K=QQ, var='q') :
    r"""
    Return the `q`-expansion of the normalized weight `k` Eisenstein series on
    `{\rm SL}_2(\ZZ)` to precision prec in the ring `K`.  (The normalization
    chosen here is the one that forces the coefficient of `q` to be 1.)

    INPUT:

    - ``k`` - an even positive integer

    - ``prec`` - (default: 10) a nonnegative integer

    - ``K`` - (default: `\QQ`) a ring in which the denominator of `B_k / 2k` is invertible

    - ``var`` - (default: 'q') variable name to use for q-expansion

    ALGORITHM:

        We know `E_k = \text{constant} + \sum_n \sigma_{k-1}(n) q^n`. So we
        compute all the `\sigma_{k-1}(n)` simultaneously, using the fact that
        `\sigma` is multiplicative.

    EXAMPLES::

        sage: eisenstein_series_qexp(2,5)
        -1/24 + q + 3*q^2 + 4*q^3 + 7*q^4 + O(q^5)
        sage: eisenstein_series_qexp(2,0)
        O(q^0)
        sage: eisenstein_series_qexp(2,5,GF(7))
        2 + q + 3*q^2 + 4*q^3 + O(q^5)
        sage: eisenstein_series_qexp(2,5,GF(7),var='T')
        2 + T + 3*T^2 + 4*T^3 + O(T^5)

        sage: eisenstein_series_qexp(10, 30, GF(17))
        15 + q + 3*q^2 + 15*q^3 + 7*q^4 + 13*q^5 + 11*q^6 + 11*q^7 + 15*q^8 + 7*q^9 + 5*q^10 + 7*q^11 + 3*q^12 + 14*q^13 + 16*q^14 + 8*q^15 + 14*q^16 + q^17 + 4*q^18 + 3*q^19 + 6*q^20 + 12*q^21 + 4*q^22 + 12*q^23 + 4*q^24 + 4*q^25 + 8*q^26 + 14*q^27 + 9*q^28 + 6*q^29 + O(q^30)

    TESTS:

    This shows that the bug reported at trac 8291 is fixed::

        sage: eisenstein_series_qexp(26, 10, GF(13))
        7 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + 12*q^6 + 8*q^7 + 2*q^8 + O(q^10)

    AUTHORS:

    - William Stein: original implementation

    - Craig Citro (2007-06-01): rewrote for massive speedup

    - Martin Raum (2009-08-02): port to cython for speedup

    - David Loeffler (2010-04-07): work around an integer overflow when k is large
    """
    ## we use this to prevent computation if it would fail anyway.
    if k <= 0 or k % 2 == 1 :
        raise ValueError, "k must be positive and even"

    a0 = - bernoulli(k) / (2*k)
    a0den = a0.denominator()
    try:
        a0fac = K(1/a0den)
    except ZeroDivisionError:
        raise ValueError, "The denominator of -B_k/(2*k) (=%s) must be invertible in the ring %s"%(a0den, K)


    R = PowerSeriesRing(K, var)
    if K is QQ :
        return a0fac*R(eisenstein_series_poly(k, prec).list(), prec=prec, check=False)
    else:
        # this is a temporary fix due to a change in the
        # polynomial constructor over finite fields; this
        # is a notable speed regression, to be fixed soon.
        return a0fac*R(eisenstein_series_poly(k, prec).list(), prec=prec, check=True)

######################################################################

def __common_minimal_basering(chi, psi):
    """
    Find the smallest basering over which chi and psi are valued, and
    return new chi and psi valued in that ring.

    EXAMPLES::

        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(1).0, DirichletGroup(1).0)
        (Dirichlet character modulo 1 of conductor 1 mapping 0 |--> 1, Dirichlet character modulo 1 of conductor 1 mapping 0 |--> 1)

        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(3).0, DirichletGroup(5).0)
        (Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4)

        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(12).0, DirichletGroup(36).0)
        (Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1, Dirichlet character modulo 36 of conductor 4 mapping 19 |--> -1, 29 |--> 1)
    """
    chi = chi.minimize_base_ring()
    psi = psi.minimize_base_ring()
    n = LCM(chi.base_ring().zeta().multiplicative_order(),\
                  psi.base_ring().zeta().multiplicative_order())
    if n <= 2:
        K = QQ
    else:
        K = CyclotomicField(n)
    chi = chi.change_ring(K)
    psi = psi.change_ring(K)
    return chi, psi

#def prim(eps):
#    print "making eps with modulus %s primitive"%eps.modulus()
#    return eps.primitive_character()

def __find_eisen_chars(character, k):
    """
    Find all characters chi such that (chi,k) gives rise to
    Eisenstein series for the group determined by character.

    EXAMPLES::

        sage: sage.modular.modform.eis_series.__find_eisen_chars(DirichletGroup(36).0, 4)
        []

        sage: pars =  sage.modular.modform.eis_series.__find_eisen_chars(DirichletGroup(36).0, 5)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (-1, 1), 1),
        ((1, 1), (-1, 1), 3),
        ((1, 1), (-1, 1), 9),
        ((1, -1), (-1, -1), 1),
        ((-1, 1), (1, 1), 1),
        ((-1, 1), (1, 1), 3),
        ((-1, 1), (1, 1), 9),
        ((-1, -1), (1, -1), 1)]
    """
    N = character.modulus()
    if character.is_trivial():
        if k%2 != 0:
            return []
        char_inv = ~character
        V = [(character, char_inv, t) for t in divisors(N) if t>1]
        if k != 2:
            V.insert(0,(character, char_inv, 1))
        if is_squarefree(N):
            return V
        # Now include all pairs (chi,chi^(-1)) such that cond(chi)^2 divides N:
        # TODO: Optimize -- this is presumably way too hard work below.
        G = dirichlet.DirichletGroup(N)
        for chi in G:
            if not chi.is_trivial():
                f = chi.conductor()
                if N % (f**2) == 0:
                    chi = chi.minimize_base_ring()
                    chi_inv = ~chi
                    for t in divisors(N//(f**2)):
                        V.insert(0, (chi, chi_inv, t))
        return V


    eps = character
    if eps(-1) != (-1)**k:
        return []
    eps = eps.maximize_base_ring()
    G = eps.parent()

    # Find all pairs chi, psi such that:
    #
    #  (1) cond(chi)*cond(psi) divides the level, and
    #
    #  (2) chi*psi == eps, where eps is the nebentypus character of self.
    #
    # See [Miyake, Modular Forms] Lemma 7.1.1.

    K = G.base_ring()
    C = {}

    t0 = misc.cputime()

    for e in G:
        m = Integer(e.conductor())
        if C.has_key(m):
            C[m].append(e)
        else:
            C[m] = [e]

    misc.verbose("Enumeration with conductors.",t0)

    params = []
    for L in divisors(N):
        misc.verbose("divisor %s"%L)
        if not C.has_key(L):
            continue
        GL = C[L]
        for R in divisors(N/L):
            if not C.has_key(R):
                continue
            GR = C[R]
            for chi in GL:
                for psi in GR:
                    if chi*psi == eps:
                        chi0, psi0 = __common_minimal_basering(chi, psi)
                        for t in divisors(N//(R*L)):
                            if k != 1 or ((psi0, chi0, t) not in params):
                                params.append( (chi0,psi0,t) )
    return params


def __find_eisen_chars_gamma1(N, k):
    """
    Find all characters chi such that (chi,k) gives rise to
    Eisenstein series for Gamma1(N).

    EXAMPLES::

        sage: pars = sage.modular.modform.eis_series.__find_eisen_chars_gamma1(12, 4)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (1, 1), 1),
        ((1, 1), (1, 1), 2),
        ((1, 1), (1, 1), 3),
        ((1, 1), (1, 1), 4),
        ((1, 1), (1, 1), 6),
        ((1, 1), (1, 1), 12),
        ((1, 1), (-1, -1), 1),
        ((-1, -1), (1, 1), 1),
        ((-1, 1), (1, -1), 1),
        ((1, -1), (-1, 1), 1)]

        sage: pars =  sage.modular.modform.eis_series.__find_eisen_chars_gamma1(12, 5)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (-1, 1), 1),
        ((1, 1), (-1, 1), 3),
        ((-1, 1), (1, 1), 1),
        ((-1, 1), (1, 1), 3),
        ((1, 1), (1, -1), 1),
        ((1, 1), (1, -1), 2),
        ((1, 1), (1, -1), 4),
        ((1, -1), (1, 1), 1),
        ((1, -1), (1, 1), 2),
        ((1, -1), (1, 1), 4)]
    """
    pairs = []
    s = (-1)**k
    G = dirichlet.DirichletGroup(N)
    E = list(G)
    parity = [c(-1) for c in E]
    for i in range(len(E)):
        for j in range(i,len(E)):
            if parity[i]*parity[j] == s and N % (E[i].conductor()*E[j].conductor()) == 0:
                chi, psi = __common_minimal_basering(E[i], E[j])
                if k != 1:
                    pairs.append((chi, psi))
                    if i!=j: pairs.append((psi,chi))
                else:
                    # if weight is 1 then (chi, psi) and (chi, psi) are the
                    # same form
                    if psi.is_trivial() and not chi.is_trivial():
                        # need to put the trivial character first to get the L-value right
                        pairs.append((psi, chi))
                    else:
                        pairs.append((chi, psi))
        #end fors
    #end if

    triples = []
    D = divisors(N)
    for chi, psi in pairs:
        c_chi = chi.conductor()
        c_psi = psi.conductor()
        D = divisors(N/(c_chi * c_psi))
        if (k==2 and chi.is_trivial() and psi.is_trivial()):
            D.remove(1)
        chi, psi = __common_minimal_basering(chi, psi)
        for t in D:
            triples.append((chi, psi, t))
    return triples

def eisenstein_series_lseries(weight, prec=53,
               max_imaginary_part=0,
               max_asymp_coeffs=40):
    r"""
    Return the L-series of the weight `2k` Eisenstein series
    on `\mathrm{SL}_2(\ZZ)`.

    This actually returns an interface to Tim Dokchitser's program
    for computing with the L-series of the Eisenstein series

    INPUT:

    - ``weight`` - even integer

    - ``prec`` - integer (bits precision)

    - ``max_imaginary_part`` - real number

    - ``max_asymp_coeffs`` - integer

    OUTPUT:

    The L-series of the Eisenstein series.

    EXAMPLES:

    We compute with the L-series of `E_{16}` and then `E_{20}`::

       sage: L = eisenstein_series_lseries(16)
       sage: L(1)
       -0.291657724743873
       sage: L = eisenstein_series_lseries(20)
       sage: L(2)
       -5.02355351645987
    """
    f = eisenstein_series_qexp(weight,prec)
    from sage.lfunctions.all import Dokchitser
    from sage.symbolic.constants import pi
    key = (prec, max_imaginary_part, max_asymp_coeffs)
    j = weight
    L = Dokchitser(conductor = 1,
                   gammaV = [0,1],
                   weight = j,
                   eps = (-1)**Integer((j/2)),
                   poles = [j],
                   residues = [(-1)**Integer((j/2))*(float(pi))**(0.5)*bernoulli(j)/j],
                   prec = prec)

    s = 'coeff = %s;'%f.list()
    L.init_coeffs('coeff[k+1]',pari_precode = s,
                  max_imaginary_part=max_imaginary_part,
                  max_asymp_coeffs=max_asymp_coeffs)
    L.check_functional_equation()
    L.rename('L-series associated to the weight %s Eisenstein series %s on SL_2(Z)'%(j,f))
    return L

def compute_eisenstein_params(character, k):
    r"""
    Compute and return a list of all parameters `(\chi,\psi,t)` that
    define the Eisenstein series with given character and weight `k`.

    Only the parity of `k` is relevant (unless k = 1, which is a slightly different case).

    If character is an integer `N`, then the parameters for
    `\Gamma_1(N)` are computed instead.  Then the condition is that
    `\chi(-1)*\psi(-1) =(-1)^k`.

    EXAMPLES::

        sage: sage.modular.modform.eis_series.compute_eisenstein_params(DirichletGroup(30)(1), 3)
        []

        sage: pars =  sage.modular.modform.eis_series.compute_eisenstein_params(DirichletGroup(30)(1), 4)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (1, 1), 1),
        ((1, 1), (1, 1), 2),
        ((1, 1), (1, 1), 3),
        ((1, 1), (1, 1), 5),
        ((1, 1), (1, 1), 6),
        ((1, 1), (1, 1), 10),
        ((1, 1), (1, 1), 15),
        ((1, 1), (1, 1), 30)]

        sage: pars = sage.modular.modform.eis_series.compute_eisenstein_params(15, 1)
        sage: [(x[0].values_on_gens(), x[1].values_on_gens(), x[2]) for x in pars]
        [((1, 1), (-1, 1), 1),
        ((1, 1), (-1, 1), 5),
        ((1, 1), (1, zeta4), 1),
        ((1, 1), (1, zeta4), 3),
        ((1, 1), (-1, -1), 1),
        ((1, 1), (1, -zeta4), 1),
        ((1, 1), (1, -zeta4), 3),
        ((-1, 1), (1, -1), 1)]

        sage: sage.modular.modform.eis_series.compute_eisenstein_params(DirichletGroup(15).0, 1)
        [(Dirichlet character modulo 15 of conductor 1 mapping 11 |--> 1, 7 |--> 1, Dirichlet character modulo 15 of conductor 3 mapping 11 |--> -1, 7 |--> 1, 1),
        (Dirichlet character modulo 15 of conductor 1 mapping 11 |--> 1, 7 |--> 1, Dirichlet character modulo 15 of conductor 3 mapping 11 |--> -1, 7 |--> 1, 5)]
    """
    if isinstance(character, (int,long,Integer)):
        N = character
        character = None
    else:
        N = character.modulus()

    if character != None:
        return __find_eisen_chars(character, k)
    else:
        return __find_eisen_chars_gamma1(N, k)
