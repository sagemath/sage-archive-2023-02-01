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

from sage.functions.constants import pi

from sage.rings.all import (bernoulli, CyclotomicField,
                            prime_range, QQ, Integer, divisors,
                            LCM, is_squarefree)

def eisenstein_series_qexp(k, prec=10, K=QQ):
    r"""
    Return the $q$-expansion of the normalized weight $k$ Eisenstein
    series to precision prec in the ring $K$.  (The normalization
    chosen here is the one that forces the coefficient of $q$ to be
    1.)

    Here's a rough description of how the algorithm works: we know
    $E_k = const + \sum_n sigma(n,k-1) q^n$. Now, we basically just
    compute all the $\sigma(n,k-1)$ simultaneously, as $\sigma$ is
    multiplicative.

    INPUT:
        k -- even positive integer
        prec -- nonnegative integer
        K -- a ring in which -(2*k)/B_k is invertible

    EXAMPLES:
        sage: eisenstein_series_qexp(2,5)
        -1/24 + q + 3*q^2 + 4*q^3 + 7*q^4 + O(q^5)
        sage: eisenstein_series_qexp(2,0)
        O(q^0)
        sage: eisenstein_series_qexp(2,5,GF(7))
        2 + q + 3*q^2 + 4*q^3 + O(q^5)

    AUTHORS:
        -- William Stein: original implementation
        -- Craig Citro (2007-06-01): rewrote for massive speedup
    """
    k = Integer(k)
    if k%2 or k < 2:
        raise ValueError, "k (=%s) must be an even positive integer"%k
    prec = int(prec)
    if prec < 0:
        raise ValueError, "prec (=%s) must be an even nonnegative integer"%prec
    if (prec == 0):
        R = K[['q']]
        return R(0).add_bigoh(0)

    one = Integer(1)
    val = [one] * prec
    expt = k - one

    try:
        a0inv = - (2*k) / bernoulli(k)
        a0 = K(1/a0inv)
    except ZeroDivisionError:
        raise ValueError, "-(2*k)/B_k (=%s) must be invertible in the %r"%(a0inv, K)

    for p in prime_range(1,prec):

        int_p = int(p)

        ppow = int_p
        mult = p**expt
        term = mult*mult
        last = mult

        while (ppow < prec):
            ind = ppow
            while (ind < prec):
                val[ind] = (val[ind]*(term-one)).divide_knowing_divisible_by(last - one)
                ind += ppow
            ppow *= int_p
            last = term
            term *= mult

    val[0] = a0
    R = K[['q']]
    return R(val, prec=prec, check=False)

######################################################################

def __common_minimal_basering(chi, psi):
    """
    Find the smallest basering over which chi and psi are valued, and
    return new chi and psi valued in that ring.

    EXAMPLES:
        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(1).0, DirichletGroup(1).0)
        ([1], [1])

        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(3).0, DirichletGroup(5).0)
        ([-1], [zeta4])

        sage: sage.modular.modform.eis_series.__common_minimal_basering(DirichletGroup(12).0, DirichletGroup(36).0)
        ([-1, 1], [-1, 1])
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

    EXAMPLES:
        sage: sage.modular.modform.eis_series.__find_eisen_chars(DirichletGroup(36).0, 4)
        []

        sage: sage.modular.modform.eis_series.__find_eisen_chars(DirichletGroup(36).0, 5)
        [([1, 1], [-1, 1], 1),
        ([1, 1], [-1, 1], 3),
        ([1, 1], [-1, 1], 9),
        ([1, -1], [-1, -1], 1),
        ([-1, 1], [1, 1], 1),
        ([-1, 1], [1, 1], 3),
        ([-1, 1], [1, 1], 9),
        ([-1, -1], [1, -1], 1)]
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
    #  (2) chi == eps*psi, where eps is the nebentypus character of self.
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
                            params.append( (chi0,psi0,t) )
    return params


def __find_eisen_chars_gamma1(N, k):
    """
    Find all characters chi such that (chi,k) gives rise to
    Eisenstein series for Gamma1(N).

    EXAMPLES:
        sage: sage.modular.modform.eis_series.__find_eisen_chars_gamma1(12, 4)
        [([1, 1], [1, 1], 1),
        ([1, 1], [1, 1], 2),
        ([1, 1], [1, 1], 3),
        ([1, 1], [1, 1], 4),
        ([1, 1], [1, 1], 6),
        ([1, 1], [1, 1], 12),
        ([1, 1], [-1, -1], 1),
        ([-1, -1], [1, 1], 1),
        ([-1, 1], [1, -1], 1),
        ([1, -1], [-1, 1], 1)]

        sage: sage.modular.modform.eis_series.__find_eisen_chars_gamma1(12, 5)
        [([1, 1], [-1, 1], 1),
        ([1, 1], [-1, 1], 3),
        ([-1, 1], [1, 1], 1),
        ([-1, 1], [1, 1], 3),
        ([1, 1], [1, -1], 1),
        ([1, 1], [1, -1], 2),
        ([1, 1], [1, -1], 4),
        ([1, -1], [1, 1], 1),
        ([1, -1], [1, 1], 2),
        ([1, -1], [1, 1], 4)]
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
                pairs.append((chi, psi))
                if i!=j: pairs.append((psi,chi))
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

def eisenstein_series_Lseries(weight, prec=53,
               max_imaginary_part=0,
               max_asymp_coeffs=40):
    r"""
    Return the L-series of the weight $2k$ Eisenstein series
    on $\SL_2(\Z)$.

    This actually returns an interface to Tim Dokchitser's program
    for computing with the L-series of the Eisenstein series

    INPUT:
       weight -- even integer
       prec -- integer (bits precision)
       max_imaginary_part -- real number
       max_asymp_coeffs -- integer

    OUTPUT:
       The L-series of the Eisenstein series.

    EXAMPLES:
    We compute with the L-series of $E_{16}$ and then $E_{20}$:
       sage: L = eisenstein_series_Lseries(16)
       sage: L(1)
       -0.291657724743873
       sage: L = eisenstein_series_Lseries(20)
       sage: L(2)
       -5.02355351645987
    """
    f = eisenstein_series_qexp(weight,prec)
    from sage.lfunctions.all import Dokchitser
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
    Compute and return a list of all parameters $(\chi,\psi,t)$ that
    define the Eisenstein series with given character and weight $k$.

    Only the parity of $k$ is relevant.

    If character is an integer $N$, then the parameters for
    $\Gamma_1(N)$ are computed instead.  Then the condition is that
    $\chi(-1)*\psi(-1) =(-1)^k$.

    EXAMPLES:
        sage: sage.modular.modform.eis_series.compute_eisenstein_params(DirichletGroup(30).0, 3)
        []

        sage: sage.modular.modform.eis_series.compute_eisenstein_params(DirichletGroup(30).0, 4)
        [([1, 1, 1], [1, 1, 1], 1),
        ([1, 1, 1], [1, 1, 1], 2),
        ([1, 1, 1], [1, 1, 1], 3),
        ([1, 1, 1], [1, 1, 1], 5),
        ([1, 1, 1], [1, 1, 1], 6),
        ([1, 1, 1], [1, 1, 1], 10),
        ([1, 1, 1], [1, 1, 1], 15),
        ([1, 1, 1], [1, 1, 1], 30)]
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
