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

from sage.rings.all import (bernoulli, CyclotomicField,
                            prime_range, QQ, Integer, divisors,
                            LCM, is_squarefree)

def eisenstein_series_qexp(k, prec=10, K=QQ):
    """
    Return the q-expansion of the weight k Eisenstein series
    to precision prec in the field K.

    Here's a rough description of how the algorithm works:
    we know E_k = const + sum_n sigma(n,k-1) * q^n. Now, we
    basically just compute all the sigma(n,k-1) simultaneously,
    as sigma is multiplicative.

    INPUT:
        k -- even positive integer
        prec -- nonnegative integer
        K -- a ring in which B_k/(2*k) is invertible

    EXAMPLES:
        sage: eisenstein_series_qexp(2,5)
        -1/24 + q + 3*q^2 + 4*q^3 + 7*q^4 + O(q^5)
        sage: eisenstein_series_qexp(2,0)
        O(q^0)

    AUTHORS:
        -- William Stein: original implementation
        -- Craig Citro (2007-06-01): rewrote for massive speedup
    """
    k = Integer(k)
    if k%2 or k < 2:
        raise ValueError, "k (=%s) must be an even positive integer"%k
    precision = Integer(prec)
    if precision < 0:
        raise ValueError, "prec (=%s) must an even nonnegative integer"%prec
    R = K[['q']]

    one = Integer(1)
    val = [one] * (prec + 1)

    pow = 0
    ind = 0
    term = 0

    expt = k - one

    for p in prime_range(1,prec+1):

        int_p = int(p)

        ppow = int_p
        mult = p**expt
        term = mult*mult
        last = mult

        while (ppow <= prec):
            ind = ppow
            while (ind <= prec):
                val[ind] = val[ind] * (term - one) // (last - one)
                ind += ppow
            ppow *= int_p
            last = term
            term *= mult

    val[0] = [-bernoulli(k) / (2*k)]

    return R(val, precision)

######################################################################

def __common_minimal_basering(chi, psi):
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
                        chi, psi = __common_minimal_basering(chi, psi)
                        for t in divisors(N/(R*L)):
                            params.append( (chi,psi,t) )
    return params


def __find_eisen_chars_gamma1(N, k):
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



def compute_eisenstein_params(character, k):
    r"""
    Compute and return a list of all parameters $(\chi,\psi,t)$ that
    define the Eisenstein series with given character and weight $k$.

    Only the parity of $k$ is relevant.

    If character is an integer $N$, then the parameters for
    $\Gamma_1(N)$ are computed instead.  Then the condition is that
    $\chi(-1)*\psi(-1) =(-1)^k$.
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


