r"""
Dimensions of spaces of modular forms

AUTHORS:
    -- William Stein
    -- Jordi Quer

ACKNOWLEDGEMENT:
    The dimension formulas and implementations in this module grew out
    of a program that Bruce Kaskel wrote (around 1996) in PARI, which
    Kevin Buzzard subsequently extended.  I (William Stein) then
    implemented it in C++ for Hecke.  I also implemented it in Magma.
    Also, the functions for dimensions of spaces with nontrivial
    character are based on a paper (that has no proofs) by Cohen and
    Oesterle (Springer Lecture notes in math, volume 627, pages
    69--78).  The formulas for GammaH(N) were found and implemented by
    Jordi Quer.

The formulas here are more complete than in Hecke or Magma.
"""

##########################################################################
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##########################################################################

import math

from sage.rings.arith import (factor, euler_phi as phi, divisors, is_prime,
                              valuation, kronecker_symbol, gcd, euler_phi, lcm)
import sage.modular.congroup as congroup
from sage.misc.misc import mul
from sage.rings.all import Mod, Integer, IntegerRing, IntegerModRing, ZZ
from sage.rings.rational_field import frac
import dirichlet
Z = ZZ  # useful abbreviation.

def mu0(n):
    r"""
    Return value of the combinatorial function
    $$
       \mu_0(n) = \prod_{p^r|n} (p+1)p^{r-1},
    $$
    where the product is over the maximal prime power
    divisors of $n$.

    INPUT:
        n -- an integer

    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.mu0(n) for n in [1..19]]
        [1, 3, 4, 6, 6, 12, 8, 12, 12, 18, 12, 24, 14, 24, 24, 24, 18, 36, 20]
    """
    return mul([(p+1)*(p**(r-1)) for p, r in factor(n)])

def mu20(n):
    r"""
    Return value of the arithmetic function $\mu_{2,0}(n)$.

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.mu20(n) for n in [1..19]]
        [1, 1, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 2, 0, 0]
    """
    if n%4 == 0:
        return Integer(0)
    return mul([1 + kronecker_symbol(-4,p) for p, _ in factor(n)])

def mu30(n):
    """
    Return value of the arithmetic function $\mu_{3,0}(n)$.

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.mu30(n) for n in [1..19]]
        [1, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2]
    """
    if n%2==0 or n%9==0:
        return Integer(0)
    return mul([1 + kronecker_symbol(-3,p) for p, _ in factor(n)])

def c0(n):
    """
    Return the number of cusps of the modular curve $X_0(n)$.

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.c0(n) for n in [1..19]]
        [1, 2, 2, 3, 2, 4, 2, 4, 4, 4, 2, 6, 2, 4, 4, 6, 2, 8, 2]
        sage: [sage.modular.dims.c0(n) for n in prime_range(2,100)]
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    """
    return sum([phi(gcd(d,n//d)) for d in divisors(n)])

def g0(n):
    """
    Return the genus of the modular curve $X_0(n)$.


    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.g0(n) for n in [1..23]]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 2, 2]
        sage: [n for n in [1..200] if sage.modular.dims.g0(n) == 1]
        [11, 14, 15, 17, 19, 20, 21, 24, 27, 32, 36, 49]
    """
    return Integer(1 + frac(mu0(n),12) - frac(mu20(n),4) - \
               frac(mu30(n),3) - frac(c0(n),2))

def mu1(n):
    r"""
    Return value of the combinatorial function $\mu_1(n)$.

    INPUT:
        n -- an integer
    OUTPUT:
        Integer

    EXAMPLES:
        sage: [sage.modular.dims.mu1(n) for n in [1..16]]
        [1, 3, 4, 6, 12, 12, 24, 24, 36, 36, 60, 48, 84, 72, 96, 96]
    """
    if n <= 2:
        return mu0(n)
    return (phi(n)*mu0(n))/2

def mu21(n):
    r"""
    Return value of $\mu_{2,1}(n)$.  This equals $\mu_{2,0}(n)$ if
    $n<4$ and equals 0 otherwise.

    EXAMPLES:
    """
    if n<4:
        return mu20(n)
    return Integer(0)

def mu31(n):
    r"""
    Return value of $\mu_{3,1}(n)$.  This equals $\mu_{3,0}(n)$ if
    $n<4$ and equals 0 otherwise.

    EXAMPLES:
        sage: [sage.modular.dims.mu31(n) for n in [1..10]]
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    """
    if n<4:
        return mu30(n)
    return Integer(0)

def c1(n):
    """
    Return value of $c_1(n)$, which is the number of cusps on the
    modular curve $X_1(n)$.

    EXAMPLES:
        sage: [sage.modular.dims.c1(n) for n in prime_range(2,100)]
        [2, 2, 4, 6, 10, 12, 16, 18, 22, 28, 30, 36, 40, 42, 46, 52, 58, 60, 66, 70, 72, 78, 82, 88, 96]
    """
    if n<3:
        return c0(n)
    if n==4:
        return Integer(3)
    return Integer(sum([frac(phi(d)*phi(n/d),2) \
                    for d in divisors(n)]))

def g1(n):
    """
    Return the genus of the modular curve $X_1(n)$.

    EXAMPLES:
        sage: [sage.modular.dims.g1(n) for n in prime_range(2,100)]
        [0, 0, 0, 0, 1, 2, 5, 7, 12, 22, 26, 40, 51, 57, 70, 92, 117, 126, 155, 176, 187, 222, 247, 287, 345]
    """
    return Integer(1+frac(mu1(n),12)-frac(mu21(n),4)-frac(mu31(n),3) - frac(c1(n),2))

def ss0(n,p):
    """

    INPUT:
        n -- integer
        p -- a prime

    OUTPUT:
        Integer

    """
    assert is_prime(p), "p must be prime"
    assert n%p==0, "n must be divisible by p"
    return g0(n*p) - 2*g0(n) + 1

def muXNp(n,p):
    return mu1(n)*mu0(p)

def mu2XNp(n,p):
    return 0

def mu3XNp(n,p):
    return 0

def cXNp(n):
    return 2*c1(n)

def gXNp(n,p):
    if n<4:
        return g0(n*p)
    return int(1+frac(muXNp(n,p),12)-frac(mu2XNp(n,p),4) \
               - frac(mu3XNp(n,p),3) - frac(cXNp(n),2))

def ss1(n,p):
    assert is_prime(p) and (n%p != 0)
    return gXNp(n,p) - 2*g1(n) + 1

def eisen(p):
   assert is_prime(p)
   return frac(p-1,12).numerator()

def S0(n,k):
    n = int(n); k = int(k)
    assert n>0
    if k<=0 or k%2!=0:
        return 0
    if k==2:
        return g0(n)
    return int((k-1)*(g0(n)-1) + \
               (frac(k,2)-1)*c0(n)+mu20(n)*(k//4)+mu30(n)*(k//3))  # // is floor div

def S1(n,k):
    n = int(n); k = int(k)
    assert n>0
    if k<=0 or (n<=2 and k%2!=0):
        return 0
    assert k!=1, "weight 1 dimension not programmed"
    if k==2:
        return g1(n)
    if n<=2:
        return S0(n,k)
    a = (k-1)*(g1(n)-1)+(frac(k,2)-1)*c1(n)
    if n == 4 and k%2!=0:
        a += frac(1,2)
    elif n == 3:
        a += k//3  # // is floor div
    return int(a)

def idxG0(N):
    r"""
    The index $[\Gamma_0(N):\SL_2(\Z)]$.
    """
    return mul([(p+1)*p**(e-1) for p, e in factor(N)])

def idxG1(N):
    r"""
    The index $[\Gamma_1(N):\SL_2(\Z)]$.
    """
    return phi(N)*idxG0(N)

#    Formula of Cohen-Oesterle for dim S_k(Gamma_1(N),eps).
#    REF: Springer Lecture notes in math, volume 627, pages 69--78.
#    The functions CO_delta and CO_nu, which were first written by Kevin Buzzard,
#    are used only by the function CohenOesterle.

def CO_delta(r,p,N,eps):
    assert is_prime(p)
    K = eps.base_ring()
    if p%4 == 3:
        return K(0)
    if p==2:
        if r==1:
            return K(1)
        return K(0)
    # interesting case: p=1(mod 4).
    # omega is a primitive 4th root of unity mod p.
    omega = (IntegerModRing(p).unit_gens()[0])**((p-1)//4)
    # this n is within a p-power root of a "local" 4th root of 1 modulo p.
    n = Mod(int(omega.crt(Mod(1,N//(p**r)))),N)
    n = n**(p**(r-1))   # this is correct now
    t = eps(n)
    if t==K(1):
        return K(2)
    if t==K(-1):
        return K(-2)
    return K(0)

def CO_nu(r, p, N, eps):
    K = eps.base_ring()
    if p%3==2:
        return K(0)
    if p==3:
        if r==1:
            return K(1)
        return K(0)
    # interesting case: p=1(mod 3)
    # omega is a cube root of 1 mod p.
    omega = (IntegerModRing(p).unit_gens()[0])**((p-1)//3)
    n = Mod(omega.crt(Mod(1,N//(p**r))), N)  # within a p-power root of a "local" cube root of 1 mod p.
    n = n**(p**(r-1))  # this is right now
    t = eps(n)
    if t==K(1):
        return K(2)
    return K(-1)

#todo: I had the following comment in my magma code.  check it.
# Kevin's clever function has a bug, so I'm not using it now:
#  K := CyclotomicField(3);
#  eps := DirichletGroup(7*43,K).1^2;
#  CuspidalSubspace(ModularForms([eps],2));
#  boom!

def CohenOesterle(eps, k):
    N    = eps.modulus()
    facN = factor(N)
    f    = eps.conductor()
    gamma_k = 0
    if k%4==2:
        gamma_k = frac(-1,4)
    elif k%4==0:
        gamma_k = frac(1,4)
    mu_k = 0
    if k%3==2:
        mu_k = frac(-1,3)
    elif k%3==0:
        mu_k = frac(1,3)
    def _lambda(r,s,p):
        if 2*s<=r:
            if r%2==0:
                return p**(r//2) + p**((r//2)-1)
            return 2*p**((r-1)//2)
        return 2*(p**(r-s))
    #end def of lambda
    K = eps.base_ring()
    return K(frac(-1,2) * mul([_lambda(r,valuation(f,p),p) for p, r in facN]) + \
               gamma_k * mul([CO_delta(r,p,N,eps)         for p, r in facN]) + \
                mu_k    * mul([CO_nu(r,p,N,eps)            for p, r in facN]))

def dimension_cusp_forms_eps(eps, k=2):
    """
    The dimension of the space of cusp forms of weight k and character
    eps.

    INPUT:
        eps -- a Dirichlet character
        k -- integer, a weight >= 2.

    OUTPUT:
        integer -- the dimension
    """
    if isinstance(eps, (int,long) ):
        return dimension_cusp_forms_gamma0(eps,k)

    if k < 0:
        return Z(0)
    if eps.is_even():
        if k % 2 == 1:
            return Z(0)
    else:  # odd
        if k % 2 == 0:
            return Z(0)
    if k == 0:
        return Z(0)
    elif k == 1:
        raise NotImplementedError, "Computation of dimensions of spaces of weight 1 modular forms not implemented."

    N = eps.modulus()
    if eps.is_trivial():
        return Z(S0(N,k))
    if (eps.is_odd() and k%2==0) or (eps.is_even() and k%2==1):
        return Z(0)
    K = eps.base_ring()
    return Z ( K(idxG0(N)*frac(k-1,12)) + CohenOesterle(eps,k) )


def dimension_eis_eps(eps, k=2):
    r"""
    The dimension of the space of eisenstein series of weight $k$ and
    character $\varepsilon$.

    INPUT:
        eps -- a Dirichlet character
        k -- integer, a weight >= 2.

    OUTPUT:
        integer -- the dimension
    """
    if isinstance(eps, (int,long) ):
        return dimension_eis_gamma0(eps,k)

    if k < 0:
        return Z(0)
    if eps.is_even():
        if k % 2 == 1:
            return Z(0)
    else:  # odd
        if k % 2 == 0:
            return Z(0)
    if k == 0:
        if eps.is_trivial():
            return Z(1)    # the constants
        else:
            return Z(0)

    N = eps.modulus()
    if (eps.is_odd() and k%2==0) or (eps.is_even() and k%2==1):
        return Z(0)
    if eps.is_trivial():
        return dimension_eis(congroup.Gamma0(N),k)
    K = eps.base_ring()
    j = 2-k
    # We use the Cohen-Oesterle formula in a subtle way to
    # compute dim M_k(N,eps) (see Ch. 6 of my book on
    # computing with modular forms).
    alpha = -Z ( K(idxG0(N)*frac(j-1,12)) + CohenOesterle(eps,j) )
    if k == 1:
        return alpha
    else:
        return alpha - dimension_cusp_forms_eps(eps, k)


def dimension_cusp_forms_gamma0(N,k=2):
    r"""
    The dimension of the space $S_k(\Gamma_0(N))$ of cusp forms.

    INPUT:
        N -- integer
        k -- integer, weight >= 2

    OUTPUT:
        integer -- the dimension
    """
    if N <= 0:
        raise ArithmeticError, "the level N (=%s) must be positive"%N
    if k < 0:
        return Z(0)
    elif k == 0:
        return Z(0)
    elif k%2 == 1:
        return Z(0)
    return Z(S0(N,k))

def dimension_cusp_forms_gamma1(N,k=2):
    r"""
    The dimension of the space $S_k(\Gamma_1(N))$ of cusp forms.

    INPUT:
        N -- integer
        k -- integer, weight >= 2

    OUTPUT:
        integer -- the dimension

    EXAMPLES:
    """
    if N <= 0:
        raise ArithmeticError, "the level N (=%s) must be positive"%N
    if k < 0:
        return Z(0)
    elif k == 0:
        return Z(0)
    if k == 1:
        raise NotImplementedError, "computation of dimensions of spaces of weight 1 modular forms not implemented in general."
    return Z(S1(N,k))


def mumu(N):
    assert N>=1
    p = 1
    for _,r in factor(N):
        if r > 2:
            return Z(0)
        elif r == 1:
            p *= -2
    return Z(p)

def dimension_new_cusp_forms_gamma0(N, k=2, p=0):
    r"""
    Dimension of the p-new subspace of $S_k(\Gamma_0(N))$.
    If $p=0$, dimension of the new subspace.

    INPUT:
        N -- a positive integer
        k -- an integer
        p -- a prime number

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.dimension_new_cusp_forms_gamma0(100, 2, 5)
        5

    Independently compute the dimension 5 above:
        sage: m = ModularSymbols(100, 2,sign=1).cuspidal_subspace()
        sage: m.new_subspace(5)
        Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 18 for Gamma_0(100) of weight 2 with sign 1 over Rational Field
    """
    if N <= 0:
        raise ArithmeticError, "the level N (=%s) must be positive"%N
    if k < 0:
        return Z(0)
    elif k == 0:
        return Z(0)
    elif k%2 == 1:
        return Z(0)
    if p==0 or N%p!=0:
        return sum([dimension_cusp_forms_gamma0(M,k)*mumu(N//M) \
                    for M in divisors(N)])
    return dimension_cusp_forms_gamma0(N,k) - \
           2*dimension_new_cusp_forms_gamma0(N//p,k)

def dimension_new_cusp_forms_gamma1(N,k=2,p=0):
    r"""
    Return the dimension of the $p$-new subspace of
    $S_k(\Gamma_1(N))$.  If $p=0$, return the dimension of the new
    subspace.

    INPUT:
        N -- a positive integer
        k -- an integer
        p -- a prime number

    OUTPUT:
        Integer

    EXAMPLES:
        sage: sage.modular.dims.dimension_new_cusp_forms_gamma1(4*25, 2, 5)
        225
    """
    if N <= 0:
        raise ArithmeticError, "the level N (=%s) must be positive"%N
    if k < 0:
        return Z(0)
    elif k == 0:
        return Z(0)
    elif k == 1:
        raise NotImplementedError, "Computation of dimensions of spaces of weight 1 modular forms not implemented."

    if p==0 or N%p!=0:
        return sum([dimension_cusp_forms_gamma1(M,k)*mumu(N//M) \
                    for M in divisors(N)])
    return dimension_cusp_forms_gamma1(N,k) - \
           2*dimension_new_cusp_forms_gamma1(N//p,k)

def dimension_new_cusp_forms_group(group, k=2, p=0):
    """
    Return the dimension of the new space of cusp forms for the
    congruence subgroup group.
    """
    assert isinstance(group, congroup.CongruenceSubgroup), \
           "Argument 1 must be a congruence subgroup."
    if isinstance(group, congroup.Gamma0):
        return dimension_new_cusp_forms_gamma0(group.level(), k, p)
    elif isinstance(group, congroup.Gamma1):
        return dimension_new_cusp_forms_gamma1(group.level(), k, p)
    elif congroup.is_GammaH(group):
        return dimension_new_cusp_forms_H(group, k, p)
    else:
        raise NotImplementedError, "Computing of dimensions for congruence subgroups besides \
        Gamma0 and Gamma1 is not yet implemented."


def dimension_new_cusp_forms_eps(eps, k=2, p=0):
    """
    Dimension of the new subspace (or p-new subspace) of cusp forms of
    weight k and character eps.
    """
    if not isinstance(eps, dirichlet.DirichletCharacter):
        raise TypeError, "eps = (%s) must be a DirichletCharacter"%eps
    if k < 0:
        return Z(0)
    if eps.is_even():
        if k % 2 == 1:
            return Z(0)
    else:  # odd
        if k % 2 == 0:
            return Z(0)
    if k == 0:
        return Z(0)

    elif k == 1:
        raise NotImplementedError, "Computation of dimensions of spaces of weight 1 modular forms not implemented."


    N = eps.modulus()
    if p == 0 or N%p != 0 or valuation(eps.conductor(),p) == valuation(N,p):
        D = [eps.conductor()*d for d in divisors(N//eps.conductor())]
        return sum([dimension_cusp_forms_eps(eps.restrict(M), k)*mumu(N//M) for M in D])
    eps_p = eps.restrict(N//p)
    old = dimension_cusp_forms(eps_p, k)
    return dimension_new_cusp_forms(eps, k) - 2*old




######################################################################
# Computing dimensions of modular forms spaces for Gamma_H.
# Algorithms found and implemented by Jordi Quer.
######################################################################
# degree of the covering $X_H(N)->X$
def muH(N,H):
    lenHpm = len(H)
    if N-Integer(1) not in H: lenHpm*=Integer(2)
    return mul([(p**Integer(2)-Integer(1))*(p**(Integer(2)*r-Integer(2))) for p, r in factor(N)])//lenHpm

# number of elliptic points of order 2 for the group $\Gamma_H(N)$
def nu2H(N,H):
    if N%Integer(4) == Integer(0): return Integer(0)
    for p, r in factor(N):
        if p%Integer(4) ==Integer(3): return Integer(0)
    return (euler_phi(N)//len(H))*len([x for x in H if (x**Integer(2)+Integer(1))%N == Integer(0)])

# number of elliptic points of order 3 for the group $\Gamma_H(N)$
def nu3H(N,H):
    if N%Integer(9) == Integer(0): return Integer(0)
    for p, r in factor(N):
        if p%Integer(3) == Integer(2): return Integer(0)
    lenHpm = len(H)
    if N-Integer(1) not in H: lenHpm*=Integer(2)
    return (euler_phi(N)//lenHpm)*len([x for x in H if (x**Integer(2)+x+Integer(1))%N == Integer(0)])

# number of cusps for the group $\Gamma_H(N)$
def nuinfH(N,H):
    c = Integer(0)
    for d in [d for d in divisors(N) if d**Integer(2)<=N]:
        Nd = lcm(d,N//d)
        Hd = set([x%Nd for x in H])
        lenHd = len(Hd)
        if Nd-Integer(1) not in Hd: lenHd*=Integer(2)
        sumand = euler_phi(d)*euler_phi(N//d)//lenHd
        if d**Integer(2)==N:
            c = c + sumand
        else:
            c = c + Integer(2)*sumand
    return c

# number of regular cusps for the group $\Gamma_H(N)$
def nuinfHreg(N,H):
    c = Integer(0)
    for d in [d for d in divisors(N) if d**Integer(2)<=N]:
        Nd = lcm(d,N//d)
        Hd = set([x%Nd for x in H])
        if Nd-Integer(1) not in Hd:
            sumand = euler_phi(d)*euler_phi(N//d)//(Integer(2)*len(Hd))
            if d**Integer(2)==N:
                c = c + sumand
            else:
                c = c + Integer(2)*sumand
    return c


# genus of the curve $X_H(N)$
# (formulas are multiplied by 12 so that everything is an integer)
def gH(N,H):
    return (Integer(12) + muH(N,H) - Integer(3)*nu2H(N,H) - Integer(4)*nu3H(N,H) - Integer(6)*nuinfH(N,H))//Integer(12)

def genus_H(G):
    return gH(G.level(),G._list_of_elements_in_H())

# coefficients of the numbers of elliptic points in the formulas
def lambda4(k):
    if k%Integer(2) == Integer(1):
        return Integer(0)
    elif k%Integer(4) == Integer(2):
        return -Integer(3)
    else:
        return Integer(3)

def lambda3(k):
    if k%Integer(3) == Integer(1):
        return Integer(0)
    elif k%Integer(3) == Integer(2):
        return -Integer(4)
    else:
        return Integer(4)

# dimensions of spaces of cusp and modular forms
def dimension_cusp_forms_H(G,k):
    N = G.level()
    H = G._list_of_elements_in_H()
    if k%Integer(2) == Integer(1) and N-Integer(1) in H:
        return Integer(0)
    dim = Integer(0)
    if k == Integer(2):
        dim = Integer(12)
    dim += (k-Integer(1))*muH(N,H)
    if k%Integer(2) == Integer(0):
        dim += lambda4(k)*nu2H(N,H)
    dim += lambda3(k)*nu3H(N,H)
    if k%Integer(2) == Integer(0):
        dim += -Integer(6)*nuinfH(N,H)
    else:
        dim += -Integer(6)*nuinfHreg(N,H)
    return dim//Integer(12)

def dimension_eis_H(G,k):
    N = G.level()
    H = G._list_of_elements_in_H()
    if k%Integer(2) == Integer(1) and N-Integer(1) in H:
        return Integer(0)
    if k%Integer(2) == Integer(0):
        dim = nuinfH(N,H)
    else:
        dim = nuinfHreg(N,H)
    if k == Integer(2):
        dim -= Integer(1)
    return dim

def dimension_new_cusp_forms_H(G,k, p=0):
    """
    INPUT:
        G -- group of the form Gamma_H(N)
        k -- an integer at least 2 (the weight)
        p -- integer (default: 0); if nonzero, compute the p-new subspace.

    OUTPUT:
        Integer

    EXAMPLES:
        sage: from sage.modular.dims import *
        sage: dimension_new_cusp_forms_H(GammaH(33,[2]), 2)
        3
    """
    N = G.level()
    if p==0 or N%p!=0:
        return sum([dimension_cusp_forms_H(H,k) * mumu(N // H.level()) \
                for H in G.divisor_subgroups()])
    else:
        return dimension_new_cusp_forms_H(G,k) - \
               2*dimension_new_cusp_forms_H(G.restrict(N//p),k)


#######

def multgroup(N):
    return [x for x in range(N) if gcd(x,N) == Integer(1)]

def dimension_cusp_forms_fromH(chi,k):
    N = chi.modulus()
    n = chi.order()
    dim = Integer(0)
    for d in divisors(n):
        G = GammaH(N,ker(chi**d))
        dim = dim + moebius(d)*dimension_cusp_formsH(G,k)
    return dim//euler_phi(n)

def dimension_modular_forms_fromH(chi,k):
    N = chi.modulus()
    n = chi.order()
    dim = Integer(0)
    for d in divisors(n):
        G = GammaH(N,ker(chi**d))
        dim = dim + moebius(d)*dimension_modular_formsH(G,k)
    return dim//euler_phi(n)

####################################################################
# Exported Functions
####################################################################

def dimension_new_cusp_forms(X, k=2, p=0):
    """
    Return the dimension of the new (or p-new) subspace of
    cusp forms for the character or group X.

    INPUT:
        X -- integer, congruence subgroup or Dirichlet character
        k -- weight (integer)
        p -- 0 or a prime

    EXAMPLES:
        sage: dimension_new_cusp_forms(100,2)
        1

        sage: dimension_new_cusp_forms(Gamma0(100),2)
        1
        sage: dimension_new_cusp_forms(Gamma0(100),4)
        5

        sage: dimension_new_cusp_forms(Gamma1(100),2)
        141
        sage: dimension_new_cusp_forms(Gamma1(100),4)
        463

        sage: dimension_new_cusp_forms(DirichletGroup(100).1^2,2)
        2
        sage: dimension_new_cusp_forms(DirichletGroup(100).1^2,4)
        8

        sage: sum(dimension_new_cusp_forms(e,3) for e in DirichletGroup(30))
        12
        sage: dimension_new_cusp_forms(Gamma1(30),3)
        12
    """
    if isinstance(X, congroup.CongruenceSubgroup):
        return dimension_new_cusp_forms_group(X,k,p)
    elif isinstance(X, dirichlet.DirichletCharacter):
        return dimension_new_cusp_forms_eps(X,k,p)
    elif isinstance(X, (int,long,Integer)):
        return dimension_new_cusp_forms_gamma0(X,k,p)
    else:
        raise TypeError, "X (=%s) must be an integer, congruence subgroup or Diirichlet character"%X

def dimension_cusp_forms(X, k=2):
    r"""
    The dimension of the space of cusp forms for the given congruence
    subgroup or Dirichlet character.

    INPUT:
        X -- congruence subgroup or Dirichlet character or integer
        k -- weight (integer)

    EXAMPLES:
        sage: dimension_cusp_forms(5,4)
        1

        sage: dimension_cusp_forms(Gamma0(11),2)
        1
        sage: dimension_cusp_forms(Gamma1(13),2)
        2

        sage: dimension_cusp_forms(DirichletGroup(13).0^2,2)
        1
        sage: dimension_cusp_forms(DirichletGroup(13).0,3)
        1

        sage: dimension_cusp_forms(Gamma0(11),2)
        1
        sage: dimension_cusp_forms(Gamma0(11),0)
        0
        sage: dimension_cusp_forms(Gamma0(1),12)
        1
        sage: dimension_cusp_forms(Gamma0(1),2)
        0
        sage: dimension_cusp_forms(Gamma0(1),4)
        0

        sage: dimension_cusp_forms(Gamma0(389),2)
        32
        sage: dimension_cusp_forms(Gamma0(389),4)
        97
        sage: dimension_cusp_forms(Gamma0(2005),2)
        199
        sage: dimension_cusp_forms(Gamma0(11),1)
        0

        sage: dimension_cusp_forms(Gamma1(11),2)
        1
        sage: dimension_cusp_forms(Gamma1(1),12)
        1
        sage: dimension_cusp_forms(Gamma1(1),2)
        0
        sage: dimension_cusp_forms(Gamma1(1),4)
        0

        sage: dimension_cusp_forms(Gamma1(389),2)
        6112
        sage: dimension_cusp_forms(Gamma1(389),4)
        18721
        sage: dimension_cusp_forms(Gamma1(2005),2)
        159201

        sage: dimension_cusp_forms(Gamma1(11),1)
        Traceback (most recent call last):
        ...
        NotImplementedError: computation of dimensions of spaces of weight 1 modular forms not implemented in general.

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_cusp_forms(e,2)
        0
        sage: dimension_cusp_forms(e^2,2)
        1
    """
    if isinstance(X, dirichlet.DirichletCharacter):
        return dimension_cusp_forms_eps(X, k)
    elif isinstance(X, congroup.Gamma0):
        return dimension_cusp_forms_gamma0(X.level(),k)
    elif isinstance(X, (Integer,int,long)):
        return dimension_cusp_forms_gamma0(Integer(X),k)
    elif isinstance(X, congroup.Gamma1):
        return dimension_cusp_forms_gamma1(X.level(),k)
    elif congroup.is_GammaH(X):
        return dimension_cusp_forms_H(X,k)
    elif not isinstance(X, congroup.CongruenceSubgroup):
        raise TypeError, "Argument 1 must be a congruence subgroup or Dirichlet character"
    else:
        raise NotImplementedError, "Computing of dimensions for congruence subgroups besides \
        Gamma0, Gamma1, and GammaH is not yet implemented."

def dimension_eis(X, k=2):
    """
    The dimension of the space of eisenstein series for the given
    congruence subgroup.

    INPUT:
        X -- congruence subgroup or Dirichlet character or integer
        k -- weight (integer)

    EXAMPLES:
        sage: dimension_eis(5,4)
        2

        sage: dimension_eis(Gamma0(11),2)
        1
        sage: dimension_eis(Gamma1(13),2)
        11
        sage: dimension_eis(Gamma1(2006),2)
        3711

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_eis(e,2)
        0
        sage: dimension_eis(e^2,2)
        2

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_eis(e,2)
        0
        sage: dimension_eis(e^2,2)
        2
        sage: dimension_eis(e,13)
        2

        sage: G = DirichletGroup(20)
        sage: dimension_eis(G.0,3)
        4
        sage: dimension_eis(G.1,3)
        6
        sage: dimension_eis(G.1^2,2)
        6

        sage: G = DirichletGroup(200)
        sage: e = prod(G.gens(), G(1))
        sage: e.conductor()
        200
        sage: dimension_eis(e,2)
        4

        sage: dimension_modular_forms(Gamma1(4), 11)
        6
    """
    if k <= 1:
        # TODO
        raise NotImplementedError, "Dimension of weight <= 1 Eisenstein series not yet implemented."
    if isinstance(X, (int,long,Integer)):
        if k%2 == 1: return 0
        d = c0(Integer(X))
        if k==2: d -= 1
        return Z(d)
    if isinstance(X, dirichlet.DirichletCharacter):
        return dimension_eis_eps(X, k)
    if isinstance(X, congroup.Gamma0):
        if k%2 == 1: return 0
        d = c0(X.level())
        if k==2: d -= 1
        return Z(d)
    elif isinstance(X, congroup.Gamma1):
        N = X.level()
        if N == 2 and k%2 == 1:
            d = 0  # level Gamma1(2) and odd weight is a special case.
        elif N == 4 and k%2 == 1:
            d = 2  # level Gamma1(4) and odd weight is a special case.
        else:
            d = c1(N)
        if k==2: d -= 1
        return Z(d)
    elif congroup.is_GammaH(X):
        return dimension_eis_H(X, k)
    elif isinstance(X, congroup.CongruenceSubgroup):
        raise NotImplementedError, "Computation of dimensions for congruence subgroups besides " + \
              "Gamma0 and Gamma1 is not yet implemented."
    else:
        raise TypeError


def dimension_modular_forms(X, k=2):
    r"""
    The dimension of the space of cusp forms for the given congruence
    subgroup (either $\Gamma_0(N)$, $\Gamma_1(N)$, or $\Gamma_H(N)$)
    or Dirichlet character.

    INPUT:
        X -- congruence subgroup or Dirichlet character
        k -- weight (integer)

    EXAMPLES:
        sage: dimension_modular_forms(Gamma0(11),2)
        2
        sage: dimension_modular_forms(Gamma0(11),0)
        1
        sage: dimension_modular_forms(Gamma1(13),2)
        13

        sage: e = DirichletGroup(20).1
        sage: dimension_modular_forms(e,3)
        9
        sage: dimension_cusp_forms(e,3)
        3
        sage: dimension_eis(e,3)
        6
        sage: dimension_modular_forms(11,2)
        2
    """
    if isinstance(X, (int, long, Integer)):
        X = congroup.Gamma0(X)
    elif not isinstance(X, congroup.CongruenceSubgroup) and \
         not isinstance(X, dirichlet.DirichletCharacter):
        raise TypeError, "Argument 1 must be a congruence subgroup or Dirichlet character."
    if k == 0:
        return 1
    if congroup.is_GammaH(X):
        return dimension_modular_forms_H(X, k)
    return dimension_cusp_forms(X, k) + dimension_eis(X, k)

def sturm_bound(level, weight):
    r"""
    Returns the Sturm bound for modules forms with given level and
    weight.

    EXAMPLES:
        sage: sturm_bound(11,2)
        2
        sage: sturm_bound(389,2)
        65
        sage: sturm_bound(1,12)
        1
        sage: sturm_bound(100,2)
        30
        sage: sturm_bound(1,36)
        3

    FURTHER DETAILS: Returns a positive integer~$n$ such that the
    Hecke operators $T_1,\ldots, T_n$ acting on \emph{cusp forms}
    generate the Hecke algebra as a $\Z$-module when the character is
    trivial or quadratic.  Otherwise, $T_1,\ldots,T_n$ generate the
    Hecke algebra at least as a $\Z[\eps]$-module, where $\Z[\eps]$ is
    the ring generated by the values of the Dirichlet character
    $\eps$.  Alternatively, this is a bound such that if two cusp
    forms associated to this space of modular symbols are congruent
    modulo $(\lambda, q^n)$, then they are congruent modulo $\lambda$.

    REFERENCES:
    See the Agashe-Stein appendix to Lario and Schoof, \emph{Some
    computations with Hecke rings and deformation rings},
    Experimental Math., 11 (2002), no. 2, 303-311.  This result
    originated in the paper Sturm, \emph{On the congruence of
    modular forms}, Springer LNM 1240, 275--280, 1987.

    REMARK:
    Kevin Buzzard pointed out to me (William Stein) in Fall 2002
    that the above bound is fine for $\Gamma_1(N)$ with character,
    as one sees by taking a power of $f$.  More precisely, if $f
    \con 0 \pmod{p}$ for first $s$ coefficients, then $f^r \con 0
    \pmod{p}$ for first $sr$ coefficents.  Since the weight of
    $f^r$ is $r\cdot k(f)$, it follows that if $s \geq b$, where
    $b$ is the Sturm bound for $\Gamma_0(N)$ at weight $k(f)$, then
    $f^r$ has valuation large enough to be forced to be $0$ at
    $r*k(f)$ by Sturm bound (which is valid if we choose $r$
    correctly).  Thus $f \con 0 \pmod{p}$.  Conclusion: For
    $\Gamma_1(N)$ with fixed character, the Sturm bound is
    \emph{exactly} the same as for $\Gamma_0(N)$.

    A key point is that we are finding $\Z[\eps]$ generators for the
    Hecke algebra here, not $\Z$-generators.  So if one wants
    generators for the Hecke algebra over $\Z$, this bound must be
    suitably modified (and I'm not sure what the modification is).

    AUTHOR:
        -- William Stein
    """
    return Integer(int(math.ceil((weight * idxG0(level) / ZZ(12)))))

