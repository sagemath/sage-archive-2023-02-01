r"""
Solve S-unit equation x + y = 1

Inspired by work of Tzanakis--de Weger, Baker--Wustholz and Smart, we use the LLL methods in Sage to implement an algorithm that returns all S-unit solutions to the equation $x + y = 1$.

REFERENCES:

- [MR2016]_

- [Sma1995]_

- [Sma1998]_

AUTHORS:

- Alejandra Alvarado, Angelos Koutsianas, Beth Malmskog, Christopher Rasmussen, Christelle Vincent, Mckenzie West (2018-04-25): original version

    
"""


#*****************************************************************************
#       Copyright (C) 2017 Alejandra Alvarado <aalvarado2 at eiu.edu>
#                          Angelos Koutsianas <koutsis.jr at gmail.com>
#                          Beth Malmskog <beth.malmskog at gmail.com>
#                          Christopher Rasmussen <crasmussen at wesleyan.edu>
#                          Christelle Vincent <christelle.vincent at uvm.edu>
#                          Mckenzie West <mckenzierwest at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

from sage.rings.all import Infinity
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RealField, RR
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import is_real_place, refine_embedding
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.finite_rings.integer_mod import mod
from sage.rings.padics.factory import Qp
from sage.combinat.combination import Combinations
from sage.misc.all import prod
from sage.arith.all import factorial
from sage.matrix.constructor import matrix, identity_matrix, vector, block_matrix, zero_matrix
from sage.modules.free_module_element import zero_vector
from itertools import combinations_with_replacement
from sage.functions.log import log
from sage.functions.other import sqrt
from copy import copy
from sage.misc.functional import round

def column_Log(SUK, iota, U, prec=106):
    r"""
    Return the log vector of ``iota``; i.e., the logs of all the valuations

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``iota`` -- an element of ``K``
    - ``U`` -- a list of places (finite or infinite) of ``K``
    - ``prec`` -- (default: 106) the precision of the real field

    OUTPUT:

    The log vector as a list of real numbers

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import column_Log
        sage: K.<xi> = NumberField(x^3-3)
        sage: S = tuple(K.primes_above(3))
        sage: SUK = UnitGroup(K, S=S)
        sage: phi_complex = K.places()[1]
        sage: v_fin = S[0]
        sage: U = [phi_complex, v_fin]
        sage: column_Log(SUK, xi^2, U) # abs tol 1e-29
        [1.464816384890812968648768625966, -2.197224577336219382790490473845]

    REFERENCES:

    - [Sma1995]_ p. 823
    """
    R = RealField(prec)

    return [ R(SUK.number_field().abs_val(v, iota, prec)).log() for v in U]

def c3_func(SUK, prec=106):
    r"""
    Return the constant `c_3` from Smart's 1995 TCDF paper, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``prec`` -- (default: 106) the precision of the real field

    OUTPUT:

    The constant ``c3``, as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c3_func
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))

        sage: c3_func(SUK) # abs tol 1e-29
        0.4257859134798034746197327286726

    .. NOTE::

        The numerator should be as close to 1 as possible, especially as the rank of the `S`-units grows large

    REFERENCES:

    - [Sma1995]_ p. 823

    """

    R = RealField(prec)

    all_places = list(SUK.primes()) + SUK.number_field().places(prec)
    Possible_U = Combinations(all_places, SUK.rank())
    c1 = R(0)
    for U in Possible_U:
        # first, build the matrix C_{i,U}
        columns_of_C = []
        for unit in SUK.fundamental_units():
            columns_of_C.append(column_Log(SUK, unit, U, prec))
        C = matrix(SUK.rank(), SUK.rank(), columns_of_C)
        # Is it invertible?
        if abs(C.determinant()) > 10**(-10):
            poss_c1 = C.inverse().apply_map(abs).norm(Infinity)
            c1 = R(max(poss_c1, c1))
    return R(0.9999999) / (c1*SUK.rank())

def c4_func(SUK, v, A, prec=106):
    r"""
    Return the constant `c_4` from Smart's TCDF paper, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of ``SUK.number_field().places(prec)``)
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` -- (default: 106) the precision of the real field

    OUTPUT:

    The constant ``c4``, as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c4_func
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c4_func(SUK,phi_real,A)
        1.000000000000000000000000000000

        sage: c4_func(SUK,phi_complex,A)
        1.000000000000000000000000000000

        sage: c4_func(SUK,v_fin,A)
        1.000000000000000000000000000000

    REFERENCES:

    - [Sma1995]_ p. 824
    """
    return max(SUK.number_field().abs_val(v, alpha, prec) for alpha in A)

def beta_k(betas_and_ns):
    r"""
    Return a pair `[\beta_k,|beta_k|_v]`, where `\beta_k` has the smallest nonzero valuation in absolute value of the list ``betas_and_ns``

    INPUT:

    - ``betas_and_ns`` -- a list of pairs ``[beta,val_v(beta)]`` outputted from the function where ``beta`` is an element of ``SUK.fundamental_units()``

    OUTPUT:

    The pair ``[beta_k,v(beta_k)]``, where ``beta_k`` is an element of ``K`` and ``val_v(beta_k)`` is a integer

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import beta_k
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: betas = [ [beta, beta.valuation(v_fin)] for beta in SUK.fundamental_units() ]
        sage: beta_k(betas)
        [xi, 1]

    REFERENCES:

    - [Sma1995]_ pp. 824-825
    """
    for pair in betas_and_ns:
        if pair[1].abs() != 0:
            good_pair = pair
            break
    for pair in betas_and_ns:
        if (pair[1].abs()!=0 and pair[1].abs()<good_pair[1].abs()):
            good_pair = pair
    return good_pair

def mus(SUK, v):
    r"""
    Return a list `[\mu]`, for `\mu` defined on pp. 824-825 of TCDF, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a finite place of ``K``

    OUTPUT:

    A list ``[mus]`` where each ``mu`` is an element of ``K``

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import mus
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: mus(SUK, v_fin)
        [xi^2 - 2]

    REFERENCES:

    - [Sma1995]_ pp. 824-825

    """
    betas = SUK.fundamental_units()
    beta_and_ns = [[beta,beta.valuation(v)] for beta in betas]
    if all(pair[1]==0 for pair in beta_and_ns):
        return betas
    else:
        good_pair = beta_k(beta_and_ns)
        temp = [(beta[0]**good_pair[1])*(good_pair[0]**(-beta[1])) for beta in beta_and_ns]
        temp.remove(1)
        return temp

def possible_mu0s(SUK, v):
    r"""
    Return a list `[\mu_0]` of all possible `\mu_0` values defined on pp. 824-825 of TCDF, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a finite place of ``K``

    OUTPUT:

    A list ``[mu0s]`` where each ``mu0`` is an element of ``K``

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import possible_mu0s
        sage: K.<xi> = NumberField(x^3-3)
        sage: S = tuple(K.primes_above(3))
        sage: SUK = UnitGroup(K, S=S)
        sage: v_fin = S[0]

        sage: possible_mu0s(SUK,v_fin)
        [-1, 1]

    .. NOTE::

        `n_0` is the valuation of the coefficient `\alpha_d` of the `S`-unit equation such that `|\alpha_d \tau_d|_v = 1`
        We have set `n_0 = 0` here since the coefficients are roots of unity
        `\alpha_0` is not defined in the paper, we set it to be 1

    REFERENCES:

    - [Sma1995]_ pp. 824-825, but we modify the definition of ``sigma`` (``sigma_tilde``) to make it easier to code

    """
    beta_and_ns = [[beta,beta.valuation(v)] for beta in SUK.fundamental_units()]
    betak, nk = beta_k(beta_and_ns)
    ns = [beta[1] for beta in beta_and_ns if beta[0] != betak]
    betas = [beta[0] for beta in beta_and_ns if beta[0] != betak]
    mu0s = []
    for rs in combinations_with_replacement(range(nk.abs()), len(betas)):
        # n_0 = valuation_v of one of the coefficients of the equation = 0 for x + y = 1 p. 824
        n_rs = zip(ns, rs)
        sigma_tilde = -(sum([n_r[0]*n_r[1] for n_r in n_rs]))
        if sigma_tilde % nk == 0:
            beta_rs = zip(betas, rs)
            temp_prod = prod([beta_r[0]**beta_r[1] for beta_r in beta_rs]) * betak**(sigma_tilde/nk)
            for alpha0 in SUK.roots_of_unity():
                if alpha0*temp_prod not in mu0s:
                    mu0s.append(alpha0*temp_prod)
    return mu0s

def c8_c9_func(SUK, v, A, prec=106):
    r"""
    Return the constants `c_8` and `c_9` from Smart's TCDF paper, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a finite place of ``K`` (a fractional ideal)
    - ``A`` -- the set of the product of the coefficients of the `S`-unit equation with each root of unity of ``K``
    - ``prec`` -- (default: 106) the precision of the real field

    OUTPUT:

    The constants ``c8`` and ``c9``, as real numbers

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c8_c9_func
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: v_fin = K.primes_above(3)[0]
        sage: A = K.roots_of_unity()

        sage: c8_c9_func(SUK, v_fin,A) # abs tol 1e-29
        (4.524941291354698258804956696127e15, 1.621521281297160786545580368612e16)

    REFERENCES:

    - [Sma1995]_ p. 825
    - [Sma1998]_ p. 226, Theorem A.2 for the local constants

    """
    R = RealField(prec)
    num_mus = len(mus(SUK, v))+1
    p = v.smallest_integer()
    f_p = v.residue_class_degree()
    d = SUK.number_field().degree()
    if p == 2:
        local_c2 = Integer(197142) * Integer(36)**num_mus
    elif p%4 == 1:
        local_c2 = Integer(35009) * (Integer(45)/Integer(2))**num_mus
    else:
        local_c2 = Integer(30760) * Integer(25)**num_mus
    x = polygen(SUK.number_field())
    if (p > 2 and not ((x**2+1).is_irreducible())) or (p==2 and not ((x**2+3).is_irreducible())):
        D = d
    else:
        D = 2*d
    l_c3 = (num_mus+1)**(2*num_mus+4) * p**(D * f_p/d) * (f_p*R(p).log())**(-num_mus-1) * D**(num_mus+2)

    def modified_height(SUK, v, D, b):
        #[Sma1998]_ p. 226
        max_log_b = max([phi(b).log().abs() for phi in SUK.number_field().places(prec)])
        return R(max([b.global_height(), max_log_b/(2*R.pi()*D), f_p*R(p).log()/d]))

    mus_prod = prod([modified_height(SUK,v,D,b) for b in mus(SUK,v)])
    local_c3 = R(max([mus_prod*modified_height(SUK, v, D, mu0) for mu0 in possible_mu0s(SUK, v)]))

    l_c3 *= local_c3
    H = max([modified_height(SUK, v, D, alpha) for alpha in mus(SUK, v)+possible_mu0s(SUK, v)])
    if p == 2:
        local_c4 = R(3 * 2**10 * (num_mus+1)**2 * D**2 * H).log()
    else:
        local_c4 = R(2**11 * (num_mus+1)**2 * D**2 * H).log()
    local_c5 = 2 * R(D).log()
    return R(local_c2*l_c3*local_c4), R(local_c2*l_c3*local_c4*local_c5)

def c11_func(SUK, v, A, prec=106):
    r"""
    Return the constant `c_{11}` from Smart's TCDF paper, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of ``SUK.number_field().places(prec)``)
    - ``A`` -- the set of the product of the coefficients of the `S`-unit equation with each root of unity of ``K``
    - ``prec`` -- (default: 106) the precision of the real field

    OUTPUT:

    The constant ``c11``, a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c11_func
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c11_func(SUK, phi_real, A) # abs tol 1e-29
        3.255848343572896153455615423662

        sage: c11_func(SUK, phi_complex, A) # abs tol 1e-29
        6.511696687145792306911230847323

    REFERENCES:

    - [Sma1995]_ p. 825
    """
    R = RealField(prec)
    if is_real_place(v):
        return R(4*c4_func(SUK, v, A, prec)).log() / c3_func(SUK, prec)
    else:
        return 2*R(4*(c4_func(SUK, v, A, prec)).sqrt()).log() / c3_func(SUK, prec)

def c13_func(SUK, v, prec=106):
    r"""
    Return the constant `c_{13}` from Smart's TCDF paper, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- an infinite place of ``K`` (element of ``SUK.number_field().places(prec)``)
    - ``prec`` -- (default: 106) the precision of the real field

    OUTPUT:

    The constant ``c13``, as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import c13_func
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]

        sage: c13_func(SUK, phi_real) # abs tol 1e-29
        0.4257859134798034746197327286726

        sage: c13_func(SUK, phi_complex) # abs tol 1e-29
        0.2128929567399017373098663643363

    It is an error to input a finite place. ::

        sage: phi_finite = K.primes_above(3)[0]
        sage: c13_func(SUK, phi_finite)
        Traceback (most recent call last):
        ...
        TypeError: Place must be infinite


    REFERENCES:

    - [Sma1995]_ p. 825
    """
    try:
        v.codomain()
    except AttributeError:
        raise TypeError('Place must be infinite')
    if is_real_place(v):
        return c3_func(SUK, prec)
    else:
        return c3_func(SUK, prec)/2

def K0_func(SUK, A, prec=106):
    r"""
    Return the constant `K_0` from Smart's TCDF paper, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``A`` -- the set of the products of the coefficients of the `S`-unit equation with each root of unity of ``K``
    - ``prec`` -- (default: 106) the precision of the real field

    OUTPUT:

    The constant ``K0``, a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import K0_func
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: A = K.roots_of_unity()

        sage: K0_func(SUK, A) # abs tol 1e-29
        9.475576673109275443280257946929e17

    REFERENCES:

    - [Sma1995]_ p. 824
    """
    R = RealField(prec)

    def c5_func(SUK, v, R):
        return c3_func(SUK, R.precision()) / (v.residue_class_degree()*R(v.smallest_integer()).log()*v.ramification_index())

    def c6_func(SUK, v, A, R):
        return c4_func(SUK, v, A, R.precision()).log() / (v.residue_class_degree()*R(v.smallest_integer()).log()*v.ramification_index())

    def c7_func(SUK, v, A, R):
        return (c4_func(SUK, v, A, R.precision())).log() / c3_func(SUK, R.precision())

    def c10_func(SUK, v, A, R):
        # [Sma1995]_ p. 824
        e_h = v.ramification_index()
        c_8, c_9 = c8_c9_func(SUK, v, A, R.precision())
        return (2/(e_h*c5_func(SUK, v, R))) * (e_h*c6_func(SUK, v, A, R)+c_9+c_8*(c_8/(e_h*c5_func(SUK, v, R))).log())

    return R(max([c10_func(SUK,v, A, R) for v in SUK.primes()] + [c7_func(SUK,v,A,R) for v in SUK.primes()]))

def K1_func(SUK, v, A, prec=106):
    r"""
    Return the constant `K_1` from Smart's TCDF paper, [Sma1995]_

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``v`` -- an infinite place of ``K`` (element of ``SUK.number_field().places(prec)``)
    - ``A`` -- a list of all products of each potential ``a``, ``b`` in the $S$-unit equation ``ax + by + 1 = 0`` with each root of unity of ``K``
    - ``prec`` -- (default: 106) the precision of the real field

    OUTPUT:

    The constant ``K1,`` a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import K1_func
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K, S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: K1_func(SUK, phi_real, A)
        4.396386097852707394927181864635e16

        sage: K1_func(SUK, phi_complex, A)
        2.034870098399844430207420286581e17

    REFERENCES:

    - [Sma1995]_ p. 825

    """
    R = RealField(prec)

    #[Sma1995]_ p. 825
    if is_real_place(v):
        c11 = R(4*c4_func(SUK, v, A, prec)).log() / c3_func(SUK, prec)
    else:
        c11 = 2*( R(4*(c4_func(SUK,v, A, prec)).sqrt()).log() ) / c3_func(SUK, prec)

    #[Sma1995]_ p. 825
    if is_real_place(v):
        c12 = R(2 * c4_func(SUK, v, A, prec))
    else:
        c12 = R(2 * c4_func(SUK, v, A, prec).sqrt())

    #[Sma1998]_ p. 225, Theorem A.1
    d = SUK.number_field().degree()
    t = SUK.rank()
    Baker_C = R(18 * factorial(t+2) * (t+1)**(t+2) * (32*d)**(t+3) * R(2*(t+1) * d).log())

    def hprime(SUK, alpha, v):
        #[Sma1998]_ p. 225
        return R(max(alpha.global_height(), 1/SUK.number_field().degree(), v(alpha).log().abs()/SUK.number_field().degree()))

    #[Sma1995]_ p. 825 and [Sma1998]_ p. 225, Theorem A.1
    c14 = Baker_C * prod([hprime(SUK, alpha, v) for alpha in SUK.gens_values()])

    #[Sma1995]_ p. 825
    c15 = 2 * (c12.log()+c14*((SUK.rank()+1)*c14/c13_func(SUK, v, prec)).log()) / c13_func(SUK, v, prec)

    return max([c11, c15])

def minimal_vector(A, y, prec=106):
    r"""
    INPUT:

    - ``A`` : a square n by n non-singular integer matrix whose rows generate a lattice `\mathcal L`
    - ``y`` : a row (1 by n) vector with integer coordinates
    - ``prec`` : (default: 106) precision of real field

    OUTPUT:

    A lower bound for the square of

    .. MATH::

        \ell (\mathcal L,\vec y) =
        \begin{cases}
        \displaystyle\min_{\vec x\in\mathcal L}\Vert\vec x-\vec y\Vert &, \vec y\not\in\mathcal L. \\
        \displaystyle\min_{0\neq\vec x\in\mathcal L}\Vert\vec x\Vert &,\vec y\in\mathcal L.
        \end{cases}`

    ALGORITHM:

    The algorithm is based on V.9 and V.10 of [Sma1998]_

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import minimal_vector
        sage: B = matrix(ZZ, 2, [1,1,1,0])
        sage: y = vector(ZZ, [2,1])
        sage: minimal_vector(B, y)
        1/2

    ::

        sage: B = random_matrix(ZZ, 3)
        sage: B #random
        [-2 -1 -1]
        [ 1  1 -2]
        [ 6  1 -1]
        sage: y = vector([1, 2, 100])
        sage: minimal_vector(B, y) #random
        15/28
    """
    if A.is_singular():
        raise ValueError('The matrix A is singular')

    R = RealField(prec)

    n = len(y)
    c1 = 2**(n-1)
    ALLL = A.LLL()
    ALLLinv = ALLL.inverse()
    ybrace = [R(a-round(a)).abs() for a in y * ALLLinv if (a-round(a)) != 0]

    if len(ybrace) == 0:
        return (ALLL.rows()[0].norm())**2 / c1
    else:
        sigma = ybrace[len(ybrace)-1]
        return ((ALLL.rows()[0].norm())**2 * sigma) / c1

def reduction_step_real_case(place, B0, G, c7):
    r"""
    INPUT:

    - ``place`` -- (ring morphism) an infinite place of a number field `K`
    - ``B0`` -- the initial bound
    - ``G`` -- a set of generators of the free part of the group
    - ``c7`` -- a positive real number

    OUTPUT:

    A tuple consisting of:

    1. a new upper bound, an integer
    2. a boolean value, ``True`` if we have to increase precision, otherwise ``False``

    .. NOTE::

        The constant ``c7`` in the reference page 137

    REFERENCES:

    - [Sma1998]_

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import reduction_step_real_case
        sage: K.<a> = NumberField(x^3-2)
        sage: SK = sum([K.primes_above(p) for p in [2,3,5]],[])
        sage: G = [g for g in K.S_unit_group(S=SK).gens_values() if g.multiplicative_order()==Infinity]
        sage: p1 = K.real_places(prec=300)[0]
        sage: reduction_step_real_case(p1, 10**10, G, 2)
        (58, False)
    """
    prec = place.codomain().precision()
    n = len(G)

    def e_s_real(a, place):
        if place(a) < 0:
            return (-1) * a
        else:
            return a
    Glog = [log(place(e_s_real(g, place))) for g in G]
    if len([1 for g in G if place(e_s_real(g, place)).is_zero()]) > 0:
        return 0, True

    #We choose the initial value of C such that the vector v not to have 0 everywhere
    C = round(max([1/abs(l) for l in Glog if l != 0]) + 1)

    #if the precision we have is not high enough we have to increase it and evaluate c7 again
    if prec < log(C)/log(2):
        return 0, True

    S = (n-1) * (B0)**2
    T = (1 + n*B0) / 2
    finish = False
    while not finish:
        A = copy(identity_matrix(ZZ, n))
        v = vector([round(g*C) for g in Glog])

        if v[n-1] == 0: #we replace the last element of v with an other non zero
            k = [i for i,a in enumerate(v) if not a.is_zero()][0]
            v[n-1] = v[k]
            v[k] = 0
        A[n-1] = v

        #We have to work with rows because of the .LLL() function

        A = A.transpose()
        y = copy(zero_vector(ZZ, n))
        l = minimal_vector(A, y)

        #On the following lines I apply Lemma VI.1 from Smart's book page 83

        if l < T**2 + S:
            C = 2 * C
            #Again if our precision is not high enough
            if prec < log(C)/log(2):
                return 0, True
        else:
            if sqrt(l-S) - T > 0:
                return round((log(C*2)-log(sqrt(l-S)-T)) / c7), False
            else:
                return B0, False

def reduction_step_complex_case(place,B0,G,g0,c7):
    r"""
    INPUT:

    - ``place`` -- (ring morphism) an infinite place of a number field `K`
    - ``B0`` -- the initial bound
    - ``G`` -- a set of generators of the free part of the group
    - ``g0`` -- an element of the torsion part of the group
    - ``c7`` -- a positive real number

    OUTPUT:

    A tuple consisting of:

    1. a new upper bound, an integer
    2. a boolean value, ``True`` if we have to increase precision, otherwise ``False``

    .. NOTE::

        The constant ``c7`` in the reference page 138

    REFERENCES:

    See [Sma1998]_.

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import reduction_step_complex_case
        sage: K.<a> = NumberField([x^3-2])
        sage: SK = sum([K.primes_above(p) for p in [2,3,5]],[])
        sage: G = [g for g in K.S_unit_group(S=SK).gens_values() if g.multiplicative_order()==Infinity]
        sage: p1 = K.places(prec=100)[1]
        sage: reduction_step_complex_case(p1, 10^5, G, -1, 2)
        (17, False)
    """
    prec = place.codomain().precision()
    R = RealField(prec)
    n = len(G)
    Glog_imag = [place(g).log().imag_part() for g in G]
    Glog_real = [place(g).log().real_part() for g in G]
    Glog_imag = Glog_imag + [2 * R.pi()]
    Glog_real = Glog_real + [0]
    a0log_imag = place(-g0).log().imag_part()
    a0log_real = place(-g0).log().real_part()

    #the case when the real part is 0 for all log(a_i)

    pl = refine_embedding(place)
    if len([g for g in G if (pl(g).abs()-1).abs() > 2**(-place.codomain().precision())])  ==  0:

        #we have only imaginary numbers and we are in case 2 as Smart's book says on page 84

        C = ZZ(1) #round(min((B0**n/100),max([1/l.abs() for l in Glog_imag if l != 0]+[0])))+1
        S = n * B0**2

        #if the precision we have is not high enough we have to increase it and evaluate c7 again
        # if precision < log(C)/log(2):
        #     return 0,True

        T = ((n+1)*B0 + 1) / 2
        finish = False
        while not finish:
            A = copy(identity_matrix(ZZ, n+1))
            v = vector([round(g * C) for g in Glog_imag])

            if v[n] == 0:
                #we replace the last element of v with an other non zero

                k = [i for i,a in enumerate(v) if not a.is_zero()][0]
                v[n] = v[k]
                v[k] = 0
            A[n] = v

            if A.is_singular():
                C = ZZ(2*C)
            else:
                #We have to work with rows because of the .LLL() function

                A = A.transpose()
                y = copy(zero_vector(ZZ, n+1))
                y[n] = (-1) * round(a0log_imag*C)
                l = minimal_vector(A, y)

                #On the following lines I apply Lemma VI.1 of the reference page 83

                if l < T**2+S:
                    C = ZZ(2*C)

                    #The same as above if for the new C the precision is low
                    if prec < log(C)/log(2):
                        return 0, True
                else:
                    Bnew = round((log(C * 2)-log(sqrt(l-S)-T)) / c7)
                    finish = True
                    if mod(y[n], A[n,n]) == 0:
                        return max(Bnew,(y[n]/A[n,n]).abs()), False
                    else:
                        return Bnew, False

    else:

        #the case when the real part is not 0 for all log(a_i)
        C = R(1)
        S = (n-1) * B0**2
        T = ((n+1)*B0+1) / R(2).sqrt()
        finish = False

        #we are relabeling the Glog_real and Glog_imag s.t. the condition Real(a_n)*Im(a_(n-1))-Real(a_(n-1))*Im(a_n)!=0 to be satisfied. See page 84 of the reference.

        k = [i for i in range(len(Glog_real)) if Glog_real[i] != 0][0]
        a = Glog_real[k]
        Glog_real[k] = Glog_real[n-1]
        Glog_real[n-1] = a

        a = Glog_imag[k]
        Glog_imag[k] = Glog_imag[n-1]
        Glog_imag[n-1] = a

        while not finish:

            A = copy(identity_matrix(ZZ, n+1))
            #return [g * C for g in Glog_imag]
            A[n-1] = vector([round(g*C) for g in Glog_real])
            A[n] = vector([round(g*C) for g in Glog_imag])

            if A.is_singular():
                C *= 2
            else:
                #On the following lines I apply Lemma VI.2 of the reference page 85

                A = A.transpose()
                y = copy(zero_vector(ZZ, n+1))
                y[n] = (-1) * round(a0log_imag*C)
                y[n-1] = (-1) * round(a0log_real*C)
                l = minimal_vector(A,y)


                if l <= T**2 + S:
                    C *= 2
                    #The same as above if for the new C the precision is low
                    if prec < log(C)/log(2):
                        return 0, True
                else:
                    Bnew = round( ((C * 2).log() - ((l-S).sqrt()-T).log()) / c7 )

                    #we take into account the second case of the theorem VI.2 of the reference page 85

                    M = matrix(ZZ, 2, [A[n-1,n-1],A[n-1,n],A[n,n-1],A[n,n]])
                    b = vector(ZZ, 2, [-y[n-1],-y[n]])
                    if M.determinant() == 1 or M.determinant() == -1:
                        x = M.inverse() * b
                        return max(Bnew, x[0].abs(), x[1].abs()), False
                    else:
                        return Bnew, False

def cx_LLL_bound(SUK, A, prec=106):
    r"""
    Return the maximum of all of the `K_1`'s as they are LLL-optimized for each infinite place `v`

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``A`` -- a list of all products of each potential ``a``, ``b`` in the `S`-unit equation ``ax + by + 1 = 0`` with each root of unity of ``K``
    - ``prec`` -- precision of real field (default 106)

    OUTPUT:

    A bound for the exponents at the infinite place, as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import cx_LLL_bound
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: A = K.roots_of_unity()

        sage: cx_LLL_bound(SUK,A)
        22
    """
    cx_LLL = 0
    #initialize a bound, a bad guess, as we iterate over the places of the number field, we will replace its value with the largest complex LLL bound we've found across the places
    for v in SUK.number_field().places(prec=prec):
        prec_v = prec
        c13_LLL = c13_func(SUK, v, prec_v)
        cx_bound = K1_func(SUK, v, A, prec_v)
        #cx_bound is the LLL bound according to this place, it will be replaced as LLL gives us smaller bounds
        if is_real_place(v):
            new_bound,inc_prec = reduction_step_real_case(v, cx_bound, SUK.fundamental_units(), c13_LLL)
            while inc_prec:
                v = refine_embedding(v)
                c13_LLL = c13_func(SUK, v, prec_v)
                cx_bound = K1_func(SUK, v, A, prec_v)
                new_bound, inc_prec = reduction_step_real_case(v, cx_bound, SUK.fundamental_units(), c13_LLL)
            counter = 0
            while (cx_bound-new_bound).abs() > .01*cx_bound and counter < 15:
                #We fear a loop that is not convergent, this is the purpose of the counter
                #Repeat complex LLL until we get essentially no change from it
                cx_bound = min(cx_bound,new_bound)
                new_bound, inc_prec = reduction_step_real_case(v, cx_bound, SUK.fundamental_units(),c13_LLL)
                while inc_prec:
                    v = refine_embedding(v)
                    c13_LLL = c13_func(SUK,v,prec_v)
                    new_bound, inc_prec = reduction_step_real_case(v, cx_bound, SUK.fundamental_units(), c13_LLL)
                counter += 1
        else:
            prec_v = prec
            new_bound, inc_prec = reduction_step_complex_case(v, cx_bound, SUK.fundamental_units(), SUK.zeta(), c13_LLL)
            while inc_prec:
                v = refine_embedding(v)
                c13_LLL = c13_func(SUK, v, prec_v)
                cx_bound = K1_func(SUK, v, A, prec_v)
                new_bound, inc_prec = reduction_step_complex_case(v, cx_bound, SUK.fundamental_units(), SUK.zeta(), c13_LLL)
            counter = 0
            while (cx_bound-new_bound).abs() > .01*cx_bound and counter < 15:
                #We fear a loop that is not convergent, this is the purpose of the counter
                #Repeat complex LLL until we get essentially no change from it
                cx_bound = min(cx_bound, new_bound)
                new_bound, inc_prec = reduction_step_complex_case(v, cx_bound, SUK.fundamental_units(), SUK.zeta(), c13_LLL)
                while inc_prec:
                    v = refine_embedding(v)
                    c13_LLL = c13_func(SUK, v, prec_v)
                    new_bound, inc_prec = reduction_step_complex_case(v, cx_bound, SUK.fundamental_units(), SUK.zeta(), c13_LLL)
                counter += 1

        cx_bound = min(cx_bound, new_bound)
        #for this place the complex LLL bound is cx_bound
        cx_LLL = max(cx_bound, cx_LLL)
        #compare this value with the complex LLL bounds we have found for the previous places, if it is bigger, replace that bound
    return cx_LLL


def log_p(a, prime, prec):
    r"""
    INPUT:

    - ``a`` -- an element of a number field `K`
    - ``prime`` -- a prime ideal of the number field `K`
    - ``prec`` -- a positive integer

    OUPUT:

    An element of `K` which is congruent to the ``prime``-adic logarithm of ``a`` with respect to ``prime`` modulo ``p^prec``, where ``p`` is the rational prime below ``prime``

    .. NOTE::

        Here we take into account the other primes in `K` above `p` in order to get coefficients with small values

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import log_p
        sage: K.<a> = NumberField(x^2+14)
        sage: p1 = K.primes_above(3)[0]
        sage: p1
        Fractional ideal (3, a + 1)
        sage: log_p(a+2, p1, 20)
        8255385638/3*a + 15567609440/3

    ::

        sage: K.<a> = NumberField(x^4+14)
        sage: p1 = K.primes_above(5)[0]
        sage: p1
        Fractional ideal (5, a + 1)
        sage: log_p(1/(a^2-4), p1, 30)
        -42392683853751591352946/25*a^3 - 113099841599709611260219/25*a^2 -
        8496494127064033599196/5*a - 18774052619501226990432/25
    """
    if a == 0:
        raise ValueError('a is the zero element')

    if a.valuation(prime) != 0:
        raise ValueError('The valuation of a with respect to prime is not zero')

    K = prime.ring()
    p = prime.smallest_integer()

    #In order to get an approximation with small coefficients we have to take into account the other primes above p
    #with negative valuation.  For example, say prime2 is another (principal ideal) prime above p, and a=(unit)(prime2)^(-k) for some unit and k
    #a postive integer, and let tilde(a):=a(prime2)^k.  Then log_p(a)=log_p(tilde(a))-k(log_p(prime2)), where the series representations
    #of these two logs will have smaller coeffiecients.

    primes = [(-(a.valuation(pr)),pr) for pr in K.primes_above(p) if a.valuation(pr) < 0]
    local_terms = []

    for (val, pr) in primes:
        #for its pair in primes we find an element in K such that it is divisible only by pr and not by any other ideal above p. Then we take this element in the correct exponent

        if pr.is_principal():
            local_terms.append(pr.gens_reduced()[0]**val)
        else:
            local_terms.append(pr.gens()[1]**val)

    return log_p_series_part(a*prod(local_terms), prime, prec) - sum([log_p_series_part(b, prime, prec) for b in local_terms])

def log_p_series_part(a, prime, prec):
    r"""
    INPUT:

    - ``a`` -- an element of a number field `K`
    - ``prime`` -- a prime ideal of the number field `K`
    - ``prec`` -- a positive integer

    OUPUT:

    The ``prime``-adic logarithm of ``a`` and accuracy ``p^prec``, where ``p`` is the rational prime below ``prime``

    ALGORITHM:

    The algorithm is based on the algorithm on page 30 of [Sma1998]_

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import log_p_series_part
        sage: K.<a> = NumberField(x^2-5)
        sage: p1 = K.primes_above(3)[0]
        sage: p1
        Fractional ideal (3)
        sage: log_p_series_part(a^2-a+1, p1, 30)
        120042736778562*a + 263389019530092

    ::

        sage: K.<a> = NumberField(x^4+14)
        sage: p1 = K.primes_above(5)[0]
        sage: p1
        Fractional ideal (5, a + 1)
        sage: log_p_series_part(1/(a^2-4), p1, 30)
        5628940883264585369224688048459896543498793204839654215019548600621221950915106576555819252366183605504671859902129729380543157757424169844382836287443485157589362653561119898762509175000557196963413830027960725069496503331353532893643983455103456070939403472988282153160667807627271637196608813155377280943180966078/1846595723557147156151786152499366687569722744011302407020455809280594038056223852568951718462474153951672335866715654153523843955513167531739386582686114545823305161128297234887329119860255600972561534713008376312342295724191173957260256352612807316114669486939448006523889489471912384033203125*a^2 + 2351432413692022254066438266577100183514828004415905040437326602004946930635942233146528817325416948515797296867947688356616798913401046136899081536181084767344346480810627200495531180794326634382675252631839139904967037478184840941275812058242995052383261849064340050686841429735092777331963400618255005895650200107/1846595723557147156151786152499366687569722744011302407020455809280594038056223852568951718462474153951672335866715654153523843955513167531739386582686114545823305161128297234887329119860255600972561534713008376312342295724191173957260256352612807316114669486939448006523889489471912384033203125
    """
    if a.valuation(prime) != 0:
        raise ValueError('The valuation of a with respect to prime is not zero')
    K = prime.ring()
    g = K.gen()
    p = prime.smallest_integer()
    f = prime.residue_class_degree()
    e = prime.absolute_ramification_index()
    q = p**f - 1
    R = RealField(prec)

    divisor = q.divisors()
    order = min([d for d in divisor if (a**d - 1).valuation(prime) > 0])
    gamma= a**order
    t = 0
    while (gamma-1).valuation(prime) <= e:
        t += 1
        gamma = gamma**p
    prec += t
    #since later we divide by p^t, we must increase the precision by t at this point.
    m = (gamma-1).valuation(prime) / e
    n = Integer(1)
    step = 10 ** (R(prec).log()/R(10).log()).floor()
    while n < (R(n).log()/R(p).log() + prec)/m:
        n += step
    #could use smaller stepsize to get actual smallest integer n, however this seems to run faster.
    w = R(R(prec).log()/R(p).log()).floor()
    gamma = sum([ZZ(gi%(p**(prec+w)))* g**i if gi.valuation(p) >= 0 else ZZ((gi * p**(-gi.valuation(p)))%(p**(prec+w-gi.valuation(p)))) * p**(gi.valuation(p)) * g**i for i,gi in enumerate(gamma) if gi != 0])


    beta = 0
    delta = 1-gamma
    for i in range(1, n+1):
        beta -= delta / i
        delta *= (1 - gamma)
        delta = sum([ZZ(di%(p**(prec+w)))* g**e if di.valuation(p) >= 0 else ZZ((di * p**(-di.valuation(p)))%(p**(prec+w-di.valuation(p)))) * p**(di.valuation(p)) * g**e for e,di in enumerate(delta) if di != 0],0)
    beta = beta/(order*p**t)

    #we try to make the coefficients small

    logp = 0
    for i,b in enumerate(beta.list()):
        val = b.valuation(p)
        if val < 0:
            t = b * p**(-val)
            t = ZZ(mod(t, p**(prec-val)))
            t = t * p**val
        else:
            t = ZZ(mod(b, p**prec))
        logp = logp + t * g**i

    return logp

def defining_polynomial_for_Kp(prime, prec=106):
    r"""
    INPUT:

    - ``prime`` -- a prime ideal of a number field `K`
    - ``prec`` -- a positive natural number (default 106)

    OUTPUT:

    A polynomial with integer coefficients that is equivalent ``mod p^prec`` to the defining polynomial of the completion of `K` associate to the defining polynomial of `K`

    .. NOTE::

        `K` has to be an absolute extension

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import defining_polynomial_for_Kp
        sage: K.<a> = QuadraticField(2)
        sage: p2 = K.prime_above(7); p2
        Fractional ideal (-2*a + 1)
        sage: defining_polynomial_for_Kp(p2, 10)
        x + 266983762

    ::

        sage: K.<a> = QuadraticField(-6)
        sage: p2 = K.prime_above(2); p2
        Fractional ideal (2, a)
        sage: defining_polynomial_for_Kp(p2, 100)
        x^2 + 6
        sage: p5 = K.prime_above(5); p5
        Fractional ideal (5, a + 2)
        sage: defining_polynomial_for_Kp(p5, 100)
        x + 3408332191958133385114942613351834100964285496304040728906961917542037
    """
    K = prime.ring()
    if not K.is_absolute():
        raise ValueError('The number field is not an absolute extension')

    theta = K.gen()
    f = K.defining_polynomial()
    p = prime.smallest_integer()
    e = prime.absolute_ramification_index()

    find = False
    N = prec
    while find == False:
        RQp = Qp(p, prec=N, type='capped-rel', print_mode='series')

        #We factor f in Integers(p**(precision)) using the factorization in Qp

        g = f.change_ring(RQp)
        factors = g.factor()

        #We are going to find which factor of f is related to the prime ideal 'prime'

        L = [factors[i][0].change_ring(ZZ) for i in range(len(factors))]
        A = [g for g in L if (g(theta)).valuation(prime) >= e*N/2];

        if len(A) == 1:
            return A[0].change_ring(Integers(p**prec)).change_ring(ZZ)
        else:
            N += 1

def embedding_to_Kp(a, prime, prec):
    r"""
    INPUT:

    - ``a`` -- an element of a number field `K`
    - ``prime`` -- a prime ideal of `K`
    - ``prec`` -- a positive natural number

    OUTPUT:

    An element of `K` that is equivalent to ``a`` modulo ``p^(prec)`` and the generator of `K` appears with exponent less than `e \cdot f`, where ``p`` is the rational prime below ``prime`` and `e,f` are the ramification index and residue degree, respectively.

    .. NOTE::

        `K` has to be an absolute extension

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import embedding_to_Kp
        sage: K.<a> = QuadraticField(17)
        sage: p = K.prime_above(13); p
        Fractional ideal (-a + 2)
        sage: embedding_to_Kp(a-3, p, 15)
        -20542890112375827

    ::

        sage: K.<a> = NumberField(x^4-2)
        sage: p = K.prime_above(7); p
        Fractional ideal (-a^2 + a - 1)
        sage: embedding_to_Kp(a^3-3, p, 15)
        -1261985118949117459462968282807202378
    """
    K = prime.ring()
    if not K.is_absolute():
        raise ValueError('K has to be an absolute extension')

    g = defining_polynomial_for_Kp(prime, prec)
    gen = K.gen()
    g = g.change_ring(QQ)
    f = K(a).lift()

    return K( sum([b*gen**j for j,b in enumerate(f.mod(g))]) )

def p_adic_LLL_bound_one_prime(prime, B0, M, M_logp, m0, c3, prec=106):
    r"""
    INPUT:

    - ``prime`` -- a prime ideal of a number field `K`
    - ``B0`` -- the initial bound
    - ``M`` -- a list of elements of `K`, the `\mu_i`'s from Lemma IX.3 of [Sma1998]_
    - ``M_logp`` -- the p-adic logarithm of elements in `M`
    - ``m0`` -- an element of `K`, this is `\mu_0` from Lemma IX.3 of [Sma1998]_
    - ``c3`` -- a positive real constant
    - ``prec`` -- the precision of the calculations (default 106)

    OUTPUT:

    A pair consisting of:

    1. a new upper bound, an integer
    2. a boolean value, ``True`` if we have to increase precision, otherwise ``False``

    .. NOTE::

        The constant `c_5` is the constant `c_5` at the page 89 of [Sma1998]_ which is equal to the constant `c_{10}` at the page 139 of [Sma1995]_.
        In this function, the `c_i` constants are in line with [Sma1998]_, but generally differ from the constants in [Sma1995]_ and other parts of this code.

    EXAMPLES:

    This example indictes a case where we must increase precision::

        sage: from sage.rings.number_field.S_unit_solver import p_adic_LLL_bound_one_prime
        sage: prec = 50
        sage: K.<a> = NumberField(x^3-3)
        sage: S = tuple(K.primes_above(3))
        sage: SUK = UnitGroup(K, S=S)
        sage: v = S[0]
        sage: A = SUK.roots_of_unity()
        sage: K0_old = 9.4755766731093e17
        sage: Mus = [a^2 - 2]
        sage: Log_p_Mus = [185056824593551109742400*a^2 + 1389583284398773572269676*a + 717897987691852588770249]
        sage: mu0 = K(-1)
        sage: c3_value = 0.42578591347980
        sage: m0_Kv_new, increase_precision = p_adic_LLL_bound_one_prime(v, K0_old, Mus, Log_p_Mus, mu0, c3_value, prec)
        sage: m0_Kv_new
        0
        sage: increase_precision
        True

    And now we increase the precision to make it all work::

        sage: prec = 106
        sage: K0_old = 9.475576673109275443280257946930e17
        sage: Log_p_Mus = [1029563604390986737334686387890424583658678662701816*a^2 + 661450700156368458475507052066889190195530948403866*a]
        sage: c3_value = 0.4257859134798034746197327286726
        sage: m0_Kv_new, increase_precision = p_adic_LLL_bound_one_prime(v, K0_old, Mus, Log_p_Mus, mu0, c3_value, prec)
        sage: m0_Kv_new
        476
        sage: increase_precision
        False
    """
    if any(g.valuation(prime) != 0 for g in M+[m0]):
        raise ValueError('There is an element with non zero valuation at prime')

    K = prime.ring()
    p = prime.smallest_integer()
    f = prime.residue_class_degree()
    e = prime.absolute_ramification_index()
    R = RealField(prec)
    c5 = c3 / (f*e*R(p).log())
    theta = K.gen()

    #if M is empty then it is easy to give an upper bound
    if len(M) == 0:
        if m0 != 1:
            return R(max(R(p).log()*f*(m0-1).valuation(prime)/c3, 0)).floor(), False
        else:
            return 0, False
    #we evaluate the p-adic logarithms of m0 and we embed it in the completion of K with respect to prime

    m0_logp = log_p(m0, prime, prec)
    m0_logp = embedding_to_Kp(m0_logp, prime, prec)
    n = len(M_logp)
    #Below we implement paragraph VI.4.2 of [Smart], pages 89-93

    #we evaluate the order of discriminant of theta

    Theta = [theta**i for i in range(K.absolute_degree())]
    ordp_Disc = (K.disc(Theta)).valuation(p)
    #Let's check the mathematics here
    #We evaluate lambda

    c8 = min([min([a.valuation(p) for a in g]) for g in M_logp])
    lam = p**c8

    #we apply lemma VI.5 of [Smart] page 90
    #c6 is 0 here because we seek to solve the equation x+y=1, so our set A
    #is contained in the roots of unity of K

    low_bound = round(1/c5)
    for a in m0_logp:
        if a != 0:
            if c8 > a.valuation(p):
                B1 = (c8 + ordp_Disc/2) / c5
                if B1 > low_bound:
                    return RR(B1).floor(), False
                else:
                    return low_bound, False

    c8 = min([a.valuation(p) for a in m0_logp] + [c8])
    B = [g/lam for g in M_logp]
    b0 = m0_logp / lam
    c9 = c8 + ordp_Disc/2

    #We evaluate 'u' and we construct the matrix A

    m = e * f
    u = 1
    while True:
        if u > (prec * R(2).log()) / R(p).log():
            return 0, True

        #We construct the matrix A as a block matrix

        A11 = copy(identity_matrix(ZZ, n))
        A12 = copy(zero_matrix(ZZ, n, m))
        A21 = copy(zero_matrix(ZZ, n, m))
        A22 = p**u * copy(identity_matrix(ZZ, m))
        for i,b in enumerate(B):
            A21[i] = vector([mod(b[j],p**u) for j in range(m)])
        A = block_matrix( [[A11,A12], [A21.transpose(),A22]] )

        y = copy(zero_vector(ZZ, n+m))
        for i in range(m):
            y[i+n] = -mod(b0[i], p**u)
        #This refers to c10 from Smart
        c10squared = minimal_vector(A.transpose(), y)
        if c10squared > n * B0**2:
            B2 = (u+c9) / c5
            if B2 > low_bound:
                return R(B2).floor(),False
            else:
                return low_bound,False
        else:
            u += 1

def p_adic_LLL_bound(SUK, A, prec=106):
    r"""
    Return the maximum of all of the `K_0`'s as they are LLL-optimized for each finite place `v`

    INPUT:

    - ``SUK`` -- a group of `S`-units
    - ``A`` -- a list of all products of each potential ``a``, ``b`` in the `S`-unit equation ``ax + by + 1 = 0`` with each root of unity of ``K``
    - ``prec``-- precision for p-adic LLL calculations (default: 106)

    OUTPUT:

    A bound for the max of exponents in the case that extremal place is finite (see [Sma1995]_) as a real number

    EXAMPLES::

        sage: from sage.rings.number_field.S_unit_solver import p_adic_LLL_bound
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: A = SUK.roots_of_unity()
        sage: prec = 100
        sage: p_adic_LLL_bound(SUK,A, prec)
        89
    """
    S = SUK.primes()
    K0_old = K0_func(SUK, A, prec)
    LLL_K0_by_finite_place = []
    for i,v in enumerate(S):
        #Kv_old = K0_by_finite_place[0]
        Mus0 = possible_mu0s(SUK, v)
        Mus = mus(SUK, v)
        Log_p_Mus = [log_p(a, v, prec) for a in Mus]
        local_prec = prec
        val = 0
        for m0 in Mus0:
            m0_Kv_old = K0_old
            m0_Kv_new, increase_precision = p_adic_LLL_bound_one_prime(v, m0_Kv_old, Mus, Log_p_Mus, m0, c3_func(SUK, local_prec), local_prec)
            while increase_precision:
                local_prec *= 2
                Log_p_Mus = [log_p(a, v, local_prec) for a in Mus]
                Log_p_Mus = [embedding_to_Kp(a, v, prec) for a in Log_p_Mus]
                m0_Kv_new, increase_precision = p_adic_LLL_bound_one_prime(v, m0_Kv_old, Mus, Log_p_Mus, m0, c3_func(SUK, local_prec), local_prec)

            while m0_Kv_new < m0_Kv_old:
                m0_Kv_old = m0_Kv_new
                m0_Kv_new, increase_precision = p_adic_LLL_bound_one_prime(v, m0_Kv_old, Mus, Log_p_Mus, m0, c3_func(SUK,local_prec), local_prec)
                while increase_precision:
                    local_prec *= 2
                    Log_p_Mus = [log_p(a, v, local_prec) for a in Mus]
                    m0_Kv_new,increase_precision = p_adic_LLL_bound_one_prime(v, m0_Kv_old, Mus, Log_p_Mus, m0, c3_func(SUK, local_prec), local_prec)

            if m0_Kv_old > val:
                val = m0_Kv_old

        LLL_K0_by_finite_place.append(val)
    return max(LLL_K0_by_finite_place)
