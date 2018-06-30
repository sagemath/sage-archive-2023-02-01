r"""
Solve S-unit equation x + y = 1

Inspired by work of Tzanakis--de Weger, Baker--W{\"{u}}stholz and Smart, we use the LLL methods in Sage to implement an algorithm that returns all S-unit solutions to the equation x + y = 1.
This is a preliminary patch and does not have all functions to solve the equation.

REFERENCES:

- [MR2016]_

- [Sma1995]_

- [Sma1998]_

AUTHORS:

- Alejandra Alvarado, Angelos Koutsianas, Beth Malmskog, Christopher Rasmussen, Christelle Vincent, Mckenzie West (2017-01-10): original version

    
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
from sage.rings.real_mpfr import RealField
from sage.rings.number_field.number_field import is_real_place
from sage.combinat.combination import Combinations
from sage.misc.all import prod
from sage.arith.all import factorial
from sage.matrix.constructor import Matrix
from itertools import combinations_with_replacement

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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: U = [phi_complex, v_fin]
        sage: column_Log(SUK, xi^2, U) # abs tol 1e-29
        [1.464816384890812968648768625966, -2.197224577336219382790490473845]

    REFERENCES:

    - [Sma1995]_ p. 823
    """
    R = RealField(prec)

    return [ R(SUK.number_field().abs_val(v,iota,prec)).log() for v in U]

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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))

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
    c1 = 0
    for U in Possible_U:
        # first, build the matrix C_{i,U}
        columns_of_C = []
        for unit in SUK.fundamental_units():
            columns_of_C.append( column_Log(SUK, unit, U, prec) )
        C = Matrix(SUK.rank(), SUK.rank(), columns_of_C)
        # Is it invertible?
        if abs(C.determinant()) > 10**(-10):
            poss_c1 = C.inverse().apply_map(abs).norm(Infinity)
            c1 = max(poss_c1,c1)
    return R(R(0.9999999)/(R(c1)*SUK.rank()))

def c4_func(SUK,v, A, prec=106):
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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
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
    R = RealField(prec)

    return R(max(SUK.number_field().abs_val(v, alpha, prec) for alpha in A))

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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: betas = [[beta,beta.valuation(v_fin)] for beta in SUK.fundamental_units()]
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
        if (pair[1].abs() != 0 and pair[1].abs() < good_pair[1].abs()):
            good_pair = pair
    return good_pair

def mus(SUK,v):
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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: mus(SUK,v_fin)
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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

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
    for rs in combinations_with_replacement(range(nk.abs()),len(betas)):
        # n_0 = valuation_v of one of the coefficients of the equation = 0 for x + y = 1 p. 824
        n_rs = zip(ns,rs)
        sigma_tilde = -(sum([n_r[0]*n_r[1] for n_r in n_rs]))
        if sigma_tilde % nk == 0:
            beta_rs = zip(betas,rs)
            temp_prod = prod([beta_r[0]**beta_r[1] for beta_r in beta_rs])*betak**(sigma_tilde/nk)
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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = K.primes_above(3)[0]
        sage: A = K.roots_of_unity()

        sage: c8_c9_func(SUK,v_fin,A) # abs tol 1e-29
        (4.524941291354698258804956696127e15, 1.621521281297160786545580368612e16)

    REFERENCES:

    - [Sma1995]_ p. 825
    - [Sma1998]_ p. 226, Theorem A.2 for the local constants

    """
    R = RealField(prec)
    num_mus = len(mus(SUK,v))+1
    p = v.smallest_integer()
    f_p = v.residue_class_degree()
    d = SUK.number_field().degree()
    if p == 2:
        local_c2 = Integer(197142)*Integer(36)**num_mus
    elif p%4 == 1:
        local_c2 = Integer(35009)*(Integer(45)/Integer(2))**num_mus
    else:
        local_c2 = Integer(30760)*Integer(25)**num_mus
    x = polygen(SUK.number_field())
    if ( p > 2 and not ((x**2+1).is_irreducible()) ) or ( p==2 and not ((x**2+3).is_irreducible()) ):
        D = d
    else:
        D = 2*d
    l_c3 = (num_mus+1)**(2*num_mus+4)*p**(D * f_p/d)*(f_p*R(p).log())**(-num_mus-1)*D**(num_mus+2)

    def modified_height(SUK,v,D,b):
        #[Sma1998]_ p. 226
        max_log_b = max([phi(b).log().abs() for phi in SUK.number_field().places(prec)])
        return R(max([b.global_height(),max_log_b/(2*R.pi()*D),f_p*R(p).log()/d]))

    mus_prod = prod([modified_height(SUK,v,D,b) for b in mus(SUK,v)])
    local_c3 = R(max([mus_prod*modified_height(SUK,v,D,mu0) for mu0 in possible_mu0s(SUK,v)]))

    l_c3 *= local_c3
    H = max([modified_height(SUK,v,D,alpha) for alpha in mus(SUK,v)+possible_mu0s(SUK,v)])
    if p == 2:
        local_c4 = R(3*2**10*(num_mus+1)**2*D**2*H).log()
    else:
        local_c4 = R(2**11*(num_mus+1)**2*D**2*H).log()
    local_c5 = 2*R(D).log()
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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c11_func(SUK,phi_real,A) # abs tol 1e-29
        3.255848343572896153455615423662

        sage: c11_func(SUK,phi_complex,A) # abs tol 1e-29
        6.511696687145792306911230847323

    REFERENCES:

    - [Sma1995]_ p. 825
    """
    R = RealField(prec)
    if is_real_place(v):
        return R(R(4*c4_func(SUK, v, A, prec)).log()/(c3_func(SUK, prec)))
    else:
        return R(2*(R(4*(c4_func(SUK,v, A, prec)).sqrt()).log())/(c3_func(SUK, prec)))

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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]

        sage: c13_func(SUK,phi_real) # abs tol 1e-29
        0.4257859134798034746197327286726

        sage: c13_func(SUK,phi_complex) # abs tol 1e-29
        0.2128929567399017373098663643363

    It is an error to input a finite place

    ::

        sage: phi_finite = K.primes_above(3)[0]
        sage: c13_func(SUK,phi_finite)
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
        return c3_func(SUK,prec)
    else:
        return c3_func(SUK,prec)/2

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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: A = K.roots_of_unity()

        sage: K0_func(SUK,A) # abs tol 1e-29
        9.475576673109275443280257946929e17

    REFERENCES:

    - [Sma1995]_ p. 824
    """
    R = RealField(prec)

    def c5_func(SUK, v, R):
        return R(c3_func(SUK, R.precision())/(v.residue_class_degree()*R(v.smallest_integer()).log()*v.ramification_index()))

    def c6_func(SUK, v, A, R):
        return R(R(c4_func(SUK, v, A, R.precision())).log()/(v.residue_class_degree()*R(v.smallest_integer()).log()*v.ramification_index()))

    def c7_func(SUK, v, A, R):
        return R(R(c4_func(SUK, v, A, R.precision())).log()/c3_func(SUK, R.precision()))

    def c10_func(SUK, v, A, R):
        # [Sma1995]_ p. 824
        e_h = v.ramification_index()
        c_8, c_9 = c8_c9_func(SUK, v, A, R.precision())
        return R((2/(e_h*c5_func(SUK, v, R)))*(e_h*c6_func(SUK, v, A, R) + c_9 + c_8 * R( c_8/(e_h*c5_func(SUK, v, R))).log()))

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
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: K1_func(SUK,phi_real,A)
        4.396386097852707394927181864635e16

        sage: K1_func(SUK,phi_complex,A)
        2.034870098399844430207420286581e17

    REFERENCES:

    - [Sma1995]_ p. 825

    """
    R = RealField(prec)

    #[Sma1995]_ p. 825
    if is_real_place(v):
        c11 = R(R(4*c4_func(SUK, v, A, prec)).log()/(c3_func(SUK, prec)))
    else:
        c11 = R(2*(R(4*(c4_func(SUK,v, A, prec)).sqrt()).log())/(c3_func(SUK, prec)))

    #[Sma1995]_ p. 825
    if is_real_place(v):
        c12 = R(2*c4_func(SUK, v, A, prec))
    else:
        c12 = R(2*(c4_func(SUK,v, A, prec)).sqrt())

    #[Sma1998]_ p. 225, Theorem A.1
    d = SUK.number_field().degree()
    t = SUK.rank()
    Baker_C = R( 18 * factorial(t+2) * (t+1)**(t+2) * (32*d)**(t + 3) * R( 2*(t+1) * d).log() )

    def hprime(SUK, alpha, v):
        #[Sma1998]_ p. 225
        return R(max(alpha.global_height(), 1/SUK.number_field().degree(), (v(alpha)).log().abs()/SUK.number_field().degree()))

    #[Sma1995]_ p. 825 and [Sma1998]_ p. 225, Theorem A.1
    c14 = Baker_C * prod([hprime(SUK, alpha, v) for alpha in SUK.gens_values()])

    #[Sma1995]_ p. 825
    c15 = R(2*((c12).log()+c14*R((SUK.rank()+1)*c14/c13_func(SUK, v, prec)).log())/c13_func(SUK, v, prec))

    return max([c11, c15])
