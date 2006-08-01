"""
Creation of modular symbols spaces

EXAMPLES:
    sage: C = HeckeModules(RationalField()); C
    Category of Hecke modules over Rational Field
    sage: M = ModularSymbols(11)
    sage: M.category()
    Category of Hecke modules over Rational Field
    sage: M in C
    True
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import weakref

import ambient
import sage.modular.congroup as congroup
import sage.modular.dirichlet as dirichlet
import sage.rings.rational_field as rational_field
import sage.rings.all as rings

import element

def canonical_parameters(group, weight, sign, base_ring):
    sign = rings.Integer(sign)
    if not (sign in [-1,0,1]):
        raise ValueError, "sign must be -1, 0, or 1"

    weight = rings.Integer(weight)
    if weight <= 1:
        raise ValueError, "the weight must be at least 2"

    if isinstance(group, (int, rings.Integer)):
        group = congroup.Gamma0(group)

    elif isinstance(group, dirichlet.DirichletCharacter):
        try:
            group = group.minimize_base_ring()
        except NotImplementedError:  # todo -- implement minimize_base_ring over finite fields
            pass

    if not rings.is_CommutativeRing(base_ring):
        raise TypeError, "base_ring (=%s) must be a commutative ring"%base_ring

    return group, weight, sign, base_ring

_cache = {}

def ModularSymbols_clear_cache():
    global _cache
    _cache = {}

def ModularSymbols(group  = 1,
                   weight = 2,
                   sign   = 0,
                   base_ring = rational_field.RationalField(),
                   use_cache = True):
    r"""
    Create an ambient space of modular symbols.

    INPUT:
        group -- A congruence subgroup or a Dirichlet
                 character eps.

        weight -- int, the weight, which must be >= 2.

        sign -- int, The sign of the involution on modular symbols
                induced by complex conjugation.  The default is 0,
                which means"no sign", i.e., take the
                whole space.

        base_ring -- the base ring.  This is ignored if group
                is a Dirichlet character.

    EXAMPLES:
    First we create some spaces with trivial character:
        sage: ModularSymbols(Gamma0(11),2).dimension()
        3
        sage: ModularSymbols(Gamma0(1),12).dimension()
        3

    If we give an integer N for the congruence subgroup, it defaults
    to $\Gamma_0(N)$:
        sage: ModularSymbols(1,12,-1).dimension()
        1
        sage: ModularSymbols(11,4, sign=1)
        Modular Symbols space of dimension 4 for Gamma_0(11) of weight 4 with sign 1 over Rational Field

    We create some spaces for $\Gamma_1(N)$.
        sage: ModularSymbols(Gamma1(13),2)
        Modular Symbols space of dimension 15 for Gamma_1(13) of weight 2 with sign 0 and over Rational Field
        sage: ModularSymbols(Gamma1(13),2, sign=1).dimension()
        13
        sage: ModularSymbols(Gamma1(13),2, sign=-1).dimension()
        2
        sage: [ModularSymbols(Gamma1(7),k).dimension() for k in [2,3,4,5]]
        [5, 8, 12, 16]
        sage: ModularSymbols(Gamma1(5),11).dimension()
        20

    We create a space with character:
        sage: e = (DirichletGroup(13).0)^2
        sage: e.order()
        6
        sage: M = ModularSymbols(e, 2); M
        Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
        sage: f = M.T(2).charpoly(); f
        x^4 + (-zeta6 - 1)*x^3 + (-8*zeta6)*x^2 + (10*zeta6 - 5)*x + 21*zeta6 - 21
        sage: f.factor()
        (x + -2*zeta6 - 1) * (x + -zeta6 - 2) * (x + zeta6 + 1)^2

    More examples of spaces with character:
        sage: e = DirichletGroup(5, RationalField()).gen(); e
        [-1]
        sage: m = ModularSymbols(e, 2); m
        Modular Symbols space of dimension 2 and level 5, weight 2, character [-1], sign 0, over Rational Field

        sage: m.T(2).charpoly()
        x^2 - 1
        sage: m = ModularSymbols(e, 6); m.dimension()
        6
        sage: m.T(2).charpoly()
        x^6 - 873*x^4 - 82632*x^2 - 1860496

    We create a space of modular symbols with nontrivial character in characteristic 2.
        sage: G = DirichletGroup(13,GF(4)); G
        Group of Dirichlet characters of modulus 13 over Finite Field in a of size 2^2
        sage: e = G.list()[2]; e
        [a + 1]
        sage: M = ModularSymbols(e,4); M
        Modular Symbols space of dimension 8 and level 13, weight 4, character [a + 1], sign 0, over Finite Field in a of size 2^2
        sage: M.basis()
        ([X*Y,(1,0)], [X*Y,(1,5)], [X*Y,(1,10)], [X*Y,(1,11)], [X^2,(0,1)], [X^2,(1,10)], [X^2,(1,11)], [X^2,(1,12)])
        sage: M.T(2).matrix()
        [    0     0     0     0     0     0     1     1]
        [    0     0     0     0     0     0     0     0]
        [    0     0     0     0     0 a + 1     1     a]
        [    0     0     0     0     0     1 a + 1     a]
        [    0     0     0     0 a + 1     0     1     1]
        [    0     0     0     0     0     a     1     a]
        [    0     0     0     0     0     0 a + 1     a]
        [    0     0     0     0     0     0     1     0]
    """
    key = canonical_parameters(group, weight, sign, base_ring)

    if use_cache and _cache.has_key(key):
         M = _cache[key]()
         if not (M is None): return M

    (group, weight, sign, base_ring) = key

    M = None
    if isinstance(group, congroup.Gamma0):
            if weight == 2:
                M = ambient.ModularSymbolsAmbient_wt2_g0(
                    group.level(),sign, base_ring)
            else:
                M = ambient.ModularSymbolsAmbient_wtk_g0(
                    group.level(), weight, sign, base_ring)

    elif isinstance(group, congroup.Gamma1):

        M = ambient.ModularSymbolsAmbient_wtk_g1(group.level(), weight, sign, base_ring)

    elif isinstance(group, dirichlet.DirichletCharacter):

        eps = group
        M = ambient.ModularSymbolsAmbient_wtk_eps(eps, weight, sign)

    if M is None:
        raise NotImplementedError, "computation of requested space of modular symbols not defined or implemented"

    _cache[key] = weakref.ref(M)
    return M

