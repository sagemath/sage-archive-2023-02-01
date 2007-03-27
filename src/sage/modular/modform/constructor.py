"""
Creating Spaces of Modular Forms

EXAMPLES:

"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import weakref


import sage.modular.congroup as congroup
import sage.modular.dirichlet as dirichlet
import sage.rings.all as rings

import ambient
import ambient_eps
import ambient_g0
import ambient_g1
import ambient_R
import defaults
import element


def canonical_parameters(group, weight, base_ring):
    weight = rings.Integer(weight)
    if weight <= 1:
        raise NotImplementedError, "weight must be at least 2"

    if isinstance(group, (int, rings.Integer)):
        group = congroup.Gamma0(group)

    elif isinstance(group, dirichlet.DirichletCharacter):
        group = group.minimize_base_ring()

    if isinstance(group, congroup.SL2Z) or \
       isinstance(group, congroup.Gamma1) and group.level() == 1:
        group = congroup.Gamma0(1)

    if not rings.is_CommutativeRing(base_ring):
        raise TypeError, "base_ring (=%s) must be a commutative ring"%base_ring

    # it is *very* important to include the level as part of the data
    # that defines the key, since dirichlet characters of different
    # levels can compare equal, but define much different modular
    # forms spaces.
    return group.level(), group, weight, base_ring

_cache = {}

def ModularForms_clear_cache():
    global _cache
    _cache = {}

def ModularForms(group  = 1,
                 weight = 2,
                 base_ring = None,
                 use_cache = True,
                 prec = defaults.DEFAULT_PRECISION):
    r"""
    Create an ambient space of modular forms.

    INPUT:
        group -- A congruence subgroup or a Dirichlet character eps.
        weight -- int, the weight, which must be an integer >= 1.
        base_ring -- the base ring (ignored if group is a Dirichlet character)

    Create using the command
        ModularForms(group, weight, character)


    EXAMPLES:
    First we create some spaces with trivial character:
        sage: ModularForms(Gamma0(11),2).dimension()
        2
        sage: ModularForms(Gamma0(1),12).dimension()
        2

    If we give an integer N for the congruence subgroup, it defaults
    to $\Gamma_0(N)$:
        sage: ModularForms(1,12).dimension()
        2
        sage: ModularForms(11,4)
        Modular Forms space of dimension 4 for Congruence Subgroup Gamma0(11) of weight 4 over Rational Field

    We create some spaces for $\Gamma_1(N)$.
        sage: ModularForms(Gamma1(13),2)
        Modular Forms space of dimension 13 for Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
        sage: ModularForms(Gamma1(13),2).dimension()
        13
        sage: [ModularForms(Gamma1(7),k).dimension() for k in [2,3,4,5]]
        [5, 7, 9, 11]
        sage: ModularForms(Gamma1(5),11).dimension()
        12

    We create a space with character:
        sage: e = (DirichletGroup(13).0)^2
        sage: e.order()
        6
        sage: M = ModularForms(e, 2); M
        Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2
        sage: f = M.T(2).charpoly('x'); f
        x^3 + (-2*zeta6 - 2)*x^2 + (-2*zeta6)*x + 14*zeta6 - 7
        sage: f.factor()
        (x + -2*zeta6 - 1) * (x + -zeta6 - 2) * (x + zeta6 + 1)

    More examples of spaces with character:
        sage: e = DirichletGroup(5, RationalField()).gen(); e
        [-1]
        sage: m = ModularForms(e, 2); m
        Modular Forms space of dimension 2, character [-1] and weight 2 over Rational Field
        sage: m.T(2).charpoly('x')
        x^2 - 1
        sage: m = ModularForms(e, 6); m.dimension()
        4
        sage: m.T(2).charpoly('x')
        x^4 - 917*x^2 - 42284
    """
    if isinstance(group, dirichlet.DirichletCharacter):
        if base_ring is None:
            base_ring = group.minimize_base_ring().base_ring()
    if base_ring is None:
        base_ring = rings.QQ

    key = canonical_parameters(group, weight, base_ring)

    if use_cache and _cache.has_key(key):
         M = _cache[key]()
         if not (M is None):
             M.set_precision(prec)
             return M

    (level, group, weight, base_ring) = key

    M = None
    if isinstance(group, congroup.Gamma0):
        M = ambient_g0.ModularFormsAmbient_g0_Q(group.level(), weight)
        if base_ring != rings.QQ:
            M = ambient_R.ModularFormsAmbient_R(M, base_ring)

    elif isinstance(group, congroup.Gamma1):
        M = ambient_g1.ModularFormsAmbient_g1_Q(group.level(), weight)
        if base_ring != rings.QQ:
            M = ambient_R.ModularFormsAmbient_R(M, base_ring)

    elif isinstance(group, dirichlet.DirichletCharacter):
        eps = group
        eps = eps.minimize_base_ring()
        if eps.is_trivial():
            return ModularForms(eps.modulus(), weight, base_ring,
                                use_cache = use_cache,
                                prec = prec)
        M = ambient_eps.ModularFormsAmbient_eps(eps, weight)
        if base_ring != eps.base_ring():
            M = ambient_R.ModularFormsAmbient_R(M, base_ring)

    if M is None:
        raise NotImplementedError, \
           "computation of requested space of modular forms not defined or implemented"

    M.set_precision(prec)
    _cache[key] = weakref.ref(M)
    return M


def CuspForms(group  = 1,
              weight = 2,
              base_ring = None,
              use_cache = True,
              prec = defaults.DEFAULT_PRECISION):
    """
    Create a space of cuspidal modular forms.

    See the documentation for the ModularForms command for a
    description of the input parameters.
    """
    return ModularForms(group, weight, base_ring,
                        use_cache=use_cache, prec=prec).cuspidal_submodule()


def EisensteinForms(group  = 1,
              weight = 2,
              base_ring = None,
              use_cache = True,
              prec = defaults.DEFAULT_PRECISION):
    """
    Create a space of eisenstein modular forms.

    See the documentation for the ModularForms command for a
    description of the input parameters.
    """
    return ModularForms(group, weight, base_ring,
                        use_cache=use_cache, prec=prec).eisenstein_submodule()
