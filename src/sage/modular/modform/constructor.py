"""
Creating Spaces of Modular Forms

EXAMPLES:

"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import weakref

import ambient

import sage.modular.congroup as congroup
import sage.modular.dirichlet as dirichlet
import sage.rings.all as rings

import element

def canonical_parameters(group, weight, base_ring):
    weight = rings.Integer(weight)
    if weight <= 1:
        raise ValueError, "the weight must be at least 2"

    if isinstance(group, (int, rings.Integer)):
        group = congroup.Gamma0(group)

    elif isinstance(group, dirichlet.DirichletCharacter):
        group = group.minimize_base_ring()

    if not rings.is_CommutativeRing(base_ring):
        raise TypeError, "base_ring (=%s) must be a commutative ring"%base_ring

    return group, weight, sign, base_ring

_cache = {}

def ModularForms_clear_cache():
    global _cache
    _cache = {}

def ModularForms(group  = 1,
                 weight = 2,
                 base_ring = rings.RationalField(),
                 use_cache = True):
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
        ???

    We create some spaces for $\Gamma_1(N)$.
        sage: ModularForms(Gamma1(13),2)
        ???
        sage: ModularForms(Gamma1(13),2).dimension()
        ???
        sage: ModularForms(Gamma1(13),2).dimension()
        ???
        sage: [ModularForms(Gamma1(7),k).dimension() for k in [2,3,4,5]]
        [???]
        sage: ModularForms(Gamma1(5),11).dimension()
        ???

    We create a space with character:
        sage: e = (DirichletGroup(13).0)^2
        sage: e.order()
        6
        sage: M = ModularForms(e, 2); M
        ???
        sage: f = M.T(2).charpoly(); f
        ???
        sage: f.factor()
        ???

    More examples of spaces with character:
        sage: e = DirichletGroup(5, RationalField()).gen(); e
        [-1]
        sage: m = ModularForms(e, 2); m
        ???
        sage: m.T(2).charpoly()
        x^2 - 1
        sage: m = ModularForms(e, 6); m.dimension()
        6
        sage: m.T(2).charpoly()
        x^6 - 873*x^4 - 82632*x^2 - 1860496
    """
    key = canonical_parameters(group, weight, base_ring)

    if use_cache and _cache.has_key(key):
         M = _cache[key]()
         if not (M is None): return M

    (group, weight, base_ring) = key

    M = None
    if isinstance(group, congroup.Gamma0):
        if base_ring is rings.RationalField():
            M = ambient.ModularFormsAmbient_g0_Q(group.level(), weight)

    elif isinstance(group, congroup.Gamma1):
        if base_ring is rings.RationalField():
            M = ambient.ModularFormsAmbient_g1_Q(group.level(), weight)

    elif isinstance(group, dirichlet.DirichletCharacter):
        eps = group
        M = ambient.ModularFormsAmbient_eps(eps, weight)

    if M is None:
        raise NotImplementedError, \
  "computation of requested space of modular symbols not defined or implemented"

    _cache[key] = weakref.ref(M)
    return M



def CuspForms(group  = 1,
              weight = 2,
              base_ring = rings.RationalField(),
              use_cache = True):
    """
    Create a space of cuspidal modular forms.

    See the documentation for the ModularForms command for a
    description of the input parameters.
    """
    return ModularForms(group, weight, base_ring,
                        use_cache=use_cache).cuspidal_submodule()


def EisensteinForms(group  = 1,
              weight = 2,
              base_ring = rings.RationalField(),
              use_cache = True):
    """
    Create a space of eisenstein modular forms.

    See the documentation for the ModularForms command for a
    description of the input parameters.
    """
    return ModularForms(group, weight, base_ring,
                        use_cache=use_cache).eisenstein_submodule()
