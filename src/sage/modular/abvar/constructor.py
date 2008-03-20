"""
Constructors for certain modular abelian varieties.

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

import weakref

from sage.rings.integer import Integer

from sage.modular.congroup import is_CongruenceSubgroup
from sage.modular.modsym.space import is_ModularSymbolsSpace
from sage.misc.misc import prod
from abvar_newform import ModularAbelianVariety_newform
import sage.modular.modform.element

_cache = {}

def _get(key):
    if _cache.has_key(key):
        z = _cache[key]()
        if z is not None:
            return z
    raise ValueError

def _saved(key, J):
    _cache[key] = weakref.ref(J)
    return J


def J0(N):
    """
    Return the Jacobian $J_0(N)$ of the modular curve $X_0(N)$.

    EXAMPLES:
        sage: J0(389)
        Abelian variety J0(389) of dimension 32

    The result is cached:
        sage: J0(33) is J0(33)
        True
    """
    key = 'J0(%s)'%N
    try:
        return _get(key)
    except ValueError:
        from sage.modular.congroup import Gamma0
        J = Gamma0(N).modular_abelian_variety()
        return _saved(key, J)

def J1(N):
    """
    Return the Jacobian $J_1(N)$ of the modular curve $X_1(N)$.

    EXAMPLES:
        sage: J1(389)
        Abelian variety J1(389) of dimension 6112
    """
    key = 'J1(%s)'%N
    try:
        return _get(key)
    except ValueError:
        from sage.modular.congroup import Gamma1
        return _saved(key, Gamma1(N).modular_abelian_variety())

def JH(N, H):
    """
    Return the Jacobian $J_H(N)$ of the modular curve $X_H(N)$.

    EXAMPLES:
        sage: JH(389,[2])
        Abelian variety JH(389,[2]) of dimension 32
    """
    key = 'JH(%s,%s)'%(N,H)
    try:
        return _get(key)
    except ValueError:
        from sage.modular.congroup import GammaH
        return _saved(key, GammaH(N, H).modular_abelian_variety())

def AbelianVariety(groups=None, lattice=None, modsym=None, base_field=None):
    """
    Create the abelian variety corresponding to the given definining data.
    """

    if is_CongruenceSubgroup(groups):
        groups = [groups]

    if isinstance(groups, str):
        from sage.modular.modform.constructor import Newform
        f = Newform(groups, names='a')
        return ModularAbelianVariety_newform(f, internal_name=True)

    elif isinstance(f, sage.modular.modform.element.Newform):
        return ModularAbelianVariety_newform(f)

    if groups is not None and all([is_CongruenceSubgroup(G) for G in groups]):
        return prod([G.modular_abelian_variety() for G in groups])

    if modsym is not None and is_ModularSymbolsSpace(modsym):
        return modsym.modular_abelian_variety()

    raise NotImplementedError, "arguments to AbelianVariety not recognized"
