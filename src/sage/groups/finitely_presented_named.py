r"""
Named Finitely Presented Groups

Construct groups of small order and "named" groups as quotients of free groups. These
groups are available through tab completion by typing ``groups.presentation.<tab>``
or by importing the required methods. Tab completion is made available through
Sage's :ref:`group catalog <sage.groups.groups_catalog>`. Some examples are engineered
from entries in [THOMAS-WOODS]_.

Groups available as finite presentations:

- Cyclic group, `C_n` of order `n` --
  :func:`groups.presentation.Cyclic <sage.groups.finitely_presented_named.CyclicPresentation>`

- Dicyclic group, nonabelian groups of order `4n` with a unique element of
  order 2 --
  :func:`groups.presentation.DiCyclic <sage.groups.finitely_presented_named.DiCyclicPresentation>`

- Dihedral group, `D_n` of order `2n` --
  :func:`groups.presentation.Dihedral <sage.groups.finitely_presented_named.DihedralPresentation>`

- Klein four group, `C_2 \times C_2` --
  :func:`groups.presentation.KleinFour <sage.groups.finitely_presented_named.KleinFourPresentation>`

AUTHORS:

- Davis Shurbert (2013-06-21): initial version

EXAMPLES::

    sage: groups.presentation.Cyclic(4)
    Finitely presented group < a | a^4 >

You can also import the desired functions::

    sage: from sage.groups.finitely_presented_named import CyclicPresentation
    sage: CyclicPresentation(4)
    Finitely presented group < a | a^4 >
"""
#*****************************************************************************
#       Copyright (C) 2013 Davis Shurbert <davis.sprout@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all      import Integer
from sage.groups.free_group import FreeGroup
from sage.groups.finitely_presented import FinitelyPresentedGroup

def CyclicPresentation(n):
    r"""
    Build cyclic group of order `n` as a finitely presented group.

    INPUT:

    - ``n`` -- The order of the cyclic presentation to be returned.

    OUTPUT:

    The cyclic group of order `n` as finite presentation.

    EXAMPLES::

        sage: groups.presentation.Cyclic(10)
        Finitely presented group < a | a^10 >
        sage: n = 8; C = groups.presentation.Cyclic(n)
        sage: C.as_permutation_group().is_isomorphic(CyclicPermutationGroup(n))
        True

    TESTS::

        sage: groups.presentation.Cyclic(0)
        Traceback (most recent call last):
        ...
        ValueError: finitely presented group order must be positive
    """
    n = Integer(n)
    if n < 1:
        raise ValueError('finitely presented group order must be positive')
    F = FreeGroup( 'a' )
    rls = F([1])**n,
    return FinitelyPresentedGroup( F, rls )

def DihedralPresentation(n):
    r"""
    Build the Dihedral group of order `2n` as a finitely presented group.

    INPUT:

    - ``n`` -- The size of the set that `D_n` is acting on.

    OUTPUT:

    Dihedral group of order `2n`.

    EXAMPLES::

        sage: D = groups.presentation.Dihedral(7); D
        Finitely presented group < a, b | a^7, b^2, (a*b)^2 >
        sage: D.as_permutation_group().is_isomorphic(DihedralGroup(7))
        True

    TESTS::

        sage: n = 9
        sage: D = groups.presentation.Dihedral(n)
        sage: D.ngens() == 2
        True
        sage: groups.presentation.Dihedral(0)
        Traceback (most recent call last):
        ...
        ValueError: finitely presented group order must be positive
    """
    n = Integer( n )
    if n < 1:
        raise ValueError('finitely presented group order must be positive')
    F = FreeGroup([ 'a', 'b' ])
    rls = F([1])**n, F([2])**2, (F([1])*F([2]))**2
    return FinitelyPresentedGroup( F, rls )

def DiCyclicPresentation(n):
    r"""
    Build the dicyclic group of order `4n`, for `n \geq 2`, as a finitely
    presented group.

    INPUT:

    - ``n`` -- positive integer, 2 or greater, determining the order of
      the group (`4n`).

    OUTPUT:

    The dicyclic group of order `4n` is defined by the presentation

    .. MATH::

        \langle a, x \mid a^{2n}=1, x^{2}=a^{n}, x^{-1}ax=a^{-1} \rangle

    .. NOTE::

        This group is also available as a permutation group via
        :class:`groups.permutation.DiCyclic <sage.groups.perm_gps.permgroup_named.DiCyclicGroup>`.

    EXAMPLES::

        sage: D = groups.presentation.DiCyclic(9); D
        Finitely presented group < a, b | a^18, b^2*a^-9, b^-1*a*b*a >
        sage: D.as_permutation_group().is_isomorphic(groups.permutation.DiCyclic(9))
        True

    TESTS::

        sage: Q = groups.presentation.DiCyclic(2)
        sage: Q.as_permutation_group().is_isomorphic(QuaternionGroup())
        True
        sage: all([groups.presentation.DiCyclic(i).as_permutation_group(
        ....: ).is_isomorphic(groups.permutation.DiCyclic(i)) for i in [5,8,12,2^5]])
        True
        sage: groups.presentation.DiCyclic(1)
        Traceback (most recent call last):
        ...
        ValueError: input integer must be greater than 1
    """
    n = Integer(n)
    if n < 2:
        raise ValueError('input integer must be greater than 1')

    F = FreeGroup(['a','b'])
    rls =  F([1])**(2*n), F([2,2])*F([-1])**n, F([-2,1,2,1])
    return FinitelyPresentedGroup(F, rls)

def KleinFourPresentation():
    r"""
    Build the Klein group of order `4` as a finitely presented group.

    OUTPUT:

    Klein four group (`C_2 \times C_2`) as a finitely presented group.

    EXAMPLES::

        sage: K = groups.presentation.KleinFour(); K
        Finitely presented group < a, b | a^2, b^2, a^-1*b^-1*a*b >
    """
    F = FreeGroup(['a','b'])
    rls = F([1])**2, F([2])**2, F([-1])*F([-2])*F([1])*F([2])
    return FinitelyPresentedGroup(F, rls)
