r"""
Congruence arithmetic subgroups of `{\rm SL}_2(\ZZ)`

Sage can compute extensively with the standard congruence subgroups
`\Gamma_0(N)`, `\Gamma_1(N)`, and `\Gamma_H(N)`.

AUTHORS:

- William Stein
- David Loeffler (2009, 10) -- modifications to work with more general arithmetic subgroups
"""

################################################################################
#
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
################################################################################

from sage.rings.all import QQ, ZZ
from sage.rings.arith import gcd
#import sage.modular.cusps # circular!

from arithgroup_element import ArithmeticSubgroupElement
from arithgroup_generic import ArithmeticSubgroup

def is_CongruenceSubgroup(x):
    """
    Return True if x is of type CongruenceSubgroup.

    EXAMPLES::

        sage: from sage.modular.arithgroup.congroup_generic import is_CongruenceSubgroup
        sage: is_CongruenceSubgroup(SL2Z)
        True
        sage: is_CongruenceSubgroup(Gamma0(13))
        True
        sage: is_CongruenceSubgroup(Gamma1(6))
        True
        sage: is_CongruenceSubgroup(GammaH(11, [3]))
        True
        sage: is_CongruenceSubgroup(SymmetricGroup(3))
        False
    """
    return isinstance(x, CongruenceSubgroup)

class CongruenceSubgroup(ArithmeticSubgroup):
    def __init__(self, level):
        """
        Create a congruence subgroup with given level.

        EXAMPLES::

            sage: Gamma0(500)
            Congruence Subgroup Gamma0(500)
        """
        level = ZZ(level)
        if level <= 0:
            raise ArithmeticError, "Congruence groups only defined for positive levels."
        self.__level = level

    def _repr_(self):
        """
        Return the string representation of self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5)._repr_()
            'Generic congruence subgroup'
        """
        return "Generic congruence subgroup"

    def is_congruence(self):
        r"""
        Return True, since this is a congruence subgroup.

        EXAMPLE::

            sage: Gamma0(7).is_congruence()
            True
        """

        return True

    def modular_symbols(self, sign=0, weight=2, base_ring=QQ):
        """
        Return the space of modular symbols of the specified weight and sign
        on the congruence subgroup self.

        EXAMPLES::

            sage: G = Gamma0(23)
            sage: G.modular_symbols()
            Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: G.modular_symbols(weight=4)
            Modular Symbols space of dimension 12 for Gamma_0(23) of weight 4 with sign 0 over Rational Field
            sage: G.modular_symbols(base_ring=GF(7))
            Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Finite Field of size 7
            sage: G.modular_symbols(sign=1)
            Modular Symbols space of dimension 3 for Gamma_0(23) of weight 2 with sign 1 over Rational Field
        """
        from sage.modular.modsym.modsym import ModularSymbols
        return ModularSymbols(self, sign=sign, weight=weight, base_ring=base_ring)

    def modular_abelian_variety(self):
        """
        Return the modular abelian variety corresponding to the congruence
        subgroup self.

        EXAMPLES::

            sage: Gamma0(11).modular_abelian_variety()
            Abelian variety J0(11) of dimension 1
            sage: Gamma1(11).modular_abelian_variety()
            Abelian variety J1(11) of dimension 1
            sage: GammaH(11,[3]).modular_abelian_variety()
            Abelian variety JH(11,[3]) of dimension 1
        """
        from sage.modular.abvar.abvar_ambient_jacobian import ModAbVar_ambient_jacobian
        return ModAbVar_ambient_jacobian(self)

    def level(self):
        """
        Return the level of this congruence subgroup.

        EXAMPLES::

            sage: SL2Z.level()
            1
            sage: Gamma0(20).level()
            20
            sage: Gamma1(11).level()
            11
            sage: GammaH(14, [5]).level()
            14
        """
        return self.__level

    def _new_group_from_level(self, level):
        r"""
        Return a new group of the same type (Gamma0, Gamma1, or
        GammaH) as self of the given level. In the case that self is
        of type GammaH, we take the largest H inside
        $(Z/\text{level}Z)^\times$ which maps to H, namely its inverse
        image under the natural reduction map.

        EXAMPLES::

            sage: G = Gamma0(20)
            sage: G._new_group_from_level(4)
            Congruence Subgroup Gamma0(4)
            sage: G._new_group_from_level(40)
            Congruence Subgroup Gamma0(40)

            sage: G = Gamma1(10)
            sage: G._new_group_from_level(6)
            Traceback (most recent call last):
            ...
            ValueError: one level must divide the other

            sage: G = GammaH(50,[3,37])
            sage: G
            Congruence Subgroup Gamma_H(50) with H generated by [3, 37]
            sage: G._new_group_from_level(25)
            Congruence Subgroup Gamma_H(25) with H generated by [3, 12]
            sage: G._new_group_from_level(100)
            Congruence Subgroup Gamma_H(100) with H generated by [3, 37, 53, 87]
            sage: G._new_group_from_level(150)
            Congruence Subgroup Gamma_H(150) with H generated by [37, 53, 103, 137]
        """
        from congroup_gamma0 import is_Gamma0
        from congroup_gamma1 import is_Gamma1
        from congroup_gammaH import is_GammaH
        from all import Gamma0, Gamma1, GammaH
        N = self.level()
        if (level%N) and (N%level):
            raise ValueError, "one level must divide the other"
        if is_Gamma0(self):
            return Gamma0(level)
        elif is_Gamma1(self):
            return Gamma1(level)
        elif is_GammaH(self):
            H = self._generators_for_H()
            if level > N:
                d = level // N
                diffs = [ N*i for i in range(d) ]
                newH = [ h + diff for h in H for diff in diffs ]
                return GammaH(level, [x for x in newH if gcd(level, x) == 1])
            else:
                return GammaH(level, [ h%level for h in H ])
        else:
            raise NotImplementedError

    def __call__(self, x, check=True):
        """
        Coerce x into self.

        NOTE: This function should be overridden by any subclass the
        user will interact with directly.

        EXAMPLES::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5).__call__(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if isinstance(x, ArithmeticSubgroupElement) and x.parent() == self:
            return x
        raise NotImplementedError
