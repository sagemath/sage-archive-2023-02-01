"""
Optimal abelian variety quotients of modular jacobians.

AUTHOR:
    -- William Stein (2007-03)

EXAMPLES:
    sage: A = J0(389)[0]; A
    Abelian variety factor of dimension 1 of J0(389)
    sage: loads(dumps(A)) == A
    True
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from abvar             import ModularAbelianVariety_modsym, is_ModularAbelianVariety
from sage.rings.all    import QQ

class ModAbVar_modsym_factor(ModularAbelianVariety_modsym):
    r"""
    A factor of a modular abelian variety.

    A \emph{factor} is an abelian subvariety (attached to a modular
    symbols space) of an ambient modular abelian variety $J$ attached
    to a congruence subgroup.
    """
    def __init__(self, ambient, modsym):
        """
        Create a modular abelian variety factor attached to a modular
        symbols space.

        EXAMPLES:
            sage: A = J0(389)[0]; A
            Abelian variety quotient of dimension 1 and level 389
            sage: type(A)
            <class 'sage.modular.abvar.abvar_modsym_factor.ModAbVar_modsym_factor'>
        """
        self._ambient = ambient
        self._modsym = {1:modsym}
        assert modsym.base_ring() == QQ
        assert modsym.sign() == 1
        ModularAbelianVariety_modsym.__init__(self, level = modsym.level(), base_field = QQ)

    def ambient_variety(self):
        """
        Return the ambient variety of which this modular abelian
        variety is a *quotient*.

        EXAMPLES:
            sage: A = J0(33)[1]; A
            Modular abelian variety quotient of dimension 2 and level 33
            sage: A.ambient_variety()
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(33)
        """
        return self._ambient

    def factor_number(self):
        """
        Return the factor number of this modular abelian variety.
        This is the position of this factor in the sorted list
        of factors of the ambient variety.

        EXAMPLES:
            sage: J = J0(37)

        Set a and b to be the two factors of $J_0(37)$:
            sage: a,b = J

        Now get their factor numbers.
            sage: a.factor_number()
            0
            sage: b.factor_number()
            1
        """
        return self.modular_symbols(sign=1).factor_number()

    def modular_symbols(self, sign=0):
        """
        Return space of modular symbols (with given sign) associated
        to this modular abelian variety.

        INPUT:
            sign -- integer, either -1, 0 or 1 (default: 0)

        EXAMPLES:
        """
        sign = int(sign)
        try:
            return self._modsym[sign]
        except KeyError:
            pass
        M = self._modsym[1]
        A = M.modular_symbols_of_sign(sign)
        self._modsym[sign] = A
        return A

    def _repr_(self):
        """
        Return string representation of this factor.

        EXAMPLES:
            sage: A = J0(389)[0]
            sage: A._repr_()
            'Modular abelian variety quotient of dimension 1 and level 389'
        """
        return "Modular abelian variety quotient of dimension %s and level %s"%(\
            self.dimension(), self.level())
