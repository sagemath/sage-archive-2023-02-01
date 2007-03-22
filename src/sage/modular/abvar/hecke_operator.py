"""
Hecke operators on modular abelian varieties

\sage can compute with Hecke operators on modular abelian varieties.  A
Hecke operator is defined by given a modular abelian variety and an index.
Given a Hecke operator, \sage can compute the characteristic polynomial,
and the action of the Hecke operator on various homology groups.

TODO: Compute kernel, image, etc. of Hecke operators.

AUTHOR:
    -- William Stein (2007-03)

EXAMPLES:
    sage: A = J0(54)
    sage: t5 = A.hecke_operator(5); t5
    Hecke operator T_5 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(54)
    sage: t5.charpoly().factor()
    (x - 3)^2 * (x + 3)^2 * x^4
    sage: B = A.new_quotient(); B
    Modular abelian quotient variety of dimension 2 and level 54
    sage: t5 = B.hecke_operator(5); t5
    Hecke operator T_5 on Modular abelian quotient variety of dimension 2 and level 54
    sage: t5.charpoly().factor()
    (x - 3)^2 * (x + 3)^2
    sage: t5.action_on_homology().matrix()
    [ 0  3  3 -3]
    [-3  3  3  0]
    [ 3  3  0 -3]
    [-3  6  3 -3]
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.rings.all import ZZ
from morphism import Morphism

class HeckeOperator(Morphism):
    def __init__(self, abvar, n):
        n = ZZ(n)
        if n <= 0:
            raise ValueError, "n must be positive"
        self.__abvar = abvar
        self._n = n

    def _repr_(self):
        return "Hecke operator T_%s on %s"%(self._n, self.__abvar)

    def index(self):
        """
        Return the index of this Hecke operator.  I.e., if this is
        the operator $T_n$, then the index is the integer $n$.

        EXAMPLES:
            sage: J = J1(12345)
            sage: t = J.hecke_operator(997)
            sage: t
            Hecke operator T_997 on Jacobian of the modular curve associated to the congruence subgroup Gamma1(12345)
            sage: t.index()
            997
        """
        return self._n

    def characteristic_polynomial(self, var='x'):
        """
        Return the characteristic polynomial of this Hecke operator
        in the given variable.

        INPUT:
            var -- a string (default: 'x')
        OUTPUT:
            a polynomial in x

        EXAMPLES:
            sage: A = J0(43)[1]; A
            Modular abelian quotient variety of dimension 2 and level 43
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Modular abelian quotient variety of dimension 2 and level 43
            sage: f = t2.characteristic_polynomial(); f
            x^4 - 4*x^2 + 4
            sage: f.factor()
            (x^2 - 2)^2
            sage: t2.characteristic_polynomial('y')
            y^4 - 4*y^2 + 4
        """
        return self.__abvar.rational_homology().hecke_polynomial(self._n, var)

    def charpoly(self, var='x'):
        """
        Synonym for \code{self.characteristic_polynomial(var)}.

        EXAMPLES:
            sage: A = J1(13)
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Jacobian of the modular curve associated to the congruence subgroup Gamma1(13)
            sage: f = t2.charpoly(); f
            x^4 + 6*x^3 + 15*x^2 + 18*x + 9
            sage: f.factor()
            (x^2 + 3*x + 3)^2
            sage: t2.charpoly('y')
            y^4 + 6*y^3 + 15*y^2 + 18*y + 9
        """
        return self.characteristic_polynomial(var)

    def action_on_homology(self, R=ZZ):
        """
        Return the action of this Hecke operator on the homology
        $H_1(A; R)$ of this abelian variety with coefficients in $R$.

        EXAMPLES:
            sage: A = J0(43)
            sage: t2 = A.hecke_operator(2); t2
            Hecke operator T_2 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
            sage: h2 = t2.action_on_homology(); h2
            Hecke operator T_2 on Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
            sage: h2.matrix()
            [-2  1  0  0  0  0]
            [-1  1  1  0 -1  0]
            [-1  0 -1  2 -1  1]
            [-1  0  1  1 -1  1]
            [ 0 -2  0  2 -2  1]
            [ 0 -1  0  1  0 -1]
            sage: h2 = t2.action_on_homology(GF(2)); h2
            Hecke operator T_2 on Homology with coefficients in Finite Field of size 2 of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
            sage: h2.matrix()
            [0 1 0 0 0 0]
            [1 1 1 0 1 0]
            [1 0 1 0 1 1]
            [1 0 1 1 1 1]
            [0 0 0 0 0 1]
            [0 1 0 1 0 1]
        """
        return self.__abvar.homology(R).hecke_operator(self.index())
