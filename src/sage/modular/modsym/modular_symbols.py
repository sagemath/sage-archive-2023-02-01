"""
Modular symbols
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

import sage.modular.cusps as cusps
import sage.modular.modsym.manin_symbols
import sage.structure.formal_sum as formal_sum
import sage.rings.arith as arith
import sage.rings.rational_field as rational_field
from sage.misc.latex import latex

_C = cusps.Cusps()

class ModularSymbol:
    r"""
    A Modular symbol $X^i\cdot Y^{k-2-i} \{\alpha, \beta\}$.
    """
    def __init__(self, space, i, alpha, beta):
        self.__space = space
        self.__i = i
        self.__alpha = _C(alpha)
        self.__beta = _C(beta)

    def __repr__(self):
        polypart = sage.modular.modsym.manin_symbols._print_polypart(self.__i, self.weight()-2-self.__i)
        if len(polypart) > 0:
            polypart = polypart + "*"
        return "%s{%s,%s}"%(polypart, self.__alpha, self.__beta)

    def _latex_(self):
        polypart = sage.modular.modsym.manin_symbols._print_polypart(self.__i, self.weight()-2-self.__i)
        return "%s\\left\\{%s, %s\\right\\}"%(latex(polypart),
                  latex(self.__alpha), latex(self.__beta))

    def __cmp__(self, other):
        if not isinstance(other, ModularSymbol):
            return -1
        if self.__space == other.__space and self.__i == other.__i \
             and self.__alpha == other.__alpha and self.__beta == other.__beta:
            return 0
        return 1

    def space(self):
        return self.__space

    def i(self):
        return self.__i

    def weight(self):
        return self.__space.weight()

    def alpha(self):
        return self.__alpha

    def beta(self):
        return self.__beta

    def apply(self, g):
        r"""
        INPUT:
            g -- a list [a,b,c,d] where we view [a,b,c,d] as defining a 2x2 matrix in GL_2(Q).
        OUTPUT:
            list -- a list of tuples (coef_i, x_i), where coef_i is a scalar and x_i
                    is a ModularSymbol, such that the sum coef_i*x_i is the image
                    of this symbol under the action of g.  No reduction is performed
                    modulo the relations that hold in self.space().

        The action of $g$ on symbols is by
        $$
           P(X,Y)\{\alpha, \beta\} \mapsto  P(dX-bY, -cx+aY) \{g(\alpha), g(\beta)\}.
        $$

        Note that for us we have $P=X^i Y^{k-2-i}$, which simplifies computation
        of the polynomial part slightly.
        """
        space = self.__space
        i = self.__i
        k = space.weight()
        a,b,c,d = tuple(g)
        coeffs = sage.modular.modsym.manin_symbols.apply_to_monomial(i, k-2, d, -b, -c, a)
        g_alpha = self.__alpha.apply(g)
        g_beta = self.__beta.apply(g)
        return formal_sum.FormalSum([(coeffs[j], ModularSymbol(space, j, g_alpha, g_beta)) \
                                     for j in reversed(range(k-1)) if coeffs[j] != 0])

    def __manin_symbol_rep(self, alpha):
        """
        Return representation of X^i*Y^(k-2-i){0,alpha}.
        """
        space = self.__space
        i = self.__i
        k = space.weight()
        v = [(0,1), (1,0)]
        if not alpha.is_infinity():
            v += [(x.numerator(), x.denominator()) for x in arith.convergents(alpha._rational_())]
        sign = 1
        apply = sage.modular.modsym.manin_symbols.apply_to_monomial
        mansym = sage.modular.modsym.manin_symbols.ManinSymbol
        z = formal_sum.FormalSum(0)
        for j in range(1,len(v)):
            c = sign*v[j][1]
            d = v[j-1][1]
            coeffs = apply(i, k-2, sign*v[j][0], v[j-1][0], sign*v[j][1], v[j-1][1])
            w = [(coeffs[j], mansym(space, (j, c, d))) \
                       for j in range(k-1) if coeffs[j] != 0]
            z += formal_sum.FormalSum(w)
            sign *= -1
        return z

    def manin_symbol_rep(self):
        """
        Returns a representation of self as a formal sum of modular symbols.
        (The result is not cached.)
        """

        alpha = self.__alpha
        beta = self.__beta
        return -1*self.__manin_symbol_rep(alpha) + self.__manin_symbol_rep(beta)

