"""
Elementary symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
import sfa, multiplicative, classical

###################################
#                                 #
# Elementary Symmetric Functions  #
#                                 #
###################################
class SymmetricFunctionAlgebra_elementary(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, R):
        """
        TESTS::

            sage: e = SFAElementary(QQ)
            sage: e == loads(dumps(e))
            True
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, R, "elementary", 'e')


    class Element(classical.SymmetricFunctionAlgebra_classical.Element):
        # TODO: fix indentation
      def omega(self):
        """
        Returns the image of self under the Frobenius / omega
        automorphism.

        EXAMPLES::

            sage: e = SFAElementary(QQ)
            sage: a = e([2,1]); a
            e[2, 1]
            sage: a.omega()
            e[1, 1, 1] - e[2, 1]

        ::

            sage: h = SFAHomogeneous(QQ)
            sage: h(e([2,1]).omega())
            h[2, 1]
        """
        e = self.parent()
        h = sfa.SFAHomogeneous(e.base_ring())
        return e( h._from_element(self) )


      def expand(self, n, alphabet='x'):
        """
        Expands the symmetric function as a symmetric polynomial in `n`
        variables.

        EXAMPLES::

            sage: e = SFAElementary(QQ)
            sage: e([2,1]).expand(3)
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 3*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
            sage: e([1,1,1]).expand(2)
            x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + x1^3
            sage: e([3]).expand(2)
            0
            sage: e([2]).expand(3)
            x0*x1 + x0*x2 + x1*x2
            sage: e([3]).expand(4,alphabet='x,y,z,t')
            x*y*z + x*y*t + x*z*t + y*z*t
            sage: e([3]).expand(4,alphabet='y')
            y0*y1*y2 + y0*y1*y3 + y0*y2*y3 + y1*y2*y3
        """
        condition = lambda part: max(part) > n
        return self._expand(condition, n, alphabet)

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.elementary', 'SymmetricFunctionAlgebraElement_elementary',  SymmetricFunctionAlgebra_elementary.Element)
