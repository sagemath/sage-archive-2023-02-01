"""
Elementary symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#                     2012 Anne Schilling <anne@math.ucdavis.edu>
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
import sage
from sage.combinat.partition import Partition


###################################
#                                 #
# Elementary Symmetric Functions  #
#                                 #
###################################
class SymmetricFunctionAlgebra_elementary(multiplicative.SymmetricFunctionAlgebra_multiplicative):
    def __init__(self, Sym):
        """
        A class for methods for the elementary basis of the symmetric functions.

        INPUT:

        - ``self`` -- an elementary basis of the symmetric functions
        - ``Sym`` -- an instance of the ring of symmetric functions

        TESTS::

            sage: e = SymmetricFunctions(QQ).e()
            sage: e == loads(dumps(e))
            True
            sage: TestSuite(e).run(skip=['_test_associativity', '_test_distributivity', '_test_prod'])
            sage: TestSuite(e).run(elements = [e[1,1]+e[2], e[1]+2*e[1,1]])
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, Sym, "elementary", 'e')

    def _dual_basis_default(self):
        """
        Returns the default value for ``self.dual_basis()``

        This method returns the dual basis to the elementary basis
        with respect to the standard scalar product, that is the
        forgotten basis.

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).e()
            sage: e.dual_basis()
            Symmetric Functions over Rational Field in the forgotten basis

        TESTS::

            sage: e._dual_basis_default() is e.dual_basis()
            True
        """
        return self.dual_basis(scalar = None, prefix="f", basis_name = "forgotten")

    def coproduct_on_generators(self, i):
        r"""
        Returns the coproduct on ``self[i]``.

        INPUT:

        - ``self`` -- an elementary basis of the symmetric functions
        - ``i`` -- a nonnegative integer

        OUTPUT:

        - returns the coproduct on the elementary generator `e(i)`

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: e = Sym.elementary()
            sage: e.coproduct_on_generators(2)
            e[] # e[2] + e[1] # e[1] + e[2] # e[]
            sage: e.coproduct_on_generators(0)
            e[] # e[]
        """
        def P(i): return Partition([i]) if i else Partition([])
        T = self.tensor_square()
        return T.sum_of_monomials( (P(j), P(i-j)) for j in range(i+1) )

    class Element(classical.SymmetricFunctionAlgebra_classical.Element):
        def omega(self):
            """
            Returns the image of self under the Frobenius / omega automorphism.

            INPUT:

            - ``self`` -- an element of the elementary basis of the symmetric functions

            OUTPUT:

            - returns the image of the element under the omega automorphism

            EXAMPLES::

                sage: e = SymmetricFunctions(QQ).e()
                sage: a = e([2,1]); a
                e[2, 1]
                sage: a.omega()
                e[1, 1, 1] - e[2, 1]

            ::

                sage: h = SymmetricFunctions(QQ).h()
                sage: h(e([2,1]).omega())
                h[2, 1]
            """
            e = self.parent()
            h = sage.combinat.sf.sf.SymmetricFunctions(e.base_ring()).h()
            return e( h._from_element(self) )


        def expand(self, n, alphabet='x'):
            """
            Expands the symmetric function as a symmetric polynomial in `n` variables.

            INPUT:

            - ``self`` -- an element of the elementary basis of the symmetric functions
            - ``n`` -- a positive integer
            - ``alphabet`` -- a variable for the expansion (default: `x`)

            OUTPUT:

            - a monomial expansion of an instance of ``self`` in `n` variables

            EXAMPLES::

                sage: e = SymmetricFunctions(QQ).e()
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
                sage: e([]).expand(2)
                1
            """
            condition = lambda part: max(part) > n
            return self._expand(condition, n, alphabet)

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.elementary', 'SymmetricFunctionAlgebraElement_elementary',  SymmetricFunctionAlgebra_elementary.Element)
