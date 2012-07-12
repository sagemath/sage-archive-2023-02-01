"""
Monomial symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2010 Anne Schilling <anne at math.ucdavis.edu> (addition)
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
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

import classical, sfa, dual
import sage.libs.symmetrica.all as symmetrica
from sage.rings.integer import Integer
from sage.combinat.partition import Partition
import sage

class SymmetricFunctionAlgebra_monomial(classical.SymmetricFunctionAlgebra_classical):
    def __init__(self, Sym):
        """
        A class for methods related to monomial symmetric functions

        INPUT:

        - ``self`` -- a monomial symmetric function basis
        - ``Sym`` -- an instance of the ring of the symmetric functions

        TESTS::

            sage: m = SymmetricFunctions(QQ).m()
            sage: m == loads(dumps(m))
            True
            sage: TestSuite(m).run(skip=['_test_associativity', '_test_distributivity', '_test_prod'])
            sage: TestSuite(m).run(elements = [m[1,1]+m[2], m[1]+2*m[1,1]])
        """
        classical.SymmetricFunctionAlgebra_classical.__init__(self, Sym, "monomial", 'm')

    def dual_basis(self, scalar = None, scalar_name="",  prefix=None):
        """
        Returns the dual basis to ``self``.

        INPUT:

        - ``self`` -- a monomial symmetric function basis
        - ``scalar`` -- optional input which specifies a function ``zee`` on partitions. The function
          ``zee`` determines the scalar product on the power sum basis
          with normalization `<p_\mu, p_\mu> = zee(mu)`.
          (default: uses standard ``zee`` function)
        - ``scalar_name`` -- specifies the name of the scalar function (optional)
        - ``prefix`` -- optional input, specifies the prefix to be used to display the basis.

        OUTPUT:

        - This method returns the dual basis of the monomial basis with respect to the
          standard scalar product, which is the homogeneous basis. If ``zee`` is specified,
          the modified scalar product is used.

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: h = SymmetricFunctions(QQ).h()
            sage: m.dual_basis() == h
            True

            sage: zee = lambda x : x.centralizer_size()
            sage: dm = m.dual_basis(zee)
            sage: dm[3,1].scalar(m[2,1,1])
            0
            sage: m[2,1,1].scalar(dm[3,1])
            0
        """
        if scalar is None:
            return sage.combinat.sf.sf.SymmetricFunctions(self.base_ring()).h()
        else:
            return dual.SymmetricFunctionAlgebra_dual(self, scalar, scalar_name, prefix)

    def _multiply(self, left, right):
        """
        Returns the product of ``left`` and ``right``.

        - ``self`` -- a monomial symmetric function basis
        - ``left``, ``right`` -- instances of the Schur basis ``self``.

        OUPUT:

        - an element of the Schur basis, the product of ``left`` and ``right``

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: a = m([2,1])
            sage: a^2
            4*m[2, 2, 1, 1] + 6*m[2, 2, 2] + 2*m[3, 2, 1] + 2*m[3, 3] + 2*m[4, 1, 1] + m[4, 2]

        ::

            sage: QQx.<x> = QQ['x']
            sage: m = SymmetricFunctions(QQx).m()
            sage: a = m([2,1])+x
            sage: 2*a # indirect doctest
            2*x*m[] + 2*m[2, 1]
            sage: a^2
            x^2*m[] + 2*x*m[2, 1] + 4*m[2, 2, 1, 1] + 6*m[2, 2, 2] + 2*m[3, 2, 1] + 2*m[3, 3] + 2*m[4, 1, 1] + m[4, 2]
        """
        #Use symmetrica to do the multiplication
        A = left.parent()
        R = A.base_ring()

        #Hack due to symmetrica crashing when both of the
        #partitions are the empty partition
        #if  R is ZZ or R is QQ:
        #    return symmetrica.mult_monomial_monomial(left, right)

        z_elt = {}
        for (left_m, left_c) in left._monomial_coefficients.iteritems():
            for (right_m, right_c) in right._monomial_coefficients.iteritems():

                #Hack due to symmetrica crashing when both of the
                #partitions are the empty partition
                if left_m == [] and right_m == []:
                    z_elt[ left_m ] = left_c*right_c
                    continue

                d = symmetrica.mult_monomial_monomial({left_m:Integer(1)}, {right_m:Integer(1)}).monomial_coefficients()
                for m in d:
                    if m in z_elt:
                        z_elt[ m ] = z_elt[m] + left_c * right_c * d[m]
                    else:
                        z_elt[ m ] = left_c * right_c * d[m]
        return z_elt

    def from_polynomial(self, f, check=True):
        """
        Returns the symmetric function in the monomial basis corresponding to the polynomial ``f``.

        INPUT:

        - ``self`` -- a monomial symmetric function basis
        - ``f`` -- a polynomial in finitely many variables over the same base ring as ``self``.
          It is assumed that this polynomial is symmetric.
        - ``check`` -- boolean (default: True), checks whether the polynomial is indeed symmetric

        OUTPUT:

        - This function converts a symmetric polynomial `f` in a polynomial ring in finitely
          many variables to a symmetric function in the monomial
          basis of the ring of symmetric functions over the same base ring.

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: P = PolynomialRing(QQ, 'x', 3)
            sage: x = P.gens()
            sage: f = x[0] + x[1] + x[2]
            sage: m.from_polynomial(f)
            m[1]
            sage: f = x[0]**2+x[1]**2+x[2]**2
            sage: m.from_polynomial(f)
            m[2]
            sage: f=x[0]^2+x[1]
            sage: m.from_polynomial(f)
            Traceback (most recent call last):
            ...
            ValueError: x0^2 + x1 is not a symmetric polynomial
            sage: f = (m[2,1]+m[1,1]).expand(3)
            sage: m.from_polynomial(f)
            m[1, 1] + m[2, 1]
            sage: f = (2*m[2,1]+m[1,1]+3*m[3]).expand(3)
            sage: m.from_polynomial(f)
            m[1, 1] + 2*m[2, 1] + 3*m[3]
        """
        assert self.base_ring() == f.base_ring()
        out = self.sum_of_terms((Partition(e), c)
                                for (e,c) in f.dict().iteritems()
                                if tuple(sorted(e)) == tuple(reversed(e)))
        if check and out.expand(f.parent().ngens(),f.parent().gens()) <> f:
            raise ValueError, "%s is not a symmetric polynomial"%f
        return out

    def from_polynomial_exp(self, p):
        r"""
        Conversion from polynomial in exponential notation

        INPUT:

        - ``self`` -- a monomial symmetric function basis
        - ``p`` -- a multivariate polynomial over the same base ring as ``self``

        OUTPUT:

        - This returns a symmetric function by mapping each monomial of
          `p` with exponents ``exp`` into `m_\lambda` where `\lambda` is
          the partition with exponential notation ``exp``.

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: P = PolynomialRing(QQ,'x',5)
            sage: x = P.gens()

        The exponential notation of the partition `(5,5,5,3,1,1)` is::

            sage: Partition([5,5,5,3,1,1]).to_exp()
            [2, 0, 1, 0, 3]

        Therefore, the monomial::

            sage: f = x[0]^2 * x[2] * x[4]^3

        is mapped to::

            sage: m.from_polynomial_exp(f)
            m[5, 5, 5, 3, 1, 1]

        Furthermore, this function is linear::

            sage: f = 3 * x[3] + 2 * x[0]^2 * x[2] * x[4]^3
            sage: m.from_polynomial_exp(f)
            3*m[4] + 2*m[5, 5, 5, 3, 1, 1]

        ..SEEALSO:: :func:`Partition`, :meth:`Partition.to_exp`
        """
        assert self.base_ring() == p.parent().base_ring()
        return self.sum_of_terms((Partition(exp=monomial), coeff)
                                 for (monomial, coeff) in p.dict().iteritems())

    def coproduct_by_coercion(self, elt):
        r"""
        Returns the coproduct of the element ``elt`` by coercion to the Schur basis.

        INPUT:

        - ``self`` -- a monomial symmetric function basis
        - ``elt`` -- an instance of monomial basis

        OUTPUT:

        - The coproduct acting on ``elt``, the result is an element of the
          tensor squared of the monomial symmetric function basis


        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: m[3,1,1].coproduct()
            m[] # m[3, 1, 1] + m[1] # m[3, 1] + m[1, 1] # m[3] + m[3] # m[1, 1] + m[3, 1] # m[1] + m[3, 1, 1] # m[]
            sage: m.coproduct_by_coercion(m[2,1])
            m[] # m[2, 1] + m[1] # m[2] + m[2] # m[1] + m[2, 1] # m[]
            sage: m.coproduct_by_coercion(m[2,1]) == m([2,1]).coproduct()
            True
        """
        from sage.categories.tensor import tensor
        s = self.realization_of().schur()
        return self.tensor_square().sum(coeff * tensor([self(s[x]), self(s[y])])
                                        for ((x,y), coeff) in s(elt).coproduct())


    class Element(classical.SymmetricFunctionAlgebra_classical.Element):
        def expand(self, n, alphabet='x'):
            """
            Expands the symmetric function as a symmetric polynomial in `n` variables.

            INPUT:

            - ``self`` -- an element of the monomial symmetric function basis
            - ``n`` -- a positive integer
            - ``alphabet`` -- a variable for the expansion (default: `x`)

            OUTPUT: a monomial expansion of an instance of ``self`` in `n` variables

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).m()
                sage: m([2,1]).expand(3)
                x0^2*x1 + x0*x1^2 + x0^2*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2
                sage: m([1,1,1]).expand(2)
                0
                sage: m([2,1]).expand(3,alphabet='z')
                z0^2*z1 + z0*z1^2 + z0^2*z2 + z1^2*z2 + z0*z2^2 + z1*z2^2
                sage: m([2,1]).expand(3,alphabet='x,y,z')
                x^2*y + x*y^2 + x^2*z + y^2*z + x*z^2 + y*z^2
            """
            condition = lambda part: len(part) > n
            return self._expand(condition, n, alphabet)

# Backward compatibility for unpickling
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.sf.monomial', 'SymmetricFunctionAlgebraElement_monomial',  SymmetricFunctionAlgebra_monomial.Element)
