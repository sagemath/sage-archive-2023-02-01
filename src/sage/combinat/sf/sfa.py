r"""
Symmetric Functions

For a comprehensive tutorial on how to use symmetric functions in sage

.. SEEALSO:: :func:`SymmetricFunctions`

We define the algebra of symmetric functions in the Schur and elementary bases::

    sage: s = SymmetricFunctionAlgebra(QQ, basis='schur')
    sage: e = SymmetricFunctionAlgebra(QQ, basis='elementary')

Each is actually is a graded Hopf algebra whose basis is indexed by
integer partitions::

    sage: s.category()
    Category of bases of Symmetric Functions over Rational Field
    sage: s.basis().keys()
    Partitions

Let us compute with some elements::

    sage: f1 = s([2,1]); f1
    s[2, 1]
    sage: f2 = e(f1); f2
    e[2, 1] - e[3]
    sage: f1 == f2
    True
    sage: f1.expand(3, alphabet=['x','y','z'])
    x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2
    sage: f2.expand(3, alphabet=['x','y','z'])
    x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2

::

    sage: m = SymmetricFunctions(QQ).monomial()
    sage: m([3,1])
    m[3, 1]
    sage: m(4)
    4*m[]
    sage: m([4])
    m[4]
    sage: 3*m([3,1])-1/2*m([4])
    3*m[3, 1] - 1/2*m[4]

We can also do computations with various bases::

    sage: p = SymmetricFunctions(QQ).power()
    sage: f = p(3)
    sage: f
    3*p[]
    sage: f.parent()
    Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
    sage: f + p([3,2])
    3*p[] + p[3, 2]

One can convert symmetric functions to symmetric polynomials and vice versa::

    sage: Sym = SymmetricFunctions(QQ)
    sage: p = Sym.powersum()
    sage: h = Sym.homogeneous()
    sage: f = h[2,1] + 2*p[3,1]
    sage: poly = f.expand(3); poly
    2*x0^4 + 2*x0^3*x1 + 2*x0*x1^3 + 2*x1^4 + 2*x0^3*x2 + 2*x1^3*x2 + 2*x0*x2^3 + 2*x1*x2^3 + 2*x2^4
    + x0^3 + 2*x0^2*x1 + 2*x0*x1^2 + x1^3 + 2*x0^2*x2 + 3*x0*x1*x2 + 2*x1^2*x2 + 2*x0*x2^2 + 2*x1*x2^2 + x2^3
    sage: Sym.from_polynomial(poly)
    3*m[1, 1, 1] + 2*m[2, 1] + m[3] + 2*m[3, 1] + 2*m[4]
    sage: Sym.from_polynomial(poly) == f
    True

::

    sage: Sym = SymmetricFunctions(QQ)
    sage: s = Sym.s()
    sage: h = Sym.h()
    sage: p = Sym.p()
    sage: e = Sym.e()
    sage: m = Sym.m()
    sage: a = s([3,1])
    sage: s(a)
    s[3, 1]
    sage: h(a)
    h[3, 1] - h[4]
    sage: p(a)
    1/8*p[1, 1, 1, 1] + 1/4*p[2, 1, 1] - 1/8*p[2, 2] - 1/4*p[4]
    sage: e(a)
    e[2, 1, 1] - e[2, 2] - e[3, 1] + e[4]
    sage: m(a)
    3*m[1, 1, 1, 1] + 2*m[2, 1, 1] + m[2, 2] + m[3, 1]
    sage: a.expand(4)
    x0^3*x1 + x0^2*x1^2 + x0*x1^3 + x0^3*x2 + 2*x0^2*x1*x2 + 2*x0*x1^2*x2 + x1^3*x2 + x0^2*x2^2 + 2*x0*x1*x2^2 + x1^2*x2^2 + x0*x2^3 + x1*x2^3 + x0^3*x3 + 2*x0^2*x1*x3 + 2*x0*x1^2*x3 + x1^3*x3 + 2*x0^2*x2*x3 + 3*x0*x1*x2*x3 + 2*x1^2*x2*x3 + 2*x0*x2^2*x3 + 2*x1*x2^2*x3 + x2^3*x3 + x0^2*x3^2 + 2*x0*x1*x3^2 + x1^2*x3^2 + 2*x0*x2*x3^2 + 2*x1*x2*x3^2 + x2^2*x3^2 + x0*x3^3 + x1*x3^3 + x2*x3^3

Here are further examples::

    sage: h(m([1]))
    h[1]
    sage: h( m([2]) +m([1,1]) )
    h[2]
    sage: h( m([3]) + m([2,1]) + m([1,1,1]) )
    h[3]
    sage: h( m([4]) + m([3,1]) + m([2,2]) + m([2,1,1]) + m([1,1,1,1]) )
    h[4]
    sage: k = 5
    sage: h( sum([ m(part) for part in Partitions(k)]) )
    h[5]
    sage: k = 10
    sage: h( sum([ m(part) for part in Partitions(k)]) )
    h[10]

::

    sage: P3 = Partitions(3)
    sage: P3.list()
    [[3], [2, 1], [1, 1, 1]]
    sage: m = SymmetricFunctions(QQ).monomial()
    sage: f = sum([m(p) for p in P3])
    sage: m.get_print_style()
    'lex'
    sage: f
    m[1, 1, 1] + m[2, 1] + m[3]
    sage: m.set_print_style('length')
    sage: f
    m[3] + m[2, 1] + m[1, 1, 1]
    sage: m.set_print_style('maximal_part')
    sage: f
    m[1, 1, 1] + m[2, 1] + m[3]
    sage: m.set_print_style('lex')

::

    sage: Sym = SymmetricFunctions(QQ)
    sage: s = Sym.s()
    sage: m = Sym.m()
    sage: m([3])*s([2,1])
    2*m[3, 1, 1, 1] + m[3, 2, 1] + 2*m[4, 1, 1] + m[4, 2] + m[5, 1]
    sage: s(m([3])*s([2,1]))
    s[2, 1, 1, 1, 1] - s[2, 2, 2] - s[3, 3] + s[5, 1]
    sage: s(s([2,1])*m([3]))
    s[2, 1, 1, 1, 1] - s[2, 2, 2] - s[3, 3] + s[5, 1]
    sage: e = Sym.e()
    sage: e([4])*e([3])*e([1])
    e[4, 3, 1]

::

    sage: s = SymmetricFunctions(QQ).s()
    sage: z = s([2,1]) + s([1,1,1])
    sage: z.coefficient([2,1])
    1
    sage: z.length()
    2
    sage: z.support()
    [[1, 1, 1], [2, 1]]
    sage: z.degree()
    3

RECENT BACKWARD INCOMPATIBLE CHANGES:

The symmetric functions code has been refactored to take
advantate of the coercion systems. This introduced a couple glitches:

- On some bases changes, coefficients in Jack polynomials are not normalized

- Except in a few cases, conversions and coercions are only defined
  between symmetric functions over the same coefficient ring. E.g.
  the following does not work anymore::

      sage: s  = SymmetricFunctions(QQ)
      sage: s2 = SymmetricFunctions(QQ['t'])
      sage: s([1]) + s2([2]) # todo: not implemented

  This feature will probably come back at some point through
  improvements to the Sage coercion system.

Backward compatibility should be essentially retained.

AUTHORS:

- Mike Hansen (2007-06-15)
- Nicolas M. Thiery (partial refactoring)
- Mike Zabrocki, Anne Schilling (2012)

"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Anne Schilling <anne at math.ucdavis.edu>
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
from sage.misc.cachefunc import cached_method
from sage.rings.all import Integer, PolynomialRing, is_Polynomial, is_MPolynomial, ZZ, QQ
import sage.combinat.partition
from sage.combinat.partition import Partitions
import sage.libs.symmetrica.all as symmetrica  # used in eval()
from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.constructor import matrix
from sage.misc.misc import repr_lincomb, prod, uniq
from functools import partial
from copy import copy


def SymmetricFunctionAlgebra(R, basis="schur"):
    r"""
    Return the free algebra over the ring `R` on `n`
    generators with given names.

    INPUT:

    -  ``R`` -- ring with identity basis
    -  ``basis`` -- a string for the name of the basis, must be one of
       'schur', 'elementary', 'homogeneous', 'power', 'monomial' or their
       abbreviations 's', 'e', 'h', 'p', 'm'

    OUTPUT: A SymmetricFunctionAlgebra

    EXAMPLES::

        sage: SymmetricFunctionAlgebra(QQ)
        Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis

    ::

        sage: SymmetricFunctionAlgebra(QQ, basis='m')
        Symmetric Function Algebra over Rational Field, Monomial symmetric functions as basis

    ::

        sage: SymmetricFunctionAlgebra(QQ, basis='power')
        Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
    """
    # Todo: this is a backward compatibility function, and should be deprecated
    from sage.combinat.sf.sf import SymmetricFunctions
    Sym = SymmetricFunctions(R)
    if basis == 'schur' or basis == 's':
        return Sym.s()
    elif basis == "elementary" or  basis ==  'e':
        return Sym.e()
    elif basis == "homogeneous" or basis ==  'h':
        return Sym.h()
    elif basis == 'power' or basis ==  'p':
        return Sym.p()
    elif basis == 'monomial' or basis ==  'm':
        return Sym.m()
    else:
        raise ValueError, "unknown basis (= %s)"%basis

def SFAPower(R):
    """
    Returns the symmetric function algebra over ``R`` with the power-sum
    symmetric functions as the basis.

    This function is deprecated.  Use instead:
    SymmetricFunctions(R).power()

    EXAMPLES::

        sage: SymmetricFunctions(QQ).power()
        Symmetric Function Algebra over Rational Field, Power symmetric functions as basis

    TESTS::

        sage: p = SFAPower(QQ)
        doctest:1: DeprecationWarning: Deprecation warning: In the future use SymmetricFunctions(R).power()
        See http://trac.sagemath.org/5457 for details.
    """
    sage.misc.superseded.deprecation(5457,"Deprecation warning: In the future use SymmetricFunctions(R).power()")
    return SymmetricFunctionAlgebra(R, basis='power')


def SFAElementary(R):
    """
    Returns the symmetric function algebra over ``R`` with the elementary
    symmetric functions as the basis.

    This function is deprecated.  Use instead:
    SymmetricFunctions(R).elementary()

    EXAMPLES::

        sage: SymmetricFunctions(QQ).elementary()
        Symmetric Function Algebra over Rational Field, Elementary symmetric functions as basis

    TESTS::

        sage: p = SFAElementary(QQ)
        doctest:1: DeprecationWarning: Deprecation warning: In the future use SymmetricFunctions(R).elementary()
        See http://trac.sagemath.org/5457 for details.
    """
    sage.misc.superseded.deprecation(5457,"Deprecation warning: In the future use SymmetricFunctions(R).elementary()")
    return SymmetricFunctionAlgebra(R, basis='elementary')


def SFAHomogeneous(R):
    """
    Returns the symmetric function algebra over R with the Homogeneous
    symmetric functions as the basis.

    This function is deprecated.  Use instead:
    SymmetricFunctions(R).homogeneous()

    EXAMPLES::

        sage: SymmetricFunctions(QQ).homogeneous()
        Symmetric Function Algebra over Rational Field, Homogeneous symmetric functions as basis

    TESTS::

        sage: p = SFAHomogeneous(QQ)
        doctest:1: DeprecationWarning: Deprecation warning: In the future use SymmetricFunctions(R).homogeneous()
        See http://trac.sagemath.org/5457 for details.
    """
    sage.misc.superseded.deprecation(5457,"Deprecation warning: In the future use SymmetricFunctions(R).homogeneous()")
    return SymmetricFunctionAlgebra(R, basis='homogeneous')


def SFASchur(R):
    """
    Returns the symmetric function algebra over R with the Schur
    symmetric functions as the basis.

    This function is deprecated.  Use instead:
    SymmetricFunctions(R).schur()

    EXAMPLES::

        sage: SymmetricFunctions(QQ).schur()
        Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis

    TESTS::

        sage: p = SFASchur(QQ)
        doctest:1: DeprecationWarning: Deprecation warning: In the future use SymmetricFunctions(R).schur()
        See http://trac.sagemath.org/5457 for details.
    """
    sage.misc.superseded.deprecation(5457,"Deprecation warning: In the future use SymmetricFunctions(R).schur()")
    return SymmetricFunctionAlgebra(R, basis='schur')


def SFAMonomial(R):
    """
    Returns the symmetric function algebra over R with the monomial
    symmetric functions as the basis.

    This function is deprecated.  Use instead:
    SymmetricFunctions(R).homogeneous()

    EXAMPLES::

        sage: SymmetricFunctions(QQ).monomial()
        Symmetric Function Algebra over Rational Field, Monomial symmetric functions as basis

    TESTS::

        sage: p = SFAMonomial(QQ)
        doctest:1: DeprecationWarning: Deprecation warning: In the future use SymmetricFunctions(R).monomial()
        See http://trac.sagemath.org/5457 for details.
    """
    sage.misc.superseded.deprecation(5457, "Deprecation warning: In the future use SymmetricFunctions(R).monomial()")
    return SymmetricFunctionAlgebra(R, basis='monomial')

def is_SymmetricFunctionAlgebra(x):
    """
    Checks whether ``x`` is a symmetric function algebra.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import is_SymmetricFunctionAlgebra
        sage: is_SymmetricFunctionAlgebra(5)
        False
        sage: is_SymmetricFunctionAlgebra(ZZ)
        False
        sage: is_SymmetricFunctionAlgebra(SymmetricFunctionAlgebra(ZZ,'schur'))
        True
        sage: is_SymmetricFunctionAlgebra(SymmetricFunctions(QQ).e())
        True
        sage: is_SymmetricFunctionAlgebra(SymmetricFunctions(QQ).macdonald(q=1,t=1).P())
        True
        sage: is_SymmetricFunctionAlgebra(SymmetricFunctions(FractionField(QQ['q','t'])).macdonald().P())
        True
    """
    return isinstance(x, SymmetricFunctionAlgebra_generic)

def zee(part):
    r"""
    Returns the size of the centralizer of permutations of cycle type ``part``.

    Note that the size of the centralizer is the inner product between `p(part)` and
    itself where `p` is the power-sum symmetric functions.

    INPUT:

    -  ``part`` -- an integer partition (for example, [2,1,1])

    OUTPUT:

    - the integer `\prod_{i} i^{m_i(part)} m_i(part)!` where `m_i(part)` is
      the number of parts in the partition ``part`` equal to `i`

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import zee
        sage: zee([2,1,1])
        4
    """
    if not isinstance(part, sage.combinat.partition.Partition_class):
        part = sage.combinat.partition.Partition_class(part)
    return part.centralizer_size()


def is_SymmetricFunction(x):
    r"""
    Checks whether ``x`` is a symmetric function.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import is_SymmetricFunction
        sage: s = SymmetricFunctions(QQ).s()
        sage: is_SymmetricFunction(2)
        False
        sage: is_SymmetricFunction(s(2))
        True
        sage: is_SymmetricFunction(s([2,1]))
        True
    """
    return isinstance(x, SymmetricFunctionAlgebra_generic.Element)

from sage.categories.realizations import Realizations, Category_realization_of_parent
class SymmetricFunctionsBases(Category_realization_of_parent):
    r"""
    The category of bases of symmetric functions.
    """
    def __init__(self, base):
        r"""
        Initializes the bases of the ring of symmetric functions.

        INPUT:

        - ``self`` -- a category of bases for the symmetric functions
        - ``base`` -- ring of symmetric functions

        TESTS::

            sage: from sage.combinat.sf.sfa import SymmetricFunctionsBases
            sage: Sym = SymmetricFunctions(QQ)
            sage: bases = SymmetricFunctionsBases(Sym); bases
            Category of bases of Symmetric Functions over Rational Field
            sage: Sym.schur() in bases
            True
        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Returns the representation of ``self``.

        INPUT:

        - ``self`` -- a category of bases for the symmetric functions

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import SymmetricFunctionsBases
            sage: Sym = SymmetricFunctions(QQ)
            sage: bases = SymmetricFunctionsBases(Sym)
            sage: bases._repr_()
            'Category of bases of Symmetric Functions over Rational Field'
        """
        return "Category of bases of %s" % self.base()

    def super_categories(self):
        r"""
        The super categories of ``self``.

        INPUT:

        - ``self`` -- a category of bases for the symmetric functions

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import SymmetricFunctionsBases
            sage: Sym = SymmetricFunctions(QQ)
            sage: bases = SymmetricFunctionsBases(Sym)
            sage: bases.super_categories()
            [Category of graded hopf algebras with basis over Rational Field, Category of realizations of Symmetric Functions over Rational Field, Category of commutative rings]
        """
        from sage.categories.all import CommutativeRings, GradedHopfAlgebrasWithBasis
        return [GradedHopfAlgebrasWithBasis(self.base().base_ring()),
                Realizations(self.base()),
                CommutativeRings()]

    class ParentMethods:
        def antipode_by_coercion(self, element):
            r"""
            The antipode of ``element`` via coercion to and from the powersum basis.

            INPUT:

            - ``self`` -- a basis of the ring of symmetric functions
            - ``element`` -- element in a basis of the ring of symmetric functions

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ)
                sage: p = Sym.p()
                sage: s = Sym.s()
                sage: h = Sym.h()
                sage: (h([]) + h([1])).antipode() # indirect doctest
                h[] - h[1]
                sage: (s([]) + s([1]) + s[2]).antipode()
                s[] - s[1] + s[1, 1]
                sage: (p([2]) + p([3])).antipode()
                -p[2] - p[3]
            """
            p = self.realization_of().powersum()
            return self(p.antipode(p(element)))

        def counit(self, element):
            r"""
            Returns the counit of ``element``.

            INPUT:

            - ``self`` -- a basis of the ring of symmetric functions
            - ``element`` -- element in a basis of the ring of symmetric functions

            The counit is the constant terms of ``element``.

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ)
                sage: m = Sym.monomial()
                sage: f = 2*m[2,1] + 3*m[[]]
                sage: f.counit()
                3
            """
            return element.degree_zero_coefficient()

    class ElementMethods:

        def degree_zero_coefficient(self):
            r"""
            Returns the degree zero coefficient of ``self``.

            INPUT:

            - ``self`` -- an element of the symmetric functions

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ)
                sage: m = Sym.monomial()
                sage: f = 2*m[2,1] + 3*m[[]]
                sage: f.degree_zero_coefficient()
                3
            """
            return self.coefficient([])

class SymmetricFunctionAlgebra_generic(CombinatorialFreeModule):
    r"""
    TESTS::

        sage: s = SymmetricFunctions(QQ).s()
        sage: m = SymmetricFunctions(ZZ).m()
        sage: s(m([2,1]))
        -2*s[1, 1, 1] + s[2, 1]
    """

    def __init__(self, Sym, basis_name = None, prefix = None):
        r"""
        Initializes the symmetric function algebra.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``Sym`` -- the ring of symmetric functions
        - ``basis_name`` -- name of basis (default: None)
        - ``prefix`` -- prefix used to display basis

        TESTS::

            sage: from sage.combinat.sf.classical import SymmetricFunctionAlgebra_classical
            sage: s = SymmetricFunctions(QQ).s()
            sage: isinstance(s, SymmetricFunctionAlgebra_classical)
            True
            sage: TestSuite(s).run()
        """
        R = Sym.base_ring()
        from sage.categories.all import CommutativeRings
        if R not in CommutativeRings():
            raise TypeError, "Argument R must be a commutative ring."
        try:
            z = R(Integer(1))
        except StandardError:
            raise ValueError, "R must have a unit element"

        if basis_name is not None:
            self._basis = basis_name
        if prefix is not None:
            self._prefix = prefix
        self._sym = Sym
        CombinatorialFreeModule.__init__(self, Sym.base_ring(), sage.combinat.partition.Partitions(),
                                         category = SymmetricFunctionsBases(Sym),
                                         bracket = "", prefix = prefix)

    @cached_method
    def one_basis(self):
        r"""
        Returns the empty partition, as per ``AlgebrasWithBasis.ParentMethods.one_basis``

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions

        EXAMPLES::

            sage: SymmetricFunctions(QQ).s().one_basis()
            []

        TESTS::

            sage: SymmetricFunctions(QQ).s().one_basis() == Partition([])
            True
        """
        return sage.combinat.partition.Partitions()([])

    _print_style = 'lex'

    # Todo: share this with ncsf and over algebras with basis indexed by word-like elements
    def __getitem__(self, c, *rest):
        r"""
        This method implements the abuses of notations `p[2,1]`, `p[[2,1]]`, `p[Partition([2,1])]`

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``c`` -- list, list of list, or partition

        Todo: should call super.term so as not to interfer with the standard notation p['x,y,z']

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s[2,1]
            s[2, 1]
            sage: s[[2,1]]
            s[2, 1]
            sage: s[Partition([2,1])]
            s[2, 1]
        """
        C = self.basis().keys()
        if isinstance(c, C.element_class):
            assert len(rest) == 0
        else:
            if len(rest) > 0 or type(c) is int or type(c) is Integer:
                c = C([c]+list(rest))
            else:
                c = C(list(c))
        return self.monomial(c)

    def _change_by_proportionality(self, x, function):
        r"""
        Return the symmetric function obtained by scaling each basis
        element corresponding to the partition `part` by ``function``(`part`).

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``x`` -- a symmetric function
        - ``function`` -- a function which takes in a partition and returns a scalar

        OUTPUT: a symmetric function in ``self`` which is a scaled version of ``x``

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([3])+s([2,1])+s([1,1,1]); a
            s[1, 1, 1] + s[2, 1] + s[3]
            sage: f = lambda part: len(part)
            sage: s._change_by_proportionality(a, f)
            3*s[1, 1, 1] + 2*s[2, 1] + s[3]
        """
        BR = self.base_ring()
        z_elt = {}
        for m, c in x._monomial_coefficients.iteritems():
            coeff = function(m)
            z_elt[m] = BR( c*coeff )
        return self._from_dict(z_elt)

    def _change_by_plethysm(self, x, expr, deg_one):
        r"""
        Returns the plethysm of ``x`` by ``expr``.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``x` -- symmetric function
        - ``expr`` -- an expression used in the plethysm
        - ``deg_one`` -- specifies the degree one terms

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: a = m([2,1])
            sage: a.omega()
            -m[2, 1] - 2*m[3]
            sage: m._change_by_plethysm(-a,-1,[])
            -m[2, 1] - 2*m[3]

        ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([3])
            sage: s._change_by_plethysm(-a,-1,[])
            s[1, 1, 1]
        """
        #Covert to the power sum
        p = sage.combinat.sf.sf.SymmetricFunctions(self.base_ring()).power()
        p_x = p(x)
        expr_k = lambda k: expr.subs(**dict([(str(x),x**k) for x in deg_one]))
        f = lambda m,c: (m, c*prod([expr_k(k) for k in m]))
        return self(p_x.map_item(f))

    # TODO:
    #  - lift to combinatorial_module
    #  - rename to _apply_bimodule_morphism or generalize to true multi_module
    #  - generalization with a "neighbor function" that given x says
    #    for which y one has f(x,y) != 0
    #  - add option orthonormal
    def _apply_multi_module_morphism(self, x, y, f, orthogonal=False):
        r"""
        Applies morphism specified by ``f`` on (``x``,``y``).

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``x`` -- an element of ``self``
        - ``y`` -- an element of ``self``
        - ``f`` -- a function that takes in two partitions
                   (basis elements) and returns an element of the target domain
        - ``orthogonal`` -- if orthogonal is set to True,
                             then `f(part1, part2)` is assumed to be 0 if part1 != part2.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([2,1])+s([1,1,1])
            sage: b = s([3])+s([2,1])
            sage: f1 = lambda p1, p2: len(p1)*len(p2)
            sage: f2 = lambda p1, p2: len(p1)+len(p2)
            sage: s._apply_multi_module_morphism(a,b,f1,orthogonal=False) #(2+3)*(2+1)
            15
            sage: s._apply_multi_module_morphism(a,b,f1,orthogonal=True)  #(2)*(2)
            4
            sage: s._apply_multi_module_morphism(a,b,f2,orthogonal=False) #2*(2+3+2+1)
            16
            sage: s._apply_multi_module_morphism(a,b,f2,orthogonal=True)  #2+2
            4
        """
        # broken for most coeff ring
        res = 0
        if orthogonal:
            # could check which of x and y has less terms
            # for mx, cx in x:
            for mx, cx in x._monomial_coefficients.iteritems():
                # if not y.has_key(mx):
                if mx not in y._monomial_coefficients:
                    continue
                else:
                    # cy = y[mx]
                    cy = y._monomial_coefficients[mx]
                # might as well call f(mx)
                res += cx*cy*f(mx, mx)
            return res
        else:
            for mx, cx in x._monomial_coefficients.iteritems():
                for my, cy in y._monomial_coefficients.iteritems():
                    res += cx*cy*f(mx,my)
            return res

    def _from_element(self, x):
        r"""
        Return the element of ``self`` with the same 'internal structure' as ``x``.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``x`` -- a symmetric function

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).e()
            sage: s = SymmetricFunctions(QQ).s()
            sage: a = e([2,1]) + e([1,1,1]); a
            e[1, 1, 1] + e[2, 1]
            sage: s._from_element(a)
            s[1, 1, 1] + s[2, 1]
        """
        return self._from_dict(x.monomial_coefficients())

    def _from_cache(self, element, cache_function, cache_dict, **subs_dict):
        r"""
        Return an ``element`` of ``self`` from cache.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions

        -  ``element`` -- a symmetric function

        -  ``cache_function`` -- a function which accepts an
           integer n as its input and creates the cache for that homogenous
           component

        -  ``cache_dict`` -- the dictionary storing the cache;
           it is indexed by the positive integers n, and it's values are
           dictionaries indexed by partitions size n. The values of those
           dictionaries are in dictionaries indexed by partitions of size n.

        -  ``subs_dict`` -- a dictionary for any substitutions
           to make after the value is extracted from cache_dict.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: Sym = SymmetricFunctions(R)
            sage: s = Sym.s()
            sage: p21 = Partition([2,1])
            sage: a = s(p21)
            sage: e = Sym.e()
            sage: cache_dict = {}
            sage: cache_dict[3] = {}
            sage: cache_dict[3][p21] = {}
            sage: cache_dict[3][p21][p21] = x^2
            sage: cache_dict[3][p21][Partition([1,1,1])] = 3*x
            sage: cache_function = lambda n: 0 #do nothing
            sage: e._from_cache(a, cache_function, cache_dict)
            3*x*e[1, 1, 1] + x^2*e[2, 1]
            sage: e._from_cache(a, cache_function, cache_dict, x=2)
            6*e[1, 1, 1] + 4*e[2, 1]
        """
        #Convert x to the monomial basis
        BR = self.base_ring()
        zero = BR(0)
        z_elt = {}
        for part,c in element.monomial_coefficients().iteritems():
            if sum(part) not in cache_dict:
                cache_function(sum(part))
            for part2, c2 in cache_dict[sum(part)][part].iteritems():
                if hasattr(c2,'subs'): # c3 may be in the base ring
                    c3 = c*BR(c2.subs(**subs_dict))
                else:
                    c3 = c*BR(c2)
#                c3 = c*c2
#                if hasattr(c3,'subs'): # c3 may be in the base ring
#                    c3 = c3.subs(**subs_dict)
                z_elt[ part2 ] = z_elt.get(part2, zero) + BR(c3)
        return self._from_dict(z_elt)

    def _invert_morphism(self, n, base_ring, self_to_other_cache, other_to_self_cache,\
                         to_other_function=None, to_self_function=None, \
                         upper_triangular=False, lower_triangular=False, \
                         ones_on_diagonal=False):
        r"""
        Returns the inverse of a morphism between ``self`` and ``other``.

        In order to use this, you must be able compute the morphism in one
        direction. This method assumes that the morphism is indeed
        invertible.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions

        -  ``n`` -- an integer, the homogeneous component of
           symmetric functions for which we want to a morphism's inverse

        -  ``base_ring`` -- the base ring being worked over

        -  ``self_to_other_cache`` -- a dictionary which
           stores the transition from ``self`` to ``other``

        -  ``other_to_self_cache`` -- a dictionary which
           stores the transition from ``other`` to ``self``

        -  ``to_other_function`` -- a function which takes in
           a partition and returns a function which gives the coefficients of
           self(part) in the other basis

        -  ``to_self_function`` -- a function which takes in a
           partition and returns a function which gives the coefficients of
           other(part) in self

        -  ``upper_triangular`` -- a boolean, if True, the
           inverse will be computed by back substitution

        -  ``lower_triangular`` -- a boolean, if True, the
           inverse will be computed by forward substitution

        -  ``ones_on_diagonal`` -- a boolean, if True, the
           entries on the diagonal of the morphism (and inverse) matrix are
           assumed to be one. This is used to remove divisions from the
           forward and back substitute algorithms.


        EXAMPLES: First, we will do an example of inverting the morphism
        which sends a Schur function to its conjugate Schur function. Note
        that this is an involution.

        ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: conj = lambda p1: lambda p2: QQ(1) if p2 == p1.conjugate() else QQ(0)
            sage: c1 = {}
            sage: c2 = {}
            sage: s._invert_morphism(4, QQ, c1, c2, to_other_function = conj)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(c1[4])
            [([1, 1, 1, 1], [([4], 1)]),
             ([2, 1, 1], [([3, 1], 1)]),
             ([2, 2], [([2, 2], 1)]),
             ([3, 1], [([2, 1, 1], 1)]),
             ([4], [([1, 1, 1, 1], 1)])]
            sage: l(c2[4])
            [([1, 1, 1, 1], [([4], 1)]),
             ([2, 1, 1], [([3, 1], 1)]),
             ([2, 2], [([2, 2], 1)]),
             ([3, 1], [([2, 1, 1], 1)]),
             ([4], [([1, 1, 1, 1], 1)])]
            sage: c2 == c1
            True

        We can check that we get the same results if specify
        to_self_function = conj.

        ::

            sage: d1 = {}
            sage: d2 = {}
            sage: s._invert_morphism(4, QQ, d1, d2, to_self_function = conj)
            sage: d1 == c1
            True
            sage: d2 == c2
            True

        Now we do an example of upper triangularity and check that we get
        the same thing whether or not we specify ones_on_diagonal.

        ::

            sage: f = lambda p1: lambda p2: QQ(1) if p2 <= p1 else QQ(0)
            sage: c1 = {}
            sage: c2 = {}
            sage: s._invert_morphism(3, QQ, c1, c2, to_other_function = f, upper_triangular=True)
            sage: l(c1[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], 1), ([2, 1], 1)]),
             ([3], [([1, 1, 1], 1), ([2, 1], 1), ([3], 1)])]
            sage: l(c2[3])
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], -1), ([2, 1], 1)]),
             ([3], [([2, 1], -1), ([3], 1)])]

        ::

            sage: d1 = {}
            sage: d2 = {}
            sage: s._invert_morphism(3, QQ, d1, d2, to_other_function = f,upper_triangular=True, ones_on_diagonal=True)
            sage: c1 == d1
            True
            sage: c2 == d2
            True

        Finally, we do the same thing for lower triangular matrices.

        ::

            sage: f = lambda p1: lambda p2: QQ(1) if p2 >= p1 else QQ(0)
            sage: c1 = {}
            sage: c2 = {}
            sage: s._invert_morphism(3, QQ, c1, c2, to_other_function = f, lower_triangular=True)
            sage: l(c1[3])
            [([1, 1, 1], [([1, 1, 1], 1), ([2, 1], 1), ([3], 1)]),
             ([2, 1], [([2, 1], 1), ([3], 1)]),
             ([3], [([3], 1)])]

        ::

            sage: l(c2[3])
            [([1, 1, 1], [([1, 1, 1], 1), ([2, 1], -1)]),
             ([2, 1], [([2, 1], 1), ([3], -1)]),
             ([3], [([3], 1)])]

        ::

            sage: d1 = {}
            sage: d2 = {}
            sage: s._invert_morphism(3, QQ, d1, d2, to_other_function = f,lower_triangular=True, ones_on_diagonal=True)
            sage: c1 == d1
            True
            sage: c2 == d2
            True
        """
        #Decide whether we know how to go from self to other or
        #from other to self
        if to_other_function is not None:
            known_cache    = self_to_other_cache  #the known direction
            unknown_cache  = other_to_self_cache  #the unknown direction
            known_function = to_other_function
        else:
            unknown_cache  = self_to_other_cache  #the known direction
            known_cache    = other_to_self_cache  #the unknown direction
            known_function = to_self_function

        #Do nothing if we've already computed the inverse
        #for degree n.
        if n in known_cache and n in unknown_cache:
            return

        #Univariate polynomial arithmetic is faster
        #over ZZ.  Since that is all we need to compute
        #the transition matrices between S and P, we
        #should use that.
        #Zt = ZZ['t']
        #t = Zt.gen()
        one  = base_ring(1)
        zero = base_ring(0)

        #Get and store the list of partition we'll need
        pn = sage.combinat.partition.Partitions_n(n).list()
        len_pn = len(pn)

        #Create the initial cache dictionaries
        known_cache_n = {}
        known_matrix_n = matrix(base_ring, len_pn, len_pn)
        unknown_cache_n = {}
        for i in range(len_pn):
            known_cache_part = {}
            f = known_function(pn[i])
            for j in range(len_pn):
                if lower_triangular and j>i:
                    break
                if upper_triangular and i>j:
                    continue
                value = f(pn[j])
                if value != zero:
                    known_cache_part[ pn[ j ] ] = value
                    known_matrix_n[i,j] = value
            known_cache_n[ pn[i] ] = known_cache_part

            unknown_cache_n[ pn[i] ] = {}

        #Compute the inverse of the matrix
        if upper_triangular is not False and lower_triangular is not False:
            raise ValueError, "only one of upper_triangular and lower_triangular can be specified"
        elif upper_triangular is not False:
            #Compute the inverse of by using back
            #substitution.  We solve a len(pn) systems of
            #equations known_matrix_n*x = b_i for x, where e_i
            #is the ith standard basis vector
            inverse = copy(known_matrix_n.parent().zero_matrix())

            delta = lambda i: lambda j: one if i == j else zero

            for column in range(len_pn):
                e = delta(column)
                x = [0]*len_pn
                for i in range(len_pn-1,-1,-1):
                    value = e(i)
                    if not ones_on_diagonal:
                        value /= known_matrix_n[i,i]
                    for j in range(i+1,len_pn):
                        if ones_on_diagonal:
                            value -= known_matrix_n[i,j]*x[j]
                        else:
                            value -= known_matrix_n[i,j]*x[j]/known_matrix_n[i,i]
                    x[i] = value

                for j in range(column+1):
                    if x[j] != zero:
                        inverse[j,column] = x[j]

        elif lower_triangular is not False:
            #Compute the inverse of by using forward
            #substitution.  We solve a len(pn) systems of
            #equations known_matrix_n*x = b_i for x, where e_i
            #is the ith standard basis vector
            inverse = copy(known_matrix_n.parent().zero_matrix())


            delta = lambda i: lambda j: one if i == j else zero

            for column in range(len_pn):
                e = delta(column)
                x = []
                for i in range(len_pn):
                    value = e(i)
                    if not ones_on_diagonal:
                        value /= known_matrix_n[i,i]
                    for j in range(len(x)):
                        if ones_on_diagonal:
                            value -= known_matrix_n[i,j]*x[j]
                        else:
                            value -= known_matrix_n[i,j]*x[j]/known_matrix_n[i,i]
                    x.append(value)
                for j in range(column,len(x)):
                    if x[j] != zero:
                        inverse[j,column] = x[j]

        else:
            inverse = ~known_matrix_n

        for i in range(len_pn):
            for j in range(len_pn):
                if inverse[i,j] != zero:
                    if hasattr(self, '_normalize_coefficients'):
                        unknown_cache_n[ pn[i] ] [ pn[j] ] = self._normalize_coefficients(inverse[i,j])
                    else:
                        unknown_cache_n[ pn[i] ] [ pn[j] ] = inverse[i,j]

        known_cache[ n ]   = known_cache_n
        unknown_cache[ n ] = unknown_cache_n

    def symmetric_function_ring(self):
        r"""
        Returns the family of symmetric functions associated to the basis ``self``

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions

        OUTPUT:

        - returns an instance of the ring of symmetric functions

        EXAMPLES::

            sage: schur = SymmetricFunctions(QQ).schur()
            sage: schur.symmetric_function_ring()
            Symmetric Functions over Rational Field
            sage: power = SymmetricFunctions(QQ['t']).power()
            sage: power.symmetric_function_ring()
            Symmetric Functions over Univariate Polynomial Ring in t over Rational Field
        """
        return self.realization_of()

    def prefix(self):
        r"""
        Returns the prefix on the elements of ``self``.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions

        EXAMPLES::

            sage: schur = SymmetricFunctions(QQ).schur()
            sage: schur([3,2,1])
            s[3, 2, 1]
            sage: schur.prefix()
            's'
        """
        return self._prefix

    def transition_matrix(self, basis, n):
        r"""
        Returns the transitions matrix between ``self`` and ``basis`` for the
        homogenous component of degree ``n``.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``basis`` -- a basis of the ring of symmetric functions
        - ``n`` -- a nonnegative integer

        OUTPUT:

        - a matrix of coefficients giving the expansion of elements of ``self`` in the expansion of elements of ``basis``

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: m = SymmetricFunctions(QQ).m()
            sage: s.transition_matrix(m,5)
            [1 1 1 1 1 1 1]
            [0 1 1 2 2 3 4]
            [0 0 1 1 2 3 5]
            [0 0 0 1 1 3 6]
            [0 0 0 0 1 2 5]
            [0 0 0 0 0 1 4]
            [0 0 0 0 0 0 1]
            sage: s.transition_matrix(m,1)
            [1]
            sage: s.transition_matrix(m,0)
            [1]

        ::

            sage: p = SymmetricFunctions(QQ).p()
            sage: s.transition_matrix(p, 4)
            [ 1/4  1/3  1/8  1/4 1/24]
            [-1/4    0 -1/8  1/4  1/8]
            [   0 -1/3  1/4    0 1/12]
            [ 1/4    0 -1/8 -1/4  1/8]
            [-1/4  1/3  1/8 -1/4 1/24]
            sage: StoP = s.transition_matrix(p,4)
            sage: a = s([3,1])+5*s([1,1,1,1])-s([4])
            sage: a
            5*s[1, 1, 1, 1] + s[3, 1] - s[4]
            sage: mon = a.support()
            sage: coeffs = a.coefficients()
            sage: coeffs
            [5, 1, -1]
            sage: mon
            [[1, 1, 1, 1], [3, 1], [4]]
            sage: cm = matrix([[-1,1,0,0,5]])
            sage: cm * StoP
            [-7/4  4/3  3/8 -5/4 7/24]
            sage: p(a)
            7/24*p[1, 1, 1, 1] - 5/4*p[2, 1, 1] + 3/8*p[2, 2] + 4/3*p[3, 1] - 7/4*p[4]

        ::

            sage: h = SymmetricFunctions(QQ).h()
            sage: e = SymmetricFunctions(QQ).e()
            sage: s.transition_matrix(m,7) == h.transition_matrix(s,7).transpose()
            True

        ::

            sage: h.transition_matrix(m, 7) == h.transition_matrix(m, 7).transpose()
            True

        ::

            sage: h.transition_matrix(e, 7) == e.transition_matrix(h, 7)
            True

        ::

            sage: p.transition_matrix(s, 5)
            [ 1 -1  0  1  0 -1  1]
            [ 1  0 -1  0  1  0 -1]
            [ 1 -1  1  0 -1  1 -1]
            [ 1  1 -1  0 -1  1  1]
            [ 1  0  1 -2  1  0  1]
            [ 1  2  1  0 -1 -2 -1]
            [ 1  4  5  6  5  4  1]

        ::

            sage: e.transition_matrix(m,7) == e.transition_matrix(m,7).transpose()
            True
        """
        P = sage.combinat.partition.Partitions_n(n)
        Plist = P.list()
        m = []
        for row_part in Plist:
            z = basis(self(row_part))
            m.append( map( lambda col_part: z.coefficient(col_part), Plist ) )
        return matrix(m)


    def _gram_schmidt(self, n, source, scalar, cache, leading_coeff=None, upper_triangular=True):
        r"""
        Returns Gram-Schmidt on ``source`` with respect to the scalar product ``scalar`` for all partitions of ``n`

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``n`` -- nonnegative integer which specifies the size of the partitions
        - ``source`` -- a basis of the ring of symmetric functions
        - ``scalar`` -- a function ``zee`` from partitions to the base ring which specifies the
                        scalar product by `\langle p_\lambda, p_\lambda \rangle = zee(\lambda)`.
        - ``cache`` -- a cache function
        - ``leading_coeff`` -- specifies the leading coefficients for Gram-Schmidt (default: None)
        - ``upper_triangular`` -- boolean, indicates whether the transition is upper triangular or not

        EXAMPLES::

            sage: cache = {}
            sage: from sage.combinat.sf.sfa import zee
            sage: s = SymmetricFunctions(QQ).s()
            sage: m = SymmetricFunctions(QQ).m()
            sage: s._gram_schmidt(3, m, zee, cache)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(cache)
            [([1, 1, 1], [([1, 1, 1], 1)]),
             ([2, 1], [([1, 1, 1], 2), ([2, 1], 1)]),
             ([3], [([1, 1, 1], 1), ([2, 1], 1), ([3], 1)])]
        """
        BR = self.base_ring(); one = BR(1)
        p = sage.combinat.sf.sf.SymmetricFunctions(BR).power()

        #Create a function which converts x and y to the power-sum basis and applies
        #the scalar product.
        pscalar = lambda x,y: p._apply_multi_module_morphism(p(x), p(y), lambda a,b:scalar(a), orthogonal=True)

        if leading_coeff is None:
            leading_coeff = lambda x: one

        #We are going to be doing everything like we are in the upper-triangular case
        #We list the partitions in "decreasing order" and work from the beginning forward.
        #If we are in the lower-triangular case, then we shouldn't reverse the list
        l = sage.combinat.partition.Partitions(n).list()
        if upper_triangular is True:
            l.reverse()

        #precomputed elements
        precomputed_elements = []

        #Handle the initial case
        cache[l[0]] = { l[0]: leading_coeff(l[0]) }
        precomputed_elements.append(leading_coeff( l[0] )*source(l[0]))

        for i in range(1, len(l)):
            start = leading_coeff( l[i] )*source(l[i])
            sub = 0
            for j in range(i):
                sub += pscalar( start, precomputed_elements[j] ) / pscalar(precomputed_elements[j], precomputed_elements[j])*precomputed_elements[j]
            res = start - sub

            if hasattr(self, '_normalize_coefficients'):
                res = res.map_coefficients(self._normalize_coefficients)
            precomputed_elements.append(res)
            cache[l[i]] = {}
            for j in range(i+1):
                cache[l[i]][l[j]] = res.coefficient(l[j])


    def dual_basis(self, scalar=None, scalar_name="", prefix=None):
        r"""
        Returns the dual basis of ``self`` with respect to the scalar product ``scalar``.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``scalar`` -- a function ``zee`` from partitions to the base ring which specifies the
          scalar product by `\langle p_\lambda, p_\lambda \rangle = zee(\lambda)`.
          If ``scalar`` is None, then the standard (Hall) scalar product is used.
        - ``scalar_name`` -- name of the scalar function
        - ``prefix`` -- prefix used to display the basis

        EXAMPLES: The duals of the elementary symmetric functions with
        respect to the Hall scalar product are the forgotten symmetric
        functions.

        ::

            sage: e = SymmetricFunctions(QQ).e()
            sage: f = e.dual_basis(prefix='f'); f
            Dual basis to Symmetric Function Algebra over Rational Field, Elementary symmetric functions as basis with respect to the Hall scalar product
            sage: f([2,1])^2
            4*f[2, 2, 1, 1] + 6*f[2, 2, 2] + 2*f[3, 2, 1] + 2*f[3, 3] + 2*f[4, 1, 1] + f[4, 2]
            sage: f([2,1]).scalar(e([2,1]))
            1
            sage: f([2,1]).scalar(e([1,1,1]))
            0

        Since the power-sum symmetric functions are orthogonal, their duals
        with respect to the Hall scalar product are scalar multiples of
        themselves.

        ::

            sage: p = SymmetricFunctions(QQ).p()
            sage: q = p.dual_basis(prefix='q'); q
            Dual basis to Symmetric Function Algebra over Rational Field, Power symmetric functions as basis with respect to the Hall scalar product
            sage: q([2,1])^2
            4*q[2, 2, 1, 1]
            sage: p([2,1]).scalar(q([2,1]))
            1
            sage: p([2,1]).scalar(q([1,1,1]))
            0
        """
        import dual
        if scalar is None:
            scalar = zee
            scalar_name = "Hall scalar product"
        return dual.SymmetricFunctionAlgebra_dual(self, scalar, scalar_name, prefix)


    def basis_name(self):
        r"""
        Returns the name of the basis of ``self``.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: s.basis_name()
            'schur'
            sage: p = Sym.p()
            sage: p.basis_name()
            'power'
            sage: h = Sym.h()
            sage: h.basis_name()
            'homogeneous'
            sage: e = Sym.e()
            sage: e.basis_name()
            'elementary'
            sage: m = Sym.m()
            sage: m.basis_name()
            'monomial'
        """
        try:
            return self._basis
        except AttributeError:
            return self._prefix.lower()

    def get_print_style(self):
        r"""
        Returns the value of the current print style for ``self``.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s.get_print_style()
            'lex'
            sage: s.set_print_style('length')
            sage: s.get_print_style()
            'length'
            sage: s.set_print_style('lex')
        """
        return self._print_style

    def set_print_style(self, ps):
        r"""
        Set the value of the current print style to ``ps``.

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``ps`` -- a string specifying the printing style

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s.get_print_style()
            'lex'
            sage: s.set_print_style('length')
            sage: s.get_print_style()
            'length'
            sage: s.set_print_style('lex')
        """
        if ps == 'lex':
            self.print_options(monomial_cmp = lambda x,y: cmp(x,y))
        elif ps == 'length':
            self.print_options(monomial_cmp = lambda x,y: cmp(len(x), len(y)))
        elif ps == 'maximal_part':
            self.print_options(monomial_cmp= lambda x,y: cmp(_lmax(x), _lmax(y)))
        else:
            raise ValueError, "the print style must be one of lex,length, or maximal_part "
        self._print_style = ps

    def _latex_term(self, m):
        r"""
        Latex terms (i.e. partitions) as plain lists (and not as ferrers diagrams)

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``m`` -- a partition or list

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: m._latex_term(Partition([3,2,1]))
            'm_{3,2,1}'
            sage: f = sum([m(p) for p in Partitions(3)])
            sage: m.set_print_style('lex')
            sage: latex(f)
            m_{1,1,1} + m_{2,1} + m_{3}
            sage: m.set_print_style('length')
            sage: latex(f)
            m_{3} + m_{2,1} + m_{1,1,1}
            sage: m.set_print_style('maximal_part')
            sage: latex(f)
            m_{1,1,1} + m_{2,1} + m_{3}
        """
        return super(SymmetricFunctionAlgebra_generic, self)._latex_term(','.join(str(i) for i in m))

    def from_polynomial(self, poly, check=True):
        r"""
        Convert polynomial to a symmetric function in the monomial basis and then to the basis ``self``

        INPUT:

        - ``self`` -- a basis of the ring of symmetric functions
        - ``poly`` -- a symmetric polynomial
        - ``check`` -- boolean, specifies whether the computation checks that the polynomial is indeed symmetric (default: True)

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: h = Sym.homogeneous()
            sage: f = (h([]) + h([2,1]) + h([3])).expand(3)
            sage: h.from_polynomial(f)
            h[] + h[2, 1] + h[3]
            sage: s = Sym.s()
            sage: g = (s([]) + s([2,1])).expand(3); g
            x0^2*x1 + x0*x1^2 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2 + x0*x2^2 + x1*x2^2 + 1
            sage: s.from_polynomial(g)
            s[] + s[2, 1]
        """
        m = sage.combinat.sf.sf.SymmetricFunctions(self.base_ring()).m()
        return self(m.from_polynomial(poly, check=check))

class SymmetricFunctionAlgebra_generic_Element(CombinatorialFreeModule.Element):
    r"""
    Class of generic elements for the symmetric function algebra.

    TESTS::

        sage: m = SymmetricFunctions(QQ).m()
        sage: f = sum([m(p) for p in Partitions(3)])
        sage: m.set_print_style('lex')
        sage: f
        m[1, 1, 1] + m[2, 1] + m[3]
        sage: m.set_print_style('length')
        sage: f
        m[3] + m[2, 1] + m[1, 1, 1]
        sage: m.set_print_style('maximal_part')
        sage: f
        m[1, 1, 1] + m[2, 1] + m[3]
        sage: m.set_print_style('lex')
    """

    def plethysm(self, x, include=None, exclude=None):
        r"""
        Returns the outer plethysm of ``self`` with ``x``.

        By default, the degree one elements are the generators for the
        self's base ring.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        -  ``x`` -- a symmetric function

        -  ``include`` -- a list of variables to be treated as
           degree one elements instead of the default degree one elements

        -  ``exclude`` -- a list of variables to be excluded
           from the default degree one elements

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: h = Sym.h()
            sage: s ( h([3])( h([2]) ) )
            s[2, 2, 2] + s[4, 2] + s[6]
            sage: p = Sym.p()
            sage: p([3])( s([2,1]) )
            1/3*p[3, 3, 3] - 1/3*p[9]
            sage: e = Sym.e()
            sage: e([3])( e([2]) )
            e[3, 3] + e[4, 1, 1] - 2*e[4, 2] - e[5, 1] + e[6]

        ::

            sage: R.<t> = QQ[]
            sage: s = SymmetricFunctions(R).s()
            sage: a = s([3])
            sage: f = t*s([2])
            sage: a(f)
            t^3*s[2, 2, 2] + t^3*s[4, 2] + t^3*s[6]
            sage: f(a)
            t*s[4, 2] + t*s[6]
            sage: s(0).plethysm(s[1])
            0
            sage: s(1).plethysm(s[1])
            s[]
        """
        if not is_SymmetricFunction(x):
            raise TypeError, "only know how to compute plethysms between symmetric functions"
        parent = self.parent()
        p = parent.realization_of().power()
        R = parent.base_ring()
        p_x = p(x)
        if self == parent.zero():
            return self

        #Handle degree one elements
        if include is not None and exclude is not None:
            raise RuntimeError, "include and exclude cannot both be specified"

        try:
            degree_one = [R(g) for g in R.variable_names_recursive()]
        except AttributeError:
            try:
                degree_one = R.gens()
            except NotImplementedError:
                degree_one= []

        if include:
            degree_one = [R(g) for g in include]
        if exclude:
            degree_one = [g for g in degree_one if g not in exclude]


        #Takes in n, and returns a function which takes in a partition and
        #scales all of the parts of that partition by n
        scale_part = lambda n: lambda m: m.__class__([i*n for i in m])

        raise_c = lambda n: lambda c: c.subs(**dict((str(g),g**n) for g in degree_one if g != 1))

        #Takes n an symmetric function f, and an n and returns the
        #symmetric function with all of its basis partitions scaled
        #by n
        pn_pleth = lambda f, n: f.map_support(scale_part(n))

        #Takes in a partition and applies
        f = lambda part: prod( pn_pleth(p_x.map_coefficients(raise_c(i)), i) for i in part )
        return parent(p._apply_module_morphism(p(self),f))

    __call__ = plethysm


    def _inner_plethysm_pk_g(self, k, g, cache):
        r"""
        Returns the inner plethysm between `p_k` and ``g``.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        -  ``k`` -- a positive integer

        -  ``g`` -- a symmetric function in the power sum basis

        -  ``cache`` -- a dictionary whose keys are (k, g) pairs
           and values are the cached output of this function

        EXAMPLES::

            sage: p = SymmetricFunctions(QQ).p()
            sage: _inner_plethysm_pk_g = p(0)._inner_plethysm_pk_g
            sage: _inner_plethysm_pk_g(2, p([1,1,1]), {})
            p[1, 1, 1] + 3*p[2, 1]
            sage: _inner_plethysm_pk_g(5, p([2,2,1,1,1]), {})
            p[2, 2, 1, 1, 1]
        """
        try:
            return cache[(k,g)]
        except KeyError:
            pass

        p = sage.combinat.sf.sf.SymmetricFunctions(self.base_ring()).p()
        res = 0
        degrees = uniq([ sum(m) for m in g.support() ])
        for d in degrees:
            for mu in sage.combinat.partition.Partitions(d):
                mu_k = mu.power(k)
                if mu_k in g:
                    res += g.coefficient(mu_k)*mu_k.centralizer_size()/mu.centralizer_size()*p(mu)

        cache[(k,g)] = res
        return res

    def _inner_plethysm_pnu_g(self, p_x, cache, nu):
        r"""
        Returns the inner plethysm of `p_\nu` with another symmetric function
        `p_x` in the power-sum basis.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        - ``p_x`` -- a symmetric function in the power sum basis

        - ``cache`` -- a cache function

        - ``nu`` -- a partition

        Note that the order of the arguments is somewhat strange in order
        to facilitate partial function application.

        EXAMPLES::

            sage: p = SymmetricFunctions(QQ).p()
            sage: s = SymmetricFunctions(QQ).s()
            sage: _inner_plethysm_pnu_g = p(0)._inner_plethysm_pnu_g
            sage: _inner_plethysm_pnu_g( p([1,1,1]), {}, Partition([2,1]))
            6*p[1, 1, 1]
            sage: _inner_plethysm_pnu_g( p([1,1,1]), {}, Partition([]))
            1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]
            sage: s(_)
            s[3]
        """
        #We handle the constant term case separately.  It should be
        #the case that p([]).inner_tensor(s(mu)) = s([ mu.size() ]).
        #Here, we get the degrees of the homogeneous pieces of
        if len(nu) == 0:
            s = sage.combinat.sf.sf.SymmetricFunctions(self.base_ring()).s()
            p = sage.combinat.sf.sf.SymmetricFunctions(self.base_ring()).p()
            degrees = [ part.size() for part in p_x.support() ]
            degrees = uniq(degrees)
            if 0 in degrees:
                ext = p([])
            else:
                ext = 0
            return ext + p(sum([s([n]) for n in degrees if n!=0]))

        #For each k in nu, we compute the inner plethysm of
        #p_k with p_x
        res = [self._inner_plethysm_pk_g(k, p_x, cache) for k in nu]

        #To get the final answer, we compute the inner tensor product
        #of all the symmetric functions in res
        return reduce(lambda x, y: x.itensor(y), res)

    def inner_plethysm(self, x):
        r"""
        Returns the inner plethysm of ``self`` with ``x``.

        The result of ``f.inner_plethysm(g)`` is linear in `f` and linear in
        'homogeneous pieces' of `g`. So, to describe this function, we assume
        without loss that `f` is some Schur function `s_\lambda` and `g` is a
        homogeneous symmetric function of degree `n`. The function `g` can be
        thought of as the character of an irreducible representation, rho,
        of the symmetric group `S_n`. Let `N` be the dimension of
        this representation. If the number of parts of `\lambda` is greater than
        `N`, then ``f.inner_plethysm(g)`` `= 0` by definition. Otherwise, we can
        interpret `f` as the character of an irreducible `GL_N`
        representation, call it `\sigma`. Now
        `\sigma \circ \rho` is an `S_n` representation
        and, by definition, the character of this representation is
        ``f.inner_plethysm(g)``.

        REFERENCES:

            .. [King] King, R.
               Branching rules for `GL_m \supset \Sigma_n` and the evaluation of inner plethysms.
               J. Math. Phys. 15, 258 (1974)

        INPUT:

        - ``self``, ``x`` -- elements of the ring of symmetric functions

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.schur()
            sage: p = Sym.power()
            sage: h = Sym.complete()
            sage: s([2,1]).inner_plethysm(s([1,1,1]))
            0
            sage: s([2]).inner_plethysm(s([2,1]))
            s[2, 1] + s[3]
            sage: s([1,1]).inner_plethysm(s([2,1]))
            s[1, 1, 1]
            sage: s[2,1].inner_tensor(s[2,1])
            s[1, 1, 1] + s[2, 1] + s[3]
            sage: s(0).inner_plethysm(s(0))
            0
            sage: s(1).inner_plethysm(s(0))
            0
            sage: s(0).inner_plethysm(s(1))
            0
            sage: s(1).inner_plethysm(s(1))
            s[]

        ::

            sage: f = s([2,1]) + 2*s([3,1])
            sage: f.itensor(f)
            s[1, 1, 1] + s[2, 1] + 4*s[2, 1, 1] + 4*s[2, 2] + s[3] + 4*s[3, 1] + 4*s[4]
            sage: s( h([1,1]).inner_plethysm(f) )
            s[1, 1, 1] + s[2, 1] + 4*s[2, 1, 1] + 4*s[2, 2] + s[3] + 4*s[3, 1] + 4*s[4]

        ::

            sage: s([]).inner_plethysm(s([1,1]) + 2*s([2,1])+s([3]))
            s[2] + s[3]
            sage: [s([]).inner_plethysm(s(p)) for p in Partitions(4)]
            [s[4], s[4], s[4], s[4], s[4]]
        """
        parent = self.parent()
        if self == parent.zero():
            return self
        p = parent.realization_of().power()

        p_x = p(x)

        cache = {}
        f = partial(self._inner_plethysm_pnu_g, p_x, cache)

        return parent(p._apply_module_morphism(p(self), f))

    def omega(self):
        r"""
        Returns the image of ``self`` under the Frobenius / omega automorphism.

        The default implementation converts to the Schurs performs the
        automorphism and changes back.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        EXAMPLES::

            sage: J = SymmetricFunctions(QQ).jack(t=1).P()
            sage: a = J([2,1]) + J([1,1,1])
            sage: a.omega()
            JackP[2, 1] + JackP[3]
            sage: J(0).omega()
            0
            sage: J(1).omega()
            JackP[]
        """
        parent = self.parent()
        s = parent.realization_of().schur()
        return parent(s(self).omega())

    def theta(self,a):
        r"""
        Returns the image of ``self`` under the theta automorphism which sends
        `p[k]` to `a*p[k]`.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions
        - ``a`` -- an element of the base ring

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s([2,1]).theta(2)
            2*s[1, 1, 1] + 6*s[2, 1] + 2*s[3]
            sage: p = SymmetricFunctions(QQ).p()
            sage: p([2]).theta(2)
            2*p[2]
            sage: p(0).theta(2)
            0
            sage: p(1).theta(2)
            p[]
        """
        p = self.parent().realization_of().power()
        p_self = p(self)
        res = p_self.map_item(lambda m,c: (m, c*a**len(m)))
        return self.parent()(res)

    def theta_qt(self,q=None,t=None):
        r"""
        Returns the image of ``self`` under the theta automorphism which sends
        `p[k]` to `(1-q^k)/(1-t^k)*p[k]`.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions
        - ``q``, ``t`` -- parameters (default: None, in which case 'q' and 't' are used)

        EXAMPLES::

            sage: QQqt = QQ['q,t'].fraction_field()
            sage: q,t = QQqt.gens()
            sage: p = SymmetricFunctions(QQqt).p()
            sage: p([2]).theta_qt(q,t)
            ((-q^2+1)/(-t^2+1))*p[2]
            sage: p([2,1]).theta_qt(q,t)
            ((q^3-q^2-q+1)/(t^3-t^2-t+1))*p[2, 1]
            sage: p(0).theta_qt(q=1,t=3)
            0
            sage: p([2,1]).theta_qt(q=2,t=3)
            3/16*p[2, 1]
            sage: s = p.realization_of().schur()
            sage: s([3]).theta_qt(q=0)*(1-t)*(1-t^2)*(1-t^3)
            t^3*s[1, 1, 1] + (t^2+t)*s[2, 1] + s[3]
            sage: p(1).theta_qt()
            p[]
        """
        parent = self.parent()
        BR = parent.base_ring()
        p = parent.realization_of().power()
        p_self = p(self)
        if t is None:
            if hasattr(parent,"t"):
                t = parent.t
            else:
                t = BR(QQ['t'].gen())
        if q is None:
            if hasattr(parent,"q"):
                q = parent.q
            else:
                q = BR(QQ['q'].gen())
        res = p_self.map_item(lambda m,c: (m, BR(prod([(1-q**k)/(1-t**k) for k in m])*c)))
        return parent(res)

    def omega_qt(self,q = None,t = None):
        r"""
        Returns the image of ``self`` under the theta automorphism which sends
        `p[k]` to `(-1)^{k-1}*(1-q^k)/(1-t^k)*p[k]`.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions
        - ``q``, ``t`` -- parameters (default: None, in which case 'q' and 't' are used)

        EXAMPLES::

            sage: QQqt = QQ['q,t'].fraction_field()
            sage: q,t = QQqt.gens()
            sage: p = SymmetricFunctions(QQqt).p()
            sage: p[5].omega_qt()
            ((-q^5+1)/(-t^5+1))*p[5]
            sage: p[5].omega_qt(q,t)
            ((-q^5+1)/(-t^5+1))*p[5]
            sage: p([2]).omega_qt(q,t)
            ((q^2-1)/(-t^2+1))*p[2]
            sage: p([2,1]).omega_qt(q,t)
            ((-q^3+q^2+q-1)/(t^3-t^2-t+1))*p[2, 1]
            sage: p([3,2]).omega_qt(5,q)
            ((-2976)/(q^5-q^3-q^2+1))*p[3, 2]
            sage: p(0).omega_qt()
            0
            sage: p(1).omega_qt()
            p[]
            sage: H = SymmetricFunctions(QQqt).macdonald().H()
            sage: H([1,1]).omega_qt()
            ((2*q^2-2*q*t-2*q+2*t)/(t^3-t^2-t+1))*McdH[1, 1] + ((q-1)/(t-1))*McdH[2]
            sage: H([1,1]).omega_qt(q,t)
            ((2*q^2-2*q*t-2*q+2*t)/(t^3-t^2-t+1))*McdH[1, 1] + ((q-1)/(t-1))*McdH[2]
            sage: H([1,1]).omega_qt(t,q)
            ((t^3-t^2-t+1)/(q^3-q^2-q+1))*McdH[2]
            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: S = Sym.macdonald().S()
            sage: S([1,1]).omega_qt()
            ((q^2-q*t-q+t)/(t^3-t^2-t+1))*McdS[1, 1] + ((-q^2*t+q*t+q-1)/(-t^3+t^2+t-1))*McdS[2]
            sage: s = Sym.schur()
            sage: s(S([1,1]).omega_qt())
            s[2]
        """
        parent = self.parent()
        BR = parent.base_ring()
        p = parent.realization_of().power()
        p_self = p(self)
        if t is None:
            if hasattr(parent,"t"):
                t = parent.t
            else:
                t = BR(QQ['t'].gen())
        if q is None:
            if hasattr(parent,"q"):
                q = parent.q
            else:
                q = BR(QQ['q'].gen())
        f = lambda part: prod([(-1)**(i-1)*(1-q**i)/(1-t**i) for i in part])
        res = p_self.map_item(lambda m,c: (m, BR(f(m)*c)))
        return parent(res)

    def itensor(self, x):
        r"""
        Returns the inner tensor product of ``self`` and ``x`` in the basis of ``self``.

        The inner tensor product (also known as the Kronecker product) can be defined as the linear extension of the
        definition on power sums `p_\lambda \otimes p_\mu = \delta_{\lambda,\mu} z_\lambda p_\lambda`, where
        `z_\lambda = (1^{r_1} r_1!) (2^{r_2} r_2!) \cdots` for `\lambda = (1^{r_1} 2^{r_2} \cdots )`.

        INPUT:

        - ``self``, ``x`` -- elements of the ring of symmetric functions of the same degree

        OUTPUT:

        - symmetric function of the same degree as ``self`` and ``x``

        The methods :meth:`itensor`, :meth:`kronecker_product`, :meth:`inner_tensor` are all
        synonyms.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([2,1])
            sage: b = s([3])
            sage: a.itensor(b)
            s[2, 1]
            sage: c = s([3,2,1])
            sage: c.itensor(c)
            s[1, 1, 1, 1, 1, 1] + 2*s[2, 1, 1, 1, 1] + 3*s[2, 2, 1, 1] + 2*s[2, 2, 2] + 4*s[3, 1, 1, 1] + 5*s[3, 2, 1] + 2*s[3, 3] + 4*s[4, 1, 1] + 3*s[4, 2] + 2*s[5, 1] + s[6]

        TESTS::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([8,8])
            sage: a.itensor(a)
            s[4, 4, 4, 4] + s[5, 5, 3, 3] + s[5, 5, 5, 1] + s[6, 4, 4, 2] + s[6, 6, 2, 2] + s[6, 6, 4] + s[7, 3, 3, 3] + s[7, 5, 3, 1] + s[7, 7, 1, 1] + s[8, 4, 2, 2] + s[8, 4, 4] + s[8, 6, 2] + s[8, 8] + s[9, 3, 3, 1] + s[9, 5, 1, 1] + s[10, 2, 2, 2] + s[10, 4, 2] + s[10, 6] + s[11, 3, 1, 1] + s[12, 2, 2] + s[12, 4] + s[13, 1, 1, 1] + s[14, 2] + s[16]
            sage: s[8].itensor(s[7])
            0
            sage: s(0).itensor(s(0))
            0
            sage: s(1).itensor(s(0))
            0
            sage: s(0).itensor(s(1))
            0
            sage: s(1).itensor(s(1))
            s[]
        """
        #Convert both self and x to the p basis
        parent = self.parent()
        p = parent.realization_of().power()
        f = lambda part1, part2: zee(part1)*p(part1)
        return parent(p._apply_multi_module_morphism(p(self),p(x),f,orthogonal=True))


    internal_product = itensor
    kronecker_product = itensor
    inner_tensor = itensor

    def nabla(self, q=None, t=None, power=1):
        r"""
        Returns the value of the nabla operator applied to ``self``.

        The eigenvectors of the nabla operator are the Macdonald polynomials in
        the Ht basis.

        If the parameter ``power`` is an integer then it calculates
        nabla to that integer.  The default value of ``power`` is 1.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions
        - ``q``, ``t`` -- optional parameters (default: None in which case `q` and `t` are used)
        - ``power`` -- an integer indicating how many times to apply the
          operator `\nabla` (default : 1).  Negative values of ``power`` indicate
          powers of `\nabla^{-1}`.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
            sage: p = Sym.power()
            sage: p([1,1]).nabla()
            (-1/2*q*t+1/2*q+1/2*t+1/2)*p[1, 1] + (1/2*q*t-1/2*q-1/2*t+1/2)*p[2]
            sage: p([2,1]).nabla(q=1)
            (-t-1)*p[1, 1, 1] + t*p[2, 1]
            sage: p([2]).nabla(q=1)*p([1]).nabla(q=1)
            (-t-1)*p[1, 1, 1] + t*p[2, 1]
            sage: s = Sym.schur()
            sage: s([2,1]).nabla()
            (-q^3*t-q^2*t^2-q*t^3)*s[1, 1, 1] + (-q^2*t-q*t^2)*s[2, 1]
            sage: s([1,1,1]).nabla()
            (q^3+q^2*t+q*t^2+t^3+q*t)*s[1, 1, 1] + (q^2+q*t+t^2+q+t)*s[2, 1] + s[3]
            sage: s([1,1,1]).nabla(t=1)
            (q^3+q^2+2*q+1)*s[1, 1, 1] + (q^2+2*q+2)*s[2, 1] + s[3]
            sage: s(0).nabla()
            0
            sage: s(1).nabla()
            s[]
            sage: s([2,1]).nabla(power=-1)
            ((-q-t)/(q^2*t^2))*s[2, 1] + ((-q^2-q*t-t^2)/(q^3*t^3))*s[3]
            sage: (s([2])+s([3])).nabla()
            (-q*t)*s[1, 1] + (q^3*t^2+q^2*t^3)*s[1, 1, 1] + q^2*t^2*s[2, 1]
        """
        parent = self.parent()
        BR = parent.base_ring()
        if q is None:
            if hasattr(parent,"q"):
                q = parent.q
            else:
                q = BR(QQ['q'].gen())
        if t is None:
            if hasattr(parent,"t"):
                t = parent.t
            else:
                t = BR(QQ['t'].gen())
        Ht = parent.realization_of().macdonald(q=q,t=t).Ht()
        return parent(Ht(self).nabla(power=power))

    def scalar(self, x, zee=None):
        r"""
        Returns standard scalar product between ``self`` and ``x``.

        INPUT:

        - ``self``, ``x`` -- elements of the ring of symmetric functions

        - ``zee`` -- an optional function on partitions giving
          the value for the scalar product between `p_\mu` and `p_\mu`
          (default is to use the standard zee function)

        This is the default implementation that converts both ``self`` and ``x``
        into Schur functions and performs the scalar product that basis.

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).e()
            sage: h = SymmetricFunctions(QQ).h()
            sage: m = SymmetricFunctions(QQ).m()
            sage: p4 = Partitions(4)
            sage: matrix([ [e(a).scalar(h(b)) for a in p4] for b in p4])
            [ 0  0  0  0  1]
            [ 0  0  0  1  4]
            [ 0  0  1  2  6]
            [ 0  1  2  5 12]
            [ 1  4  6 12 24]
            sage: matrix([ [h(a).scalar(e(b)) for a in p4] for b in p4])
            [ 0  0  0  0  1]
            [ 0  0  0  1  4]
            [ 0  0  1  2  6]
            [ 0  1  2  5 12]
            [ 1  4  6 12 24]
            sage: matrix([ [m(a).scalar(e(b)) for a in p4] for b in p4])
            [-1  2  1 -3  1]
            [ 0  1  0 -2  1]
            [ 0  0  1 -2  1]
            [ 0  0  0 -1  1]
            [ 0  0  0  0  1]
            sage: matrix([ [m(a).scalar(h(b)) for a in p4] for b in p4])
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]

            sage: p = SymmetricFunctions(QQ).p()
            sage: m(p[3,2]).scalar(p[3,2], zee=lambda mu: 2**mu.length())
            4
            sage: m(p[3,2]).scalar(p[2,2,1], lambda mu: 1)
            0
            sage: m[3,2].scalar(h[3,2], zee=lambda mu: 2**mu.length())
            2/3

      TESTS::

            sage: m(1).scalar(h(1))
            1
            sage: m(0).scalar(h(1))
            0
            sage: m(1).scalar(h(0))
            0
            sage: m(0).scalar(h(0))
            0
        """
        if zee is None:
            s = self.parent().realization_of().schur()
            s_self = s(self)
            s_x = s(x)
            return s_self.scalar(s_x)
        else:
            p = self.parent().realization_of().power()
            p_self = p(self)
            p_x = p(x)
            return sum(zee(mu)*p_x.coefficient(mu)*p_self.coefficient(mu) for mu in p_self.support())

    def scalar_qt(self, x, q = None, t = None):
        r"""
        Returns the standard Hall-Littlewood scalar product of ``self`` and ``x``.

        INPUT:

        - ``self``, ``x`` -- elements of the ring of symmetric functions

        - ``q``, ``t`` -- parameters (default: None in which case `q` and `t` are used)

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([2,1])
            sage: sp = a.scalar_qt(a); factor(sp)
            (t - 1)^-3 * (q - 1) * (t^2 + t + 1)^-1 * (q^2*t^2 - q*t^2 + q^2 - 2*q*t + t^2 - q + 1)
            sage: sp.parent()
            Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
            sage: a.scalar_qt(a,q=0)
            (-t^2 - 1)/(t^5 - 2*t^4 + t^3 - t^2 + 2*t - 1)
            sage: a.scalar_qt(a,t=0)
            -q^3 + 2*q^2 - 2*q + 1
            sage: a.scalar_qt(a,5,7) # q=5 and t=7
            490/1539
            sage: (x,y) = var('x,y')
            sage: a.scalar_qt(a,q=x,t=y)
            2/3*(x - 1)^3/(y - 1)^3 + 1/3*(x^3 - 1)/(y^3 - 1)
            sage: Rn = QQ['q','t','y','z'].fraction_field()
            sage: (q,t,y,z) = Rn.gens()
            sage: Mac = SymmetricFunctions(Rn).macdonald(q=y,t=z)
            sage: a = Mac._sym.schur()([2,1])
            sage: factor(Mac.P()(a).scalar_qt(Mac.Q()(a),q,t))
            (t - 1)^-3 * (q - 1) * (t^2 + t + 1)^-1 * (q^2*t^2 - q*t^2 + q^2 - 2*q*t + t^2 - q + 1)
            sage: factor(Mac.P()(a).scalar_qt(Mac.Q()(a)))
            (z - 1)^-3 * (y - 1) * (z^2 + z + 1)^-1 * (y^2*z^2 - y*z^2 + y^2 - 2*y*z + z^2 - y + 1)
        """
        parent = self.parent()
        p = parent.realization_of().power()
        if t is None:
            if hasattr(parent,"t"):
                t = self.parent().t
            else:
                if q is None:
                    t = QQ['q','t'].gens()[1]
                else:
                    t = QQ['t'].gen()
        if q is None:
            if hasattr(parent,"q"):
                q = parent.q
            else:
                q = QQ['q','t'].gens()[0]
        f = lambda part1, part2: part1.centralizer_size(t = t, q = q)
        return p._apply_multi_module_morphism(p(self), p(x), f, orthogonal=True)

    def scalar_t(self, x, t = None):
        r"""
        Returns the standard Hall-Littlewood scalar product of ``self`` and ``x``.

        INPUT:

        - ``self``, ``x`` -- elements of the ring of symmetric functions

        - ``t`` -- parameter (default: None in which case `q` and `t` are used)

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([2,1])
            sage: sp = a.scalar_t(a); sp
            (-t^2 - 1)/(t^5 - 2*t^4 + t^3 - t^2 + 2*t - 1)
            sage: sp.parent()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self.scalar_qt( x, q=self.base_ring()(0), t=t )

    scalar_hl = scalar_t

    def scalar_jack(self, x, t=None):
        r"""
        Returns the Jack-scalar product beween ``self`` and ``x``

        This scalar product is defined so that the power sum elements
        `p_\mu` are orthogonal and `\langle p_\mu, p_\mu \rangle = z_\mu t^{length(\mu)}`

        INPUT:

        - ``self``, ``x`` -- elements of the ring of symmetric functions
        - ``t`` -- an optional parameter (default: None in which case `t` is used)

        EXAMPLES ::

            sage: p = SymmetricFunctions(QQ['t']).power()
            sage: matrix([[p(mu).scalar_jack(p(nu)) for nu in Partitions(4)] for mu in Partitions(4)])
            [   4*t      0      0      0      0]
            [     0  3*t^2      0      0      0]
            [     0      0  8*t^2      0      0]
            [     0      0      0  4*t^3      0]
            [     0      0      0      0 24*t^4]
            sage: matrix([[p(mu).scalar_jack(p(nu),2) for nu in Partitions(4)] for mu in Partitions(4)])
            [  8   0   0   0   0]
            [  0  12   0   0   0]
            [  0   0  32   0   0]
            [  0   0   0  32   0]
            [  0   0   0   0 384]
            sage: JQ = SymmetricFunctions(QQ['t'].fraction_field()).jack().Q()
            sage: matrix([[JQ(mu).scalar_jack(JQ(nu)) for nu in Partitions(3)] for mu in Partitions(3)])
            [(2*t^2 + 3*t + 1)/(6*t^3)                         0                         0]
            [                        0     (t + 2)/(2*t^3 + t^2)                         0]
            [                        0                         0     6/(t^3 + 3*t^2 + 2*t)]
        """
        parent = self.parent()
        p = parent.realization_of().power()
        if t is None:
            if hasattr(parent,"t"):
                t = self.parent().t
            else:
                t = QQ['t'].gen()
        zee = lambda part: part.centralizer_size()*t**part.length()
        return self.scalar(x, zee)

    def derivative_with_respect_to_p1(self, n=1):
        r"""
        Returns the symmetric function obtained by taking the derivative of
        ``self`` with respect to the power-sum symmetric function `p([1])` when
        the expansion of ``self`` in the power-sum basis is considered as a
        polynomial in `p([1])`'s.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions
        - ``n`` -- positive integer which determines which power of the derivative is taken (default: `n=1`)

        EXAMPLES::

            sage: p = SymmetricFunctions(QQ).p()
            sage: a = p([1,1,1])
            sage: a.derivative_with_respect_to_p1()
            3*p[1, 1]
            sage: a.derivative_with_respect_to_p1(1)
            3*p[1, 1]
            sage: a.derivative_with_respect_to_p1(2)
            6*p[1]
            sage: a.derivative_with_respect_to_p1(3)
            6*p[]

        ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s([3]).derivative_with_respect_to_p1()
            s[2]
            sage: s([2,1]).derivative_with_respect_to_p1()
            s[1, 1] + s[2]
            sage: s([1,1,1]).derivative_with_respect_to_p1()
            s[1, 1]
            sage: s(0).derivative_with_respect_to_p1()
            0
            sage: s(1).derivative_with_respect_to_p1()
            0
            sage: s([1]).derivative_with_respect_to_p1()
            s[]
        """
        p = self.parent().realization_of().power()
        res = p(self)
        for i in range(n):
            res = res._derivative_with_respect_to_p1()
        return self.parent()(res)


    def _expand(self, condition, n, alphabet = 'x'):
        r"""
        Expands the symmetric function as a symmetric polynomial in ``n`` variables.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        - ``condition`` -- a function on partitions with a boolean output, selecting only certain terms

        - ``n`` -- a positive integer

        - ``alphabet`` -- a variable for the expansion (default: `x`)

        OUTPUT: a monomial expansion of an instance of ``self`` in `n` variables

        EXAMPLES::

            sage: p = SymmetricFunctions(QQ).p()
            sage: a = p([2])+p([3])
            sage: a._expand(lambda part: False, 3)
            x0^3 + x1^3 + x2^3 + x0^2 + x1^2 + x2^2
            sage: a._expand(lambda part: max(part)>2, 3)
            x0^2 + x1^2 + x2^2
            sage: p(0).expand(3)
            0
            sage: p([]).expand(3)
            1
        """
        import classical
        parent = self.parent()
        resPR = PolynomialRing(parent.base_ring(), n, alphabet)
        if self==parent.zero():
            return resPR.zero()
        e = eval('symmetrica.compute_' + str(classical.translate[parent.basis_name()]).lower() + '_with_alphabet')
        def f(part):
            if part==[]:
                return resPR(1)
            else:
                return resPR(0) if condition(part) else resPR(e(part, n, alphabet))
        return parent._apply_module_morphism(self, f)

    def is_schur_positive(self):
        r"""
        Returns True if and only if ``self`` is Schur positive.

        If s is the space of Schur functions over ``self``'s base ring, then this is the
        same as self._is_positive(s).

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        EXAMPLES ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([2,1]) + s([3])
            sage: a.is_schur_positive()
            True
            sage: a = s([2,1]) - s([3])
            sage: a.is_schur_positive()
            False

        ::

            sage: QQx = QQ['x']
            sage: s = SymmetricFunctions(QQx).s()
            sage: x = QQx.gen()
            sage: a = (1+x)*s([2,1])
            sage: a.is_schur_positive()
            True
            sage: a = (1-x)*s([2,1])
            sage: a.is_schur_positive()
            False
            sage: s(0).is_schur_positive()
            True
            sage: s(1+x).is_schur_positive()
            True
        """
        return self._is_positive( self.parent().realization_of().schur() )


    def _is_positive(self, s):
        r"""
        Returns True if and only if ``self`` has nonnegative coefficients in the basis s.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions
        - ``s`` -- a basis of the ring of symmetric functions

        EXAMPLES ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([2,1]) + s([3])
            sage: a._is_positive(s)
            True
            sage: a = s([2,1]) - s([3])
            sage: a._is_positive(s)
            False

            sage: m = SymmetricFunctions(QQ).m()
            sage: a = s([2,1]) + s([3])
            sage: a._is_positive(m)
            True
            sage: a = -s[2,1]
            sage: a._is_positive(m)
            False

            sage: (s[2,1] - s[1,1,1])._is_positive(s)
            False
            sage: (s[2,1] - s[1,1,1])._is_positive(m)
            True
        """
        s_self = s(self)
        return all([ _nonnegative_coefficients(c) for c in s_self.coefficients()])

    def degree(self):
        r"""
        Returns the degree of ``self``

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        EXAMPLES ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1]) + 3
            sage: z.degree()
            4
            sage: s(1).degree()
            0
            sage: s(0).degree()
            0
        """
        return max( map( sum, self._monomial_coefficients ) + [0] )

    def restrict_degree(self, d, exact = True):
        r"""
        Returns the degree ``d`` terms of ``self``.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        - ``d`` -- positive integer, degree of the terms to be returned

        - ``exact`` -- boolean, if True, returns the terms of degree exactly ``d``, otherwise returns all terms of degree less than or equal to ``d``

        OUTPUT:

        - returns the homogeneous component of ``self`` of degree ``d``

        EXAMPLES ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_degree(2)
            0
            sage: z.restrict_degree(1)
            s[1]
            sage: z.restrict_degree(3)
            s[1, 1, 1] + s[2, 1]
            sage: z.restrict_degree(3, exact=False)
            s[1] + s[1, 1, 1] + s[2, 1]
            sage: z.restrict_degree(0)
            0
        """
        if exact:
            res = dict( filter( lambda x: sum(x[0]) == d, self._monomial_coefficients.items()) )
        else:
            res = dict( filter( lambda x: sum(x[0]) <= d, self._monomial_coefficients.items()) )
        return self.parent()._from_dict(res)

    def restrict_partition_lengths(self, l, exact = True):
        r"""
        Returns the terms of ``self`` labelled by partitions of length ``l``

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        - ``l`` -- positive integer, length of the partitions of the terms to be returned

        - ``exact`` -- boolean, if True, returns the terms of degree exactly ``d``, otherwise returns all terms of degree less than or equal to ``d``

        EXAMPLES ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_partition_lengths(2)
            s[2, 1]
            sage: z.restrict_partition_lengths(0)
            0
            sage: z.restrict_partition_lengths(2, exact = False)
            s[1] + s[2, 1] + s[4]
        """
        if exact:
            res = dict( filter( lambda x: len(x[0]) == l, self._monomial_coefficients.items()) )
        else:
            res = dict( filter( lambda x: len(x[0]) <= l, self._monomial_coefficients.items()) )
        return self.parent()._from_dict(res)

    def restrict_parts(self, n):
        r"""
        Returns the terms of ``self`` labelled by partitions `\lambda` with `\lambda_1 \le n`.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        - ``n`` -- positive integer, restrict the parts of the partitions of the terms to be returned

        EXAMPLES ::

            sage: s = SymmetricFunctions(QQ).s()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_parts(2)
            s[1] + s[1, 1, 1] + s[2, 1]
            sage: z.restrict_parts(1)
            s[1] + s[1, 1, 1]
        """
        res = dict( filter( lambda x: _lmax(x[0]) <= n, self._monomial_coefficients.items()) )
        return self.parent()._from_dict(res)

    def expand(self, n, alphabet = 'x'):
        r"""
        Expands the symmetric function as a symmetric polynomial in ``n`` variables.

        INPUT:

        - ``self`` -- an element of the ring of symmetric functions

        - ``n`` -- a positive integer

        - ``alphabet`` -- a variable for the expansion (default: `x`)

        OUTPUT: a monomial expansion of an instance of ``self`` in `n` variables

        EXAMPLES ::

            sage: J = SymmetricFunctions(QQ).jack(t=2).J()
            sage: J([2,1]).expand(3)
            4*x0^2*x1 + 4*x0*x1^2 + 4*x0^2*x2 + 6*x0*x1*x2 + 4*x1^2*x2 + 4*x0*x2^2 + 4*x1*x2^2
        """
        s = self.parent().realization_of().schur()
        condition = lambda part: len(part) > n
        return s(self)._expand(condition, n, alphabet)

    def skew_by(self, x):
        r"""
        Returns the element whose result is the dual to multiplication by ``x`` applied to ``self``

        INPUT:

        - ``self``, ``x`` -- elements of the ring of symmetric functions

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s([3,2]).skew_by(s([2]))
            s[2, 1] + s[3]
            sage: s([3,2]).skew_by(s([1,1,1]))
            0
            sage: s([3,2,1]).skew_by(s([2,1]))
            s[1, 1, 1] + 2*s[2, 1] + s[3]

        ::

            sage: p = SymmetricFunctions(QQ).powersum()
            sage: p([4,3,3,2,2,1]).skew_by(p([2,1]))
            4*p[4, 3, 3, 2]
            sage: zee = sage.combinat.sf.sfa.zee
            sage: zee([4,3,3,2,2,1])/zee([4,3,3,2])
            4
            sage: s(0).skew_by(s([1]))
            0
            sage: s(1).skew_by(s([1]))
            0
            sage: s([]).skew_by(s([]))
            s[]
            sage: s([]).skew_by(s[1])
            0

        TESTS::

            sage: f=s[3,2]
            sage: f.skew_by([1])
            Traceback (most recent call last):
            ...
            ValueError: x needs to be a symmetric function
        """
        if x not in self.parent().realization_of():
            raise ValueError("x needs to be a symmetric function")
        s = self.parent().realization_of().schur()
        f = lambda part1, part2: s([part1,part2]) if part1.contains(part2) else 0
        return self.parent()(s._apply_multi_module_morphism(s(self),s(x),f))

    def hl_creation_operator(self, nu, t = None):
        r"""
        This is the vertex operator that generalizes Jing's operator.

        It is a linear operator that raises the degree by
        sum(nu). This creation operator is a t-analogue of
        multiplication by s(nu)

        .. SEEALSO:: Proposition 5 in [SZ2001]_.

        INPUT:

        -  ``nu`` -- a partition

        - ``t`` -- a parameter (default: None, in this case `t` is used)

        REFERENCES:

        .. [SZ2001] M. Shimozono, M. Zabrocki,
           Hall-Littlewood vertex operators and generalized Kostka polynomials.
           Adv. Math. 158 (2001), no. 1, 66-85.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ['t']).s()
            sage: s([2]).hl_creation_operator([3,2])
            s[3, 2, 2] + t*s[3, 3, 1] + t*s[4, 2, 1] + t^2*s[4, 3] + t^2*s[5, 2]

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: HLQp = Sym.hall_littlewood().Qp()
            sage: s = Sym.s()
            sage: HLQp(s([2]).hl_creation_operator([2]).hl_creation_operator([3]))
            HLQp[3, 2, 2]
            sage: s([2,2]).hl_creation_operator([2,1])
            t*s[2, 2, 2, 1] + t^2*s[3, 2, 1, 1] + t^2*s[3, 2, 2] + t^3*s[3, 3, 1] + t^3*s[4, 2, 1] + t^4*s[4, 3]
            sage: s(1).hl_creation_operator([2,1,1])
            s[2, 1, 1]
            sage: s(0).hl_creation_operator([2,1,1])
            0
            sage: s([3,2]).hl_creation_operator([2,1,1])
            (t^2-t)*s[2, 2, 2, 2, 1] + t^3*s[3, 2, 2, 1, 1] + (t^3-t^2)*s[3, 2, 2, 2] + t^3*s[3, 3, 1, 1, 1] + t^4*s[3, 3, 2, 1] + t^3*s[4, 2, 1, 1, 1] + t^4*s[4, 2, 2, 1] + 2*t^4*s[4, 3, 1, 1] + t^5*s[4, 3, 2] + t^5*s[4, 4, 1] + t^4*s[5, 2, 1, 1] + t^5*s[5, 3, 1]

        TESTS::

            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: ks = Sym.kschur(4)
            sage: ks([3,1,1]).lift().hl_creation_operator([1])
            (t-1)*s[2, 2, 1, 1] + t^2*s[3, 1, 1, 1] + (t^3+t^2-t)*s[3, 2, 1] + (t^3-t^2)*s[3, 3] + (t^4+t^3)*s[4, 1, 1] + t^4*s[4, 2] + t^5*s[5, 1]
            sage: ks(ks([3,1,1]).lift().hl_creation_operator([1]))
            (t-1)*ks4[2, 2, 1, 1] + t^2*ks4[3, 1, 1, 1] + t^3*ks4[3, 2, 1] + (t^3-t^2)*ks4[3, 3] + t^4*ks4[4, 1, 1]
            sage: s(0).hl_creation_operator([1])
            0
        """
        s = self.parent().realization_of().schur()
        if t is None:
            if hasattr(self.parent(),"t"):
                t = self.parent().t
            else:
                t = QQ['t'].gen()
        P = self.parent()
        self = s(self)
        return P(self*s(nu) +
                 s.sum( s.sum_of_terms( (lam,c) for lam, c in s(mu)*s(nu) if len(lam) <= len(nu) ) *
                        self.skew_by(s(mu).plethysm((t-1)*s([1])))
                        for d in range(self.degree())
                        for mu in Partitions(d+1, max_length=len(nu)) )
                )

SymmetricFunctionAlgebra_generic.Element = SymmetricFunctionAlgebra_generic_Element


###################
def _lmax(x):
    r"""
    Returns the max of ``x`` where ``x`` is a list.

    If ``x`` is the empty list, ``_lmax`` returns 0.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import _lmax
        sage: _lmax([3,2,1])
        3
        sage: _lmax([])
        0
    """
    if x == []:
        return 0
    else:
        return max(x)


def _nonnegative_coefficients(x):
    r"""
    Returns True if ``x`` has nonnegative coefficients.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import _nonnegative_coefficients
        sage: _nonnegative_coefficients(2)
        True
        sage: _nonnegative_coefficients(-2)
        False
        sage: R.<x> = ZZ[]
        sage: _nonnegative_coefficients(x^2+4)
        True
        sage: _nonnegative_coefficients(x^2-4)
        False
    """
    if is_Polynomial(x) or is_MPolynomial(x):
        return all([ c >= 0 for c in x.coeffs() ])
    else:
        return x >= 0
