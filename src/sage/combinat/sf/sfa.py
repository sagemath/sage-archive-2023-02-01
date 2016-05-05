r"""
Symmetric Functions

For a comprehensive tutorial on how to use symmetric functions in Sage

.. SEEALSO:: :func:`SymmetricFunctions`

We define the algebra of symmetric functions in the Schur and elementary bases::

    sage: s = SymmetricFunctions(QQ).schur()
    sage: e = SymmetricFunctions(QQ).elementary()

Each is actually a graded Hopf algebra whose basis is indexed by
integer partitions::

    sage: s.category()
    Category of graded bases of Symmetric Functions over Rational Field
    sage: s.basis().keys()
    Partitions

Let us compute with some elements in different bases::

    sage: f1 = s([2,1]); f1
    s[2, 1]
    sage: f2 = e(f1); f2 # basis conversion
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
    sage: m(4) # This is the constant 4, not the partition 4.
    4*m[]
    sage: m([4]) # This is the partition 4.
    m[4]
    sage: 3*m([3,1])-1/2*m([4])
    3*m[3, 1] - 1/2*m[4]

::

    sage: p = SymmetricFunctions(QQ).power()
    sage: f = p(3)
    sage: f
    3*p[]
    sage: f.parent()
    Symmetric Functions over Rational Field in the powersum basis
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
    sage: g = h[1,1,1,1]
    sage: poly = g.expand(3)
    sage: Sym.from_polynomial(poly) == g
    False

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

TESTS:

Check that we can handle large integers properly (:trac:`13413`)::

    sage: s = SymmetricFunctions(QQ).s()
    sage: p = SymmetricFunctions(QQ).p()
    sage: max(s(p([1]*36)).coefficients())  # long time (4s on sage.math, 2013)
    40971642983700000000

BACKWARD INCOMPATIBLE CHANGES (:trac:`5457`):

The symmetric functions code has been refactored to take
advantage of the coercion systems. This introduced a couple of glitches:

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
- Darij Grinberg (2013) Sym over rings that are not characteristic 0

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
from sage.rings.all import Integer, PolynomialRing, QQ, ZZ
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.rings.polynomial.multi_polynomial import is_MPolynomial
from sage.combinat.partition import _Partitions, Partitions, Partitions_n, Partition
from sage.categories.algebras import Algebras
from sage.categories.hopf_algebras import HopfAlgebras
import sage.libs.symmetrica.all as symmetrica  # used in eval()
from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.constructor import matrix
from sage.misc.all import prod, uniq
from copy import copy
from functools import reduce


def is_SymmetricFunctionAlgebra(x):
    """
    Checks whether ``x`` is a symmetric function algebra.

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import is_SymmetricFunctionAlgebra
        sage: is_SymmetricFunctionAlgebra(5)
        False
        sage: is_SymmetricFunctionAlgebra(ZZ)
        False
        sage: is_SymmetricFunctionAlgebra(SymmetricFunctions(ZZ).schur())
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
    Return the size of the centralizer of any permutation of cycle type
    ``part``.

    Note that the size of the centralizer is the inner product between
    ``p(part)`` and itself, where `p` is the power-sum symmetric
    functions.

    INPUT:

    - ``part`` -- an integer partition (for example, ``[2,1,1]``)

    OUTPUT:

    - the integer `\prod_{i} i^{m_i(part)} m_i(part)!` where `m_i(part)` is
      the number of parts in the partition ``part`` equal to `i`

    EXAMPLES::

        sage: from sage.combinat.sf.sfa import zee
        sage: zee([2,1,1])
        4
    """
    if not isinstance(part, Partition):
        part = _Partitions(part)
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

#####################################################################
## Bases categories

from sage.categories.realizations import Category_realization_of_parent
class SymmetricFunctionsBases(Category_realization_of_parent):
    r"""
    The category of bases of the ring of symmetric functions.

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
    def _repr_(self):
        r"""
        Return the representation of ``self``.

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

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import SymmetricFunctionsBases
            sage: Sym = SymmetricFunctions(QQ)
            sage: bases = SymmetricFunctionsBases(Sym)
            sage: bases.super_categories()
            [Category of realizations of Symmetric Functions over Rational Field,
             Category of commutative hopf algebras with basis over Rational Field,
             Join of Category of realizations of hopf algebras over Rational Field
                 and Category of graded algebras over Rational Field]
        """
        # FIXME: The last one should also be commutative, but this triggers a
        #   KeyError when doing the C3 algorithm!!!
        cat = HopfAlgebras(self.base().base_ring())
        return [self.base().Realizations(),
                cat.Commutative().WithBasis(),
                cat.Graded().Realizations()]

    class ParentMethods:

        def is_integral_domain(self, proof=True):
            """
            Return whether ``self`` is an integral domain. (It is if
            and only if the base ring is an integral domain.)

            INPUT:

            - ``self`` -- a basis of the symmetric functions
            - ``proof`` -- an optional argument (default value: ``True``)

            EXAMPLES::

                sage: s = SymmetricFunctions(QQ).s()
                sage: s.is_integral_domain()
                True

            The following doctest is disabled pending :trac:`15475`::

                sage: s = SymmetricFunctions(Zmod(14)).s() # not tested
                sage: s.is_integral_domain() # not tested
                False
            """
            return self.base_ring().is_integral_domain()

        def is_field(self, proof=True):
            """
            Return whether ``self`` is a field. (It is not.)

            INPUT:

            - ``self`` -- a basis of the symmetric functions
            - ``proof`` -- an optional argument (default value: ``True``)

            EXAMPLES::

                sage: s = SymmetricFunctions(QQ).s()
                sage: s.is_field()
                False
            """
            return False

        def is_commutative(self):
            """
            Returns whether this symmetric function algebra is commutative.

            INPUT:

            - ``self`` -- a basis of the symmetric functions

            EXAMPLES::

                sage: s = SymmetricFunctions(QQ).s()
                sage: s.is_commutative()
                True
            """
            return self.base_ring().is_commutative()

        def _repr_(self):
            """
            Text representation of this basis of symmetric functions

            INPUT:

            - ``self`` -- a basis of the symmetric functions

            EXAMPLES::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q,t'])); Sym
                Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
                sage: Sym.p()
                Symmetric Functions over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field in the powersum basis

            In the following examples, we rename {{{Sym}}} for brevity::

                sage: Sym.rename("Sym"); Sym
                Sym

            Classical bases::

                sage: Sym.s()
                Sym in the Schur basis
                sage: Sym.p()
                Sym in the powersum basis
                sage: Sym.m()
                Sym in the monomial basis
                sage: Sym.e()
                Sym in the elementary basis
                sage: Sym.h()
                Sym in the homogeneous basis
                sage: Sym.f()
                Sym in the forgotten basis

            Macdonald polynomials::

                sage: Sym.macdonald().P()
                Sym in the Macdonald P basis
                sage: Sym.macdonald().Q()
                Sym in the Macdonald Q basis
                sage: Sym.macdonald().J()
                Sym in the Macdonald J basis
                sage: Sym.macdonald().H()
                Sym in the Macdonald H basis
                sage: Sym.macdonald().Ht()
                Sym in the Macdonald Ht basis
                sage: Sym.macdonald().S()
                Sym in the Macdonald S basis

            Macdonald polynomials, with specialized parameters::

                sage: Sym.macdonald(q=1).S()
                Sym in the Macdonald S with q=1 basis
                sage: Sym.macdonald(q=1,t=3).P()
                Sym in the Macdonald P with q=1 and t=3 basis

            Hall-Littlewood polynomials:

                sage: Sym.hall_littlewood().P()
                Sym in the Hall-Littlewood P basis
                sage: Sym.hall_littlewood().Q()
                Sym in the Hall-Littlewood Q basis
                sage: Sym.hall_littlewood().Qp()
                Sym in the Hall-Littlewood Qp basis

            Hall-Littlewood polynomials, with specialized parameter::

                sage: Sym.hall_littlewood(t=1).P()
                Sym in the Hall-Littlewood P with t=1 basis

            Jack polynomials::

                sage: Sym.jack().J()
                Sym in the Jack J basis
                sage: Sym.jack().P()
                Sym in the Jack P basis
                sage: Sym.jack().Q()
                Sym in the Jack Q basis
                sage: Sym.jack().Qp()
                Sym in the Jack Qp basis

            Jack polynomials, with specialized parameter::

                sage: Sym.jack(t=1).J()
                Sym in the Jack J with t=1 basis

            Zonal polynomials::

                sage: Sym.zonal()
                Sym in the zonal basis

            LLT polynomials::

                sage: Sym.llt(3).hspin()
                Sym in the level 3 LLT spin basis
                sage: Sym.llt(3).hcospin()
                Sym in the level 3 LLT cospin basis

            LLT polynomials, with specialized parameter::

                sage: Sym.llt(3, t=1).hspin()
                Sym in the level 3 LLT spin with t=1 basis
                sage: Sym.llt(3, t=1).hcospin()
                Sym in the level 3 LLT cospin with t=1 basis

            TESTS::

                sage: Sym.s()._repr_()
                'Sym in the Schur basis'
                sage: Sym.s()._repr_.__module__
                'sage.combinat.sf.sfa'

            ::

                sage: Sym.rename()
            """
            return "%s in the %s basis"%(self.realization_of(), self.basis_name())

        @cached_method
        def one_basis(self):
            r"""
            Returns the empty partition, as per ``AlgebrasWithBasis.ParentMethods.one_basis``

            INPUT:

            - ``self`` -- a basis of the ring of symmetric functions

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
                sage: s = Sym.s()
                sage: s.one_basis()
                []
                sage: Q = Sym.hall_littlewood().Q()
                sage: Q.one_basis()
                []

            .. TODO:: generalize to Modules.Graded.Connected.ParentMethods
            """
            return _Partitions([])

        def degree_on_basis(self, b):
            r"""
            Return the degree of the basis element indexed by ``b``.

            INPUT:

            - ``self`` -- a basis of the symmetric functions
            - ``b`` -- a partition

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ['q,t'].fraction_field())
                sage: m = Sym.monomial()
                sage: m.degree_on_basis(Partition([3,2]))
                5
                sage: P = Sym.macdonald().P()
                sage: P.degree_on_basis(Partition([]))
                0
            """
            return sum(b)

        def corresponding_basis_over(self, R):
            r"""
            Return the realization of symmetric functions corresponding to
            ``self`` but over the base ring ``R``. Only works when ``self``
            is one of the classical bases, not one of the `q,t`-dependent
            ones. In the latter case, ``None`` is returned instead.

            INPUT:

            - ``R`` -- a commutative ring

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ)
                sage: m = Sym.monomial()
                sage: m.corresponding_basis_over(ZZ)
                Symmetric Functions over Integer Ring in the monomial basis

                sage: Sym = SymmetricFunctions(CyclotomicField())
                sage: s = Sym.schur()
                sage: s.corresponding_basis_over(Integers(13))
                Symmetric Functions over Ring of integers modulo 13 in the Schur basis

                sage: P = ZZ['q','t']
                sage: Sym = SymmetricFunctions(P)
                sage: mj = Sym.macdonald().J()
                sage: mj.corresponding_basis_over(Integers(13))

            TESTS:

            Let's check that this handles each of the bases properly::

                sage: P = QQ['q','t']
                sage: Sym = SymmetricFunctions(P)
                sage: Q = CyclotomicField()['q','t']
                sage: Sym.s().corresponding_basis_over(CyclotomicField())
                Symmetric Functions over Universal Cyclotomic Field in the Schur basis
                sage: Sym.p().corresponding_basis_over(CyclotomicField())
                Symmetric Functions over Universal Cyclotomic Field in the powersum basis
                sage: Sym.m().corresponding_basis_over(CyclotomicField())
                Symmetric Functions over Universal Cyclotomic Field in the monomial basis
                sage: Sym.e().corresponding_basis_over(CyclotomicField())
                Symmetric Functions over Universal Cyclotomic Field in the elementary basis
                sage: Sym.h().corresponding_basis_over(CyclotomicField())
                Symmetric Functions over Universal Cyclotomic Field in the homogeneous basis
                sage: Sym.f().corresponding_basis_over(CyclotomicField())
                Symmetric Functions over Universal Cyclotomic Field in the forgotten basis
                sage: Sym.w().corresponding_basis_over(CyclotomicField())
                Symmetric Functions over Universal Cyclotomic Field in the Witt basis
                sage: Sym.macdonald().P().corresponding_basis_over(CyclotomicField())
                sage: Sym.macdonald().Q().corresponding_basis_over(CyclotomicField())
                sage: Sym.macdonald().J().corresponding_basis_over(CyclotomicField())
                sage: Sym.macdonald().H().corresponding_basis_over(CyclotomicField())
                sage: Sym.macdonald().Ht().corresponding_basis_over(CyclotomicField())
                sage: Sym.macdonald().S().corresponding_basis_over(CyclotomicField())
                sage: Sym.macdonald(q=1).S().corresponding_basis_over(CyclotomicField())
                sage: Sym.macdonald(q=1,t=3).P().corresponding_basis_over(CyclotomicField())
                sage: Sym.hall_littlewood().P().corresponding_basis_over(CyclotomicField())
                sage: Sym.hall_littlewood().Q().corresponding_basis_over(CyclotomicField())
                sage: Sym.hall_littlewood().Qp().corresponding_basis_over(CyclotomicField())
                sage: Sym.hall_littlewood(t=1).P().corresponding_basis_over(CyclotomicField())
                sage: Sym.jack().J().corresponding_basis_over(CyclotomicField())
                sage: Sym.jack().P().corresponding_basis_over(CyclotomicField())
                sage: Sym.jack().Q().corresponding_basis_over(CyclotomicField())
                sage: Sym.jack().Qp().corresponding_basis_over(CyclotomicField())
                sage: Sym.jack(t=1).J().corresponding_basis_over(CyclotomicField())
                sage: Sym.zonal().corresponding_basis_over(CyclotomicField())
                Symmetric Functions over Universal Cyclotomic Field in the zonal basis
                sage: Sym.llt(3).hspin().corresponding_basis_over(CyclotomicField())
                sage: Sym.llt(3).hcospin().corresponding_basis_over(CyclotomicField())
                sage: Sym.llt(3, t=1).hspin().corresponding_basis_over(CyclotomicField())
                sage: Sym.llt(3, t=1).hcospin().corresponding_basis_over(CyclotomicField())

            .. TODO::

                This function is an ugly hack using strings. It should be
                rewritten as soon as the bases of ``SymmetricFunctions`` are
                put on a more robust and systematic footing.
            """
            from sage.combinat.sf.sf import SymmetricFunctions
            from sage.misc.misc import attrcall
            try:
                return attrcall(self._basis)(SymmetricFunctions(R))
            except AttributeError:   # or except (AttributeError, ValueError):
                return None
            #Alternative code proposed by Florent Hivert, which sadly fails for the
            #forgotten basis (which reduces differently than the other ones):
            #try:
            #    parentred1 = self._reduction
            #    parentred2 = parentred1[1][0]._reduction
            #    parentred2prime = tuple([parentred2[0], tuple([R]), parentred2[2]])
            #    from sage.structure.unique_representation import unreduce
            #    parent2 = unreduce(*parentred2prime)
            #    parentred1prime = tuple([parentred1[0], tuple([parent2]), parentred1[2]])
            #    return unreduce(*parentred1prime)
            #except (AttributeError, ValueError):
            #    return None
            #This code relied heavily on the construction of bases of
            #``SymmetricFunctions`` and on their reduction.

        def skew_schur(self, x):
            """
            Return the skew Schur function indexed by ``x`` in ``self``.

            INPUT:

            - ``x`` -- a skew partition

            EXAMPLES::

                sage: sp = SkewPartition([[5,3,3,1], [3,2,1]])
                sage: s = SymmetricFunctions(QQ).s()
                sage: s.skew_schur(sp)
                s[2, 2, 1, 1] + s[2, 2, 2] + s[3, 1, 1, 1] + 3*s[3, 2, 1]
                 + s[3, 3] + 2*s[4, 1, 1] + 2*s[4, 2] + s[5, 1]

                sage: e = SymmetricFunctions(QQ).e()
                sage: ess = e.skew_schur(sp); ess
                e[2, 1, 1, 1, 1] - e[2, 2, 1, 1] - e[3, 1, 1, 1] + e[3, 2, 1]
                sage: ess == e(s.skew_schur(sp))
                True

            TESTS::

                sage: s.skew_schur([[2,1], [1]])
                s[1, 1] + s[2]

                sage: s.skew_schur([[2,1], [3]])
                Traceback (most recent call last):
                ...
                ValueError: not a valid skew partition
            """
            from sage.combinat.skew_partition import SkewPartitions
            if x not in SkewPartitions():
                raise ValueError("not a valid skew partition")
            import sage.libs.lrcalc.lrcalc as lrcalc
            s = self.realization_of().schur()
            skewschur = lrcalc.skew(x[0], x[1])
            return self(s._from_dict(skewschur))

        def Eulerian(self, n, j, k=None):
            """
            Return the Eulerian symmetric function `Q_{n,j}` (with `n`
            either an integer or a partition) or `Q_{n,j,k}` (if the
            optional argument ``k`` is specified) in terms of the basis
            ``self``.

            It is known that the Eulerian quasisymmetric functions are
            in fact symmetric functions [SW2010]_. For more information,
            see :meth:`QuasiSymmetricFunctions.Fundamental.Eulerian()`,
            which accepts the same syntax as this method.

            INPUT:

            - ``n`` -- the nonnegative integer `n` or a partition
            - ``j`` -- the number of excedances
            - ``k`` -- (optional) if specified, determines the number of fixed
              points of the permutations which are being summed over

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ)
                sage: m = Sym.m()
                sage: m.Eulerian(3, 1)
                4*m[1, 1, 1] + 3*m[2, 1] + 2*m[3]
                sage: h = Sym.h()
                sage: h.Eulerian(4, 2)
                h[2, 2] + h[3, 1] + h[4]
                sage: s = Sym.s()
                sage: s.Eulerian(5, 2)
                s[2, 2, 1] + s[3, 1, 1] + 5*s[3, 2] + 6*s[4, 1] + 6*s[5]
                sage: s.Eulerian([2,2,1], 2)
                s[2, 2, 1] + s[3, 2] + s[4, 1] + s[5]
                sage: s.Eulerian(5, 2, 2)
                s[3, 2] + s[4, 1] + s[5]

            We check Equation (5.4) in [SW2010]_::

                sage: h.Eulerian([6], 3)
                h[3, 2, 1] - h[4, 1, 1] + 2*h[4, 2] + h[5, 1]
                sage: s.Eulerian([6], 3)
                s[3, 2, 1] + s[3, 3] + 3*s[4, 2] + 3*s[5, 1] + 3*s[6]
            """
            from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
            F = QuasiSymmetricFunctions(self.base_ring()).F()
            if n in _Partitions:
                n = _Partitions(n)
            return self(F.Eulerian(n, j, k).to_symmetric_function())

        def gessel_reutenauer(self, lam):
            r"""
            Return the Gessel-Reutenauer symmetric function
            corresponding to the partition ``lam`` written in the basis
            ``self``.

            Let `\lambda` be a partition. The *Gessel-Reutenauer
            symmetric function* `\mathbf{GR}_\lambda` corresponding to
            `\lambda` is the symmetric function denoted `L_\lambda` in
            [GR1993]_ and in Exercise 7.89 of [STA]_. It can be defined
            in several ways:

            - It is the sum of the monomials `\mathbf{x}_w` over all
              words `w` over the alphabet
              `\left\{ 1, 2, 3, \ldots \right\}` which have CFL type
              `\lambda`. Here, the monomial `\mathbf{x}_w` for a word
              `w = \left(w_1, w_2, \ldots, w_k\right)` is defined as
              `x_{w_1} x_{w_2} \cdots x_{w_k}`, and the *CFL type* of
              a word `w` is defined as the partition obtained by
              sorting (in decreasing order) the lengths of the factors
              in the Lyndon factorization
              (:meth:`~sage.combinat.words.finite_word.FiniteWord_class.lyndon_factorization`)
              of `w`. The fact that this power series
              `\mathbf{GR}_\lambda` is symmetric is not obvious.

            - It is the sum of the fundamental quasisymmetric
              functions `F_{\operatorname{Des} \sigma}` over all
              permutations `\sigma` which have cycle type `\lambda`. See
              :class:`sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Fundamental`
              for the definition of fundamental quasisymmetric functions,
              and :meth:`~sage.combinat.permutation.Permutation.cycle_type`
              for that of cycle type. For a permutation `\sigma`, we use
              `\operatorname{Des} \sigma` to denote the descent composition
              (:meth:`~sage.combinat.permutation.Permutation.descents_composition`)
              of `\sigma`. Again, this definition makes the symmetry
              of `\mathbf{GR}_\lambda` far from obvious.

            - For every positive integer `n`, we have

              .. MATH::

                  \mathbf{GR}_{\left(n\right)}
                  = \frac{1}{n} \sum_{d \mid n} \mu(d) p_d^{n/d},

              where `p_d` denotes the `d`-th power-sum symmetric
              function. This `\mathbf{GR}_{\left(n\right)}` is also
              denoted by `L_n`. Now, `\mathbf{GR}_\lambda` is defined
              as the product:

              .. MATH::

                  h_{m_1} \left[L_1\right] \cdot h_{m_2} \left[L_2\right]
                  \cdot h_{m_3} \left[L_3\right] \cdots,

              where `m_i` denotes the multiplicity of the part `i` in
              `\lambda`, and where the square brackets stand for
              plethysm (:meth:`plethysm`). This definition makes
              the symmetry (but not the integrality!) of
              `\mathbf{GR}_\lambda` obvious.

            The equivalences of these three definitions are proven in
            [GR1993]_ Sections 2-3.

            INPUT:

            - ``lam`` -- a partition or a positive integer (in the latter
              case, it is understood to mean the partition ``[lam]``)

            OUTPUT:

            The Gessel-Reutenauer symmetric function
            `\mathbf{GR}_\lambda`, where `\lambda` is ``lam``,
            expanded in the basis ``self``.

            REFERENCES:

            .. [GR1993] Ira M. Gessel, Christophe Reutenauer.
               *Counting Permutations with Given Cycle Structure
               and Descent Set*.
               Journal of Combinatorial Theory, Series A, 64 (1993),
               pp. 189--215.

            EXAMPLES:

            The first few values of `\mathbf{GR}_{(n)} = L_n`::

                sage: Sym = SymmetricFunctions(ZZ)
                sage: h = Sym.h()
                sage: h.gessel_reutenauer(1)
                h[1]
                sage: h.gessel_reutenauer(2)
                h[1, 1] - h[2]
                sage: h.gessel_reutenauer(3)
                h[2, 1] - h[3]
                sage: h.gessel_reutenauer(4)
                h[2, 1, 1] - h[2, 2]
                sage: h.gessel_reutenauer(5)
                h[2, 1, 1, 1] - h[2, 2, 1] - h[3, 1, 1] + h[3, 2] + h[4, 1] - h[5]
                sage: h.gessel_reutenauer(6)
                h[2, 1, 1, 1, 1] - h[2, 2, 1, 1] - h[2, 2, 2]
                 - 2*h[3, 1, 1, 1] + 5*h[3, 2, 1] - 2*h[3, 3] + h[4, 1, 1]
                 - h[4, 2] - h[5, 1] + h[6]

            Gessel-Reutenauer functions indexed by partitions::

                sage: h.gessel_reutenauer([2, 1])
                h[1, 1, 1] - h[2, 1]
                sage: h.gessel_reutenauer([2, 2])
                h[1, 1, 1, 1] - 3*h[2, 1, 1] + 2*h[2, 2] + h[3, 1] - h[4]

            The Gessel-Reutenauer functions are Schur-positive::

                sage: s = Sym.s()
                sage: s.gessel_reutenauer([2, 1])
                s[1, 1, 1] + s[2, 1]
                sage: s.gessel_reutenauer([2, 2, 1])
                s[1, 1, 1, 1, 1] + s[2, 1, 1, 1] + s[2, 2, 1] + s[3, 2]

            They do not form a basis, as the following example (from
            [GR1993]_ p. 201) shows::

                sage: s.gessel_reutenauer([4]) == s.gessel_reutenauer([2, 1, 1])
                True

            Of the above three equivalent definitions of
            `\mathbf{GR}_\lambda`, we use the third one for
            computations. Let us check that the second one gives the
            same results::

                sage: QSym = QuasiSymmetricFunctions(ZZ)
                sage: F = QSym.F() # fundamental basis
                sage: def GR_def2(lam): # `\mathbf{GR}_\lambda`
                ....:     n = lam.size()
                ....:     r = F.sum_of_monomials([sigma.descents_composition()
                ....:                             for sigma in Permutations(n)
                ....:                             if sigma.cycle_type() == lam])
                ....:     return r.to_symmetric_function()
                sage: all( GR_def2(lam) == h.gessel_reutenauer(lam)
                ....:      for n in range(5) for lam in Partitions(n) )
                True

            And the first one, too (assuming symmetry)::

                sage: m = Sym.m()
                sage: def GR_def1(lam): # `\mathbf{GR}_\lambda`
                ....:     n = lam.size()
                ....:     Permus_mset = sage.combinat.permutation.Permutations_mset
                ....:     def coeff_of_m_mu_in_result(mu):
                ....:         words_to_check = Permus_mset([i for (i, l) in enumerate(mu)
                ....:                                       for _ in range(l)])
                ....:         return sum((1 for w in words_to_check if
                ....:                     Partition(list(reversed(sorted([len(v) for v in Word(w).lyndon_factorization()]))))
                ....:                     == lam))
                ....:     r = m.sum_of_terms([(mu, coeff_of_m_mu_in_result(mu))
                ....:                         for mu in Partitions(n)],
                ....:                        distinct=True)
                ....:     return r
                sage: all( GR_def1(lam) == h.gessel_reutenauer(lam)
                ....:      for n in range(5) for lam in Partitions(n) )
                True

            TESTS:

            This works fine over other base rings::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
                sage: P = Sym.macdonald().P()
                sage: h = Sym.h()
                sage: P.gessel_reutenauer(3) == P(h.gessel_reutenauer(3))
                True

            .. NOTE::

                The currently existing implementation of this function is
                technically unsatisfactory. It distinguishes the case when the
                base ring is a `\QQ`-algebra from the case
                where it isn't. In the latter, it does a computation using
                universal coefficients, again distinguishing the case when it is
                able to compute the "corresponding" basis of the symmetric function
                algebra over `\QQ` (using the ``corresponding_basis_over`` hack)
                from the case when it isn't (in which case it transforms everything
                into the Schur basis, which is slow).
            """
            if lam in ZZ:
                lam = [lam]
            lam = _Partitions(lam)
            R = self.base_ring()
            # We use [GR1993]_ Theorem 3.6 and work over `\QQ` to
            # compute the Gessel-Reutenauer symmetric function.
            if self.has_coerce_map_from(QQ):
                # [GR1993]_ Theorem 3.6
                m = lam.to_exp_dict() # == {i: m_i | i occurs in lam}
                p = self.realization_of().power()
                h = self.realization_of().complete()
                from sage.arith.all import Moebius, squarefree_divisors
                mu = Moebius()
                def component(i, g): # == h_g[L_i]
                    L_i = p.sum_of_terms([(_Partitions([d] * (i//d)), R(mu(d)))
                                          for d in squarefree_divisors(i)],
                                         distinct=True) / i
                    return p(h[g]).plethysm(L_i)
                return self( p.prod(component(i, g) for i, g in m.items()) )

            # The base ring does not coerce into `\QQ`

            # comp_parent is the parent that is going to be used for
            # computations. In most cases it will just be self.
            comp_parent = self
            # Now let's try to find out what basis self is in, and
            # construct the corresponding basis of symmetric functions
            # over QQ.
            corresponding_parent_over_QQ = self.corresponding_basis_over(QQ)
            if corresponding_parent_over_QQ is None:
                # This is the case where the corresponding basis
                # over QQ cannot be found. This can have two reasons:
                # Either the basis depends on variables (like the
                # Macdonald symmetric functions), or its basis_name()
                # is not identical to the name of the method on
                # SymmetricFunctions(QQ) that builds it. Either way,
                # give up looking for the corresponding parent, and
                # transform everything into the Schur basis (very
                # slow!) instead.
                comp_parent = self.realization_of().schur()
                from sage.combinat.sf.sf import SymmetricFunctions
                corresponding_parent_over_QQ = SymmetricFunctions(QQ).schur()
            corresponding_result = corresponding_parent_over_QQ.gessel_reutenauer(lam)
            comp_base_ring = comp_parent.base_ring()
            result = comp_parent.sum_of_terms((nu, comp_base_ring(c))
                                               for nu, c in corresponding_result)
            return self(result)    # just in case comp_parent != self.

        def carlitz_shareshian_wachs(self, n, d, s, comparison=None):
            r"""
            Return the Carlitz-Shareshian-Wachs symmetric function
            `X_{n, d, s}` (if ``comparison`` is ``None``), or
            `U_{n, d, s}` (if ``comparison`` is ``-1``), or
            `V_{n, d, s}` (if ``comparison`` is ``0``), or
            `W_{n, d, s}` (if ``comparison`` is ``1``) written in the
            basis ``self``. These functions are defined below.

            The Carlitz-Shareshian-Wachs symmetric functions have been
            introduced in [GriRei2014]_, Exercise 2.84, as
            refinements of a certain particular case of chromatic
            quasisymmetric functions defined by Shareshian and Wachs.
            Their definitions are as follows:

            Let `n`, `d` and `s` be three nonnegative integers. Let
            `W(n, d, s)` denote the set of all `n`-tuples
            `(w_1, w_2, \ldots, w_n)` of positive integers having the
            property that there exist precisely `d` elements `i`
            of `\left\{ 1, 2, \ldots, n-1 \right\}` satisfying
            `w_i > w_{i+1}`, and precisely `s` elements `i` of
            `\left\{ 1, 2, \ldots, n-1 \right\}` satisfying
            `w_i = w_{i+1}`. For every
            `w = (w_1, w_2, \ldots, w_n) \in W(n, d, s)`, let `x_w`
            be the monomial `x_{w_1} x_{w_2} \cdots x_{w_n}`. We then
            define the power series `X_{n, d, s}` by

            .. MATH::

                X_{n, d, s} = \sum_{w \in W(n, d, s)} x_w .

            This is a symmetric function (according to
            [GriRei2014]_, Exercise 2.84(b)), and for `s = 0` equals
            the `t^d`-coefficient of the descent enumerator of Smirnov
            words of length `n` (an example of a chromatic
            quasisymmetric function which happens to be symmetric --
            see [ShaWach2014]_, Example 2.5).

            Assume that `n > 0`. Then, we can define three further
            power series as follows:

            .. MATH::

                U_{n, d, s} = \sum_{w_1 < w_n} x_w ; \qquad
                V_{n, d, s} = \sum_{w_1 = w_n} x_w ; \qquad
                W_{n, d, s} = \sum_{w_1 > w_n} x_w ,

            where all three sums range over
            `w = (w_1, w_2, \ldots, w_n) \in W(n, d, s)`. These
            three power series `U_{n, d, s}`, `V_{n, d, s}` and
            `W_{n, d, s}` are symmetric functions as well
            ([GriRei2014]_, Exercise 2.84(c)). Their sum is
            `X_{n, d, s}`.

            REFERENCES:

            .. [ShaWach2014] John Shareshian, Michelle L. Wachs.
               *Chromatic quasisymmetric functions*.
               :arxiv:`1405.4629v1`.

            INPUT:

            - ``n`` -- a nonnegative integer

            - ``d`` -- a nonnegative integer

            - ``s`` -- a nonnegative integer

            - ``comparison`` (default: ``None``) -- a variable
              which can take the forms ``None``, ``-1``, ``0``
              and ``1``

            OUTPUT:

            The Carlitz-Shareshian-Wachs symmetric function
            `X_{n, d, s}` (if ``comparison`` is ``None``), or
            `U_{n, d, s}` (if ``comparison`` is ``-1``), or
            `V_{n, d, s}` (if ``comparison`` is ``0``), or
            `W_{n, d, s}` (if ``comparison`` is ``1``) written in the
            basis ``self``.

            EXAMPLES:

            The power series `X_{n, d, s}`::

                sage: Sym = SymmetricFunctions(ZZ)
                sage: m = Sym.m()
                sage: m.carlitz_shareshian_wachs(3, 2, 1)
                0
                sage: m.carlitz_shareshian_wachs(3, 1, 1)
                m[2, 1]
                sage: m.carlitz_shareshian_wachs(3, 2, 0)
                m[1, 1, 1]
                sage: m.carlitz_shareshian_wachs(3, 0, 2)
                m[3]
                sage: m.carlitz_shareshian_wachs(3, 1, 0)
                4*m[1, 1, 1] + m[2, 1]
                sage: m.carlitz_shareshian_wachs(3, 0, 1)
                m[2, 1]
                sage: m.carlitz_shareshian_wachs(3, 0, 0)
                m[1, 1, 1]
                sage: m.carlitz_shareshian_wachs(5, 2, 2)
                m[2, 2, 1] + m[3, 1, 1]
                sage: m.carlitz_shareshian_wachs(1, 0, 0)
                m[1]
                sage: m.carlitz_shareshian_wachs(0, 0, 0)
                m[]

            The power series `U_{n, d, s}`::

                sage: m.carlitz_shareshian_wachs(3, 2, 1, comparison=-1)
                0
                sage: m.carlitz_shareshian_wachs(3, 1, 1, comparison=-1)
                0
                sage: m.carlitz_shareshian_wachs(3, 2, 0, comparison=-1)
                0
                sage: m.carlitz_shareshian_wachs(3, 0, 2, comparison=-1)
                0
                sage: m.carlitz_shareshian_wachs(3, 1, 0, comparison=-1)
                2*m[1, 1, 1]
                sage: m.carlitz_shareshian_wachs(3, 0, 1, comparison=-1)
                m[2, 1]
                sage: m.carlitz_shareshian_wachs(3, 0, 0, comparison=-1)
                m[1, 1, 1]
                sage: m.carlitz_shareshian_wachs(5, 2, 2, comparison=-1)
                0
                sage: m.carlitz_shareshian_wachs(4, 2, 0, comparison=-1)
                3*m[1, 1, 1, 1]
                sage: m.carlitz_shareshian_wachs(1, 0, 0, comparison=-1)
                0

            The power series `V_{n, d, s}`::

                sage: m.carlitz_shareshian_wachs(3, 2, 1, comparison=0)
                0
                sage: m.carlitz_shareshian_wachs(3, 1, 1, comparison=0)
                0
                sage: m.carlitz_shareshian_wachs(3, 2, 0, comparison=0)
                0
                sage: m.carlitz_shareshian_wachs(3, 0, 2, comparison=0)
                m[3]
                sage: m.carlitz_shareshian_wachs(3, 1, 0, comparison=0)
                m[2, 1]
                sage: m.carlitz_shareshian_wachs(3, 0, 1, comparison=0)
                0
                sage: m.carlitz_shareshian_wachs(3, 0, 0, comparison=0)
                0
                sage: m.carlitz_shareshian_wachs(5, 2, 2, comparison=0)
                0
                sage: m.carlitz_shareshian_wachs(4, 2, 0, comparison=0)
                m[2, 1, 1]
                sage: m.carlitz_shareshian_wachs(1, 0, 0, comparison=0)
                m[1]

            The power series `W_{n, d, s}`::

                sage: m.carlitz_shareshian_wachs(3, 2, 1, comparison=1)
                0
                sage: m.carlitz_shareshian_wachs(3, 1, 1, comparison=1)
                m[2, 1]
                sage: m.carlitz_shareshian_wachs(3, 2, 0, comparison=1)
                m[1, 1, 1]
                sage: m.carlitz_shareshian_wachs(3, 0, 2, comparison=1)
                0
                sage: m.carlitz_shareshian_wachs(3, 1, 0, comparison=1)
                2*m[1, 1, 1]
                sage: m.carlitz_shareshian_wachs(3, 0, 1, comparison=1)
                0
                sage: m.carlitz_shareshian_wachs(3, 0, 0, comparison=1)
                0
                sage: m.carlitz_shareshian_wachs(5, 2, 2, comparison=1)
                m[2, 2, 1] + m[3, 1, 1]
                sage: m.carlitz_shareshian_wachs(4, 2, 0, comparison=1)
                8*m[1, 1, 1, 1] + 2*m[2, 1, 1] + m[2, 2]
                sage: m.carlitz_shareshian_wachs(1, 0, 0, comparison=1)
                0

            TESTS:

            This works fine over other base rings::

                sage: Sym = SymmetricFunctions(FractionField(QQ['q','t']))
                sage: P = Sym.macdonald().P()
                sage: m = Sym.m()
                sage: m.carlitz_shareshian_wachs(4, 1, 1)
                4*m[2, 1, 1] + 2*m[2, 2] + 2*m[3, 1]
                sage: P.carlitz_shareshian_wachs(4, 1, 1) == P(m.carlitz_shareshian_wachs(4, 1, 1))
                True
            """
            # Stupid implementation.
            R = self.base_ring()
            m = self.realization_of().m()
            from sage.combinat.permutation import Permutations_mset
            # Defining a ``check_word`` function. This function will be used
            # to check if an `n`-tuple `w` of positive integers belongs to
            # `W(n, d, s)` and satisfies the additional requirement
            # determined by ``comparison``.
            # The ``comparison`` check has been factored out so that
            # ``comparison`` needs not be called a myriad of times. Might
            # be folly.
            if comparison is None:
                def check_word(w):
                    if sum((1 for i in range(n-1) if w[i] > w[i+1])) != d:
                        return False
                    if sum((1 for i in range(n-1) if w[i] == w[i+1])) != s:
                        return False
                    return True
            elif comparison == -1:
                def check_word(w):
                    if sum((1 for i in range(n-1) if w[i] > w[i+1])) != d:
                        return False
                    if sum((1 for i in range(n-1) if w[i] == w[i+1])) != s:
                        return False
                    return (w[0] < w[-1])
            elif comparison == 0:
                def check_word(w):
                    if sum((1 for i in range(n-1) if w[i] > w[i+1])) != d:
                        return False
                    if sum((1 for i in range(n-1) if w[i] == w[i+1])) != s:
                        return False
                    return (w[0] == w[-1])
            elif comparison == 1:
                def check_word(w):
                    if sum((1 for i in range(n-1) if w[i] > w[i+1])) != d:
                        return False
                    if sum((1 for i in range(n-1) if w[i] == w[i+1])) != s:
                        return False
                    return (w[0] > w[-1])

            def coeff_of_m_mu_in_result(mu):
                # Compute the coefficient of the monomial symmetric
                # function ``m[mu]`` in the result.
                words_to_check = Permutations_mset([i for (i, l) in enumerate(mu)
                                                    for _ in range(l)])
                return R( sum(1 for w in words_to_check if check_word(w)) )

            from sage.combinat.partition import Partitions_n
            r = m.sum_of_terms([(mu, coeff_of_m_mu_in_result(mu))
                                for mu in Partitions_n(n)],
                               distinct=True)
            return self(r)

class FilteredSymmetricFunctionsBases(Category_realization_of_parent):
    r"""
    The category of filtered bases of the ring of symmetric functions.

    TESTS::

        sage: from sage.combinat.sf.sfa import FilteredSymmetricFunctionsBases
        sage: Sym = SymmetricFunctions(QQ)
        sage: bases = FilteredSymmetricFunctionsBases(Sym); bases
        Category of filtered bases of Symmetric Functions over Rational Field
        sage: Sym.schur() in bases
        True
        sage: Sym.sp() in bases
        True
    """
    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import FilteredSymmetricFunctionsBases
            sage: Sym = SymmetricFunctions(QQ)
            sage: bases = FilteredSymmetricFunctionsBases(Sym)
            sage: bases._repr_()
            'Category of filtered bases of Symmetric Functions over Rational Field'
        """
        return "Category of filtered bases of %s" % self.base()

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import FilteredSymmetricFunctionsBases
            sage: Sym = SymmetricFunctions(QQ)
            sage: bases = FilteredSymmetricFunctionsBases(Sym)
            sage: bases.super_categories()
            [Category of bases of Symmetric Functions over Rational Field,
             Join of Category of hopf algebras with basis over Rational Field
                 and Category of filtered algebras with basis over Rational Field
                 and Category of commutative algebras over Rational Field]
        """
        cat = HopfAlgebras(self.base().base_ring()).Commutative().WithBasis().Filtered()
        return [SymmetricFunctionsBases(self.base()), cat]

class GradedSymmetricFunctionsBases(Category_realization_of_parent):
    r"""
    The category of graded bases of the ring of symmetric functions.

    These are further required to have the property that the basis element
    indexed by the empty partition is `1`.

    TESTS::

        sage: from sage.combinat.sf.sfa import GradedSymmetricFunctionsBases
        sage: Sym = SymmetricFunctions(QQ)
        sage: bases = GradedSymmetricFunctionsBases(Sym); bases
        Category of graded bases of Symmetric Functions over Rational Field
        sage: Sym.schur() in bases
        True
        sage: Sym.sp() in bases
        False
    """
    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import GradedSymmetricFunctionsBases
            sage: Sym = SymmetricFunctions(QQ)
            sage: bases = GradedSymmetricFunctionsBases(Sym)
            sage: bases._repr_()
            'Category of graded bases of Symmetric Functions over Rational Field'
        """
        return "Category of graded bases of %s" % self.base()

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.sf.sfa import GradedSymmetricFunctionsBases
            sage: Sym = SymmetricFunctions(QQ)
            sage: bases = GradedSymmetricFunctionsBases(Sym)
            sage: bases.super_categories()
            [Category of filtered bases of Symmetric Functions over Rational Field,
             Category of commutative graded hopf algebras with basis over Rational Field]
        """
        cat = HopfAlgebras(self.base().base_ring()).Commutative().WithBasis().Graded()
        return [FilteredSymmetricFunctionsBases(self.base()), cat]

    class ParentMethods:
        def antipode_by_coercion(self, element):
            r"""
            The antipode of ``element``.

            INPUT:

            - ``element`` -- element in a basis of the ring of symmetric functions

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ)
                sage: p = Sym.p()
                sage: s = Sym.s()
                sage: e = Sym.e()
                sage: h = Sym.h()
                sage: (h([]) + h([1])).antipode() # indirect doctest
                h[] - h[1]
                sage: (s([]) + s([1]) + s[2]).antipode()
                s[] - s[1] + s[1, 1]
                sage: (p([2]) + p([3])).antipode()
                -p[2] - p[3]
                sage: (e([2]) + e([3])).antipode()
                e[1, 1] - e[1, 1, 1] - e[2] + 2*e[2, 1] - e[3]
                sage: f = Sym.f()
                sage: f([3,2,1]).antipode()
                -f[3, 2, 1] - 4*f[3, 3] - 2*f[4, 2] - 2*f[5, 1] - 6*f[6]

            The antipode is an involution::

                sage: Sym = SymmetricFunctions(ZZ)
                sage: s = Sym.s()
                sage: all( s[u].antipode().antipode() == s[u] for u in Partitions(4) )
                True

            The antipode is an algebra homomorphism::

                sage: Sym = SymmetricFunctions(FiniteField(23))
                sage: h = Sym.h()
                sage: all( all( (s[u] * s[v]).antipode() == s[u].antipode() * s[v].antipode()
                ....:           for u in Partitions(3) )
                ....:      for v in Partitions(3) )
                True

            TESTS:

            Everything works over `\ZZ`::

                sage: Sym = SymmetricFunctions(ZZ)
                sage: p = Sym.p()
                sage: s = Sym.s()
                sage: e = Sym.e()
                sage: h = Sym.h()
                sage: (h([]) + h([1])).antipode() # indirect doctest
                h[] - h[1]
                sage: (s([]) + s([1]) + s[2]).antipode()
                s[] - s[1] + s[1, 1]
                sage: (p([2]) + p([3])).antipode()
                -p[2] - p[3]
                sage: (e([2]) + e([3])).antipode()
                e[1, 1] - e[1, 1, 1] - e[2] + 2*e[2, 1] - e[3]
            """
            return self.degree_negation(element.omega())

        def counit(self, element):
            r"""
            Return the counit of ``element``.

            The counit is the constant term of ``element``.

            INPUT:

            - ``element`` -- element in a basis of the ring of symmetric functions

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ)
                sage: m = Sym.monomial()
                sage: f = 2*m[2,1] + 3*m[[]]
                sage: f.counit()
                3
            """
            return element.degree_zero_coefficient()

        def degree_negation(self, element):
            r"""
            Return the image of ``element`` under the degree negation
            automorphism of the ring of symmetric functions.

            The degree negation is the automorphism which scales every
            homogeneous element of degree `k` by `(-1)^k` (for all `k`).

            INPUT:

            - ``element`` -- symmetric function written in ``self``

            EXAMPLES::

                sage: Sym = SymmetricFunctions(ZZ)
                sage: m = Sym.monomial()
                sage: f = 2*m[2,1] + 4*m[1,1] - 5*m[1] - 3*m[[]]
                sage: m.degree_negation(f)
                -3*m[] + 5*m[1] + 4*m[1, 1] - 2*m[2, 1]

            TESTS:

            Using :meth:`degree_negation` on an element of a different
            basis works correctly::

                sage: e = Sym.elementary()
                sage: m.degree_negation(e[3])
                -m[1, 1, 1]
                sage: m.degree_negation(m(e[3]))
                -m[1, 1, 1]
            """
            return self.sum_of_terms([ (lam, (-1)**(sum(lam)%2) * a)
                                       for lam, a in self(element) ])

    class ElementMethods:
        def degree_negation(self):
            r"""
            Return the image of ``self`` under the degree negation
            automorphism of the ring of symmetric functions.

            The degree negation is the automorphism which scales every
            homogeneous element of degree `k` by `(-1)^k` (for all `k`).

            Calling ``degree_negation(self)`` is equivalent to calling
            ``self.parent().degree_negation(self)``.

            EXAMPLES::

                sage: Sym = SymmetricFunctions(ZZ)
                sage: m = Sym.monomial()
                sage: f = 2*m[2,1] + 4*m[1,1] - 5*m[1] - 3*m[[]]
                sage: f.degree_negation()
                -3*m[] + 5*m[1] + 4*m[1, 1] - 2*m[2, 1]
                sage: x = m.zero().degree_negation(); x
                0
                sage: parent(x) is m
                True
            """
            return self.parent().sum_of_terms([ (lam, (-1)**(sum(lam)%2) * a)
                                                for lam, a in self ])

        def degree_zero_coefficient(self):
            r"""
            Returns the degree zero coefficient of ``self``.

            EXAMPLES::

                sage: Sym = SymmetricFunctions(QQ)
                sage: m = Sym.monomial()
                sage: f = 2*m[2,1] + 3*m[[]]
                sage: f.degree_zero_coefficient()
                3
            """
            return self.coefficient([])

#SymmetricFunctionsBases.Filtered = FilteredSymmetricFunctionsBases
#SymmetricFunctionsBases.Graded = GradedSymmetricFunctionsBases

#####################################################################
## ABC for bases of the symmetric functions

class SymmetricFunctionAlgebra_generic(CombinatorialFreeModule):
    r"""
    Abstract base class for symmetric function algebras.

    .. TODO::

        Most of the methods in this class are generic (manipulations of
        morphisms, ...) and should be generalized (or removed)

    TESTS::

        sage: s = SymmetricFunctions(QQ).s()
        sage: m = SymmetricFunctions(ZZ).m()
        sage: s(m([2,1]))
        -2*s[1, 1, 1] + s[2, 1]
    """
    def __init__(self, Sym, basis_name=None, prefix=None, graded=True):
        r"""
        Initializes the symmetric function algebra.

        INPUT:

        - ``Sym`` -- the ring of symmetric functions
        - ``basis_name`` -- name of basis (default: ``None``)
        - ``prefix`` -- prefix used to display basis
        - ``graded`` -- (default: ``True``) if ``True``, then the basis is
          considered to be graded, otherwise the basis is filtered

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
            raise TypeError("Argument R must be a commutative ring.")
        try:
            R(Integer(1))
        except Exception:
            raise ValueError("R must have a unit element")

        if basis_name is not None:
            self._basis = basis_name
        if prefix is not None:
            self._prefix = prefix
        self._sym = Sym
        if graded:
            cat = GradedSymmetricFunctionsBases(Sym)
        else: # Right now, there are no non-filted bases
            cat = FilteredSymmetricFunctionsBases(Sym)
        CombinatorialFreeModule.__init__(self, Sym.base_ring(), _Partitions,
                                         category=cat,
                                         bracket="", prefix=prefix)

    _print_style = 'lex'

    # Todo: share this with ncsf and over algebras with basis indexed by word-like elements
    def __getitem__(self, c, *rest):
        r"""
        This method implements the abuses of notations ``p[2,1]``,
        ``p[[2,1]]``, ``p[Partition([2,1])]``.

        INPUT:

        - ``c`` -- a list, list of lists, or partition

        .. TODO::

            Should call ``super.term`` so as not to interfere with the
            standard notation ``p['x,y,z']``.

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
            if len(rest) != 0:
                raise ValueError("invalid number of arguments")
        else:
            if len(rest) > 0 or isinstance(c, int) or isinstance(c, Integer):
                c = C([c]+list(rest))
            else:
                c = C(list(c))
        return self.monomial(c)

    def _change_by_proportionality(self, x, function):
        r"""
        Return the symmetric function obtained from ``x`` by scaling
        each basis element corresponding to the partition `\lambda` by
        the value of ``function`` on `\lambda`.

        INPUT:

        - ``x`` -- a symmetric function
        - ``function`` -- a function which takes in a partition
          and returns a scalar

        OUTPUT:

        A symmetric function in ``self`` which is a scaled version of ``x``.

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
        Return the plethysm of ``x`` by ``expr``.

        INPUT:

        - ``x`` -- a symmetric function
        - ``expr`` -- an expression used in the plethysm
        - ``deg_one`` -- a list (or iterable) specifying the degree one
          variables (that is, the terms to be treated as degree-one
          elements when encountered in ``x``; they will be taken to the
          appropriate powers when computing the plethysm)

        OUTPUT:

        The plethysm of ``x`` by ``expr``.

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
        p = self.realization_of().power()
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
        Applies morphism specified by ``f`` on (``x``, ``y``).

        INPUT:

        - ``x`` -- an element of ``self``
        - ``y`` -- an element of ``self``
        - ``f`` -- a function that takes in two partitions
          (basis elements) and returns an element of the target domain
        - ``orthogonal`` -- if orthogonal is set to ``True``, then
          ``f(part1, part2)`` is assumed to be 0 if ``part1 != part2``.

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
        Return the element of ``self`` with the same 'internal structure' as
        ``x``. This means the element whose coefficients in the basis ``self``
        are the respective coefficients of ``x`` in the basis of ``x``.

        INPUT:

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
        Return the image of an element ``element`` of some realization `M`
        of the ring of symmetric functions under a linear map from `M` to
        ``self`` whose matrix representation is cached (with ``cache_dict``
        serving as cache, and ``cache_function`` as the function which
        precomputes this cache).

        INPUT:

        -  ``element`` -- an element of a realization `M` of the ring of
           symmetric functions. Note that `M` can be a different realization
           than the one in which ``self`` is written, and does not have to
           be specified. It is assumed that the basis of ``self`` is indexed
           by partitions, and the degree of a basis element is the size of
           the partition indexing it.

        -  ``cache_function`` -- a function which accepts an
           integer `n` as its input and creates the cache for that homogeneous
           component (saving it in ``cache_dict``).

        -  ``cache_dict`` -- a dictionary storing a cache.
           It should be indexed by the positive integers `n`. Its values
           are dictionaries indexed by the partitions of size `n`. The values
           of those latter dictionaries are, again, dictionaries indexed by
           partitions of size `n`. Altogether, ``cache_dict`` should be
           understood to encode a graded linear map from `M` to the
           realization ``self`` of the ring of symmetric functions; the
           encoding is done in such a way that, for any `n` and any partitions
           ``lam`` and ``mu`` of `n`, the ``self[mu]``-coordinate of the image
           of ``M[lam]`` under this linear map (in the basis ``self``) is
           ``cache_dict[lam][mu]``.

        -  ``subs_dict`` -- (optional) a dictionary for any substitutions
           to make after the value is extracted from ``cache_dict``.

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
        # Convert x to the monomial basis
        BR = self.base_ring()
        zero = BR.zero()
        z_elt = {}
        for part, c in element.monomial_coefficients().iteritems():
            if sum(part) not in cache_dict:
                cache_function(sum(part))
            # Make sure it is a partition (for #13605), this is
            #   needed for the old kschur functions - TCS
            part = _Partitions(part)
            for part2, c2 in cache_dict[sum(part)][part].iteritems():
                if hasattr(c2,'subs'): # c3 may be in the base ring
                    c3 = c*BR(c2.subs(**subs_dict))
                else:
                    c3 = c*BR(c2)
                # c3 = c*c2
                # if hasattr(c3,'subs'): # c3 may be in the base ring
                #     c3 = c3.subs(**subs_dict)
                z_elt[ part2 ] = z_elt.get(part2, zero) + BR(c3)
        return self._from_dict(z_elt)

    def _invert_morphism(self, n, base_ring, self_to_other_cache, other_to_self_cache,\
                         to_other_function=None, to_self_function=None, \
                         upper_triangular=False, lower_triangular=False, \
                         ones_on_diagonal=False):
        r"""
        Compute the inverse of a morphism between ``self`` and ``other``
        (more precisely, its `n`-th graded component).

        In order to use this, you must be able to compute the morphism in
        one direction. This method assumes that the morphism is indeed
        invertible.

        INPUT:

        -  ``n`` -- an integer, the homogeneous component of
           symmetric functions for which we want to a morphism's inverse

        -  ``base_ring`` -- the base ring being worked over

        -  ``self_to_other_cache`` -- a dictionary which
           stores the transition from ``self`` to ``other``

        -  ``other_to_self_cache`` -- a dictionary which
           stores the transition from ``other`` to ``self``

        -  ``to_other_function`` -- a function which takes in
           a partition and returns a function which gives the coefficients of
           ``self(part)`` in the ``other`` basis

        -  ``to_self_function`` -- a function which takes in a
           partition and returns a function which gives the coefficients of
           ``other(part)`` in ``self``

        -  ``upper_triangular`` -- a boolean, if ``True``, the
           inverse will be computed by back substitution

        -  ``lower_triangular`` -- a boolean, if ``True``, the
           inverse will be computed by forward substitution

        -  ``ones_on_diagonal`` -- a boolean, if ``True``, the
           entries on the diagonal of the morphism (and inverse) matrix are
           assumed to be ones. This is used to remove divisions from the
           forward and back substitute algorithms.

        OUTPUT:

        Nothing is returned, but the caches ``self_to_other_cache``
        and ``other_to_self_cache`` are updated with the `n`-th degree
        components of the respective transition matrices.

        EXAMPLES:

        First, we will do an example of inverting the morphism
        which sends a Schur function to its conjugate Schur function. Note
        that this is an involution. ::

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

        We can check that we get the same results if we specify
        ``to_self_function = conj``::

            sage: d1 = {}
            sage: d2 = {}
            sage: s._invert_morphism(4, QQ, d1, d2, to_self_function = conj)
            sage: d1 == c1
            True
            sage: d2 == c2
            True

        Now we do an example of upper triangularity and check that we get
        the same thing whether or not we specify ``ones_on_diagonal``::

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

        Finally, we do the same thing for lower triangular matrices::

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
        one  = base_ring.one()
        zero = base_ring.zero()

        #Get and store the list of partitions we'll need
        pn = Partitions_n(n).list()
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
            raise ValueError("only one of upper_triangular and lower_triangular can be specified")
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
        Return the family of symmetric functions associated to the
        basis ``self``.

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
        Return the prefix on the elements of ``self``.

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
        Return the transition matrix between ``self`` and ``basis`` for the
        homogeneous component of degree ``n``.

        INPUT:

        - ``basis`` -- a basis of the ring of symmetric functions
        - ``n`` -- a nonnegative integer

        OUTPUT:

        - a matrix of coefficients giving the expansion of the
          homogeneous degree-`n` elements of ``self`` in the
          degree-`n` elements of ``basis``

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
            sage: mon = sorted(a.support())
            sage: coeffs = [a[i] for i in mon]
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
        P = Partitions_n(n)
        Plist = P.list()
        m = []
        for row_part in Plist:
            z = basis(self(row_part))
            m.append( [z.coefficient(col_part) for col_part in Plist] )
        return matrix(m)


    def _gram_schmidt(self, n, source, scalar, cache, leading_coeff=None, upper_triangular=True):
        r"""
        Apply Gram-Schmidt to ``source`` with respect to the scalar product
        ``scalar`` for all partitions of `n`. The scalar product is supposed
        to make the power-sum basis orthogonal. The Gram-Schmidt algorithm
        computes an orthogonal basis (with respect to the scalar product
        given by ``scalar``) of the `n`-th homogeneous component of the
        ring of symmetric functions such that the transition matrix from
        the basis ``source`` to this orthogonal basis is triangular.

        The result is not returned, but instead, the coefficients of the
        elements of the orthogonal basis with respect to the basis
        ``source`` are stored in the cache ``cache``.

        The implementation uses the powersum basis, so this function
        shouldn't be used unless the base ring is a `\QQ`-algebra
        (or ``self`` and ``source`` are both the powersum basis).

        INPUT:

        - ``n`` -- nonnegative integer which specifies the size of
          the partitions
        - ``source`` -- a basis of the ring of symmetric functions
        - ``scalar`` -- a function ``zee`` from partitions to the base ring
          which specifies the scalar product by `\langle p_{\lambda},
          p_{\lambda} \rangle = \mathrm{zee}(\lambda)`.
        - ``cache`` -- a cache function
        - ``leading_coeff`` -- (default: ``None``) specifies the leading
          coefficients for Gram-Schmidt
        - ``upper_triangular`` -- (defaults to ``True``) boolean, indicates
          whether the transition is upper triangular or not

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
        BR = self.base_ring(); one = BR.one()
        p = self.realization_of().p()

        # Create a function which converts x and y to the power-sum basis and applies
        # the scalar product.
        pscalar = lambda x,y: p._apply_multi_module_morphism(p(x), p(y), lambda a,b:scalar(a), orthogonal=True)

        if leading_coeff is None:
            leading_coeff = lambda x: one

        # We are going to be doing everything like we are in the upper-triangular case
        # We list the partitions in "decreasing order" and work from the beginning forward.
        # If we are in the lower-triangular case, then we shouldn't reverse the list
        l = Partitions_n(n).list()
        if upper_triangular:
            l.reverse()

        # precomputed elements
        precomputed_elements = []

        # Handle the initial case
        cache[l[0]] = { l[0]: leading_coeff(l[0]) }
        precomputed_elements.append(leading_coeff( l[0] )*source(l[0]))

        for i in range(1, len(l)):
            start = leading_coeff( l[i] )*source(l[i])
            sub = 0
            for j in range(i):
                sub += pscalar( start, precomputed_elements[j] ) / pscalar(precomputed_elements[j], precomputed_elements[j]) * precomputed_elements[j]
            res = start - sub

            if hasattr(self, '_normalize_coefficients'):
                res = res.map_coefficients(self._normalize_coefficients)
            precomputed_elements.append(res)
            # Now, res == precomputed_elements[i]
            cache[l[i]] = {}
            for j in range(i+1):
                cache[l[i]][l[j]] = res.coefficient(l[j])


    def _inner_plethysm_pk_g(self, k, g, cache):
        r"""
        Return the inner plethysm between the power-sum symmetric
        function `p_k` and the symmetric function ``g``.

        See :meth:`inner_plethysm` for the definition of inner
        plethysm.

        .. WARNING::

            The function ``g`` *must* be given in the power-sum
            basis for this method to return a correct result.

        ALGORITHM:

        Express ``g`` in the power sum basis as
        `g = \sum_\mu c_\mu p_\mu/z_\mu`
        (where `z_\mu` is the size of the centralizer of any
        permutation with cycle type `\mu`). Then, the inner plethysm
        is calculated as

        .. MATH::

            p_k \{ g \} = \sum_\mu c_\mu p_k \{ p_\mu/z_\mu \}~.

        The inner plethysm `p_k \{ p_mu/z_\mu \}` is given by the formula

        .. MATH::

            p_k \{ p_\mu/z_\mu \} = \sum_{\nu : \nu^k = \mu } p_{\nu}/z_{\nu}~,

        where `\nu^k` is the `k`-th power of `nu` (see
        :~sage.combinat.partition.partition_power`).

        .. SEEALSO:: :func:`~sage.combinat.partition.partition_power`,
            :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.inner_plethysm`

        INPUT:

        -  ``k`` -- a positive integer

        -  ``g`` -- a symmetric function in the power sum basis

        -  ``cache`` -- a dictionary whose keys are (k, g) pairs
           and values are the cached output of this function

        EXAMPLES::

            sage: p = SymmetricFunctions(QQ).p()
            sage: p._inner_plethysm_pk_g(2, p([1,1,1]), {})
            p[1, 1, 1] + 3*p[2, 1]
            sage: p._inner_plethysm_pk_g(5, p([2,2,1,1,1]), {})
            p[2, 2, 1, 1, 1]
        """
        try:
            return cache[(k,g)]
        except KeyError:
            pass

        p = self.realization_of().p()
        res = 0
        degrees = uniq([ sum(m) for m in g.support() ])
        for d in degrees:
            for mu in Partitions_n(d):
                mu_k = mu.power(k)
                if mu_k in g:
                    res += g.coefficient(mu_k)*mu_k.centralizer_size()/mu.centralizer_size()*p(mu)

        cache[(k,g)] = res
        return res

    def _inner_plethysm_pnu_g(self, p_x, cache, nu):
        r"""
        Return the inner plethysm of the power-sum symmetric function
        `p_\nu` with another symmetric function ``p_x`` in the
        power-sum basis.

        See :meth:`inner_plethysm` for the definition of inner
        plethysm.

        .. WARNING::

            The function ``p_x`` *must* be given in the power-sum
            basis for this method to return a correct result.

        The computation uses the inner plethysm of `p_k` and ``p_x``
        and the identity

        .. MATH::

            (f \cdot g) \{ h \} = (f \{ h \}) \ast (g \{ h \})~.

        .. SEEALSO:: :meth:`_inner_plethysm_pk_g`, 
            :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.itensor`,
            :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.inner_plethysm`

        INPUT:

        - ``p_x`` -- a symmetric function in the power sum basis

        - ``cache`` -- a cache function

        - ``nu`` -- a partition

        Note that the order of the arguments is somewhat strange in order
        to facilitate partial function application.

        OUTPUT:

        - an element of the basis ``self``

        EXAMPLES::

            sage: p = SymmetricFunctions(QQ).p()
            sage: s = SymmetricFunctions(QQ).s()
            sage: p._inner_plethysm_pnu_g( p([1,1,1]), {}, Partition([2,1]))
            6*p[1, 1, 1]
            sage: p._inner_plethysm_pnu_g( p([1,1,1]), {}, Partition([]))
            1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]
            sage: s(_)
            s[3]
        """
        #We handle the constant term case separately.  It should be
        #the case that p([]).inner_tensor(s(mu)) = s([ mu.size() ]).
        #Here, we get the degrees of the homogeneous pieces of
        if not nu._list:
            s = self.realization_of().s()
            degrees = [ part.size() for part in p_x.support() ]
            degrees = uniq(degrees)
            if 0 in degrees:
                ext = self([])
            else:
                ext = 0
            return ext + self(sum([s([n]) for n in degrees if n!=0]))

        #For each k in nu, we compute the inner plethysm of
        #p_k with p_x
        res = [self._inner_plethysm_pk_g(k, p_x, cache) for k in nu]

        #To get the final answer, we compute the inner tensor product
        #of all the symmetric functions in res
        return self(reduce(lambda x, y: 0 if x==0 else x.itensor(y), res))


    def _dual_basis_default(self):
        """
        Returns the default value for ``self.dual_basis()``

        .. SEEALSO:: :meth:`dual_basis`

        EXAMPLES:

        This default implementation constructs the dual basis using
        the standard (Hall) scalar product::

            sage: Sym = SymmetricFunctions(QQ)
            sage: Sym.p()._dual_basis_default()
            Dual basis to Symmetric Functions over Rational Field in the powersum basis with respect to the Hall scalar product

        This is meant to be overriden by subclasses for which an
        explicit dual basis is known::

            sage: Sym.s()._dual_basis_default()
            Symmetric Functions over Rational Field in the Schur basis
            sage: Sym.h()._dual_basis_default()
            Symmetric Functions over Rational Field in the monomial basis
            sage: Sym.m()._dual_basis_default()
            Symmetric Functions over Rational Field in the homogeneous basis
            sage: Sym.f()._dual_basis_default()
            Symmetric Functions over Rational Field in the elementary basis
            sage: Sym.e()._dual_basis_default()
            Symmetric Functions over Rational Field in the forgotten basis
            sage: Sym.f()._dual_basis_default()
            Symmetric Functions over Rational Field in the elementary basis
        """
        return self.dual_basis(scalar=zee, scalar_name = "Hall scalar product")


    def dual_basis(self, scalar=None, scalar_name="", basis_name=None, prefix=None):
        r"""
        Return the dual basis of ``self`` with respect to the scalar
        product ``scalar``.

        INPUT:

        - ``scalar`` -- A function ``zee`` from partitions to the base ring
          which specifies the scalar product by `\langle p_{\lambda},
          p_{\lambda} \rangle = \mathrm{zee}(\lambda)`. (Independently on the
          function chosen, the power sum basis will always be orthogonal;
          the function ``scalar`` only determines the norms of the basis
          elements.) If ``scalar`` is None, then the standard (Hall) scalar
          product is used.
        - ``scalar_name`` -- name of the scalar function
        - ``prefix`` -- prefix used to display the basis

        EXAMPLES:

        The duals of the elementary symmetric functions with respect to the
        Hall scalar product are the forgotten symmetric functions.

        ::

            sage: e = SymmetricFunctions(QQ).e()
            sage: f = e.dual_basis(prefix='f'); f
            Dual basis to Symmetric Functions over Rational Field in the elementary basis with respect to the Hall scalar product
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
            Dual basis to Symmetric Functions over Rational Field in the powersum basis with respect to the Hall scalar product
            sage: q([2,1])^2
            4*q[2, 2, 1, 1]
            sage: p([2,1]).scalar(q([2,1]))
            1
            sage: p([2,1]).scalar(q([1,1,1]))
            0
        """
        import dual
        if scalar is None:
            if basis_name is None and prefix is None:
                return self._dual_basis_default()
            scalar = zee
            scalar_name = "Hall scalar product"
        return dual.SymmetricFunctionAlgebra_dual(self, scalar, scalar_name,
                                                  basis_name = basis_name,
                                                  prefix = prefix)

    def basis_name(self):
        r"""
        Return the name of the basis of ``self``.

        This is used for output and, for the classical bases of
        symmetric functions, to connect this basis with Symmetrica.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: s.basis_name()
            'Schur'
            sage: p = Sym.p()
            sage: p.basis_name()
            'powersum'
            sage: h = Sym.h()
            sage: h.basis_name()
            'homogeneous'
            sage: e = Sym.e()
            sage: e.basis_name()
            'elementary'
            sage: m = Sym.m()
            sage: m.basis_name()
            'monomial'
            sage: f = Sym.f()
            sage: f.basis_name()
            'forgotten'
        """
        return self._basis

    def get_print_style(self):
        r"""
        Return the value of the current print style for ``self``.

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
            self.print_options(generator_cmp = lambda x,y: cmp(x,y))
        elif ps == 'length':
            self.print_options(generator_cmp = lambda x,y: cmp(len(x), len(y)))
        elif ps == 'maximal_part':
            self.print_options(generator_cmp = lambda x,y: cmp(_lmax(x), _lmax(y)))
        else:
            raise ValueError("the print style must be one of lex, length, or maximal_part ")
        self._print_style = ps

    def _latex_term(self, m):
        r"""
        Latex terms (i.e. partitions) as plain lists (and not as
        ferrers diagrams).

        INPUT:

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
        Convert polynomial to a symmetric function in the monomial basis
        and then to the basis ``self``.

        INPUT:

        - ``poly`` -- a symmetric polynomial
        - ``check`` -- (default: ``True``) boolean, specifies whether
          the computation checks that the polynomial is indeed symmetric

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
        m = self.realization_of().m()
        return self(m.from_polynomial(poly, check=check))

    def product_by_coercion(self, left, right):
        r"""
        Return the product of elements ``left`` and ``right`` by coercion to
        the Schur basis.

        INPUT:

        - ``left``, ``right`` -- instances of this basis

        OUTPUT:

        - the product of ``left`` and ``right`` expressed in the basis ``self``

        EXAMPLES::

            sage: p = SymmetricFunctions(QQ).p()
            sage: p.product_by_coercion(p[3,1,1], p[2,2])
            p[3, 2, 2, 1, 1]
            sage: m = SymmetricFunctions(QQ).m()
            sage: m.product_by_coercion(m[2,1],m[1,1]) == m[2,1]*m[1,1]
            True
        """
        s = self.realization_of().schur()
        return self(s.product(s(left),s(right)))

    def coproduct_by_coercion(self, elt):
        r"""
        Return the coproduct of the element ``elt`` by coercion to
        the Schur basis.

        INPUT:

        - ``elt`` -- an instance of this basis

        OUTPUT:

        - The image of ``elt`` under the comultiplication (=coproduct)
          of the coalgebra of symmetric functions. The result is an
          element of the tensor squared of the basis ``self``.

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).m()
            sage: m[3,1,1].coproduct()
            m[] # m[3, 1, 1] + m[1] # m[3, 1] + m[1, 1] # m[3] + m[3] # m[1, 1] + m[3, 1] # m[1] + m[3, 1, 1] # m[]
            sage: m.coproduct_by_coercion(m[2,1])
            m[] # m[2, 1] + m[1] # m[2] + m[2] # m[1] + m[2, 1] # m[]
            sage: m.coproduct_by_coercion(m[2,1]) == m([2,1]).coproduct()
            True
            sage: McdH = SymmetricFunctions(QQ['q','t'].fraction_field()).macdonald().H()
            sage: McdH[2,1].coproduct()
            McdH[] # McdH[2, 1] + ((q^2*t-1)/(q*t-1))*McdH[1] # McdH[1, 1] + ((q*t^2-1)/(q*t-1))*McdH[1] # McdH[2] + ((q^2*t-1)/(q*t-1))*McdH[1, 1] # McdH[1] + ((q*t^2-1)/(q*t-1))*McdH[2] # McdH[1] + McdH[2, 1] # McdH[]
            sage: HLQp = SymmetricFunctions(QQ['t'].fraction_field()).hall_littlewood().Qp()
            sage: HLQp[2,1].coproduct()
            HLQp[] # HLQp[2, 1] + HLQp[1] # HLQp[1, 1] + HLQp[1] # HLQp[2] + HLQp[1, 1] # HLQp[1] + HLQp[2] # HLQp[1] + HLQp[2, 1] # HLQp[]
            sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
            sage: LLT = Sym.llt(3)
            sage: LLT.cospin([3,2,1]).coproduct()
            (t+1)*m[] # m[1, 1] + m[] # m[2] + (t+1)*m[1] # m[1] + (t+1)*m[1, 1] # m[] + m[2] # m[]
            sage: f = SymmetricFunctions(ZZ).f()
            sage: f[3].coproduct()
            f[] # f[3] + f[3] # f[]
            sage: f[3,2,1].coproduct()
            f[] # f[3, 2, 1] + f[1] # f[3, 2] + f[2] # f[3, 1] + f[2, 1] # f[3] + f[3] # f[2, 1] + f[3, 1] # f[2] + f[3, 2] # f[1] + f[3, 2, 1] # f[]
        """
        from sage.categories.tensor import tensor
        s = self.realization_of().schur()
        return self.tensor_square().sum(coeff * tensor([self(s[x]), self(s[y])])
                                        for ((x,y), coeff) in s(elt).coproduct())


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
        Return the outer plethysm of ``self`` with ``x``. This is
        implemented only over base rings which are `\QQ`-algebras.
        (To compute outer plethysms over general binomial rings, change
        bases to the fraction field.)

        The outer plethysm of `f` with `g` is commonly denoted by
        `f \left[ g \right]` or by `f \circ g`. It is an algebra map
        in `f`, but not (generally) in `g`.

        By default, the degree one elements are taken to be the
        generators for the ``self``'s base ring. This setting can be
        modified by specifying the ``include`` and ``exclude`` keywords.

        INPUT:

        -  ``x`` -- a symmetric function over the same base ring as
           ``self``

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
            sage: s(1).plethysm(s(0))
            s[]

        .. SEEALSO::

            :meth:`frobenius`
        """
        if not is_SymmetricFunction(x):
            raise TypeError("only know how to compute plethysms between symmetric functions")
        parent = self.parent()
        p = parent.realization_of().power()
        R = parent.base_ring()
        p_x = p(x)
        if self == parent.zero():
            return self

        #Handle degree one elements
        if include is not None and exclude is not None:
            raise RuntimeError("include and exclude cannot both be specified")

        try:
            degree_one = [R(g) for g in R.variable_names_recursive()]
        except AttributeError:
            try:
                degree_one = R.gens()
            except NotImplementedError:
                degree_one = []

        if include:
            degree_one = [R(g) for g in include]
        if exclude:
            degree_one = [g for g in degree_one if g not in exclude]


        #Takes in n, and returns a function which takes in a partition and
        #scales all of the parts of that partition by n
        scale_part = lambda n: lambda m: m.__class__(m.parent(), [i*n for i in m])

        raise_c = lambda n: lambda c: c.subs(**dict((str(g),g**n) for g in degree_one if g != 1))

        #Takes n an symmetric function f, and an n and returns the
        #symmetric function with all of its basis partitions scaled
        #by n
        pn_pleth = lambda f, n: f.map_support(scale_part(n))

        #Takes in a partition and applies
        f = lambda part: prod( pn_pleth(p_x.map_coefficients(raise_c(i)), i) for i in part )
        return parent(p._apply_module_morphism(p(self),f))

    __call__ = plethysm


    def inner_plethysm(self, x):
        r"""
        Return the inner plethysm of ``self`` with ``x``.

        Whenever `R` is a `\QQ`-algebra, and `f` and `g` are two
        symmetric functions over `R` such that the constant term of `f`
        is zero, the inner plethysm of `f` with `g` is a symmetric
        function over `R`, and the degree of this symmetric function is
        the same as the degree of `g`. We will denote the inner plethysm
        of `f` with `g` by `f \{ g \}` (in contrast to the notation of
        outer plethysm which is generally denoted `f [ g ]`); in Sage
        syntax, it is ``f.inner_plethysm(g)``.

        First we describe the axiomatic definition of the operation; see
        below for a representation-theoretic interpretation.
        In the following equations, we denote the outer product
        (i.e., the standard product on the ring of symmetric functions,
        :meth:`~sage.categories.algebras_with_basis.AlgebrasWithBasis.ParentMethods.product`)
        by `\cdot` and the Kronecker product (:meth:`itensor`) by `\ast`).

        .. MATH::

            (f + g) \{ h \} = f \{ h \} + g \{ h \}

            (f \cdot g) \{ h \} = (f \{ h \}) \ast (g \{ h \})

            p_k \{ f + g \} = p_k \{ f \} + p_k \{ g \}

        where `p_k` is the `k`-th power-sum symmetric function for every
        `k > 0`.

        Let `\sigma` be a permutation of cycle type `\mu` and let `\mu^k`
        be the cycle type of `\sigma^k`. Then,

        .. MATH::

            p_k \{ p_\mu/z_\mu \} = \sum_{\nu : \nu^k = \mu } p_{\nu}/z_{\nu}

        Since `(p_\mu/z_\mu)_{\mu}` is a basis for the symmetric
        functions, these four formulas define the symmetric function
        operation `f \{ g \}` for any symmetric functions `f` and `g`
        (where `f` has constant term `0`) by expanding `f` in the
        power sum basis and `g` in the dual basis `p_\mu/z_\mu`.

        .. SEEALSO:: :meth:`itensor`, :func:`~sage.combinat.partition.partition_power`,
            :meth:`plethysm`

        This operation admits a representation-theoretic interpretation
        in the case where `f` is a Schur function `s_\lambda` and
        `g` is a homogeneous degree `n` symmetric function with
        nonnegative integral coefficients in the Schur basis.
        The symmetric function `f \{ g \}` is the Frobenius
        image of the `S_n`-representation constructed as follows.

        The assumptions on `g` imply that `g` is the Frobenius image of a
        representation `\rho` of the symmetric group `S_n`:

        .. MATH::

            \rho : S_n \to GL_N.

        If the degree `N` of this representation is greater than or equal
        to the number of parts of `\lambda`, then `f`, which denotes `s_\lambda`,
        corresponds to the character of some irreducible `GL_N`-representation, say

        .. MATH::

            \sigma : GL_N \to GL_M.

        The composition `\sigma \circ \rho : S_n \to GL_M` is a representation
        of `S_n` whose Frobenius image is precisely `f \{ g \}`.

        If `N` is less than the number of parts of `\lambda`,
        then `f \{ g \}` is `0` by definition.

        When `f` is a symmetric function with constant term `\neq 0`, the
        inner plethysm `f \{ g \}` isn't well-defined in the ring of
        symmetric functions. Indeed, it is not clear how to define
        `1 \{ g \}`. The most sensible way to get around this probably is
        defining it as the infinite sum `h_0 + h_1 + h_2 + \cdots` (where
        `h_i` means the `i`-th complete homogeneous symmetric function)
        in the completion of this ring with respect to its grading. This is
        how [SchaThi1994]_ defines `1 \{ g \}`. The present method,
        however, sets it to be the sum of `h_i` over all `i` for which the
        `i`-th homogeneous component of `g` is nonzero. This is rather a
        hack than a reasonable definition. Use with caution!

        .. NOTE::

            If a symmetric function `g` is written in the form
            `g = g_0 + g_1 + g_2 + \cdots` with each `g_i` homogeneous
            of degree `i`, then
            `f \{ g \} = f \{ g_0 \} + f \{ g_1 \} + f \{ g_2 \} + \cdots`
            for every `f` with constant term `0`. But in general, inner
            plethysm is not linear in the second variable.

        REFERENCES:

        .. [King] King, R. Branching rules for `GL_m \supset \Sigma_n`
           and the evaluation of inner plethysms.
           J. Math. Phys. 15, 258 (1974) :doi:`10.1063/1.1666632`

        .. [SchaThi1994] Thomas Scharf, Jean-Yves Thibon.
           *A Hopf-algebra approach to inner plethysm*.
           Advances in Mathematics 104 (1994), pp. 30-58.
           ftp://ftp.mathe2.uni-bayreuth.de/axel/papers/scharf:a_hopf_algebra_approach_to_inner_plethysm.ps.gz

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the same
          base ring as ``self``

        OUTPUT:

        - an element of symmetric functions in the parent of ``self``

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

        ::

            sage: f = s([2,1]) + 2*s([3,1])
            sage: f.itensor(f)
            s[1, 1, 1] + s[2, 1] + 4*s[2, 1, 1] + 4*s[2, 2] + s[3] + 4*s[3, 1] + 4*s[4]
            sage: s( h([1,1]).inner_plethysm(f) )
            s[1, 1, 1] + s[2, 1] + 4*s[2, 1, 1] + 4*s[2, 2] + s[3] + 4*s[3, 1] + 4*s[4]

        ::

            sage: s([]).inner_plethysm(s([1,1]) + 2*s([2,1])+s([3]))
            s[2] + s[3]
            sage: [s([]).inner_plethysm(s(la)) for la in Partitions(4)]
            [s[4], s[4], s[4], s[4], s[4]]
            sage: s([3]).inner_plethysm(s([]))
            s[]
            sage: s[1,1,1,1].inner_plethysm(s[2,1])
            0
            sage: s[1,1,1,1].inner_plethysm(2*s[2,1])
            s[3]

        ::

            sage: p[3].inner_plethysm(p[3])
            0
            sage: p[3,3].inner_plethysm(p[3])
            0
            sage: p[3].inner_plethysm(p[1,1,1])
            p[1, 1, 1] + 2*p[3]
            sage: p[4].inner_plethysm(p[1,1,1,1]/24)
            1/24*p[1, 1, 1, 1] + 1/4*p[2, 1, 1] + 1/8*p[2, 2] + 1/4*p[4]
            sage: p[3,3].inner_plethysm(p[1,1,1])
            6*p[1, 1, 1] + 12*p[3]

        TESTS::

            sage: s(0).inner_plethysm(s(0))
            0
            sage: s(1).inner_plethysm(s(0))
            0
            sage: s(0).inner_plethysm(s(1))
            0
            sage: s(1).inner_plethysm(s(1))
            s[]
            sage: s(2).inner_plethysm(s(1))
            2*s[]
            sage: s(1).inner_plethysm(s(2))
            s[]
        """
        parent = self.parent()
        if self == parent.zero():
            return self
        p = parent.realization_of().power()
        cache = {}
        ip_pnu_g = parent._inner_plethysm_pnu_g
        return parent.sum( c*ip_pnu_g(p(x), cache, nu)
                           for (nu, c) in p(self).monomial_coefficients().iteritems() )


    def omega(self):
        r"""
        Return the image of ``self`` under the omega automorphism.

        The *omega automorphism* is defined to be the unique algebra
        endomorphism `\omega` of the ring of symmetric functions that
        satisfies `\omega(e_k) = h_k` for all positive integers `k`
        (where `e_k` stands for the `k`-th elementary symmetric
        function, and `h_k` stands for the `k`-th complete homogeneous
        symmetric function). It furthermore is a Hopf algebra
        endomorphism and an involution, and it is also known as the
        *omega involution*. It sends the power-sum symmetric function
        `p_k` to `(-1)^{k-1} p_k` for every positive integer `k`.

        The images of some bases under the omega automorphism are given by

        .. MATH::

            \omega(e_{\lambda}) = h_{\lambda}, \qquad
            \omega(h_{\lambda}) = e_{\lambda}, \qquad
            \omega(p_{\lambda}) = (-1)^{|\lambda| - \ell(\lambda)}
            p_{\lambda}, \qquad
            \omega(s_{\lambda}) = s_{\lambda^{\prime}},

        where `\lambda` is any partition, where `\ell(\lambda)` denotes
        the length (:meth:`~sage.combinat.partition.Partition.length`)
        of the partition `\lambda`, where `\lambda^{\prime}` denotes the
        conjugate partition
        (:meth:`~sage.combinat.partition.Partition.conjugate`) of
        `\lambda`, and where the usual notations for bases are used
        (`e` = elementary, `h` = complete homogeneous, `p` = powersum,
        `s` = Schur).

        The default implementation converts to the Schur basis, then
        performs the automorphism and changes back.

        :meth:`omega_involution()` is a synonym for the :meth:`omega()` method.

        EXAMPLES::

            sage: J = SymmetricFunctions(QQ).jack(t=1).P()
            sage: a = J([2,1]) + J([1,1,1])
            sage: a.omega()
            JackP[2, 1] + JackP[3]
            sage: J(0).omega()
            0
            sage: J(1).omega()
            JackP[]

        The forgotten symmetric functions are the images of the monomial
        symmetric functions under omega::

            sage: Sym = SymmetricFunctions(ZZ)
            sage: m = Sym.m()
            sage: f = Sym.f()
            sage: all( f(lam) == m(lam).omega() for lam in Partitions(3) )
            True
            sage: all( m(lam) == f(lam).omega() for lam in Partitions(3) )
            True
        """
        parent = self.parent()
        s = parent.realization_of().schur()
        return parent(s(self).omega())

    omega_involution = omega

    def theta(self,a):
        r"""
        Return the image of ``self`` under the theta endomorphism which sends
        `p_k` to `a \cdot p_k` for every positive integer `k`.

        In general, this is well-defined outside of the powersum basis only
        if the base ring is a `\QQ`-algebra.

        INPUT:

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
        res = p_self.map_item(lambda m,c: (m, c * a**len(m)))
        return self.parent()(res)

    def theta_qt(self, q=None, t=None):
        r"""
        Return the image of ``self`` under the `q,t`-deformed theta
        endomorphism which sends `p_k` to `\frac{1-q^k}{1-t^k} \cdot p_k`
        for all positive integers `k`.

        In general, this is well-defined outside of the powersum basis only
        if the base ring is a `\QQ`-algebra.

        INPUT:

        - ``q``, ``t`` -- parameters (default: ``None``, in which case 'q'
          and 't' are used)

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
        one = BR.one()
        if not t:
            res = p._from_dict({m: BR(prod(one - q**k for k in m) * c)
                                for m,c in p_self})
        else:
            res = p._from_dict({m: BR(prod((one-q**k) / (one-t**k) for k in m)*c)
                                for m,c in p_self})
        return parent(res)

    def omega_qt(self, q=None, t=None):
        r"""
        Return the image of ``self`` under the `q,t`-deformed omega
        automorphism which sends `p_k` to
        `(-1)^{k-1} \cdot \frac{1-q^k}{1-t^k} \cdot p_k` for all positive
        integers `k`.

        In general, this is well-defined outside of the powersum basis only
        if the base ring is a `\QQ`-algebra.

        If `q = t`, then this is the omega automorphism (:meth:`omega`).

        INPUT:

        - ``q``, ``t`` -- parameters (default: ``None``, in which case
          ``'q'`` and ``'t'`` are used)

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
            -(2976/(q^5-q^3-q^2+1))*p[3, 2]
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
            ((-t^3+t^2+t-1)/(-q^3+q^2+q-1))*McdH[2]
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
        one = BR.one()
        if not t:
            res = p._from_dict({m: c * (-one)**(sum(m)-len(m))
                                     * BR(prod(one-q**i for i in m))
                                for m,c in p_self})
        else:
            res = p._from_dict({m: c * (-one)**(sum(m)-len(m))
                                     * BR(prod((one-q**i) / (one-t**i)
                                               for i in m))
                                for m,c in p_self})
        return parent(res)

    def itensor(self, x):
        r"""
        Return the internal (tensor) product of ``self`` and ``x`` in the
        basis of ``self``.

        The internal tensor product can be defined as the linear extension
        of the definition on power sums
        `p_{\lambda} \ast p_{\mu} = \delta_{\lambda,\mu} z_{\lambda}
        p_{\lambda}`, where `z_{\lambda} = (1^{r_1} r_1!) (2^{r_2} r_2!)
        \cdots` for `\lambda = (1^{r_1} 2^{r_2} \cdots )` and where `\ast`
        denotes the internal tensor product.
        The internal tensor product is also known as the Kronecker product,
        or as the second multiplication on the ring of symmetric functions.

        Note that the internal product of any two homogeneous symmetric
        functions of equal degrees is a homogeneous symmetric function of the
        same degree. On the other hand, the internal product of two homogeneous
        symmetric functions of distinct degrees is `0`.

        .. NOTE::

            The internal product is sometimes referred to as "inner product"
            in the literature, but unfortunately this name is shared by a
            different operation, namely the Hall inner product
            (see :meth:`scalar`).

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the
          same base ring as ``self``

        OUTPUT:

        - the internal product of ``self`` with ``x`` (an element of the
          ring of symmetric functions in the same basis as ``self``)

        The methods :meth:`itensor`, :meth:`internal_product`,
        :meth:`kronecker_product`, :meth:`inner_tensor` are all
        synonyms.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([2,1])
            sage: b = s([3])
            sage: a.itensor(b)
            s[2, 1]
            sage: c = s([3,2,1])
            sage: c.itensor(c)
            s[1, 1, 1, 1, 1, 1] + 2*s[2, 1, 1, 1, 1] + 3*s[2, 2, 1, 1] + 2*s[2, 2, 2]
             + 4*s[3, 1, 1, 1] + 5*s[3, 2, 1] + 2*s[3, 3] + 4*s[4, 1, 1]
             + 3*s[4, 2] + 2*s[5, 1] + s[6]

        There are few quantitative results pertaining to Kronecker products
        in general, which makes their computation so difficult. Let us test
        a few of them in different bases.

        The Kronecker product of any homogeneous symmetric function `f` of
        degree `n` with the `n`-th complete homogeneous symmetric function
        ``h[n]`` (a.k.a. ``s[n]``) is `f`::

            sage: h = SymmetricFunctions(ZZ).h()
            sage: all( h([5]).itensor(h(p)) == h(p) for p in Partitions(5) )
            True

        The Kronecker product of a Schur function `s_{\lambda}` with the `n`-th
        elementary symmetric function ``e[n]``, where `n = \left| \lambda
        \right|`, is `s_{\lambda'}` (where `\lambda'` is the conjugate
        partition of `\lambda`)::

            sage: F = CyclotomicField(12)
            sage: s = SymmetricFunctions(F).s()
            sage: e = SymmetricFunctions(F).e()
            sage: all( e([5]).itensor(s(p)) == s(p.conjugate()) for p in Partitions(5) )
            True

        The Kronecker product is commutative::

            sage: e = SymmetricFunctions(FiniteField(19)).e()
            sage: m = SymmetricFunctions(FiniteField(19)).m()
            sage: all( all( e(p).itensor(m(q)) == m(q).itensor(e(p)) for q in Partitions(4) )
            ....:      for p in Partitions(4) )
            True

            sage: F = FractionField(QQ['q','t'])
            sage: mq = SymmetricFunctions(F).macdonald().Q()
            sage: mh = SymmetricFunctions(F).macdonald().H()
            sage: all( all( mq(p).itensor(mh(r)) == mh(r).itensor(mq(p))   # long time
            ....:           for r in Partitions(4) )
            ....:      for p in Partitions(3) )
            True

        Let us check (on examples) Proposition 5.2 of Gelfand, Krob, Lascoux, Leclerc,
        Retakh, Thibon, "Noncommutative symmetric functions", :arXiv:`hep-th/9407124`, for
        `r = 2`::

            sage: e = SymmetricFunctions(FiniteField(29)).e()
            sage: s = SymmetricFunctions(FiniteField(29)).s()
            sage: m = SymmetricFunctions(FiniteField(29)).m()
            sage: def tensor_copr(u, v, w):  # computes \mu ((u \otimes v) * \Delta(w)) with
            ....:                            # * meaning Kronecker product and \mu meaning the
            ....:                            # usual multiplication.
            ....:     result = w.parent().zero()
            ....:     for partition_pair, coeff in w.coproduct():
            ....:         result += coeff * w.parent()(u).itensor(partition_pair[0]) * w.parent()(v).itensor(partition_pair[1])
            ....:     return result
            sage: all( all( all( tensor_copr(e[u], s[v], m[w])   # long time
            ....:                == (e[u] * s[v]).itensor(m[w])
            ....:                for w in Partitions(5) )
            ....:           for v in Partitions(2) )
            ....:      for u in Partitions(3) )
            True

        Some examples from Briand, Orellana, Rosas, "The stability of the Kronecker
        products of Schur functions." :arXiv:`0907.4652`::

            sage: s = SymmetricFunctions(ZZ).s()
            sage: s[2,2].itensor(s[2,2])
            s[1, 1, 1, 1] + s[2, 2] + s[4]
            sage: s[3,2].itensor(s[3,2])
            s[2, 1, 1, 1] + s[2, 2, 1] + s[3, 1, 1] + s[3, 2] + s[4, 1] + s[5]
            sage: s[4,2].itensor(s[4,2])
            s[2, 2, 2] + s[3, 1, 1, 1] + 2*s[3, 2, 1] + s[4, 1, 1] + 2*s[4, 2] + s[5, 1] + s[6]

        An example from p. 220 of Thibon, "Hopf algebras of symmetric functions
        and tensor products of symmetric group representations", International
        Journal of Algebra and Computation, 1991::

            sage: s = SymmetricFunctions(QQbar).s()
            sage: s[2,1].itensor(s[2,1])
            s[1, 1, 1] + s[2, 1] + s[3]

        TESTS::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([8,8])
            sage: a.itensor(a) # long time
            s[4, 4, 4, 4] + s[5, 5, 3, 3] + s[5, 5, 5, 1] + s[6, 4, 4, 2]
             + s[6, 6, 2, 2] + s[6, 6, 4] + s[7, 3, 3, 3] + s[7, 5, 3, 1]
             + s[7, 7, 1, 1] + s[8, 4, 2, 2] + s[8, 4, 4] + s[8, 6, 2]
             + s[8, 8] + s[9, 3, 3, 1] + s[9, 5, 1, 1] + s[10, 2, 2, 2]
             + s[10, 4, 2] + s[10, 6] + s[11, 3, 1, 1] + s[12, 2, 2]
             + s[12, 4] + s[13, 1, 1, 1] + s[14, 2] + s[16]
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

        Same over the ring of integers::

            sage: s = SymmetricFunctions(ZZ).s()
            sage: a = s([8,8])
            sage: a.itensor(a) # long time
            s[4, 4, 4, 4] + s[5, 5, 3, 3] + s[5, 5, 5, 1] + s[6, 4, 4, 2]
             + s[6, 6, 2, 2] + s[6, 6, 4] + s[7, 3, 3, 3] + s[7, 5, 3, 1]
             + s[7, 7, 1, 1] + s[8, 4, 2, 2] + s[8, 4, 4] + s[8, 6, 2]
             + s[8, 8] + s[9, 3, 3, 1] + s[9, 5, 1, 1] + s[10, 2, 2, 2]
             + s[10, 4, 2] + s[10, 6] + s[11, 3, 1, 1] + s[12, 2, 2]
             + s[12, 4] + s[13, 1, 1, 1] + s[14, 2] + s[16]
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

        Theorem 2.1 in Bessenrodt, van Willigenburg, :arXiv:`1105.3170v2`::

            sage: s = SymmetricFunctions(ZZ).s()
            sage: all( all( max( r[0] for r in s(p).itensor(s(q)).monomial_coefficients().keys() )
            ....:           == sum( min(p[i], q.get_part(i)) for i in range(len(p)) )
            ....:           for p in Partitions(4) )
            ....:      for q in Partitions(4) )
            True
            sage: all( all( max( len(r) for r in s(p).itensor(s(q)).monomial_coefficients().keys() )
            ....:           == sum( min(p[i], q.conjugate().get_part(i)) for i in range(len(p)) )
            ....:           for p in Partitions(4) )
            ....:      for q in Partitions(4) )
            True

        Check that the basis and ground ring of ``self`` are preserved::

            sage: F = CyclotomicField(12)
            sage: s = SymmetricFunctions(F).s()
            sage: e = SymmetricFunctions(F).e()
            sage: e[3].itensor(s[3])
            e[3]
            sage: s[3].itensor(e[3])
            s[1, 1, 1]
            sage: parent(e[3].itensor(s[3]))
            Symmetric Functions over Cyclotomic Field of order 12 and degree 4 in the elementary basis
            sage: parent(s[3].itensor(e[3]))
            Symmetric Functions over Cyclotomic Field of order 12 and degree 4 in the Schur basis

        .. NOTE::

            The currently existing implementation of this function is
            technically unsatisfactory. It distinguishes the case when the
            base ring is a `\QQ`-algebra (in which case the Kronecker product
            can be easily computed using the power sum basis) from the case
            where it isn't. In the latter, it does a computation using
            universal coefficients, again distinguishing the case when it is
            able to compute the "corresponding" basis of the symmetric function
            algebra over `\QQ` (using the ``corresponding_basis_over`` hack)
            from the case when it isn't (in which case it transforms everything
            into the Schur basis, which is slow).
        """
        parent = self.parent()
        if parent.has_coerce_map_from(QQ):
            # Convert both self and x to the p basis
            p = parent.realization_of().power()
            f = lambda part1, part2: zee(part1)*p(part1)
            return parent(p._apply_multi_module_morphism(p(self),p(x),f,orthogonal=True))
        else:
            # comp_parent is the parent that is going to be used for
            # computations. In most cases it will just be parent.
            # Similarly for comp_self and comp_x.
            comp_parent = parent
            comp_self = self
            # Now let's try to find out what basis self is in, and
            # construct the corresponding basis of symmetric functions
            # over QQ.
            corresponding_parent_over_QQ = parent.corresponding_basis_over(QQ)
            if corresponding_parent_over_QQ is None:
                # This is the case where the corresponding basis
                # over QQ cannot be found. This can have two reasons:
                # Either the basis depends on variables (like the
                # Macdonald symmetric functions), or its basis_name()
                # is not identical to the name of the method on
                # SymmetricFunctions(QQ) that builds it. Either way,
                # give up looking for the corresponding parent, and
                # transform everything into the Schur basis (very
                # slow!) instead.
                comp_parent = parent.realization_of().schur()
                comp_self = comp_parent(self)
                from sage.combinat.sf.sf import SymmetricFunctions
                corresponding_parent_over_QQ = SymmetricFunctions(QQ).schur()
            comp_x = comp_parent(x)    # For simplicity, let self and x be in the same basis.
            result = comp_parent.zero()
            for lam, a in comp_self:
                # lam is a partition, a is an element of the base ring.
                for mu, b in comp_x:
                    # mu is a partition, b is an element of the base ring.
                    lam_star_mu = corresponding_parent_over_QQ(lam).itensor(corresponding_parent_over_QQ(mu))
                    # lam_star_mu is now a symmetric function over QQ.
                    for nu, c in lam_star_mu:
                        # nu is a partition, c is an element of QQ.
                        result += a * b * comp_parent.base_ring()(c) * comp_parent(nu)
            return parent(result)    # just in case comp_parent != parent.

    internal_product = itensor
    kronecker_product = itensor
    inner_tensor = itensor

    def reduced_kronecker_product(self, x):
        r"""
        Return the reduced Kronecker product of ``self`` and ``x`` in the
        basis of ``self``.

        The reduced Kronecker product is a bilinear map mapping two
        symmetric functions to another, not necessarily preserving degree.
        It can be defined as follows: Let `*` denote the Kronecker product
        (:meth:`itensor`) on the space of symmetric functions. For any
        partitions `\alpha`, `\beta`, `\gamma`, let
        `g^{\gamma}_{\alpha, \beta}` denote the coefficient of the Schur
        function `s_{\gamma}` in the Kronecker product
        `s_{\alpha} * s_{\beta}` (this is called a Kronecker coefficient).
        For every partition
        `\lambda = (\lambda_1, \lambda_2, \lambda_3, \ldots)`
        and every integer `n > \left| \lambda \right| + \lambda_1`, let
        `\lambda[n]` denote the `n`-completion of `\lambda` (this is the
        partition
        `(n - \left| \lambda \right|, \lambda_1, \lambda_2, \lambda_3, \ldots)`;
        see :meth:`~sage.combinat.partition.Partition.t_completion`).
        Then, Theorem 1.2 of [BOR09]_ shows that for any partitions
        `\alpha` and `\beta` and every integer
        `n \geq \left|\alpha\right| + \left|\beta\right| + \alpha_1 + \beta_1`,
        we can write the Kronecker product `s_{\alpha[n]} * s_{\beta[n]}`
        in the form

        .. MATH::

            s_{\alpha[n]} * s_{\beta[n]} = \sum_{\gamma} g^{\gamma[n]}_{\alpha[n], \beta[n]} s_{\gamma[n]}

        with `\gamma` ranging over all partitions. The
        coefficients `g^{\gamma[n]}_{\alpha[n], \beta[n]}`
        are independent on `n`. These coefficients
        `g^{\gamma[n]}_{\alpha[n], \beta[n]}` are denoted by
        `\overline{g}^{\gamma}_{\alpha, \beta}`, and the symmetric
        function

        .. MATH::

            \sum_{\gamma} \overline{g}^{\gamma}_{\alpha, \beta} s_{\gamma}

        is said to be the *reduced Kronecker product* of `s_{\alpha}` and
        `s_{\beta}`. By bilinearity, this extends to a definition of a
        reduced Kronecker product of any two symmetric functions.

        The definition of the reduced Kronecker product goes back to
        Murnaghan, and has recently been studied in [BOR09]_, [BdVO12]_
        and other places (our notation
        `\overline{g}^{\gamma}_{\alpha, \beta}` appears in these two
        sources).

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the
          same base ring as ``self``

        OUTPUT:

        - the reduced Kronecker product of ``self`` with ``x`` (an element
          of the ring of symmetric functions in the same basis as
          ``self``)

        EXAMPLES:

        The example from page 2 of [BOR09]_::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.schur()
            sage: s[2].reduced_kronecker_product(s[2])
            s[] + s[1] + s[1, 1] + s[1, 1, 1] + 2*s[2] + 2*s[2, 1] + s[2, 2] + s[3] + s[3, 1] + s[4]

        Taking the reduced Kronecker product with `1 = s_{\empty}` is the
        identity map on the ring of symmetric functions::

            sage: all( s[Partition([])].reduced_kronecker_product(s[lam])
            ....:      == s[lam] for i in range(4)
            ....:      for lam in Partitions(i) )
            True

        While reduced Kronecker products are hard to compute in general,
        there is a rule for taking reduced Kronecker products with
        `s_1`. Namely, for every partition `\lambda`, the reduced
        Kronecker product of `s_{\lambda}` with `s_1` is
        `\sum_{\mu} a_{\mu} s_{\mu}`, where the sum runs over all
        partitions `\mu`, and the coefficient `a_{\mu}` is defined as the
        number of ways to obtain `\mu` from `\lambda` by one of the
        following three operations:

        - Add an addable cell
          (:meth:`~sage.combinat.partition.Partition.addable_cells`) to
          `\lambda`.
        - Remove a removable cell
          (:meth:`~sage.combinat.partition.Partition.removable_cells`)
          from `\lambda`.
        - First remove a removable cell from `\lambda`, then add an
          addable cell to the resulting Young diagram.

        This is, in fact, Proposition 5.15 of [CO10]_ in an elementary
        wording. We check this for partitions of size `\leq 4`::

            sage: def mults1(lam):
            ....:     # Reduced Kronecker multiplication by s[1], according
            ....:     # to [CO10]_.
            ....:     res = s.zero()
            ....:     for mu in lam.up_list():
            ....:         res += s(mu)
            ....:     for mu in lam.down_list():
            ....:         res += s(mu)
            ....:         for nu in mu.up_list():
            ....:             res += s(nu)
            ....:     return res
            sage: all( mults1(lam) == s[1].reduced_kronecker_product(s[lam])
            ....:      for i in range(5) for lam in Partitions(i) )
            True

        Here is the example on page 3 of Christian Gutschwager's
        :arxiv:`0912.4411v3`::

            sage: s[1,1].reduced_kronecker_product(s[2])
            s[1] + 2*s[1, 1] + s[1, 1, 1] + s[2] + 2*s[2, 1] + s[2, 1, 1] + s[3] + s[3, 1]

        Example 39 from F. D. Murnaghan, "The analysis of the Kronecker
        product of irreducible representations of the symmetric group",
        American Journal of Mathematics, Vol. 60, No. 3, Jul. 1938::

            sage: s[3].reduced_kronecker_product(s[2,1])
            s[1] + 2*s[1, 1] + 2*s[1, 1, 1] + s[1, 1, 1, 1] + 2*s[2] + 5*s[2, 1] + 4*s[2, 1, 1]
            + s[2, 1, 1, 1] + 3*s[2, 2] + 2*s[2, 2, 1] + 2*s[3] + 5*s[3, 1] + 3*s[3, 1, 1]
            + 3*s[3, 2] + s[3, 2, 1] + 2*s[4] + 3*s[4, 1] + s[4, 1, 1] + s[4, 2] + s[5]
            + s[5, 1]

        TESTS::

            sage: h = SymmetricFunctions(QQ).h()
            sage: (2*h([])).reduced_kronecker_product(3*h([]))
            6*h[]

        Different bases and base rings::

            sage: h = SymmetricFunctions(ZZ).h()
            sage: e = SymmetricFunctions(ZZ).e()
            sage: h(e[2].reduced_kronecker_product(h[2]))
            h[1] + 2*h[1, 1] + h[1, 1, 1] - h[2] + h[2, 1, 1] - h[2, 2]

            sage: F = CyclotomicField(12)
            sage: s = SymmetricFunctions(F).s()
            sage: e = SymmetricFunctions(F).e()
            sage: v = e[2].reduced_kronecker_product(e[2]); v
            e[] + e[1] + 2*e[1, 1] + e[1, 1, 1] + (-1)*e[2] + e[2, 2]
            sage: parent(v)
            Symmetric Functions over Cyclotomic Field of order 12 and degree 4 in the elementary basis

            sage: s = SymmetricFunctions(ZZ).s()
            sage: v = s[1].reduced_kronecker_product(s[1]); parent(v)
            Symmetric Functions over Integer Ring in the Schur basis

        .. TODO::

            This implementation of the reduced Kronecker product is
            painfully slow.
        """
        parent = self.parent()
        comp_parent = parent.realization_of().schur()
        comp_self = comp_parent(self)
        comp_x = comp_parent(x)
        # Now, comp_self and comp_x are the same as self and x, but in the
        # Schur basis, which we call comp_parent.
        schur_Q = comp_parent.corresponding_basis_over(QQ)
        # schur_Q is the Schur basis of the symmetric functions over QQ.
        result = comp_parent.zero()
        for lam, a in comp_self:
            # lam is a partition, a is an element of the base ring.
            lam_list = lam._list
            if not lam_list:
                # Special handling for the empty partition. The reduced
                # Kronecker product of 1 with any symmetric function f is
                # f.
                result += a * comp_x
                continue
            sum_lam = sum(lam_list)
            for mu, b in comp_x:
                # mu is a partition, b is an element of the base ring.
                mu_list = mu._list
                if not mu_list:
                    # Special handling for the empty partition.
                    result += a * b * comp_parent(lam)
                    continue
                # Now, both lam and mu are nonempty.
                sum_mu = sum(mu_list)
                stab = lam_list[0] + mu_list[0] + sum_lam + sum_mu
                s_lam_stabilized = schur_Q(_Partitions([stab - sum_lam] + lam_list))
                s_mu_stabilized = schur_Q(_Partitions([stab - sum_mu] + mu_list))
                lam_star_mu = s_lam_stabilized.itensor(s_mu_stabilized)
                # lam_star_mu is now a symmetric function over QQ.
                for nu, c in lam_star_mu:
                    # nu is a partition of the integer stab, c is an element of QQ.
                    nu_unstabilized = _Partitions(nu[1:])
                    result += a * b * comp_parent.base_ring()(c) \
                                    * comp_parent(nu_unstabilized)
        return parent(result)

    def left_padded_kronecker_product(self, x):
        r"""
        Return the left-padded Kronecker product of ``self`` and ``x`` in
        the basis of ``self``.

        The left-padded Kronecker product is a bilinear map mapping two
        symmetric functions to another, not necessarily preserving degree.
        It can be defined as follows: Let `*` denote the Kronecker product
        (:meth:`itensor`) on the space of symmetric functions. For any
        partitions `\alpha`, `\beta`, `\gamma`, let
        `g^{\gamma}_{\alpha, \beta}` denote the coefficient of the
        complete homogeneous symmetric function `h_{\gamma}` in the
        Kronecker product `h_{\alpha} * h_{\beta}`.
        For every partition
        `\lambda = (\lambda_1, \lambda_2, \lambda_3, \ldots)`
        and every integer `n > \left| \lambda \right| + \lambda_1`, let
        `\lambda[n]` denote the `n`-completion of `\lambda` (this is the
        partition
        `(n - \left| \lambda \right|, \lambda_1, \lambda_2, \lambda_3, \ldots)`;
        see :meth:`~sage.combinat.partition.Partition.t_completion`).
        Then, for any partitions `\alpha` and `\beta` and every integer
        `n \geq \left|\alpha\right| + \left|\beta\right| + \alpha_1 + \beta_1`,
        we can write the Kronecker product `h_{\alpha[n]} * h_{\beta[n]}`
        in the form

        .. MATH::

            h_{\alpha[n]} * h_{\beta[n]} = \sum_{\gamma}
            g^{\gamma[n]}_{\alpha[n], \beta[n]} h_{\gamma[n]}

        with `\gamma` ranging over all partitions. The
        coefficients `g^{\gamma[n]}_{\alpha[n], \beta[n]}`
        are independent on `n`. These coefficients
        `g^{\gamma[n]}_{\alpha[n], \beta[n]}` are denoted by
        `\overline{g}^{\gamma}_{\alpha, \beta}`, and the symmetric
        function

        .. MATH::

            \sum_{\gamma} \overline{g}^{\gamma}_{\alpha, \beta} h_{\gamma}

        is said to be the *left-padded Kronecker product* of `h_{\alpha}`
        and `h_{\beta}`. By bilinearity, this extends to a definition of a
        left-padded Kronecker product of any two symmetric functions.

        This notion of left-padded Kronecker product can be lifted to the
        non-commutative symmetric functions
        (:meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmeticFunctions.Bases.ElementMethods.left_padded_kronecker_product`).

        .. WARNING::

            Don't mistake this product for the reduced Kronecker product
            (:meth:`reduced_kronecker_product`), which uses the Schur
            functions instead of the complete homogeneous functions in
            its definition.

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the
          same base ring as ``self``

        OUTPUT:

        - the left-padded Kronecker product of ``self`` with ``x`` (an
          element of the ring of symmetric functions in the same basis
          as ``self``)

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: h = Sym.h()
            sage: h[2,1].left_padded_kronecker_product(h[3])
            h[1, 1, 1, 1] + h[2, 1] + h[2, 1, 1] + h[2, 1, 1, 1] + h[2, 2, 1] + h[3, 2, 1]
            sage: h[2,1].left_padded_kronecker_product(h[1])
            h[1, 1, 1] + h[2, 1] + h[2, 1, 1]
            sage: h[1].left_padded_kronecker_product(h[2,1])
            h[1, 1, 1] + h[2, 1] + h[2, 1, 1]
            sage: h[1,1].left_padded_kronecker_product(h[2])
            h[1, 1] + 2*h[1, 1, 1] + h[2, 1, 1]
            sage: h[1].left_padded_kronecker_product(h[2,1,1])
            h[1, 1, 1, 1] + 2*h[2, 1, 1] + h[2, 1, 1, 1]
            sage: h[2].left_padded_kronecker_product(h[3])
            h[2, 1] + h[2, 1, 1] + h[3, 2]

        Taking the left-padded Kronecker product with `1 = h_{\empty}` is
        the identity map on the ring of symmetric functions::

            sage: all( h[Partition([])].left_padded_kronecker_product(h[lam])
            ....:      == h[lam] for i in range(4)
            ....:      for lam in Partitions(i) )
            True

        Here is a rule for the left-padded Kronecker product of `h_1`
        (this is the same as `h_{(1)}`) with any complete homogeneous
        function: Let `\lambda` be a partition. Then, the left-padded
        Kronecker product of `h_1` and `h_{\lambda}` is
        `\sum_{\mu} a_{\mu} h_{\mu}`, where the sum runs over all
        partitions `\mu`, and the coefficient `a_{\mu}` is defined as the
        number of ways to obtain `\mu` from `\lambda` by one of the
        following two operations:

        - Insert a `1` into `\lambda`.
        - Subtract `1` from one of the entries of `\lambda` (and remove
          the entry if it thus becomes `0`), and insert a `1` into
          `\lambda`.

        We check this for partitions of size `\leq 4`::

            sage: def mults1(I):
            ....:     # Left-padded Kronecker multiplication by h[1].
            ....:     res = h[I[:] + [1]]
            ....:     for k in range(len(I)):
            ....:         I2 = I[:]
            ....:         if I2[k] == 1:
            ....:             I2 = I2[:k] + I2[k+1:]
            ....:         else:
            ....:             I2[k] -= 1
            ....:         res += h[sorted(I2 + [1], reverse=True)]
            ....:     return res
            sage: all( mults1(I) == h[1].left_padded_kronecker_product(h[I])
            ....:                == h[I].left_padded_kronecker_product(h[1])
            ....:      for i in range(5) for I in Partitions(i) )
            True

        The left-padded Kronecker product is commutative::

            sage: all( h[lam].left_padded_kronecker_product(h[mu])
            ....:      == h[mu].left_padded_kronecker_product(h[lam])
            ....:      for lam in Partitions(3) for mu in Partitions(3) )
            True

        TESTS::

            sage: h = SymmetricFunctions(QQ).h()
            sage: (2*h([])).left_padded_kronecker_product(3*h([]))
            6*h[]

        Different bases and base rings::

            sage: h = SymmetricFunctions(ZZ).h()
            sage: e = SymmetricFunctions(ZZ).e()
            sage: h(e[2].left_padded_kronecker_product(h[2]))
            h[1, 1] + h[1, 1, 1] - h[2] + h[2, 1, 1] - h[2, 2]

            sage: F = CyclotomicField(12)
            sage: s = SymmetricFunctions(F).s()
            sage: e = SymmetricFunctions(F).e()
            sage: v = e[2].left_padded_kronecker_product(e[2]); v
            e[1, 1] + e[1, 1, 1] + (-1)*e[2] + e[2, 2]
            sage: parent(v)
            Symmetric Functions over Cyclotomic Field of order 12 and degree 4 in the elementary basis

            sage: s = SymmetricFunctions(ZZ).s()
            sage: v = s[1].left_padded_kronecker_product(s[1]); parent(v)
            Symmetric Functions over Integer Ring in the Schur basis
        """
        from sage.combinat.composition import Compositions
        _Compositions = Compositions()
        parent = self.parent()
        h = parent.realization_of().h()
        h_self = h(self)
        h_x = h(x)
        # Now, h_self and h_x are the same as self and x, but in the
        # h (=complete homogeneous) basis, which we call h.
        R = self.base_ring()
        from sage.combinat.ncsf_qsym.ncsf import NonCommutativeSymmetricFunctions
        # We lift to the noncommutative symmetric functions.
        S = NonCommutativeSymmetricFunctions(R).S()
        result = h.zero()
        for lam, a in h_self:
            # lam is a partition, a is an element of the base ring.
            if not lam._list:
                # Special handling for the empty partition. The reduced
                # Kronecker product of 1 with any symmetric function f is f.
                result += a * h_x
                continue
            c_lam = _Compositions(lam)
            for mu, b in h_x:
                # mu is a partition, b is an element of the base ring.
                if not mu._list:
                    # Special handling for the empty partition.
                    result += a * b * h(lam)
                    continue
                # Now, both lam and mu are nonempty.
                c_mu = _Compositions(mu)
                result += a * b * S[c_lam].left_padded_kronecker_product(S[c_mu]).to_symmetric_function()
        return parent(result)

    def internal_coproduct(self):
        r"""
        Return the inner coproduct of ``self`` in the basis of ``self``.

        The inner coproduct (also known as the Kronecker coproduct, as the
        internal coproduct, or as the second comultiplication on the ring of
        symmetric functions) is a ring homomorphism `\Delta^\times` from the
        ring of symmetric functions to the tensor product (over the base
        ring) of this ring with itself. It is uniquely characterized by the
        formula

        .. MATH::

            \Delta^{\times}(h_n) = \sum_{\lambda \vdash n} s_{\lambda}
            \otimes s_{\lambda} = \sum_{\lambda \vdash n} h_{\lambda} \otimes
            m_{\lambda} = \sum_{\lambda \vdash n} m_{\lambda} \otimes
            h_{\lambda},

        where `\lambda \vdash n` means `\lambda` is a partition of `n`, and
        `n` is any nonnegative integer. It also satisfies

        .. MATH::

            \Delta^\times (p_n) = p_n \otimes p_n

        for any positive integer `n`. If the base ring is a `\QQ`-algebra, it
        also satisfies

        .. MATH::

            \Delta^{\times}(h_n) = \sum_{\lambda \vdash n} z_{\lambda}^{-1}
            p_{\lambda} \otimes p_{\lambda},

        where

        .. MATH::

            z_{\lambda} = \prod_{i=1}^\infty i^{m_i(\lambda)} m_i(\lambda)!

        with `m_i(\lambda)` meaning the number of appearances of `i`
        in `\lambda` (see :meth:`~sage.combinat.sf.sfa.zee`).

        The method :meth:`kronecker_coproduct` is a synonym of
        :meth:`internal_coproduct`.

        EXAMPLES::

            sage: s = SymmetricFunctions(ZZ).s()
            sage: a = s([2,1])
            sage: a.internal_coproduct()
            s[1, 1, 1] # s[2, 1] + s[2, 1] # s[1, 1, 1] + s[2, 1] # s[2, 1] + s[2, 1] # s[3] + s[3] # s[2, 1]

            sage: e = SymmetricFunctions(QQ).e()
            sage: b = e([2])
            sage: b.internal_coproduct()
            e[1, 1] # e[2] + e[2] # e[1, 1] - 2*e[2] # e[2]

        The internal coproduct is adjoint to the internal product with respect
        to the Hall inner product: Any three symmetric functions `f`, `g` and
        `h` satisfy `\langle f * g, h \rangle = \sum_i \langle f, h^{\prime}_i
        \rangle \langle g, h^{\prime\prime}_i \rangle`, where we write
        `\Delta^{\times}(h)` as `\sum_i h^{\prime}_i \otimes
        h^{\prime\prime}_i`. Let us check this in degree `4`::

            sage: e = SymmetricFunctions(FiniteField(29)).e()
            sage: s = SymmetricFunctions(FiniteField(29)).s()
            sage: m = SymmetricFunctions(FiniteField(29)).m()
            sage: def tensor_incopr(f, g, h):  # computes \sum_i \left< f, h'_i \right> \left< g, h''_i \right>
            ....:     result = h.base_ring().zero()
            ....:     for partition_pair, coeff in h.internal_coproduct():
            ....:         result += coeff * h.parent()(f).scalar(partition_pair[0]) * h.parent()(g).scalar(partition_pair[1])
            ....:     return result
            sage: all( all( all( tensor_incopr(e[u], s[v], m[w]) == (e[u].itensor(s[v])).scalar(m[w])  # long time (10s on sage.math, 2013)
            ....:                for w in Partitions(5) )
            ....:           for v in Partitions(2) )
            ....:      for u in Partitions(3) )
            True

        Let us check the formulas for `\Delta^{\times}(h_n)` and
        `\Delta^{\times}(p_n)` given in the description of this method::

            sage: e = SymmetricFunctions(QQ).e()
            sage: p = SymmetricFunctions(QQ).p()
            sage: h = SymmetricFunctions(QQ).h()
            sage: s = SymmetricFunctions(QQ).s()
            sage: all( s(h([n])).internal_coproduct() == sum([tensor([s(lam), s(lam)]) for lam in Partitions(n)])
            ....:      for n in range(6) )
            True
            sage: all( h([n]).internal_coproduct() == sum([tensor([h(lam), h(m(lam))]) for lam in Partitions(n)])
            ....:      for n in range(6) )
            True
            sage: all( factorial(n) * h([n]).internal_coproduct()
            ....:      == sum([lam.conjugacy_class_size() * tensor([h(p(lam)), h(p(lam))])
            ....:              for lam in Partitions(n)])
            ....:      for n in range(6) )
            True

        TESTS::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s([]).internal_coproduct()
            s[] # s[]
        """
        parent = self.parent()
        h = parent.realization_of().homogeneous()
        s = parent.realization_of().schur()
        from sage.categories.tensor import tensor
        result = tensor([parent.zero(), parent.zero()])
        result_parent = result.parent()
        from sage.misc.cachefunc import cached_function
        @cached_function
        def hnimage(n):
            return result_parent.sum((tensor([parent(s(lam)), parent(s(lam))]) for lam in Partitions(n)))
        for lam, a in h(self):
            result += a * prod((hnimage(i) for i in lam))
        return result

    kronecker_coproduct = internal_coproduct

    def arithmetic_product(self, x):
        r"""
        Return the arithmetic product of ``self`` and ``x`` in the
        basis of ``self``.

        The arithmetic product is a binary operation `\boxdot` on the
        ring of symmetric functions which is bilinear in its two
        arguments and satisfies

        .. MATH::

            p_{\lambda} \boxdot p_{\mu} = \prod\limits_{i \geq 1, j \geq 1}
            p_{\mathrm{lcm}(\lambda_i, \mu_j)}^{\mathrm{gcd}(\lambda_i, \mu_j)}

        for any two partitions `\lambda = (\lambda_1, \lambda_2, \lambda_3,
        \dots )` and `\mu = (\mu_1, \mu_2, \mu_3, \dots )` (where `p_{\nu}`
        denotes the power-sum symmetric function indexed by the partition
        `\nu`, and `p_i` denotes the `i`-th power-sum symmetric function).
        This is enough to define the arithmetic product if the base ring
        is torsion-free as a `\ZZ`-module; for all other cases the
        arithmetic product is uniquely determined by requiring it to be
        functorial in the base ring. See
        http://mathoverflow.net/questions/138148/ for a discussion of
        this arithmetic product.

        If `f` and `g` are two symmetric functions which are homogeneous
        of degrees `a` and `b`, respectively, then `f \boxdot g` is
        homogeneous of degree `ab`.

        The arithmetic product is commutative and associative and has
        unity `e_1 = p_1 = h_1`.

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the
          same base ring as ``self``

        OUTPUT:

        Arithmetic product of ``self`` with ``x``; this is a symmetric
        function over the same base ring as ``self``.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s([2]).arithmetic_product(s([2]))
            s[1, 1, 1, 1] + 2*s[2, 2] + s[4]
            sage: s([2]).arithmetic_product(s([1,1]))
            s[2, 1, 1] + s[3, 1]

        The symmetric function ``e[1]`` is the unity for the arithmetic
        product::

            sage: e = SymmetricFunctions(ZZ).e()
            sage: all( e([1]).arithmetic_product(e(q)) == e(q) for q in Partitions(4) )
            True

        The arithmetic product is commutative::

            sage: e = SymmetricFunctions(FiniteField(19)).e()
            sage: m = SymmetricFunctions(FiniteField(19)).m()
            sage: all( all( e(p).arithmetic_product(m(q)) == m(q).arithmetic_product(e(p))  # long time (26s on sage.math, 2013)
            ....:           for q in Partitions(4) )
            ....:      for p in Partitions(4) )
            True

        .. NOTE::

            The currently existing implementation of this function is
            technically unsatisfactory. It distinguishes the case when the
            base ring is a `\QQ`-algebra (in which case the arithmetic product
            can be easily computed using the power sum basis) from the case
            where it isn't. In the latter, it does a computation using
            universal coefficients, again distinguishing the case when it is
            able to compute the "corresponding" basis of the symmetric function
            algebra over `\QQ` (using the ``corresponding_basis_over`` hack)
            from the case when it isn't (in which case it transforms everything
            into the Schur basis, which is slow).
        """
        # The following code is analogous to the code of itensor, so comments
        # have been removed for brevity.
        parent = self.parent()
        if parent.has_coerce_map_from(QQ):
            from sage.combinat.partition import Partition
            from sage.arith.all import gcd, lcm
            from itertools import product, repeat, chain
            p = parent.realization_of().power()
            def f(lam, mu):
                # This is the map sending two partitions lam and mu to the
                # arithmetic product p[lam] \boxdot p[mu].
                # Code shamelessly stolen from Andrew Gainer-Dewar, trac #14542.
                term_iterable = chain.from_iterable( repeat(lcm(pair), times=gcd(pair))
                                                     for pair in product(lam, mu) )
                term_list = sorted(term_iterable, reverse=True)
                res = Partition(term_list)
                return p(res)
            return parent(p._apply_multi_module_morphism(p(self),p(x),f))
        comp_parent = parent
        comp_self = self
        corresponding_parent_over_QQ = parent.corresponding_basis_over(QQ)
        if corresponding_parent_over_QQ is None:
            comp_parent = parent.realization_of().schur()
            comp_self = comp_parent(self)
            from sage.combinat.sf.sf import SymmetricFunctions
            corresponding_parent_over_QQ = SymmetricFunctions(QQ).schur()
        comp_x = comp_parent(x)
        result = comp_parent.zero()
        for lam, a in comp_self:
            for mu, b in comp_x:
                lam_star_mu = corresponding_parent_over_QQ(lam).arithmetic_product(corresponding_parent_over_QQ(mu))
                for nu, c in lam_star_mu:
                    result += a * b * comp_parent.base_ring()(c) * comp_parent(nu)
        return parent(result)

    def nabla(self, q=None, t=None, power=1):
        r"""
        Return the value of the nabla operator applied to ``self``.

        The eigenvectors of the nabla operator are the Macdonald polynomials in
        the Ht basis.

        If the parameter ``power`` is an integer then it calculates
        nabla to that integer.  The default value of ``power`` is 1.

        INPUT:

        - ``q``, ``t`` -- optional parameters (default: ``None``, in which
          case ``q`` and ``t`` are used)
        - ``power`` -- (default: ``1``) an integer indicating how many times to
          apply the operator `\nabla`.  Negative values of ``power``
          indicate powers of `\nabla^{-1}`.

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
            ((-q-t)/(q^2*t^2))*s[2, 1] + ((q^2+q*t+t^2)/(-q^3*t^3))*s[3]
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
        Return the standard scalar product between ``self`` and ``x``.

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the
          same base ring as ``self``

        - ``zee`` -- an optional function on partitions giving
          the value for the scalar product between `p_{\mu}` and `p_{\mu}`
          (default is to use the standard :meth:`~sage.combinat.sf.sfa.zee` function)

        This is the default implementation that converts both ``self`` and
        ``x`` into either Schur functions (if ``zee`` is not specified) or
        power-sum functions (if ``zee`` is specified) and performs the scalar
        product in that basis.

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

        Over the integers, too (as long as ``zee`` is not set)::

            sage: Sym = SymmetricFunctions(ZZ)
            sage: m = Sym.m()
            sage: m([2]).scalar(m([2]))
            2
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
        Returns the `q,t`-deformed standard Hall-Littlewood scalar product of
        ``self`` and ``x``.

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the same
          base ring as ``self``

        - ``q``, ``t`` -- parameters (default: ``None`` in which case ``q``
          and ``t`` are used)

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
            1/3*(x^3 - 1)/(y^3 - 1) + 2/3*(x - 1)^3/(y - 1)^3
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
        Return the `t`-deformed standard Hall-Littlewood scalar product of
        ``self`` and ``x``.

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the same
          base ring as ``self``

        - ``t`` -- parameter (default: ``None``, in which case ``t`` is used)

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: a = s([2,1])
            sage: sp = a.scalar_t(a); sp
            (-t^2 - 1)/(t^5 - 2*t^4 + t^3 - t^2 + 2*t - 1)
            sage: sp.parent()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        return self.scalar_qt( x, q=self.base_ring().zero(), t=t )

    scalar_hl = scalar_t

    def scalar_jack(self, x, t=None):
        r"""
        Return the Jack-scalar product beween ``self`` and ``x``.

        This scalar product is defined so that the power sum elements
        `p_{\mu}` are orthogonal and `\langle p_{\mu}, p_{\mu} \rangle =
        z_{\mu} t^{\ell(\mu)}`, where `\ell(\mu)` denotes the length of
        `\mu`.

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the
          same base ring as ``self``
        - ``t`` -- an optional parameter (default: ``None`` in which
          case ``t`` is used)

        EXAMPLES::

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
        if t is None:
            if hasattr(parent,"t"):
                t = self.parent().t
            else:
                t = QQ['t'].gen()
        zee = lambda part: part.centralizer_size()*t**part.length()
        return self.scalar(x, zee)

    def derivative_with_respect_to_p1(self, n=1):
        r"""
        Return the symmetric function obtained by taking the derivative of
        ``self`` with respect to the power-sum symmetric function `p_1`
        when the expansion of ``self`` in the power-sum basis is considered
        as a polynomial in `p_k`'s (with `k \geq 1`).

        This is the same as skewing ``self`` by the first power-sum symmetric
        function `p_1`.

        INPUT:

        - ``n`` -- (default: 1) nonnegative integer which determines
          which power of the derivative is taken

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

        Let us check that taking the derivative with respect to ``p[1]``
        is equivalent to skewing by ``p[1]``::

            sage: p1 = s([1])
            sage: all( s(lam).derivative_with_respect_to_p1()
            ....:      == s(lam).skew_by(p1) for lam in Partitions(4) )
            True
        """
        p = self.parent().realization_of().power()
        res = p(self)
        for i in range(n):
            res = res._derivative_with_respect_to_p1()
        return self.parent()(res)

    def frobenius(self, n):
        r"""
        Return the image of the symmetric function ``self`` under the
        `n`-th Frobenius operator.

        The `n`-th Frobenius operator `\mathbf{f}_n` is defined to be the
        map from the ring of symmetric functions to itself that sends
        every symmetric function `P(x_1, x_2, x_3, \ldots)` to
        `P(x_1^n, x_2^n, x_3^n, \ldots)`. This operator `\mathbf{f}_n`
        is a Hopf algebra endomorphism, and satisfies

        .. MATH::

            \mathbf{f}_n m_{(\lambda_1, \lambda_2, \lambda_3, \ldots)} =
            m_{(n\lambda_1, n\lambda_2, n\lambda_3, \ldots)}

        for every partition `(\lambda_1, \lambda_2, \lambda_3, \ldots)`
        (where `m` means the monomial basis). Moreover,
        `\mathbf{f}_n (p_r) = p_{nr}` for every positive integer `r` (where
        `p_k` denotes the `k`-th powersum symmetric function).

        The `n`-th Frobenius operator is also called the `n`-th
        Frobenius endomorphism. It is not related to the Frobenius map
        which connects the ring of symmetric functions with the
        representation theory of the symmetric group.

        The `n`-th Frobenius operator is also the `n`-th Adams operator
        of the `\Lambda`-ring of symmetric functions over the integers.

        The `n`-th Frobenius operator can also be described via plethysm:
        Every symmetric function `P` satisfies
        `\mathbf{f}_n(P) = p_n \circ P = P \circ p_n`,
        where `p_n` is the `n`-th powersum symmetric function, and `\circ`
        denotes (outer) plethysm.

        INPUT:

        - ``n`` -- a positive integer

        OUTPUT:

        The result of applying the `n`-th Frobenius operator (on the ring of
        symmetric functions) to ``self``.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(ZZ)
            sage: p = Sym.p()
            sage: h = Sym.h()
            sage: s = Sym.s()
            sage: m = Sym.m()
            sage: s[3].frobenius(2)
            -s[3, 3] + s[4, 2] - s[5, 1] + s[6]
            sage: m[4,2,1].frobenius(3)
            m[12, 6, 3]
            sage: p[4,2,1].frobenius(3)
            p[12, 6, 3]
            sage: h[4].frobenius(2)
            h[4, 4] - 2*h[5, 3] + 2*h[6, 2] - 2*h[7, 1] + 2*h[8]

        The Frobenius endomorphisms are multiplicative::

            sage: all( all( s(lam).frobenius(3) * s(mu).frobenius(3) # long time
            ....:           == (s(lam) * s(mu)).frobenius(3)
            ....:           for mu in Partitions(3) )
            ....:      for lam in Partitions(3) )
            True
            sage: all( all( m(lam).frobenius(2) * m(mu).frobenius(2)
            ....:           == (m(lam) * m(mu)).frobenius(2)
            ....:           for mu in Partitions(4) )
            ....:      for lam in Partitions(4) )
            True
            sage: all( all( p(lam).frobenius(2) * p(mu).frobenius(2)
            ....:           == (p(lam) * p(mu)).frobenius(2)
            ....:           for mu in Partitions(3) )
            ....:      for lam in Partitions(4) )
            True

        Being Hopf algebra endomorphisms, the Frobenius operators
        commute with the antipode::

            sage: all( p(lam).frobenius(4).antipode()
            ....:      == p(lam).antipode().frobenius(4)
            ....:      for lam in Partitions(3) )
            True

        Testing the `\mathbf{f}_n(P) = p_n \circ P = P \circ p_n`
        equality (over `\QQ`, since plethysm is currently not
        defined over `\ZZ` in Sage)::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: p = Sym.p()
            sage: all( s(lam).frobenius(3) == s(lam).plethysm(p[3])
            ....:      == s(p[3].plethysm(s(lam)))
            ....:      for lam in Partitions(4) )
            True

        By Exercise 7.61 in Stanley's EC2 [STA]_ (see the errata on his
        website), `\mathbf{f}_n(h_m)` is a linear combination of
        Schur polynomials (of straight shapes) using coefficients `0`,
        `1` and `-1` only; moreover, all partitions whose Schur
        polynomials occur with coefficient `\neq 0` in this
        combination have empty `n`-cores. Let us check this on
        examples::

            sage: all( all( all( (coeff == -1 or coeff == 1)
            ....:                and lam.core(n) == Partition([])
            ....:                for lam, coeff in s([m]).frobenius(n) )
            ....:           for n in range(2, 4) )
            ....:      for m in range(4) )
            True

        .. SEEALSO::

            :meth:`plethysm`

        .. TODO::

            This method is fast on the monomial and the powersum
            bases, while all other bases get converted to the
            monomial basis. For most bases, this is probably the
            quickest way to do, but at least the Schur basis should
            have a better option. (Quoting from Stanley's EC2 [STA]_:
            "D. G. Duncan, J. London Math. Soc. 27 (1952), 235-236,
            or Y. M. Chen, A. M. Garsia, and J. B. Remmel, Contemp.
            Math. 34 (1984), 109-153".)
        """
        # Convert to the monomial basis, there apply Frobenius componentwise,
        # then convert back.
        parent = self.parent()
        m = parent.realization_of().monomial()
        from sage.combinat.partition import Partition
        dct = {Partition([n * i for i in lam]): coeff
               for (lam, coeff) in m(self)}
        result_in_m_basis = m._from_dict(dct)
        return parent(result_in_m_basis)

    def adams_operation(self, *args, **opts):
        from sage.misc.superseded import deprecation
        deprecation(19255, "Do not use this method! Please use `frobenius` or `adams_operator` methods following what you expect.")

    def verschiebung(self, n):
        r"""
        Return the image of the symmetric function ``self`` under the
        `n`-th Verschiebung operator.

        The `n`-th Verschiebung operator `\mathbf{V}_n` is defined to be
        the unique algebra endomorphism `V` of the ring of symmetric
        functions that satisfies `V(h_r) = h_{r/n}` for every positive
        integer `r` divisible by `n`, and satisfies `V(h_r) = 0` for
        every positive integer `r` not divisible by `n`. This operator
        `\mathbf{V}_n` is a Hopf algebra endomorphism. For every
        nonnegative integer `r` with `n \mid r`, it satisfies

        .. MATH::

            \mathbf{V}_n(h_r) = h_{r/n},
            \quad \mathbf{V}_n(p_r) = n p_{r/n},
            \quad \mathbf{V}_n(e_r) = (-1)^{r - r/n} e_{r/n}

        (where `h` is the complete homogeneous basis, `p` is the
        powersum basis, and `e` is the elementary basis). For every
        nonnegative integer `r` with `n \nmid r`, it satisfes

        .. MATH::

            \mathbf{V}_n(h_r) = \mathbf{V}_n(p_r) = \mathbf{V}_n(e_r) = 0.

        The `n`-th Verschiebung operator is also called the `n`-th
        Verschiebung endomorphism. Its name derives from the Verschiebung
        (German for "shift") endomorphism of the Witt vectors.

        The `n`-th Verschiebung operator is adjoint to the `n`-th
        Frobenius operator (see :meth:`frobenius` for its definition)
        with respect to the Hall scalar product (:meth:`scalar`).

        The action of the `n`-th Verschiebung operator on the Schur basis
        can also be computed explicitly. The following (probably clumsier
        than necessary) description can be obtained by solving exercise
        7.61 in Stanley's [STA]_.

        Let `\lambda` be a partition. Let `n` be a positive integer. If
        the `n`-core of `\lambda` is nonempty, then
        `\mathbf{V}_n(s_\lambda) = 0`. Otherwise, the following method
        computes `\mathbf{V}_n(s_\lambda)`: Write the partition `\lambda`
        in the form `(\lambda_1, \lambda_2, \ldots, \lambda_{ns})` for some
        nonnegative integer `s`. (If `n` does not divide the length of
        `\lambda`, then this is achieved by adding trailing zeroes to
        `\lambda`.) Set `\beta_i = \lambda_i + ns - i` for every
        `s \in \{ 1, 2, \ldots, ns \}`. Then,
        `(\beta_1, \beta_2, \ldots, \beta_{ns})` is a strictly decreasing
        sequence of nonnegative integers. Stably sort the list
        `(1, 2, \ldots, ns)` in order of (weakly) increasing remainder of
        `-1 - \beta_i` modulo `n`. Let `\xi` be the sign of the
        permutation that is used for this sorting. Let `\psi` be the sign
        of the permutation that is used to stably sort the list
        `(1, 2, \ldots, ns)` in order of (weakly) increasing remainder of
        `i - 1` modulo `n`. (Notice that `\psi = (-1)^{n(n-1)s(s-1)/4}`.)
        Then, `\mathbf{V}_n(s_\lambda) = \xi \psi \prod_{i = 0}^{n - 1}
        s_{\lambda^{(i)}}`, where
        `(\lambda^{(0)}, \lambda^{(1)}, \ldots, \lambda^{(n - 1)})`
        is the `n`-quotient of `\lambda`.

        INPUT:

        - ``n`` -- a positive integer

        OUTPUT:

        The result of applying the `n`-th Verschiebung operator (on the ring of
        symmetric functions) to ``self``.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(ZZ)
            sage: p = Sym.p()
            sage: h = Sym.h()
            sage: s = Sym.s()
            sage: m = Sym.m()
            sage: s[3].verschiebung(2)
            0
            sage: s[3].verschiebung(3)
            s[1]
            sage: p[3].verschiebung(3)
            3*p[1]
            sage: m[3,2,1].verschiebung(3)
            -18*m[1, 1] - 3*m[2]
            sage: p[3,2,1].verschiebung(3)
            0
            sage: h[4].verschiebung(2)
            h[2]
            sage: p[2].verschiebung(2)
            2*p[1]
            sage: m[3,2,1].verschiebung(6)
            12*m[1]

        The Verschiebung endomorphisms are multiplicative::

            sage: all( all( s(lam).verschiebung(2) * s(mu).verschiebung(2)
            ....:           == (s(lam) * s(mu)).verschiebung(2)
            ....:           for mu in Partitions(4) )
            ....:      for lam in Partitions(4) )
            True

        Being Hopf algebra endomorphisms, the Verschiebung operators
        commute with the antipode::

            sage: all( p(lam).verschiebung(3).antipode()
            ....:      == p(lam).antipode().verschiebung(3)
            ....:      for lam in Partitions(6) )
            True

        Testing the adjointness between the Frobenius operators
        `\mathbf{f}_n` and the Verschiebung operators
        `\mathbf{V}_n`::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.s()
            sage: p = Sym.p()
            sage: all( all( s(lam).verschiebung(2).scalar(p(mu))
            ....:           == s(lam).scalar(p(mu).frobenius(2))
            ....:           for mu in Partitions(3) )
            ....:      for lam in Partitions(6) )
            True
        """
        # Convert to the complete homogenenous basis, there apply
        # Verschiebung componentwise, then convert back.
        parent = self.parent()
        h = parent.realization_of().homogeneous()
        from sage.combinat.partition import Partition
        dct = {Partition([i // n for i in lam]): coeff
               for (lam, coeff) in h(self)
               if all( i % n == 0 for i in lam )}
        result_in_h_basis = h._from_dict(dct)
        return parent(result_in_h_basis)

    def bernstein_creation_operator(self, n):
        r"""
        Return the image of ``self`` under the `n`-th Bernstein creation
        operator.

        Let `n` be an integer. The `n`-th Bernstein creation operator
        `\mathbf{B}_n` is defined as the endomorphism of the space
        `Sym` of symmetric functions which sends every `f` to

        .. MATH::

            \sum_{i \geq 0} (-1)^i h_{n+i} e_i^\perp,

        where usual notations are in place (`h` stands for the complete
        homogeneous symmetric functions, `e` for the elementary ones,
        and `e_i^\perp` means skewing (:meth:`skew_by`) by `e_i`).

        This has been studied in [BBSSZ2012]_, section 2.2, where the
        following rule is given for computing `\mathbf{B}_n` on a
        Schur function: If `(\alpha_1, \alpha_2, \ldots, \alpha_n)` is
        an `n`-tuple of integers (positive or not), then

        .. MATH::

            \mathbf{B}_n s_{(\alpha_1, \alpha_2, \ldots, \alpha_n)}
            = s_{(n, \alpha_1, \alpha_2, \ldots, \alpha_n)}.

        Here, `s_{(\alpha_1, \alpha_2, \ldots, \alpha_n)}` is the
        "Schur function" associated to the `n`-tuple
        `(\alpha_1, \alpha_2, \ldots, \alpha_n)`, and defined by
        literally applying the Jacobi-Trudi identity, i.e., by

        .. MATH::

            s_{(\alpha_1, \alpha_2, \ldots, \alpha_n)}
            = \det \left( (h_{\alpha_i - i + j})_{i, j = 1, 2, \ldots, n} \right).

        This notion of a Schur function clearly extends the classical
        notion of Schur function corresponding to a partition, but is
        easily reduced to the latter (in fact, for any `n`-tuple
        `\alpha` of integers, one easily sees that `s_\alpha` is
        either `0` or minus-plus a Schur function corresponding to a
        partition; and it is easy to determine which of these is the
        case and find the partition by a combinatorial algorithm).

        EXAMPLES:

        Let us check that what this method computes agrees with the
        definition::

            sage: Sym = SymmetricFunctions(ZZ)
            sage: e = Sym.e()
            sage: h = Sym.h()
            sage: s = Sym.s()
            sage: def bernstein_creation_by_def(n, f):
            ....:     # `n`-th Bernstein creation operator applied to `f`
            ....:     # computed according to its definition.
            ....:     res = f.parent().zero()
            ....:     if not f:
            ....:         return res
            ....:     max_degree = max(sum(m) for m, c in f)
            ....:     for i in range(max_degree + 1):
            ....:         if n + i >= 0:
            ....:             res += (-1) ** i * h[n + i] * f.skew_by(e[i])
            ....:     return res
            sage: all( bernstein_creation_by_def(n, s[l]) == s[l].bernstein_creation_operator(n)
            ....:      for n in range(-2, 3) for l in Partitions(4) )
            True
            sage: all( bernstein_creation_by_def(n, s[l]) == s[l].bernstein_creation_operator(n)
            ....:      for n in range(-3, 4) for l in Partitions(3) )
            True
            sage: all( bernstein_creation_by_def(n, e[l]) == e[l].bernstein_creation_operator(n)
            ....:      for n in range(-3, 4) for k in range(3) for l in Partitions(k) )
            True

        Some examples::

            sage: s[3,2].bernstein_creation_operator(3)
            s[3, 3, 2]
            sage: s[3,2].bernstein_creation_operator(1)
            -s[2, 2, 2]
            sage: h[3,2].bernstein_creation_operator(-2)
            h[2, 1]
            sage: h[3,2].bernstein_creation_operator(-1)
            h[2, 1, 1] - h[2, 2] - h[3, 1]
            sage: h[3,2].bernstein_creation_operator(0)
            -h[3, 1, 1] + h[3, 2]
            sage: h[3,2].bernstein_creation_operator(1)
            -h[2, 2, 2] + h[3, 2, 1]
            sage: h[3,2].bernstein_creation_operator(2)
            -h[3, 3, 1] + h[4, 2, 1]
        """
        # We use the formula for the Bernstein creation operator on
        # a Schur function given in the docstring.
        from sage.combinat.partition import _Partitions
        parent = self.parent()
        s = parent.realization_of().schur()
        res = s.zero()
        for m, c in s(self): # m = monomial (= corresponding partition), c = coefficient
            # Add ``c * s[m].bernstein_creation_operator()`` to ``res``.
            # There is a simple combinatorial algorithm for this (using
            # the Jacobi-Trudi formula), which returns either 0 or
            # minus-plus a single Schur function.
            for j, p in enumerate(m + [0]):
                # The "+ [0]" is important and corresponds to moving the ``n``
                # to the very end!
                if n == p - j - 1:
                    break
                if n > p - j - 1:
                    if n + j < 0:
                        break
                    m_new = [k - 1 for k in m[:j]] + [n + j] + m[j:]
                    m_new = _Partitions(m_new)
                    res += (-1) ** j * c * s[m_new]
                    break
        return parent(res)

    def _expand(self, condition, n, alphabet = 'x'):
        r"""
        Expand the symmetric function as a symmetric polynomial in ``n``
        variables.

        INPUT:

        - ``condition`` -- a function on partitions with a boolean output,
          selecting only certain terms (namely, only the items failing
          the condition are being expanded)

        - ``n`` -- a nonnegative integer

        - ``alphabet`` -- (default: ``'x'``) a variable for the expansion

        OUTPUT:

        A monomial expansion of ``self`` in the `n` variables labelled
        by ``alphabet``.

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

        .. NOTE::

            The term corresponding to the empty partition is always
            selected, even if ``condition`` returns ``False`` or an
            error when applied to the empty partition. This is in
            order to simplify using the ``_expand`` method with
            conditions like ``lambda part: max(part) < 3`` which
            would require extra work to handle the empty partition.
        """
        import classical
        parent = self.parent()
        resPR = PolynomialRing(parent.base_ring(), n, alphabet)
        if self == parent.zero():
            return resPR.zero()
        e = eval('symmetrica.compute_' + str(classical.translate[parent.basis_name()]).lower() + '_with_alphabet')
        def f(part):
            if not part:
                return resPR.one()
            else:
                return resPR.zero() if condition(part) else resPR(e(part, n, alphabet))
        return parent._apply_module_morphism(self, f)

    def is_schur_positive(self):
        r"""
        Return ``True`` if and only if ``self`` is Schur positive.

        If `s` is the space of Schur functions over ``self``'s base ring, then
        this is the same as ``self._is_positive(s)``.

        EXAMPLES::

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
        Return ``True`` if and only if ``self`` has nonnegative coefficients
        in the basis `s`.

        INPUT:

        - ``s`` -- a basis of the ring of symmetric functions

        EXAMPLES::

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
        return all(( _nonnegative_coefficients(c) for c in s_self.coefficients() ))

    def degree(self):
        r"""
        Return the degree of ``self`` (which is defined to be `0`
        for the zero element).

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1]) + 3
            sage: z.degree()
            4
            sage: s(1).degree()
            0
            sage: s(0).degree()
            0
        """
        return max( [sum(_) for _ in self._monomial_coefficients] + [0] )

    def restrict_degree(self, d, exact = True):
        r"""
        Return the degree ``d`` component of ``self``.

        INPUT:

        - ``d`` -- positive integer, degree of the terms to be returned

        - ``exact`` -- boolean, if ``True``, returns the terms of degree
          exactly ``d``, otherwise returns all terms of degree less than
          or equal to ``d``

        OUTPUT:

        - the homogeneous component of ``self`` of degree ``d``

        EXAMPLES::

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
            res = dict(x for x in self._monomial_coefficients.items() if sum(x[0]) == d)
        else:
            res = dict(x for x in self._monomial_coefficients.items() if sum(x[0]) <= d)
        return self.parent()._from_dict(res)

    def restrict_partition_lengths(self, l, exact = True):
        r"""
        Return the terms of ``self`` labelled by partitions of length ``l``.

        INPUT:

        - ``l`` -- nonnegative integer

        - ``exact`` -- boolean, defaulting to ``True``

        OUTPUT:

        - if ``True``, returns the terms labelled by
          partitions of length precisely ``l``; otherwise returns all terms
          labelled by partitions of length less than or equal to ``l``

        EXAMPLES::

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
            res = dict(x for x in self._monomial_coefficients.items() if len(x[0]) == l)
        else:
            res = dict(x for x in self._monomial_coefficients.items() if len(x[0]) <= l)
        return self.parent()._from_dict(res)

    def restrict_parts(self, n):
        r"""
        Return the terms of ``self`` labelled by partitions `\lambda` with
        `\lambda_1 \leq n`.

        INPUT:

        - ``n`` -- positive integer, to restrict the parts of the partitions
          of the terms to be returned

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.restrict_parts(2)
            s[1] + s[1, 1, 1] + s[2, 1]
            sage: z.restrict_parts(1)
            s[1] + s[1, 1, 1]
        """
        res = dict(x for x in self._monomial_coefficients.items() if _lmax(x[0]) <= n)
        return self.parent()._from_dict(res)

    def expand(self, n, alphabet = 'x'):
        r"""
        Expand the symmetric function ``self`` as a symmetric polynomial
        in ``n`` variables.

        INPUT:

        - ``n`` -- a nonnegative integer

        - ``alphabet`` -- (default: ``'x'``) a variable for the expansion

        OUTPUT:

        A monomial expansion of ``self`` in the `n` variables
        labelled ``x0``, ``x1``, ..., ``x{n-1}`` (or just ``x``
        if `n = 1`), where ``x`` is ``alphabet``.

        EXAMPLES::

            sage: J = SymmetricFunctions(QQ).jack(t=2).J()
            sage: J([2,1]).expand(3)
            4*x0^2*x1 + 4*x0*x1^2 + 4*x0^2*x2 + 6*x0*x1*x2 + 4*x1^2*x2 + 4*x0*x2^2 + 4*x1*x2^2
            sage: (2*J([2])).expand(0)
            0
            sage: (3*J([])).expand(0)
            3
        """
        s = self.parent().realization_of().schur()
        condition = lambda part: len(part) > n
        return s(self)._expand(condition, n, alphabet)

    def skew_by(self, x):
        r"""
        Return the result of skewing ``self`` by ``x``. (Skewing by ``x`` is
        the endomorphism (as additive group) of the ring of symmetric
        functions adjoint to multiplication by ``x`` with respect to the
        Hall inner product.)

        INPUT:

        - ``x`` -- element of the ring of symmetric functions over the same
          base ring as ``self``

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
        parent = self.parent()
        Sym = parent.realization_of()
        if x not in Sym:
            raise ValueError("x needs to be a symmetric function")
        s = Sym.schur()
        zero = s.zero()
        f = lambda part1, part2: s([part1,part2]) if part1.contains(part2) else zero
        return parent(s._apply_multi_module_morphism(s(self), s(x), f))

    def hl_creation_operator(self, nu, t = None):
        r"""
        This is the vertex operator that generalizes Jing's operator.

        It is a linear operator that raises the degree by
        `|\nu|`. This creation operator is a t-analogue of
        multiplication by ``s(nu)`` .

        .. SEEALSO:: Proposition 5 in [SZ2001]_.

        INPUT:

        - ``nu`` -- a partition

        - ``t`` -- (default: ``None``, in which case ``t`` is used) a parameter

        REFERENCES:

        .. [SZ2001] \M. Shimozono, M. Zabrocki,
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

    def eval_at_permutation_roots(self, rho):
        r"""
        Evaluate at eigenvalues of a permutation matrix.

        Evaluate a symmetric function at the eigenvalues of a permutation
        matrix whose cycle structure is ``rho``.  This computation is
        computed by coercing to the power sum basis where the value may
        be computed on the generators.

        This function evaluates an element at the roots of unity

        .. MATH::

            \Xi_{\rho_1},\Xi_{\rho_2},\ldots,\Xi_{\rho_\ell}

        where

        .. MATH::

            \Xi_{m} = 1,\zeta_m,\zeta_m^2,\ldots,\zeta_m^{m-1}

        and `\zeta_m` is an `m` root of unity.
        These roots of unity represent the eigenvalues of permutation
        matrix with cycle structure `\rho`.

        INPUT:

        - ``rho`` -- a partition or a list of non-negative integers

        OUTPUT:

        - an element of the base ring

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s([3,3]).eval_at_permutation_roots([6])
            0
            sage: s([3,3]).eval_at_permutation_roots([3])
            1
            sage: s([3,3]).eval_at_permutation_roots([1])
            0
            sage: s([3,3]).eval_at_permutation_roots([3,3])
            4
            sage: s([3,3]).eval_at_permutation_roots([1,1,1,1,1])
            175
            sage: (s[1]+s[2]+s[3]).eval_at_permutation_roots([3,2])
            2
        """
        p = self.parent().symmetric_function_ring().p()
        return p(self).eval_at_permutation_roots(rho)

    def character_to_frobenius_image(self, n):
        r"""
        Interpret ``self`` as a `Gl_n` character and then take the Frobenius
        image of this character of the permutation matrices `S_n` which
        naturally sit inside of `Gl_n`.

        To know the value of this character at a permutation of cycle structure
        `\rho` the symmetric function ``self`` is evaluated at the
        eigenvalues of of a permutation of cycle structure `\rho`.  The
        Frobenius image is then defined as
        `\sum_{\rho \vdash n} f[ \Xi_\rho ] p_\rho/z_\rho`.

        .. SEEALSO::
            :meth:`eval_at_permutation_roots`

        INPUT:

        - ``n`` -- a non-negative integer to interpret ``self`` as
          a character of `Gl_n`

        OUTPUT:

        - a symmetric function of degree ``n``

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: s([1,1]).character_to_frobenius_image(5)
            s[3, 1, 1] + s[4, 1]
            sage: s([2,1]).character_to_frobenius_image(5)
            s[2, 2, 1] + 2*s[3, 1, 1] + 2*s[3, 2] + 3*s[4, 1] + s[5]
            sage: s([2,2,2]).character_to_frobenius_image(3)
            s[3]
            sage: s([2,2,2]).character_to_frobenius_image(4)
            s[2, 2] + 2*s[3, 1] + 2*s[4]
            sage: s([2,2,2]).character_to_frobenius_image(5)
            2*s[2, 2, 1] + s[3, 1, 1] + 4*s[3, 2] + 3*s[4, 1] + 2*s[5]
        """
        p = self.parent().symmetric_function_ring().p()
        return self.parent()(p.sum(self.eval_at_permutation_roots(rho) \
            *p(rho)/rho.centralizer_size() for rho in Partitions(n)))

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
    return max(x) if x else 0

def _nonnegative_coefficients(x):
    r"""
    Returns ``True`` if ``x`` has nonnegative coefficients.

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
        return all([ c >= 0 for c in x.coefficients(sparse=False) ])
    else:
        return x >= 0

