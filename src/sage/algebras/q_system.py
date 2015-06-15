r"""
Q-Systems

AUTHORS:

- Travis Scrimshaw (2013-10-08): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.algebras import Algebras
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.rings.all import ZZ, QQ
from sage.rings.infinity import infinity
from sage.rings.arith import LCM
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.combinat.cartesian_product import CartesianProduct
from sage.combinat.free_module import CombinatorialFreeModule
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid, IndexedMonoid
from sage.matrix.constructor import matrix
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.cartan_matrix import CartanMatrix

class QSystem(Parent, UniqueRepresentation):
    r"""
    A Q-system.

    Let `\mathfrak{g}` be a symmetrizable Kac-Moody algebra with index
    set `I` over a field `k`. Follow the presentation given in [HKOTY99]_,
    an unrestricted Q-system is a `k`-algebra in infinitely many variables
    `Q_m^{(a)}`, where `m \in \ZZ_{>0}` and `a \in I`, which satisifies
    the relations

    .. MATH::

        (Q_m^{(a)})^2 = Q_{m+1}^{(a)} Q_{m-1}^{(a)} +
        \prod_{b \sim a} \prod_{k=0}^{-C_{ab} - 1}
        Q^{(b)}_{\lfloor \frac{m C_{ba} - k}{C_{ab}} \rfloor},

    with `Q_0^{(a)} := 1`. Q-systems can be considered as T-systems where
    we forget the spectral parameter `u` and for `\mathfrak{g}` of finite
    type, have a solution given by the characters of Kirillov-Reshetikhin
    modules (again without the spectral parameter) for an affine Kac-Moody
    algebra `\widehat{\mathfrak{g}}` with `\mathfrak{g}` as its classical
    subalgebra.

    Q-systems have two natural bases:

    - :class:`~QSystem.Fundamental` -- given by polynomials of the
      fundamental representations `Q_1^{(a)}`
    - :class:`~QSystem.SqaureFree` -- given by square-free terms

    There is also a level `\ell` restricted Q-system (with unit boundary
    condition) given by setting `Q_{d_a \ell}^{(a)}` = 1`, where `d_a`
    are the entries of the symmetrizing matrix for the dual type of
    `\mathfrak{g}`.

    EXAMPLES:

    We begin by constructing a Q-system and the two bases::

        sage: Q = QSystem(['A', 4])
        sage: F = Q.Fundamental()
        sage: SF = Q.SquareFree()

    REFERENCES:

    .. [HKOTY99] G. Hatayama, A. Kuniba, M. Okado, T. Tagaki, and Y. Yamada.
       *Remarks on the fermionic formula*. Contemp. Math., **248** (1999).
    """
    @staticmethod
    def __classcall__(cls, base_ring, cartan_type, level=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: Q1 = QSystem(QQ, ['A',4])
            sage: Q2 = QSystem(QQ, 'A4')
            sage: Q1 is Q2
            True
        """
        cartan_type = CartanType(cartan_type)
        # TODO: Check for tamely laced!!!
        return super(QSystem, cls).__classcall__(cls, base_ring, cartan_type, level)

    def __init__(self, base_ring, cartan_type, level):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: TestSuite(Q).run()
        """
        category = Algebras(base_ring).Commutative()
        self._cartan_type = cartan_type
        self._level = level
        Parent.__init__(self, base_ring, category=category.WithRealizations())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: QSystem(QQ, ['A',4])
            Q-system of type ['A', 4] over Rational Field
        """
        if self._level is not None:
            res = "Restricted level {} ".format(self._level)
        else:
            res = ''
        return "{}Q-system of type {} over {}".format(res, self._cartan_type, self.base_ring())

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: Q.cartan_type()
            ['A', 4]
        """
        return self._cartan_type

    def level(self):
        """
        Return the restriction level of ``self`` or ``None`` if
        the system is unrestricted.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: Q.level()

            sage: Q = QSystem(QQ, ['A',4], 5)
            sage: Q.level()
            5
        """
        return self._level

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the fundamental basis).

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: Q.a_realization()
            Q-system of type ['A', 4] over Rational Field
             in the Fundamental basis
        """
        return self.Fundamental()

    class _BasesCategory(Category_realization_of_parent):
        r"""
        The category of bases of a Iwahori-Hecke algebra.
        """
        def super_categories(self):
            r"""
            The super categories of ``self``.

            EXAMPLES::

                sage: Q = QSystem(QQ, ['A',4])
                sage: bases = Q._BasesCategory()
                sage: bases.super_categories()
                [Category of realizations of Q-system of type ['A', 4] over Rational Field,
                 Category of commutative algebras with basis over quotient fields]
            """
            cat = Algebras(self.base().base_ring().category()).Commutative().WithBasis()
            return [Realizations(self.base()), cat]

        def _repr_(self):
            r"""
            Return the representation of ``self``.

            EXAMPLES::

                sage: Q = QSystem(QQ, ['A',4])
                sage: Q._BasesCategory()
                Category of bases of Q-system of type ['A', 4] over Rational Field
            """
            return "Category of bases of {}".format(self.base())

        class ParentMethods:
            r"""
            This class collects code common to all the various bases.
            """
            def _repr_(self):
                """
                Text representation of this basis of Iwahori-Hecke algebra.

                EXAMPLES::

                    sage: Q = QSystem(QQ, ['A',4])
                    sage: Q.Fundamental()
                    Q-system of type ['A', 4] over Rational Field
                     in the Fundamental basis
                    sage: Q.SquareFree()
                    Q-system of type ['A', 4] over Rational Field
                     in the Square-Free basis
                """
                return "{} in the {} basis".format(self.realization_of(), self._basis_name)

            def level(self):
                """
                Return the restriction level of ``self`` or ``None`` if
                the system is unrestricted.

                EXAMPLES::

                    sage: F = QSystem(QQ, ['A',4]).Fundamental()
                    sage: F.level()

                    sage: F = QSystem(QQ, ['A',4], 5).Fundamental()
                    sage: F.level()
                    5
                """
                return self._level

            def cartan_type(self):
                """
                Return the Cartan type of ``self``.

                EXAMPLES::

                    sage: F = QSystem(QQ, ['A',4]).Fundamental()
                    sage: F.cartan_type()
                    ['A', 4]
                """
                return self._cartan_type

            def index_set(self):
                """
                Return the index set of ``self``.

                EXAMPLES::

                    sage: F = QSystem(QQ, ['A',4]).Fundamental()
                    sage: F.index_set()
                    (1, 2, 3, 4)
                """
                return self._cartan_type.index_set()

            @cached_method
            def one_basis(self):
                """
                Return the basis element indexing `1`.

                EXAMPLES::

                    sage: F = QSystem(QQ, ['A',4]).Fundamental()
                    sage: F.one_basis()
                    1
                    sage: F.one_basis().parent() is F._indices
                    True
                """
                return self._indices.one()

            @cached_method
            def algebra_generators(self):
                """
                Return the algebra generators of ``self``.

                EXAMPLES::

                    sage: F = QSystem(QQ, ['A',4]).Fundamental()
                    sage: F.algebra_generators()
                    Finite family {1: Q^(1)[1], 2: Q^(2)[1], 3: Q^(3)[1], 4: Q^(4)[1]}
                """
                return Family({a: self.gen(a, 1) for a in self.index_set()})

            def dimension(self):
                """
                Return the dimension of ``self``, which is `\infty`.

                EXAMPLES::

                    sage: F = QSystem(QQ, ['A',4]).Fundamental()
                    sage: F.dimension()
                    +Infinity
                """
                return infinity

    class _BasisAbstract(CombinatorialFreeModule, BindableClass):
        """
        Abstract base class for a basis of a Q-system.
        """
        def __init__(self, Q_system, indices):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: F = QSystem(QQ, ['A',2]).Fundamental()
                sage: TestSuite(F).run()
            """
            self._cartan_type = Q_system._cartan_type
            self._level = Q_system._level
            category = Q_system._BasesCategory()
            CombinatorialFreeModule.__init__(self, Q_system.base_ring(), indices,
                                             prefix='Q', category=category)

        def _repr_term(self, t):
            """
            Return a string representation of the basis element indexed by ``a``
            with all `m = 1`.

            EXAMPLES::

                sage: F = QSystem(QQ, ['A',4]).Fundamental()
                sage: I = F._indices
                sage: F._repr_term( I.gen((1,1)) * I.gen((4,1)) )
                'Q^(1)[1]*Q^(4)[1]'
            """
            if len(t) == 0:
                return '1'
            def repr_gen(x):
                ret = 'Q^({})[{}]'.format(*(x[0]))
                if x[1] > 1:
                    ret += '^{}'.format(x[1])
                return ret
            return '*'.join(repr_gen(x) for x in t._sorted_items())

        def _latex_term(self, t):
            r"""
            Return a `\LaTeX` representation of the basis element indexed
            by ``m``.

            EXAMPLES::

                sage: F = QSystem(QQ, ['A',4]).Fundamental()
                sage: I = F._indices
                sage: F._latex_term( I.gen((3,1)) * I.gen((4,1)) )
                'Q^{(3)}_{1} Q^{(4)}_{1}'
            """
            if len(t) == 0:
                return '1'
            def repr_gen(x):
                ret = 'Q^{{({})}}_{{{}}}'.format(*(x[0]))
                if x[1] > 1:
                    ret = '\\bigl(' + ret + '\\bigr)^{{{}}}'.format(x[1])
                return ret
            return ' '.join(repr_gen(x) for x in t._sorted_items())

        # TODO: Once T-systems are implemented
        #def _coerce_map_from_(self, M):
        #    """
        #    Return a morphism if there is a coercion map from ``M``
        #    into ``self``.
        #    """
        #    if isinstance(M, TSystem):
        #        if (M.cartan_type() == self.cartan_type()
        #            and self.base_ring().has_coerce_map_from(M.base_ring())):
        #            return M.module_morphism(M.restriction_on_basis, codomain=self)
        #    return super(QSystem._BasisAbstract, self)._coerce_map_from_(M)

    class Fundamental(_BasisAbstract):
        """
        The basis of a Q-system given by polynomials of the
        fundamental representations.
        """
        _basis_name = "Fundamental"

        def __init__(self, Q_system):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: F = QSystem(QQ, ['A',4]).Fundamental()
                sage: TestSuite(F).run()
            """
            cartan_type = Q_system._cartan_type
            indices = CartesianProduct(cartan_type.index_set(), [1])
            basis = IndexedFreeAbelianMonoid(indices, prefix='Q', bracket=False)
            # This is used to do the reductions
            self._poly = PolynomialRing(ZZ, ['q'+str(i) for i in cartan_type.index_set()])
            QSystem._BasisAbstract.__init__(self, Q_system, basis)

        def gen(self, a, m):
            """
            Return the generator `Q^{(a)}_i` of ``self``.

            EXAMPLES::

                sage: F = QSystem(QQ, ['A',8]).Fundamental()
                sage: F.gen(2, 1)
                Q^(2)[1]
                sage: F.gen(6, 2)
                Q^(6)[1]^2 - Q^(5)[1]*Q^(7)[1]
                sage: F.gen(7, 3)
                -2*Q^(6)[1]*Q^(7)[1]*Q^(8)[1] + Q^(6)[1]^2
                 + Q^(5)[1]*Q^(8)[1]^2 - Q^(5)[1]*Q^(7)[1] + Q^(7)[1]^3
                sage: F.gen(1, 0)
                1
            """
            if a not in self.index_set():
                raise ValueError("a is not in the index set")
            if m == 0:
                return self.one()
            if self._level:
                t = self._cartan_type.dual().cartan_matrix().symmetrizer()
                if m == t[a] * self._level:
                    return self.one()
            if m == 1:
                return self.monomial( self._indices.gen((a,1)) )
            #if self._cartan_type.type() == 'A' and self._level is None:
            #    return self._jacobi_trudy(a, m)
            I = self._cartan_type.index_set()
            p = self._Q_poly(a, m)
            return p.subs({ g: self.gen(I[i], 1) for i,g in enumerate(self._poly.gens()) })

        @cached_method
        def _Q_poly(self, a, m):
            r"""
            Return the element `Q^{(a)}_m` as a polynomial.

            We start with the relation

            .. MATH::

                Q^{(a)}_{m-1}^2 = Q^{(a)}_m Q^{(a)}_{m-2} + \mathcal{Q}_{a,m-1},

            which implies

            .. MATH::

                Q^{(a)}_m = \frac{Q^{(a)}_{m-1}^2 - \mathcal{Q}_{a,m-1}}{
                Q^{(a)}_{m-2}}.

            .. NOTE::

                This helper method is defined in order to use the
                division implemented in polynomial rings.

            EXAMPLES::

                sage: F = QSystem(QQ, ['A',8]).Fundamental()
                sage: F._Q_poly(1, 2)
                q1^2 - q2
                sage: F._Q_poly(3, 2)
                q3^2 - q2*q4
                sage: F._Q_poly(6, 3)
                q6^3 - 2*q5*q6*q7 + q4*q7^2 + q5^2*q8 - q4*q6*q8
            """
            if m == 0 or m == self._level:
                return self._poly.one()
            if m == 1:
                return self._poly.gen(self._cartan_type.index_set().index(a))

            cm = CartanMatrix(self._cartan_type)
            I = self._cartan_type.index_set()
            m -= 1 # So we don't have to do it everywhere

            cur = self._Q_poly(a, m)**2
            i = I.index(a)
            ret = self._poly.one()
            for b in self._cartan_type.dynkin_diagram().neighbors(a):
                j = I.index(b)
                for k in range(-cm[i,j]):
                    ret *= self._Q_poly(b, (m * cm[j,i] - k) // cm[i,j])
            cur -= ret
            if m > 1:
                cur //= self._Q_poly(a, m-1)
            return cur

        class Element(CombinatorialFreeModule.Element):
            def _mul_(self, x):
                """
                Return the product of ``self`` and ``x``.

                EXAMPLES::

                    sage: F = QSystem(QQ, ['A',8]).Fundamental()
                    sage: x = F.gen(1, 2)
                    sage: y = F.gen(3, 2)
                    sage: x * y
                    -Q^(1)[1]^2*Q^(2)[1]*Q^(4)[1] + Q^(2)[1]^2*Q^(4)[1]
                     - Q^(2)[1]*Q^(3)[1]^2 + Q^(1)[1]^2*Q^(3)[1]^2
                """
                return self.parent().sum_of_terms((tl*tr, cl*cr)
                                                  for tl,cl in self for tr,cr in x)

    class SquareFree(_BasisAbstract):
        """
        The Q-system basis given in square-free monomials.
        """
        _basis_name = "Square-Free"

        def __init__(self, Q_system):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: SF = QSystem(QQ, ['A',2]).SquareFree()
                sage: TestSuite(SF).run()
            """
            basis = SquareFreeMonomials(Q_system._cartan_type)
            QSystem._BasisAbstract.__init__(self, Q_system, basis)

        def gen(self, a, m):
            """
            Return the generator `Q^{(a)}_i` of ``self``.

            EXAMPLES::

                sage: SF = QSystem(QQ, ['A',20]).SquareFree()
                sage: SF.gen(2, 1)
                Q^(2)[1]
                sage: SF.gen(12, 2)
                Q^(12)[2]
                sage: SF.gen(1, 0)
                1
            """
            if a not in self.index_set():
                raise ValueError("a is not in the index set")
            if m == 0 or m == self._level:
                return self.one()
            return self.monomial(self._indices.gen((a,m)))

        @lazy_attribute
        def _diagonal(self):
            """
            Return the diagonal of the symmetrized Cartan matrix of ``self``.

            EXAMPLES::

                sage: SF = QSystem(QQ, ['B',4]).SquareFree()
                sage: SF._diagonal
                Finite family {1: 2, 2: 2, 3: 2, 4: 1}
            """
            d = self._cartan_type.cartan_matrix().is_symmetrizable(True)
            lcd = LCM([QQ(x).denominator() for x in d])
            return Family({a: ZZ(d[i]*lcd) for i,a in enumerate(self.index_set())})

        @lazy_attribute
        def _lcm_diagonal(self):
            """
            Return the least common multiple of the diagonal.

            EXAMPLES::

                sage: SF = QSystem(QQ, ['B',4]).SquareFree()
                sage: SF._lcm_diagonal
                2
            """
            return LCM(list(self._diagonal))

        def _reduce_square_free(self, a, m):
            r"""
            Return the element `(Q^{(a)}_m)^2` reduced to a square-free
            polynomial.

            Return the square reduced by the formula

            .. MATH::

                \left( Q^{(a)}_m \right)^2 = Q^{(a)}_{m-1} Q^{(a)}_{m+1}
                + \prod_{a \sim b} Q(b)

            where if `d_a > 1` then `Q(b) = Q_{\beta}^{(b)}` where
            `\beta = d_a \left\lfloor \frac{m-1}{d_b} \right\rfloor`, otherwise

            .. MATH::

                Q(b) = \prod_{k=1}^{d_b} Q^{(b)}_{\beta_k}

            where `\beta_k = 1 + \left\lfloor \frac{m-1-k}{d_b} \right\rfloor`.


            EXAMPLES::

                sage: SF = QSystem(QQ, ['A', 6]).SquareFree()
                sage: SF._reduce_square_free(3, 7)
                Q^(2)[7]*Q^(4)[7] + Q^(3)[6]*Q^(3)[8]

                sage: SF = QSystem(QQ, ['B', 4]).SquareFree()
                sage: SF._reduce_square_free(3, 4)
                Q^(2)[4]*Q^(4)[8] + Q^(3)[3]*Q^(3)[5]
            """
            d = self._diagonal
            da = d[a]

            cur = self.gen(a,m-1) * self.gen(a,m+1)
            N = self._cartan_type.dynkin_diagram().neighbors(a)
            if da > 1:
                cur += self.prod(self.gen(b, (da * m) // d[b])
                                 for b in N)
            else:
                cur += self.prod(self.gen(b, 1 + (m-k) // d[b])
                                 for b in N for k in range(1, d[b]+1))
            return cur

        class Element(CombinatorialFreeModule.Element):
            def _mul_(self, x):
                """
                Return the product of ``self`` and ``x``.

                EXAMPLES::

                    sage: SF = QSystem(QQ, ['A', 4]).SquareFree()
                    sage: x = SF.an_element() + SF.gen(4, 2) + SF.gen(3,5); x
                    Q^(3)[5] + Q^(4)[2] + Q^(1)[1]*Q^(1)[2]*Q^(1)[3]
                    sage: y = SF.gen(1, 4) + SF.gen(2,2)
                    sage: x * y
                    Q^(2)[2]*Q^(4)[2] + Q^(1)[4]*Q^(4)[2] + Q^(1)[4]*Q^(3)[5]
                     + Q^(1)[1]*Q^(1)[2]*Q^(1)[3]*Q^(2)[2]
                     + Q^(1)[1]*Q^(1)[2]*Q^(1)[3]*Q^(1)[4] + Q^(2)[2]*Q^(3)[5]
                """
                P = self.parent()
                cur = P.sum_of_terms((tl*tr, cl*cr) for tl,cl in self for tr,cr in x)
                is_square_free = False
                while not is_square_free:
                    is_square_free = True
                    for m,c in cur:
                        ret = []
                        for a,exp in m._sorted_items():
                            if exp > 1:
                                assert exp == 2, "larger than square"
                                ret.append(P._reduce_square_free(*a))
                                is_square_free = False
                            else:
                                ret.append(P.gen(*a))
                        if not is_square_free:
                            cur += c * P.prod(ret) - P._from_dict({m:c}, remove_zeros=False)
                            break
                return cur


###########################################################
## Helper classes

class SquareFreeMonomials(IndexedFreeAbelianMonoid):
    """
    Parent representing the set of square-free monomials.
    """
    @staticmethod
    def __classcall__(cls, ct):
        """
        Normalize input.

        EXAMPLES::

            sage: from sage.algebras.q_system import SquareFreeMonomials
            sage: S1 = SquareFreeMonomials(['A',3])
            sage: S2 = SquareFreeMonomials('A3')
            sage: S1 is S2
            True
        """
        return super(IndexedMonoid, cls).__classcall__(cls, CartanType(ct))

    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.q_system import SquareFreeMonomials
            sage: SFM = SquareFreeMonomials(['A',2])
            sage: TestSuite(SFM).run()
        """
        self._cartan_type = cartan_type
        indices = CartesianProduct(cartan_type.index_set(), PositiveIntegers())
        IndexedFreeAbelianMonoid.__init__(self, indices, prefix='Q')

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.q_system import SquareFreeMonomials
            sage: SquareFreeMonomials(['A',2])
            Square-free monomials of the Q-system of type ['A', 2]
        """
        return "Square-free monomials of the Q-system of type {}".format(self._cartan_type)

    def _an_element_(self):
        """
        Return an element in ``self``.

        EXAMPLES::

            sage: from sage.algebras.q_system import SquareFreeMonomials
            sage: SFM = SquareFreeMonomials(['A',2])
            sage: SFM._an_element_()
            Q[[1, 1]]*Q[[1, 2]]*Q[[1, 3]]
        """
        x = self.one()
        I = self._indices
        try:
            g = iter(self._indices)
            for c in range(1,4):
                x *= self.gen(next(g))
        except Exception:
            pass
        return x

