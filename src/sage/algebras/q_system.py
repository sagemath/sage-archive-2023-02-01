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

import itertools
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

from sage.categories.algebras import Algebras
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.rings.all import ZZ, QQ
from sage.rings.infinity import infinity
from sage.rings.arith import LCM
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.combinat.free_module import CombinatorialFreeModule
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid, IndexedMonoid
from sage.matrix.constructor import matrix
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.cartan_matrix import CartanMatrix

class QSystem(CombinatorialFreeModule):
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

    Q-systems have a natural bases given by polynomials of the
    fundamental representations `Q_1^{(a)}`

    There is also a level `\ell` restricted Q-system (with unit boundary
    condition) given by setting `Q_{d_a \ell}^{(a)}` = 1`, where `d_a`
    are the entries of the symmetrizing matrix for the dual type of
    `\mathfrak{g}`.

    EXAMPLES:

    We begin by constructing a Q-system and doing some basic computations
    in type `A_4`::

        sage: Q = QSystem(QQ, ['A', 4])
        sage: Q.gen(3,1)
        Q^(3)[1]
        sage: Q.gen(1,2)
        Q^(1)[1]^2 - Q^(2)[1]
        sage: Q.gen(3,3)
        -Q^(1)[1]*Q^(3)[1] + Q^(1)[1]*Q^(4)[1]^2 + Q^(2)[1]^2
         - 2*Q^(2)[1]*Q^(3)[1]*Q^(4)[1] + Q^(3)[1]^3
        sage: x = Q.gen(1,1) + Q.gen(2,1); x
        Q^(1)[1] + Q^(2)[1]
        sage: x * x
        Q^(1)[1]^2 + 2*Q^(1)[1]*Q^(2)[1] + Q^(2)[1]^2

    Next we do some basic computations in type `C_4`::

        sage: Q = QSystem(QQ, ['C', 4])
        sage: Q.gen(4,1)
        Q^(4)[1]
        sage: Q.gen(1,2)
        Q^(1)[1]^2 - Q^(2)[1]
        sage: Q.gen(3,3)
        Q^(1)[1]*Q^(4)[1]^2 - 2*Q^(2)[1]*Q^(3)[1]*Q^(4)[1] + Q^(3)[1]^3

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

            sage: Q = QSystem(QQ, ['A',2])
            sage: TestSuite(Q).run()
        """
        self._cartan_type = cartan_type
        self._level = level
        indices = tuple(itertools.product(cartan_type.index_set(), [1]))
        basis = IndexedFreeAbelianMonoid(indices, prefix='Q', bracket=False)
        # This is used to do the reductions
        self._poly = PolynomialRing(ZZ, ['q'+str(i) for i in cartan_type.index_set()])

        category = Algebras(base_ring).Commutative().WithBasis()
        CombinatorialFreeModule.__init__(self, base_ring, basis,
                                         prefix='Q', category=category)

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

    def _repr_term(self, t):
        """
        Return a string representation of the basis element indexed by ``a``
        with all `m = 1`.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: I = Q._indices
            sage: Q._repr_term( I.gen((1,1)) * I.gen((4,1)) )
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

            sage: Q = QSystem(QQ, ['A',4])
            sage: I = Q._indices
            sage: Q._latex_term( I.gen((3,1)) * I.gen((4,1)) )
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

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: Q.cartan_type()
            ['A', 4]
        """
        return self._cartan_type

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: Q.index_set()
            (1, 2, 3, 4)
        """
        return self._cartan_type.index_set()

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

    @cached_method
    def one_basis(self):
        """
        Return the basis element indexing `1`.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: Q.one_basis()
            1
            sage: Q.one_basis().parent() is Q._indices
            True
        """
        return self._indices.one()

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',4])
            sage: Q.algebra_generators()
            Finite family {1: Q^(1)[1], 2: Q^(2)[1], 3: Q^(3)[1], 4: Q^(4)[1]}
        """
        return Family({a: self.gen(a, 1) for a in self._cartan_type.index_set()})

    def dimension(self):
        """
        Return the dimension of ``self``, which is `\infty`.

        EXAMPLES::

            sage: F = QSystem(QQ, ['A',4])
            sage: F.dimension()
            +Infinity
        """
        return infinity

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

    def gen(self, a, m):
        """
        Return the generator `Q^{(a)}_i` of ``self``.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',8])
            sage: Q.gen(2, 1)
            Q^(2)[1]
            sage: Q.gen(6, 2)
            -Q^(5)[1]*Q^(7)[1] + Q^(6)[1]^2
            sage: Q.gen(7, 3)
            -Q^(5)[1]*Q^(7)[1] + Q^(5)[1]*Q^(8)[1]^2 + Q^(6)[1]^2
             - 2*Q^(6)[1]*Q^(7)[1]*Q^(8)[1] + Q^(7)[1]^3
            sage: Q.gen(1, 0)
            1
        """
        if a not in self._cartan_type.index_set():
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

        This becomes our relation used for reducing the Q-system to the
        fundamental representations.

        .. NOTE::

            This helper method is defined in order to use the
            division implemented in polynomial rings.

        EXAMPLES::

            sage: Q = QSystem(QQ, ['A',8])
            sage: Q._Q_poly(1, 2)
            q1^2 - q2
            sage: Q._Q_poly(3, 2)
            q3^2 - q2*q4
            sage: Q._Q_poly(6, 3)
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

                sage: Q = QSystem(QQ, ['A',8])
                sage: x = Q.gen(1, 2)
                sage: y = Q.gen(3, 2)
                sage: x * y
                -Q^(1)[1]^2*Q^(2)[1]*Q^(4)[1] + Q^(1)[1]^2*Q^(3)[1]^2
                 + Q^(2)[1]^2*Q^(4)[1] - Q^(2)[1]*Q^(3)[1]^2
            """
            return self.parent().sum_of_terms((tl*tr, cl*cr)
                                              for tl,cl in self for tr,cr in x)

