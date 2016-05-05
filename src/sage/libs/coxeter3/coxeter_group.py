"""
Coxeter Groups implemented with Coxeter3
"""
#*****************************************************************************
#       Copyright (C) 2009-2013 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.coxeter3.coxeter import get_CoxGroup, CoxGroupElement
from sage.misc.cachefunc import cached_method

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.all import CoxeterGroups
from sage.structure.parent import Parent

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class CoxeterGroup(UniqueRepresentation, Parent):
    @staticmethod
    def __classcall__(cls, cartan_type, *args, **options):
        """
        TESTS::

            sage: from sage.libs.coxeter3.coxeter_group import CoxeterGroup # optional - coxeter3
            sage: CoxeterGroup(['B',2])                                     # optional - coxeter3
            Coxeter group of type ['B', 2] implemented by Coxeter3

        """
        from sage.combinat.all import CartanType
        ct = CartanType(cartan_type)
        return super(CoxeterGroup, cls).__classcall__(cls, ct, *args, **options)

    def __init__(self, cartan_type):
        """
        TESTS::

            sage: from sage.libs.coxeter3.coxeter_group import CoxeterGroup  # optional - coxeter3
            sage: CoxeterGroup(['A',2])                                     # optional - coxeter3
            Coxeter group of type ['A', 2] implemented by Coxeter3
            sage: TestSuite(CoxeterGroup(['A',2])).run()                    # optional - coxeter3
        """
        category = CoxeterGroups()
        if cartan_type.is_finite():
            category = category.Finite()
        Parent.__init__(self, category=category)
        self._coxgroup = get_CoxGroup(cartan_type)
        self._cartan_type = cartan_type

    def _repr_(self):
        """
        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3'); W      # optional - coxeter3 # indirect doctest
            Coxeter group of type ['A', 3] implemented by Coxeter3
            sage: W = CoxeterGroup(['A', 3, 1], implementation='coxeter3'); W   # optional - coxeter3
            Coxeter group of type ['A', 3, 1] implemented by Coxeter3
        """
        return "Coxeter group of type %s implemented by Coxeter3"%(self.cartan_type())

    def __iter__(self):
        """
        EXAMPLES::

            sage: W = CoxeterGroup(['A', 2], implementation='coxeter3')    # optional - coxeter3
            sage: list(W)                                                  # optional - coxeter3
            [[], [1], [2], [1, 2], [2, 1], [1, 2, 1]]
        """
        for x in self._coxgroup:
            yield CoxeterGroup.Element(self, x)

    def cartan_type(self):
        """
        Return the Cartan type for this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: W.cartan_type()                                         # optional - coxeter3
            ['A', 3]
        """
        return self._cartan_type

    def index_set(self):
        """
        Return the index set for the generators of this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')  # optional - coxeter3
            sage: W.index_set()                                          # optional - coxeter3
            (1, 2, 3)
            sage: C = CoxeterGroup(['A', 3,1], implementation='coxeter3') # optional - coxeter3
            sage: C.index_set()                                           # optional - coxeter3
            (0, 1, 2, 3)
        """
        return self.cartan_type().index_set()
        #return range(1, self.rank()+1)

    def bruhat_interval(self, u, v):
        """
        Return the Bruhat interval between ``u`` and ``v``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')    # optional - coxeter3
            sage: W.bruhat_interval([1],[3,1,2,3])                         # optional - coxeter3
            [[1], [1, 2], [1, 3], [1, 2, 3], [1, 3, 2], [1, 2, 3, 2]]
        """
        u, v = self(u), self(v)
        return self._coxgroup.bruhat_interval(u.value, v.value)

    def cardinality(self):
        """
        Return the cardinality of this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')    # optional - coxeter3
            sage: W.cardinality()                                          # optional - coxeter3
            24
        """
        return self._coxgroup.order()

    def one(self):
        """
        Return the identity element of this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: W.one()                                                 # optional - coxeter3
            []

        """
        return self.element_class(self, [])

    def simple_reflections(self):
        """
        Return the family of generators for this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: s = W.simple_reflections()                                            # optional - coxeter3
            sage: s[2]*s[1]*s[2]                                          # optional - coxeter3
            [2, 1, 2]
        """
        from sage.combinat.family import Family
        return Family(self.index_set(), lambda i: self.element_class(self, [i]))

    gens = simple_reflections

    def rank(self):
        """
        Return the rank of this Coxeter group, that is, the number of generators.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: W.rank()                                                # optional - coxeter3
            3
        """
        return self._coxgroup.rank()

    def is_finite(self):
        """
        Return True if this is a finite Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: W.is_finite()                                           # optional - coxeter3
            True
        """
        return self._coxgroup.is_finite()

    def length(self, x):
        """
        Return the length of an element ``x`` in this Coxeter group.
        This is just the length of a reduced word for ``x``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: W.length(W([1,2]))                                      # optional - coxeter3
            2
            sage: W.length(W([1,1]))                                      # optional - coxeter3
            0

        """
        return x.length()

    @cached_method
    def coxeter_matrix(self):
        """
        Return the Coxeter matrix for this Coxeter group.

        The columns and rows are ordered according to the result of
        :meth:`index_set`.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')    # optional - coxeter3
            sage: W.coxeter_matrix()                                       # optional - coxeter3
            [1 3 2]
            [3 1 3]
            [2 3 1]

        """
        return self._coxgroup.coxeter_matrix()

    def root_system(self):
        """
        Return the root system associated with this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: R = W.root_system(); R                                  # optional - coxeter3
            Root system of type ['A', 3]
            sage: alpha = R.root_space().basis()                          # optional - coxeter3
            sage: alpha[2] + alpha[3]                                     # optional - coxeter3
            alpha[2] + alpha[3]
        """
        return self.cartan_type().root_system()

    def _an_element_(self):
        """
        Return an element of this Coxeter group.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: W._an_element_()                                        # optional - coxeter3
            []

        """
        return self.element_class(self, [])

    def m(self, i, j):
        """
        Return the entry in the Coxeter matrix between the generator
        with label ``i`` and the generator with label ``j``.

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: W.m(1,1)                                                # optional - coxeter3
            1
            sage: W.m(1,0)                                                # optional - coxeter3
            2
        """
        return self.coxeter_matrix()[i-1,j-1]

    def kazhdan_lusztig_polynomial(self, u, v, constant_term_one=True):
        r"""
        Return the Kazhdan-Lusztig polynomial `P_{u,v}`.

        INPUT:

        - ``u``, ``v`` -- elements of the underlying Coxeter group
        - ``constant_term_one`` -- (default: True) True uses the constant equals one convention,
           False uses the Leclerc-Thibon convention

        .. SEEALSO::

            - :class:`~sage.combinat.kazhdan_lusztig.KazhdanLusztigPolynomial`
            - :meth:`parabolic_kazhdan_lusztig_polynomial`

        EXAMPLES::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
            sage: W.kazhdan_lusztig_polynomial([], [1,2, 1])              # optional - coxeter3
            1
            sage: W.kazhdan_lusztig_polynomial([1],[3,2])                 # optional - coxeter3
            0
            sage: W = CoxeterGroup(['A',3],implementation='coxeter3')     # optional - coxeter3
            sage: W.kazhdan_lusztig_polynomial([2],[2,1,3,2])             # optional - coxeter3
            q + 1

        .. NOTE::

            Coxeter3, as well as Sage's native implementation in
            :class:`~sage.combinat.kazhdan_lusztig.KazhdanLusztigPolynomial`
            use the convention under which Kazhdan-Lusztig
            polynomials give the change of basis from the `(C_w)_{w\in W}`
            basis to the `(T_w)_{w\in W}` of the Hecke algebra of `W` with
            parameters `q` and `q^{-1}`:

                .. MATH:: C_w = \sum_u  P_{u,w} T_w

            In particular, `P_{u,u}=1`::

                sage: all(W.kazhdan_lusztig_polynomial(u,u) == 1 for u in W) # optional - coxeter3
                True

            This convention differs from Theorem 2.7 in [LeclercThibon1998]_ by:

            .. MATH::

                {}^{LT} P_{y,w}(q) = q^{\ell(w)-\ell(y)} P_{y,w}(q^{-2})

            To access the Leclerc-Thibon convention use::

                sage: W = CoxeterGroup(['A',3],implementation='coxeter3')                         # optional - coxeter3
                sage: W.kazhdan_lusztig_polynomial([2],[2,1,3,2],constant_term_one=False)         # optional - coxeter3
                q^3 + q

        TESTS:

        We check that Coxeter3 and Sage's implementation give the same results::

            sage: C = CoxeterGroup(['B', 3], implementation='coxeter3')                           # optional - coxeter3
            sage: W = WeylGroup("B3",prefix="s")
            sage: [s1,s2,s3] = W.simple_reflections()
            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: KL = KazhdanLusztigPolynomial(W,q)
            sage: all(KL.P(1,w) == C.kazhdan_lusztig_polynomial([],w.reduced_word()) for w in W)  # optional - coxeter3  # long (15s)
            True
        """
        u, v = self(u), self(v)
        p = u.value.kazhdan_lusztig_polynomial(v.value)
        if constant_term_one:
            return p
        ZZq = PolynomialRing(ZZ, 'q', sparse=True)
        # This is the same as q**len_diff * p(q**(-2))
        len_diff = v.length()-u.length()
        d = {-2*deg+len_diff: coeff for deg,coeff in enumerate(p) if coeff != 0}
        return ZZq(d)

    def parabolic_kazhdan_lusztig_polynomial(self, u, v, J, constant_term_one=True):
        """
        Return the parabolic Kazhdan-Lusztig polynomial `P_{u,v}^{-,J}`.

        INPUT:

        - ``u``, ``v`` -- minimal length coset representatives of `W/W_J` for this Coxeter group `W`
        - ``J`` -- a subset of the index set of ``self`` specifying the parabolic subgroup

        This method implements the parabolic Kazhdan-Lusztig polynomials
        `P^{-,J}_{u,v}` of [Deodhar1987]_, which are defined as
        `P^{-,J}_{u,v} = \sum_{z\in W_J} (-1)^{\ell(z)} P_{yz,w}(q)`
        with the conventions in Sage.
        As for :meth:`kazhdan_lusztig_polynomial` the convention
        differs from Theorem 2.7 in [LeclercThibon1998]_ by:

        .. MATH::

            {}^{LT} P_{y,w}^{-,J}(q) = q^{\ell(w)-\ell(y)} P_{y,w}^{-,J}(q^{-2})

        REFERENCES:

            .. [Deodhar1987] \V.V. Deodhar, On some geometric aspects of Bruhat orderings II. The parabolic analogue of Kazhdan-Lusztig polynomials, J. Alg. 111 (1987) 483-506.

            .. [LeclercThibon1998] \B. Leclerc, J.-Y. Thibon, Littlewood-Richardson coefficients and Kazhdan-Lusztig polynomials, http://front.math.ucdavis.edu/9809.5122

        EXAMPLES::

            sage: W = CoxeterGroup(['A',3], implementation='coxeter3')                      # optional - coxeter3
            sage: W.parabolic_kazhdan_lusztig_polynomial([],[3,2],[1,3])                    # optional - coxeter3
            0
            sage: W.parabolic_kazhdan_lusztig_polynomial([2],[2,1,3,2],[1,3])               # optional - coxeter3
            q

            sage: C = CoxeterGroup(['A',3,1], implementation='coxeter3')                    # optional - coxeter3
            sage: C.parabolic_kazhdan_lusztig_polynomial([],[1],[0])                        # optional - coxeter3
            1
            sage: C.parabolic_kazhdan_lusztig_polynomial([],[1,2,1],[0])                    # optional - coxeter3
            1
            sage: C.parabolic_kazhdan_lusztig_polynomial([],[0,1,0,1,2,1],[0])              # optional - coxeter3
            q
            sage: w=[1, 2, 1, 3, 0, 2, 1, 0, 3, 0, 2]
            sage: v=[1, 2, 1, 3, 0, 1, 2, 1, 0, 3, 0, 2, 1, 0, 3, 0, 2]
            sage: C.parabolic_kazhdan_lusztig_polynomial(w,v,[1,3])                         # optional - coxeter3
            q^2 + q
            sage: C.parabolic_kazhdan_lusztig_polynomial(w,v,[1,3],constant_term_one=False) # optional - coxeter3
            q^4 + q^2

        TESTS::

            sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')                     # optional - coxeter3
            sage: type(W.parabolic_kazhdan_lusztig_polynomial([2],[],[1]))                  # optional - coxeter3
            <type 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
        """
        u = self(u)
        v = self(v)
        if any(d in J for d in u.descents()) or any(d in J for d in v.descents()):
            raise ValueError("u and v have to be minimal coset representatives")
        J_set = set(J)
        WOI = self.weak_order_ideal(lambda x: J_set.issuperset(x.descents()))
        if constant_term_one:
            P = PolynomialRing(ZZ, 'q')
            return P.sum((-1)**(z.length()) * self.kazhdan_lusztig_polynomial(u*z,v)
                         for z in WOI if (u*z).bruhat_le(v))
        P = PolynomialRing(ZZ, 'q', sparse=True)
        return P.sum((-1)**(z.length()) * self.kazhdan_lusztig_polynomial(u*z,v, constant_term_one=False).shift(z.length())
                     for z in WOI if (u*z).bruhat_le(v))


    class Element(ElementWrapper):
        wrapped_class = CoxGroupElement

        def __init__(self, parent, x):
            """
            TESTS::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')    # optional - coxeter3
                sage: W([2,1,2])                                               # optional - coxeter3
                [1, 2, 1]
            """
            if not isinstance(x, CoxGroupElement):
                x = CoxGroupElement(parent._coxgroup, x).reduced()
            ElementWrapper.__init__(self, parent, x)

        def __iter__(self):
            """
            Return an iterator for the elements in the reduced word.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
                sage: w = W([1,2,1])                                          # optional - coxeter3
                sage: list(iter(w))                                           # optional - coxeter3
                [1, 2, 1]
            """
            return iter(self.value)

        def coatoms(self):
            """
            Return the coatoms (or co-covers) of this element in the Bruhat order.

            EXAMPLES::

                sage: W = CoxeterGroup(['B', 3], implementation='coxeter3')  # optional - coxeter3
                sage: w = W([1,2,3])                                         # optional - coxeter3
                sage: w.coatoms()                                            # optional - coxeter3
                [[2, 3], [1, 3], [1, 2]]
            """
            W = self.parent()
            return [W(w) for w in self.value.coatoms()]

        def __cmp__(self, other):
            """
            Return lexicographic comparison of ``self`` and ``other``.

            EXAMPLES::

                sage: W = CoxeterGroup(['B', 3], implementation='coxeter3')   # optional - coxeter3
                sage: w = W([1,2,3])                                          # optional - coxeter3
                sage: v = W([3,1,2])                                          # optional - coxeter3
                sage: w.__cmp__(v)                                            # optional - coxeter3
                -1
                sage: v.__cmp__(w)                                            # optional - coxeter3
                1
            """
            if type(self) is not type(other):
                return cmp(type(self), type(other))
            return cmp(list(self), list(other))

        def reduced_word(self):
            """
            Return the reduced word of ``self``.

            EXAMPLES::

                sage: W = CoxeterGroup(['B', 3], implementation='coxeter3')  # optional - coxeter3
                sage: w = W([1,2,3])                                         # optional - coxeter3
                sage: w.reduced_word()                                       # optional - coxeter3
                [1, 2, 3]
            """
            return list(self)

        def __invert__(self):
            """
            Return the inverse of this Coxeter group element.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
                sage: w = W([1,2,3])                                          # optional - coxeter3
                sage: ~w                                                      # optional - coxeter3
                [3, 2, 1]
            """
            return self.__class__(self.parent(), ~self.value)

        inverse = __invert__

        def __getitem__(self, i):
            """
            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
                sage: w0 = W([1,2,1])                                         # optional - coxeter3
                sage: w0[0]                                                   # optional - coxeter3
                1
                sage: w0[1]                                                   # optional - coxeter3
                2

            """
            # Allow the error message to be raised by the underlying element
            return self.value[i]

        def _mul_(self, y):
            """
            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
                sage: s = W.gens()                                            # optional - coxeter3
                sage: s[1]._mul_(s[1])                                        # optional - coxeter3
                []
                sage: s[1]*s[2]*s[1]                                          # optional - coxeter3
                [1, 2, 1]
            """
            return self.__class__(self.parent(), self.value * y.value)

        def __eq__(self, y):
            """
            Return whether this Coxeter group element is equal to ``y``.

            This is computed by computing the normal form of both elements.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')    # optional - coxeter3
                sage: W([1,2,1]) == W([2,1,2])                                 # optional - coxeter3
                True
                sage: W([1,2,1]) == W([2,1])                                   # optional - coxeter3
                False

            """
            if not isinstance(y, self.parent().element_class):
                return False

            return list(self) == list(y)

        def __len__(self):
            """
            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')    # optional - coxeter3
                sage: w = W([1,2,1])                                  # optional - coxeter3
                sage: w.length()                                      # optional - coxeter3
                3
                sage: len(w)                                          # optional - coxeter3
                3
            """
            return len(self.value)

        length = __len__

        def bruhat_le(self, v):
            """
            Return whether ``self`` `\le` ``v`` in Bruhat order.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
                sage: W([]).bruhat_le([1,2,1])                                # optional - coxeter3
                True
            """
            v = self.parent()(v)
            return self.value.bruhat_le(v.value)

        def poincare_polynomial(self):
            """
            Return the Poincare polynomial associated with this element.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 2], implementation='coxeter3')      # optional - coxeter3
                sage: W.long_element().poincare_polynomial()                     # optional - coxeter3
                t^3 + 2*t^2 + 2*t + 1
                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')      # optional - coxeter3
                sage: W([2,1,3,2]).poincare_polynomial()                         # optional - coxeter3
                t^4 + 4*t^3 + 5*t^2 + 3*t + 1
                sage: W([1,2,3,2,1]).poincare_polynomial()                       # optional - coxeter3
                t^5 + 4*t^4 + 6*t^3 + 5*t^2 + 3*t + 1

                sage: rw = sage.combinat.permutation.from_reduced_word           # optional - coxeter3
                sage: p = map(attrcall('poincare_polynomial'), W)                # optional - coxeter3
                sage: [rw(w.reduced_word()) for i,w in enumerate(W) if p[i] != p[i].reverse()] # optional - coxeter3
                [[3, 4, 1, 2], [4, 2, 3, 1]]
            """
            return self.value.poincare_polynomial()

        def has_right_descent(self, i):
            """
            Return whether ``i`` is a right descent of this element.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 4], implementation='coxeter3')   # optional - coxeter3
                sage: W([1,2]).has_right_descent(1)                           # optional - coxeter3
                False
                sage: W([1,2]).has_right_descent(2)                           # optional - coxeter3
                True
            """
            return i in self.value.right_descents()

        def has_left_descent(self, i):
            """
            Return True if ``i`` is a left descent of this element.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 4], implementation='coxeter3')   # optional - coxeter3
                sage: W([1,2]).has_left_descent(1)                            # optional - coxeter3
                True
                sage: W([1,2]).has_left_descent(2)                            # optional - coxeter3
                False
            """
            return i in self.value.left_descents()

        def action(self, v):
            """
            Return the action of of this Coxeter group element on the root space.

            INPUT:

            - ``v`` -- an element of the root space associated with the Coxeter group for ``self``

            EXAMPLES::

                sage: W = CoxeterGroup(['B', 3], implementation='coxeter3')   # optional - coxeter3
                sage: R = W.root_system().root_space()                        # optional - coxeter3
                sage: v = R.an_element(); v                                   # optional - coxeter3
                2*alpha[1] + 2*alpha[2] + 3*alpha[3]
                sage: w = W([1,2,3])                                          # optional - coxeter3
                sage: w.action(v)                                             # optional - coxeter3
                -alpha[1] + alpha[2] + alpha[3]
            """
            #TODO: Find a better way to do this
            W = self.parent().root_system().root_space().weyl_group()
            w = W.from_reduced_word(list(self))
            return w.action(v)

        def action_on_rational_function(self, f):
            r"""
            Return the natural action of this Coxeter group element on a
            polynomial considered as an element of `S(\mathfrak{h}^*)`.

            .. NOTE::

               Note that the number of variables in the polynomial
               ring must correspond to the rank of this Coxeter
               group. The ordering of the variables is assumed to
               coincide with the result of :meth:`index_set`.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3], implementation='coxeter3')   # optional - coxeter3
                sage: S = PolynomialRing(QQ, 'x,y,z').fraction_field()        # optional - coxeter3
                sage: x,y,z = S.gens()                                        # optional - coxeter3
                sage: W([1]).action_on_rational_function(x+y+z)               # optional - coxeter3
                (x^2*y + x*z + 1)/x
                sage: W([2]).action_on_rational_function(x+y+z)               # optional - coxeter3
                (x*y^2 + y^2*z + 1)/y
                sage: W([3]).action_on_rational_function(x+y+z)               # optional - coxeter3
                (y*z^2 + x*z + 1)/z
            """
            Q = f.parent()
            Q_gens = Q.gens()
            W = self.parent()
            R = W.root_system().root_space()
            alpha = R.basis()
            n = W.rank()

            if Q.ngens() != n:
                raise ValueError("the number of generators for the polynomial ring must be the same as the rank of the root system")

            basis_elements = [alpha[i] for i in W.index_set()]
            basis_to_order = {s: i for i, s in enumerate(W.index_set())}

            results = []
            for poly in [f.numerator(), f.denominator()]:
                result = 0
                exponents = poly.exponents()

                for exponent in exponents:
                    #Construct something in the root lattice from the exponent vector
                    exponent = sum(e*b for e, b in zip(exponent, basis_elements))
                    exponent = self.action(exponent)

                    monomial = 1
                    for s, c in exponent.monomial_coefficients().iteritems():
                        monomial *= Q_gens[basis_to_order[s]]**int(c)

                    result += monomial

                results.append(result)

            numerator, denominator = results
            return numerator / denominator

