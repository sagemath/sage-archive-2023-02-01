from sage.misc.cachefunc import cached_method
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement

class QuiverAlgebra(CombinatorialFreeModule):
    """
    Creates the Quiver Algebra of a Quiver over a given field.

    Given a Quiver Q and a field k the quiver algebra kQ is defined as follows.  As
    a vector space it has basis the set of all paths in Q.  Multiplication is
    defined on this basis and extended bilinearly.  If p is a path with terminal
    vertex t and q is a path with initial vertex i then the product p*q is defined
    to be the composition of the paths p and q if t = i and 0 otherwise.

    INPUT:

    - ``k`` - field, the base field of the quiver algebra.

    - ``Q`` - Quiver or DiGraph, the quiver of the quiver algebra.

    OUTPUT:

        - QuiverAlgebra

    EXAMPLES::

        sage: Q = Quiver({1:{2:['a']}, 2:{3:['b']}})
        sage: A = Q.algebra(GF(7))
        sage: A
        Algebra of Quiver on 3 vertices over Finite Field of size 7
        sage: A.variable_names()
        ('e_1', 'e_2', 'e_3', 'a', 'b')

    Note that QuiverAlgebras are uniquely defined by their Quiver and field::

        sage: A is Q.algebra(GF(7))
        True
        sage: A is Q.algebra(RR)
        False
        sage: A is Quiver({1:{2:['a']}}).algebra(GF(7))
        False

    The QuiverAlgebra of an acyclic quiver has a finite basis::

        sage: A.dimension()
        6
        sage: list(A.basis())
        [e_1, e_2, e_3, a, b, a*b]

    The QuiverAlgebra can create elements from QuiverPaths or from elements of the
    base ring::

        sage: A(5)
        5*e_1 + 5*e_2 + 5*e_3
        sage: S = A.semigroup()
        sage: S
        Free small category of Quiver on 3 vertices
        sage: p = S((1, 2, 'a'))
        sage: r = S((2, 3, 'b'))
        sage: e2 = S((2, 2))
        sage: x = A(p) + A(e2)
        sage: x
        a + e_2
        sage: y = A(p) + A(r)
        sage: y
        a + b

    QuiverAlgebras are graded algebras.  The grading is given by assigning to each
    basis element the length of the path corresponding to that basis element::

        sage: x.is_homogeneous()
        False
        sage: x.degree()
        Traceback (most recent call last):
        ...
        ValueError: Element is not homogeneous.
        sage: y.is_homogeneous()
        True
        sage: y.degree()
        1
        sage: A[1]
        Free module spanned by [a, b] over Finite Field of size 7
        sage: A[2]
        Free module spanned by [a*b] over Finite Field of size 7

    TESTS::

        sage: TestSuite(A).run()

    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, k, P):
        """
        Creates a QuiverAlgebra object.  Type QuiverAlgebra? for more information.

        INPUT:

        - ``k``, a commutagive ring
        - ``P``, the partial semigroup formed by the paths of a quiver.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}}).path_semigroup()
            sage: Q.algebra(GF(5))
            Algebra of Quiver on 4 vertices over Finite Field of size 5
        """
        # The following hidden methods are relevant:
        #
        # - _base
        #       The base ring of the quiver algebra.
        # - _basis_keys
        #       Finite enumerated set containing the QuiverPaths that form the
        #       basis.
        # - _quiver
        #       The quiver of the quiver algebra
        # - _free_small_category
        #       Shortcut for _quiver.free_small_category()

        from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
        self._quiver = P.quiver()
        self._semigroup = P
        super(QuiverAlgebra, self).__init__(k, self._semigroup,
                                            prefix='',
                                            element_class=self.Element,
                                            category=GradedAlgebrasWithBasis(k),
                                            bracket=False)
        self._assign_names(self._semigroup.variable_names())

    @cached_method
    def gens(self):
        """
        Generators of this algebra (idempotents and arrows).

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}})
            sage: A = Q.algebra(GF(5))
            sage: A.variable_names()
            ('e_1', 'e_2', 'e_3', 'e_4', 'a', 'b', 'c')
            sage: A.gens()
            (e_1, e_2, e_3, e_4, a, b, c)

        """
        return tuple(self._from_dict( {index: self.base_ring().one()}, remove_zeros = False ) for index in self._semigroup.gens())

    @cached_method
    def arrows(self):
        """
        Arrows of this algebra (corresponding to edges of the underlying quiver).

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}})
            sage: A = Q.algebra(GF(5))
            sage: A.arrows()
            (a, b, c)

        """
        return tuple(self._from_dict( {index: self.base_ring().one()}, remove_zeros = False ) for index in self._semigroup.arrows())

    @cached_method
    def idempotents(self):
        """
        Idempotents of this algebra (corresponding to vertices of the underlying quiver).

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}})
            sage: A = Q.algebra(GF(5))
            sage: A.idempotents()
            (e_1, e_2, e_3, e_4)

        """
        return tuple(self._from_dict( {index: self.base_ring().one()}, remove_zeros = False ) for index in self._semigroup.idempotents())


    def gen(self, i):
        """
        `i`-th generator of this algebra.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}})
            sage: A = Q.algebra(GF(5))
            sage: A.gens()
            (e_1, e_2, e_3, e_4, a, b, c)
            sage: A.gen(2)
            e_3
            sage: A.gen(5)
            b

        """
        return self._from_dict( {self._semigroup.gen(i): self.base_ring().one()}, remove_zeros = False )

    def ngens(self):
        """
        Number of generators of this algebra.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}})
            sage: A = Q.algebra(GF(5))
            sage: A.ngens()
            7

        """
        return self._semigroup.ngens()

    def _element_constructor_(self, x):
        """
        Attempts to construct an element of self from x.

        TESTS::

            sage: A = Quiver({1:{2:['a']}}).algebra(QQ)
            sage: B = Quiver({1:{2:['a']}, 2:{3:['b']}}).algebra(QQ)
            sage: x = A((1, 2, 'a')) + 1 # indirect doctest
            sage: x
            a + e_1 + e_2
            sage: B(x) # indirect doctest
            a + e_1 + e_2
            sage: A(1) # indirect doctest
            e_1 + e_2
        """
        from sage.quivers.paths import QuiverPath
        # If it's an element of another quiver algebra, do a linear combination
        # of the basis
        if isinstance(x, CombinatorialFreeModuleElement) and isinstance(x.parent(), QuiverAlgebra):
            coeffs = x.monomial_coefficients()
            result = self.zero()
            for key in coeffs:
                result += coeffs[key]*self.monomial(key)
            return result

        # If it's a QuiverPath return the associated basis element
        if isinstance(x, QuiverPath):
            if x:
                return self.monomial(x)
            else:
                return self.zero()

        # If it's a tuple or a list try and create a QuiverPath from it and
        # then return the associated basis element
        if isinstance(x, tuple) or isinstance(x, list):
            return self.monomial(self._semigroup(x))

        # Otherwise let CombinatorialFreeModule try
        return super(QuiverAlgebra, self)._element_constructor_(x)

    def _coerce_map_from_(self, other):
        """
        True if there is a coercion from other to self.

        The algebras that coerce into a quiver algebra are rings k or quiver algebras
        kQ such that k has a coercion into the base ring of self and Q is a subquiver
        of the quiver of self.

        In addition, the free algebra of a subquiver coerces into the algebra.

        TESTS::

            sage: Q1 = Quiver({1:{2:['a']}})
            sage: Q2 = Quiver({1:{2:['a','b']}})
            sage: A1 = Q1.algebra(GF(3))
            sage: A2 = Q2.algebra(GF(3))
            sage: A1.coerce_map_from(A2) # indirect doctest
            sage: A2.coerce_map_from(A1) # indirect doctest
            Conversion map:
                  From: Algebra of Quiver on 2 vertices over Finite Field of size 3
                  To:   Algebra of Quiver on 2 vertices over Finite Field of size 3
            sage: A1.coerce_map_from(ZZ) # indirect doctest
            Composite map:
                  From: Integer Ring
                  To:   Algebra of Quiver on 2 vertices over Finite Field of size 3
                  Defn:   Natural morphism:
                          From: Integer Ring
                          To:   Finite Field of size 3
                        then
                          Generic morphism:
                          From: Finite Field of size 3
                          To:   Algebra of Quiver on 2 vertices over Finite Field of size 3
            sage: A1.coerce_map_from(QQ) # indirect doctest
            sage: A1.coerce_map_from(ZZ)
            Composite map:
              From: Integer Ring
              To:   Algebra of Quiver on 2 vertices over Finite Field of size 3
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Finite Field of size 3
                    then
                      Generic morphism:
                      From: Finite Field of size 3
                      To:   Algebra of Quiver on 2 vertices over Finite Field of size 3

        ::

            sage: A2.coerce_map_from(Q1.path_semigroup())
            Conversion map:
              From: Free small category of Quiver on 2 vertices
              To:   Algebra of Quiver on 2 vertices over Finite Field of size 3
            sage: a = Q1.path_semigroup()(Q1.path_semigroup().arrows()[0]); a
            a
            sage: A2.one() * a == a     # indirect doctest
            True

        """

        if (isinstance(other, QuiverAlgebra) and
                self._base.has_coerce_map_from(other._base) and
                other._quiver <= self._quiver):
            return True
        if self._semigroup.has_coerce_map_from(other):
            return True
        return self._base.has_coerce_map_from(other)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b']}})
            sage: Q.algebra(RR) # indirect doctest
            Algebra of Quiver on 3 vertices over Real Field with 53 bits of precision
        """

        return "Algebra of {0} over {1}".format(self._quiver, self._base)

    ###########################################################################
    #                                                                         #
    # CATEGORY METHODS                                                        #
    #    These functions are used by the category to implement the algebra    #
    #    structure.                                                           #
    #                                                                         #
    ###########################################################################

    def _monomial(self, index):
        """
        This method makes sure that the invalid path evaluates as zero.

        TESTS::

            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: inv = Q([(1,3,'x')])
            sage: inv
            invalid path
            sage: Q.algebra(ZZ)(inv)          # indirect doctest
            0

        """
        if index:
            return self._from_dict( {index: self.base_ring().one()}, remove_zeros = False )
        return self.zero()

    def product_on_basis(self, p1, p2):
        """
        Returns the element corresponding to p1*p2 in the quiver algebra.

        INPUT:

        - ``p1``, ``p2`` - QuiverPaths

        OUTPUT:

        - CombinatorialFreeModuleElement

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p1 = Q((1, 2, 'a'))
            sage: p2 = Q([(2, 3, 'b'), (3, 4, 'c')])
            sage: A = Q.algebra(QQ)
            sage: A.product_on_basis(p1, p2)
            a*b*c
            sage: A.product_on_basis(p2, p1)
            0

        """
        FSC = self._semigroup
        p = FSC(p1)*FSC(p2)
        if p:
            return self.basis()[p]
        else:
            return self.zero()

    def degree_on_basis(self, p):
        """
        Return the degree of the monomial specified by p.

        EXAMPLES::

            sage: A = Quiver({1:{2:['a']}, 2:{3:['b']}}).algebra(QQ)
            sage: A.degree_on_basis((1, 1))
            0
            sage: A.degree_on_basis((1, 2, 'a'))
            1
            sage: A.degree_on_basis([(1, 2, 'a'), (2, 3, 'b')])
            2
        """

        return len(self._semigroup(p))

    def one(self):
        """
        Return the multiplicative identity element.

        The multiplicative identity of a quiver algebra is the sum of the basis
        elements corresponding to the trivial paths at each vertex.

        EXAMPLES::

            sage: A = Quiver({1:{2:['a']}, 2:{3:['b']}}).algebra(QQ)
            sage: A.one()
            e_1 + e_2 + e_3
        """

        x = self.zero()
        B = self.basis()
        for v in self._quiver:
            x += B[self._semigroup([(v, v)], check=False)]
        return x

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data and subspaces of the quiver algebra.     #
    #                                                                         #
    ###########################################################################

    def quiver(self):
        """
        Return the quiver of this algebra.

        OUTPUT:

        - Quiver

        EXAMPLES:

            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: A = Q.algebra(GF(3))
            sage: A.quiver() is Q
            True
        """
        return self._quiver

    def semigroup(self):
        """
        Return the semigroup of this algebra.

        NOTE:

        The semigroup is the same as the free small category that is
        associated with the underlying quiver.

        OUTPUT:

        - Quiver

        EXAMPLES:

            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: A = Q.algebra(GF(3))
            sage: A.semigroup() is Q.path_semigroup()
            True
        """
        return self._semigroup

    def homogeneous_component(self, n):
        """
        Return the nth homogeneous piece of the quiver algebra.

        INPUT:

        - ``n`` - integer

        OUTPUT:

        - CombinatorialFreeModule, module spanned by the paths of lenth n in the
          quiver.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{4:['c']}, 3:{4:['d']}})
            sage: A = Q.algebra(GF(7))
            sage: A.homogeneous_component(2)
            Free module spanned by [a*c, b*d] over Finite Field of size 7
        """

        basis = [p for p in self._basis_keys if len(p) == n]
        M = CombinatorialFreeModule(self._base, basis, prefix='', bracket=False)
        M._name = "Free module spanned by {0}".format(basis)
        return M

    __getitem__ = homogeneous_component

    ###########################################################################
    #                                                                         #
    # ELEMENT CLASS                                                           #
    #    The class of elements of the quiver algebra.                         #
    #                                                                         #
    ###########################################################################

    class Element(CombinatorialFreeModuleElement):
        def is_homogeneous(self):
            """
            Return True if and only if this element is homogeneous.

            EXAMPLES::

                sage: A = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}}).algebra(QQ)
                sage: (A((1, 2, 'a')) + A((1, 2, 'b'))).is_homogeneous()
                True
                sage: (A((1, 1)) + A((1, 2, 'a'))).is_homogeneous()
                False
            """

            # Get the support, the zero element is homogeneous
            paths = self.support()
            if not paths:
                return True

            # Compare the rest of the paths, they must be the same length
            for p in paths[1:]:
                if len(p) != len(paths[0]):
                    return False

            return True

        def degree(self):
            """
            The degree of self.

            EXAMPLES::

                sage: A = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}}).algebra(QQ)
                sage: A((1, 1)).degree()
                0
                sage: (A((1, 2, 'a')) + A((1, 2, 'b'))).degree()
                1

            An error is raised if the element is not homogeneous::

                sage: (A((1, 1)) + A((1, 2, 'a'))).degree()
                Traceback (most recent call last):
                ...
                ValueError: Element is not homogeneous.
            """

            # Deal with zero
            paths = self.support()
            if not paths:
                raise ValueError("The zero element does not have a well-defined degree.")

            # Check that the element is homogeneous
            for p in paths[1:]:
                if len(p) != len(paths[0]):
                    raise ValueError("Element is not homogeneous.")

            return len(paths[0])

