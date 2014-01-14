r"""
Quotients of finite rank free modules over a field.
"""

####################################################################################
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
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
####################################################################################

from free_module import FreeModule_ambient_field

class FreeModule_ambient_field_quotient(FreeModule_ambient_field):
    r"""
    A quotient V/W of two vector spaces as a vector space.

    To obtain V or W use \code{self.V()} and \code{self.W()}.

    EXAMPLES:
        sage: k.<i> = QuadraticField(-1)
        sage: A = k^3; V = A.span([[1,0,i], [2,i,0]])
        sage: W = A.span([[3,i,i]])
        sage: U = V/W; U
        Vector space quotient V/W of dimension 1 over Number Field in i with defining polynomial x^2 + 1 where
        V: Vector space of degree 3 and dimension 2 over Number Field in i with defining polynomial x^2 + 1
        Basis matrix:
        [ 1  0  i]
        [ 0  1 -2]
        W: Vector space of degree 3 and dimension 1 over Number Field in i with defining polynomial x^2 + 1
        Basis matrix:
        [    1 1/3*i 1/3*i]
        sage: U.V()
        Vector space of degree 3 and dimension 2 over Number Field in i with defining polynomial x^2 + 1
        Basis matrix:
        [ 1  0  i]
        [ 0  1 -2]
        sage: U.W()
        Vector space of degree 3 and dimension 1 over Number Field in i with defining polynomial x^2 + 1
        Basis matrix:
        [    1 1/3*i 1/3*i]
        sage: U.quotient_map()
        Vector space morphism represented by the matrix:
        [  1]
        [3*i]
        Domain: Vector space of degree 3 and dimension 2 over Number Field in i with defining polynomial x^2 + 1
        Basis matrix:
        [ 1  0  i]
        [ 0  1 -2]
        Codomain: Vector space quotient V/W of dimension 1 over Number Field in i with defining polynomial x^2 + 1 where
        V: Vector space of degree 3 and dimension 2 over Number Field in i with defining polynomial x^2 + 1
        Basis matrix:
        [ 1  0  i]
        [ 0  1 -2]
        W: Vector space of degree 3 and dimension 1 over Number Field in i with defining polynomial x^2 + 1
        Basis matrix:
        [    1 1/3*i 1/3*i]
        sage: Z = V.quotient(W)
        sage: Z == U
        True

    TESTS::

        sage: A = QQ^0; V = A.span([]) # corner case
        sage: W = A.span([])
        sage: U = V/W

        sage: loads(dumps(U)) == U
        True
        sage: type(loads(dumps(U)) )
        <class 'sage.modules.quotient_module.FreeModule_ambient_field_quotient_with_category'>
    """
    def __init__(self, domain, sub, quotient_matrix, lift_matrix, inner_product_matrix = None):
        """
        Create this quotient space, from the given domain, sub-module, and quotient_matrix.

        EXAMPLES:
            sage: A = QQ^5; V = A.span_of_basis([[1,0,-1,1,1], [1,-1,0,2/3,3/4]]); V
            Vector space of degree 5 and dimension 2 over Rational Field
            User basis matrix:
            [  1   0  -1   1   1]
            [  1  -1   0 2/3 3/4]
            sage: W = V.span_of_basis([V.0 - 2/3*V.1]); W
            Vector space of degree 5 and dimension 1 over Rational Field
            User basis matrix:
            [1/3 2/3  -1 5/9 1/2]

        This creates a quotient vector space, which calls the init method:
            sage: Q = V / W

        Behold the type of Q:
            sage: type(Q)
            <class 'sage.modules.quotient_module.FreeModule_ambient_field_quotient_with_category'>

        We do some consistency checks on the extra quotient and
        lifting structure of Q.
            sage: Q(V.0)
            (1)
            sage: Q( V.0 - 2/3*V.1 )
            (0)
            sage: v = Q.lift(Q.0); v
            (1, 0, -1, 1, 1)
            sage: Q( v )
            (1)
        """
        base_field = domain.base_field()
        dimension = quotient_matrix.ncols()
        sparse = domain.is_sparse()
        self.__sub = sub
        self.__domain = domain
        self.__hash = hash((domain, sub))
        FreeModule_ambient_field.__init__(self, base_field, dimension, sparse)
        self.__quo_map = domain.Hom(self)(quotient_matrix)
        self.__lift_map = self.Hom(domain)(lift_matrix)

    def _repr_(self):
        """
        Return the rather verbose string representation of this quotient space V/W.

        EXAMPLES:
        We create a quotient vector space over a finite field:
            sage: k.<a> = GF(9); A = k^3; V = A.span_of_basis([[1,0,a], [a,a,1]]); W = V.span([V.1])
            sage: Q = V/W

        Note the type:
            sage: type(Q)
            <class 'sage.modules.quotient_module.FreeModule_ambient_field_quotient_with_category'>

        The string representation mentions that this is a quotient V/W, that the quotient
        has dimension 1 and is over a finite field, and also describes V and W:
            sage: Q._repr_()
            'Vector space quotient V/W of dimension 1 over Finite Field in a of size 3^2 where\nV: Vector space of degree 3 and dimension 2 over Finite Field in a of size 3^2\nUser basis matrix:\n[1 0 a]\n[a a 1]\nW: Vector space of degree 3 and dimension 1 over Finite Field in a of size 3^2\nBasis matrix:\n[    1     1 a + 2]'
        """
        return "%s space quotient V/W of dimension %s over %s where\nV: %s\nW: %s"%(
            "Sparse vector" if self.is_sparse() else "Vector",
            self.dimension(), self.base_ring(),
            self.V(), self.W())

    def __hash__(self):
        """
        Return hash of this quotient space V/W, which is by definition the hash of
        the tuple (V,W).

        EXAMPLES:
        We compute the hash of a certain 0-dimension quotient vector space:
            sage: A = QQ^2; V = A.span_of_basis([[1,0], [1,1]]); W = V.span([V.1, V.0])
            sage: Q = V/W; Q.dimension()
            0
            sage: hash(Q)
            954887582               # 32-bit
            -5856620741060301410    # 64-bit

        The hash is just got by hashing both V and W.
            sage: hash((V, W))
            954887582             # 32-bit
            -5856620741060301410  # 64-bit
        """
        return self.__hash

    def __cmp__(self, other):
        """
        Compare self and other.

        If other is not a quotient of vector spaces, returns
        comparison of the underlying types.  If it is, return
        comparison of the pair (V,W) so that self is V/W for each of
        self and other.

        EXAMPLES:
        We create three quotient spaces and compare them:
            sage: A = QQ^2; V = A.span_of_basis([[1,0], [1,1]]);
            sage: W0 = V.span([V.1, V.0]); W1 = V.span([V.1]); W2 = V.span([V.1])
            sage: Q0 = V/W0; Q1 = V/W1; Q2 = V/W2
            sage: cmp(Q0, Q1)
            1
            sage: cmp(Q1, Q0)
            -1
            sage: cmp(Q1, Q2)
            0
            sage: cmp(Q1, 5) != 0
            True
        """
        if not isinstance(other, FreeModule_ambient_field_quotient):
            return cmp(type(self), type(other))
        return cmp((self.V(), self.W()), (other.V(), other.W()))

    def __call__(self, x):
        """
        Coerce an element into this quotient space V/W if there is a way to make
        sense of it.

        An element coerces in if it can be coerced into V, or if not at least if
        if it can be made sense of as a list of length the dimension of self.

        EXAMPLES:
        We create a 2-dimensional quotient of a 3-dimension ambient vector space.
            sage: M = QQ^3 / [[1,2,3]]

        A list of length 3 coerces into QQ^3, so it coerces into M.
            sage: M([1,2,4])
            (-1/3, -2/3)
            sage: M([1,2,3])
            (0, 0)

        A list of length 2 at least coerces into M, where here it just gives
        the corresponding linear combination of the basis for M.
            sage: M([1,2])
            (1, 2)
            sage: M.0 + 2*M.1
            (1, 2)
        """
        try:
            if x.parent() is self:
                return x
        except AttributeError:
            pass
        try:
            return FreeModule_ambient_field.__call__(self, x)
        except TypeError:
            pass
        return self._coerce_impl(self.__domain(x))

    def _coerce_impl(self, x):
        """
        Canonical coercion into this quotient space V/W.

        Elements canonically coerce into self if they canonically
        coerce into V.

        EXAMPLES:
            sage: V = QQ^3; W = V.span([[1,0,0]]); Q = V/W
            sage: Q._coerce_impl(V.0)
            (0, 0)
            sage: Q._coerce_impl(0)
            (0, 0)
            sage: Q._coerce_impl(V.1)
            (1, 0)
            sage: Q.0
            (1, 0)
            sage: Q.0 + V.1
            (2, 0)

        Here we coerce in something that is over ZZ, so it canonically coerce into V hence self.
            sage: Q._coerce_impl((ZZ^3)([1,2,3]))
            (2, 3)

        Here there is no canonical coercion:
            sage: Q._coerce_impl((GF(17)^3)([1,2,3]))
            Traceback (most recent call last):
            ...
            TypeError: Automatic coercion supported only for vectors or 0.
        """
        try:
            if x.parent() is self:
                return x
        except AttributeError:
            pass
        return self.__quo_map(self.__domain._coerce_impl(x))

    def quotient_map(self):
        """
        Given this quotient space $Q = V/W$, return the natural quotient map from V to Q.

        EXAMPLES:
            sage: M = QQ^3 / [[1,2,3]]
            sage: M.quotient_map()
            Vector space morphism represented by the matrix:
            [   1    0]
            [   0    1]
            [-1/3 -2/3]
            Domain: Vector space of dimension 3 over Rational Field
            Codomain: Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]

            sage: M.quotient_map()( (QQ^3)([1,2,3]) )
            (0, 0)
        """
        return self.__quo_map

    def lift_map(self):
        """
        Given this quotient space $Q = V/W$, return a fixed choice of linear homomorphism
        (a section) from Q to V.

        EXAMPLES:
            sage: M = QQ^3 / [[1,2,3]]
            sage: M.lift_map()
            Vector space morphism represented by the matrix:
            [1 0 0]
            [0 1 0]
            Domain: Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of dimension 3 over Rational Field
            W: Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
            Codomain: Vector space of dimension 3 over Rational Field
        """
        return self.__lift_map

    def lift(self, x):
        """
        Lift element of this quotient V/W to V by applying the fixed lift homomorphism.

        The lift is a fixed homomorphism.

        EXAMPLES:
            sage: M = QQ^3 / [[1,2,3]]
            sage: M.lift(M.0)
            (1, 0, 0)
            sage: M.lift(M.1)
            (0, 1, 0)
            sage: M.lift(M.0 - 2*M.1)
            (1, -2, 0)
        """
        return self.__lift_map(x)

    def W(self):
        """
        Given this quotient space $Q = V/W$, return W.

        EXAMPLES:
            sage: M = QQ^10 / [range(10), range(2,12)]
            sage: M.W()
            Vector space of degree 10 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1 -2 -3 -4 -5 -6 -7 -8]
            [ 0  1  2  3  4  5  6  7  8  9]
        """
        return self.__sub

    def V(self):
        """
        Given this quotient space $Q = V/W$, return $V$.

        EXAMPLES:
            sage: M = QQ^10 / [range(10), range(2,12)]
            sage: M.V()
            Vector space of dimension 10 over Rational Field
        """
        return self.__domain

    def cover(self):
        """
        Given this quotient space $Q = V/W$, return $V$.  This is the same as self.V().

        EXAMPLES:
            sage: M = QQ^10 / [range(10), range(2,12)]
            sage: M.cover()
            Vector space of dimension 10 over Rational Field
        """
        return self.V()

    def relations(self):
        """
        Given this quotient space $Q = V/W$, return $W$.  This is the same as self.W().

        EXAMPLES:
            sage: M = QQ^10 / [range(10), range(2,12)]
            sage: M.relations()
            Vector space of degree 10 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1 -2 -3 -4 -5 -6 -7 -8]
            [ 0  1  2  3  4  5  6  7  8  9]
        """
        return self.W()
