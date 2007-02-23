r"""
Free modules

SAGE supports computation with free modules over an arbitrary
commutative ring.  Nontrivial functionality is available over $\Z$ and
fields.  All free modules over an integral domain are equipped with an
embedding in an ambient vector space and an inner product, which you
can specify and change.

Create the free module of rank $n$ over an arbitrary commutative
ring $R$ using the command \code{FreeModule(R,n)}.  Equivalently,
\code{R\^n} also creates that free module.

The following example illustrates the creation of both a vector spaces
and a free module over the integers and a submodule of it.  Use the functions
\code{FreeModule}, \code{span} and member functions of free modules
to create free modules.  \emph{Do not use the FreeModule\_xxx constructors
directly.}

EXAMPLES:
    sage: V = VectorSpace(QQ,3)
    sage: W = V.subspace([[1,2,7], [1,1,0]])
    sage: W
    Vector space of degree 3 and dimension 2 over Rational Field
    Basis matrix:
    [ 1  0 -7]
    [ 0  1  7]
    sage: C = VectorSpaces(FiniteField(7))
    sage: C
    Category of vector spaces over Finite Field of size 7
    sage: C(W)
    Vector space of degree 3 and dimension 2 over Finite Field of size 7
    Basis matrix:
    [1 0 0]
    [0 1 0]

    sage: M = ZZ^3
    sage: C = VectorSpaces(FiniteField(7))
    sage: C(M)
    Vector space of dimension 3 over Finite Field of size 7
    sage: W = M.submodule([[1,2,7], [8,8,0]])
    sage: C(W)
    Vector space of degree 3 and dimension 2 over Finite Field of size 7
    Basis matrix:
    [1 0 0]
    [0 1 0]

We illustrate the exponent notation for creation of free modules.
    sage: ZZ^4
    Ambient free module of rank 4 over the principal ideal domain Integer Ring
    sage: QQ^2
    Vector space of dimension 2 over Rational Field
    sage: RR^3
    Vector space of dimension 3 over Real Field with 53 bits of precision

Base ring:
    sage: R.<x,y> = MPolynomialRing(QQ,2)
    sage: M = FreeModule(R,2)
    sage: M.base_ring()
    Polynomial Ring in x, y over Rational Field

    sage: VectorSpace(QQ, 10).base_ring()
    Rational Field
"""

####################################################################################
#       Copyright (C) 2005, 2007 William Stein <wstein@gmail.com>
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

# python imports
import random
import weakref


# SAGE imports
import free_module_element

import module

import sage.matrix.matrix_space

import sage.misc.latex as latex

import sage.rings.commutative_ring as commutative_ring
import sage.rings.principal_ideal_domain as principal_ideal_domain
import sage.rings.field as field
import sage.rings.integral_domain as integral_domain
import sage.rings.ring as ring
import sage.rings.integer_ring
import sage.rings.infinity
import sage.rings.integer
import sage.structure.parent_gens as gens
import sage.modules.real_double_vector
import sage.modules.complex_double_vector
from sage.structure.sequence import Sequence

from sage.structure.parent_gens import ParentWithGens

###############################################################################
#
# Constructor functions
#
###############################################################################
_cache = {}
def FreeModule(base_ring, rank, sparse=False, inner_product_matrix=None):
    r"""
    Create the free module over the given commutative ring of the given rank.

    INPUT:
        base_ring -- a commutative ring
        rank -- a nonnegative integer
        sparse -- bool; (default False)
        inner_product_matrix -- the inner product matrix (default None)

    OUTPUT:
        a free module

    \note{In \sage it is {\em not} the case that there is only one
    free ambient module of rank $n$ over $R$.  If you create $R^n$ twice
    \sage creates two separate objects.  This is because one can
    change the inner product on an ambient free module at any time.}

    EXAMPLES:
    First we illustrate creating free modules over various base
    fields.  The base field affects the free module that is created.
    For example, free modules over a field are vector spaces, and free
    modules over a principal ideal domain are special in that more
    functionality is available for them than for completely general
    free modules.

        sage: FreeModule(Integers(8),10)
        Ambient free module of rank 10 over Ring of integers modulo 8
        sage: FreeModule(QQ,10)
        Vector space of dimension 10 over Rational Field
        sage: FreeModule(ZZ,10)
        Ambient free module of rank 10 over the principal ideal domain Integer Ring
        sage: FreeModule(FiniteField(5),10)
        Vector space of dimension 10 over Finite Field of size 5
        sage: FreeModule(Integers(7),10)
        Vector space of dimension 10 over Ring of integers modulo 7
        sage: FreeModule(PolynomialRing(QQ,'x'),5)
        Ambient free module of rank 5 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
        sage: FreeModule(PolynomialRing(ZZ,'x'),5)
        Ambient free module of rank 5 over the integral domain Univariate Polynomial Ring in x over Integer Ring

    Of course we can make rank 0 free modules:

        sage: FreeModule(RealField(100),0)
        Vector space of dimension 0 over Real Field with 100 bits of precision

    Next we create a free module with sparse representation of
    elements.  Functionality with sparse modules is \emph{identical} to
    dense modules, but they may use less memory and arithmetic may
    be faster (or slower!).

        sage: M = FreeModule(ZZ,200,sparse=True)
        sage: M.is_sparse()
        True
        sage: type(M.0)
        <type 'sage.modules.free_module_element.FreeModuleElement_generic_sparse'>

    The default is dense.
        sage: M = ZZ^200
        sage: type(M.0)
        <type 'sage.modules.vector_integer_dense.Vector_integer_dense'>


    Note that matrices associated in some way to sparse free modules
    are sparse by default:
        sage: M = FreeModule(Integers(8), 2)
        sage: A = M.basis_matrix()
        sage: A.is_sparse()
        False
        sage: Ms = FreeModule(Integers(8), 2, sparse=True)
        sage: M == Ms  # as mathematical objects they are equal
        True
        sage: Ms.basis_matrix().is_sparse()
        True


    We can also specify an inner product matrix, which is used when
    computing inner products of elements.
        sage: A = MatrixSpace(ZZ,2)([[1,0],[0,-1]])
        sage: M = FreeModule(ZZ,2,inner_product_matrix=A)
        sage: v, w = M.gens()
        sage: v.inner_product(w)
        0
        sage: v.inner_product(v)
        1
        sage: w.inner_product(w)
        -1
        sage: (v+2*w).inner_product(w)
        -2

    You can also specify the inner product matrix by giving anything that
    coerces to an appropriate matrix.   This is only useful if the inner
    product matrix takes values in the base ring.
        sage: FreeModule(ZZ,2,inner_product_matrix=1).inner_product_matrix()
        [1 0]
        [0 1]
        sage: FreeModule(ZZ,2,inner_product_matrix=[1,2,3,4]).inner_product_matrix()
        [1 2]
        [3 4]
        sage: FreeModule(ZZ,2,inner_product_matrix=[[1,2],[3,4]]).inner_product_matrix()
        [1 2]
        [3 4]
    """
    global _cache
    rank = int(rank)
    if not (inner_product_matrix is None):
        inner_product_matrix = sage.matrix.matrix_space.MatrixSpace(base_ring, rank)(inner_product_matrix)
        inner_product_matrix.set_immutable()

    key = (base_ring, rank, sparse, inner_product_matrix)
    if _cache.has_key(key):
        M = _cache[key]()
        if not (M is None):
            return M

    if not base_ring.is_commutative():
        raise TypeError, "base_ring must be a commutative ring"

    if not sparse and isinstance(base_ring,sage.rings.real_double.RealDoubleField_class):
        M = RealDoubleVectorSpace_class(rank)

    elif not sparse and isinstance(base_ring,sage.rings.complex_double.ComplexDoubleField_class):
        M=ComplexDoubleVectorSpace_class(rank)

    elif base_ring.is_field():
        M = FreeModule_ambient_field(base_ring, rank,
                                        sparse=sparse, inner_product_matrix=inner_product_matrix)

    elif isinstance(base_ring, principal_ideal_domain.PrincipalIdealDomain):
        M = FreeModule_ambient_pid(base_ring, rank,
                                   sparse=sparse, inner_product_matrix=inner_product_matrix)

    elif isinstance(base_ring, integral_domain.IntegralDomain):
        M = FreeModule_ambient_domain(base_ring, rank,
                                         sparse=sparse, inner_product_matrix=inner_product_matrix)


    else:
        M = FreeModule_ambient(base_ring, rank,
                                  sparse=sparse, inner_product_matrix=inner_product_matrix)

    _cache[key] = weakref.ref(M)
    return M


def VectorSpace(K, dimension, sparse=False,  inner_product_matrix=None):
    """
    EXAMPLES:
    The base can be complicated, as long as it is a field.
        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x over Integer Ring
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1)
        ]

    The base must be a field or a \code{TypeError} is raised.
        sage: VectorSpace(ZZ,5)
        Traceback (most recent call last):
        ...
        TypeError: K must be a field
    """

    if not field.is_Field(K):
        raise TypeError, "K must be a field"
    return FreeModule(K, dimension, sparse, inner_product_matrix)

def span(R, gens, check=True, already_echelonized=False):
    """
    Return the $R$-span of gens.

    EXAMPLES:
        sage: V = span(QQ, [[1,2,5], [2,2,2]]); V
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -3]
        [ 0  1  4]
        sage: span(QuadraticField(-7,'a'), [V.gen(0)])
        Vector space of degree 3 and dimension 1 over Number Field in a with defining polynomial x^2 + 7
        Basis matrix:
        [ 1  0 -3]
        sage: span(GF(2), [[1,2,3], [2,2,2], [1,2,5]])
        Vector space of degree 3 and dimension 1 over Finite Field of size 2
        Basis matrix:
        [1 0 1]
    """
    if len(gens) == 0:
        return FreeModule(R, 0)
    else:
        x = gens[0]
        M = FreeModule(R, len(x))
        return M.span(gens, check=check, already_echelonized=already_echelonized)


###############################################################################
#
# Base class for all free modules
#
###############################################################################

def is_FreeModule(x):
    return isinstance(x, FreeModule_generic)

class FreeModule_generic(module.Module):
    """
    Base class for all free modules.
    """
    def __init__(self, base_ring, rank, degree, sparse=False, inner_product_matrix=None):
        """
        Create the free module of given rank over the given base_ring.

        INPUT:
            base_ring -- a commutative ring
            rank -- a non-negative integer
        """
        if not isinstance(base_ring, commutative_ring.CommutativeRing):
            raise TypeError, "base_ring (=%s) must be a commutative ring"%base_ring
        rank = sage.rings.integer.Integer(rank)
        if rank < 0:
            raise ValueError, "rank (=%s) must be nonnegative"%rank
        degree = sage.rings.integer.Integer(degree)
        if degree < 0:
            raise ValueError, "degree (=%s) must be nonnegative"%degree

        ParentWithGens.__init__(self, base_ring)     # names aren't used anywhere.
        self.__uses_ambient_inner_product = True
        self.__rank = rank
        self.__degree = degree
        self.__is_sparse = sparse
        self._inner_product_matrix = inner_product_matrix
        self.element_class()

    def _an_element_impl(self):
        return self.zero_vector()

    def element_class(self):
        try:
            return self._element_class
        except AttributeError:
            pass
        C = element_class(self.base_ring(), self.is_sparse())
        self._element_class = C
        return C

    def __call__(self, x, coerce=True, copy=True, check=True):
        if isinstance(x, (int, long, sage.rings.integer.Integer)) and x==0:
            return self.zero_vector()
        elif isinstance(x, free_module_element.FreeModuleElement):
            if x.parent() is self:
                if copy:
                    return x.__copy__()
                else:
                    return x
            x = x.list()
        w = self._element_class(self, x, coerce, copy)
        if check:
            if isinstance(self, FreeModule_ambient):
                return w
            self.coordinates(w)
        return w

    def _has_coerce_map_from_space(self, V):
        """
        Return True if V canonically coerces to self.
        """
        try:
            return self.__has_coerce_map_from_space[V]
        except AttributeError:
            self.__has_coerce_map_from_space = {}
        except KeyError:
            pass
        if self.base_ring() is V.base_ring():
            h = V.is_submodule(self)
        elif not self.base_ring().has_coerce_map_from(V.base_ring()):
            self.__has_coerce_map_from_space[V] = False
            return False
        else:
            h = V.base_extend(self.base_ring()).is_submodule(self)
        self.__has_coerce_map_from_space[V] = h
        return h

    def _coerce_impl(self, x):
        """
        Canonical coercion of x into this free module.
            (0, 4/3, 8/3, 4, 16/3)
        """
        if isinstance(x, (int, long, sage.rings.integer.Integer)) and x==0:
            return self.zero_vector()
        if isinstance(x, free_module_element.FreeModuleElement):
            # determining if the map exists is expensive the first time,
            # so we cache it.
            if self._has_coerce_map_from_space(x.parent()):
                return self(x)
        raise TypeError

    def __cmp__(self, right):
        raise NotImplementedError


    def __contains__(self, v):
        r"""
        EXAMPLES:

        We create the module $\Z^3$, and the submodule generated by
        one vector $(1,1,0)$, and check whether certain elements are
        in the submodule.

            sage: R = FreeModule(ZZ, 3)
            sage: V = R.submodule([R.gen(0) + R.gen(1)])
            sage: R.gen(0) + R.gen(1) in V
            True
            sage: R.gen(0) + 2*R.gen(1) in V
            False

            sage: Q = QQ
            sage: w = Q('1/2')*(R.gen(0) + R.gen(1))
            sage: w
            (1/2, 1/2, 0)
            sage: w.parent()
            Vector space of dimension 3 over Rational Field
            sage: w in V
            False
            sage: V.coordinates(w)
            [1/2]

        """
        if not isinstance(v, free_module_element.FreeModuleElement):
            return False
        if v.parent() is self:
            return True
        try:
            c = self.coordinates(v)
        except (ArithmeticError, TypeError):
            return False
        # Finally, check that each coordinate lies in the base ring.
        R = self.base_ring()
        if not self.base_ring().is_field():
            for a in c:
                try:
                    b = R(a)
                except TypeError:
                    return False
        return True

    def __iter__(self):
        """
        Return iterator over the elements of this free module.

        EXAMPLES:
            sage: V = VectorSpace(GF(4,'a'),2)
            sage: [x for x in V]
            [(0, 0), (a, 0), (a + 1, 0), (1, 0), (0, a), (a, a), (a + 1, a), (1, a), (0, a + 1), (a, a + 1), (a + 1, a + 1), (1, a + 1), (0, 1), (a, 1), (a + 1, 1), (1, 1)]

            sage: W = V.subspace([V([1,1])])
            sage: print [x for x in W]
            [(0, 0), (a, a), (a + 1, a + 1), (1, 1)]
        """
        G = self.gens()
        if len(G) == 0:
            yield self(0)
            return
        R     = self.base_ring()
        iters = [iter(R) for _ in range(len(G))]
        for x in iters: x.next()     # put at 0
        zero  = R(0)
        v = [zero for _ in range(len(G))]
        n = 0
        z = self(0)
        yield z
        while n < len(G):
            try:
                v[n] = iters[n].next()
                yield self.linear_combination_of_basis(v)
                n = 0
            except StopIteration:
                iters[n] = iter(R)  # reset
                iters[n].next()     # put at 0
                v[n] = zero
                n += 1


    def __len__(self):
        return len(self.base_ring())**self.rank()

    def ambient_module(self):
        """
        Return the ambient module associated to this module.

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: M = FreeModule(R,2)
            sage: M.ambient_module()
            Ambient free module of rank 2 over the integral domain Polynomial Ring in x, y over Rational Field

            sage: V = FreeModule(QQ, 4).span([[1,2,3,4], [1,0,0,0]]); V
            Vector space of degree 4 and dimension 2 over Rational Field
            Basis matrix:
            [  1   0   0   0]
            [  0   1 3/2   2]
            sage: V.ambient_module()
            Vector space of dimension 4 over Rational Field
        """
        return FreeModule(self.base_ring(), self.degree())

    def base_extend(self, R):
        r"""
        Return the base extension of self to R.  This
        is the same as \code{self.change_ring(R)} except
        that a TypeError is raised if there is no canonical
        coerce map from the base ring of self to R.

        INPUT:
            R -- ring
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        raise TypeError, "base extension of self (over '%s') to ring '%s' not defined."%(self.base_ring(),R)

    def basis(self):
        """
        Return the basis of this module.

        EXAMPLES:
            sage: FreeModule(Integers(12),3).basis()
            [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
        """
        raise NotImplementedError

    def basis_matrix(self):
        """
        Return the matrix whose rows are the basis for this free module.

        EXAMPLES:
            sage: FreeModule(Integers(12),3).basis_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: M = FreeModule(GF(7),3).span([[2,3,4],[1,1,1]]); M
            Vector space of degree 3 and dimension 2 over Finite Field of size 7
            Basis matrix:
            [1 0 6]
            [0 1 2]
            sage: M.basis_matrix()
            [1 0 6]
            [0 1 2]

            sage: M = FreeModule(GF(7),3).span_of_basis([[2,3,4],[1,1,1]]);
            sage: M.basis_matrix()
            [2 3 4]
            [1 1 1]
        """
        try:
            return self.__basis_matrix
        except AttributeError:
            MAT = sage.matrix.matrix_space.MatrixSpace(self.base_ring(),
                            len(self.basis()), self.degree(),
                            sparse = self.is_sparse())
            self.__basis_matrix = MAT(self.basis())
            return self.__basis_matrix

    def echelonized_basis_matrix(self):
        raise NotImplementedError

    def category(self):
        """
        Return the category to which this free module belongs.  This is the
        category of all free modules over the base ring.

        EXAMPLES:
            sage: FreeModule(GF(7),3).category()
            Category of vector spaces over Finite Field of size 7
        """
        import sage.categories.all
        return sage.categories.all.FreeModules(self.base_ring())

    def matrix(self):
        """
        Return the basis matrix of this module, which is the matrix
        whose rows are a basis for this module.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 2)
            sage: M.matrix()
            [1 0]
            [0 1]
            sage: M.submodule([M.gen(0) + M.gen(1), M.gen(0) - 2*M.gen(1)]).matrix()
            [1 1]
            [0 3]
        """
        return self.basis_matrix()

    def coordinates(self, v, check=True):
        """
        Write $v$ in terms of the basis for self.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list

        Returns a list $c$ such that if $B$ is the basis for self, then
        $$
                \sum c_i B_i = v.
        $$
        If $v$ is not in self, raises an \code{ArithmeticError} exception.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0-M1)
            [2, -1]
        """
        return self.coordinate_vector(v, check=check).list()

    def coordinate_vector(self, v, check=True):
        """
        Return the a vector whose cofficients give $v$ as a linear combination
        of the basis for self.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list

        EXAMPLES:
            sage: M = FreeModule(ZZ, 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinate_vector(2*M0 - M1)
            (2, -1)
        """
        raise NotImplementedError

    def degree(self):
        """
        Return the degree of this free module.  This is the dimension of the
        ambient vector space in which it is embedded.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 10)
            sage: W = M.submodule([M.gen(0), 2*M.gen(3) - M.gen(0), M.gen(0) + M.gen(3)])
            sage: W.degree()
            10
            sage: W.rank()
            2
        """
        return self.__degree

    def dimension(self):
        """
        Return the dimension of this free module.

        EXAMPLES:
            sage: M = FreeModule(FiniteField(19), 100)
            sage: W = M.submodule([M.gen(50)])
            sage: W.dimension()
            1
        """
        return self.rank()

    def discriminant(self):
        """
        Return the discriminant of this free module.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 3)
            sage: M.discriminant()
            1
            sage: W = M.span([[1,2,3]])
            sage: W.discriminant()
            14
            sage: W2 = M.span([[1,2,3], [1,1,1]])
            sage: W2.discriminant()
            6
        """
        return self.inner_product_matrix().determinant()

    def free_module(self):
        """
        Return this free module.  (This is used by the \code{FreeModule} functor,
        and simply returns self.)

        EXAMPLES:
            sage: M = FreeModule(ZZ, 3)
            sage: M.free_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        return self

    def gen(self, i=0):
        if i < 0 or i >= self.rank():
            raise ValueError, "Generator %s not defined."%i
        return self.basis()[int(i)]

    def has_user_basis(self):
        """
        Return \code{True} if the basis of this free module is specified by the user,
        as opposed to being the default echelon form.

        EXAMPLES:
            sage: V = QQ^3
            sage: W = V.subspace([[2,'1/2', 1]])
            sage: W.has_user_basis()
            False
            sage: W = V.subspace_with_basis([[2,'1/2',1]])
            sage: W.has_user_basis()
            True
        """
        return False

    def inner_product_matrix(self):
        """
        Return the inner product matrix associated to this free module.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 3)
            sage: M.inner_product_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: FreeModule(ZZ,2,inner_product_matrix=[[1,-1],[2,5]]).inner_product_matrix()
            [ 1 -1]
            [ 2  5]
        """
        if self._inner_product_matrix is None:
            # Not set, so induced by the ambient space or the identity matrix
            # if this space is ambient.
            if isinstance(self, FreeModule_ambient):
                M = sage.matrix.matrix_space.MatrixSpace(
                     self.base_ring(), self.rank(), sparse=True)(1)
            else:
                # TODO: probably much better way to do this using matrix multiplication!
                A = self.ambient_vector_space()
                I = [A(v).inner_product(A(w)) for v in self.basis() for w in self.basis()]
                M = sage.matrix.matrix_space.MatrixSpace( \
                      self.base_field(), self.rank(), sparse=False)(I)

            self._inner_product_matrix = M
        return self._inner_product_matrix

    def _inner_product_is_dot_product(self):
        """
        Return whether or not the inner product on this module is induced by
        the dot product on the ambient vector space.  This is used internally
        by the inner_product function for optimization.

        EXAMPLES:
            sage: FreeModule(ZZ, 3)._inner_product_is_dot_product()
            True
            sage: FreeModule(ZZ, 3, inner_product_matrix=1)._inner_product_is_dot_product()
            True
            sage: FreeModule(ZZ, 2, inner_product_matrix=[1,0,-1,0])._inner_product_is_dot_product()
            False

            sage: M = FreeModule(QQ, 3)
            sage: M2 = M.span([[1,2,3]])
            sage: M2._inner_product_is_dot_product()
            True
        """
        try:
            return self._inner_product_is_dot_product_cache
        except AttributeError:
            if self.is_ambient():
                if self._inner_product_matrix is None or self._inner_product_matrix == sage.rings.integer.Integer(1):
                    self._inner_product_is_dot_product_cache = True
                else:
                    self._inner_product_is_dot_product_cache = False
            else:
                self._inner_product_is_dot_product_cache = self.ambient_module()._inner_product_is_dot_product()
            return self._inner_product_is_dot_product_cache

# I commented out both setting and unsetting the inner product matrix
# since being able to change them contradicts that unique objects are
# immutable.
##     def unset_inner_product_matrix(self):
##         r"""
##         If an inner product was set on this module using
##         \code{self.set_inner_product_matrix(...)}, this function
##         unsets that inner product, thus reverting to the inner
##         product induced from that on the ambient module.

##         EXAMPLES:
##             sage: M = FreeModule(QQ, 2)
##             sage: M.set_inner_product_matrix([-1,0,0,-1])
##             sage: (M.0).inner_product(M.0)
##             -1

##         We set and unset an inner product matrix on a submodule.  Note
##         that unsetting the inner product on the submodule switches back
##         to using the ambient inner product.
##             sage: M2 = M.submodule([[1,1]])
##             sage: M2.gen(0).inner_product(M2.gen(0))
##             -2
##             sage: M2.set_inner_product_matrix([1])
##             sage: M2.gen(0).inner_product(M2.gen(0))
##             1
##             sage: M2.unset_inner_product_matrix()
##             sage: M2.gen(0).inner_product(M2.gen(0))
##             -2

##         Changing the ambient inner product changes that
##         on the submodule, since the submodule uses the ambient
##         product after calling \code{unset_inner_product_matrix}.
##             sage: M.unset_inner_product_matrix()
##             sage: M2.gen(0).inner_product(M2.gen(0))
##             2
##         """
##         self.__uses_ambient_inner_product = True
##         try:
##             del self._inner_product_is_dot_product_cache
##             self._inner_product_matrix = None
##         except AttributeError:
##             pass

##     def set_inner_product_matrix(self, A):
##         """
##         Sets the inner product matrix of this module to the matrix A.

##         EXAMPLES:
##         We change the inner product matrix over an ambient free module.
##             sage: M = FreeModule(ZZ, 2)
##             sage: M.inner_product_matrix()
##             [1 0]
##             [0 1]
##             sage: M.set_inner_product_matrix([0,-1,-1,0])
##             sage: M.inner_product_matrix()
##             [ 0 -1]
##             [-1  0]
##             sage: (M.0).inner_product(M.1)
##             -1

##         We can also set the inner product matrix of a submodule:
##             sage: W = M.submodule([[1,2]])
##             sage: W.set_inner_product_matrix([2])
##             sage: W
##             Free module of degree 2 and rank 1 over Integer Ring
##             Echelon basis matrix:
##             [1 2]
##             sage: W.inner_product_matrix()
##             [2]
##             sage: v = W.gen(0)
##             sage: v.inner_product(v)
##             2
##         """
##         if not isinstance(A, sage.matrix.matrix_generic.Matrix):
##             A = sage.matrix.matrix_space.MatrixSpace(self.base_ring(), self.rank())(A)
##         else:
##             if A.nrows() != self.rank():
##                 raise ArithmeticError, "A must have %s rows but it has %s rows"%(self.rank(),
##                                                                                 A.nrows())
##             if A.ncols() != self.rank():
##                 raise ArithmeticError, "A must have %s columns but it has %s rows"%(self.rank(),
##                                                                                 A.ncols())
##         self._inner_product_matrix = A
##         self.__uses_ambient_inner_product = False
##         try:
##             del self._inner_product_is_dot_product_cache
##         except AttributeError:
##             pass

    def is_ambient(self):
        """
        Returns False sense this is not an ambient free module.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 3).span([[1,2,3]]); M
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            sage: M.is_ambient()
            False
            sage: M = (ZZ^2).span([[1,0], [0,1]])
            sage: M
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 1]
            sage: M.is_ambient()
            False
            sage: M == M.ambient_module()
            True
        """
        return False

    def is_dense(self):
        """
        Return \code{True} if the underlying representation of this module uses dense vectors,
        and False otherwise.

        EXAMPLES:
            sage: FreeModule(ZZ, 2).is_dense()
            True
            sage: FreeModule(ZZ, 2, sparse=True).is_dense()
            False
        """
        return not self.__is_sparse

    def is_full(self):
        """
        Return \code{True} if the rank of this module equals its degree.

        EXAMPLES:
            sage: FreeModule(ZZ, 2).is_full()
            True
            sage: M = FreeModule(ZZ, 2).span([[1,2]])
            sage: M.is_full()
            False
        """
        return self.rank() == self.degree()

    def is_finite(self):
        """
        Returns True if the underlying set of this free module is finite.

        EXAMPLES:
            sage: FreeModule(ZZ, 2).is_finite()
            False
            sage: FreeModule(Integers(8), 2).is_finite()
            True
            sage: FreeModule(ZZ, 0).is_finite()
            True
        """
        return self.base_ring().is_finite() or self.rank() == 0

    def is_sparse(self):
        """
        Return \code{True} if the underlying representation of this module uses sparse vectors,
        and False otherwise.

        EXAMPLES:
            sage: FreeModule(ZZ, 2).is_sparse()
            False
            sage: FreeModule(ZZ, 2, sparse=True).is_sparse()
            True
        """
        return self.__is_sparse

    def ngens(self):
        """
        Returns the number of basis elements of this free module.

        EXAMPLES:
            sage: FreeModule(ZZ, 2).ngens()
            2
            sage: FreeModule(ZZ, 0).ngens()
            0
            sage: FreeModule(ZZ, 2).span([[1,1]]).ngens()
            1
        """
        try:
            return self.__ngens
        except AttributeError:
            self.__ngens = self.rank()
        return self.__ngens

    def nonembedded_free_module(self):
        """
        Returns an ambient free module that is isomorphic to this free module.

        Thus if this free module is of rank $n$ over a ring $R$, then this function
        returns $R^n$, as an ambient free module.

        EXAMPLES:
            sage: FreeModule(ZZ, 2).span([[1,1]]).nonembedded_free_module()
            Ambient free module of rank 1 over the principal ideal domain Integer Ring
        """
        return FreeModule(self.base_ring(), self.rank())

    def random_element(self, bound=2, prob=1.0):
        """
        Returns a random element of self.

        INPUT:
            bound -- integer; coefficients of linear combination of basis are chosen
                     from the closed interval [-bound,bound]
            prob -- float; probability that given coefficient is nonzero.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 2).span([[1,1]])
            sage: x = M.random_element()
            sage.: x
            (1, 1)
            sage.: M.random_element()
            (-2, -2)
            sage.: M.random_element()
            (-1, -1)
        """
        R = self.base_ring()
        v = self(0)
        prob = float(prob)
        for i in range(self.rank()):
            if random.random() <= prob:
                v += self.gen(i) * R.random_element(bound)
        return v

    def rank(self):
        """
        Return the rank of this free module.

        EXAMPLES:
            sage: FreeModule(Integers(6), 10000000).rank()
            10000000
            sage: FreeModule(ZZ, 2).span([[1,1], [2,2], [3,4]]).rank()
            2
        """
        return self.__rank

##     def uses_ambient_inner_product(self):
##         """
##         Return \code{True} if the inner product on this module is the one induced
##         by the ambient inner product.  This is True exactly if
##         self.set_inner_product_matrix(...) has not been called (or
##         self.unset_inner_product_matrix() was subsequently called).

##         EXAMPLES:
##             sage: M = FreeModule(ZZ, 2)
##             sage: W = M.submodule([[1,2]])
##             sage: W.uses_ambient_inner_product()
##             True
##             sage: W.inner_product_matrix()
##             [5]
##             sage: W = FreeModule(ZZ, 2, inner_product_matrix = [2])
##             sage: W.uses_ambient_inner_product()
##             False
##             sage: W.inner_product_matrix()
##             [2]
##         """
##         return self.__uses_ambient_inner_product

    def zero_vector(self):
        """
        Returns the zero vector in this free module.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 2)
            sage: M.zero_vector()
            (0, 0)
            sage: M(0)
            (0, 0)
            sage: M.span([[1,1]]).zero_vector()
            (0, 0)
            sage: M.zero_submodule().zero_vector()
            (0, 0)
        """
        # Do *not* cache this -- it must be computed fresh each time, since
        # it is is used by __call__ to make a new copy of the 0 element.

        return self._element_class(self, 0)

class FreeModule_generic_pid(FreeModule_generic):
    """
    Base class for all free modules over a PID.
    """
    def __init__(self, base_ring, rank, degree, sparse=False, inner_product_matrix=None):
        """
        Create a free module over a PID.

        EXAMPLES:
            sage: FreeModule(ZZ, 2)
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: FreeModule(PolynomialRing(GF(7),'x'), 2)
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Finite Field of size 7
        """
        if not isinstance(base_ring, principal_ideal_domain.PrincipalIdealDomain):
            raise TypeError, "base_ring must be a principal ideal domain"
        FreeModule_generic.__init__(self, base_ring, rank, degree, sparse, inner_product_matrix)

    def scale(self, other):
        """
        Return the product of this module by the number other, which
        is the module spanned by other times each basis vector.

        EXAMPLES:
            sage: M = FreeModule(ZZ, 3)
            sage: M.scale(2)
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [2 0 0]
            [0 2 0]
            [0 0 2]

            sage: a = QQ('1/3')
            sage: M.scale(a)
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1/3   0   0]
            [  0 1/3   0]
            [  0   0 1/3]
        """
        if other == 0:
            return self.zero_submodule()
        if other == 1 or other == -1:
            return self
        return self.span([v*other for v in self.basis()])

    def __add__(self, other):
        r"""
        Return the sum of self and other, where both self and other
        must be submodules of the ambient vector space.

        EXAMPLES:
        We add two vector spaces:

            sage: V  = VectorSpace(QQ, 3)
            sage: W  = V.subspace([V([1,1,0])])
            sage: W2 = V.subspace([V([1,-1,0])])
            sage: W + W2
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]

        We add two free $\Z$-modules.
            sage: M = FreeModule(ZZ, 3)
            sage: W = M.submodule([M([1,0,2])])
            sage: W2 = M.submodule([M([2,0,-4])])
            sage: W + W2
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0 2]
            [0 0 8]

        We can also add free $\Z$-modules embedded non-integrally
        into an ambient space.

            sage: V = VectorSpace(QQ, 3)
            sage: W = M.span([1/2*V.0 - 1/3*V.1])

        Here the command \code{M.span(...)} creates the span of the
        indicated vectors over the base ring of $M$.

            sage: W2 = M.span([1/3*V.0 + V.1])
            sage: W + W2
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1/6  7/3    0]
            [   0 11/3    0]

        We add two modules over $\Z$:
            sage: A = Matrix(ZZ, 3, 3, [3, 0, -1, 0, -2, 0, 0, 0, -2])
            sage: V = (A+2).kernel()
            sage: W = (A-3).kernel()
            sage: V+W
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [5 0 0]
            [0 1 0]
            [0 0 1]
        """
        if not isinstance(other, FreeModule_generic):
            raise TypeError, "other (=%s) must be a free module"%other
        if not (self.ambient_vector_space() == other.ambient_vector_space()):
            raise TypeError, "ambient vector spaces must be equal"
        return self.span(self.basis() + other.basis())


    def base_field(self):
        """
        Return the base field, which is the fraction field of the base
        ring of this module.

        EXAMPLES:
            sage: FreeModule(GF(3), 2).base_field()
            Finite Field of size 3
            sage: FreeModule(ZZ, 2).base_field()
            Rational Field
            sage: FreeModule(PolynomialRing(GF(7),'x'), 2).base_field()
            Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 7
        """
        return self.base_ring().fraction_field()

    def basis_matrix(self):
        """
        Return the matrix whose rows are the basis for this free module.

        EXAMPLES:
            sage: M = FreeModule(QQ,2).span_of_basis([[1,-1],[1,0]]); M
            Vector space of degree 2 and dimension 2 over Rational Field
            User basis matrix:
            [ 1 -1]
            [ 1  0]
            sage: M.basis_matrix()
            [ 1 -1]
            [ 1  0]
        """
        try:
            return self.__basis_matrix
        except AttributeError:
            MAT = sage.matrix.matrix_space.MatrixSpace(self.base_field(),
                            len(self.basis()), self.degree(),
                            sparse = self.is_sparse())
            self.__basis_matrix = MAT(self.basis())
            return self.__basis_matrix

    def index(self, other):
        """
        Return the lattice index [self : other] which is an element
        of the base field.  When other is contained in self, the lattice
        index is the usual index.  If the index is infinite, then this
        function returns infinity.

        EXAMPLES:
            sage: L1 = span(ZZ, [[1,2]])
            sage: L2 = span(ZZ, [[3,6]])
            sage: L1.index(L2)
            3

        Note that the free modules being compared need not be integral.
            sage: L1 = span(ZZ, [['1/2','1/3'], [4,5]])
            sage: L2 = span(ZZ, [[1,2], [3,4]])
            sage: L1.index(L2)
            12/7
            sage: L2.index(L1)
            7/12
            sage: L1.discriminant() / L2.discriminant()
            49/144

        The index of a lattice of infinite index is infinite.
            sage: L1 = FreeModule(ZZ, 2)
            sage: L2 = span(ZZ, [[1,2]])
            sage: L1.index(L2)
            Infinity
        """
        if not isinstance(other, FreeModule_generic):
            raise TypeError, "other must be a free module"

        if self.ambient_vector_space() != other.ambient_vector_space():
            raise ArithmeticError, "self and other must be embedded in the same ambient space."

        if self.base_ring() != other.base_ring():
            raise NotImplementedError, "lattice index only defined for modules over the same base ring."

        if self.base_ring().is_field():
            if self == other:
                return sage.rings.integer.Integer(1)
            else:
                if other.is_subspace(self):
                    return sage.rings.infinity.infinity
            raise ArithmeticError, "other must be contained in the vector space spanned by self."

        try:
            C = [self.coordinates(b) for b in other.basis()]
        except ArithmeticError:
            raise ArithmeticError, "other must be contained in the vector space spanned by self."

        if other.rank() < self.rank():
            return sage.rings.infinity.infinity

        A = sage.matrix.matrix_space.MatrixSpace(self.base_field(), self.rank())(C)
        return A.determinant()


    def intersection(self, other):
        r"""
        Return the intersection of self and other.

        EXAMPLES:
        We intersect two submodules one of which is clearly contained in
        the other.
            sage: A = ZZ^2
            sage: M1 = A.span([[1,1]])
            sage: M2 = A.span([[3,3]])
            sage: M1.intersection(M2)
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [3 3]

        We intersection two submodules of $\Z^3$ of rank $2$, whose intersection
        has rank $1$.
            sage: A = ZZ^3
            sage: M1 = A.span([[1,1,1], [1,2,3]])
            sage: M2 = A.span([[2,2,2], [1,0,0]])
            sage: M1.intersection(M2)
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [2 2 2]

        We compute an intersection of two $\Z$-modules that are not submodules
        of $\Z^2$.
            sage: A = ZZ^2
            sage: M1 = A.span([[1,2]]).scale(1/6)
            sage: M2 = A.span([[1,2]]).scale(1/15)
            sage: M1.intersection(M2)
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1/3 2/3]

        We intersect a $\Z$-module with a $\Q$-vector space.
            sage: A = ZZ^3
            sage: L = ZZ^3
            sage: V = QQ^3
            sage: W = L.span([[1/2,0,1/2]])
            sage: K = V.span([[1,0,1], [0,0,1]])
            sage: W.intersection(K)
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1/2   0 1/2]
            sage: K.intersection(W)
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1/2   0 1/2]
        """
        if not isinstance(other, FreeModule_generic):
            raise TypeError, "other must be a free module"

        if self.ambient_vector_space() != other.ambient_vector_space():
            raise ArithmeticError, "self and other must be embedded in the same ambient space."

        if self.base_ring() != other.base_ring():
            if other.base_ring().is_field():
                return other.intersection(self)
            raise NotImplementedError, "intersection of modules over different base rings (neither a field) is not implemented."

        # dispense with the three easy cases
        if self == self.ambient_vector_space():
            return other
        elif other == other.ambient_vector_space():
            return self
        elif self.rank() == 0 or other.rank() == 0:
            return self.zero_submodule()

        # standard algorithm for computing intersection of general submodule
        if self.dimension() <= other.dimension():
            V1 = self; V2 = other
        else:
            V1 = other; V2 = self
        A1 = V1.basis_matrix()
        A2 = V2.basis_matrix()
        S  = A1.stack(A2)
        K  = S.integer_kernel()
        n  = int(V1.dimension())
        B = [A1.linear_combination_of_rows(v.list()[:n]) for v in K.basis()]
        return self.span(B, check=False)

    def is_submodule(self, other):
        """
        True if this module is a submodule of other.

        EXAMPLES:
            sage: M = FreeModule(ZZ,2)
            sage: M.is_submodule(M)
            True
            sage: N = M.scale(2)
            sage: N.is_submodule(M)
            True
            sage: M.is_submodule(N)
            False
            sage: N = M.scale(1/2)
            sage: N.is_submodule(M)
            False
            sage: M.is_submodule(N)
            True
        """
        if not isinstance(other, FreeModule_generic):
            return False
        if self.ambient_vector_space() != other.ambient_vector_space():
            return False
        if other == other.ambient_vector_space():
            return True
        if other.rank() < self.rank():
            return False
        if self.base_ring() != other.base_ring():
            try:
                if not self.base_ring().is_subring(other.base_ring()):
                    return False
            except NotImplementedError:
                return False
        for b in self.basis():
            if not (b in other):
                return False
        return True

    def zero_submodule(self):
        """
        Return the zero submodule of this module.

        EXAMPLES:
            sage: V = FreeModule(ZZ,2)
            sage: V.zero_submodule()
            Free module of degree 2 and rank 0 over Integer Ring
            Echelon basis matrix:
            []
        """
        return FreeModule_submodule_pid(self.ambient_module(),[], check=False)

    def index_in_saturation(self):
        r"""
        Return the index of this module in its saturation, i.e., its
        intersection with $R^n$.

        EXAMPLES:
            sage: W = span(ZZ,[[2,4,6]])
            sage: W.index_in_saturation()
            2
            sage: W = span(ZZ,[[1/2,1/3]])
            sage: W.index_in_saturation()
            1/6
        """
        # TODO: There is probably a much faster algorithm in this case.
        return self.saturation().index(self)

    def saturation(self):
        r"""
        Return the saturated submodule of $R^n$ that spans the same
        vector space as self.

        EXAMPLES:
        We create a 1-dimensional lattice that is obviously not
        saturated and saturate it.
            sage: L = span(ZZ, [[9,9,6]]); L
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [9 9 6]
            sage: L.saturation()
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [3 3 2]

        We create a lattice spanned by two vectors, and saturate.
        Comptuation of discriminants shows that the index of lattice
        in its saturation is $3$, which is a prime of congruence between
        the two generating vectors.
            sage: L = span(ZZ, [[1,2,3], [4,5,6]])
            sage: L.saturation()
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: L.discriminant()
            54
            sage: L.saturation().discriminant()
            6

        Notice that the saturation of a non-integral lattice $L$ is defined,
        but the result is integral hence does not contain $L$:
            sage: L = span(ZZ, [['1/2',1,3]])
            sage: L.saturation()
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1 2 6]

        """
        R = self.base_ring()
        if R.is_field():
            return self
        V = self.vector_space()
        A = self.ambient_module()
        return V.intersection(A)

    def span(self, gens, check=True, already_echelonized=False):
        """
        Return the R-span of the given list of gens, where R
        is the base ring of self.  Note that this span need not
        be a submodule of self, nor even of the ambient space.
        It must, however, be contained in the ambient vector space, i.e.,
        the ambient space tensored with the fraction field of R.

        EXAMPLES:
            sage: V = FreeModule(ZZ,3)
            sage: W = V.submodule([V.gen(0)])
            sage: W.span([V.gen(1)])
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1 0]
            sage: W.submodule([V.gen(1)])
            Traceback (most recent call last):
            ...
            ArithmeticError: gens does not generate a submodule of self
        """
        return FreeModule_submodule_pid(self.ambient_module(), gens, check=check,
                                        already_echelonized=already_echelonized)

    def submodule(self, gens, check=True, already_echelonized=False):
        r"""
        Create the R-submodule of the ambient vector space with given
        generators, where R is the base ring of self.

        INPUT:
            gens  -- a list of free module elements or a free module
            check -- (default: True) whether or not to verify
                      that the gens are in self.

        OUTPUT:
            FreeModule -- the submodule spanned by the vectors in the list gens.
            The basis for the subspace is always put in reduced row echelon form.

        EXAMPLES:
        We create a submodule of $\ZZ^3$:
            sage: M = FreeModule(ZZ, 3)
            sage: B = M.basis()
            sage: W = M.submodule([B[0]+B[1], 2*B[1]-B[2]])
            sage: W
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  1  0]
            [ 0  2 -1]

        We create a submodule of a submodule.
            sage: W.submodule([3*B[0] + 3*B[1]])
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [3 3 0]

        We try to create a submodule that isn't really a submodule,
        which results in an ArithmeticError exception:
            sage: W.submodule([B[0] - B[1]])
            Traceback (most recent call last):
            ...
            ArithmeticError: gens does not generate a submodule of self
        """
# need general HNF for this example:
#
#        Next we create a submodule of a free module over the principal ideal
#        domain $\Q[x]$:
#            sage: R = PolynomialRing(QQ, 'x'); x = R.gen()
#            sage: M = FreeModule(R, 3)
#            sage: B = M.basis()
#            sage: W = M.submodule([x*B[0], 2*B[1]- x*B[2]])
#            sage: W
#            Free module of degree 3 and rank 2 over Integer Ring
#            Echelon basis matrix:
#            [1 1 0]
#            [2 0 1]
#       """
        if is_FreeModule(gens):
            gens = gens.gens()
        if not isinstance(gens, (list, tuple)):
            raise TypeError, "gens (=%s) must be a list or tuple"%gens
        V = self.span(gens, check=check, already_echelonized=already_echelonized)
        if check:
            if not V.is_submodule(self):
                raise ArithmeticError, "gens does not generate a submodule of self"
        return V

    def span_of_basis(self, basis, check=True, already_echelonized=False):
        r"""
        Return the free R-module with the given basis, where R
        is the base ring of self.  Note that this R-module need not
        be a submodule of self, nor even of the ambient space.  It
        must, however, be contained in the ambient vector space, i.e.,
        the ambient space tensored with the fraction field of R.

        EXAMPLES:
            sage: M = FreeModule(ZZ,3)
            sage: W = M.span_of_basis([M([1,2,3])])

        Next we create two free $\Z$-modules, neither of
        which is a submodule of $W$.
            sage: W.span_of_basis([M([2,4,0])])
            Free module of degree 3 and rank 1 over Integer Ring
            User basis matrix:
            [2 4 0]

        The following module isn't even in the ambient space.
            sage: Q = QQ
            sage: W.span_of_basis([ Q('1/5')*M([1,2,0]), Q('1/7')*M([1,1,0]) ])
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1/5 2/5   0]
            [1/7 1/7   0]

        Of course the input basis vectors must be linearly independent.
            sage: W.span_of_basis([ [1,2,0], [2,4,0] ])
            Traceback (most recent call last):
            ...
            ArithmeticError: basis vectors must be linearly independent.
        """
        return FreeModule_submodule_with_basis_pid(
                                    self.ambient_module(), basis, check=check,
                                    already_echelonized=already_echelonized)

    def submodule_with_basis(self, basis, check=True, already_echelonized=False):
        """
        Create the R-submodule of the ambient vector space with given
        basis, where R is the base ring of self.

        INPUT:
            basis -- a list of linearly independent vectors
            check -- whether or not to verify that each gen is in
                     the ambient vector space

        OUTPUT:
            FreeModule -- the R-submodule with given basis

        EXAMPLES:
        First we create a submodule of $\Z^3$:
            sage: M = FreeModule(ZZ, 3)
            sage: B = M.basis()
            sage: W = M.submodule_with_basis([B[0]+B[1], 2*B[1]-B[2]])
            sage: W
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [ 1  1  0]
            [ 0  2 -1]

        """
# need general HNF for this:
#
#        Next we create a submodule of a free module over the principal ideal
#        domain $\Q[x]$:
#            sage: R = PolynomialRing(QQ, 'x'); x = R.gen()
#            sage: M = FreeModule(R, 3)
#            sage: B = M.basis()
#            sage: W = M.submodule_with_basis([x*B[0], 2*B[1]- x*B[2]])
#            sage: W
#            ???
#        """
        V = self.span_of_basis(basis, check=True, already_echelonized=already_echelonized)
        if check:
            if not V.is_submodule(self):
                raise ArithmeticError, "basis does not span a submodule of self"
        return V

    def vector_space_span(self, gens, check=True):
        r"""
        Create the vector subspace of the ambient vector space with
        given generators.

        INPUT:
            gens  -- a list of vector in self
            check -- whether or not to verify that each gen is in
                     the ambient vector space

        OUTPUT:
            a vector subspace

        EXAMPLES:
        We create a $2$-dimensional subspace of a $\Q^3$.
            sage: V = VectorSpace(QQ, 3)
            sage: B = V.basis()
            sage: W = V.vector_space_span([B[0]+B[1], 2*B[1]-B[2]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0  1/2]
            [   0    1 -1/2]

        We create a subspace of a vector space over $\Q(i)$.
            sage: R.<x> = QQ[]
            sage: K = NumberField(x^2 + 1, 'a'); a = K.gen()
            sage: V = VectorSpace(K, 3)
            sage: V.vector_space_span([2*V.gen(0) + 3*V.gen(2)])
            Vector space of degree 3 and dimension 1 over Number Field in a with defining polynomial x^2 + 1
            Basis matrix:
            [  1   0 3/2]

        We use the \code{vector_space_span} command to create a vector subspace of the
        ambient vector space of a submodule of $\Z^3$.
            sage: M = FreeModule(ZZ,3)
            sage: W = M.submodule([M([1,2,3])])
            sage: W.vector_space_span([M([2,3,4])])
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [  1 3/2   2]
        """
        if is_FreeModule(gens):
            gens = gens.gens()
        if not isinstance(gens, (list, tuple)):
            raise TypeError, "gens must be a list"
        return FreeModule_submodule_field(self.ambient_vector_space(),
                                        gens, check=check)

    def vector_space_span_of_basis(self, basis, check=True):
        """
        Create the vector subspace of the ambient vector space with
        given basis.

        INPUT:
            basis -- a list of linearly independent vectors
            check -- whether or not to verify that each gen is in
                     the ambient vector space

        OUTPUT:
            a vector subspace with user-specified basis

        EXAMPLES:
            sage: V = VectorSpace(QQ, 3)
            sage: B = V.basis()
            sage: W = V.vector_space_span_of_basis([B[0]+B[1], 2*B[1]-B[2]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [ 1  1  0]
            [ 0  2 -1]
        """
        return FreeModule_submodule_with_basis_field(
                    self.ambient_vector_space(), basis, check=check)



class FreeModule_generic_field(FreeModule_generic_pid):
    """
    Base class for all free modules over fields.
    """
    def __init__(self, base_field, dimension, degree, sparse=False, inner_product_matrix=None):
        if not isinstance(base_field, field.Field):
            raise TypeError, "base_field (=%s) must be a field"%base_field
        FreeModule_generic_pid.__init__(self, base_ring, dimension, degree, sparse, inner_product_matrix)

    def scale(self, other):
        """
        Return the product of self by the number other, which is the module
        spanned by other times each basis vector.  Since self is a vector
        space this product equals self if other is nonzero, and is the zero
        vector space if other is 0.

        EXAMPLES:
            sage: V = QQ^4
            sage: V.scale(5)
            Vector space of dimension 4 over Rational Field
            sage: V.scale(0)
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []

            sage: W = V.span([[1,1,1,1]])
            sage: W.scale(2)
            Vector space of degree 4 and dimension 1 over Rational Field
            Basis matrix:
            [1 1 1 1]
            sage: W.scale(0)
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []

            sage: V = QQ^4; V
            Vector space of dimension 4 over Rational Field
            sage: V.scale(3)
            Vector space of dimension 4 over Rational Field
            sage: V.scale(0)
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        if other == 0:
            return self.zero_submodule()
        return self

    def __add__(self, other):
        """
        Return the sum of self and other.

        EXAMPLES:
            sage: V = VectorSpace(QQ,3)
            sage: V0 = V.span([V.gen(0)])
            sage: V2 = V.span([V.gen(2)])
            sage: V0 + V2
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 0 1]
        """
        if not isinstance(other, FreeModule_generic_field):
            raise TypeError, "other must be a Vector Space"
        V = self.ambient_vector_space()
        if V != other.ambient_vector_space():
            raise ArithmeticError, "self and other must have the same ambient space"
        return V.span(self.basis() + other.basis())

    def category(self):
        """
        Return the category to which this vector space belongs.

        EXAMPLES:
            sage: V = QQ^4; V.category()
            Category of vector spaces over Rational Field
            sage: V = GF(5)**20; V.category()
            Category of vector spaces over Finite Field of size 5
        """
        import sage.categories.all
        return sage.categories.all.VectorSpaces(self.base_field())


    def echelonized_basis_matrix(self):
        """
        Return basis matrix for self in row echelon form.

        EXAMPLES:
            sage: V = FreeModule(QQ, 3).span_of_basis([[1,2,3],[4,5,6]])
            sage: V.basis_matrix()
            [1 2 3]
            [4 5 6]
            sage: V.echelonized_basis_matrix()
            [ 1  0 -1]
            [ 0  1  2]
        """
        try:
            return self.__echelonized_basis_matrix
        except AttributeError:
            pass
        self.__echelonized_basis_matrix = self.basis_matrix().echelon_form()
        return self.__echelonized_basis_matrix

    def intersection(self, other):
        """
        Return the intersection of self and other, which must be
        R-submodules of a common ambient vector space.

        EXAMPLES:
            sage: V  = VectorSpace(QQ,3)
            sage: W1 = V.submodule([V.gen(0), V.gen(0) + V.gen(1)])
            sage: W2 = V.submodule([V.gen(1), V.gen(2)])
            sage: W1.intersection(W2)
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
            sage: W2.intersection(W1)
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0]
            sage: V.intersection(W1)
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            sage: W1.intersection(V)
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            sage: Z = V.submodule([])
            sage: W1.intersection(Z)
            Vector space of degree 3 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        if not isinstance(other, FreeModule_generic):
            raise TypeError, "other must be a free module"

        if self.ambient_vector_space() != other.ambient_vector_space():
            raise ArithmeticError, "self and other must have the same ambient space."

        if self.base_ring() != other.base_ring():
            # Now other is over a ring R whose fraction field K is the base field of V = self.
            # We compute the intersection using the following algorithm:
            # 1. By explicitly computing the nullspace of the matrix whose rows
            #    are a basis for self, we obtain the matrix over a linear map
            #         phi:  K^n ----> W
            #    with kernel equal to V = self.
            # 2. Compute the kernel over R of Phi restricted to other.  Do this
            #    by clearing denominators, computing the kernel of a matrix with
            #    entries in R, then restoring denominators to the answer.
            K = self.base_ring()
            R = other.base_ring()
            B = self.basis_matrix().transpose()
            W = B.kernel()
            phi = W.basis_matrix().transpose()

            # To restrict phi to other, we multiply the basis matrix for other
            # by phi, thus computing the image of each basis vector.
            X = other.basis_matrix()
            psi = X * phi

            # Now psi is a matrix that defines an R-module morphism from other to some
            # R-module, whose kernel defines the long sought for intersection of self and other.
            L = psi.integer_kernel()

            # Finally the kernel of the intersection has basis the linear combinations of
            # the basis of other given by a basis for L.
            G = L.basis_matrix() * other.basis_matrix()
            return other.span(G.rows())

        # dispense with the three easy cases
        if self == self.ambient_vector_space():
            return other
        elif other == other.ambient_vector_space():
            return self
        elif self.dimension() == 0 or other.dimension() == 0:
            return self.zero_submodule()

        # standard algorithm for computing intersection of general subspaces
        if self.dimension() <= other.dimension():
            V1 = self; V2 = other
        else:
            V1 = other; V2 = self
        A1 = V1.basis_matrix()
        A2 = V2.basis_matrix()
        S  = A1.stack(A2)
        K  = S.kernel()
        n = int(V1.dimension())
        B = [A1.linear_combination_of_rows(v.list()[:n]) for v in K.basis()]
        return self.ambient_vector_space().submodule(B, check=False)

    def is_subspace(self, other):
        """
        True if this vector space is a subspace of other.

        EXAMPLES:
            sage: V = VectorSpace(QQ,3)
            sage: W = V.subspace([V.gen(0), V.gen(0) + V.gen(1)])
            sage: W2 = V.subspace([V.gen(1)])
            sage: W.is_subspace(V)
            True
            sage: W2.is_subspace(V)
            True
            sage: W.is_subspace(W2)
            False
            sage: W2.is_subspace(W)
            True
        """
        return self.is_submodule(other)

    def span(self, gens, check=True, already_echelonized=False):
        """
        Return the K-span of the given list of gens, where K is the
        base field of self.  Note that this span is a subspace of the
        ambient vector space, but need not be a suspace of self.

        INPUT:
            gens -- list of vectors
            check -- bool (default: True): whether or not to coerce entries of gens
                                           into base field
            already_echelonized -- bool (default: False): set this if you know the gens
                                   are already in echelon form

        EXAMPLES:
            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace([[2,3,4]]); W
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 5 2]
            sage: W.span([[1,1,1]])
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 1 1]
        """
        if is_FreeModule(gens):
            gens = gens.gens()
        if not isinstance(gens, (list, tuple)):
            raise TypeError, "gens (=%s) must be a list or tuple"%gens

        return FreeModule_submodule_field(self.ambient_module(), gens, check=check,
                                          already_echelonized=already_echelonized)

    def span_of_basis(self, basis, check=True, already_echelonized=False):
        r"""
        Return the free K-module with the given basis, where K
        is the base field of self.  Note that this span is
        a subspace of the ambient vector space, but need
        not be a suspace of self.

        INPUT:
            basis -- list of vectors
            check -- bool (default: True): whether or not to coerce entries of gens
                                           into base field
            already_echelonized -- bool (default: False): set this if you know the gens
                                   are already in echelon form

        EXAMPLES:
            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace([[2,3,4]]); W
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 5 2]
            sage: W.span_of_basis([[2,2,2], [3,3,0]])
            Vector space of degree 3 and dimension 2 over Finite Field of size 7
            User basis matrix:
            [2 2 2]
            [3 3 0]

        The basis vectors must be linearly independent or an ArithmeticError exception
        is raised.
            sage: W.span_of_basis([[2,2,2], [3,3,3]])
            Traceback (most recent call last):
            ...
            ArithmeticError: basis vectors must be linearly independent.
        """
        return FreeModule_submodule_with_basis_field(
                                    self.ambient_module(), basis, check=check,
                                    already_echelonized=already_echelonized)

    def subspace(self, gens, check=True, already_echelonized=False):
        """
        Return the subspace of self spanned by the elements of gens.

        INPUT:
            gens -- list of vectors
            check -- bool (default: True) verify that gens are all in self.
            already_echelonized -- bool (default: False) set to True if you know the
                                   gens are in Echelon form.

        EXAMPLES:
        First we create a 1-dimensional vector subspace of an ambient $3$-dimensional
        space over the finite field of order $7$.
            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace([[2,3,4]]); W
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 5 2]

        Next we create an invalid subspace, but it's allowed since \code{check=False}.
        This is just equivalent to computing the span of the element.
            sage: W.subspace([[1,1,0]], check=False)
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 1 0]

        With \code{check=True} (the default) the mistake is correctly detected
        and reported with an \code{ArithmeticError} exception.
            sage: W.subspace([[1,1,0]], check=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: gens does not generate a submodule of self
        """
        return self.submodule(gens, check=check, already_echelonized=already_echelonized)

    def subspace_with_basis(self, gens, check=True, already_echelonized=False):
        """
        Same as \code{self.submodule_with_basis(...)}.

        EXAMPLES:
        We create a subspace with a user-defined basis.
            sage: V = VectorSpace(GF(7), 3)
            sage: W = V.subspace_with_basis([[2,2,2], [1,2,3]]); W
            Vector space of degree 3 and dimension 2 over Finite Field of size 7
            User basis matrix:
            [2 2 2]
            [1 2 3]

        We then create a subspace of the subspace with user-defined basis.
            sage: W1 = W.subspace_with_basis([[3,4,5]]); W1
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            User basis matrix:
            [3 4 5]

        Notice how the basis for the same subspace is different if we merely
        use the \code{subspace} command.
            sage: W2 = W.subspace([[3,4,5]]); W2
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 6 4]

        Nonetheless the two subspaces are equal (as mathematical objects):
            sage: W1 == W2
            True
        """
        return self.submodule_with_basis(gens, check=check, already_echelonized=already_echelonized)

    def vector_space(self, base_field=None):
        """
        Return the vector space associated to self.  Since self is a vector
        space this function simply returns self, unless the base field
        is different.

        EXAMPLES:
            sage: V = span(QQ, [[1,2,3]]); V
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
            sage: V.vector_space()
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 2 3]
        """
        if base_field is None:
            return self
        return self.change_ring(base_field)

    def zero_submodule(self):
        """
        Return the zero submodule of self.

        EXAMPLES:
            sage: (QQ^4).zero_submodule()
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        return self.zero_subspace()

    def zero_subspace(self):
        """
        Return the zero subspace of self.

        EXAMPLES:
            sage: (QQ^4).zero_subspace()
            Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        return FreeModule_submodule_field(self.ambient_vector_space(),[], check=False)




###############################################################################
#
# Generic ambient free modules, i.e., of the form R^n for some commutative ring R.
#
###############################################################################

class FreeModule_ambient(FreeModule_generic):
    """
    Ambient free module over a commutative ring.
    """
    def __init__(self, base_ring, rank, sparse=False, inner_product_matrix=None):
        """
        The free module of given rank over the given base_ring.

        INPUT:
            base_ring -- a commutative ring
            rank -- a non-negative integer

        EXAMPLES:
            sage: FreeModule(ZZ, 4)
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
        """
        FreeModule_generic.__init__(self, base_ring, rank,
                                    rank, sparse, inner_product_matrix)


    def echelonized_basis_matrix(self):
        return self.basis_matrix()

    def __cmp__(self, other):
        """
        Compare self and other.

        Modules are ordered by their ambient spaces, then by
        dimension, then in order by their echelon matrices.

        EXAMPLES:
        We compare rank three free modules over the integers and rationals:
            sage: QQ^3 < CC^3
            True
            sage: CC^3 < QQ^3
            False
            sage: CC^3 > QQ^3
            True

            sage: Q = QQ; Z = ZZ
            sage: Q^3 > Z^3
            True
            sage: Q^3 < Z^3
            False
            sage: Z^3 < Q^3
            True
            sage: Z^3 > Q^3
            False
            sage: Q^3 == Z^3
            False
            sage: Q^3 == Q^3
            True

            sage: V = span(QQ, [[1,2,3], [5,6,7], [8,9,10]])
            sage: V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: A = QQ^3
            sage: V < A
            True
            sage: A < V
            False
        """
        if self is other:
            return 0
        if not isinstance(other, FreeModule_generic):
            return cmp(type(self), type(other))
        if isinstance(other, FreeModule_ambient):
            c = cmp(self.rank(), other.rank())
            if c: return c
            c = cmp(self.base_ring(), other.base_ring())
            if not c:
                return c
            try:
                if self.base_ring().is_subring(other.base_ring()):
                    return -1
                elif other.base_ring().is_subring(self.base_ring()):
                    return 1
            except NotImplementedError:
                pass
            return c
        else:  # now other is not ambient; it knows how to do the comparison.
            return -other.__cmp__(self)

    def _repr_(self):
        if self.is_sparse():
            return "Ambient sparse free module of rank %s over %s"%(self.rank(), self.base_ring())
        else:
            return "Ambient free module of rank %s over %s"%(self.rank(), self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of this ambient free module.

        EXAMPLES:
            sage: latex(QQ^3)
            \mathbf{Q}^{3}

            sage: A = GF(5)^20; latex(A)
            \mathbf{F}_{5}^{20}

            sage: A = MPolynomialRing(QQ,3,'x') ^ 20; latex(A)
            (\mathbf{Q}[x_{0}, x_{1}, x_{2}])^{20}
        """
        t = "%s"%latex.latex(self.base_ring())
        if t.find(" ") != -1:
            t = "(%s)"%t
        return "%s^{%s}"%(t, self.rank())

    def is_ambient(self):
        """
        Return \code{True} since this module is an ambient module.

        EXAMPLES:
            sage: A = QQ^5; A.is_ambient()
            True
            sage: A = (QQ^5).span([[1,2,3,4,5]]); A.is_ambient()
            False
        """
        return True



    def ambient_module(self):
        """
        Return self, since self is ambient.

        EXAMPLES:
            sage: A = QQ^5; A.ambient_module()
            Vector space of dimension 5 over Rational Field
            sage: A = ZZ^5; A.ambient_module()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
        """
        return self

    def basis(self):
        """
        Return a basis for this ambient free module.

        OUTPUT:
            Sequence -- an immutable sequence with universe this ambient free module

        EXAMPLES:
            sage: A = ZZ^3; B = A.basis(); B
            [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: B.universe()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        try:
            return self.__basis
        except  AttributeError:
            ZERO = self(0)
            one = self.base_ring()(1)
            w = []
            for n in range(self.rank()):
                v = ZERO.__copy__()
                v[n] = one
                w.append(v)
            self.__basis = basis_seq(self, w)
            return self.__basis

    def echelonized_basis(self):
        """
        Return a basis for this ambient free module in echelon form.

        EXAMPLES:
            sage: A = ZZ^3; A.echelonized_basis()
            [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
        """
        return self.basis()

    def change_ring(self, R):
        """
        Return the ambient free module over R of the same rank as self.

        EXAMPLES:
            sage: A = ZZ^3; A.change_ring(QQ)
            Vector space of dimension 3 over Rational Field
            sage: A = ZZ^3; A.change_ring(GF(5))
            Vector space of dimension 3 over Finite Field of size 5

        For ambient modules any change of rings is defined.
            sage: A = GF(5)**3; A.change_ring(QQ)
            Vector space of dimension 3 over Rational Field
        """
        if self.base_ring() == R:
            return self
        return FreeModule(R, self.rank())


    def linear_combination_of_basis(self, v):
        """
        Return the linear combination of the basis for self
        obtained from the elements of the list v.

        EXAMPLES:
            sage: V = span(ZZ, [[1,2,3], [4,5,6]])
            sage: V
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            [0 3 6]
            sage: V.linear_combination_of_basis([1,1])
            (1, 5, 9)
        """
        return self(v)

    def coordinate_vector(self, v, check=True):
        """
        Write $v$ in terms of the standard basis for self and return the
        resulting coeffcients in a vector over the fraction field of the
        base ring.

        Returns a vector c such that if B is the basis for self, then
                sum c[i] B[i] = v
        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: V = Integers(16)^3
            sage: v = V.coordinate_vector([1,5,9]); v
            (1, 5, 9)
            sage: v.parent()
            Ambient free module of rank 3 over Ring of integers modulo 16
        """
        return self(v)

    def echelon_coordinate_vector(self, v, check=True):
        r"""
        Same as \code{self.coordinate_vector(v)}, since self is an ambient free module.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list
        """
        return self.coordinate_vector(v, check=check)

    def echelon_coordinates(self, v, check=True):
        return self.coordinate_vector(v, check=check)


###############################################################################
#
# Ambient free modules over an integral domain.
#
###############################################################################

class FreeModule_ambient_domain(FreeModule_ambient):
    """
    Ambient free module over an integral domain.
    """
    def __init__(self, base_ring, rank, sparse=False, inner_product_matrix=None):
        """
        EXAMPLES:
            sage: FreeModule(PolynomialRing(GF(5),'x'), 3)
            Ambient free module of rank 3 over the principal ideal domain
            Univariate Polynomial Ring in x over Finite Field of size 5
        """
        FreeModule_ambient.__init__(self, base_ring, rank, sparse, inner_product_matrix)

    def _repr_(self):
        if self.is_sparse():
            return "Ambient sparse free module of rank %s over the integral domain %s"%(
                self.rank(), self.base_ring())
        else:
            return "Ambient free module of rank %s over the integral domain %s"%(
                self.rank(), self.base_ring())

    def base_field(self):
        """
        Return the fraction field of the base ring of self.

        EXAMPLES:
            sage: M = ZZ^3;  M.base_field()
            Rational Field
            sage: M = PolynomialRing(GF(5),'x')^3;  M.base_field()
            Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 5
        """
        return self.base_ring().fraction_field()

    def ambient_vector_space(self):
        """
        Returns the ambient vector space, which is this free module tensored
        with its fraction field.

        EXAMPLES:
            sage: M = ZZ^3;  M.ambient_vector_space()
            Vector space of dimension 3 over Rational Field
        """
        try:
            return self.__ambient_vector_space
        except AttributeError:
            if self._inner_product_is_dot_product():
                mat = None
            else:
                mat = self.inner_product_matrix()
            self.__ambient_vector_space = FreeModule(self.base_field(), self.rank(),
                                             sparse = self.is_sparse(),
                                             inner_product_matrix = mat)
            return self.__ambient_vector_space

    def coordinate_vector(self, v, check=True):
        """
        Write $v$ in terms of the standard basis for self and return the
        resulting coeffcients in a vector over the fraction field of the
        base ring.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list

        Returns a vector c such that if B is the basis for self, then
                sum c[i] B[i] = v
        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: V = ZZ^3
            sage: v = V.coordinate_vector([1,5,9]); v
            (1, 5, 9)
            sage: v.parent()
            Vector space of dimension 3 over Rational Field
        """
        return self.ambient_vector_space()(v)

    def vector_space(self, base_field=None):
        """
        Returns the vector space obtained from self by tensoring with
        the fraction field of the base ring and extending to the field.

        EXAMPLES:
            sage: M = ZZ^3;  M.vector_space()
            Vector space of dimension 3 over Rational Field
        """
        if base_field is None:
            R = self.base_ring()
            return self.change_ring(R.fraction_field())
        else:
            return self.change_ring(base_field)

###############################################################################
#
# Ambient free modules over a principal ideal domain.
#
###############################################################################

class FreeModule_ambient_pid(FreeModule_generic_pid, FreeModule_ambient_domain):
    """
    Ambient free module over a principal ideal domain.
    """
    def __init__(self, base_ring, rank, sparse=False, inner_product_matrix=None):
        """
        Create the ambient free module of given rank over the given
        principal ideal domain

        INPUT:
            base_ring -- a principal ideal domain
            rank -- a non-negative integer
            sparse -- bool (default: False)
            inner_product_matrix -- bool (default: None)

        EXAMPLES:
            sage: ZZ^3
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        FreeModule_ambient_domain.__init__(self, base_ring, rank, sparse, inner_product_matrix)

    def _repr_(self):
        if self.is_sparse():
            return "Ambient sparse free module of rank %s over the principal ideal domain %s"%(
                self.rank(), self.base_ring())
        else:
            return "Ambient free module of rank %s over the principal ideal domain %s"%(
                self.rank(), self.base_ring())


###############################################################################
#
# Ambient free modules over a field (i.e., a vector space).
#
###############################################################################

class FreeModule_ambient_field(FreeModule_generic_field, FreeModule_ambient_pid):
    def __init__(self, base_field, dimension, sparse=False, inner_product_matrix=None):
        FreeModule_ambient_pid.__init__(self, base_field, dimension, sparse, inner_product_matrix)

    def _repr_(self):
        if self.is_sparse():
            return "Sparse vector space of dimension %s over %s"%(self.dimension(), self.base_ring())
        else:
            return "Vector space of dimension %s over %s"%(self.dimension(), self.base_ring())

    def ambient_vector_space(self):
        """
        Returns the ambient vector space.

        EXAMPLES:
            sage: M = QQ^3
            sage: M.ambient_vector_space()
            Vector space of dimension 3 over Rational Field
        """
        return self

    def base_field(self):
        """
        Returns the base field of this vector space.

        EXAMPLES:
            sage: M = QQ^3
            sage: M.base_field()
            Rational Field
        """
        return self.base_ring()


###############################################################################
#
# R-Submodule of K^n where K is the fraction field of a principal ideal domain $R$.
#
###############################################################################


class FreeModule_submodule_with_basis_pid(FreeModule_generic_pid):
    """
    An $R$-submodule of $K^n$ with distinguished basis, where $K$ is
    the fraction field of a principal ideal domain $R$.
    """
    def __init__(self, ambient, basis,
                 check=True,
                 echelonize=False,
                 inner_product_matrix=None,
                 echelonized_basis=None,
                 already_echelonized=False):
        """
        Create a free module with basis over a PID.

        EXAMPLES:
            sage: M = ZZ^3
            sage: W = M.span_of_basis([[1,2,3],[4,5,6]]); W
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1 2 3]
            [4 5 6]
        """
        if not isinstance(ambient, FreeModule_ambient_pid):
            raise TypeError, "ambient (=%s) must be ambient."%ambient
        if not isinstance(basis, (list, tuple)):
            raise TypeError, "basis (=%s) must be a list"%basis
        if check:
            V = ambient.vector_space()
            basis = [V(x) for x in basis]

        self.__ambient_module = ambient

        V = ambient.ambient_vector_space()
        if check:
            for i in range(len(basis)):
                x = basis[i]
                if not (x in V):
                    try:
                        basis[i] = V(x)
                    except TypeError:
                        raise TypeError, "each element of basis must be in the ambient vector space"
        self.__ambient_vector_space = V
        R = ambient.base_ring()

        if echelonize and not already_echelonized:
            self.__basis = basis
            basis = self._echelonized_basis(ambient, basis)

        FreeModule_generic.__init__(self, R, len(basis), ambient.degree(),
                            ambient.is_sparse(), inner_product_matrix=inner_product_matrix)

        C = self.element_class()
        try:
            w = [C(self, x.list(),
                          coerce=False, copy=True) for x in basis]
        except TypeError:
            C = element_class(R.fraction_field(), self.is_sparse())
            w = [C(self, x.list(),
                          coerce=False, copy=True) for x in basis]

        self.__basis = basis_seq(self, w)

        if echelonized_basis != None:

            self.__echelonized_basis = basis_seq(self, echelonized_basis)

        else:

            if echelonize or already_echelonized:
                self.__echelonized_basis = self.__basis
            else:
                w = self._echelonized_basis(ambient, basis)
                self.__echelonized_basis = basis_seq(self, w)

        if check and len(basis) != len(self.__echelonized_basis):
            raise ArithmeticError, "basis vectors must be linearly independent."

    def echelonized_basis_matrix(self):
        """
        Return basis matrix for self in row echelon form.

        EXAMPLES:
            sage: V = FreeModule(ZZ, 3).span_of_basis([[1,2,3],[4,5,6]])
            sage: V.basis_matrix()
            [1 2 3]
            [4 5 6]
            sage: V.echelonized_basis_matrix()
            [1 2 3]
            [0 3 6]
        """
        try:
            return self.__echelonized_basis_matrix
        except AttributeError:
            pass
        self._echelonized_basis(self.ambient_module(), self.__basis)
        return self.__echelonized_basis_matrix

    def _echelonized_basis(self, ambient, basis):
        # Return the first rank rows (i.e., the nonzero rows).
        d = self._denominator(basis)
        MAT = sage.matrix.matrix_space.MatrixSpace(ambient.base_ring(),
                        len(basis), ambient.degree(),
                        sparse = ambient.is_sparse())
        if d != 1:
            basis = [x*d for x in basis]
        A = MAT(basis)
        E = A.echelon_form()
        if d != 1:
            E = E.matrix_over_field()*(~d)   # divide out denominator
        r = E.rank()
        if r < E.nrows():
            E = E.matrix_from_rows(range(r))
        self.__echelonized_basis_matrix = E
        return E.rows()

    def __cmp__(self, other):
        r"""
        Compare self and other.

        Modules are ordered by their ambient spaces, then by
        dimension, then in order by their echelon matrices.

        NOTE: Use the \code{is_submodule} to determine if one module
        is a submodule of another.

        EXAMPLES:
        First we compare two equal vector spaces.
            sage: V = span(QQ, [[1,2,3], [5,6,7], [8,9,10]])
            sage: W = span(QQ, [[5,6,7], [8,9,10]])
            sage: V == W
            True

        Next we compare a one dimensional space to the two dimensional space
        defined above.
            sage: M = span(QQ, [[5,6,7]])
            sage: V == M
            False
            sage: M < V
            True
            sage: V < M
            False

        We compare a $\Z$-module to the one-dimensional space above.
            sage: V = span(ZZ, [[5,6,7]]).scale(1/11);  V
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [5/11 6/11 7/11]
            sage: V < M
            True
            sage: M < V
            False
        """
        if self is other:
            return 0
        if not isinstance(other, FreeModule_generic):
            return cmp(type(self), type(other))
        c = cmp(self.ambient_vector_space(), other.ambient_vector_space())
        if c: return c
        c = cmp(self.dimension(), other.dimension())
        if c: return c
        # We use self.echelonized_basis_matrix() == other.echelonized_basis_matrix()
        # with the matrix to avoid a circular reference.
        return cmp(self.echelonized_basis_matrix(), other.echelonized_basis_matrix())

##         if self.base_ring() == other.base_ring() and \
##                self.echelonized_basis_matrix() == other.echelonized_basis_matrix():
##             return 0
##         try:
##             if self.is_submodule(other):
##                 return -1
##             else:
##                 return 1
##         except NotImplementedError:
##             return 1

    def _denominator(self, B):  # for internal use only!
        if len(B) == 0:
            return 1
        d = B[0].denominator()
        for x in B[1:]:
            d = d.lcm(x.denominator())
        return d

    def _repr_(self):
        if self.is_sparse():
            s = "Sparse free module of degree %s and rank %s over %s\n"%(
                self.degree(), self.rank(), self.base_ring()) + \
                "User basis matrix:\n%s"%self.basis_matrix()
        else:
            s = "Free module of degree %s and rank %s over %s\n"%(
                self.degree(), self.rank(), self.base_ring()) + \
                "User basis matrix:\n%s"%self.basis_matrix()
        return s

    def _latex_(self):
        r"""
        Return latex representation of this free module.

        EXAMPLES:
            sage: A = ZZ^3
            sage: M = A.span_of_basis([[1,2,3],[4,5,6]])
            sage: M._latex_()
            '\\mbox{\\rm RowSpan}_{\\mathbf{Z}}\\left(\\begin{array}{rrr}\n1&2&3\\\\\n4&5&6\n\\end{array}\\right)'
        """
        return "\\mbox{\\rm RowSpan}_{%s}%s"%(latex.latex(self.base_ring()), latex.latex(self.basis_matrix()))
        #return "\\mbox{\\rm Span}_{%s}(%s)"%(latex.latex(self.base_ring()), latex.latex(self.basis()))
#        return "\\text{Free module over %s spanned by the rows of:}%s."%(
#            latex.latex(self.base_ring()),latex.latex(self.basis_matrix()))


    def ambient_module(self):
        """
        Return the ambient module related to the $R$-module self, which
        was used when creating this module, and is of the form $R^n$.   Note
        that self need not be contained in the ambient module, though self
        will be contained in the ambient vector space.

        EXAMPLES:
            sage: A = ZZ^3
            sage: M = A.span_of_basis([[1,2,'3/7'],[4,5,6]])
            sage: M
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [  1   2 3/7]
            [  4   5   6]
            sage: M.ambient_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: M.is_submodule(M.ambient_module())
            False
        """
        return self.__ambient_module

    def echelon_coordinates(self, v, check=True):
        r"""
        Write $v$ in terms of the echelonized basis for self.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list

        Returns a list $c$ such that if $B$ is the basis for self, then
        $$
                \sum c_i B_i = v.
        $$
        If $v$ is not in self, raises an \code{ArithmeticError} exception.

        EXAMPLES:
            sage: A = ZZ^3
            sage: M = A.span_of_basis([[1,2,'3/7'],[4,5,6]])
            sage: M.coordinates([8,10,12])
            [0, 2]
            sage: M.echelon_coordinates([8,10,12])
            [8, -2]
            sage: B = M.echelonized_basis(); B
            [
            (1, 2, 3/7),
            (0, 3, -30/7)
            ]
            sage: 8*B[0] - 2*B[1]
            (8, 10, 12)

        We do an example with a sparse vector space:
            sage: V = VectorSpace(QQ,5, sparse=True)
            sage: W = V.subspace_with_basis([[0,1,2,0,0], [0,-1,0,0,-1/2]])
            sage: W.echelonized_basis()
            [
            (0, 1, 0, 0, 1/2),
            (0, 0, 1, 0, -1/4)
            ]
            sage: W.echelon_coordinates([0,0,2,0,-1/2])
            [0, 2]
        """
        if not isinstance(v, free_module_element.FreeModuleElement):
            v = self.ambient_vector_space()(v)
        elif v.degree() != self.degree():
            raise ArithmeticError, "vector is not in free module"
        # Find coordinates of v with respect to rref basis.
        E = self.echelonized_basis_matrix()
        P = E.pivots()
        w = v.list_from_positions(P)
        # Next use the transformation matrix from the rref basis
        # to the echelon basis.
        T = self._rref_to_echelon_matrix()
        x = T.linear_combination_of_rows(w).list(copy=False)
        if not check:
            return x
        if v.parent() is self:
            return x
        lc = E.linear_combination_of_rows(x)
        if lc != v:
            raise ArithmeticError, "vector is not in free module"
        return x

    def user_to_echelon_matrix(self):
        """
        Return matrix that transforms a vector written with respect to the user basis
        of self to one written with respect to the echelon basis.  The matrix acts
        from the right, as is usual in SAGE.

        EXAMPLES:
            sage: A = ZZ^3
            sage: M = A.span_of_basis([[1,2,3],[4,5,6]])
            sage: M.echelonized_basis()
            [
            (1, 2, 3),
            (0, 3, 6)
            ]
            sage: M.user_to_echelon_matrix()
            [ 1  0]
            [ 4 -1]

        The vector $v=(5,7,9)$ in $M$ is $(1,1)$ with respect to the user basis.
        Multiplying the above matrix on the right by this vector yields $(5,-1)$,
        which has components the coordinates of $v$ with respect to the echelon basis.
            sage: v0,v1 = M.basis(); v = v0+v1
            sage: e0,e1 = M.echelonized_basis()
            sage: v
            (5, 7, 9)
            sage: 5*e0 + (-1)*e1
            (5, 7, 9)
        """
        try:
            return self.__user_to_echelon_matrix
        except AttributeError:
            if self.base_ring().is_field():
                self.__user_to_echelon_matrix = self._user_to_rref_matrix()
            else:
                rows = sum([self.echelon_coordinates(b,check=False) for b in self.basis()], [])
                M = sage.matrix.matrix_space.MatrixSpace(self.base_ring().fraction_field(),
                                             self.dimension(),
                                             sparse = self.is_sparse())
                self.__user_to_echelon_matrix = M(rows)
        return self.__user_to_echelon_matrix


    def echelon_to_user_matrix(self):
        """
        Return matrix that transforms the echelon basis to the user basis of self.
        This is a matrix $A$ such that if $v$ is a vector written with respect to the echelon
        basis for self then $vA$ is that vector written with respect to the user
        basis of self.

        EXAMPLES:
            sage: V = QQ^3
            sage: W = V.span_of_basis([[1,2,3],[4,5,6]])
            sage: W.echelonized_basis()
            [
            (1, 0, -1),
            (0, 1, 2)
            ]
            sage: A = W.echelon_to_user_matrix(); A
            [-5/3  2/3]
            [ 4/3 -1/3]

        The vector $(1,1,1)$ has coordinates $v=(1,1)$ with respect to the echelonized
        basis for self.  Multiplying $vA$ we find the coordinates of this vector with
        respect to the user basis.
            sage: v = vector(QQ, [1,1]); v
            (1, 1)
            sage: v * A
            (-1/3, 1/3)
            sage: u0, u1 = W.basis()
            sage: (-u0 + u1)/3
            (1, 1, 1)
        """
        try:
            return self.__echelon_to_user_matrix
        except AttributeError:
            self.__echelon_to_user_matrix = ~self.user_to_echelon_matrix()
            return self.__echelon_to_user_matrix

    ## rref means the echelon form but over a field.
    def _user_to_rref_matrix(self):
        try:
            return self.__user_to_rref_matrix
        except AttributeError:
            A = self.basis_matrix()
            P = self.echelonized_basis_matrix().pivots()
            T = A.matrix_from_columns(P)
            self.__user_to_rref_matrix = T
        return self.__user_to_rref_matrix

    def _rref_to_user_matrix(self):
        try:
            return self.__rref_to_user_matrix
        except AttributeError:
            self.__rref_to_user_matrix = ~self._user_to_rref_matrix()
            return self.__rref_to_user_matrix

    def _echelon_to_rref_matrix(self):
        try:
            return self.__echelon_to_rref_matrix
        except AttributeError:
            A = self.echelonized_basis_matrix()
            T = A.matrix_from_columns(A.pivots())
            self.__echelon_to_rref_matrix = T
        return self.__echelon_to_rref_matrix

    def _rref_to_echelon_matrix(self):
        try:
            return self.__rref_to_echelon_matrix
        except AttributeError:
            self.__rref_to_echelon_matrix = ~self._echelon_to_rref_matrix()
            return self.__rref_to_echelon_matrix


    def vector_space(self, base_field=None):
        """
        Return the vector spaces associated to this free module via tensor product
        with the fraction field of the base ring.

        EXAMPLES:
            sage: A = ZZ^3; A
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: A.vector_space()
            Vector space of dimension 3 over Rational Field
            sage: M = A.span_of_basis([['1/3',2,'3/7'],[4,5,6]]); M
            Free module of degree 3 and rank 2 over Integer Ring
            User basis matrix:
            [1/3   2 3/7]
            [  4   5   6]
            sage: M.vector_space()
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1/3   2 3/7]
            [  4   5   6]
        """
        if base_field is None:
            return self.change_ring(self.base_ring().fraction_field())
        return self.change_ring(base_field)

    def ambient_vector_space(self):
        """
        Return the ambient vector space in which this free module is embedded.

        EXAMPLES:
            sage: V = ZZ^3
            sage: M = V.span_of_basis([[1,2,'1/5']])
            sage: M
            Free module of degree 3 and rank 1 over Integer Ring
            User basis matrix:
            [  1   2 1/5]
            sage: M.ambient_vector_space()
            Vector space of dimension 3 over Rational Field
        """
        return self.__ambient_vector_space

    def basis(self):
        """
        Return the user basis for this free module.

        EXAMPLES:
            sage: V = ZZ^3
            sage: V.basis()
            [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
            sage: M = V.span_of_basis([['1/8',2,1]])
            sage: M.basis()
            [
            (1/8, 2, 1)
            ]
        """
        return self.__basis

    def change_ring(self, R):
        """
        Return the free module over R obtained by coercing each
        element of self into a vector over the fraction field of R,
        then taking the resulting R-module.  Raises a TypeError
        if coercion is not possible.

        INPUT:
            R -- a principal ideal domain

        EXAMPLES:
            sage: V = QQ^3
            sage: W = V.subspace([[2,'1/2', 1]])
            sage: W.change_ring(GF(7))
            Vector space of degree 3 and dimension 1 over Finite Field of size 7
            Basis matrix:
            [1 2 4]
        """
        if self.base_ring() == R:
            return self
        K = R.fraction_field()
        V = VectorSpace(K, self.degree())
        B = [V(b) for b in self.basis()]
        M = FreeModule(R, self.degree())
        if self.has_user_basis():
            return M.span_of_basis(B)
        else:
            return M.span(B)

    def coordinate_vector(self, v, check=True):
        """
        Write $v$ in terms of the user basis for self.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list

        Returns a vector c such that if B is the basis for self, then
                sum c[i] B[i] = v
        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: V = ZZ^3
            sage: M = V.span_of_basis([['1/8',2,1]])
            sage: M.coordinate_vector([1,16,8])
            (8)
        """
        # First find the coordinates of v wrt echelon basis.
        w = self.echelon_coordinate_vector(v, check=check)
        # Next use transformation matrix from echelon basis to
        # user basis.
        T = self.echelon_to_user_matrix()
        return T.linear_combination_of_rows(w)

    def echelonized_basis(self):
        """
        Return the basis for self in echelon form.

        EXAMPLES:
            sage: V = ZZ^3
            sage: M = V.span_of_basis([['1/2',3,1], [0,'1/6',0]])
            sage: M.basis()
            [
            (1/2, 3, 1),
            (0, 1/6, 0)
            ]
            sage: B = M.echelonized_basis(); B
            [
            (1/2, 0, 1),
            (0, 1/6, 0)
            ]
            sage: V.span(B) == M
            True
        """
        return self.__echelonized_basis

    def echelon_coordinate_vector(self, v, check=True):
        """
        Write $v$ in terms of the user basis for self.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.

        Returns a vector c such that if B is the echelonized basis for self, then
                sum c[i] B[i] = v
        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: V = ZZ^3
            sage: M = V.span_of_basis([['1/2',3,1], [0,'1/6',0]])
            sage: B = M.echelonized_basis(); B
            [
            (1/2, 0, 1),
            (0, 1/6, 0)
            ]
            sage: M.echelon_coordinate_vector(['1/2', 3, 1])
            (1, 18)
        """
        return FreeModule(self.base_ring().fraction_field(), self.rank())(self.echelon_coordinates(v, check=check))

    def has_user_basis(self):
        """
        Return \code{True} if the basis of this free module is specified by the user,
        as opposed to being the default echelon form.

        EXAMPLES:
            sage: V = ZZ^3; V.has_user_basis()
            False
            sage: M = V.span_of_basis([[1,3,1]]); M.has_user_basis()
            True
            sage: M = V.span([[1,3,1]]); M.has_user_basis()
            False
        """
        return True

    def linear_combination_of_basis(self, v):
        """
        Return the linear combination of the basis for self
        obtained from the coordinates of v.

        EXAMPLES:
            sage: V = span(ZZ, [[1,2,3], [4,5,6]]); V
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            [0 3 6]
            sage: V.linear_combination_of_basis([1,1])
            (1, 5, 9)
        """
        return self.basis_matrix().linear_combination_of_rows(v)


class FreeModule_submodule_pid(FreeModule_submodule_with_basis_pid):
    """
    An $R$-submodule of $K^n$ where $K$ is the fraction field of a
    principal ideal domain $R$.

    EXAMPLES:
        sage: M = ZZ^3
        sage: W = M.span_of_basis([[1,2,3],[4,5,19]]); W
        Free module of degree 3 and rank 2 over Integer Ring
        User basis matrix:
        [ 1  2  3]
        [ 4  5 19]

    We can save and load submodules and elements.

        sage: loads(W.dumps()) == W
        True
        sage: v = W.0 + W.1
        sage: loads(v.dumps()) == v
        True


    """
    def __init__(self, ambient, gens, check=True, inner_product_matrix=None,
                 already_echelonized=False):
        """
        Create an embedded free module over a PID.

        EXAMPLES:
            sage: V = ZZ^3
            sage: W = V.span([[1,2,3],[4,5,6]])
            sage: W
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            [0 3 6]
        """
        FreeModule_submodule_with_basis_pid.__init__(self, ambient,
                                                     gens, echelonize=True,
                                                     inner_product_matrix=inner_product_matrix,
                                                     already_echelonized=already_echelonized)

    def _repr_(self):
        if self.is_sparse():
            s = "Sparse free module of degree %s and rank %s over %s\n"%(
                self.degree(), self.rank(), self.base_ring()) + \
                "Echelon basis matrix:\n%s"%self.basis_matrix()
        else:
            s = "Free module of degree %s and rank %s over %s\n"%(
                self.degree(), self.rank(), self.base_ring()) + \
                "Echelon basis matrix:\n%s"%self.basis_matrix()
        return s

    def coordinate_vector(self, v, check=True):
        """
        Write $v$ in terms of the user basis for self.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list

        Returns a vector c such that if B is the basis for self, then
                sum c[i] B[i] = v
        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: V = ZZ^3
            sage: W = V.span_of_basis([[1,2,3],[4,5,6]])
            sage: W.coordinate_vector([1,5,9])
            (5, -1)
        """
        return self.echelon_coordinate_vector(v, check=check)

    def has_user_basis(self):
        r"""
        Return \code{True} if the basis of this free module is specified by the user,
        as opposed to being the default echelon form.

        EXAMPLES:
            sage: A = ZZ^3; A
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: A.has_user_basis()
            False
            sage: W = A.span_of_basis([[2,'1/2',1]])
            sage: W.has_user_basis()
            True
            sage: W = A.span([[2,'1/2',1]])
            sage: W.has_user_basis()
            False
        """
        return False


class FreeModule_submodule_with_basis_field(FreeModule_generic_field, FreeModule_submodule_with_basis_pid):
    """
    A vector subspace with distinguished basis.

    EXAMPLES:
        sage: M = QQ^3; W = M.span([[1,2,3],[4,5,19]]); W
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [   1    0 23/3]
        [   0    1 -7/3]

    We can load and save submodules:

        sage: loads(W.dumps()) == W
        True

        sage: K.<x> = FractionField(PolynomialRing(QQ,'x'))
        sage: M = K^3; W = M.span_of_basis([[1,1,x]])
        sage: loads(W.dumps()) == W
        True
    """
    def __init__(self, ambient, basis, check=True, echelonize=False,
                 inner_product_matrix=None, echelonized_basis=None,
                 already_echelonized=False):
        """
        Create a vector space with given basis.

        EXAMPLES:
            sage: V = QQ^3
            sage: W = V.span_of_basis([[1,2,3],[4,5,6]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            User basis matrix:
            [1 2 3]
            [4 5 6]
        """
        FreeModule_submodule_with_basis_pid.__init__(self, ambient,
                                                     basis,
                                                     check,
                                                     echelonize=echelonize,
                                                     inner_product_matrix=inner_product_matrix,
                                                     echelonized_basis=echelonized_basis,
                                                     already_echelonized=already_echelonized)

    def _repr_(self):
        if self.is_sparse():
            return "Sparse vector space of degree %s and dimension %s over %s\n"%(
                    self.degree(), self.dimension(), self.base_field()) + \
                    "User basis matrix:\n%s"%self.basis_matrix()
        else:
            return "Vector space of degree %s and dimension %s over %s\n"%(
                    self.degree(), self.dimension(), self.base_field()) + \
                    "User basis matrix:\n%s"%self.basis_matrix()

    def _denominator(self, B):
        # for internal use only
        return 1

    def _echelonized_basis(self, ambient, basis):
        MAT = sage.matrix.matrix_space.MatrixSpace(ambient.base_ring(),
                        len(basis), ambient.degree(),
                        sparse = ambient.is_sparse())
        A = MAT(basis)
        E = A.echelon_form()
        # Return the first rank rows (i.e., the nonzero rows).
        return E.rows()[:E.rank()]

    def is_ambient(self):
        """
        Return False since this is not an ambient module.

        EXAMPLES:
            sage: V = QQ^3
            sage: V.is_ambient()
            True
            sage: W = V.span_of_basis([[1,2,3],[4,5,6]])
            sage: W.is_ambient()
            False
        """
        return False



class FreeModule_submodule_field(FreeModule_submodule_with_basis_field):
    """
    An embedded vector subspace with echelonized basis.
    """
    def __init__(self, ambient, gens, check=True,
                 inner_product_matrix=None, already_echelonized=False):
        """
        Create an embedded vector subspace with echelonized basis.

        EXAMPLES:
            sage: V = QQ^3
            sage: W = V.span([[1,2,3],[4,5,6]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
        """
        if is_FreeModule(gens):
            gens = gens.gens()
        if not isinstance(gens, (list, tuple)):
            raise TypeError, "gens (=%s) must be a list or tuple"%gens
        FreeModule_submodule_with_basis_field.__init__(self, ambient, gens, check,
                                                       echelonize = not already_echelonized,
                                                       inner_product_matrix=inner_product_matrix,
                                                       already_echelonized = already_echelonized)

    def _repr_(self):
        if self.is_sparse():
            return "Sparse vector space of degree %s and dimension %s over %s\n"%(
                self.degree(), self.dimension(), self.base_field()) + \
                "Basis matrix:\n%s"%self.basis_matrix()
        else:
            return "Vector space of degree %s and dimension %s over %s\n"%(
                    self.degree(), self.dimension(), self.base_field()) + \
                    "Basis matrix:\n%s"%self.basis_matrix()

    def echelon_coordinates(self, v, check=True):
        """
        Write $v$ in terms of the user basis for self.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list

        Returns a list $c$ such that if $B$ is the basis for self, then
        $$
                \sum c_i B_i = v.
        $$
        If $v$ is not in self, raises an \code{ArithmeticError} exception.

            An embedded vector subspace with echelonized basis.

        EXAMPLES:
            sage: V = QQ^3
            sage: W = V.span([[1,2,3],[4,5,6]])
            sage: W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: W.echelon_coordinates([1,5,9])
            [1, 5]
        """
        if not isinstance(v, free_module_element.FreeModuleElement):
            v = self.ambient_vector_space()(v)
        if v.degree() != self.degree():
            raise ArithmeticError, "v (=%s) is not in self"%v
        E = self.echelonized_basis_matrix()
        P = E.pivots()
        if len(P) == 0:
            if check and v != 0:
                raise ArithmeticError, "vector is not in free module"
            return []
        w = v.list_from_positions(P)
        if not check:
            # It's really really easy.
            return w
        if v.parent() is self:   # obvious that v is really in here.
            return w
        # the "linear_combination_of_rows" call dominates the runtime
        # of this function, in the check==False case when the parent
        # of v is not self.
        lc = E.linear_combination_of_rows(w)
        if lc != v:
            raise ArithmeticError, "vector is not in free module"
        return w

    def coordinate_vector(self, v, check=True):
        """
        Write $v$ in terms of the user basis for self.

        INPUT:
            v -- vector
            check -- bool (default: True); if True, also verify that v is really in self.
        OUTPUT:
            list

        Returns a vector c such that if B is the basis for self, then
                sum c[i] B[i] = v
        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: V = QQ^3
            sage: W = V.span([[1,2,3],[4,5,6]]); W
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: W.coordinate_vector([1,5,9])
            (1, 5)

            sage: V = VectorSpace(QQ,5, sparse=True)
            sage: W = V.subspace([[0,1,2,0,0], [0,-1,0,0,-1/2]])
            sage: W.coordinate_vector([0,0,2,0,-1/2])
            (0, 2)
        """
        return self.echelon_coordinate_vector(v, check=check)

    def has_user_basis(self):
        """
        Return \code{True} if the basis of this free module is specified by the user,
        as opposed to being the default echelon form.

        EXAMPLES:
            sage: V = QQ^3
            sage: W = V.subspace([[2,'1/2', 1]])
            sage: W.has_user_basis()
            False
            sage: W = V.subspace_with_basis([[2,'1/2',1]])
            sage: W.has_user_basis()
            True
        """
        return False


def basis_seq(V, w):
    return Sequence(w, universe=V, check = False, immutable=True, cr=True)






class RealDoubleVectorSpace_class(FreeModule_ambient_field):
    def __init__(self,n):
        FreeModule_ambient_field.__init__(self,sage.rings.real_double.RDF,n)
        self._element_class = sage.modules.real_double_vector.RealDoubleVectorSpaceElement

    def coordinates(self,v):
        return v

class ComplexDoubleVectorSpace_class(FreeModule_ambient_field):
    def __init__(self,n):
        FreeModule_ambient_field.__init__(self,sage.rings.complex_double.CDF,n)
        self._element_class = sage.modules.complex_double_vector.ComplexDoubleVectorSpaceElement

    def coordinates(self,v):
        return v



######################################################

def element_class(R, is_sparse):
    if sage.rings.integer_ring.is_IntegerRing(R) and not is_sparse:
        from vector_integer_dense import Vector_integer_dense
        return Vector_integer_dense
    else:
        if is_sparse:
            return free_module_element.FreeModuleElement_generic_sparse
        else:
            return free_module_element.FreeModuleElement_generic_dense
    raise NotImplementedError
