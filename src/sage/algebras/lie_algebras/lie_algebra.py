"""
Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.indexed_generators import IndexedGenerators
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper

from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.rings import Rings
from sage.categories.morphism import Morphism, SetMorphism
from sage.categories.map import Map
from sage.categories.homset import Hom

from sage.algebras.free_algebra import FreeAlgebra, is_FreeAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import (LieAlgebraElement,
    LieAlgebraElementWrapper)
from sage.rings.all import ZZ
from sage.rings.ring import Ring
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.infinity import infinity
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.modules.free_module import FreeModule, span
from sage.combinat.root_system.cartan_type import CartanType, CartanType_abstract
from sage.sets.family import Family, AbstractFamily
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet

class LieAlgebra(Parent, UniqueRepresentation): # IndexedGenerators):
    r"""
    A Lie algebra `L` over a base ring `R`.

    A Lie algebra is an algebra with a bilinear operation called Lie bracket
    `[\cdot, \cdot] : L \times L \to L` such that `[x, x] = 0` and
    the following relation holds:

    .. MATH::

        \bigl[x, [y, z] \bigr] + \bigl[ y, [z, x] \bigr]
        + \bigl[ z, [x, y] \bigr] = 0.

    This relation is known as the *Jacobi identity* (or sometimes the Jacobi
    relation). We note that from `[x, x] = 0`, we have `[x + y, x + y] = 0`.
    Next from bilinearity, we see that

    .. MATH::

        0 = [x + y, x + y] = [x, x] + [x, y] + [y, x] + [y, y]
        = [x, y] + [y, x],

    thus `[x, y] = -[y, x]` and the Lie bracket is antisymmetric.

    Lie algebras are closely related to Lie groups. Let `G` be a Lie group
    and fix some `g \in G`, we can construct the Lie algebra `L` of `G` by
    considering the tangent space at `g`. We can also (partially) recover `G`
    from `L` by using what is known as the exponential map.

    Given any associative algebra `A`, we can construct a Lie algebra `L` by
    defining the Lie bracket to be the commutator `[a, b] = ab - ba`. We call
    an associative algebra `A` which contains `L` in this fashion an
    *enveloping algebra*. We can the embedding which sends the Lie bracket to
    the commutator a Lie embedding. Now if we are given a Lie algebra `L`, we
    can construct an enveloping algebra `U_L` with Lie embedding `h : L \to
    U_L` which has the following universal property: for any enveloping
    algebra `A` with Lie embedding `f : L \to A`, there exists a unique unital
    algebra homomorphism `g : U_L \to A` such that `f = g \circ h`. The
    algebra `U_L` is known as the *universal enveloping algebra*.

    EXAMPLES:

    We can also create abelian Lie algebras using the ``abelian`` keyword::

        sage: L.<x,y,z> = LieAlgebra(QQ, abelian=True); L
        Abelian Lie algebra on 3 generators (x, y, z) over Rational Field

    We can also input a set of structure coefficients. For example, we want
    to create the Lie algebra of `\QQ^3` under the Lie bracket of `\times`
    (cross-product)::

        sage: d = {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}}
        sage: L.<x,y,z> = LieAlgebra(QQ, d)
        sage: L
        Lie algebra on 3 generators (x, y, z) over Rational Field

    To compute the Lie backet of two objects, you cannot use the ``*``.
    This will automatically lift up to the universal enveloping algebra.
    To get elements in the Lie algebra, you must use :meth:`bracket`::

        sage: L = LieAlgebra(QQ, {('e','h'): {'e':-2}, ('f','h'): {'f':2},
        ....:                     ('e','f'): {'h':1}}, names='e,f,h')
        sage: e,f,h = L.lie_algebra_generators()
        sage: L.bracket(h, e)
        2*e
        sage: elt = h*e; elt
        e*h + 2*e
        sage: elt.parent()
        Noncommutative Multivariate Polynomial Ring in e, f, h over Rational Field,
         nc-relations: {f*e: e*f - h, h*f: f*h - 2*f, h*e: e*h + 2*e}

    For convienence, there is are two shorthand notations for computing
    Lie backets::

        sage: L([h,e])
        2*e
        sage: L([h,[e,f]])
        0
        sage: L([[h,e],[e,f]])
        -4*e
        sage: L[h, e]
        2*e
        sage: L[h, L[e, f]]
        0

    .. WARNING::

        Because this is a modified (abused) version of python syntax, it
        does **NOT** work with addition. For example ``L([e + [h, f], h])``
        or ``L[e + [h, f], h]`` will both raise errors. Instead you must
        use ``L[e + L[h, f], h]``.

    Now we construct a free Lie algebra in a few different ways. There are
    two primary representations, as brackets and as polynomials::

        sage: L = LieAlgebra(QQ, 'x,y,z'); L # not tested #16823
        Free Lie algebra generated by (x, y, z) over Rational Field
        sage: P.<a,b,c> = LieAlgebra(QQ, representation="polynomial"); P
        Lie algebra generated by (a, c, b) in
         Free Algebra on 3 generators (a, b, c) over Rational Field

    We currently (:trac:`16823`) have the free Lie algebra given in the
    polynomial representation, which is the Lie algebra of the Free algebra.
    So the generators of the free Lie algebra are the generators of the
    free algebra and the Lie bracket is the commutator::

        sage: P.bracket(a, b) + P.bracket(a - c, b + 3*c)
        2*a*b + 3*a*c - 2*b*a + b*c - 3*c*a - c*b

    REFERENCES:

    .. [deGraaf] Willem A. de Graaf. *Lie Algebras: Theory and Algorithms*.
       North-Holland Mathemtaical Library. (2000). Elsevier Science B.V.

    - Victor Kac. *Infinite Dimensional Lie Algebras*.

    - :wikipedia:`Lie_algebra`
    """
    # This works because it is an abstract base class and this
    #    __classcall_private__ will only be called when calling LieAlgebra
    @staticmethod
    def __classcall_private__(cls, R=None, arg0=None, arg1=None, names=None,
                              index_set=None, abelian=False, **kwds):
        """
        Select the correct parent based upon input.

        TESTS::

            sage: LieAlgebra(QQ, abelian=True, names='x,y,z')
            Abelian Lie algebra on 3 generators (x, y, z) over Rational Field
            sage: LieAlgebra(QQ, {('e','h'): {'e':-2}, ('f','h'): {'f':2},
            ....:                 ('e','f'): {'h':1}}, names='e,f,h')
            Lie algebra on 3 generators (e, f, h) over Rational Field
        """
        # Parse associative algebra input
        # -----

        assoc = kwds.get("associative", None)
        if assoc is not None:
            return LieAlgebraFromAssociative(assoc, names=names, index_set=index_set)

        # Parse the remaining arguments
        # -----

        if R is None:
            raise ValueError("invalid arguments")

        check_assoc = lambda A: (isinstance(A, (Ring, MatrixSpace))
                                 or A in Rings()
                                 or A in Algebras(R).Associative())
        if arg0 in ZZ or check_assoc(arg1):
            # Check if we need to swap the arguments
            arg0, arg1 = arg1, arg0

        # Parse the first argument
        # -----

        if isinstance(arg0, dict):
            if not arg0:
                from sage.algebras.lie_algebras.structure_coefficients import AbelianLieAlgebra
                return AbelianLieAlgebra(R, names, index_set)
            elif isinstance(arg0.keys()[0], (list,tuple)):
                # We assume it is some structure coefficients
                arg1, arg0 = arg0, arg1

        if isinstance(arg0, (list, tuple)):
            if all(isinstance(x, str) for x in arg0):
                # If they are all strings, then it is a list of variables
                names = tuple(arg0)

        if isinstance(arg0, str):
            names = tuple(arg0.split(','))
        elif isinstance(names, str):
            names = tuple(names.split(','))

        # Parse the second argument

        if isinstance(arg1, dict):
            # Assume it is some structure coefficients
            from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
            return LieAlgebraWithStructureCoefficients(R, arg1, names, index_set, **kwds)

        # Otherwise it must be either a free or abelian Lie algebra

        if arg1 in ZZ:
            if isinstance(arg0, str):
                names = arg0
            if names is None:
                index_set = range(arg1)
            else:
                if isinstance(names, str):
                    names = tuple(names.split(','))
                    if arg1 != 1 and len(names) == 1:
                        names = tuple('{}{}'.format(names[0],i) for i in xrange(arg1))
                if arg1 != len(names):
                    raise ValueError("the number of names must equal the"
                                     " number of generators")

        if abelian:
            from sage.algebras.lie_algebras.structure_coefficients import AbelianLieAlgebra
            return AbelianLieAlgebra(R, names, index_set)

        # Otherwise it is the free Lie algebra
        rep = kwds.get("representation", "bracket")
        if rep == "polynomial":
            # Construct the free Lie algebra from polynomials in the
            #   free (associative unital) algebra
            # TODO: Change this to accept an index set once FreeAlgebra accepts one
            F = FreeAlgebra(R, names)
            return LieAlgebraFromAssociative(F, F.gens(), names=names, index_set=index_set)

        raise NotImplementedError("the free Lie algebra has only been implemented using polynomials in the free algebra, see trac ticket #16823")

    @staticmethod
    def _standardize_names_index_set(names=None, index_set=None, ngens=None):
        """
        Standardize the ``names`` and ``index_set`` for a Lie algebra.

        .. TODO::

            This function could likely be generalized for any parent
            inheriting from :class:`IndexedGenerators` and (potentially)
            having ``names``. Should this method be moved to
            :class:`IndexedGenerators`?
        """
        if index_set is None:
            if names is None:
                raise ValueError("either the names of the generators"
                                 " or the index set must be specified")
            # If only the names are specified, then we make the indexing set
            #   be the names
            index_set = tuple(names)

        if isinstance(names, str):
            names = tuple(names.split(','))
        elif names is not None:
            names = tuple(names)

        if isinstance(index_set, (tuple, list)):
            index_set = FiniteEnumeratedSet(index_set)

        if names is not None:
            if index_set is None:
                index_set = names
            elif len(names) != index_set.cardinality():
                raise ValueError("the number of names must equal"
                                 " the size of the indexing set")
            if ngens is not None and len(names) != ngens:
                raise ValueError("the number of names must equal the number of generators")
        elif ngens is not None and index_set.cardinality() != ngens:
            raise ValueError("the size of the indexing set must equal"
                             " the number of generators")

        return names, index_set

    # TODO: Should this inherit from IndexedGenerators or should this
    #   be a subclass?
    def __init__(self, R, names=None, index_set=None, category=None):
        """
        The Lie algebra.

        INPUT:

        - ``R`` -- the base ring

        - ``names`` -- (optional) the names of the generators

        - ``index_set`` -- (optional) the indexing set

        - ``category`` -- the category of the Lie algebra; the default is the
          category of Lie algebras over ``R``

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.category()
            Category of finite dimensional lie algebras with basis over Rational Field
        """
        category = LieAlgebras(R).or_subcategory(category)

        self._indices = index_set
        Parent.__init__(self, base=R, names=names, category=category)

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation="polynomial")
            sage: elt = L([x, y]); elt
            x*y - y*x
            sage: elt.parent() is L
            True
        """
        if isinstance(x, list) and len(x) == 2:
            return self(x[0])._bracket_(self(x[1]))

        if x in self.base_ring():
            if x != 0:
                raise ValueError("can only convert the scalar 0 into a Lie algebra element")
            return self.zero()

        return self.element_class(self, x)

    def __getitem__(self, x):
        """
        If `x` is a pair `(a, b)`, return the Lie bracket `(a, b)`. Otherwise
        try to return the `x`-th element of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation="polynomial")
            sage: L[x, [y, x]]
            -x^2*y + 2*x*y*x - y*x^2
        """
        if isinstance(x, tuple) and len(x) == 2:
            return self(x[0])._bracket_(self(x[1]))
        return super(LieAlgebra, self).__getitem__(x)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion from ``R`` into ``self`` and
        ``False`` otherwise.

        The things that coerce into ``self`` are:

        - Lie algebras in the same variables over a base with a coercion
          map into ``self.base_ring()``.

        - A module which coerces into the base vector space of ``self``.

        TESTS::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L._coerce_map_from_(L.free_module())
            True
            sage: L._coerce_map_from_(FreeModule(ZZ, 2))
            True
        """
        if not isinstance(R, LieAlgebra):
            # Should be moved to LieAlgebrasWithBasis somehow since it is a generic coercion
            if self.free_module is not NotImplemented:
                return self.free_module().has_coerce_map_from(R)
            return False

        # We check if it is a subalgebra of something that can coerce into ``self``
        #from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
        #if isinstance(R, LieSubalgebra) and self.has_coerce_map_from(R._ambient):
        #    return R.ambient_lift

        # Lie algebras in the same indices over any base that coerces in
        if R._indices != self._indices:
            return False

        return self.base_ring().has_coerce_map_from(R.base_ring())

    @cached_method
    def zero(self):
        """
        Return the element `0`.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation="polynomial")
            sage: L.zero()
            0
        """
        return self.element_class(self, {})

    # TODO: Find a better place for this method?
    # TODO: Use IndexedGenerators?
    def indices(self):
        """
        Return the indices of the basis of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation="polynomial")
            sage: L.indices()
            {'x', 'y'}
        """
        return self._indices

    # The following methods should belong to ModulesWithBasis?
    def _from_dict(self, d, coerce=False, remove_zeros=True):
        """
        Construct an element of ``self`` from an ``{index: coefficient}``
        dictionary.

        INPUT:

        - ``d`` -- a dictionary ``{index: coeff}`` where each ``index`` is the
          index of a basis element and each ``coeff`` belongs to the
          coefficient ring ``self.base_ring()``

        - ``coerce`` -- a boolean (default: ``False``), whether to coerce the
          ``coeff``s to the coefficient ring

        - ``remove_zeros`` -- a boolean (default: ``True``), if some
          ``coeff``s may be zero and should therefore be removed

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: d = {('p', 1): 4, ('q', 3): 1/2, 'z': -2}
            sage: L._from_dict(d)
            -2*z + 4*p1 + 1/2*q3
        """
        assert isinstance(d, dict)
        if coerce:
            R = self.base_ring()
            d = dict((key, R(coeff)) for key,coeff in d.iteritems())
        if remove_zeros:
            d = dict((key, coeff) for key, coeff in d.iteritems() if coeff)
        return self.element_class(self, d)

    def monomial(self, i):
        """
        Return the monomial indexed by ``i``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.monomial(('p', 1))
            p1
        """
        return self.element_class(self, {i: self.base_ring().one()})

    def term(self, i, c=None):
        """
        Return the term indexed by ``i`` with coefficient ``c``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.term(('p', 1), 4)
            4*p1
        """
        if c is None:
            c = self.base_ring().one()
        else:
            c = self.base_ring()(c)
        return self.element_class(self, {i: c})

    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation="polynomial")
            sage: L.lie_algebra_generators()
            Finite family {'y': y, 'x': x}
        """
        return Family(self._indices, self.monomial, name="monomial map")

    Element = LieAlgebraElement # Default for all Lie algebras

class FinitelyGeneratedLieAlgebra(LieAlgebra):
    """
    An fintely generated Lie algebra.

    INPUT:

    - ``R`` -- the base ring

    - ``names`` -- the names of the generators

    - ``index_set`` -- the index set of the generators

    - ``category`` -- the category of the Lie algebra
    """
    def __init__(self, R, names=None, index_set=None, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.category()
            Category of finite dimensional lie algebras with basis over Rational Field
        """
        LieAlgebra.__init__(self, R, names, index_set, category)
        self.__ngens = len(self._indices)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F.<x,y> = LieAlgebra(QQ, {('x','y'): {'x': 1}})
            sage: F
            Lie algebra on 2 generators (x, y) over Rational Field
        """
        if self.__ngens == 1:
            return "Lie algebra on the generator {0} over {1}".format(
                    self.gens()[0], self.base_ring())
        return "Lie algebra on {0} generators {1} over {2}".format(
                self.__ngens, self.gens(), self.base_ring())

    @cached_method
    def gens(self):
        """
        Return a tuple whose entries are the generators for this
        object, in some order.

        EXAMPLES::
    
            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.gens()
            (x, y)
        """
        G = self.lie_algebra_generators()
        try:
            return tuple(G[i] for i in self.variable_names())
        except (KeyError, IndexError):
            return tuple(G[i] for i in self.indices())
        except (KeyError, ValueError):
            return tuple(G)

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.gen(0)
            x
        """
        return self.gens()[i]

    @lazy_attribute
    def _ordered_indices(self):
        """
        Return the index set of ``self`` in order.

        EXAMPLES::
    
            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L._ordered_indices
            ('x', 'y')
        """
        return tuple(self._indices)

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::
    
            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.an_element()
            x + y
        """
        return self.sum(self.lie_algebra_generators())

class InfinitelyGeneratedLieAlgebra(LieAlgebra):
    r"""
    An infinitely generated Lie algebra.
    """
    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L._an_element_()
            z + p2 + q2 - 1/2*q3
        """
        return self.lie_algebra_generators()[self._indices.an_element()]

    def dimension(self):
        r"""
        Return the dimension of ``self``, which is `\infty`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.dimension()
            +Infinity
        """
        return infinity

class LieAlgebraFromAssociative(FinitelyGeneratedLieAlgebra):
    """
    A Lie algebra whose elements are from an associative algebra and whose
    bracket is the commutator.

    .. WARNING::

        The returned universal enveloping algebra is too large in general.
        To fix this, we need subalgebras implemented.

    .. TODO::

        Return the subalgebra generated by the basis
        elements of ``self`` for the universal enveloping algebra.

    EXAMPLES:

    We create a simple example with a commutative algebra as the base algebra.
    Note that the bracket of everything will be 0::

        sage: R.<a,b> = PolynomialRing(QQ)
        sage: L = LieAlgebra(associative=R)
        sage: L.bracket(x, y)
        0

    Next we use a free algebra and do some simple computations::

        sage: R.<a,b> = FreeAlgebra(QQ, 2)
        sage: L.<x,y> = LieAlgebra(associative=R)
        sage: x-y
        a - b
        sage: L.bracket(x-y, x)
        a*b - b*a
        sage: L.bracket(x-y, L.bracket(x,y))
        a^2*b - 2*a*b*a + a*b^2 + b*a^2 - 2*b*a*b + b^2*a

    We can also use a subset of the generators to use in our Lie algebra::

        sage: R.<a,b,c> = FreeAlgebra(QQ, 3)
        sage: L.<x,y> = LieAlgebra(associative=[a,b])

    Now for a more complicated example using the group ring of `S_3` as our
    base algebra::

        sage: G = SymmetricGroup(3)
        sage: S = GroupAlgebra(G, QQ)
        sage: L.<x,y> = LieAlgebra(associative=S.gens())
        sage: L.bracket(x, y)
        (2,3) - (1,3)
        sage: L.bracket(x, y-x)
        (2,3) - (1,3)
        sage: L.bracket(L.bracket(x, y), y)
        2*(1,2,3) - 2*(1,3,2)
        sage: L.bracket(x, L.bracket(x, y))
        (2,3) - 2*(1,2) + (1,3)
        sage: L.bracket(x, L.bracket(L.bracket(x, y), y))
        0

    Here is an example using matrices::

        sage: MS = MatrixSpace(QQ,2)
        sage: m1 = MS([[0, -1], [1, 0]])
        sage: m2 = MS([[-1, 4], [3, 2]])
        sage: L.<x,y> = LieAlgebra(associative=[m1, m2])
        sage: x
        [ 0 -1]
        [ 1  0]
        sage: y
        [-1  4]
        [ 3  2]
        sage: L.bracket(x,y)
        [-7 -3]
        [-3  7]
        sage: L.bracket(y,y)
        [0 0]
        [0 0]
        sage: L.bracket(y,x)
        [ 7  3]
        [ 3 -7]
        sage: L.bracket(x, L.bracket(y,x))
        [-6 14]
        [14  6]

    If we use a subset of the generators to construct our Lie algebra,
    the result of :meth:`universal_enveloping_algebra()` can be too large::

        sage: R.<a,b,c> = FreeAlgebra(QQ, 3)
        sage: L = LieAlgebra(associative=[a,b], names='x,y')
        sage: L.universal_enveloping_algebra() is R
        True
    """
    @staticmethod
    def __classcall_private__(cls, A, gens=None, names=None, index_set=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L1 = LieAlgebra(associative=S.gens(), names='x,y')
            sage: L2 = LieAlgebra(associative=[ S(G((1,2,3))), S(G((1,2))) ], names='x,y')
            sage: L1 is L2
            True

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: L1 = LieAlgebra(associative=F)
            sage: L2.<x,y,z> = LieAlgebra(associative=F.gens())
            sage: L3.<x,y,z> = LieAlgebra(QQ, representation="polynomial")
            sage: L1 is L2
            True
        """
        # If A is not a ring, then we treat it as a set of generators
        if isinstance(A, Parent) and A.category().is_subcategory(Rings()):
            # TODO: As part of #16823, this should instead construct a
            #   subclass with specialized methods for the free Lie algebra
            if is_FreeAlgebra(A):
                if gens is None:
                    gens = A.gens()
                if index_set is None:
                    index_set = A.variable_names()
            elif gens is None:
                # Parse possible generators (i.e., a basis) from the parent
                try:
                    gens = A.basis()
                except (AttributeError, NotImplementedError):
                    pass
            if names is None:
                try:
                    names = A.variable_names()
                except ValueError:
                    pass
        else:
            gens = A
            A = None

        if index_set is None:
            # See if we can get an index set from the generators
            try:
                index_set = gens.keys()
            except (AttributeError, ValueError):
                pass
        try:
            ngens = len(gens)
        except (TypeError, NotImplementedError):
            ngens = None

        names, index_set = LieAlgebra._standardize_names_index_set(names, index_set, ngens)

        is_fg = False # Is finitely generated
        if isinstance(gens, AbstractFamily):
            if gens.cardinality() < float('inf'):
                is_fg = True
                try:
                    gens = tuple([gens[i] for i in index_set])
                except KeyError:
                    gens = tuple(gens)

        elif isinstance(gens, dict):
            is_fg = True
            gens = gens.values()
        elif gens: # Assume it is list-like
            is_fg = True
            gens = tuple(gens)

        if is_fg: # If we have a finite generating set
            if A is None:
                A = gens[0].parent()
            # Make sure all the generators have the same parent of A
            gens = tuple(map(A, gens))

        if ngens:
            try:
                # Try to make things, such as matrices, immutable
                #    since we need to hash them
                for g in gens:
                    g.set_immutable()
            except AttributeError:
                pass

        return super(LieAlgebraFromAssociative, cls).__classcall__(cls,
                     A, gens, names=names, index_set=index_set)

    def __init__(self, A, gens, names=None, index_set=None, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: TestSuite(L).run()
        """
        self._assoc = A
        R = self._assoc.base_ring()

        # We strip the following axioms from the category of the assoc. algebra:
        #   FiniteDimensional and WithBasis
        category = LieAlgebras(R).or_subcategory(category)
        if 'FiniteDimensional' in self._assoc.category().axioms():
            category = category.FiniteDimensional()
        if 'WithBasis' in self._assoc.category().axioms():
            category = category.WithBasis()

        FinitelyGeneratedLieAlgebra.__init__(self, R, names, index_set, category)

        if isinstance(gens, tuple):
            gens = Family({self._indices[i]: self.element_class(self, v)
                           for i,v in enumerate(gens)})
        elif gens is not None: # It is a family
            gens = Family(self._indices, lambda i: self.element_class(self, gens[i]),
                          name="generator map")
        self._gens = gens
        # We don't need to store the original generators because we can
        #   get them from lifting this object's generators

        # We construct the lift map in order to register the coercion
        #   here since the UEA already exists
        self.lift

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: MS = MatrixSpace(QQ,2)
            sage: L.<x> = LieAlgebra(associative=[MS.one()])
            sage: L._repr_option('element_ascii_art')
            True
        """
        return self._assoc._repr_option(key)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: LieAlgebra(associative=S)
            Lie algebra generated by ((2,3), (1,3), (1,2), (1,3,2), (), (1,2,3))
             in Group algebra of group
             "Symmetric group of order 3! as a permutation group"
             over base ring Rational Field
            sage: LieAlgebra(associative=S.gens())
            Lie algebra generated by ((1,2,3), (1,2))
             in Group algebra of group
             "Symmetric group of order 3! as a permutation group"
             over base ring Rational Field
        """
        if self._gens is not None:
            return "Lie algebra generated by {} in {}".format(tuple(self._gens), self._assoc)
        return "Lie algebra of {}".format(self._assoc)

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: x,y = S.algebra_generators()
            sage: elt = L(x - y); elt
            -(1,2) + (1,2,3)
            sage: elt.parent() is L
            True
            sage: elt == L(x) - L(y)
            True
            sage: L([x, y])
            (2,3) - (1,3)
        """
        if isinstance(x, list) and len(x) == 2:
            return self(x[0])._bracket_(self(x[1]))
        return self.element_class(self, self._assoc(x))

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: L._construct_UEA() is S
            True
        """
        return self._assoc

    def lie_algebra_generators(self):
        """
        Return the Lie algebra generators of ``self``.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: L.lie_algebra_generators()
            Finite family {(2,3): (2,3), (1,3): (1,3), (1,2,3): (1,2,3),
                           (): (), (1,2): (1,2), (1,3,2): (1,3,2)}
        """
        return self._gens

    def monomial(self, i):
        """
        Return the monomial indexed by ``i``.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(QQ)
            sage: L = LieAlgebra(associative=F)
            sage: L.monomial('x')
            x
        """
        if i not in self._indices:
            #return self(self._assoc.monomial(i))
            raise ValueError("not an index")
        return self._gens[i]

    def term(self, i, c=None):
        """
        Return the term indexed by ``i`` with coefficient ``c``.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(QQ)
            sage: L = LieAlgebra(associative=F)
            sage: L.term('x', 4)
            4*x
        """
        if i not in self._indices:
            #return self(self._assoc.term(i, c))
            raise ValueError("not an index")
        return c * self._gens[i]

    @cached_method
    def zero(self):
        """
        Return `0`.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: L.zero()
            0
        """
        return self.element_class(self, self._assoc.zero())

    def is_abelian(self):
        """
        Return ``True`` if ``self`` is abelian.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 2, 'x,y')
            sage: L = LieAlgebra(associative=R)
            sage: L.is_abelian()
            False

            sage: R = PolynomialRing(QQ, 'x,y')
            sage: L = LieAlgebra(associative=R)
            sage: L.is_abelian()
            True

        An example with a Lie algebra from the group algebra::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L = LieAlgebra(associative=S)
            sage: L.is_abelian()
            False

        Now we construct a Lie algebra from commuting elements in the group
        algebra::

            sage: G = SymmetricGroup(5)
            sage: S = GroupAlgebra(G, QQ)
            sage: gens = map(S, [G((1, 2)), G((3, 4))])
            sage: L.<x,y> = LieAlgebra(associative=gens)
            sage: L.is_abelian()
            True
        """
        if self._assoc.is_commutative():
            return True
        return super(LieAlgebraFromAssociative, self).is_abelian()

    class Element(LieAlgebraElementWrapper):
        def _bracket_(self, rhs):
            """
            Return the bracket ``[self, rhs]``.

            EXAMPLES::

                sage: L.<x,y,z> = LieAlgebra(QQ, representation="polynomial")
                sage: L.bracket(x, y)
                x*y - y*x

                sage: G = SymmetricGroup(3)
                sage: S = GroupAlgebra(G, QQ)
                sage: L.<x,y> = LieAlgebra(associative=S.gens())
                sage: L.bracket(x, y)
                (2,3) - (1,3)

                sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
                sage: L.bracket(L.gen(0), L.gen(1))
                [ 1  0]
                [ 0 -1]
            """
            ret = self.value * rhs.value - rhs.value * self.value
            return self.__class__(self.parent(), ret)

        def lift(self):
            """
            Lift ``self`` to the universal enveloping algebra.

            EXAMPLES::

                sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
                sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
                sage: x.lift()
                x
                sage: x.lift().parent()
                Free Algebra on 3 generators (x, y, z) over Rational Field
            """
            return self.value

