"""
This module contains tools for computing with quiver representations.

AUTHOR:

- Jim Stark (2012-03-04): Initial implementation of acyclic quivers without
  relations.
- Simon King (2013-05, 2014-02): Split code up. Allow cyclic quivers where
  possible.

A Quiver is a directed graph used for representation theory. In our
representation theoretic code, it is assumed that

- the vertices of the quiver are labelled by integers, and

- each edge of the quiver is labelled with a nonempty string.  The label cannot
  begin with 'e_' or contain '*' and distinct edges must have distinct labels.

As far as the :class:`~sage.graphs.digraph.DiGraph` class is concerned, a path
is a finite list of vertices `v_1, ..., v_n` such that there exists an edge
from `v_i` to `v_{i + 1}`.  If there are multiple edges between the same two
vertices this does not contribute additional paths as listed by the DiGraph
class; for example only two paths are listed from 1 to 3 in Q::

    sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
    sage: Q.edges()
    [(1, 2, 'a'), (1, 2, 'b'), (1, 3, 'c'), (2, 3, 'd')]
    sage: Q.all_paths(1, 3)
    [[1, 2, 3], [1, 3]]

When listing paths in a Quiver in representation theory, it is of theoretical
importance to distinguish parallel edges between the same two vertices of a
Quiver.  Specifically we say a path is given by two vertices, ``start`` and
``end``, and a finite (possibly empty) list of edges `e_1, e_2, ..., e_n`
such that the initial vertex of `e_1` is ``start``, the final vertex of `e_i`
is the initial vertex of `e_{i + 1}`, and the final vertex of `e_n` is
``end``.  In the case where no edges are specified, we must have
``start = end`` and the path is called the trivial path at the given vertex.

Quiver paths in the sense stated above correspond to the elements of a
partial semigroup, with multiplication of paths given by concatenation. Hence,
rather than overloading the method name inherited from DiGraph or inventing a
new method name, we move this functionality to this so-called *path
semigroup*.  Note that with this definition there are three paths from 1 to 3
in our example::

    sage: Q.path_semigroup().all_paths(1, 3)
    [a*d, b*d, c]

The returned paths are of type :class:`~sage.quivers.paths.QuiverPath`, which
are elements in the path semigroup that is associated with the quiver. You can
specify a QuiverPath by giving an edge or a list of edges, passed as arguments
to the path semigroup containing this path.  Here an edge is a tuple of
the form ``(i, j, l)``, where ``i`` and ``j`` are vertices and ``l`` is the
label of an edge from i to j::

    sage: p = Q.path_semigroup()([(1, 2, 'a'), (2, 3, 'd')])
    sage: p
    a*d

Trivial paths are indicated by passing the tuple ``(vertex, vertex)``::

    sage: Q.path_semigroup()((6, 6))
    e_6

Trivial "edges" can occur in the input.  They are simply deleted if their
vertex matches the start and end vertex of adjacent edges. Here is an
alternative way to define a path::

    sage: PQ = Q.path_semigroup()
    sage: q = PQ([(1, 1), (1, 2, 'a'), (2, 2), (2, 3, 'd'), (3, 3)])
    sage: p == q
    True

If the vertex of a trivial path does not match with adjacent edges, or if two
adjacent edges do not match, an error is raised.

::

    sage: inv1 = PQ([(1, 2, 'a'), (1, 1)])
    Traceback (most recent call last):
    ...
    ValueError: Cannot interpret [(1, 2, 'a'), (1, 1)] as element of
    Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices
    sage: inv2 = PQ([(1, 2, 'a'), (1, 2, 'a')])
    Traceback (most recent call last):
    ...
    ValueError: Cannot interpret [(1, 2, 'a'), (1, 2, 'a')] as element of
    Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices
    sage: inv3 = PQ([(1, 2, 'x')])
    Traceback (most recent call last):
    ...
    ValueError: Cannot interpret [(1, 2, 'x')] as element of
    Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices

The `*` operator is concatenation of paths. If the two paths do not compose,
then the result is ``None`` (whence the "partial" in "partial semigroup").  ::

    sage: print p*q
    None

Let us now construct a larger quiver::

    sage: Qbig = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}, 3:{4:['e']}, 4:{5:['f']}, 5:{1:['g']} })
    sage: Pbig = Qbig.path_semigroup()

Since ``Q`` is a sub-digraph of ``Qbig``, we have a coercion of the associated
path semigroups::

    sage: Pbig.has_coerce_map_from(PQ)
    True

In particular, `p` is considered to be an element of ``Pbig``, and can be
composed with paths that were defined for the larger quiver::

    sage: p in Pbig
    True
    sage: p*Pbig([(3, 4, 'e')])
    a*d*e
    sage: Pbig([(4, 5, 'f'), (5, 1, 'g')])*p
    f*g*a*d

The length of a path is the number of edges in that path::

    sage: len(p)
    2
    sage: triv = PQ((1, 1))
    sage: len(triv)
    0

List index and slice notation can be used to access the edges in a path.
QuiverPaths can also be iterated over.  Trivial paths have no elements::

    sage: for x in p: print x
    (1, 2, 'a')
    (2, 3, 'd')
    sage: triv[:]
    []

There are methods giving the initial and terminal vertex of a path::

    sage: p.initial_vertex()
    1
    sage: p.terminal_vertex()
    3

QuiverPaths form the basis of the quiver algebra of a quiver.  Given a
field `k` and a Quiver `Q`, the quiver algebra `kQ` is, as a vector space,
the free `k`-vector space whose basis is the set of all paths in `Q`.
Multiplication is defined on this basis and extended bilinearly.  The
product of two basis elements is given by path composition when it
makes sense and is set to be zero otherwise.  Specifically, if the
terminal vertex of the left path equals the initial vertex of the right
path, then their product is the concatenation of the two paths, and
otherwise their product is zero. In sage, quiver algebras
are handled by the :class:`QuiverAlgebra` class::

    sage: A = PQ.algebra(GF(7))
    sage: A
    Path algebra of Multi-digraph on 3 vertices over Finite Field of size 7

Quivers have a method that creates their algebra over a given field (or,
more generally, commutative ring).  Note that QuiverAlgebras are uniquely
defined by their Quiver and field, and play nicely with coercions of the
underlying path semigroups::

    sage: A is PQ.algebra(GF(7))
    True
    sage: A is PQ.algebra(RR)
    False
    sage: Q1 = Q.copy()
    sage: Q1.add_vertex(4)
    sage: PQ1 = Q1.path_semigroup()
    sage: A is PQ1.algebra(GF(7))
    False
    sage: Pbig.algebra(GF(7)).has_coerce_map_from(A)
    True

The QuiverAlgebra can create elements from QuiverPaths or from elements of the
base ring::

    sage: A(5)
    5*e_1 + 5*e_2 + 5*e_3
    sage: r = PQ([(1, 2, 'b'), (2, 3, 'd')])
    sage: e2 = PQ((2, 2))
    sage: x = A(p) + A(e2)
    sage: x
    e_2 + a*d
    sage: y = A(p) + A(r)
    sage: y
    a*d + b*d

QuiverAlgebras are `\NN`-graded algebras.  The grading is given by
assigning to each basis element the length of the path corresponding to
that basis element::

    sage: x.is_homogeneous()
    False
    sage: x.degree()
    Traceback (most recent call last):
    ...
    ValueError: Element is not homogeneous.
    sage: y.is_homogeneous()
    True
    sage: y.degree()
    2
    sage: A[1]
    Free module spanned by [a, b, c, d] over Finite Field of size 7
    sage: A[2]
    Free module spanned by [a*d, b*d] over Finite Field of size 7

The category of right modules over a given quiver algebra is equivalent to the
category of representations of that quiver.  A quiver representation is a
diagram in the category of vector spaces whose underlying graph is the quiver.
So to each vertex of the quiver we assign a vector space and to each edge of
the quiver a linear map between the vector spaces assigned to the start and end
vertices of that edge.  To create the zero representation we just specify the
base ring and the path semigroup::

    sage: Z = Q1.path_semigroup().representation(GF(5))
    sage: Z.is_zero()
    True

To each vertex of a Quiver there is associated a simple module, an
indecomposable projective, and an indecomposable injective, and these can
be created from the Quiver::

    sage: S = PQ.S(GF(3), 1)
    sage: I = PQ.I(QQ, 2)
    sage: P = PQ.P(GF(3), 1)

Radicals, socles, tops, and quotients can all be computed and we can test if
modules are simple or semisimple, get their dimension, and test for equality.
Like Quivers, :class:`~sage.quivers.representation.QuiverRep` objects are
unique and therefore equal if and only if they are identical::

    sage: P.is_simple()
    False
    sage: P.dimension()
    6
    sage: R = P.radical()
    sage: P.socle()
    Representation with dimension vector (0, 0, 3)
    sage: (P/R).is_simple()
    True
    sage: P == R
    False
    sage: P.top() is P/R
    True

There are special methods to deal with modules that are given as right ideals
in the quiver algebra.  To create such a module pass the keyword option='paths'
along with a path or list of paths that generate the desired ideal::

    sage: M = PQ.representation(QQ, [[(1, 1)], [(1, 2, 'a')]], option='paths')
    sage: M.dimension_vector()
    (1, 2, 3)

There are also special methods to deal with modules that are given as the
linear dual of a right ideal in the quiver algebra.  To create such a
module, pass the keyword ``option='dual paths'`` to the constructor along
with a path or list of paths.  The module returned is the dual of the
ideal created in the opposite quiver by the reverses of the given paths::

    sage: D = PQ.representation(QQ, [[(1, 1)], [(1, 2, 'a')]], option='dual paths')
    sage: D.dimension_vector()
    (2, 0, 0)

For modules that are not a standard module or an ideal of the quiver algebra
:class:`~sage.quivers.representation.QuiverRep` can take as input two
dictionaries.  The first associates to each vertex a vector space or an
integer (the desired dimension of the vector space), the second associates to
each edge a map or a matrix or something from which sage can construct a map::

    sage: PQ2 = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
    sage: M2 = PQ2.representation(QQ, {1: QQ^2, 2: QQ^1}, {(1, 2, 'a'): [1, 0], (1, 2, 'b'): [0, 1]})
    sage: M.get_space(2)
    Vector space of dimension 2 over Rational Field
    sage: M2.get_space(2)
    Vector space of dimension 1 over Rational Field
    sage: M.get_map((1, 2, 'a'))
    Vector space morphism represented by the matrix:
    [1 0]
    Domain: Vector space of dimension 1 over Rational Field
    Codomain: Vector space of dimension 2 over Rational Field

A homomorphism between two quiver representations is given by homomorphisms
between the spaces assigned to the vertices of those representations such that
those homomorphisms commute with the edge maps of the representations. The
homomorphisms are created in the usual Sage syntax, the defining data given by
a dictionary associating maps to vertices::

    sage: P2 = PQ2.P(QQ, 1)
    sage: f = P2.hom({1:[1, 1], 2:[[1], [1]]}, M2)

When the domain is given as a right ideal in the quiver algebra we can also
create a homomorphism by just giving a single element in the codomain.  The map
is then induced by acting on that element::

    sage: x = P2.gens('x')[0]
    sage: x
    x_0
    sage: f == P2.hom(f(x), M2)
    True

As you can see, the above homomorphisms can be applied to elements.  Just
like elements, addition is defined via the + operator.  On elements scalar
multiplication is defined via the `*` operator but on homomorphisms `*`
defines composition, so scalar multiplication is done using a method::

    sage: g = f + f
    sage: g == f.scalar_mult(2)
    True
    sage: g == 2*f       # This multiplies the map with the scalar 2
    True
    sage: g(x) == 2*f(x) # This applies the map, then multiplies by the scalar
    True

The ``direct_sum`` method for modules returns only the resulting module by
default.  But can also return the projection and inclusion homomorphisms into
the various factors::

    sage: N2, inclusions, projections = M2.direct_sum([P2], return_maps=True)
    sage: inclusions[0].domain() is M2
    True
    sage: projections[0].codomain() is M2
    True
    sage: (projections[0]*inclusions[0]).is_isomorphism()
    True

As you see above we can determine if a given map is an isomorphism.  Testing
for injectivity and surjectivity works as well::

    sage: f.is_injective()
    False
    sage: f.is_surjective()
    False

We can create all the standard modules associated to maps::

    sage: f.kernel()
    Representation with dimension vector (0, 1)
    sage: f.cokernel()
    Representation with dimension vector (1, 0)
    sage: im = f.image()
    sage: im
    Representation with dimension vector (1, 1)

These methods, as well as the ``submodule`` and ``quotient`` methods that are
defined for representations, return only the resulting representation.  To get
the inclusion map of a submodule or the factor homomorphism of a quotient use
``coerce_map_from``::

    sage: incl = M2.coerce_map_from(im)
    sage: incl.domain() is im
    True
    sage: incl.codomain() is M2
    True
    sage: incl.is_injective()
    True

Both :class:`~sage.quivers.representation.QuiverRep` objects and
:class:`~sage.quivers.homspace.QuiverRepHom` objects have ``linear_dual`` and
``algebraic_dual`` methods.  The ``linear_dual`` method applies the functor
`Hom_k(..., k)` where `k` is the base ring of the representation, and the
``algebraic_dual`` method applies the functor `Hom_Q(..., kQ)` where `kQ`
is the quiver algebra.  Both these functors yield left modules.  A left
module is equivalent to a right module over the opposite algebra, and the
opposite of a quiver algebra is the algebra of the opposite quiver, so both
these methods yield modules and representations of the opposite quiver::

    sage: f.linear_dual()
    Homomorphism of representations of Reverse of (): Multi-digraph on 2 vertices
    sage: D = M2.algebraic_dual()
    sage: D.quiver() is PQ2.reverse().quiver()
    True

.. TODO::

    Change the wording ``Reverse of ()`` into something more meaningful.

There is a method returning the projective cover of any module.  Note that this
method returns the homomorphism; to get the module take the domain of the
homomorphism::

    sage: cov = M2.projective_cover()
    sage: cov
    Homomorphism of representations of Multi-digraph on 2 vertices
    sage: cov.domain()
    Representation with dimension vector (2, 4)

As projective covers are computable, so are the transpose and Auslander-Reiten
translates of modules::

    sage: M2.transpose()
    Representation with dimension vector (4, 3)
    sage: PQ2.I(QQ, 1).AR_translate()
    Representation with dimension vector (3, 2)

We have already used the ``gens`` method above to get an element of a quiver
representation.  An element of a quiver representation is simply a choice of
element from each of the spaces assigned to the vertices of the quiver.
Addition, subtraction, and scalar multiplication are performed pointwise and
implemented by the usual operators::

    sage: M2.dimension_vector()
    (2, 1)
    sage: x, y, z = M2.gens('xyz')
    sage: 2*x + y != x + 2*y
    True

To create a specific element of a given representation we just specify the
representation and a dictionary associating to each vertex an element of the
space associated to that vertex in the representation::

    sage: w = M2({1:(1, -1), 2:(3,)})
    sage: w.get_element(1)
    (1, -1)

The right action of a quiver algebra on an element is implemented via the `*`
operator::

    sage: A2 = x.quiver().path_semigroup().algebra(QQ)
    sage: a = A2((1, 2, 'a'))
    sage: x*a == z
    True
"""

#*****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013, 2014 Simon King <simon.king@uni-jena.de>
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
