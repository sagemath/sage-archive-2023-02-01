"""
This module contains classes for Quivers and their representations.

AUTHOR:

- Jim Stark (2012-03-04): Initial implementation of acyclic quivers without
                        relations.
- Simon King (2013-05): Split code up. Allow cyclic quivers where possible.

A Quiver is a directed graph used for representation theory.  In Sage, a Quiver
is different from a directed graph in the following ways:

- The vertices of a :class:`~sage.graphs.digraph.DiGraph` are arbitrary sage
  objects, but the vertices of a Quiver must be labeled by integers.

- The edges of a DiGraph are labeled with arbitrary sage objects or None if no
  label is specified.  Each edge of a Quiver must be labeled with a nonempty
  string.  The label cannot begin with 'e_' or contain '*' and distinct edges
  must have distinct labels.

- DiGraphs do not have a unique representation in Sage; Quivers do.

- DiGraphs are mutable, Quivers are immutable.

Quivers can be described using a dictionary of dictionaries.  The keys of this
dictionary are vertices ``u``.  The value associated to each ``u`` is a dictionary
whose keys are also vertices ``v`` and the value associated to ``v`` is a list of
strings that label the edges from ``u`` to ``v``::

    sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
    sage: Q.edges()
    [(1, 2, 'a'), (1, 2, 'b'), (1, 3, 'c'), (2, 3, 'd')]

Note that once created we cannot modify the Quiver even though methods from the
DiGraph class that do so appear to be available::

    sage: Q.add_vertex(4)
    Traceback (most recent call last):
    ...
    AttributeError: Quivers are immutable.

If we wish to modify the Quiver, for example adding a vertex, then we convert
the Quiver to a DiGraph, add the vertex to the DiGraph, and then construct a
new quiver from the DiGraph::

    sage: Q.vertices()
    [1, 2, 3]
    sage: D = Q.to_directed()
    sage: D.add_vertex(4)
    sage: Q1 = Quiver(D)
    sage: Q1.vertices()
    [1, 2, 3, 4]

Unique representation means that equal Quivers are identical::

    sage: Q2 = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}, 4:{}})
    sage: Q2 == Q1
    True
    sage: Q2 is Q1
    True

A subquiver of a Quiver Q is another Quiver Q1 such that each vertex and edge
of Q1 is also a vertex and edge of Q.  Testing for subquivers is done via <, >,
<=, and >=::

    sage: Q < Q
    False
    sage: Q <= Q
    True
    sage: Q < Q1
    True
    sage: Q > Q1
    False

Some methods, beyond what the DiGraph class implements, have been added.  We
can list the sources and sinks of a Quiver and test for these properties::

    sage: Q.sources()
    [1]
    sage: Q.sinks()
    [3]
    sage: Q.is_source(2)
    False
    sage: Q.is_sink(2)
    False

As far as the :class:`~sage.graphs.digraph.DiGraph` class is concerned, a path
is a finite list of vertices [v_1, ..., v_n] such that there exists an edge
from v_i to v_{i + 1}.  If there are multiple edges between the same two
vertices this does not contribute to additional paths as listed by the DiGraph
class, for example only two paths are listed from 1 to 3 in Q::

    sage: Q.all_paths(1, 3)
    [[1, 2, 3], [1, 3]]

When listing paths in a Quiver, it is of theoretical importance to distinguish
parallel edges between the same two vertices of a Quiver.  Specifically
we say a path is given by two vertices, start and end, and a finite (possibly
empty) list of edges e_1, e_2, ..., e_n such that the initial vertex of e_1 is
start, the final vertex of e_i is the initial vertex of e_{i + 1}, and the
final vertex of e_n is end.  In the case where no edges are specified we must
have start = end and the path is called the trivial path at the given vertex.

Paths in a Quiver are considered to be elements of the free small category (or
semigroup) formed by concatenation of paths. Hence, rather than overloading
the method name inherited from DiGraph or inventing a new method name, we move
this functionality to this so-called free small category.  Note that with this
definition there are three paths from 1 to 3 in our example::

    sage: Q.free_small_category().all_paths(1, 3)
    [b*d, a*d, c]

The all_quiver_paths method returns a list of objects of type
:class:`~sage.quivers.paths.QuiverPath`, which are elements in the free small
category that is associated with the quiver. You can specify a QuiverPath by
giving an edge or a list of edges, passed as arguments to the free small
category containing this path.  Here an edge is a tuple of the form (i, j, l),
where i and j are vertices and l is the label of an edge from i to j::

    sage: p = Q.free_small_category()([(1, 2, 'a'), (2, 3, 'd')])
    sage: p
    a*d

Trivial paths are indicated by passing the tuple (vertex, vertex)::

    sage: Q.free_small_category()((6, 6))
    e_6

Trivial edges can occur in the input.  They are simply deleted if their vertex
matches the start and end vertex of adjacent edges. Here is an alternative way
to define a path::

    sage: F = Q.free_small_category()
    sage: q = F([(1, 1), (1, 2, 'a'), (2, 2), (2, 3, 'd'), (3, 3)])
    sage: p == q
    True

If the vertex of a trivial path does not match with adjacent edges, or if two
adjacent edges do not match, no error is raised.  Instead the invalid path is
returned.  The invalid path can be detected by converting a QuiverPath to a
Boolean.  Valid paths are True, the invalid path is False::

    sage: inv1 = F([(1, 2, 'a'), (1, 1)])
    sage: print inv1
    invalid path
    sage: inv2 = F([(1, 2, 'a'), (1, 2, 'a')])
    sage: inv1 is inv2
    False
    sage: inv1 == inv2
    True
    sage: bool(p)
    True
    sage: bool(inv1)
    False

The * operator is concatenation of paths. If the two paths do not compose, then the
result is the invalid path.
::

    sage: bool(p*q)
    False

Let us now construct a larger quiver::

    sage: Qbig = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}, 3:{4:['e']}, 4:{5:['f']}, 5:{1:['g']} })
    sage: Fbig = Qbig.free_small_category()

Since ``Q`` is a sub-digraph of ``Qbig``, we have a coercion of the associated
free small categories::

    sage: Fbig.has_coerce_map_from(F)
    True

In particular, `p` is considered to be an element of ``Fbig``, and can be
composed with paths that were defined for the larger quiver::

    sage: p in Fbig
    True
    sage: p*Fbig([(3, 4, 'e')])
    a*d*e
    sage: Fbig([(4, 5, 'f'), (5, 1, 'g')])*p
    f*g*a*d

The length of a path is the number of edges in that path.  The invalid path and
trivial paths are therefore length 0::

    sage: len(p)
    2
    sage: triv = F((1, 1))
    sage: len(triv)
    0
    sage: len(inv1)
    0

List index and slice notation can be used to access the edges in a path.
QuiverPaths can also be iterated over.  Trivial paths and the invalid path have
no elements::

    sage: for x in p: print x
    (1, 2, 'a')
    (2, 3, 'd')
    sage: triv[:]
    []
    sage: inv1[0]
    Traceback (most recent call last):
    ...
    IndexError: tuple index out of range

There are methods giving the initial and terminal vertex of a path.  These
return None when called on the invalid path::

    sage: inv1.initial_vertex()
    sage: inv1.terminal_vertex()
    sage: p.initial_vertex()
    1
    sage: p.terminal_vertex()
    3

QuiverPaths form the basis of the quiver algebra of a quiver.  Given a field k
and a Quiver Q the quiver algebra kQ is, as a vector space it has basis the set
of all paths in Q.  Multiplication is defined on this basis and extended
bilinearly.  We multiplication is given as path composition when it makes sense
and is zero otherwise.  Specifically if the terminal vertex of the left path
equals the initial vertex of the right path then their product is the
composition of the two paths, otherwise it is zero. In sage quiver algebras
are handled by the QuiverAlgebra class::

    sage: A = Q.algebra(GF(7))
    sage: A
    Algebra of Quiver on 3 vertices over Finite Field of size 7

Quivers have a method that creates their algebra over a given field.  Note that
QuiverAlgebras are uniquely defined by their Quiver and field, and play nicely
with coercions of the underlying free small categories::

    sage: A is Q.algebra(GF(7))
    True
    sage: A is Q.algebra(RR)
    False
    sage: A is Q1.algebra(GF(7))
    False
    sage: Qbig.algebra(GF(7)).has_coerce_map_from(A)
    True

The QuiverAlgebra can create elements from QuiverPaths or from elements of the
base ring::

    sage: A(5)
    5*e_1 + 5*e_2 + 5*e_3
    sage: r = F([(1, 2, 'b'), (2, 3, 'd')])
    sage: e2 = F((2, 2))
    sage: x = A(p) + A(e2)
    sage: x
    a*d + e_2
    sage: y = A(p) + A(r)
    sage: y
    a*d + b*d

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
    2
    sage: A[1]
    Free module spanned by [b, a, c, d] over Finite Field of size 7
    sage: A[2]
    Free module spanned by [b*d, a*d] over Finite Field of size 7

The category of right modules over a given quiver algebra is equivalent to the
category of representations of that quiver.  A quiver representation is a
diagram in the category of vector spaces whose underlying graph is the quiver.
So to each vertex of the quiver we assign a vector space and to each edge of
the quiver a linear map between the vector spaces assigned to the start and end
vertexes of that edge.  To create the zero representation we just specify the
base ring and Quiver::


    sage: Z = Q1.representation(GF(5))
    sage: Z.is_zero()
    True

To each vertex of a Quiver there is associated a simple module, an
indecomposable projective, and an indecomposable injective and these can be
created from the Quiver::

    sage: S = Q.S(GF(3), 1)
    sage: I = Q.I(QQ, 2)
    sage: P = Q.P(GF(3), 1)

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

    sage: M = Q.representation(QQ, [[(1, 1)], [(1, 2, 'a')]], option='paths')
    sage: M.dimension_vector()
    (1, 2, 3)

There are also special methods to deal with modules that are given as the
linear dual of a right ideal in the quiver algebra.  To create such a module
pass the keyword option='dual paths' to the constructor along with a path or
list of paths.  The module returned is the dual of the ideal created in the
opposite quiver by the reverse of the given paths::

    sage: D = Q.representation(QQ, [[(1, 1)], [(1, 2, 'a')]], option='dual paths')
    sage: D.dimension_vector()
    (2, 0, 0)

For modules that are not a standard module or an ideal of the quiver algebra
:class:`~sage.quivers.representation.QuiverRep` can take as input two
dictionaries.  The first associates to each vertex a vector space or an
integer (the desired dimension of the vector space), the second associates to
each edge a map or a matrix or something from which sage can construct a map::

    sage: Q2 = Quiver({1:{2:['a', 'b']}})
    sage: M2 = Q2.representation(QQ, {1: QQ^2, 2: QQ^1}, {(1, 2, 'a'): [1, 0], (1, 2, 'b'): [0, 1]})
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
a dictionary associating maps to vertexes::

    sage: P2 = Q2.P(QQ, 1)
    sage: f = P2.hom({1:[1, 1], 2:[[1], [1]]}, M2)

When the domain is given as a right ideal in the quiver algebra we can also
create a homomorphism by just giving a single element in the codomain.  The map
is then induced by acting on that element::

    sage: x = P2.gens('x')[0]
    sage: x
    x_0
    sage: f == P2.hom(f(x), M2)
    True

As you can see above homomorphisms can be applied to elements.  Just like
elements, addition is defined via the + operator.  On elements scalar
multiplication is defined via the * operator but on homomorphisms * defines
composition, so scalar multiplication is done using a method::

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
defined for representations return only the resulting representation.  To get
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
`Hom_k(..., k)` where k is the base ring of the representation and the
``algebraic_dual`` method applies the functor `Hom_Q(..., kQ)` where kQ is the
quiver algebra.  Both these functors yeild left modules.  A left module is
equivalent to a right module over the opposite algebra and the opposite of a
quiver algebra is the algebra of the opposite quiver, so both these methods
yeild modules and representations of the opposite quiver::

    sage: f.linear_dual()
    Homomorphism of representations of Quiver on Reverse of ()
    sage: D = M2.algebraic_dual()
    sage: D.quiver() is Q2.reverse()
    True

.. TODO::

    Change the wording ``Reverse of ()`` into something more meaningful.

There is a method returning the projective cover of any module.  Note that this
method returns the homomorphism, to get the module take the domain of the
homomorphism::

    sage: cov = M2.projective_cover()
    sage: cov
    Homomorphism of representations of Quiver on 2 vertices
    sage: cov.domain()
    Representation with dimension vector (2, 4)

As projective covers are computable so are the transpose and Auslander-Reiten
translates of modules::

    sage: M2.transpose()
    Representation with dimension vector (4, 3)
    sage: Q2.I(QQ, 1).AR_translate()
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

The right action of a quiver algebra on an element is implemented via the *
operator::

    sage: A2 = x.quiver().algebra(QQ)
    sage: a = A2((1, 2, 'a'))
    sage: x*a == z
    True
"""

#*****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013 Simon King <simon.king@uni-jena.de>
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

from sage.structure.unique_representation import UniqueRepresentation, CachedRepresentation
from sage.misc.function_mangling import ArgumentFixer
from sage.misc.cachefunc import cached_method
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.free_module import CombinatorialFreeModuleElement
from sage.modules.module import Module
from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.categories.magmas import Magmas
from sage.structure.element import ModuleElement
from sage.categories.morphism import CallMorphism
from sage.graphs.digraph import DiGraph
from sage.rings.integer_ring import ZZ

class Quiver(UniqueRepresentation, DiGraph):
    """
    Generic class for a quiver.  Type Quiver? for more information.

    TESTS::

        sage: Q1 = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
    """
    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################
    @staticmethod
    def __classcall__(cls, *args, **kwds):
        """
        Ensure uniqueness of quivers.

        TESTS::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: Q is Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            True
            sage: Q is Quiver(DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}}))
            True
            sage: Q is loads(dumps(Q))
            True

        """
        cache = UniqueRepresentation.__classcall__.cache
        # For normalisation of the arguments in the cache key,
        # we use an "argument fixer" related with the init method
        # of DiGraph
        Args, KeyKwds = _DigraphArgumentFixer.fix_to_pos(cls, *args, **kwds)
        data = Args[1] # Args[0] is cls
        if isinstance(data,tuple):
            # This would normally occur during unpickling
            try:
                return cache[Args,KeyKwds]
            except KeyError:
                pass
            # It is not in the cache, and thus we need to revert it into a dict,
            # so that the DiGraph constructor can be called.
            tuple_data = data
            data = dict((v, {}) for v in tuple_data[0])
            for e in tuple_data[1:]:
                if e[1] not in data[e[0]]:
                    data[e[0]][e[1]] = []
                data[e[0]][e[1]].append(e[2])
        # We create an instance of cls and initialise it as a digraph.
        # Then we ask this digraph for vertices and edges (ordered).
        # This is used as part of a key for caching.
        Args = (cls, data)+Args[2:]
        Kwds = dict(KeyKwds)
        G = object.__new__(cls,*Args, **Kwds)
        DiGraph.__init__(G,*(Args[1:]), **Kwds)
        key = [tuple(G.vertices())]
        key.extend(G.edges(labels=True, sort=True))
        KeyArgs = (cls, tuple(key))+Args[2:]

        # For caching, we replace the argument "pos" in the
        # above normalised arguments.
        # It could be that after normalisation we are again able to find
        # another instance in the cache.
        try:
            return cache[KeyArgs,KeyKwds]
        except KeyError:
            pass
        # Now, we finish initialisation of G, take care of pickling,
        # put it into the cache, and return it
        cls.__init__(G, *(Args[1:]), **Kwds)
        if G.__class__.__reduce__ == CachedRepresentation.__reduce__:
            G._reduction = (cls, KeyArgs[1:], Kwds)
        cache[KeyArgs, KeyKwds] = G
        return G

    # boundary is (), not [], so that the arguments can be used in caching.
    # In our current applications, quivers have no boundary anyway.
    def __init__(self, data=None, pos=None, loops=True, format=None,
                 boundary=(), weighted=True, implementation='c_graph',
                 sparse=True, vertex_labels=True, name=None,
                 multiedges=True, convert_empty_dict_labels_to_None=False, **kwds):
        """
        Creates a Quiver object.  Type Quiver? for more information.

        TESTS::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
        """
        # DiGraph.__init__ is normally called during __classcall__
        # If it isn't, we do here, but printing a warning
        if not hasattr(self,'_backend'):
            print "WARNING: It seems that you worked around the cache of Quiver."
            DiGraph.__init__(self, data, pos, loops, format, boundary, weighted, implementation, sparse, vertex_labels, name, multiedges, convert_empty_dict_labels_to_None, **kwds)
        self._make_immutable()
        self._assert_labels_valid()

    def _assert_labels_valid(self):
        """
        Check whether all edge labels are pairwise distinct nonempty strings
        that do not start with "e_", and whether  all vertex labels are integers.

        NOTE:

        Admissibility of labels is tested for during initialisation. So, the
        following are indirect tests.

        EXAMPLES::

            sage: Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            Quiver on 3 vertices
            sage: Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['e_']}}) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Edge labels of Representation Quiver cannot begin with 'e_'.
            sage: Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['c']}})
            Traceback (most recent call last):
            ...
            ValueError: Edge labels of Representation Quiver must be unique.
            sage: Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['']}})
            Traceback (most recent call last):
            ...
            ValueError: Edge labels of Quivers must be nonempty strings.
            sage: Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:[2]}})
            Traceback (most recent call last):
            ...
            ValueError: Edge labels of Quivers must be nonempty strings.
            sage: Quiver({1:{2:['a','b'], 3:['c']}, 's':{3:['d']}})
            Traceback (most recent call last):
            ...
            ValueError: Vertices of Representation Quiver must be labeled by integers.

        """
        # Check that edges are labeled with nonempty strings and don't begin
        # with 'e_' or contain '*'
        for x in self.edge_labels():
            if not isinstance(x, str) or x == '':
                raise ValueError("Edge labels of Quivers must be nonempty strings.")
            if x[0:2] == 'e_':
                raise ValueError("Edge labels of Representation Quiver cannot begin with 'e_'.")
            if x.find('*') != -1:
                raise ValueError("Edge labels of Representation Quiver cannot contain '*'.")
        # Check that vertices are labeled 1,2,3,... and that edge labels are
        # unique
        for v in self:
            if v not in ZZ:
                raise ValueError("Vertices of Representation Quiver must be labeled by integers.")
        if len(set(self.edge_labels())) != self.num_edges():
            raise ValueError("Edge labels of Representation Quiver must be unique.")

    def _forbidden_method(self, *args, **kwds):
        """
        Functions that modify the quiver may not be called; quivers are unique.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}})
            sage: Q.add_vertex(3) # indirect doctest
            Traceback (most recent call last):
            ...
            AttributeError: Quivers are immutable.

        To modify a quiver you must use the to_directed method to create a DiGraph,
        modify the DiGraph, and then create a Quiver from the DiGraph::

            sage: G = Q.to_directed()
            sage: G.add_vertex(3)
            sage: Q = Quiver(G)
            sage: Q.vertices()
            [1, 2, 3]
        """
        raise AttributeError("Quivers are immutable.")

    def _make_immutable(self):
        """
        Invalidate certain mutating methods inherited from :class:`~sage.graphs.digraph.DiGraph`.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}})
            sage: Q.add_vertex(3) # indirect doctest
            Traceback (most recent call last):
            ...
            AttributeError: Quivers are immutable.

        Not that attributes of these names still exist, but calling them
        results in an attribute error::

            sage: Q.set_edge_label
            <bound method Quiver._forbidden_method of Quiver on 2 vertices>
            sage: Q.set_edge_label('bla')
            Traceback (most recent call last):
            ...
            AttributeError: Quivers are immutable.

        """
        # Forbid methods that change the vertices or edges
        self.add_cycle = self._forbidden_method
        self.add_edge = self._forbidden_method
        self.add_edges = self._forbidden_method
        self.add_path = self._forbidden_method
        self.add_vertex = self._forbidden_method
        self.add_vertices = self._forbidden_method
        self.allow_loops = self._forbidden_method
        self.allow_multiple_edges = self._forbidden_method
        self.clear = self._forbidden_method
        self.set_edge_label = self._forbidden_method
        self.subdivide_edge = self._forbidden_method
        self.subdivide_edges = self._forbidden_method

    @cached_method
    def is_acyclic_quiver(self):
        """
        Tests whether this quiver has loops or other directed closed paths.

        TESTS::

            sage: acyQ = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: acyQ.is_acyclic_quiver()
            True
            sage: cyQ1 = Quiver({1:{2:['a','b'], 1:['c']}, 2:{3:['d']}})
            sage: cyQ1.is_acyclic_quiver()
            False
            sage: cyQ2 = Quiver({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: cyQ2.is_acyclic_quiver()
            False

        """

        # Check that it's directed, acyclic, and has no loops
        return self.is_directed_acyclic() and not self.has_loops()

    def _repr_(self):
        """
        Default string representation of a quiver.

        TESTS::

            sage: Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}}) # indirect doctest
            Quiver on 3 vertices
        """
        name = self.name()
        if name:
            return "Quiver on "+name
        else:
            return "Quiver on " + str(len(self)) + " vertices"

    def __lt__(self, other):
        """
        Returns true if other is a quiver and the vertices and edges of self are
        vertices and edges of other, but self is not other.  Raises a TypeError if
        other is not a Quiver.

        OUTPUT:

        - bool

        TESTS::

            sage: Q1 = Quiver({1:{2:['a']}, 2:{3:['b']}})
            sage: Q2 = Quiver({1:{2:['a']}})
            sage: Q1 < Q1
            False
            sage: Q1 < Q2
            False
            sage: Q2 < Q1
            True
        """

        # Fail if not a Quiver
        if not isinstance(other, Quiver):
            raise TypeError("Can only compare Quivers with other Quivers")

        # Check for equality
        if self is other:
            return False

        # Check vertices and edges for containment
        for v in self:
            if v not in other:
                return False

        OE = other.edges()
        for e in self.edges():
            if e not in OE:
                return False

        return True

    def __gt__(self, other):
        """
        Returns true if other is a quiver and the vertices and edges of other are
        vertices and edges of self, but self is not other.  Raises a TypeError if
        other is not a Quiver.

        OUTPUT:

        - bool

        TESTS::

            sage: Q1 = Quiver({1:{2:['a']}, 2:{3:['b']}})
            sage: Q2 = Quiver({1:{2:['a']}})
            sage: Q1 > Q1
            False
            sage: Q1 > Q2
            True
            sage: Q2 > Q1
            False
        """

        # Fail if not a Quiver
        if not isinstance(other, Quiver):
            raise TypeError("Can only compare Quivers with other Quivers")

        # Check for equality
        if self is other:
            return False

        # Check vertices and edges for containment
        for v in other:
            if v not in self:
                return False
        SE = self.edges()
        for e in other.edges():
            if e not in SE:
                return False

        return True

    def __le__(self, other):
        """
        Returns true if other is a quiver and the vertices and edges of self are
        vertices and edges of other.  Raises a TypeError if other is not a Quiver.

        OUTPUT:

        - bool

        TESTS::

            sage: Q1 = Quiver({1:{2:['a']}, 2:{3:['b']}})
            sage: Q2 = Quiver({1:{2:['a']}})
            sage: Q1 <= Q1
            True
            sage: Q1 <= Q2
            False
            sage: Q2 <= Q1
            True
        """

        # Fail if not a Quiver
        if not isinstance(other, Quiver):
            raise TypeError("Can only compare Quivers with other Quivers")

        # Check for equality
        if self is other:
            return True

        # Check vertices and edges for containment
        for v in self:
            if v not in other:
                return False

        OE = other.edges()
        for e in self.edges():
            if e not in OE:
                return False

        return True

    def __ge__(self, other):
        """
        Returns true if other is a quiver and the vertices and edges of other are
        vertices and edges of self.  Raises a TypeError if other is not a Quiver.

        OUTPUT:

        - bool

        TESTS::

            sage: Q1 = Quiver({1:{2:['a']}, 2:{3:['b']}})
            sage: Q2 = Quiver({1:{2:['a']}})
            sage: Q1 >= Q1
            True
            sage: Q1 >= Q2
            True
            sage: Q2 >= Q1
            False
        """

        # Fail if not a Quiver
        if not isinstance(other, Quiver):
            raise TypeError("Can only compare Quivers with other Quivers")

        # Check for equality
        if self is other:
            return True

        # Check vertices and edges for containment
        for v in other:
            if v not in self:
                return False
        SE = self.edges()
        for e in other.edges():
            if e not in SE:
                return False

        return True

    ###########################################################################
    #                                                                         #
    # GRAPH THEORETIC FUNCTIONS                                               #
    #    These functions involve the graph theoretic structure of the quiver. #
    #                                                                         #
    ###########################################################################

    def to_directed(self):
        """
        Returns the underlying Digraph.

        OUTPUT:

        - a DiGraph with the same edges and vertices as the quiver

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a']}})
            sage: type(Q.to_directed()) == DiGraph
            True
        """

        return DiGraph(self)

    def is_source(self, vertex):
        """
        Tests whether the vertex is a source in the quiver.

        INPUT:

        - ``vertex`` - integer, the vertex to be tested

        OUTPUT:

        - bool, True if the vertex is not the terminus of any edge, False if there is
          an edge terminating at the vertex

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})

        There are no edges into vertex 1::

            sage: Q.is_source(1)
            True

        but there are edges into vertexes 2 and 3::

            sage: Q.is_source(2)
            False
            sage: Q.is_source(3)
            False
        """

        return len(self.neighbors_in(vertex)) == 0

    def is_sink(self, vertex):
        """
        Tests whether the vertex is a sink in the quiver.

        INPUT:

        - ``vertex`` - integer, the vertex to be tested

        OUTPUT:

        - bool, True if the vertex is not the source of any edge, False if there is an
          edge beginning at the vertex

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})

        There are edges out of vertexes 1 and 2::

            sage: Q.is_sink(1)
            False
            sage: Q.is_sink(2)
            False

        but there are no edges out of vertex 3::

            sage: Q.is_sink(3)
            True
        """

        return len(self.neighbors_out(vertex)) == 0

    def sources(self):
        """
        Returns a list of sources of the quiver.

        OUTPUT:

        - list, the vertices of the quiver that have no edges going into them

        EXAMPLES::

            sage: Q = Quiver({1:{3:['a']}, 2:{3:['b']}})
            sage: Q.sources()
            [1, 2]
            sage: T = Quiver({1:{}})
            sage: T.sources()
            [1]
        """

        return [x for x in self if self.is_source(x)]

    def sinks(self):
        """
        Returns a list of sinks of the quiver.

        OUTPUT:

        - list, the vertices of the quiver that have no edges beginning at them

        EXAMPLES::

            sage: Q = Quiver({1:{3:['a']}, 2:{3:['b']}})
            sage: Q.sinks()
            [3]
            sage: T = Quiver({1:{}})
            sage: T.sinks()
            [1]
        """

        return [x for x in self if self.is_sink(x)]

    def all_paths(self, start=None, end=None):
        """
        Returns a list of all vertex paths between a pair of vertices (start, end).

        INPUT:

        - ``start`` - integer or None (default: None), the initial vertex of the paths in
          the output.  If None is given then the initial vertex is arbitrary.

        - ``end`` - integer or None (default: None), the terminal vertex of the paths in
          the output.  If None is given then the terminal vertex is arbitrary.

        OUTPUT:

        - a list of vertex paths (also lists) in the quiver

        .. NOTE::

            This function differs from all_quiver_paths in that paths are given as
            lists of vertices instead of lists of edges.  If there are multiple edges
            between two vertices all_paths will not differentiate between them but
            all_quiver_paths will.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: Q.all_paths(1, 3)
            [[1, 2, 3], [1, 3]]

        If there are no paths an empty list is returned::

            sage: Q.all_paths(3, 1)
            []

        If start=end then a list containing only the trivial vertex path is returned.
        A trivial vertex path is just a list containing a single vertex::

            sage: Q.all_paths(2, 2)
            [[2]]

        If end=None then all vertex paths begining at start are returned, including
        trivial paths::

            sage: Q.all_paths(2)
            [[2], [2, 3]]

        If start=None then all vertex paths ending at end are returned, including
        trivial paths.  Note that even though there are two edges from vertex 1 to
        vertex 2 there is only one vertex path::

            sage: Q.all_paths(None, 2)
            [[1, 2], [2]]
            sage: Q.all_paths(end=2)
            [[1, 2], [2]]

        If start=end=None then all vertex paths are returned, including trivial paths::

            sage: Q.all_paths()
            [[1], [1, 2], [1, 2, 3], [1, 3], [2], [2, 3], [3]]

        The vertex given must be a vertex of the quiver::

            sage: Q.all_paths(1, 4)
            Traceback (most recent call last):
            ...
            ValueError: The end vertex 4 is not a vertex of the quiver.
        """
        # This function modifies the sage version by allowing start or end=None and not
        # producing an error when start=end.

        # Check that given arguments are vertices
        if start is not None and start not in self:
            raise ValueError("The start vertex " + str(start) + " is not a vertex of the quiver.")
        if end is not None and end not in self:
            raise ValueError("The end vertex " + str(end) + " is not a vertex of the quiver.")

        # Handle start=None
        if start == None:
            results = []
            for v in self:
                results += self.all_paths(v, end)
            return results

        # Handle end=None
        if end == None:
            results = []
            for v in self:
                results += self.all_paths(start, v)
            return results

        # Handle start=end
        if start == end:
            return [[start]]

        # Otherwise call the sage version
        return super(Quiver, self).all_paths(start, end)

    def reverse(self):
        """
        Returns a copy of the quiver with edges reversed in direction.

        OUTPUT:

        - Quiver

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: Qrev = Q.reverse()
            sage: Q.edges()
            [(1, 2, 'a'), (1, 2, 'b')]
            sage: Qrev.edges()
            [(2, 1, 'a'), (2, 1, 'b')]

        Reversing a quiver twice returns the original quiver::

            sage: Qrevrev = Qrev.reverse()
            sage: Qrevrev is Q
            True
        """

        return Quiver(DiGraph.reverse(self))

    ###########################################################################
    #                                                                         #
    # REPRESENTATION THEORETIC FUNCTIONS                                      #
    #    These functions involve the representation theory of quivers.        #
    #                                                                         #
    ###########################################################################

    def free_small_category(self):
        """
        The free small category formed by the paths of this quiver.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','c']}, 2:{3:['b']}})
            sage: F = Q.free_small_category(); F
            Free small category of Quiver on 3 vertices
            sage: list(F)
            [e_1, e_2, e_3, c, a, b, c*b, a*b]

        """
        from free_small_category import FreeSmallCategory
        return FreeSmallCategory(self)

    def representation(self, k, *args, **kwds):
        """
        Returns a representation of the quiver.

        For more information see the QuiverRep documentation.

        TESTS::

            sage: Q = Quiver({1:{3:['a']}, 2:{3:['b']}})
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = Q.representation(QQ, spaces, maps)
        """
        from sage.quivers.representation import QuiverRep
        return QuiverRep(k, self, *args, **kwds)

    def S(self, k, vertex):
        """
        Returns the simple module over k at the given vertex.

        INPUT:

        - ``k`` - ring, the base ring of the representation

        - ``vertex`` - integer, a vertex of the quiver

        OUTPUT:

        - QuiverRep, the simple module at vertex with base ring k

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c','d']}})
            sage: S1 = Q.S(GF(3), 1)
            sage: Q.S(ZZ, 3).dimension_vector()
            (0, 0, 1)
            sage: Q.S(ZZ, 1).dimension_vector()
            (1, 0, 0)

        The vertex given must be a vertex of the quiver::

            sage: Q.S(QQ, 4)
            Traceback (most recent call last):
            ...
            ValueError: Must specify a valid vertex of the quiver.
        """

        # Raise an error if the given vertex is not a vertex
        if vertex not in self:
            raise ValueError("Must specify a valid vertex of the quiver.")

        # This is the module with k at the given vertex and zero elsewhere.  As
        # all maps are zero we only need to specify that the given vertex has
        # dimension 1 and the constructor will zero out everything else.
        from sage.quivers.representation import QuiverRep
        return QuiverRep(k, self, {vertex: 1})

    def P(self, k, vertex):
        """
        Returns the indecomposable projective module over k at the given vertex.

        INPUT:

        - ``k`` - ring, the base ring of the representation

        - ``vertex`` - integer, a vertex of the quiver

        OUTPUT:

        - QuiverRep, the indecomposable projective module at vertex with base ring k

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c','d']}})
            sage: P2 = Q.P(GF(3), 2)
            sage: Q.P(ZZ, 3).dimension_vector()
            (0, 0, 1)
            sage: Q.P(ZZ, 1).dimension_vector()
            (1, 2, 4)

        The vertex given must be a vertex of the quiver::

            sage: Q.P(QQ, 4)
            Traceback (most recent call last):
            ...
            ValueError: Must specify a valid vertex of the quiver.
        """

        # Raise an error if the given vertex is not a vertex
        if vertex not in self:
            raise ValueError("Must specify a valid vertex of the quiver.")
        from sage.quivers.representation import QuiverRep
        return QuiverRep(k, self, [[(vertex, vertex)]], option='paths')

    def I(self, k, vertex):
        """
        Returns the indecomposable injective module over k at the given vertex.

        INPUT:

        - ``k`` - ring, the base ring of the representation

        - ``vertex`` - integer, a vertex of the quiver

        OUTPUT:

        - QuiverRep, the indecomposable injective module at vertex with base ring k

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c','d']}})
            sage: I2 = Q.I(GF(3), 2)
            sage: Q.I(ZZ, 3).dimension_vector()
            (4, 2, 1)
            sage: Q.I(ZZ, 1).dimension_vector()
            (1, 0, 0)

        The vertex given must be a vertex of the quiver::

            sage: Q.I(QQ, 4)
            Traceback (most recent call last):
            ...
            ValueError: Must specify a valid vertex of the quiver.
        """

        # Raise an error if the given vertex is not a vertex
        if vertex not in self:
            raise ValueError("Must specify a valid vertex of the quiver.")
        from sage.quivers.representation import QuiverRep
        return QuiverRep(k, self, [[(vertex, vertex)]], option='dual paths')


    def free_module(self, k):
        """
        Returns a free module of rank 1 over k.

        INPUT:

        - ``k`` - ring, the base ring of the representation.

        OUTPUT:

        - QuiverRep_with_path_basis, the path algebra considered as a right module over
          itself.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a', 'b'], 3: ['c', 'd']}, 2:{3:['e']}})
            sage: Q.free_module(GF(3)).dimension_vector()
            (1, 3, 6)
        """
        from sage.quivers.representation import QuiverRep
        return QuiverRep(k, self, [[(v, v)] for v in self], option='paths')

    def algebra(self, k):
        """
        Return the algebra of the quiver.

        INPUT:

        - ``k`` - ring, the base ring of the quiver algebra

        OUTPUT:

        - QuiverAlgebra

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{4:['c']}, 3:{4:['d']}})
            sage: Q.algebra(GF(7))
            Algebra of Quiver on 4 vertices over Finite Field of size 7
        """
        from sage.quivers.algebra import QuiverAlgebra
        return QuiverAlgebra(k, self)

#####################
# A helper, to normalise the arguments for a quiver
_DigraphArgumentFixer = ArgumentFixer(Quiver.__init__, False)
