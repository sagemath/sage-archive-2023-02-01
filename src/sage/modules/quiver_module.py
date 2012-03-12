"""
This module contains classes for Quivers and their representations.

AUTHOR:

- Jim Stark (2012-3-4): Initial implementation of acyclic quivers without
                        relations.

A Quiver is a directed graph used for representation theory.  In Sage a Quiver
is different from a directed graph in the following ways:

- The vertices of a DiGraph are arbitrary sage objects, but the vertices of a
  Quiver must be labeled by integers.

- DiGraphs can have cycles (paths that start and end at the same vertex) and
  even loops (edges whose initial and terminal vertices are equal).  In this
  implementation a Quiver must be acyclic (and can not have loops).

- The edges of a DiGraph are labeled with arbitrary sage objects or None if no
  label is specified.  Each edge of a Quiver must be labeled with a nonempty
  string.  The label cannot begin with 'e_' or contain '*' and distinct edges
  must have distinct labels.

- DiGraphs do not have a unique representation in Sage; Quivers do.

Quivers can be described using a dictionary of dictionaries.  The keys of this
dictionary are vertices ``u``.  The value associated to each ``u`` is a dictionary
whose keys are also vertices ``v`` and the value associated to ``v`` is a list of
strings that label the edges from ``u`` to ``v``::

    sage: from sage.modules.quiver_module import *
    sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
    sage: Q.edges()
    [(1, 2, 'a'), (1, 2, 'b'), (1, 3, 'c'), (2, 3, 'd')]

Note that once created we cannot modify the Quiver even though methods from the
DiGraph class that do so appear to be available::

    sage: Q.add_vertex(4)
    Traceback (most recent call last):
    ...
    AttributeError: Quivers cannot be modified.

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

As far as the DiGraph class is concerned a path is a finite list of vertices
[v_1, ..., v_n] such that there exists an edge from v_i to v_{i + 1}.  If there
are multiple edges between the same two vertices this does not contribute to
additional paths as listed by the DiGraph class, for example only two paths are
listed from 1 to 3 in Q::

    sage: Q.all_paths(1, 3)
    [[1, 2, 3], [1, 3]]

When listing paths in a Quiver it is of theoretical importance to distinguish
parallel edges between the same two vertices of a Quiver.  Specifically we say
a path is given by two vertices, start and end, and a finite (possibly empty)
list of edges e_1, e_2, ..., e_n such that the initial vertex of e_1 is start,
the final vertex of e_i is the initial vertex of e_{i + 1}, and the final
vertex of e_n is end.  In the case where no edges are specified we must have
start = end and the path is called the trivial path at the given vertex.  Note
that with this definition there are three paths from 1 to 3 in Q::

    sage: Q.all_quiver_paths(1, 3)
    [b*d, a*d, c]

The all_quiver_paths method returns a list of objects of type QuiverPath; a
type defined by this module.  You can specify a QuiverPath by giving an edge or
a list of edges.  Here an edge is a tuple of the form (i, j, l), where i and j
are vertices and l is the label of an edge from i to j::

    sage: p = QuiverPath([(1, 2, 'a'), (2, 3, 'd')])
    sage: p
    a*d

Trivial paths are indicated by passing the tuple (vertex, vertex)::

    sage: QuiverPath((6, 6))
    e_6

Trivial edges can occur in the input.  They are simply deleted if their vertex
matches the start and end vertex of adjacent edges.  Note that QuiverPaths are
unique::

    sage: q = QuiverPath([(1, 1), (1, 2, 'a'), (2, 2), (2, 3, 'd'), (3, 3)])
    sage: p is q
    True

If the vertex of a trivial path does not match with adjacent edges, or if two
adjacent edges do not match, no error is raised.  Instead the invalid path is
returned.  There is only one invalid path and it can be detected by converting
a QuiverPath to a Boolean.  Valid paths are True, the invalid path is False::

    sage: inv1 = QuiverPath([(1, 2, 'a'), (1, 1)])
    sage: print inv1

    sage: inv2 = QuiverPath([(1, 2, 'a'), (1, 2, 'a')])
    sage: inv1 is inv2
    True
    sage: bool(p)
    True
    sage: bool(inv1)
    False

The * operator is concatenation of paths.  So long as one of the two operands
is a QuiverPath the other can be a QuiverPath, an edge (tuple), or a list of
edges.  If the two paths do not compose the result is the invalid path::

    sage: bool(p*q)
    False
    sage: p*(3, 4, 'c')
    a*d*c
    sage: [(4, 5, 'e'), (5, 1, 'f')]*p
    e*f*a*d

The length of a path is the number of edges in that path.  The invalid path and
trivial paths are therefore length 0::

    sage: len(p)
    2
    sage: triv = QuiverPath((1, 1))
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

The 'in' keyword tests whether a QuiverPath is a path in the given Quiver::

    sage: p in Q
    True
    sage: QuiverPath((1, 4, 'x')) in Q
    False

If 'in' is called and the argument is not a QuiverPath then the call is passed
on to the DiGraph class, where v is in Q if and only if v is a vertex of Q::

    sage: 1 in Q
    True
    sage: (1, 2, 'a') in Q
    False

QuiverPaths form the basis of the quiver algebra of a quiver.  Given a field k
and a Quiver Q the quiver algebra kQ is, as a vector space it has basis the set
of all paths in Q.  Multiplication is defined on this basis and extended
bilinearly.  We multiplication is given as path composition when it makes sense
and is zero otherwise.  Specifically if the terminal vertex of the left path
equals the initial vertex of the right path then their product is the
composition of the two paths, otherwise it is zero.  In sage quiver algebras
are handled by the QuiverAlgebra class::

    sage: A = QuiverAlgebra(GF(7), Q)
    sage: A
    Algebra of Quiver on 3 vertices over Finite Field of size 7

Quivers have a method that creates their algebra over a given field.  Note that
QuiverAlgebras are uniquely defined by their Quiver and field::

    sage: A is Q.algebra(GF(7))
    True
    sage: A is Q.algebra(RR)
    False
    sage: A is Q1.algebra(GF(7))
    False

The QuiverAlgebra can create elements from QuiverPaths or from elements of the
base ring::

    sage: A(5)
    5*e_1 + 5*e_2 + 5*e_3
    sage: r = QuiverPath([(1, 2, 'b'), (2, 3, 'd')])
    sage: e2 = QuiverPath((2, 2))
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

    sage: Z = QuiverRep(GF(5), Q1)
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
Like Quivers, QuiverRep objects are unique and therefore equal if and only if
they are identical::

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
to the constructor along with a path or list of paths that generate the desired
ideal::

    sage: M = QuiverRep(QQ, Q, [[(1, 1)], [(1, 2, 'a')]], option='paths')
    sage: M.dimension_vector()
    (1, 2, 3)

There are also special methods to deal with modules that are given as the
linear dual of a right ideal in the quiver algebra.  To create such a module
pass the keyword option='dual paths' to the constructor along with a path or
list of paths.  The module returned is the dual of the ideal created in the
opposite quiver by the reverse of the given paths::

    sage: D = QuiverRep(QQ, Q, [[(1, 1)], [(1, 2, 'a')]], option='dual paths')
    sage: D.dimension_vector()
    (2, 0, 0)

For modules that are not a standard module or an ideal of the quiver algebra
QuiverRep can take as input two dictionaries.  The first associates to each
vertex a vector space or an integer (the desired dimension of the vector
space), the second associates to each edge a map or a matrix or something from
which sage can construct a map::

    sage: Q2 = Quiver({1:{2:['a', 'b']}})
    sage: M2 = QuiverRep(QQ, Q2, {1: QQ^2, 2: QQ^1}, {(1, 2, 'a'): [1, 0], (1, 2, 'b'): [0, 1]})
    sage: M.get_space(2)
    Vector space of dimension 2 over Rational Field
    sage: M.get_map((1, 2, 'a'))
    Vector space morphism represented by the matrix:
    [1 0]
    Domain: Vector space of dimension 1 over Rational Field
    Codomain: Vector space of dimension 2 over Rational Field

A homomorphism between two quiver representations is given by homomorphisms
between the spaces assigned to the vertices of those representations such that
those homomorphisms commute with the edge maps of the representations.  In Sage
these are created by the QuiverRepHom command.  It takes as input the domain,
codomain, and a dictionary associating maps to vertexes::

    sage: P2 = Q2.P(QQ, 1)
    sage: f = QuiverRepHom(P2, M2, {1:[1, 1], 2:[[1], [1]]})

When the domain is given as a right ideal in the quiver algebra we can also
create a homomorphism by just giving a single element in the codomain.  The map
is then induced by acting on that element::

    sage: x = P2.gens('x')[0]
    sage: x
    x_0
    sage: f == QuiverRepHom(P2, M2, f(x))
    True

As you can see above homomorphisms can be applied to elements.  Just like
elements, addition is defined via the + operator.  On elements scalar
multiplication is defined via the * operator but on homomorphisms * defines
composition, so scalar multiplication is done using a method::

    sage: g = f + f
    sage: g == f.scalar_mult(2)
    True
    sage: g(x) == 2*f(x) # This applies the map, then multiplies by a scalar
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

Both QuiverRep objects and QuiverRepHom objects have ``linear_dual`` and
``algebraic_dual`` methods.  The ``linear_dual`` method applies the functor
Hom_k(-, k) where k is the base ring of the representation and the
``algebraic_dual`` method applies the functor Hom_Q(-, kQ) where kQ is the quiver
algebra.  Both these functors yeild left modules.  A left module is equivalent
to a right module over the opposite algebra and the opposite of a quiver
algebra is the algebra of the opposite quiver, so both these methods yeild
modules and representations of the opposite quiver::

    sage: f.linear_dual()
    Homomorphism of representations of Quiver on 2 vertices
    sage: D = M2.algebraic_dual()
    sage: D.quiver() is Q2.reverse()
    True

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

    sage: w = QuiverRepElement(M2, {1:(1, -1), 2:(3,)})
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

from sage.structure.factory import UniqueFactory
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.free_module import CombinatorialFreeModuleElement
from sage.modules.module import Module
from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.structure.element import ModuleElement
from sage.categories.morphism import CallMorphism
from sage.graphs.digraph import DiGraph

###############################################################################
#
# TODO list for Quiver theory computations:
#
###############################################################################
#
# Phase 1: Finite acyclic quivers without relations
#
#   - The QuiverRep_with_path_basis class and it's dual should have a hidden
#     option where as input you can give  a dictionary instead of a list.  The
#     constructor should assume that the dictionary is valid and not change the
#     order of the entries.  Then the linear dual methods can use this option
#     to make the double dual mesh with unique representation.  This would be
#     better then the funny stuff linear_dual currently does.
#
#   - Maybe there should be default names for the basis elements of a QuiverRep
#     that way the _repr_ method of QuiverRepElement can print something useful
#     instead of essentially just giving the type.  I'm sure the 'proper' thing
#     to do is cache this basis or something like that, but I don't really know
#     and I haven't looked into it yet.
#
#   - Code for extending the base ring
#
#   - Injective envelope code.  QuiverRepHoms should be able to specify a map
#     into one of the dual classes mentioned above by specifying a single
#     element in the domain.  Not exactly sure how that code would work but
#     theoretically everything is dual to the situation with projectives so
#     it shouldn't be hard to figure out when I get around to it.
#
#   - It would be nice if Sage recognized the QuiverHomSpace as a module and
#     we could define addition and scalar multiplication through the coercion
#     system.  But when I tried to make that happen I ran into problems with
#     the * operator doing double duty as scalar multiplication and
#     composition.  I just don't know how actions work so I don't how how to
#     tell the coercion system to do something different based on the type of
#     the operands.
#
###############################################################################
#
# Phase 2: Finite acyclic quivers with relations
#
#   - Figure out the best way to represent relations internally and as
#     input to the constructor.  Much of the code for quivers without
#     relations should work as is for quivers with relations, but the
#     following functions will need to be modified/rewritten/added:
#
#       - Quiver_generic:
#           - Constructor rewritten to take relations as input
#           - _assert_valid_quiver should check that relations are valid
#           - dunder comparisons should require that relations are respected
#           - _repr_ should indicate if relations are present
#           - algebra should give the quotient algebra
#           - all_quiver_paths should only give one representative from each
#             equivalence class of paths
#           - I, P should compute with equivalence classes of paths
#           - should add impose_relation method which returns a quiver equal
#             as a graph to self but with the additional relations
#           - should add drop_relations method which returns a quiver equal
#             as a graph to self but with no relations
#
#       - QuiverAlgebra:
#             Will have to work out how to either quotient the algebra or
#             construct an algebra on equivalence classes of paths.  The first
#             option would be preferable and hopefully there would be a
#             cononical coercion from the algebra sans quotienting.
#
#       - QuiverFactory:
#             Relations should be encoded into the keys in some canonical way
#
#       - QuiverHomSpace:
#           - The constructor should translate the relations of the quiver into
#             additional equations in the system of equations that defines the
#             Hom space.
#
#       - QuiverRep_generic:
#           - _assert_valid_quiverrep should check that the maps respect the
#             relations
#           - Other then ^ I think most of the code here should work as is.  It
#             should definitely be tested though.  There could always be ussues
#             with the dual functors that I'm not anticipating.  But assuming
#             that there are reasonably easy ways to manipulate relations it
#             shouldn't be hard to change whatever needs changing.
#
#       - QuiverRep_with_path_basis:
#           - The constructor will need to associate equivalence classes of
#             paths to the basis instead of just paths.
#           - _left_edge_action will also need to work on equivalence classes
#
#       - QuiverRepHom:
#           - Constructor will need to be able to induce maps from the new
#             QuiverRep_with_path_basis class.  Specifically the option to
#             induce a map by giving a single element of the codomain.
#             Depending on what changes are made to QuiverRep_with_path_basis
#             maybe this code will work as is?
#
#   - The Quiver_generic class should be written first, the others can be
#     modified to return a NotImplementedError when relations are present.
#     Then those functions can be rewritten and tested one by one instead of
#     someone having to do it all at once.
#
###############################################################################
#
# Phase 3: Finite quivers with relations
#
#     If we remove the condition that quivers be acyclic I think a majority of
#     the Quiver and QuiverRep code will work.  The QuiverAlgebra class will
#     need significant changes though.  Do we create an infinite dimensional
#     algebra and then quotient it?  Or do we figure out the equivalence
#     classes ahead of time and use those for the basis?  We would also need
#     to change Quiver_generic.all_quiver_paths; maybe it will give paths
#     without loops and then we add a different method to enumerate loops.  Or
#     maybe it can return some type of collection object with an iterator.
#     That object should maybe be of a type that can be passed directly to
#     CombinatorialFreeModule (if such a thing is possible).
#
#
###############################################################################

def _factory_version():
    """
    Returns the version of sage that UniqueFactory expects.

    TESTS::

        sage: from sage.modules.quiver_module import _factory_version
        sage: _factory_version() #random
        (4, 8)
    """

    from sage.version import version as sage_version
    sage_version = sage_version.split('.')
    for i in range(len(sage_version)):
        try:
            sage_version[i] = int(sage_version[i])
        except ValueError:
            pass
    return tuple(sage_version)

class QuiverFactory(UniqueFactory):
    """
    A Quiver is a directed graph used for representation theory.  In Sage a Quiver
    is different from a directed graph in the following ways:

    - The vertices of a DiGraph are arbitrary sage objects, but the vertices of a
      Quiver must be labeled by integers.

    - DiGraphs can have cycles (paths that start and end at the same vertex) and
      even loops (edges whose initial and terminal vertices are equal).  In this
      implementation a Quiver must be acyclic (and can not have loops).

    - The edges of a DiGraph are labeled with arbitrary sage objects or None if no
      label is specified.  Each edge of a Quiver must be labeled with a nonempty
      string.  The label cannot begin with 'e_' or contain '*' and distinct edges
      must have distinct labels.

    - DiGraphs do not have a unique representation in Sage; Quivers do.

    .. WARNING::

        The Quiver class is derived from the DiGraph class so all the methods of
        the DiGraph class are available to Quivers.  But in order for Quivers to
        have a unique representation in Sage the methods of the DiGraph class that
        alter the structure or labeling of the graph must be disabled.  Thus
        calling the following methods will result in an AttributeError:

          - ``add_cycle``
          - ``add_edge``
          - ``add_edges``
          - ``add_path``
          - ``add_vertex``
          - ``add_vertices``
          - ``allow_loops``
          - ``allow_multiple_edges``
          - ``clear``
          - ``set_edge_label``
          - ``subdivide_edge``
          - ``subdivide_edges``

    INPUT:

    See the documentation of the DiGraph class for the full specification of
    possible inputs.  The examples in the documentation of functions relating to
    quivers will use the following specification:

    - ``data`` - dictionary of dictionaries (default: empty), the keys of this
      dictionary are vertices ``u``.  The value associated to each ``u`` is a dictionary
      whose keys are also vertices ``v`` and the value associated to ``v`` is a list of
      strings that label the edges from ``u`` to ``v``.

    OUTPUT:

    - Quiver

    EXAMPLES::

        sage: from sage.modules.quiver_module import Quiver
        sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
        sage: Q.edges()
        [(1, 2, 'a'), (1, 2, 'b'), (1, 3, 'c'), (2, 3, 'd')]

    Note that once created we cannot modify the Quiver::

        sage: Q.add_vertex(4)
        Traceback (most recent call last):
        ...
        AttributeError: Quivers cannot be modified.

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

    Unique representation means that if we create a Quiver that is equal to the
    first then it is in fact identical::

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

    When listing paths in a Quiver it is of theoretical importance to distinguish
    parallel edges between the same two vertices of a Quiver::

        sage: Q.all_quiver_paths(1, 3)
        [b*d, a*d, c]

    Compare this with the output of the DiGraph method, which lists paths by the
    vertices traversed, and hence lists only two paths, instead of three::

        sage: Q.all_paths(1, 3)
        [[1, 2, 3], [1, 3]]

    The all_quiver_paths method returns a list of objects of type QuiverPath.  See
    the documentation of the QuiverPath class for more information.
    """

    def create_key(self, *args, **kwds):
        """
        Returns a key for the specified Quiver.

        The key is a tuple.  The first entry of the tuple is a tuple containing the
        vertices of the quiver, in the order returned by calling ``DiGraph.vertices()``.
        The following entries in the tuple are the edges of the quiver, sorted as in
        ``DiGraph.edges(sort=True)``.

        INPUT:

        - See the documentation for DiGraph

        OUTPUT:

        - tuple of tuples

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Quiver.create_key({1:{2:['a', 'b'], 4:['hi']}, 2:{3:['x', 'a']}})
            ((1, 2, 3, 4), (1, 2, 'a'), (1, 2, 'b'), (1, 4, 'hi'), (2, 3, 'a'), (2, 3, 'x'))
        """

        # Create a DiGraph and construct a tuple from its edge list
        G = DiGraph(*args, **kwds)
        key = [tuple(G.vertices())]
        key.extend(G.edges(labels=True, sort=True))
        return tuple(key)

    def create_object(self, version, key, **extra_args):
        """
        Creates a Quiver_generic object from the key.

        The key is a tuple.  The first entry of the tuple is a tuple containing the
        vertices of the quiver, in the order returned by calling ``DiGraph.vertices()``.
        The following entries in the tuple are the edges of the quiver, sorted as in
        ``DiGraph.edges(sort=True)``.

        INPUT:

        - ``version`` - The version of sage, this is currently ignored.

        - ``key`` - tuple of tuples

        OUTPUT:

        - Quiver_generic object

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, _factory_version
            sage: Q = Quiver.create_object(_factory_version(), ((1, 2, 3), (1, 2, 'a'), (1, 3, 'b')))
            sage: Q.vertices()
            [1, 2, 3]
            sage: Q.edges()
            [(1, 2, 'a'), (1, 3, 'b')]
        """

        # Create the edge dictionary
        data = dict((v, {}) for v in key[0])
        for e in key[1:]:
            if e[1] not in data[e[0]]:
                data[e[0]][e[1]] = []
            data[e[0]][e[1]].append(e[2])

        # Create and return the quiver
        return Quiver_generic(data)

Quiver = QuiverFactory("Quiver")

class QuiverRepFactory(UniqueFactory):
    """
    A quiver representation is a diagram in the category of vector spaces whose
    underlying graph is the quiver.  Giving a finite dimensional representation
    is equivalent to giving a finite dimensional right module for the path
    algebra of the quiver.

    INPUT:

    The first two arguments specify the base ring and the Quiver, and they are
    always required:

    - ``k`` - ring, the base ring of the representation

    - ``Q`` - Quiver, the quiver of the representation

    Then to specify the spaces and maps associated to the Quiver there are three
    possible options.  The first is the ``values`` option, where the next two
    arguments give the data to be assigned.  The following can either be the next
    two entries in the argument list or they can be passed by keyword.  If the
    argument list is long enough the keywords are ignored; the keywords are only
    checked in the event that the argument list does not have enough entries after
    ``Q``.

    - ``spaces`` - dict (default: empty), a dictionary associating to each vertex a
      free module over the base ring k.  Not all vertices must be specified,
      unspecified vertices are automatically set to k^0.  Keys of the dictionary
      that don't correspond to vertices are ignored.

    - ``maps`` - dict (default: empty), a dictionary associating to each edge a map
      whose domain and codomain are the spaces associated to the initial and
      terminal vertex of the edge respectively.  Not all edges must be specified,
      unspecified edges are automatically set to the zero map.  Keys of the
      dictionary that don't correspond to edges are ignored.

    The second option is the ``paths`` option which creates a module by generating a
    right ideal from a list of paths.  Thus the basis elements of this module
    correspond to paths of the quiver and the maps are given by right
    multiplication by the corresponding edge.  As above this can be passed either
    as the next entry in the argument list or as a keyword.  The keyword is only
    checked if there is no entry in the argument list after ``Q``.

    - ``basis`` - list, a nonempty list of paths in the quiver Q.  Entries that do not
      represent valid paths are ignored and duplicate paths are deleted.  There
      must be at least one valid path in the list or a ValueError is raised.  The
      closure of this list under right multiplication forms the basis of the
      resulting representation.

    The third option is the ``dual paths`` option which creates the dual of a left
    ideal in the quiver algebra.  Thus the basis elements of this module correspond
    to paths of the quiver and the maps are given by deleting the corresponding
    edge from the start of the path (the edge map is zero on a path if that edge is
    not the initial edge of the path).  As above this can be passed either as the
    next entry in the argument list or as a keyword.

    - ``basis`` - list, a nonempty list of paths in the quiver Q.  Entries that do not
      represent valid paths are ignored and duplicate paths are deleted.  There
      must be at least one valid path in the list or a ValueError is raised.  The
      closure of this list under left multiplication of edges forms the basis of
      the resulting representation.

    Using the second and third options requires that the following keyword be
    passed to the constructor.  This must be passed as a keyword.

    - ``option`` - string (default: None), either 'values' or 'paths'.  None is
      equivalent to 'values'.

    OUTPUT:

    - QuiverRep

    EXAMPLES::

        sage: from sage.modules.quiver_module import Quiver, QuiverRep
        sage: Q1 = Quiver({1:{2:['a']}})

    When the ``option`` keyword is not supplied the constructor uses the 'values'
    option and expects the spaces and maps to be specified.  If no maps or spaces
    are given the zero module is created::

        sage: M = QuiverRep(GF(5), Q1)
        sage: M.is_zero()
        True

    The simple modules, indecomposable projectives, and indecomposable injectives are
    examples of quiver representations::

        sage: S = Q1.S(GF(3), 1)
        sage: I = Q1.I(QQ, 2)
        sage: P = Q1.P(GF(3), 1)

    Various standard submodules can be computed, such as radicals and socles.  We can
    also form quotients and test for certain attributes such as semisimplicity::

        sage: R = P.radical()
        sage: R.is_zero()
        False
        sage: (P/R).is_simple()
        True
        sage: P == R
        False

    With the option 'paths' the input data should be a list of QuiverPaths or
    things that QuiverPaths can be constructed from.  The resulting module is the
    submodule generated by these paths in the quiver algebra, when considered as a
    right module over itself::

        sage: P1 = QuiverRep(QQ, Q1, [[(1, 1)]], option='paths')
        sage: P1.dimension()
        2

    Notice that the 3rd and 4th paths are actually the same, the duplicate is
    removed::

        sage: N = QuiverRep(QQ, Q1, [[(1, 1)], [(2, 2)], [(1, 1), (1, 2, 'a'), (2, 2)], [(1, 2, 'a')]], option='paths')
        sage: N.dimension()
        3

    The dimension at each vertex equals the number of paths in the closed basis whose
    terminal point is that vertex::

        sage: Q2 = Quiver({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}})
        sage: M = QuiverRep(QQ, Q2, [[(2, 2)], [(1, 2, 'a')]], option='paths')
        sage: M.dimension_vector()
        (0, 2, 2)
        sage: N = QuiverRep(QQ, Q2, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
        sage: N.dimension_vector()
        (0, 1, 2)
    """

    def create_key(self, k, Q, *args, **kwds):
        """
        Returns a key for the specified module.

        The key is a tuple.  The first and second entries are the base ring ``k`` and the
        key of the Quiver ``Q``.  The third entry is the 'option' and the remaining entries
        depend on that option.  If the option is 'values' and the Quiver has n vertices
        then the next n entries are the vector spaces to be assigned to those vertices.
        After that are the matrices of the maps assigned to edges, listed in the same
        order that ``Q.edges()`` uses.  If the option is 'paths' or 'dual paths' then the
        next entry is a tuple containing a sorted list of the paths that form a basis
        of the Quiver.

        INPUT:

        - See the class documentation

        OUTPUT:

        - tuple

        EXAMPLES::

            sage: from sage.modules.quiver_module import *
            sage: Q = Quiver({1:{2:['a']}})
            sage: QuiverRep.create_key(GF(5), Q)
            (Finite Field of size 5, ((1, 2), (1, 2, 'a')), 'values', Vector space of dimension 0 over Finite Field of size 5, Vector space of dimension 0 over Finite Field of size 5, [])
        """

        key = [k, Quiver.create_key(Q)]
        if 'option' in kwds and (kwds['option'] == 'paths' or kwds['option'] == 'dual paths'):
            # Follow the 'paths' specification for the input
            key.append(kwds['option'])
            if args:
                basis = args[0]
            else:
                basis = kwds['basis']

            # Add as QuiverPaths to a set
            paths = set()
            for p in basis:
                paths.add(QuiverPath(p))

            if kwds['option'] == 'paths':
                # Close the set under right mult by edges
                edges = [QuiverPath(e) for e in Q.edges()]
                just_added = paths
                while just_added:
                    to_be_added = []
                    for e in edges:
                        for p in just_added:
                            if p*e not in paths:
                                to_be_added.append(p*e)

                    paths.update(to_be_added)
                    just_added = to_be_added

            if kwds['option'] == 'dual paths':
                # Close the set under left mult by edges
                edges = [QuiverPath(e) for e in Q.edges()]
                just_added = paths
                while just_added:
                    to_be_added = []
                    for e in edges:
                        for p in just_added:
                            if e*p not in paths:
                                to_be_added.append(e*p)

                    paths.update(to_be_added)
                    just_added = to_be_added

            # Remove the invalid path
            paths.difference_update([QuiverPath()])

            # Add to the key
            key.append(tuple(sorted(paths)))

        else:
            # Assume the input type the 'values' option
            key.append('values')
            if args:
                spaces = args[0]
            elif 'spaces' in kwds:
                spaces = kwds['spaces']
            else:
                spaces = {}
            if len(args) > 1:
                maps = args[1]
            elif 'maps' in kwds:
                maps = kwds['maps']
            else:
                maps = {}

            # If the vertex is not specified set it as a free module of rank 0, if
            # an integer is given set it as a free module of that rank, otherwise
            # assume the object is a module and assign it to the vertex.
            from sage.rings.finite_rings.integer_mod_ring import Integers
            verts = Q.vertices()
            for x in verts:
                if x not in spaces:
                    key.append(k**0)
                elif spaces[x] in Integers():
                    key.append(k**spaces[x])
                else:
                    key.append(spaces[x])

            # The preferred method of specifying an edge is as a tuple
            # (i, t, l) where i is the initial vertex, t is the terminal
            # vertex, and l is the label.  This is the form in which
            # quiver.edges() and other such functions give the edge.  But here
            # edges can be specified by giving only the two vertices or giving
            # only the edge label.
            #
            # Note that the first space is assigned to key[3] and the first
            # vertex is 1 so the space assigned to vertex v is key[2 + v]
            from sage.matrix.constructor import Matrix
            from sage.categories.morphism import is_Morphism
            for x in Q.edges():
                if x in maps:
                    e = maps[x]
                elif (x[0], x[1]) in maps:
                    e = maps[(x[0], x[1])]
                elif x[2] in maps:
                    e = maps[x[2]]
                else:
                    e = Matrix(k, key[3 + verts.index(x[0])].dimension(), key[3 + verts.index(x[1])].dimension())

                # If a morphism is specified take it's matrix.  Create one if
                # needed.  Otherwise assume the Matrix function can convert the
                # object to a Matrix.
                if is_Morphism(e):
                    if hasattr(e, 'matrix'):
                        key.append(e.matrix())
                    else:
                        gens_images = [key[3 + verts.index(x[1])].coordinate_vector(e(x))
                                        for x in key[3 + verts.index(x[0])].gens()]
                        key.append(Matrix(k, key[3 + verts.index(x[0])].dimension(),
                                          key[3 + verts.index(x[1])].dimension(), gens_images))
                else:
                    key.append(Matrix(k, key[3 + verts.index(x[0])].dimension(),
                                      key[3 + verts.index(x[1])].dimension(), e))

                # Make sure the matrix is immutable so it hashes
                key[-1].set_immutable()

        # Wrap as a tuple and return
        return tuple(key)

    def create_object(self, version, key, **extra_args):
        """
        Creates a QuiverRep_generic or QuiverRep_with_path_basis object from the key.

        The key is a tuple.  The first and second entries are the base ring ``k`` and the
        key of the Quiver ``Q``.  The third entry is the 'option' and the remaining entries
        depend on that option.  If the option is 'values' and the Quiver has n vertices
        then the next n entries are the vector spaces to be assigned to those vertices.
        After that are the matrices of the maps assigned to edges, listed in the same
        order that ``Q.edges()`` uses.  If the option is 'paths' or 'dual paths' then the
        next entry is a tuple containing a sorted list of the paths that form a basis
        of the Quiver.

        INPUT:

        - ``version`` - The version of sage, this is currently ignored.

        - ``key`` - tuple

        OUTPUT:

        - QuiverRep_generic or QuiverRep_with_path_basis object

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, _factory_version
            sage: Q = Quiver({1:{2:['a']}})
            sage: key = QuiverRep.create_key(GF(5), Q)
            sage: QuiverRep.create_object(_factory_version(), key)
            Representation with dimension vector (0, 0)
        """

        if len(key) < 4:
            raise ValueError("Invalid key used in QuiverRepFactory!")

        # Get the quiver
        Q = Quiver.get_object(version, key[1], {})

        if key[2] == 'values':
            # Get the spaces
            spaces = {}
            i = 3
            for v in Q:
                spaces[v] = key[i]
                i += 1

            # Get the maps
            maps = {}
            for e in Q.edges():
                maps[e] = key[i]
                i += 1

            # Create and return the module
            return QuiverRep_generic(key[0], Q, spaces, maps)

        elif key[2] == 'paths':
            # Create and return the module
            return QuiverRep_with_path_basis(key[0], Q, key[3])

        elif key[2] == 'dual paths':
            # Create and return the module
            return QuiverRep_with_dual_path_basis(key[0], Q, key[3])

        else:
            raise ValueError("Invalid key used in QuiverRepFactory!")

QuiverRep = QuiverRepFactory("QuiverRep")

class QuiverHomSpaceFactory(UniqueFactory):
    """
    A homomorphism of quiver representations is for each vertex of the quiver a
    homomorphism of the spaces assigned to those vertices such that these
    homomorphisms commute with the edge maps.  This class handles the set of all
    such maps, Hom_Q(M, N).

    INPUT:

    - ``domain`` - QuiverRep, the domain of the homomorphism space

    - ``codomain`` - QuiverRep, the codomain of the homomorphism space

    OUTPUT:

    - QuiverHomSpace, the homomorphism space Hom_Q(domain, codomain)

    .. NOTES::

        The quivers of the domain and codomain must be equal or a ValueError is
        raised.

    EXAMPLES::

        sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
        sage: Q = Quiver({1:{2:['a', 'b']}})
        sage: H = QuiverHomSpace(Q.S(QQ, 2), Q.P(QQ, 1))
        sage: H.dimension()
        2
        sage: H.gens()
        [Homomorphism of representations of Quiver on 2 vertices, Homomorphism of representations of Quiver on 2 vertices]
    """

    def create_key(self, domain, codomain):
        """
        Returns a key for the specified Hom space.

        The key is a tuple of length two, the two entries being the keys for the domain
        and codomain.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: QuiverHomSpace.create_key(Q.P(GF(3), 1), Q.S(GF(3), 1))
            ((Finite Field of size 3, ((1, 2), (1, 2, 'a'), (1, 2, 'b')), 'paths', (a, b, e_1)), (Finite Field of size 3, ((1, 2), (1, 2, 'a'), (1, 2, 'b')), 'values', Vector space of dimension 1 over Finite Field of size 3, Vector space of dimension 0 over Finite Field of size 3, [], []))
        """

        # Check that the quivers and base rings are the same
        if domain._quiver is not codomain._quiver:
            raise ValueError("The quivers of the domain and codomain must be identical.")
        if domain._base_ring is not codomain._base_ring:
            raise ValueError("The base rings of the domain and codomain must be identical.")

        # Create and return the key
        return tuple([domain._factory_data[2], codomain._factory_data[2]])

    def create_object(self, version, key, **extra_args):
        """
        Creates a QuiverHomSpace from the key.

        The key is a tuple of length two, the two entries being the keys for the domain
        and codomain.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace, _factory_version
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: key = QuiverHomSpace.create_key(Q.P(GF(3), 1), Q.S(GF(3), 1))
            sage: QuiverHomSpace.create_object(_factory_version(), key)
            Dimension 1 QuiverHomSpace
        """

        return QuiverHomSpace_generic(QuiverRep.get_object(version, key[0], []),
                                      QuiverRep.get_object(version, key[1], []))

QuiverHomSpace = QuiverHomSpaceFactory("QuiverHomSpace")

class Quiver_generic(DiGraph):
    """
    Generic class for a quiver.  Type Quiver? for more information.

    TESTS::

        sage: from sage.modules.quiver_module import Quiver
        sage: Q1 = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, data=None, pos=None, loops=False, format=None,
             boundary=[], weighted=True, implementation='c_graph',
             sparse=True, vertex_labels=True, multiedges=True, **kwds):
        """
        Creates a Quiver object.  Type Quiver? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
        """
        #print "boundary=",boundary
        # Use the Digraph construction
        super(Quiver_generic, self).__init__(data, pos, False, format, boundary,
             True, implementation, sparse, vertex_labels, **kwds)

        # Ensure the input is valid
        self._assert_valid_quiver()

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

    def _assert_valid_quiver(self):
        """
        Raises an exception if self does not satisfy the following:

        - Each edge must have a unique string label.

        - Labels of edges can't begin with 'e_' or contain '*'.

        - Vertices must be labeled by the integers.

        - Q must be acyclic (no loops).

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
            sage: goodQ = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}}) # indirect doctest
            sage: badQ1 = Quiver({1:{2:['a',2], 3:['c']}, 2:{3:['d']}}) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Edge labels of Representation Quiver must be nonempty strings.
            sage: badQ2 = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['a']}}) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Edge labels of Representation Quiver must be unique.
        """

        # Check that it's directed, acyclic, and has no loops
        if not self.is_directed_acyclic():
            raise ValueError("Representation Quiver must be directed acyclic.")
        if self.has_loops():
            raise NotImplementedError("Representation Quiver may not have loops.")

        # Check that vertices are labeled 1,2,3,... and that edge labels are
        # unique
        from sage.rings.finite_rings.integer_mod_ring import Integers
        for v in self:
            if v not in Integers():
                raise ValueError("Vertices of Representation Quiver must be labeled by integers.")
        if len(set(self.edge_labels())) != self.num_edges():
            raise ValueError("Edge labels of Representation Quiver must be unique.")

        # Check that edges are labeled with nonempty strings and don't begin
        # with 'e_' or contain '*'
        for x in self.edge_labels():
            if not isinstance(x, str) or x == '':
                raise ValueError("Edge labels of Representation Quiver must be nonempty strings.")
            if x[0:2] == 'e_':
                raise ValueError("Edge labels of Representation Quiver cannot begin with 'e_'.")
            if x.find('*') != -1:
                raise ValueError("Edge labels of Representation Quiver cannot contain '*'.")

    def _repr_(self):
        """
        Default string representation of a quiver.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
            sage: Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}}) # indirect doctest
            Quiver on 3 vertices
        """

        return "Quiver on " + str(len(self)) + " vertices"

    def _forbidden_method(self, *args):
        """
        Functions that modify the quiver may not be called; quivers are unique.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a']}})
            sage: Q.add_vertex(3) # indirect doctest
            Traceback (most recent call last):
            ...
            AttributeError: Quivers cannot be modified.

        To modify a quiver you must use the to_directed method to create a DiGraph,
        modify the DiGraph, and then create a Quiver from the DiGraph::

            sage: G = Q.to_directed()
            sage: G.add_vertex(3)
            sage: Q = Quiver(G)
            sage: Q.vertices()
            [1, 2, 3]
        """

        raise AttributeError("Quivers cannot be modified.")

    def __contains__(self, other):
        """
        Implements the ``in`` keyword.

        If other is a QuiverPath then ``other in self`` returns True if and only if
        each edge of other is an edge of self.  In the case of a trivial path at a
        vertex the vertex must be a vertex of the graph.  The invalid path belongs to
        every quiver.  If other is not a QuiverPath then the call is passed up to
        ``DiGraph.__contains__()`` which returns True if and only if other is a vertex.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverPath
            sage: Q = Quiver({1: {2:['a']}, 2:{3:['b']}})
            sage: x = QuiverPath([(1, 2, 'a'), (2, 3, 'b')])
            sage: y = x*(3, 4, 'c')
            sage: x in Q
            True
            sage: y in Q
            False
            sage: 1 in Q
            True
            sage: 4 in Q
            False
            sage: QuiverPath((1, 1)) in Q
            True
            sage: QuiverPath((4, 4)) in Q
            False
            sage: QuiverPath([(1, 1), (2, 2)]) in Q
            True
        """

        # Pass up if not a QuiverPath
        if not isinstance(other, QuiverPath):
            return super(Quiver_generic, self).__contains__(other)

        # Empty paths are in every quiver
        if not other._path:
            return True

        # If its a trivial path just check the vertex
        if other._path[0][0] == other._path[0][1]:
            return super(Quiver_generic, self).__contains__(other._path[0][0])

        # Otherwise check each edge
        for e in other:
            if (not super(Quiver_generic, self).__contains__(e[0]) or
                    not super(Quiver_generic, self).__contains__(e[1]) or
                    e[2] not in self.edge_label(e[0], e[1])):
                return False
        return True

    def __lt__(self, other):
        """
        Returns true if other is a quiver and the vertices and edges of self are
        vertices and edges of other, but self is not other.  Raises a TypeError if
        other is not a Quiver.

        OUTPUT:

        - bool

        TESTS::

            sage: from sage.modules.quiver_module import *
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
        if not isinstance(other, Quiver_generic):
            raise TypeError("Can only compare Quivers with other Quivers")

        # Check for equality
        if self is other:
            return False

        # Check vertices and edges for containment
        for v in self:
            if v not in other:
                return False

        for e in self.edges():
            if QuiverPath(e) not in other:
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

            sage: from sage.modules.quiver_module import *
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
        if not isinstance(other, Quiver_generic):
            raise TypeError("Can only compare Quivers with other Quivers")

        # Check for equality
        if self is other:
            return False

        # Check vertices and edges for containment
        for v in other:
            if v not in self:
                return False

        for e in other.edges():
            if QuiverPath(e) not in self:
                return False

        return True

    def __le__(self, other):
        """
        Returns true if other is a quiver and the vertices and edges of self are
        vertices and edges of other.  Raises a TypeError if other is not a Quiver.

        OUTPUT:

        - bool

        TESTS::

            sage: from sage.modules.quiver_module import *
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
        if not isinstance(other, Quiver_generic):
            raise TypeError("Can only compare Quivers with other Quivers")

        # Check for equality
        if self is other:
            return True

        # Check vertices and edges for containment
        for v in self:
            if v not in other:
                return False

        for e in self.edges():
            if QuiverPath(e) not in other:
                return False

        return True

    def __ge__(self, other):
        """
        Returns true if other is a quiver and the vertices and edges of other are
        vertices and edges of self.  Raises a TypeError if other is not a Quiver.

        OUTPUT:

        - bool

        TESTS::

            sage: from sage.modules.quiver_module import *
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
        if not isinstance(other, Quiver_generic):
            raise TypeError("Can only compare Quivers with other Quivers")

        # Check for equality
        if self is other:
            return True

        # Check vertices and edges for containment
        for v in other:
            if v not in self:
                return False

        for e in other.edges():
            if QuiverPath(e) not in self:
                return False

        return True

    ###########################################################################
    #                                                                         #
    # GRAPH THEORETIC FUNCTIONS                                               #
    #    These functions involve the graph theoretic structure of the quiver. #
    #                                                                         #
    ###########################################################################

    def all_quiver_paths(self, start=None, end=None):
        """
        Returns a list of all paths between a pair of vertices (start, end).

        INPUT:

        - ``start`` - integer or None (default: None), the initial vertex of the paths in
          the output.  If None is given then the initial vertex is arbitrary.

        - ``end`` - integer or None (default: None), the terminal vertex of the paths in
          the output.  If None is given then the terminal vertex is arbitrary.

        OUTPUT:

        - list of QuiverPaths

        .. NOTE::

            If there are multiple edges between two vertices all_paths will not
            differentiate between them but all_quiver_paths will.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: Q.all_quiver_paths(1, 3)
            [b*d, a*d, c]

        If start=end then we expect only the trivial path at that vertex::

            sage: Q.all_quiver_paths(1, 1)
            [e_1]

        The empty list is returned if there are no paths between the given vertices::

            sage: Q.all_quiver_paths(3, 1)
            []

        If end=None then all edge paths beginning at start are returned, including the
        trivial path::

            sage: Q.all_quiver_paths(2)
            [e_2, d]

        If start=None then all edge paths ending at end are returned, including the
        trivial path.  Note that the two edges from vertex 1 to vertex 2 count as two
        different edge paths::

            sage: Q.all_quiver_paths(None, 2)
            [b, a, e_2]
            sage: Q.all_quiver_paths(end=2)
            [b, a, e_2]

        If start=end=None then all edge paths are returned, including trivial paths::

            sage: Q.all_quiver_paths()
            [e_1, b, a, b*d, a*d, c, e_2, d, e_3]

        The vertex given must be a vertex of the quiver::

            sage: Q.all_quiver_paths(1, 4)
            Traceback (most recent call last):
            ...
            ValueError: The end vertex 4 is not a vertex of the quiver.
        """

        # Check that given arguments are vertices
        if start is not None and start not in self:
            raise ValueError("The start vertex " + str(start) + " is not a vertex of the quiver.")
        if end is not None and end not in self:
            raise ValueError("The end vertex " + str(end) + " is not a vertex of the quiver.")

        # Handle start=None
        if start == None:
            results = []
            for v in self:
                results += self.all_quiver_paths(v, end)
            return results

        # Handle end=None
        if end == None:
            results = []
            for v in self:
                results += self.all_quiver_paths(start, v)
            return results

        # Handle the trivial case
        if start == end:
            return [QuiverPath((start, end))]

        # This function will recursively convert a path given in terms of
        # vertices to a list of QuiverPaths.
        def _v_to_e(path):
            if len(path) == 1:
                return [QuiverPath((path[0], path[0]))]
            paths = []
            for a in self.edge_label(path[0], path[1]):
                for b in _v_to_e(path[1:]):
                    paths.append(QuiverPath((path[0], path[1], a))*b)
            return paths

        # For each vertex path we append the resulting edge paths
        result = []
        for path in self.all_paths(start, end):
            result += _v_to_e(path)

        # The result is all paths from start to end
        return result

    def to_directed(self):
        """
        Returns the underlying Digraph.

        OUTPUT:

        - a DiGraph with the same edges and vertices as the quiver

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
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

            sage: from sage.modules.quiver_module import Quiver
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

            sage: from sage.modules.quiver_module import Quiver
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

            sage: from sage.modules.quiver_module import Quiver
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

            sage: from sage.modules.quiver_module import Quiver
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

            sage: from sage.modules.quiver_module import Quiver
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
        return super(Quiver_generic, self).all_paths(start, end)

    def reverse(self):
        """
        Returns a copy of the quiver with edges reversed in direction.

        OUTPUT:

        - Quiver

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
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

    def representation(self, k, *args, **kwds):
        """
        Returns a representation of the quiver.

        For more information see the QuiverRep documentation.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep_generic
            sage: Q = Quiver({1:{3:['a']}, 2:{3:['b']}})
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = Q.representation(QQ, spaces, maps)
        """

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

            sage: from sage.modules.quiver_module import Quiver
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

            sage: from sage.modules.quiver_module import Quiver
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

            sage: from sage.modules.quiver_module import Quiver
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

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b'], 3: ['c', 'd']}, 2:{3:['e']}})
            sage: Q.free_module(GF(3)).dimension_vector()
            (1, 3, 6)
        """

        return QuiverRep(k, self, [[(v, v)] for v in self], option='paths')

    def algebra(self, k):
        """
        Return the algebra of the quiver.

        INPUT:

        - ``k`` - ring, the base ring of the quiver algebra

        OUTPUT:

        - QuiverAlgebra

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{4:['c']}, 3:{4:['d']}})
            sage: Q.algebra(GF(7))
            Algebra of Quiver on 4 vertices over Finite Field of size 7
        """

        return QuiverAlgebra(k, self)

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

        sage: from sage.modules.quiver_module import Quiver, QuiverAlgebra, QuiverPath
        sage: Q = Quiver({1:{2:['a']}, 2:{3:['b']}})
        sage: A = QuiverAlgebra(GF(7), Q)
        sage: A
        Algebra of Quiver on 3 vertices over Finite Field of size 7

    Quivers have a method that creates their algebra over a given field.  Note that
    QuiverAlgebras are uniquely defined by their Quiver and field::

        sage: A is Q.algebra(GF(7))
        True
        sage: A is Q.algebra(RR)
        False
        sage: A is Quiver({1:{2:['a']}}).algebra(GF(7))
        False

    The QuiverAlgebra can create elements from QuiverPaths or from elements of the
    base ring::

        sage: A(5)
        5*e_1 + 5*e_2 + 5*e_3
        sage: p = QuiverPath((1, 2, 'a'))
        sage: r = QuiverPath((2, 3, 'b'))
        sage: e2 = QuiverPath((2, 2))
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
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    @staticmethod
    def __classcall__(self, k, Q):
        """
        Convert the Quiver Q into its key so that unique representation works.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q1 = Quiver({4:{2:['a']}})
            sage: Q2 = Quiver(DiGraph({4:{2:['a']}}))
            sage: Q1 is Q2
            True
        """

        Qkey = Quiver.create_key(Q)
        return super(QuiverAlgebra, self).__classcall__(self, k, Qkey)

    def __init__(self, k, Qkey):
        """
        Creates a QuiverAlgebra object.  Type QuiverAlgebra? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverAlgebra
            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b', 'c']}, 4:{}})
            sage: QuiverAlgebra(GF(5), Q)
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

        from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
        self._quiver = Quiver.get_object(_factory_version(), Qkey, [])
        super(QuiverAlgebra, self).__init__(k, self._quiver.all_quiver_paths(),
                                            prefix='',
                                            element_class=self.Element,
                                            category=GradedAlgebrasWithBasis(k),
                                            bracket=False)

    def _element_constructor_(self, x):
        """
        Attempts to construct an element of self from x.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
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
            return self.monomial(x)

        # If it's a tuple or a list try and create a QuiverPath from it and
        # then return the associated basis element
        if isinstance(x, tuple) or isinstance(x, list):
            return self.monomial(QuiverPath(x))

        # Otherwise let CombinatorialFreeModule try
        return super(QuiverAlgebra, self)._element_constructor_(x)

    def _coerce_map_from_(self, other):
        """
        True if there is a coercion from other to self.

        The algebras that coerce into a quiver algebra are rings k or quiver algebras
        kQ such that k has a coercion into the base ring of self and Q is a subquiver
        of the quiver of self.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
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
        """

        if (isinstance(other, QuiverAlgebra) and
                self._base.has_coerce_map_from(other._base) and
                other._quiver <= self._quiver):
            return True
        else:
            return self._base.has_coerce_map_from(other)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverAlgebra
            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b']}})
            sage: QuiverAlgebra(RR, Q) # indirect doctest
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

    def product_on_basis(self, p1, p2):
        """
        Returns the element corresponding to p1*p2 in the quiver algebra.

        INPUT:

        - ``p1``, ``p2`` - QuiverPaths

        OUTPUT:

        - CombinatorialFreeModuleElement

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverPath, QuiverAlgebra
            sage: Q = Quiver({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}})
            sage: p1 = QuiverPath((1, 2, 'a'))
            sage: p2 = QuiverPath([(2, 3, 'b'), (3, 4, 'c')])
            sage: A = QuiverAlgebra(QQ, Q)
            sage: A.product_on_basis(p1, p2)
            a*b*c
        """

        p = QuiverPath(p1)*QuiverPath(p2)
        if p:
            return self.basis()[p]
        else:
            return self.zero()

    def degree_on_basis(self, p):
        """
        Return the degree of the monomial specified by p.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: A = Quiver({1:{2:['a']}, 2:{3:['b']}}).algebra(QQ)
            sage: A.degree_on_basis((1, 1))
            0
            sage: A.degree_on_basis((1, 2, 'a'))
            1
            sage: A.degree_on_basis([(1, 2, 'a'), (2, 3, 'b')])
            2
        """

        return len(QuiverPath(p))

    def one(self):
        """
        Return the multiplicative identity element.

        The multiplicative identity of a quiver algebra is the sum of the basis
        elements corresponding to the trivial paths at each vertex.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: A = Quiver({1:{2:['a']}, 2:{3:['b']}}).algebra(QQ)
            sage: A.one()
            e_1 + e_2 + e_3
        """

        x = self.zero()
        B = self.basis()
        for v in self._quiver:
            x += B[QuiverPath((v, v))]
        return x

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data and subspaces of the quiver algebra.     #
    #                                                                         #
    ###########################################################################

    def quiver(self):
        """
        Return the quiver of the representation.

        OUTPUT:

        - Quiver

        EXAMPLES:

            sage: from sage.modules.quiver_module import Quiver, QuiverAlgebra
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: A = QuiverAlgebra(GF(3), Q)
            sage: A.quiver() is Q
            True
        """

        return self._quiver

    def homogeneous_component(self, n):
        """
        Return the nth homogeneous piece of the quiver algebra.

        INPUT:

        - ``n`` - integer

        OUTPUT:

        - CombinatorialFreeModule, module spanned by the paths of lenth n in the
          quiver.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
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

                sage: from sage.modules.quiver_module import Quiver
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

                sage: from sage.modules.quiver_module import Quiver
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

class QuiverPath(UniqueRepresentation, SageObject):
    """
    Class for paths in a quiver.

    A path is given by two vertices, start and end, and a finite (possibly empty)
    list of edges e_1, e_2, ..., e_n such that the initial vertex of e_1 is start,
    the final vertex of e_i is the initial vertex of e_{i + 1}, and the final
    vertex of e_n is end.  In the case where no edges are specified we must have
    start = end and the path is called the trivial path at the given vertex.  This
    class can also represent an invalid path.

    INPUT:

    - ``path`` - tuple or iterable (default: empty list), if ``path`` is a tuple then it is
      assumed to be of the form (vertex, vertex) or (start, end, label).  In the
      first case the trivial path at the given vertex is created.  In the second
      case a path consisting of just the given edge is created.  If ``path`` is not a
      tuple then it is assumed to be an iterable variable giving the edges of a
      path, where each edge is in one of the two forms above.  If these edges do
      not compose to form a valid path then an invalid path is returned.

    OUTPUT:

    - QuiverPath

    EXAMPLES::

    Specify a path by giving a list of edges::

        sage: from sage.modules.quiver_module import QuiverPath
        sage: p = QuiverPath([(1, 2, 'a'), (2, 2), (2, 3, 'b')])
        sage: p
        a*b

    Paths are unique::

        sage: q = QuiverPath([(1, 1), (1, 2, 'a'), (2, 3, 'b'), (3, 3)])
        sage: p is q
        True

    In particular there is only one invalid path.  It can be detected by converting
    a path to a Boolean.  Valid paths are True, the invalid path is False::

        sage: inv1 = QuiverPath([(1, 1), (2, 2)])
        sage: inv2 = QuiverPath([(1, 2, 'a'), (1, 2, 'a')])
        sage: inv1 is inv2
        True
        sage: bool(p)
        True
        sage: bool(inv1)
        False

    The * operator is concatenation of paths.  So long as one of the two operands
    is a QuiverPath the other can be a QuiverPath, an edge (tuple), or a list of
    edges.  If the two paths do not compose the result is the invalid path::

        sage: bool(p*q)
        False
        sage: p*(3, 4, 'c')
        a*b*c
        sage: [(4, 5, 'c'), (5, 1, 'd')]*p
        c*d*a*b

    The length of a path is the number of edges in that path.  The invalid path and
    trivial paths are therefore length 0::

        sage: len(p)
        2
        sage: triv = QuiverPath((1, 1))
        sage: len(triv)
        0
        sage: len(inv1)
        0

    List index and slice notation can be used to access the edges in a path.
    QuiverPaths can also be iterated over.  Trivial paths and the invalid path have
    no elements::

        sage: for x in p: print x
        (1, 2, 'a')
        (2, 3, 'b')
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
    """

    @staticmethod
    def __classcall__(cls, path=None):
        """
        Puts path into a standard form for UniqueRepresentation.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: a = QuiverPath([(1, 1), (1, 1)])
            sage: b = QuiverPath((1, 1))
            sage: a is b
            True
            sage: QuiverPath([(1, 1), (1, 2, 'a'), (2, 2), (2, 3, 'b'), (3, 3)])._path
            ((1, 2, 'a'), (2, 3, 'b'))
            sage: QuiverPath() is QuiverPath(())
            True
            sage: QuiverPath() is QuiverPath([])
            True
            sage: QuiverPath() is QuiverPath([(1, 1), (2, 2)])
            True
        """

        # For a QuiverPath we take it's tuple and can assume it is already in a
        # canonical form
        if isinstance(path, QuiverPath):
            return super(QuiverPath, cls).__classcall__(cls, path._path)

        # A tuple is assumed to be an edge, anything else is assumed to be a
        # list of edges
        if path is None:
            new_path = []
        elif isinstance(path, tuple):
            new_path = [path]
        else:
            new_path = list(path)

        # Check that each edge in the path is valid
        from sage.rings.finite_rings.integer_mod_ring import Integers
        good = True
        for x in new_path:
            if (len(x) < 2 or x[0] not in Integers() or x[1] not in Integers()
                           or len(x) == 2 and x[0] != x[1]
                           or len(x) == 3 and not isinstance(x[2], str)
                           or len(x) > 3):
                good = False
                break
        if not good:
            new_path = []

        # Delete trivial edges, and clear the path if not valid
        i = 0
        while i + 1 < len(new_path):
            if new_path[i][1] != new_path[i + 1][0]:
                new_path = []
            elif new_path[i][0] == new_path[i][1]:
                del new_path[i]
            else:
                i += 1
        if len(new_path) > 1 and new_path[-1][0] == new_path[-1][1]:
            del new_path[-1]

        return super(QuiverPath, cls).__classcall__(cls, tuple(new_path))

    def __init__(self, path):
        """
        Creates a path object.  Type QuiverPath? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: QuiverPath([(1, 1), (1, 1)]) is QuiverPath((1, 1))
            True
            sage: QuiverPath([(1, 1), (1, 2, 'a'), (2, 2), (2, 3, 'b'), (3, 3)])._path
            ((1, 2, 'a'), (2, 3, 'b'))
        """

        self._path = path

    def _repr_(self):
        """
        Default representation of a path.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: QuiverPath([(1, 2, 'a'), (1, 2, 'a')]) # indirect doctest
            sage: QuiverPath([(1, 2, 'a'), (2, 3, 'b')]) # indirect doctest
            a*b
            sage: QuiverPath((1, 1)) # indirect doctest
            e_1
        """

        if not self._path:
            return ''
        elif self._path[0][0] == self._path[0][1]:
            return 'e_{0}'.format(self._path[0][0])
        else:
            return '*'.join([e[2] for e in self._path])

    def __len__(self):
        """
        Returns the length of the path.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: len(QuiverPath([(1, 2, 'a'), (1, 2, 'a')]))
            0
            sage: len(QuiverPath([(1, 2, 'a'), (2, 3, 'b')]))
            2
            sage: len(QuiverPath((1, 1)))
            0
            sage: len(QuiverPath((1, 2, 'a')))
            1
        """

        if not self._path:
            return 0
        elif self._path[0][0] == self._path[0][1]:
            return 0
        else:
            return len(self._path)

    def __nonzero__(self):
        """
        Implements boolean values for the object.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: a = QuiverPath((1, 2, 'a'))
            sage: b = QuiverPath((2, 3, 'b'))
            sage: bool(a*b)
            True
            sage: bool(b*a)
            False
        """

        return bool(self._path)

    def __cmp__(self, other):
        """
        Comparison for QuiverPaths.

        If other is not of type QuiverPath then the comparison is made using `id`.
        If other is a QuiverPath then equality is tested using `is`.  If the
        QuiverPaths are unequal then one of the following data (listed in order of
        preferance) is unequal and used for comparison:

        - Length of the path
        - String representation of the path
        - Vertexes of the path (from initial to final)

        .. NOTES::

            This code is used by CombinatorialFreeModule to order the monomials
            when printing elements of QuiverAlgebras.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: a = QuiverPath((1, 2, 'a'))
            sage: b = QuiverPath((1, 2, 'b'))
            sage: c = QuiverPath((2, 3, 'c'))
            sage: c2 = QuiverPath((2, 4, 'c'))
            sage: a < a*b
            True
            sage: a*b < a
            False
            sage: a < b
            True
            sage: b < a
            False
            sage: a*c < a*c2
            True
            sage: a*c2 < a*c
            False
            sage: a < a
            False
        """

        # Compare id if other is not a QuiverPath
        if not isinstance(other, QuiverPath):
            if id(self) < id(other):
                return -1
            return 1

        # Check if equal
        if self is other:
            return 0

        # Compare lengths if different
        if len(self) != len(other):
            return len(other) - len(self)

        # Compare string representations if different
        if self._repr_() != other._repr_():
            if self._repr_() < other._repr_():
                return -1
            return 1

        # Compare internal tuple
        if self._path < other._path:
            return -1
        return 1

    def __getitem__(self, *args):
        """
        Implements index notation.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: p = QuiverPath([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: p
            a*b*c
            sage: p[0]
            (1, 2, 'a')
            sage: p[-1]
            (3, 4, 'c')
            sage: p[1:]
            ((2, 3, 'b'), (3, 4, 'c'))
        """

        if self._path and self._path[0][0] == self._path[0][1]:
            return list().__getitem__(*args)
        else:
            return self._path.__getitem__(*args)

    def __iter__(self):
        """
        Implements iteration over the path.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: p = QuiverPath([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: for e in p: print e
            (1, 2, 'a')
            (2, 3, 'b')
            (3, 4, 'c')
        """

        # Return an iterator over an empty tuple for trivial paths, otherwise
        # return an iterator for _path as a list
        if self._path and self._path[0][0] == self._path[0][1]:
            return list().__iter__()
        else:
            return list(self._path).__iter__()

    def __mul__(self, other):
        """
        Composes two paths.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: x = QuiverPath([(1, 2, 'a'), (2, 3, 'b')])
            sage: y = QuiverPath([(3, 4, 'c'), (4, 5, 'd')])
            sage: y*x

            sage: x*y
            a*b*c*d
            sage: x*(3, 4, 'c')
            a*b*c
            sage: x*[(3, 4, 'c'), (4, 5, 'd')]
            a*b*c*d
            sage: x*6
            Traceback (most recent call last):
            ...
            TypeError: QuiverPath cannot be multiplied with <type 'sage.rings.integer.Integer'>
        """

        # Handle the invalid path
        if not self._path:
            return self

        if isinstance(other, QuiverPath):
            if not other._path:
                return other
            return QuiverPath(list(self._path) + list(other._path))
        if isinstance(other, tuple):
            return QuiverPath(list(self._path) + [other])
        if isinstance(other, list):
            return QuiverPath(list(self._path) + other)
        raise TypeError("QuiverPath cannot be multiplied with {0}".format(type(other)))

    def __rmul__(self, other):
        """
        Composes two paths.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: y = QuiverPath([(3, 4, 'c'), (4, 5, 'd')])
            sage: (2, 3, 'b')*y
            b*c*d
            sage: [(1, 2, 'a'), (2, 3, 'b')]*y
            a*b*c*d
            sage: 6*y
            Traceback (most recent call last):
            ...
            TypeError: QuiverPath cannot be multiplied with <type 'sage.rings.integer.Integer'>
        """

        # Handle the invalid path
        if not self._path:
            return self

        if isinstance(other, QuiverPath):
            if not other._path:
                return other
            return QuiverPath(list(other._path) + list(self._path))
        if isinstance(other, tuple):
            return QuiverPath([other] + list(self._path))
        if isinstance(other, list):
            return QuiverPath(other + list(self._path))
        raise TypeError("QuiverPath cannot be multiplied with {0}".format(type(other)))

    def __mod__(self, other):
        """
        Returns self with other deleted from the beginning.

        If other is not the beginning of self the result is the invalid path.  Deleting
        the trivial path at vertex v from a path that begins at v does nothing.
        Deleting it from a path that does not begin at v returns the invalid path.
        Deleting the invalid path returns the invalid path.

        TESTS::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: p = QuiverPath([(1, 2, 'a'), (2, 3, 'b')])
            sage: a = QuiverPath((1, 2, 'a'))
            sage: b = QuiverPath((2, 3, 'b'))
            sage: e1 = QuiverPath((1, 1))
            sage: e2 = QuiverPath((2, 2))
            sage: p % a
            b
            sage: p % b
            sage: p % e1
            a*b
            sage: p % e2
        """

        # Convert other to a QuiverPath
        oth = QuiverPath(other)

        # Handle invalid path and other == self
        if not self._path:
            return self
        if not oth._path:
            return oth
        if self is oth:
            v = self._path[-1][1]
            return QuiverPath((v, v))

        # Handle trivial paths
        if oth._path[0][0] == oth._path[0][1]:
            if self._path[0][0] == oth._path[0][0]:
                return self
            else:
                return QuiverPath()

        # If other is the beginning, return the rest
        if self._path[:len(oth._path)] == oth._path:
            return QuiverPath(list(self._path[len(oth._path):]))
        else:
            return QuiverPath()

    def initial_vertex(self):
        """
        Returns the initial vertex of the path.

        The invalid path does not have an initial vertex, so None is returned.

        OUTPUT:

        - integer or None

        EXAMPLES::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: y = QuiverPath([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.initial_vertex()
            1
            sage: print QuiverPath([(1, 1), (2, 2)]).initial_vertex()
            None
        """

        if self._path:
            return self._path[0][0]
        else:
            return None

    def terminal_vertex(self):
        """
        Returns the terminal vertex of the path.

        The invalid path does not have an terminal vertex, so None is returned.

        OUTPUT:

        - integer or None

        EXAMPLES::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: y = QuiverPath([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.terminal_vertex()
            3
            sage: print QuiverPath([(1, 1), (2, 2)]).terminal_vertex()
            None
        """

        if self._path:
            return self._path[-1][1]
        else:
            return None

    def reverse(self):
        """
        Returns the path along the same edges in the opposite quiver.

        EXAMPLES::

            sage: from sage.modules.quiver_module import QuiverPath
            sage: p = QuiverPath([(1, 2, 'a'), (2, 3, 'b')])
            sage: p.reverse()
            b*a
        """

        # Handle the invalid path and trivial paths
        if not self._path or self._path[0][0] == self._path[0][1]:
            return self

        # Reverse all the edges in the path, then reverse the path
        new_path = [(e[1], e[0], e[2]) for e in self._path]
        return QuiverPath(reversed(new_path))

class QuiverRep_generic(Module):
    """
    This function should not be called by the user.

    Call QuiverRep with option='values' (which is the default) instead.

    INPUT:

    - ``k`` - ring, the base ring of the representation

    - ``Q`` - Quiver, the quiver of the representation

    - ``spaces`` - dict (default: empty), a dictionary associating to each vertex a
      free module over the base ring k.  Not all vertices must be specified,
      unspecified vertices are automatically set to k^0.  Keys of the dictionary
      that don't correspond to vertices are ignored.

    - ``maps`` - dict (default: empty), a dictionary associating to each edge a map
      whose domain and codomain are the spaces associated to the initial and
      terminal vertex of the edge respectively.  Not all edges must be specified,
      unspecified edges are automatically set to the zero map.  Keys of the
      dictionary that don't correspond to edges are ignored.

    OUTPUT:

    - QuiverRep

    TESTS::

        sage: from sage.modules.quiver_module import Quiver, QuiverRep_generic
        sage: Q = Quiver({1:{3:['a']}, 2:{3:['b']}})
        sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
        sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
        sage: M = QuiverRep_generic(QQ, Q, spaces, maps)
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, k, Q, spaces, maps):
        """
        Type QuiverRep? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a']}})
            sage: M = QuiverRep(GF(5), Q)
        """
        # This class can handle representations over
        # an arbitrary base ring, not necessarily a field, so long as sage can
        # construct free modules over that ring.  The data of a representation is held
        # in the following private variables:
        #
        # * _quiver
        #      The quiver of the representation.
        # * _base
        #       The quiver algebra of _quiver over _base_ring
        # * _base_ring
        #      The base ring of the representation.
        # * _spaces
        #      A dictionary which associates to each vertex of the quiver a free
        #      module over the base ring.
        # * _maps
        #      A dictionary which associates to each edge of the quiver a homomorphism
        #      whose domain and codomain equal the initial and terminal vertex of the
        #      edge.

        self._quiver = Q
        self._base_ring = k
        self._spaces = {}
        self._maps = {}
        self.Element = QuiverRepElement
        self._element_constructor_ = QuiverRepElement

        # If the vertex is not specified set it as a free module of rank 0, if
        # an integer is given set it as a free module of that rank, otherwise
        # assume the object is a module and assign it to the vertex.
        from sage.rings.finite_rings.integer_mod_ring import Integers
        for x in Q:
            if x not in spaces:
                self._spaces[x] = k**0
            elif spaces[x] in Integers():
                self._spaces[x] = k**spaces[x]
            else:
                self._spaces[x] = spaces[x]

        # The preferred method of specifying an edge is as a tuple (i, t, l)
        # where i is the initial vertex, t is the terminal vertex, and l is the
        # label.  This is the form in which quiver.edges() and other such
        # functions give the edge.  But here edges can be specified by giving
        # only the two vertices or giving only the edge label.
        for x in Q.edges():
            if x in maps:
                e = maps[x]
            elif (x[0], x[1]) in maps:
                e = maps[(x[0], x[1])]
            elif x[2] in maps:
                e = maps[x[2]]
            else:
                e = self._spaces[x[0]].Hom(self._spaces[x[1]]).zero_element()

            #If a morphism is specified use it, otherwise assume the hom
            # function can convert the object to a morphism.  Matrices and the
            # zero and one of the base ring are valid inputs (one is valid only
            # when the domain and codomain are equal.
            from sage.categories.morphism import Morphism
            if isinstance(e, Morphism):
                self._maps[x] = e
            else:
                self._maps[x] = self._spaces[x[0]].hom(e, self._spaces[x[1]])

        self._assert_valid_quiverrep()

        super(QuiverRep_generic, self).__init__(QuiverAlgebra(k, Q))

    def _assert_valid_quiverrep(self):
        """
        Raises an error if the representation is not well defined.

        Specifically it checks the map assigned to each edge.  The domain and codomain
        must equal the modules assigned to the initial and terminal vertices of the
        edge.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a']}})
            sage: M = Q.P(GF(3), 2) # indirect doctest

        Due to unique representation, we will cause bugs in later code if we modify M
        to be an invalid representation.  So we make sure to store the original values
        and replace them after testing::

            sage: sv = M._spaces[1]
            sage: M._spaces[1] = 0
            sage: M._assert_valid_quiverrep()
            Traceback (most recent call last):
            ...
            ValueError: Domain of map at edge 'a' does not match.
            sage: M._spaces[1] = sv
            sage: M = Q.P(GF(3), 2)
            sage: sv = M._maps[(1, 2, 'a')]
            sage: M._maps[(1, 2, 'a')] = (QQ^2).Hom(QQ^1).zero_element()
            sage: M._assert_valid_quiverrep()
            Traceback (most recent call last):
            ...
            ValueError: Domain of map at edge 'a' does not match.
            sage: M._maps[(1, 2, 'a')] = sv
        """

        for x in self._quiver.edges():
            if self._maps[x].domain() != self._spaces[x[0]]:
                raise ValueError("Domain of map at edge '" + str(x[2]) + "' does not match.")
            if self._maps[x].codomain() != self._spaces[x[1]]:
                raise ValueError("Codomain of map at edge '" + str(x[2]) + "' does not match.")

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a']}})
            sage: Q.P(GF(3), 2) # indirect doctest
            Representation with dimension vector (0, 1)
        """

        return "Representation with dimension vector " + str(self.dimension_vector())

    def __div__(self, sub):
        """
        This and __Truediv__ below together overload the / operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a']}})
            sage: P = Q.P(GF(3), 1)
            sage: R = P.radical()
            sage: (P/R).is_simple()
            True
        """

        return self.quotient(sub)

    def __truediv__(self, sub):
        """
        This and __Truediv__ below together overload the / operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a']}})
            sage: P = Q.P(GF(3), 1)
            sage: R = P.radical()
            sage: (P/R).is_simple()
            True
        """

        return self.quotient(sub)

    def __contains__(self, element):
        """
        This overloads the 'in' operator

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a']}})
            sage: P = Q.P(GF(3), 1)
            sage: I = Q.I(QQ, 1)
            sage: P.an_element() in P
            True
            sage: I.an_element() in P
            False
        """

        # Representations only contain elements from the same quiver
        if not isinstance(element, QuiverRepElement) or self._quiver != element._quiver:
            return False

        # The element is in the representation if and only if the element at
        # each vertex is in the space assigned to that vertex
        for v in self._quiver:
            if element._elems[v] not in self._spaces[v]:
                return False

        return True

    def _submodule(self, spaces={}):
        """
        Returns the submodule specified by the data.

        This differs from self.submodule in that it assumes the data correctly
        specifies a submodule whereas self.submodule returns the smallest submodule
        containing the data.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{3:['a']}, 2:{3:['b']}})
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: v = M.an_element()
            sage: M.submodule([v]) # indirect doctest
            Representation with dimension vector (1, 1, 1)
            sage: M.submodule(spaces={1: QQ^2}) # indirect doctest
            Representation with dimension vector (2, 0, 2)
            sage: M.submodule(spaces={2: QQ^3}) # indirect doctest
            Representation with dimension vector (0, 3, 1)
            sage: v.support()
            [1, 2, 3]
            sage: M.submodule([v], {2: QQ^3}) # indirect doctest
            Representation with dimension vector (1, 3, 1)
            sage: M.submodule().is_zero() # indirect doctest
            True
            sage: M.submodule(M.gens()) is M # indirect doctest
            True
        """

        # Add zero submodules
        for v in self._quiver:
            if v not in spaces:
                spaces[v] = self._spaces[v].zero_submodule()

        # Create edge homomorphisms restricted to the new domains and codomains
        maps = {}
        for e in self._quiver.edges():
            maps[e] = self._maps[e].restrict_domain(spaces[e[0]]).restrict_codomain(spaces[e[1]])

        return QuiverRep(self._base_ring, self._quiver, spaces, maps)

    def _coerce_map_from_(self, domain):
        """
        Returns either a QuiverRepHom from domain to self, or False.

        .. NOTES::

            This function simply tries to coerce a map at each vertex and then check if
            the result is a valid homomorphism.  If it is then that homomorphism is
            returned.  If it is not or if no coercion was possible then it returns
            False.

        INPUT:

        - ``domain`` - a Sage object

        OUTPUT:

        - QuiverRepHom or bool

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = Q.P(QQ, 1)
            sage: S = M.radical()
            sage: M.coerce_map_from(S) # indirect doctest
            Homomorphism of representations of Quiver on 3 vertices
            sage: (M/S).coerce_map_from(M) # indirect doctest
            Homomorphism of representations of Quiver on 3 vertices

        Here sage coerces a map but the result is not a homomorphism::

            sage: S.coerce_map_from(M) # indirect doctest

        Here sage cannot coerce a map::

            sage: N = Q.P(QQ, 3)
            sage: N.coerce_map_from(M) # indirect doctest
        """

        # Domain must be a QuiverRep
        if not isinstance(domain, QuiverRep_generic):
            return False

        # Coerce a map at each vertex, return false if it fails
        maps = {}
        for v in self._quiver:
            maps[v] = self._spaces[v].coerce_map_from(domain._spaces[v])
            if maps[v]==None or maps[v]==False:
                return False

        # Create and return the hom, return False if it wasn't valid
        try:
            return QuiverRepHom(domain, self, maps)
        except ValueError:
            return False

    def _Hom_(self, codomain, category):
        """
        This function is used by the coercion model.

        INPUT:

        - ``codomain`` - QuiverRepHom

        - ``category`` - This input is ignored

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: from sage.categories.homset import Hom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c', 'd']}})
            sage: P = Q.P(GF(3), 2)
            sage: S = P/P.radical()
            sage: Hom(P, S) # indirect doctest
            Dimension 1 QuiverHomSpace
        """

        if not isinstance(codomain, QuiverRep_generic):
            raise TypeError("codomain must be of type QuiverRep_generic")

        return QuiverHomSpace(self, codomain)

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view the representation data.            #
    #                                                                         #
    ###########################################################################

    def get_space(self, vertex):
        """
        Returns the module associated to the given vertex.

        INPUT:

        - ``vertex`` - integer, a vertex of the quiver of the module

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a'], 3:['b']}})
            sage: Q.P(QQ, 1).get_space(1)
            Vector space of dimension 1 over Rational Field
        """

        return self._spaces[vertex]

    def get_map(self, edge):
        """
        Returns the map associated to the given edge.

        INPUT:

        - ``edge`` - tuple of the form (initial vertex, terminal vertex, label) specifying
          the edge whose map is returned

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: Q.P(ZZ, 1).get_map((1, 2, 'a'))
            Free module morphism defined by the matrix
            [1 0]
            Domain: Ambient free module of rank 1 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
        """

        return self._maps[edge]

    def quiver(self):
        """
        Returns the quiver of the representation.

        OUTPUT:

        - Quiver

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a']}})
            sage: M = QuiverRep(GF(5), Q)
            sage: M.quiver() is Q
            True
        """

        return self._quiver

    def base_ring(self):
        """
        Returns the base ring of the representation.

        OUTPUT:

        - ring

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a']}})
            sage: M = QuiverRep(GF(5), Q)
            sage: M.base_ring() is GF(5)
            True
        """

        return self._base_ring

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the representation.       #
    #                                                                         #
    ###########################################################################

    def dimension(self, vertex=None):
        """
        Returns the dimension of the space associated to the given vertex.

        INPUT:

        - ``vertex`` - integer or None (default: None), the given vertex

        OUTPUT:

        - integer, the dimension over the base ring of the space associated to the
          given vertex.  If vertex=None then the dimension over the base ring of the
          module is returned

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: P = Q.P(GF(2), 1)
            sage: P.dimension(1)
            1
            sage: P.dimension(2)
            2

        The total dimension of the module is the sum of the dimensions at each vertex::

            sage: P.dimension()
            3
        """

        if vertex == None:
            # Sum the dimensions of each vertex
            dim = 0
            for x in self._quiver:
                dim += self._spaces[x].dimension()
            return dim
        else:
            # Return the dimension of just the one vertex
            return self._spaces[vertex].dimension()

    def dimension_vector(self):
        """
        Returns the dimension vector of the representation.

        OUTPUT:

        - tuple

        .. NOTE::

            The order of the entries in the tuple matches the order given by calling
            the vertices() method on the quiver.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: P = Q.P(GF(2), 1)
            sage: P.dimension_vector()
            (1, 2)

        Each coordinate of the dimension vector is the dimension of the space
        associated to that coordinate::

            sage: P.get_space(2).dimension()
            2
        """

        return tuple(self._spaces[x].dimension() for x in self._quiver)

    def is_zero(self):
        """
        Tests whether the representation is zero.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = QuiverRep(ZZ, Q)
            sage: N = QuiverRep(ZZ, Q, {1: 1})
            sage: M
            Representation with dimension vector (0, 0)
            sage: N
            Representation with dimension vector (1, 0)
            sage: M.is_zero()
            True
            sage: N.is_zero()
            False
        """

        return self.dimension() == 0

    def is_simple(self):
        """
        Tests whether the representation is simple.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: Q.P(RR, 1).is_simple()
            False
            sage: Q.S(RR, 1).is_simple()
            True
        """

        # A module for an acyclic quiver is simple if and only if it has total
        # dimension 1.
        return self.dimension() == 1

    def is_semisimple(self):
        """
        Tests whether the representation is semisimple.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = Q.P(QQ, 1)
            sage: (M/M.radical()).is_semisimple()
            True
        """

        # A quiver representation is semisimple if and only if the zero map is
        # assigned to each edge.
        for x in self._quiver.edges():
            if not self._maps[x].is_zero():
                return False
        return True

    def an_element(self):
        """
        Returns an element of self.

        OUTPUT:

        - QuiverRepElement

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = Q.P(QQ, 1)
            sage: M.an_element()
            Element of quiver representation
        """

        # Here we just the the an_element function from each space.
        elements = dict((v, self._spaces[v].an_element()) for v in self._quiver)
        return QuiverRepElement(self, elements)

    def zero(self):
        """
        Returns the zero element.

        OUTPUT:

        - QuiverRepElement

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = Q.P(QQ, 1)
            sage: M.zero().is_zero()
            True
        """

        # If we don't specify elements this constructor automatically returns
        # the zero element.
        return QuiverRepElement(self)

    def support(self):
        """
        Returns the support of self as a list.

        OUTPUT:

        - list, the vertices of the representation that have nonzero spaces associated
          to them

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a']}, 3:{2:['b'], 4:['c']}})
            sage: M = Q.P(QQ, 3)
            sage: M
            Representation with dimension vector (0, 1, 1, 1)
            sage: M.support()
            [2, 3, 4]
        """

        sup = []
        for v in self._quiver:
            if self._spaces[v].dimension() != 0:
                sup.append(v)

        return sup

    def gens(self, names='v'):
        """
        Returns a list of generators.

        INPUT:

        - ``names`` - an iterable variable of length equal to the number of generators or a
          string (default: 'v'), gives the names of the generators either by giving a
          name to each generator or by giving a name to which an index will be appended

        OUTPUT:

        - list of QuiverRepElement objects, the linear generators of the module (over
          the base ring)

        .. NOTES::

            The generators are ordered first by vertex and then by the order given by
            the gens() method of the space associated to that vertex.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = Q.P(QQ, 1)
            sage: M.gens()
            [v_0, v_1, v_2]

        If a string is given then it is used as the name of each generator, with the
        index of the generator appended in order to differentiate them::

            sage: M.gens('generator')
            [generator_0, generator_1, generator_2]

        If a list or other iterable variable is given then each generator is named
        using the appropriate entry.  The length of the variable must equal the number
        of generators (the dimension of the module)::

            sage: M.gens(['w', 'x', 'y', 'z'])
            Traceback (most recent call last):
            ...
            TypeError: can only concatenate list (not "str") to list
            sage: M.gens(['x', 'y', 'z'])
            [x, y, z]

        Strings are iterable, so if the length of the string is equal to the number of
        generators then the characters of the string will be used as the names::

            sage: M.gens('xyz')
            [x, y, z]
        """

        # Use names as a list if and only if it is the correct length
        uselist = (len(names) == self.dimension())
        i = 0

        # Create bases for each space and construct QuiverRepElements from
        # them.  Other functions depend on the ordering of generators produced
        # by this code.  Make sure you know which they are before you change
        # anything
        basis = []
        for v in self._quiver:
            for m in self._spaces[v].gens():
                if uselist:
                    basis.append(QuiverRepElement(self, {v: m}, names[i]))
                else:
                    basis.append(QuiverRepElement(self, {v: m}, names + "_" + str(i)))
                i += 1

        return basis

    def coordinates(self, vector):
        """
        Returns the coordinates when vector is expressed in terms of the gens.

        INPUT:

        - ``vector`` - QuiverRepElement

        OUTPUT:

        - list, the coefficients when the vector is expressed as a linear combination
          of the generators of the module

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = Q.P(QQ, 1)
            sage: x, y, z = M.gens('xyz')
            sage: M.coordinates(x - y + z)
            [1, -1, 1]
            sage: M.coordinates(M.an_element())
            [1, 1, 0]
            sage: M.an_element() == x + y
            True
        """

        # Just use the coordinates functions on each space
        coords = []
        for v in self._quiver:
            coords += self._spaces[v].coordinates(vector._elems[v])

        return coords

    def linear_combination_of_basis(self, coordinates):
        """
        Return the linear combination of the basis for self given by ``coordinates``.

        INPUT:

        - ``coordinates`` - list, a list whose length is the dimension of self.  The ith
          element of this list defines the coefficient of the ith basis vector in the
          linear combination.

        OUTPUT:

        - QuiverRepElement

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = Q.P(QQ, 1)
            sage: x, y, z = M.gens('xyz')
            sage: e = x - y + z
            sage: M.coordinates(e)
            [1, -1, 1]
            sage: M.linear_combination_of_basis([1, -1, 1]) == e
            True
        """

        # Make sure the input is valid
        gens = self.gens()
        if len(gens) != len(coordinates):
            raise ValueError("The coordinates do not match the dimension of the module.")

        result = QuiverRepElement(self)
        for i in range(0, len(gens)):
            result += coordinates[i]*gens[i]

        return result

    ###########################################################################
    #                                                                         #
    # CONSTRUCTION FUNCTIONS                                                  #
    #    These functions create and return submodules and homomorphisms.      #
    #                                                                         #
    ###########################################################################

    def submodule(self, elements=[], spaces=None):
        """
        Returns the submodule generated by the data.

        INPUT:

        - ``elements`` - a collection of QuiverRepElements (default: empty list), each
          should be an element of self

        - ``spaces`` - dictionary (default: empty), this dictionary should contain entries
          of the form {v: S} where v is a vertex of the quiver and S is a subspace of
          the vector space associated to v

        OUTPUT:

        - QuiverRep, the smallest subspace of self containing the given elements and
          the given subspaces

        .. NOTE::

            This function returns only a QuiverRep object ``sub``.  The inclusion map
            of ``sub`` into ``M``=self can be obtained by calling ``M.coerce_map_from(sub)``.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{3:['a']}, 2:{3:['b']}})
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: v = M.an_element()
            sage: M.submodule([v])
            Representation with dimension vector (1, 1, 1)

        The smallest submodule containing the vector space at vertex 1 also contains
        the image of the rank 1 homomorphism associated to the edge (1, 3, 'a')::

            sage: M.submodule(spaces={1: QQ^2})
            Representation with dimension vector (2, 0, 2)

        The smallest submodule containing the vector space at vertex 2 also contains
        the entire vector space associated to vertex 3 because there is an isomorphism
        associated to the edge (2, 3, 'b')::

            sage: M.submodule(spaces={2: QQ^3})
            Representation with dimension vector (0, 3, 1)

        As v is not already contained in this submodule adding it as a generator yields
        a larger submodule::

            sage: v.support()
            [1, 2, 3]
            sage: M.submodule([v], {2: QQ^3})
            Representation with dimension vector (1, 3, 1)

        Giving no generating data yields the zero submodule::

            sage: M.submodule().is_zero()
            True

        If the given data generates all of M then the result is M::

            sage: M.submodule(M.gens()) is M
            True
        """

        if spaces is None:
            spaces = {}

        # For each vertex generate a submodule from the given data
        dim = old_dim = 0
        for v in self._quiver:
            #Start with the zero submodule if no space is specified
            if v not in spaces:
                spaces[v] = self._spaces[v].zero_submodule()

            # Sum this with the submodule generated by the given elements.
            # Note that we are only taking the part of the element at the
            # vertex v.  We can always multiply an element of a quiver
            # representation by a constant path so we don't need to worry about
            # subspaces being embedded diagonally across multiple vertices.
            spaces[v] += self._spaces[v].submodule([m._elems[v] for m in elements])
            dim += spaces[v].dimension()

        #Now to enlarge the subspace to a submodule we sum a subspace at a
        # vertex with the images of the subspaces at adjacent vertices.  The
        # dimension of the subspace will strictly increase until we generate a
        # submodule.  At that point the dimension stabilizes and we can exit
        # the loop.
        while old_dim != dim:
            old_dim, dim = dim, 0

            # First sum the subspaces
            for e in self._quiver.edges():
                spaces[e[1]] += self._maps[e](spaces[e[0]])

            # Then get the resulting dimensions
            for v in self._quiver:
                dim += spaces[v].dimension()

        # Return self if the entire module was generated, otherwise return a
        # submodule
        if dim == self.dimension():
            return self
        else:
            return self._submodule(spaces)

    def quotient(self, sub, check=True):
        """
        Returns the quotient of self by the submodule sub.

        INPUT:

        - ``sub`` - QuiverRep, this must be a submodule of self, meaning the space
          associated to each vertex v of sub is a subspace of the space associated to v
          in self and the map associated to each edge e of sub is the restriction of
          the map associated to e in self

        - ``check`` - bool, if True then sub is checked to verify that it is indeed a
          submodule of self and an error is raised if it is not

        OUTPUT:

        - QuiverRep, the quotient module self/sub

        .. NOTE::

            This function returns only a QuiverRep object ``quot``.  The inclusion map of
            ``quot`` into ``M``=self can be obtained by calling ``M.coerce_map_from(quot)``.

        EXAMPLES:

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = Q.I(GF(3), 3)
            sage: N = Q.S(GF(3), 3)
            sage: M.quotient(N)
            Representation with dimension vector (2, 1, 0)
            sage: M.quotient(M.radical())
            Representation with dimension vector (2, 0, 0)
        """

        # First form the quotient space at each vertex
        spaces = {}
        for v in self._quiver:
            spaces[v] = self._spaces[v].quotient(sub._spaces[v], check)

        # Check the maps of sub if desired
        if check:
            for e in self._quiver.edges():
                for x in sub._spaces[e[0]].gens():
                    if sub._maps[e](x) != self._maps[e](x):
                        raise ValueError("The quotient method was not passed a submodule.")

        # Then pass the edge maps to the quotient
        maps = {}
        for e in self._quiver.edges():
            # Sage can automatically coerce an element of a module to an
            # element of a quotient of that module but not the other way
            # around.  So in order to pass a map to the quotient we need to
            # construct the quotient map for the domain so that we can take
            # inverse images to lift elments.  As sage can coerce to a quotient
            # this is easy, we just send generators to themselves and set the
            # domain to be the quotient.

            # TODO: This 'if' shouldn't need to be here, but sage crashes when
            # coercing elements into a quotient that is zero.  Once Trac ticket
            # 12413 gets fixed only the else should need to execute.
            if spaces[e[1]].dimension() == 0:
                maps[e] = spaces[e[0]].Hom(spaces[e[1]]).zero_element()
            else:
                factor = self._spaces[e[0]].hom(self._spaces[e[0]].gens(), spaces[e[0]])
                # Now we create a homomorphism by specifying the images of
                # generators.  Each generator is lifted to the original domain and
                # mapped over using the original map.  The codomain is set as the
                # quotient so sage will take care of pushing the result to the
                # quotient in the codomain.
                maps[e] = spaces[e[0]].hom([self._maps[e](factor.lift(x)) for x in spaces[e[0]].gens()], spaces[e[1]])

        return QuiverRep(self._base_ring, self._quiver, spaces, maps)

    def socle(self):
        """
        The socle of self.

        OUTPUT:

        - QuiverRep, the socle

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = Q.P(QQ, 1)
            sage: M.socle()
            Representation with dimension vector (0, 0, 2)
        """

        # The socle of a representation is the intersection of the kernels of
        # all the edge maps.  The empty intersection is defined to be the
        # entire space so this is what we start with.
        spaces = self._spaces.copy()
        for e in self._quiver.edges():
            spaces[e[0]] = spaces[e[0]].intersection(self._maps[e].kernel())

        return self._submodule(spaces)

    def radical(self):
        """
        Returns the Jacobson radical of self.

        OUTPUT:

        - QuiverRep, the socle

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = Q.P(QQ, 1)
            sage: M.radical()
            Representation with dimension vector (0, 2, 2)
        """

        #The Jacobson radical of a representation is the sum of the images of
        # all of the edge maps.  The empty sum is defined to be zero so this is
        # what we start with.
        spaces = dict((v, self._spaces[v].zero_submodule()) for v in self._quiver)
        for e in self._quiver.edges():
            spaces[e[1]] += self._maps[e].image()

        return self._submodule(spaces)

    def top(self):
        """
        Returns the top of self.

        OUTPUT:

        - QuiverRep, the quotient of self by its radical

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = Q.P(QQ, 1)
            sage: M.top()
            Representation with dimension vector (1, 0, 0)
            sage: M.top() == M/M.radical()
            True
        """

        return self.quotient(self.radical())

    def zero_submodule(self):
        """
        Returns the zero submodule.

        OUTPUT:

        - QuiverRep, the quotient of self by its radical

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = Q.P(QQ, 1)
            sage: M.zero_submodule()
            Representation with dimension vector (0, 0, 0)
            sage: M.zero_submodule().is_zero()
            True
        """

        # When no data is specified this constructor automatically returns the
        # zero submodule
        return self._submodule()

    def hom(self, codomain, maps={}):
        """
        Returns a homomorphism from self to codomain defined by maps.

        For more information see the QuiverRepHom documentation.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: f = S.hom(M)
            sage: f.is_zero()
            True
            sage: maps2 = {2:[1, -1], 3:1}
            sage: g = S.hom(M, maps2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: h = S.hom(M, [x, y])
            sage: g == h
            True
            sage: Proj = Q.P(GF(7), 3)
            sage: Simp = Q.S(GF(7), 3)
            sage: im = QuiverRepElement(Simp, {3: (1,)})
            sage: Proj.hom(Simp, im).is_surjective()
            True
        """

        return QuiverRepHom(self, codomain, maps)

    def linear_dual(self):
        """
        Computes the linear dual Hom_k(M, k) of the module M=self over the base ring k.

        OUTPUT:

        - QuiverRep, the dual representation

        .. NOTES::

            If e is an edge of the quiver Q then we let (fe)(m) = f(me).  This gives
            Hom_k(M, k) a module structure over the opposite quiver Q.reverse().

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = Q.P(QQ, 1)
            sage: M.linear_dual()
            Representation with dimension vector (1, 2, 2)
            sage: M.linear_dual().quiver() == Q.reverse()
            True
        """

        # This module is formed by taking the transpose of the edge maps.
        spaces = self._spaces.copy()
        maps = dict(((e[1], e[0], e[2]), self._spaces[e[1]].hom(self._maps[e].matrix().transpose(), self._spaces[e[0]]))
                    for e in self._quiver.edges())

        # Reverse the bases if present
        if hasattr(self, '_bases'):
            bases = {}
            basis = []
            for v in self._bases:
                bases[v] = [p.reverse() for p in self._bases[v]]
                basis.extend(bases[v])

        if isinstance(self, QuiverRep_with_path_basis):
            result = QuiverRep(self._base_ring, self._quiver.reverse(), basis, option='dual paths')
            result._maps = maps
            result._bases = bases
        elif isinstance(self, QuiverRep_with_dual_path_basis):
            result = QuiverRep(self._base_ring, self._quiver.reverse(), basis, option='paths')
            result._maps = maps
            result._bases = bases
        else:
            result = QuiverRep(self._base_ring, self._quiver.reverse(), spaces, maps)

        return result

    def algebraic_dual(self, basis=False):
        """
        Computes the algebraic dual Hom_Q(M, kQ) of the module M=self.

        INPUT:

        - ``basis`` - bool, if false then only the module is returned.  If true then a
          tuple is returned.  The first element is the QuiverRep and the second element
          is a dictionary which associates to each vertex a list.  The elements of this
          list are the homomorphisms which correspond to the basis elements of that
          vertex in the module.

        OUTPUT:

        - QuiverRep or tuple

        .. NOTES::

            Here kQ is the path algebra considered as a right module over itself.  If e
            is an edge of the quiver Q then we let (fe)(m) = ef(m).  This gives
            Hom_Q(M, kQ) a module structure over the opposite quiver Q.reverse().

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b'], 3: ['c', 'd']}, 2:{3:['e']}})
            sage: Q.free_module(GF(7)).algebraic_dual().dimension_vector()
            (7, 2, 1)
        """

        return QuiverHomSpace(self, self._quiver.free_module(self._base_ring)).left_module(basis)

    def Hom(self, codomain):
        """
        Returns the hom space from self to codomain.

        For more information see the QuiverHomSpace documentation.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: Q.S(QQ, 2).Hom(Q.P(QQ, 1))
            Dimension 2 QuiverHomSpace
        """

        return QuiverHomSpace(self, codomain)

    def direct_sum(self, modules, return_maps=False):
        """
        Returns the direct sum of self with the given modules.

        The modules must be modules over the same quiver and base ring.

        INPUT:

        - ``modules`` - QuiverRep or list of QuiverReps

        - ``return_maps`` - Boolean (default: False), if False then the output is a single
          QuiverRep object which is the direct sum of self with the given module or
          modules.  If True then the output is a list ``[sum, iota, pi]``.  The
          first entry ``sum`` is the direct sum of self with the given module or
          modules.  Both ``iota`` and ``pi`` are lists of QuiverRepHoms with one entry for
          each summand; ``iota[i]`` is the inclusion map and ``pi[i]`` is the projection
          map of the ith summand.  The summands are ordered as given with self
          being the zeroth summand.

        OUTPUT:

        - QuiverRep or tuple

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{4:['c']}, 3:{4:['d']}})
            sage: P1 = Q.P(QQ, 1)
            sage: P2 = Q.P(QQ, 2)
            sage: S = P1.direct_sum(P2)
            sage: P1.dimension_vector()
            (1, 1, 1, 2)
            sage: P2.dimension_vector()
            (0, 1, 0, 1)
            sage: S.dimension_vector()
            (1, 2, 1, 3)
            sage: S, iota, pi = P1.direct_sum(P2, return_maps=True)
            sage: iota[0].domain() is P1
            True
            sage: iota[1].domain() is P2
            True
            sage: pi[0].codomain() is P1
            True
            sage: pi[1].codomain() is P2
            True
            sage: iota[0].codomain() is S
            True
            sage: iota[1].codomain() is S
            True
            sage: pi[0].domain() is S
            True
            sage: pi[1].domain() is S
            True
            sage: iota[0].get_matrix(4)
            [1 0 0]
            [0 1 0]
            sage: pi[0].get_matrix(4)
            [1 0]
            [0 1]
            [0 0]
            sage: P1prime = S/iota[1].image()
            sage: f = P1prime.coerce_map_from(S)
            sage: g = f*iota[0]
            sage: g.is_isomorphism()
            True
        """

        # Convert a single module into a length 1 list and check validity
        if isinstance(modules, QuiverRep_generic):
            mods = [self, modules]
        else:
            mods = [self] + list(modules)
        for i in range(1, len(mods)):
            if not isinstance(mods[i], QuiverRep_generic):
                raise ValueError("modules must be a QuiverRep or a list of QuiverReps")
            if self._quiver is not mods[i]._quiver:
                raise ValueError("Cannot direct sum modules of different quivers")
            if self._base_ring is not mods[i]._base_ring:
                raise ValueError("Cannot direct sum modules with different base rings")

        # Get the dimensions of all spaces at each vertex
        dims = dict((v, [x._spaces[v].dimension() for x in mods]) for v in self._quiver)

        # Create spaces of the correct dimensions
        spaces = dict((v, self._base_ring**sum(dims[v])) for v in self._quiver)

        # Take block sums of matrices to form the maps
        from sage.matrix.constructor import block_diagonal_matrix
        maps = {}
        for e in self._quiver.edges():
            maps[e] = block_diagonal_matrix([x._maps[e].matrix() for x in mods], subdivide=False)

        # Create the QuiverRep, return if the maps aren't wanted
        result = QuiverRep(self._base_ring, self._quiver, spaces, maps)
        if not return_maps:
            return result

        # Create the inclusions and projections
        from sage.matrix.constructor import block_matrix, Matrix
        from sage.rings.integer import Integer
        iota = []
        pi = []
        for i in range(0, len(mods)):
            incl_maps = {}
            proj_maps = {}
            for v in self._quiver:
                # Create the maps using block matrices
                pre_dims = sum(dims[v][:i])
                post_dims = sum(dims[v][i + 1:])
                incl_maps[v] = block_matrix(1, 3, [Matrix(dims[v][i], pre_dims),
                                                   Matrix(dims[v][i], dims[v][i], Integer(1)),
                                                   Matrix(dims[v][i], post_dims)])
                proj_maps[v] = block_matrix(3, 1, [Matrix(pre_dims, dims[v][i]),
                                                   Matrix(dims[v][i], dims[v][i], Integer(1)),
                                                   Matrix(post_dims, dims[v][i])])

            # Create and save the QuiverRepHom
            iota.append(QuiverRepHom(mods[i], result, incl_maps))
            pi.append(QuiverRepHom(result, mods[i], proj_maps))

        # Return all the data
        return [result, iota, pi]

    def projective_cover(self, return_maps=False):
        """
        Return the projective cover of self.

        EXAMPLES::

            sage: from sage.modules.quiver_module import *
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c','d']}})
            sage: S1 = Q.S(GF(3), 1)
            sage: f1 = S1.projective_cover()
            sage: f1.is_surjective()
            True
            sage: f1._domain
            Representation with dimension vector (1, 2, 4)
            sage: Q.P(GF(3), 1)
            Representation with dimension vector (1, 2, 4)
        """

        # Deal with the zero module
        if self.dimension() == 0:
            return self.coerce_map_from(self)

        # Get the projection onto the top and generators
        top = self.top()
        proj_to_top = top.coerce_map_from(self)
        gens = top.gens()

        # Lift each of the generators
        lifts = [proj_to_top.lift(x) for x in gens]

        # Get projective summands of the cover
        Ps = [self._quiver.P(self._base_ring, x.support()[0]) for x in gens]

        # Factor the maps through self
        maps = [Ps[i].hom(self, lifts[i]) for i in range(0, len(gens))]

        # Sum them and return
        return maps[0].direct_sum(maps[1:], return_maps, 'codomain')

    def transpose(self):
        """
        Returns the transpose of self.

        The transpose, Tr M, of a module M is defined as follows.  Let p: P1 -> P2 be
        the second map in a minimal projective presentation P1 -> P2 -> M -> 0 of M.
        If p^t is the algebraic dual of p then define Tr M = coker p^t.

        OUTPUT:

        - QuiverRep

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = QuiverRep(GF(3), Q, {1: 1, 2: 1}, {(1, 2, 'a'): 1})
            sage: M.transpose()
            Representation with dimension vector (1, 1)
        """

        # Get the second map in the minimal projective resolution of self.
        # Note we create this map as the projective cover of a kernel so we
        # need to alter the codomain to be the projective and not the kernel.
        p1 = self.projective_cover()
        k = p1.kernel()
        p = p1.domain().coerce_map_from(k)*k.projective_cover()

        # Return the cokernel
        return p.algebraic_dual().cokernel()

    def AR_translate(self):
        """
        Return the Auslander-Reiten translate of self.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: M = QuiverRep(GF(3), Q, {1: 1, 2: 1}, {(1, 2, 'a'): 1})
            sage: tauM = M.AR_translate()
            sage: tauM
            Representation with dimension vector (1, 1)
            sage: tauM.get_map((1, 2, 'a')).matrix()
            [1]
            sage: tauM.get_map((1, 2, 'b')).matrix()
            [0]

        The module M above is its own AR translate.  This is not always true::

            sage: Q2 = Quiver({3:{1:['b']}, 5:{3:['a']}})
            sage: Q2.S(QQ, 5).AR_translate()
            Representation with dimension vector (0, 1, 0)
        """

        return self.transpose().linear_dual()

    ###########################################################################
    #                                                                         #
    # ADDITIONAL OPERATIONS                                                   #
    #    These functions operations that are not implemented via binary       #
    #    operators.                                                           #
    #                                                                         #
    ###########################################################################

    def right_edge_action(self, element, path):
        """
        Returns the result of element*edge.

        INPUT:

        - ``element`` - QuiverRepElement, an element of self

        - ``path`` - QuiverPath or list of tuples

        OUTPUT:

        - QuiverRepElement, the result of element*path when path is considered an
          element of the path algebra of the quiver

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = Q.P(QQ, 1)
            sage: v = M.an_element()
            sage: v.support()
            [1, 2, 3]
            sage: M.right_edge_action(v, (1, 1)).support()
            [1]
            sage: M.right_edge_action(v, [(1, 1)]).support()
            [1]
            sage: M.right_edge_action(v, [(1, 1), (2, 2)]).support()
            []
            sage: M.right_edge_action(v, [(1, 1), (1, 2, 'a')]).support()
            [2]
            sage: M.right_edge_action(v, (1, 2, 'a')) == M.right_edge_action(v, [(1, 1), (1, 2, 'a'), (2, 2)])
            True
        """

        # Convert to a QuiverPath
        qpath = QuiverPath(path)

        # Cannot act by a path not in the quiver
        if qpath not in self._quiver:
            raise ValueError("{0} is not a path in the quiver.".format(path))

        # Invalid paths are zero in the quiver algebra
        result = self.zero()
        if not qpath:
            return result

        # Start with the element at the initial vertex
        x = element._elems[qpath.initial_vertex()]

        # Act by each edge
        for e in qpath:
            x = self.get_map(e)(x)

        # Assign the result to the terminal vertex and return
        result._elems[qpath.terminal_vertex()] = x
        return result

class QuiverRep_with_path_basis(QuiverRep_generic):
    """
    The basis of the module must be closed under right multiplication by an
    edge; that is, appending any edge to the end of any path in the basis must
    result in either an invalid path or a valid path also contained in the
    basis of the module.

    INPUT:

    - ``k`` - ring, the base ring of the representation

    - ``Q`` - Quiver, the quiver of the representation

    - ``basis`` - list (default: empty), should be a list of paths (also lists) in the
      quiver Q.  Entries that do not represent valid paths are ignored and duplicate
      paths are deleted.  The closure of this list under right multiplication forms the
      basis of the resulting representation.
    """
    # This class implements quiver representations whose bases correspond to
    # paths in the path algebra and whose maps are path multiplication.  The
    # main advantage to having such a basis is that a homomorphism can be
    # defined by giving a single element in the codomain.  This class derives
    # from the QuiverRep class and the following private methods and variables
    # have been added:
    #
    # * _bases
    #       A dictionary associating to each vertex a list of paths (also lists)
    #       which correspond to the basis elements of the space assigned to that
    #       vertex.
    #
    # If the right closure of the basis is also closed under left mult by an
    # edge then the object will have the following private variable:
    #
    # * _left_action_mats
    #       A dictionary of dictionaries.  If e is an edge or trivial path of
    #       the quiver then _left_action_mats[e] is a dictionary that
    #       associates a matrix to each vertex.  The left action of e on a
    #       QuiverRepElement is given by multiplying, on the left, the vector
    #       associated to vertex v in the QuiverRepElement by the matrix
    #       _left_action_mats[e][v]

    def __init__(self, k, Q, basis):
        """
        Type QuiverRep_with_path_basis? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q1 = Quiver({1:{2:['a']}})
            sage: P1 = QuiverRep(QQ, Q1, [[(1, 1)]], option='paths')
            sage: P1.dimension()
            2
            sage: kQ = QuiverRep(QQ, Q1, [[(1, 1)], [(2, 2)], [(1, 1), (1, 2, 'a'), (2, 2)], [(1, 2, 'a')]], option='paths')
            sage: kQ.dimension()
            3
            sage: Q2 = Quiver({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}})
            sage: M = QuiverRep(QQ, Q2, [[(2, 2)], [(1, 2, 'a')]], option='paths')
            sage: M.dimension_vector()
            (0, 2, 2)
            sage: N = QuiverRep(QQ, Q2, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
            sage: N.dimension_vector()
            (0, 1, 2)
        """

        self._quiver = Q
        self._base_ring = k

        # Add the paths to the basis dictionary.  The terminal vertex is the
        # key
        self._bases = dict((v, []) for v in Q)
        for path in basis:
            self._bases[path.terminal_vertex()].append(path)

        # Create the matrixes of the maps
        from sage.matrix.constructor import Matrix
        maps = {}
        for e in Q.edges():
            # Start with the zero matrix and fill in from there
            maps[e] = Matrix(self._base_ring, len(self._bases[e[0]]), len(self._bases[e[1]]))
            for i in range(0, len(self._bases[e[0]])):
                # Add an entry to the matrix coresponding to where the new path is found
                j = self._bases[e[1]].index(self._bases[e[0]][i]*e)
                maps[e][i, j] = self._base_ring(1)

        # Create the spaces and then the representation
        spaces = dict((v, len(self._bases[v])) for v in Q)
        super(QuiverRep_with_path_basis, self).__init__(k, Q, spaces, maps)

        # Try and create the matrices for the left edge action of edges.  If it
        # fails just return, there's no edge action and the construction is
        # done
        action_mats = {}
        for e in self._quiver.edges():
            action_mats[e] = {}
            for v in self._quiver:
                # Start with the zero matrix and fill in
                l = len(self._bases[v])
                action_mats[e][v] = Matrix(self._base_ring, l, l)

                for j in range(0, l):
                    if e[1] == self._bases[v][j].initial_vertex():
                        try:
                            action_mats[e][v][self._bases[v].index(e*self._bases[v][j]), j] = self._base_ring(1)
                        except ValueError:
                            # There is no left action
                            return

        # If haven't returned yet then there is a left action.  Create the
        # matrices for acting by trivial paths
        for vert in self._quiver:
            e = (vert, vert)
            action_mats[e] = {}
            for v in self._quiver:
                # Start with the zero matrix and fill in
                l = len(self._bases[v])
                action_mats[e][v] = Matrix(self._base_ring, l, l)

                # Paths not beginning at vert are sent to zero, paths beginning
                # at vert are fixed
                for i in range(0, l):
                    if self._bases[v][i].initial_vertex() == vert:
                        action_mats[e][v][i, i] = self._base_ring(1)

        # Define the method and save the matrices
        self.left_edge_action = self._left_edge_action
        self._left_action_mats = action_mats

    def _left_edge_action(self, edge, element):
        """
        Returns the result of edge*element.

        INPUT:

        - ``element`` - QuiverRepElement, an element of self

        - ``edge`` - an edge of the quiver (a tuple) or a list of edges in the quiver.
          Such a list can be empty (in which case no action is performed) and can
          contain trivial paths (tuples of the form (v, v) where v is a vertex of the
          quiver)

        OUTPUT:

        - QuiverRepElement, the result of edge*element when edge is considered an
          element of the path algebra of the quiver

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c']}})
            sage: M = QuiverRep(QQ, Q, [[(1, 1)], [(2, 2)], [(3, 3)]], option='paths')
            sage: v = M.an_element()
            sage: v.support()
            [1, 2, 3]

        The sum of all trivial paths is the identity element under this action::

            sage: x = M.left_edge_action((1, 1), v) # indirect doctest
            sage: y = M.left_edge_action((2, 2), v) # indirect doctest
            sage: z = M.left_edge_action((3, 3), v) # indirect doctest
            sage: x + y + z == v
            True

        Note that the action only depends on the element of the path algebra
        that the input specifies::

            sage: a = M.left_edge_action([(1, 1), (1, 2, 'a'), (2, 2)], v) # indirect doctest
            sage: b = M.left_edge_action((1, 2, 'a'), v) # indirect doctest
            sage: a == b
            True

        These two edges multiply to zero in the path algebra::

            sage: M.left_edge_action([(1, 2, 'a'), (1, 2, 'b')], v).is_zero() # indirect doctest
            True
        """

        # Deal with lists by calling this function recursively
        if isinstance(edge, list):
            if len(edge) == 0:
                return element;
            else:
                return self.left_edge_action(edge[:-1], self.left_edge_action(edge[-1], element))

        # Now we are just acting by a single edge
        elems = dict((v, self._left_action_mats[edge][v]*element._elems[v]) for v in self._quiver)
        return QuiverRepElement(self, elems)

    def is_left_module(self):
        """
        Tests whether the basis is closed under left multiplication.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q1 = Quiver({1:{2:['a']}})
            sage: P2 = QuiverRep(QQ, Q1, [[(2, 2)]], option='paths')
            sage: P2.is_left_module()
            False

        The supplied basis is not closed under left multiplication, but it's not closed
        under right multiplication either.  When the closure under right multiplication
        is taken the result is also closed under left multiplication and therefore
        produces a left module structure::

            sage: kQ = QuiverRep(QQ, Q1, [[(1, 1)], [(2, 2)]], option='paths')
            sage: kQ.is_left_module()
            True

        Taking the right closure of a left closed set produces another left closed set::

            sage: Q2 = Quiver({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}})
            sage: M = QuiverRep(QQ, Q2, [[(2, 2)], [(1, 2, 'a')]], option='paths')
            sage: M.is_left_module()
            True

        Note that the second path is length 2, so even though the edge (1, 2, 'a')
        appears in the input the path [(1, 2, 'a')] is not in the right closure::

            sage: N = QuiverRep(QQ, Q2, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
            sage: N.is_left_module()
            False
        """

        return hasattr(self, 'left_edge_action')

class QuiverRep_with_dual_path_basis(QuiverRep_generic):
    """
    The basis of the module must be closed under left deletion of an edge; that
    is, deleting any edge from the beginning of any path in the basis must
    result in a path also contained in the basis of the module.

    INPUT:

    - ``k`` - ring, the base ring of the representation

    - ``Q`` - Quiver, the quiver of the representation

    - ``basis`` - list (default: empty), should be a list of paths (also lists) in the
      quiver Q.  Entries that do not represent valid paths are ignored and duplicate
      paths are deleted.  The closure of this list under left deletion forms the
      basis of the resulting representation.
    """
    # This class implements quiver representations whose bases correspond to
    # paths in the path algebra and whose maps are edge deletion.  The
    # main advantage to having such a basis is that a homomorphism can be
    # defined by giving a single element in the domain.  This class derives
    # from the QuiverRep class and the following methods and private variables
    # have been added:
    #
    # * _bases
    #       A dictionary associating to each vertex a list of paths which
    #       correspond to the basis elements of the space assigned to that
    #       vertex.

    def __init__(self, k, Q, basis):
        """
        Type QuiverRep_with_dual_path_basis? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep
            sage: Q1 = Quiver({1:{2:['a']}})
            sage: I2 = QuiverRep(QQ, Q1, [(2, 2)], option='dual paths')
            sage: I2.dimension()
            2
            sage: kQdual = QuiverRep(QQ, Q1, [[(1, 1)], [(2, 2)], [(1, 1), (1, 2, 'a'), (2, 2)], [(1, 2, 'a')]], option='dual paths')
            sage: kQdual.dimension()
            3
            sage: Q2 = Quiver({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}})
            sage: M = QuiverRep(QQ, Q2, [[(1, 2, 'a'), (2, 3, 'd')], [(1, 3, 'b')]], option='dual paths')
            sage: M.dimension_vector()
            (2, 0, 0)
            sage: N = QuiverRep(QQ, Q2, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='dual paths')
            sage: N.dimension_vector()
            (2, 1, 0)
        """

        self._quiver = Q
        self._base_ring = k

        # Add the paths to the basis dictionary.  The initial vertex is the
        # key
        self._bases = dict((v, []) for v in Q)
        for path in basis:
            self._bases[path.initial_vertex()].append(path)

        # Create the matrixes of the maps
        from sage.matrix.constructor import Matrix
        maps = {}
        for e in Q.edges():
            # Start with the zero matrix and fill in from there
            maps[e] = Matrix(self._base_ring, len(self._bases[e[0]]), len(self._bases[e[1]]))
            for i in range(0, len(self._bases[e[0]])):
                # Add an entry to the matrix coresponding to where the new path is found
                if self._bases[e[0]][i] % e in self._bases[e[1]]:
                    j = self._bases[e[1]].index(self._bases[e[0]][i] % e)
                    maps[e][i, j] = self._base_ring(1)

        # Create the spaces and then the representation
        spaces = dict((v, len(self._bases[v])) for v in Q)
        super(QuiverRep_with_dual_path_basis, self).__init__(k, Q, spaces, maps)

class QuiverHomSpace_generic(Parent):
    """
    A homomorphism of quiver representations is for each vertex of the quiver a
    homomorphism of the spaces assigned to those vertices such that these
    homomorphisms commute with the edge maps.  This class handles the set of all
    such maps, Hom_Q(M, N).

    INPUT:

    - ``domain`` - QuiverRep, the domain of the homomorphism space

    - ``codomain`` - QuiverRep, the codomain of the homomorphism space

    OUTPUT:

    - QuiverHomSpace, the homomorphism space Hom_Q(domain, codomain)

    .. NOTES::

        The quivers of the domain and codomain must be equal or a ValueError is
        raised.

    EXAMPLES::

        sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
        sage: Q = Quiver({1:{2:['a', 'b']}})
        sage: H = QuiverHomSpace(Q.S(QQ, 2), Q.P(QQ, 1))
        sage: H.dimension()
        2
        sage: H.gens()
        [Homomorphism of representations of Quiver on 2 vertices, Homomorphism of representations of Quiver on 2 vertices]
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, domain, codomain):
        """
        Type QuiverHomSpace? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: H = QuiverHomSpace(Q.S(QQ, 2), Q.P(QQ, 1))
            sage: H.dimension()
            2
            sage: H.gens()
            [Homomorphism of representations of Quiver on 2 vertices, Homomorphism of representations of Quiver on 2 vertices]
        """
        # The data in the class is stored in the following private variables:
        #
        # * _base
        #      The base ring of the representations M and N.
        # * _codomain
        #      The QuiverRep object of the codomain N.
        # * _domain
        #      The QuiverRep object of the domain M.
        # * _quiver
        #      The quiver of the representations M and N.
        # * _space
        #      A free module with ambient space.
        #
        # The free module _space is the homomorphism space.  The ambient space
        # is k^n where k is the base ring and n is the sum of the dimensions of
        # the spaces of homomorphisms between the free modules attached in M
        # and N to the vertices of the quiver.  Each coordinate represents a
        # single entry in one of those matrices.

        # Get the quiver and base ring and check they they are the same for
        # both modules
        self._quiver = domain._quiver
        self._domain = domain
        self._codomain = codomain
        if self._quiver != codomain._quiver:
            raise ValueError("Representations are not over the same quiver.")
        if codomain._base_ring != domain._base_ring:
            raise ValueError("Representations are not over the same base ring.")

        # To compute the Hom Space we set up a 'generic' homomorphism where the
        # maps at each vertex are described by matrices whose entries are
        # variables.  Then the commutativity of edge diagrams gives us a
        # system of equations whose solution space is the Hom Space we're
        # looking for.  The variables will be numbered consecutively starting
        # at 0, ordered first by the vertex the matrix occurs at, then by row
        # then by column.  We'll have to keep track of which variables
        # correspond to which matrices.

        # eqs will count the number of equations in our system of equations,
        # varstart will be a list whose ith entry is the number of the variable
        # located at (0, 0) in the matrix assigned to the ith vertex.
        eqs = 0
        verts = domain._quiver.vertices()
        varstart = [0]*(len(verts) + 1)

        # First assign to varstart the dimension of the matrix assigned to the
        # previous vertex.
        for v in verts:
            varstart[verts.index(v) + 1] = domain._spaces[v].dimension()*codomain._spaces[v].dimension()
        for e in domain._quiver.edges():
            eqs += domain._spaces[e[0]].dimension()*codomain._spaces[e[1]].dimension()

        # After this cascading sum varstart[v] will be the sum of the
        # dimensions of the matrixes assigned to vertices ordered before v.
        # This is equal to the number of the first variable assigned to v.
        for i in range(2, len(varstart)):
            varstart[i] += varstart[i-1]

        # This will be the coefficient matrix for the system of equations.  We
        # start with all zeros and will fill in as we go.  We think of this
        # matrix as acting on the right so the columns correspond to equations,
        # the rows correspond to variables, and .kernel() will give a right
        # kernel as is needed.
        from sage.matrix.constructor import Matrix
        coef_mat = Matrix(codomain._base_ring, varstart[-1], eqs)

        # row keeps track of what equation we are on.  If the maps X and Y are
        # assigned to an edge e and A and B are the matrices of variables that
        # describe the generic maps between the initial and final vertices of e
        # then commutativity of the edge diagram is described by the equation
        # AY = XB, or
        #
        #          Sum_k A_ik*Y_kj - Sum_k X_ikB_kj == 0 for all i and j.
        #
        # Below we loop through these values of i,j,k and write the
        # coefficients of the equation above into the coefficient matrix.
        eqn = 0
        for e in domain._quiver.edges():
            X = domain._maps[e].matrix()
            Y = codomain._maps[e].matrix()
            for i in range(0, X.nrows()):
                for j in range(0, Y.ncols()):
                    for k in range(0, Y.nrows()):
                        coef_mat[varstart[verts.index(e[0])] + i*Y.nrows() + k, eqn] = Y[k, j]
                    for k in range(0, X.ncols()):
                        coef_mat[varstart[verts.index(e[1])] + k*Y.ncols() + j, eqn] = -X[i, k]
                    eqn += 1

        # Now we can create the hom space
        self._space = coef_mat.kernel()

        # Bind identity if domain = codomain
        if domain is codomain:
            self.identity = self._identity

        super(QuiverHomSpace_generic, self).__init__(base=codomain._base_ring,
                                                     element_constructor=self._element_constructor_)

    def _coerce_map_from_(self, other):
        """
        A coercion exists if and only if `other` is also a
        QuiverHomSpace and there is a coercion from the domain of `self`
        to the domain of `other` and from the codomain of `other` to the
        domain of `self`.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a']}})
            sage: P = Q.P(QQ, 1)
            sage: S = Q.S(QQ, 1)
            sage: H1 = P.Hom(S)
            sage: H2 = (P/P.radical()).Hom(S)
            sage: H1.coerce_map_from(H2) # indirect doctest
            Conversion map:
                  From: Dimension 1 QuiverHomSpace
                  To:   Dimension 1 QuiverHomSpace
        """

        if not isinstance(other, QuiverHomSpace_generic):
            return False
        if not other._domain.has_coerce_map_from(self._domain):
            return False
        if not self._codomain.has_coerce_map_from(other._codomain):
            return False
        return True

    def _element_constructor_(self, data=None):
        """
        A homomorphism of quiver representations is for each vertex of the quiver a
        homomorphism of the spaces assigned to those vertices such that these
        homomorphisms commute with the edge maps.  The domain and codomain of the
        homomorphism are required to be representations over the same quiver with
        the same base ring.

        INPUT:

        - ``data`` - dict, list, QuiverRepElement or QuiverRepHom (default: empty dict)
          as follows:
          - list, data can be a list of images for the generators of the domain.  An
            error will be generated if the map so defined is not equivariant with
            respect to the action of the quiver.
          - dictionary, data can be a dictionary associating to each vertex of the
            quiver either a homomorphism with domain and codomain the spaces associated
            to this vertex in the domain and codomain modules respectively, or a matrix
            defining such a homomorphism, or an object that sage can construct such a
            matrix from.  Not all vertices must be specified, unspecified vertices are
            assigned the zero map, and keys not corresponding to vertices of the quiver
            are ignored.  An error will be generated if these maps do not commute with
            the edge maps of the domain and codomain.
          - QuiverRepElement, if the domain is a QuiverRep_with_path_basis then data
            can be a single QuiverRepElement belonging to the codomain.  The map is
            then defined by sending each path, p, in the basis to data*p.  If data is
            not an element of the codomain or the domain is not a
            QuiverRep_with_path_basis then an error will be generated.
          - QuiverRepHom, the input can also be a map ``f:D -> C`` such that there is a
            coercion from the domain of self to ``D`` and from ``C`` to the codomain of
            self.  The composition of these maps is the result.

        OUTPUT:

        - QuiverRepHom

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverHomSpace, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: H = QuiverHomSpace(S, M)

        With no additional data this creates the zero map::

            sage: f = H() # indirect doctest
            sage: f.is_zero()
            True

        We must specify maps at the vertices to get a nonzero homomorphism.  Note that
        if the dimensions of the spaces assigned to the domain and codomain of a vertex
        are equal then Sage will construct the identity matrix from ``1``::

            sage: maps2 = {2:[1, -1], 3:1}
            sage: g = H(maps2) # indirect doctest

        Here we create the same map by specifying images for the generators::

            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: h = H([x, y]) # indirect doctest
            sage: g == h
            True

        If the domain is a module of type QuiverRep_with_path_basis (for example, the
        indecomposable projectives) we can create maps by specifying a single image::

            sage: Proj = Q.P(GF(7), 3)
            sage: Simp = Q.S(GF(7), 3)
            sage: im = QuiverRepElement(Simp, {3: (1,)})
            sage: H2 = QuiverHomSpace(Proj, Simp)
            sage: H2(im).is_surjective() # indirect doctest
            True
        """

        if data is None or data == 0:
            data = {}

        return QuiverRepHom(self._domain, self._codomain, data)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a']}})
            sage: Q.P(GF(3), 2).Hom(Q.S(GF(3), 2)) # indirect doctest
            Dimension 1 QuiverHomSpace
        """

        return "Dimension " + str(self._space.dimension()) + " QuiverHomSpace"

    def __contains__(self, map):
        """
        This overloads the in operator

        TESTS::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a']}})
            sage: H1 = Q.P(GF(3), 2).Hom(Q.S(GF(3), 2))
            sage: H2 = Q.P(GF(3), 2).Hom(Q.S(GF(3), 1))
            sage: H1.an_element() in H1
            True
            sage: H2.an_element() in H1
            False
        """

        # First check the type
        if not isinstance(map, QuiverRepHom):
            return False

        # Then check the quivers, domain, and codomain
        if self._quiver != map._quiver or self._domain != map._domain or self._codomain != map._codomain:
            return False

        # Finally check the vector
        return map._vector in self._space

    def _identity(self):
        """
        Returns the identity map.

        OUTPUT:

        - QuiverRepHom

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a']}})
            sage: P = Q.P(QQ, 1)
            sage: H = P.Hom(P)
            sage: f = H.identity() # indirect doctest
            sage: f.is_isomorphism()
            True
        """

        from sage.matrix.constructor import Matrix
        maps = dict((v, Matrix(self._domain._spaces[v].dimension(),
                               self._domain._spaces[v].dimension(), self._base(1)))
                               for v in self._quiver)
        return QuiverRepHom(self._domain, self._codomain, maps)

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view and modify the representation data. #
    #                                                                         #
    ###########################################################################

    def base_ring(self):
        """
        Returns the base ring of the representations.

        OUTPUT:

        - ring, the base ring of the representations

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: H = QuiverHomSpace(Q.S(QQ, 2), Q.P(QQ, 1))
            sage: H.base_ring()
            Rational Field
        """

        return self._base

    def quiver(self):
        """
        Returns the quiver of the representations.

        OUTPUT:

        - Quiver, the quiver of the representations

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: H = QuiverHomSpace(Q.S(QQ, 2), Q.P(QQ, 1))
            sage: H.quiver() is Q
            True
        """

        return self._quiver

    def domain(self):
        """
        Returns the domain of the hom space.

        OUTPUT:

        - QuiverRep, the domain of the Hom space

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: S = Q.S(QQ, 2)
            sage: H = QuiverHomSpace(S, Q.P(QQ, 1))
            sage: H.domain() is S
            True
        """

        return self._domain

    def codomain(self):
        """
        Returns the codomain of the hom space.

        OUTPUT:

        - QuiverRep, the codomain of the Hom space

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: P = Q.P(QQ, 1)
            sage: H = QuiverHomSpace(Q.S(QQ, 2), P)
            sage: H.codomain() is P
            True
        """

        return self._codomain

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the representation.       #
    #                                                                         #
    ###########################################################################

    def dimension(self):
        """
        Returns the dimension of the hom space.

        OUTPUT:

        - integer, the dimension

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: H = QuiverHomSpace(Q.S(QQ, 2), Q.P(QQ, 1))
            sage: H.dimension()
            2
        """

        return self._space.dimension()

    def gens(self):
        """
        Returns a list of generators of the hom space

        OUTPUT:

        - list of QuiverRepHom objects, the generators

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: H = QuiverHomSpace(Q.S(QQ, 2), Q.P(QQ, 1))
            sage: H.gens()
            [Homomorphism of representations of Quiver on 2 vertices, Homomorphism of representations of Quiver on 2 vertices]
        """

        return [QuiverRepHom(self._domain, self._codomain, f) for f in self._space.gens()]

    def coordinates(self, hom):
        """
        Returns the coordinates of the map when expressed in terms of gens.

        INTPUT:

        - ``hom`` - QuiverRepHom

        OUTPUT:

        - list, the coordinates of the given map when written in terms of the
          generators of the QuiverHomSpace

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: S = Q.S(QQ, 2)
            sage: P = Q.P(QQ, 1)
            sage: H = QuiverHomSpace(S, P)
            sage: f = S.hom(P, {2: [[1,-1]]})
            sage: H.coordinates(f)
            [1, -1]
        """

        #Use the coordinates function on space
        return self._space.coordinates(hom._vector)

        ###########################################################################
        #                                                                         #
        # CONSTRUCTION FUNCTIONS                                                  #
        #    These functions create and return modules and homomorphisms.         #
        #                                                                         #
        ###########################################################################

    def hom(self, maps = {}):
        """
        Creates a homomorphism.

        See the QuiverRepHom documentation for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverHomSpace, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: H1 = QuiverHomSpace(S, M)
            sage: f = H1()
            sage: f.is_zero()
            True
            sage: maps2 = {2:[1, -1], 3:1}
            sage: g = H1(maps2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: h = H1([x, y])
            sage: g == h
            True
            sage: Proj = Q.P(GF(7), 3)
            sage: Simp = Q.S(GF(7), 3)
            sage: H2 = QuiverHomSpace(Proj, Simp)
            sage: im = H2({3: (1,)})
            sage: H2(im).is_surjective()
            True
        """

        return QuiverRepHom(self._domain, self._codomain, maps)

    def an_element(self):
        """
        Returns a homomorphism in the Hom space.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverHomSpace
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: S = Q.S(QQ, 2)
            sage: P = Q.P(QQ, 1)
            sage: H = QuiverHomSpace(S, P)
            sage: H.an_element() in H
            True
        """

        return QuiverRepHom(self._domain, self._codomain, self._space.an_element())

    def left_module(self, basis=False):
        """
        Creates the QuiverRep of self as a module over the opposite quiver.

        INPUT:

        - ``basis`` - bool, if false then only the module is returned.  If true then a
          tuple is returned.  The first element is the QuiverRep and the second element
          is a dictionary which associates to each vertex a list.  The elements of this
          list a the homomorphisms which correspond to the basis elements of that
          vertex in the module.

        OUTPUT:

        - QuiverRep or tuple

        .. WARNING::

            The codomain of the Hom space must be a left module.

        .. NOTES::

            The left action of a path e on a map f is given by (ef)(m) = ef(m).  This
            gives the Hom space its structure as a left module over the path algebra.
            This is then converted to a right module over the path algebra of the
            opposite quiver ``Q.reverse()`` and returned.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b'], 3: ['c', 'd']}, 2:{3:['e']}})
            sage: P = Q.P(GF(3), 3)
            sage: A = Q.free_module(GF(3))
            sage: H = P.Hom(A)
            sage: H.dimension()
            6
            sage: M, basis_dict = H.left_module(true)
            sage: M.dimension_vector()
            (4, 1, 1)
            sage: Q.reverse().P(GF(3), 3).dimension_vector()
            (4, 1, 1)

        As lists start indexing at 0 the ith vertex corresponds to the (i-1)th entry of
        the dimension vector::

            sage: len(basis_dict[2]) == M.dimension_vector()[1]
            True
        """

        if not self._codomain.is_left_module():
            raise ValueError("The codomain must be a left module.")

        # Create the spaces
        spaces = {}
        for v in self._quiver:
            im_gens = [self.hom([self._codomain.left_edge_action((v, v), f(x))
                for x in self._domain.gens()])._vector
                for f in self.gens()]
            spaces[v] = self._space.submodule(im_gens)

        # Create the maps
        maps = {}
        for e in self._quiver.edges():
            e_op = (e[1], e[0], e[2])
            maps[e_op] = []
            for vec in spaces[e[1]].gens():
                vec_im = spaces[e_op[1]].coordinate_vector(self.hom([self._codomain.left_edge_action(e, self.hom(vec)(x))
                                                                     for x in self._domain.gens()])._vector)
                maps[e_op].append(vec_im)

        # Create and return the module (and the dict if desired)
        if basis:
            basis_dict = {}
            for v in self._quiver:
                basis_dict[v] = [QuiverRepHom(self._domain, self._codomain, vec) for vec in spaces[v].gens()]
            return (QuiverRep(self._base, self._quiver.reverse(), spaces, maps), basis_dict)
        else:
            return QuiverRep(self._base, self._quiver.reverse(), spaces, maps)

class QuiverRepHom(CallMorphism):
    """
    A homomorphism of quiver representations is for each vertex of the quiver a
    homomorphism of the spaces assigned to those vertices such that these
    homomorphisms commute with the edge maps.  The domain and codomain of the
    homomorphism are required to be representations over the same quiver with
    the same base ring.

    INPUT:

    - ``domain`` - QuiverRep, the domain of the homomorphism

    - ``codomain`` - QuiverRep, the codomain of the homomorphism

    - ``data`` - dict, list, or QuiverRepElement (default: empty dict) as follows
      - list, data can be a list of images for the generators of the domain.  An
        error will be generated if the map so defined is not equivariant with
        respect to the action of the quiver.
      - dictionary, data can be a dictionary associating to each vertex of the
        quiver either a homomorphism with domain and codomain the spaces associated
        to this vertex in the domain and codomain modules respectively, or a matrix
        defining such a homomorphism, or an object that sage can construct such a
        matrix from.  Not all vertices must be specified, unspecified vertices are
        assigned the zero map, and keys not corresponding to vertices of the quiver
        are ignored.  An error will be generated if these maps do not commute with
        the edge maps of the domain and codomain.
      - QuiverRepElement, if the domain is a QuiverRep_with_path_basis then data
        can be a single QuiverRepElement belonging to the codomain.  The map is
        then defined by sending each path, p, in the basis to data*p.  If data is
        not an element of the codomain or the domain is not a
        QuiverRep_with_path_basis then an error will be generated.
      - QuiverRepHom, the input can also be a map ``f:D -> C`` such that there is a
        coercion from the domain of self to ``D`` and from ``C`` to the codomain of
        self.  The composition of these maps is the result.

    OUTPUT:

    - QuiverRepHom

    EXAMPLES::

        sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
        sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
        sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
        sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
        sage: M = QuiverRep(QQ, Q, spaces, maps)
        sage: spaces2 = {2: QQ^1, 3: QQ^1}
        sage: S = QuiverRep(QQ, Q, spaces2)

    With no additional data this creates the zero map::

        sage: f = QuiverRepHom(S, M)
        sage: f.is_zero()
        True

    We must specify maps at the vertices to get a nonzero homomorphism.  Note that
    if the dimensions of the spaces assigned to the domain and codomain of a vertex
    are equal then Sage will construct the identity matrix from ``1``::

        sage: maps2 = {2:[1, -1], 3:1}
        sage: g = QuiverRepHom(S, M, maps2)

    Here we create the same map by specifying images for the generators::

        sage: x = QuiverRepElement(M, {2: (1, -1)})
        sage: y = QuiverRepElement(M, {3: (1,)})
        sage: h = QuiverRepHom(S, M, [x, y])
        sage: g == h
        True

    If the domain is a module of type QuiverRep_with_path_basis (for example, the
    indecomposable projectives) we can create maps by specifying a single image::

        sage: Proj = Q.P(GF(7), 3)
        sage: Simp = Q.S(GF(7), 3)
        sage: im = QuiverRepElement(Simp, {3: (1,)})
        sage: QuiverRepHom(Proj, Simp, im).is_surjective()
        True
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, domain, codomain, data={}):
        """
        Type QuiverRepHom? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: f = QuiverRepHom(S, M)
            sage: f.is_zero()
            True
            sage: maps2 = {2:[1, -1], 3:1}
            sage: g = QuiverRepHom(S, M, maps2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: h = QuiverRepHom(S, M, [x, y])
            sage: g == h
            True
            sage: Proj = Q.P(GF(7), 3)
            sage: Simp = Q.S(GF(7), 3)
            sage: im = QuiverRepElement(Simp, {3: (1,)})
            sage: QuiverRepHom(Proj, Simp, im).is_surjective()
            True
        """
        # The data of a representation is held in the following private
        # variables:
        #
        # * _quiver
        #      The quiver of the representation.
        # * _base_ring
        #      The base ring of the representation.
        # * _domain
        #      The QuiverRep object that is the domain of the homomorphism.
        # * _codomain
        #      The QuiverRep object that is the codomain of the homomorphism.
        # * _vector
        #      A vector in some free module over the base ring of a length such
        #      that each coordinate corresponds to an entry in the matrix of a
        #      homomorphism attached to a vertex.
        #
        # The variable data can also be a vector of appropriate length.  When
        # this is the case it will be loaded directly into _vector and then
        # _assert_valid_hom is called.

        self._domain = domain
        self._codomain = codomain
        self._quiver = domain._quiver
        self._base_ring = domain._base_ring

        # Check that the quiver and base ring match
        if codomain._quiver != self._quiver:
            raise ValueError("The quivers of the domain and codomain must be equal.")
        if codomain._base_ring != self._base_ring:
            raise ValueError("The base ring of the domain and codomain must be equal.")

        # Get the dimensions of the spaces
        mat_dims = {}
        domain_dims = {}
        codomain_dims = {}
        for v in self._quiver:
            domain_dims[v] = domain._spaces[v].dimension()
            codomain_dims[v] = codomain._spaces[v].dimension()
            mat_dims[v] = domain_dims[v]*codomain_dims[v]
        total_dim = sum(mat_dims.values())

        # Handle the case when data is a vector
        if data in self._base_ring**total_dim:
            self._vector = data
            self._assert_valid_hom()
            return

        # If data is not a dict, create one
        if isinstance(data, dict):
            maps_dict = data
        else:
            # If data is not a list create one, then create a dict from it
            if isinstance(data, list):
                im_list = data
            else:
                # If data is a QuiverRepHom, create a list from it
                if isinstance(data, QuiverRepHom):
                    f = data._domain.coerce_map_from(domain)
                    g = self._codomain.coerce_map_from(data._codomain)
                    im_list = [g(data(f(x))) for x in domain.gens()]

                # The only case left is that data is a QuiverRepElement
                else:
                    if not isinstance(data, QuiverRepElement):
                        raise TypeError("Input data must be dictionary, list, " +
                                        "QuiverRepElement or vector.")
                    if not isinstance(domain, QuiverRep_with_path_basis):
                        raise TypeError("If data is a QuiverRepElement then domain " +
                                        "must be a QuiverRep_with_path_basis.")
                    if data not in codomain:
                        raise ValueError("If data is a QuiverRepElement then it must " +
                                         "be an element of codomain.")
                    im_list = [codomain.right_edge_action(data, p) for v in domain._quiver for p in domain._bases[v]]

            # WARNING: This code assumes that the function QuiverRep.gens() returns
            # the generators ordered first by vertex and then by the order of the
            # gens() method of the space associated to that vertex.  In particular
            # this is the order that corresponds to how maps are represented via
            # matrices

            # Get the gens of the domain and check that im_list is the right length
            dom_gens = domain.gens()
            if len(im_list) != len(dom_gens):
                raise ValueError("Domain is dimension " + str(len(dom_gens)) + " but only " + str(len(im_list)) +
                                 " images were supplied.")

            # Get the matrices of the maps
            start_index = 0
            maps_dict = {}
            for v in self._quiver:
                maps_dict[v] = []
                dim = domain._spaces[v].dimension()
                for i in range(start_index, start_index + dim):
                    if len(im_list[i].support()) != 0 and im_list[i].support() != [v]:
                        # If the element doesn't have the correct support raise
                        # an error here, otherwise we might create a valid hom
                        # that does not map the generators to the supplied
                        # images
                        raise ValueError("Generator supported at vertex " + str(v) +
                                         " cannot map to element with support " + str(im_list[i].support()))
                    else:
                        # If the support works out add the images coordinates
                        # as a row of the matrix
                        maps_dict[v].append(codomain._spaces[v].coordinates(im_list[i]._elems[v]))

                start_index += dim

        # Get the coordinates of the vector
        from sage.categories.morphism import is_Morphism
        from sage.matrix.constructor import Matrix
        vector = []
        for v in self._quiver:
            if v in maps_dict:
                if is_Morphism(maps_dict[v]):
                    if hasattr(maps_dict[v], 'matrix'):
                        m = maps_dict[v].matrix()
                    else:
                        gens_images = [codomain._spaces[v].coordinate_vector(maps_dict[v](x))
                                       for x in domain._spaces[v].gens()]
                        m = Matrix(self._base_ring, domain_dims[v], codomain_dims[v], gens_images)
                else:
                    m = Matrix(self._base_ring, domain_dims[v], codomain_dims[v], maps_dict[v])
            else:
                m = Matrix(self._base_ring, domain_dims[v], codomain_dims[v])
            for i in range(0, domain_dims[v]):
                vector += list(m[i])

        # Wrap as a vector, check it, and return
        self._vector = (self._base_ring**total_dim)(vector)
        self._assert_valid_hom()

        super(QuiverRepHom, self).__init__(QuiverHomSpace(domain, codomain), codomain)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: QuiverRepHom(S, M) # indirect doctest
            Homomorphism of representations of Quiver on 3 vertices
        """

        return "Homomorphism of representations of " + self._quiver.__repr__()

    def __call__(self, x):
        """
        This function overloads functional notation f(x).

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: h = QuiverRepHom(S, M, [x, y])
            sage: h(S.gens()[0]) == x
            True
            sage: h(S.gens()[1]) == y
            True
        """

        # Check the input
        if not isinstance(x, QuiverRepElement):
            raise ValueError("QuiverRepHom can only be called on QuiverRepElement")

        elements = dict((v, self.get_map(v)(x._elems[v])) for v in self._quiver)
        return QuiverRepElement(self._codomain, elements)

    def __add__(left, right):
        """
        This function overloads the + operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: z = M.zero()
            sage: h = QuiverRepHom(S, M, [x, z])
            sage: g = QuiverRepHom(S, M, [z, z])
            sage: f = g + h
            sage: f(S.gens()[0]) == x
            True
            sage: f(S.gens()[1]) == z
            True
        """

        new_vector = left._vector + right._vector
        return QuiverRepHom(left._domain, left._codomain, new_vector)

    def __iadd__(self, other):
        """
        This function overloads the += operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: z = M.zero()
            sage: h = QuiverRepHom(S, M, [x, z])
            sage: g = QuiverRepHom(S, M, [z, z])
            sage: g += h
            sage: g(S.gens()[0]) == x
            True
            sage: g(S.gens()[1]) == z
            True
        """

        self._vector += other._vector

        return self

    def __sub__(left, right):
        """
        This function overloads the - operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: z = M.zero()
            sage: h = QuiverRepHom(S, M, [x, z])
            sage: g = QuiverRepHom(S, M, [z, y])
            sage: f = h - g
            sage: f(S.gens()[0]) == x
            True
            sage: f(S.gens()[1]) == -y
            True
        """

        new_vector = left._vector - right._vector
        return QuiverRepHom(left._domain, left._codomain, new_vector)

    def __isub__(self, other):
        """
        This function overloads the -= operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: z = M.zero()
            sage: h = QuiverRepHom(S, M, [x, z])
            sage: g = QuiverRepHom(S, M, [z, y])
            sage: h -= g
            sage: h(S.gens()[0]) == x
            True
            sage: h(S.gens()[1]) == -y
            True
        """

        self._vector -= other._vector

        return self

    def __neg__(self):
        """
        This function overrides the unary - operator

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: h = QuiverRepHom(S, M, [x, y])
            sage: g = -h
            sage: g(S.gens()[0]) == -x
            True
            sage: g(S.gens()[1]) == -y
            True
        """

        return QuiverRepHom(self._domain, self._codomain, -self._vector)

    def __pos__(self):
        """
        This function overrides the unary - operator

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: h = QuiverRepHom(S, M, [x, y])
            sage: g = +h
            sage: g == h
            True
        """

        return self

    def __eq__(self, other):
        """
        This function overrides the == operator

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: g = QuiverRepHom(S, M, [x, y])
            sage: h = QuiverRepHom(S, M, [x, y])
            sage: g == h
            True
        """

        # A homomorphism can only be equal to another homomorphism between the
        # same domain and codomain
        if not isinstance(other, QuiverRepHom) or self._domain != other._domain or self._codomain != other._codomain:
            return False

        # If all that holds just check the vectors
        return self._vector == other._vector

    def __ne__(self, other):
        """
        This function overrides the != operator

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = QuiverRepElement(M, {2: (1, -1)})
            sage: y = QuiverRepElement(M, {3: (1,)})
            sage: z = M.zero()
            sage: g = QuiverRepHom(S, M, [x, y])
            sage: h = QuiverRepHom(S, M, [x, z])
            sage: g != h
            True
        """

        # A homomorphism can only be equal to another homomorphism between the
        # same domain and codomain
        if not isinstance(other, QuiverRepHom) or self._domain != other._domain or self._codomain != other._codomain:
            return True

        # If all that holds just check the vectors
        return self._vector != other._vector

    def __mul__(self, other):
        """
        This function overrides the * operator

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom, QuiverRepElement
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: x = S.gens()[0]
            sage: y = S.gens()[1]
            sage: g = QuiverRepHom(S, S, [x, y])
            sage: h = QuiverRepHom(S, S)
            sage: (g*h).is_zero()
            True
        """

        maps = dict((v, other.get_matrix(v)*self.get_matrix(v)) for v in self._quiver)
        return QuiverRepHom(other._domain, self._codomain, maps)

    ###########################################################################
    #                                                                         #
    # WELL DEFINEDNESS FUNCTIONS                                              #
    #    These functions test and assert well definedness of the              #
    #    homomorphism.                                                        #
    #                                                                         #
    ###########################################################################

    def _assert_valid_hom(self):
        """
        Raises a ValueError if the homomorphism is not well defined.

        Specifically it checks that the domain and codomains of the maps are correct
        and that the edge diagrams commute.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^1, 3: QQ^1}
            sage: S = QuiverRep(QQ, Q, spaces2)
            sage: maps2 = {2:[1, -1], 3:1}
            sage: g = QuiverRepHom(S, M, maps2) # indirect doctest
            sage: f = QuiverRepHom(S, S, maps2) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: entries has the wrong length
        """

        # Check that the domain and codomains dimensions add correctly
        totaldim = 0
        for v in self._quiver:
            totaldim += self._domain._spaces[v].dimension()*self._codomain._spaces[v].dimension()
        if totaldim != len(self._vector):
            raise ValueError("Dimensions do not match domain and codomain.")

        # Check that the edge diagrams commute
        for e in self._quiver.edges():
            if self.get_matrix(e[0])*self._codomain._maps[e].matrix() != self._domain._maps[e].matrix()*self.get_matrix(e[1]):
                raise ValueError("The diagram of edge " + str(e) + " does not commute.")

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view the homomorphism data.              #
    #                                                                         #
    ###########################################################################

    def domain(self):
        """
        Returns the domain of the homomorphism.

        OUTPUT:

        - QuiverRep, the domain

        sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom
        sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
        sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
        sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
        sage: M = QuiverRep(QQ, Q, spaces, maps)
        sage: S = QuiverRep(QQ, Q)
        sage: g = QuiverRepHom(M, S)
        sage: g.domain() is M
        True
        """

        return self._domain

    def codomain(self):
        """
        Returns the codomain of the homomorphism.

        OUTPUT:

        - QuiverRep, the codomain

        sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom
        sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
        sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
        sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
        sage: M = QuiverRep(QQ, Q, spaces, maps)
        sage: S = QuiverRep(QQ, Q)
        sage: g = QuiverRepHom(S, M)
        sage: g.codomain() is M
        True
        """

        return self._codomain

    def get_matrix(self, vertex):
        """
        Returns the matrix of the homomorphism attached to vertex.

        INPUT:

        - ``vertex`` - integer, a vertex of the quiver

        OUTPUT:

        - matrix, the matrix representing the homomorphism associated to the given
          vertex

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: I = Q.I(QQ, 3)
            sage: M = I/I.radical()
            sage: f = M.coerce_map_from(I)
            sage: f.get_matrix(1)
            [1 0]
            [0 1]
        """

        # Get dimensions
        startdim = 0
        for v in self._quiver:
            if v == vertex:
                break
            startdim += self._domain._spaces[v].dimension()*self._codomain._spaces[v].dimension()
        rows = self._domain._spaces[vertex].dimension()
        cols = self._codomain._spaces[vertex].dimension()

        # Slice out the matrix and return
        from sage.matrix.constructor import Matrix
        mat = list(self._vector[startdim:startdim + rows*cols])
        return Matrix(self._base_ring, rows, cols, mat)

    def get_map(self, vertex):
        """
        Returns the homomorphism at the given vertex.

        INPUT:

        - ``vertex`` - integer, a vertex of the quiver

        OUTPUT:

        - homomorphism, the homomorphism associated to the given vertex

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: S = P/P.radical()
            sage: f = S.coerce_map_from(P)
            sage: f.get_map(1).is_bijective()
            True
        """

        return self._domain._spaces[vertex].hom(self.get_matrix(vertex), self._codomain._spaces[vertex])

    def quiver(self):
        """
        Return the quiver of the representations in the domain/codomain.

        OUTPUT:

        - Quiver, the quiver of the representations in the domain and codomain

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: f = QuiverRepHom(P, P, {1: 1, 2: 1, 3: 1})
            sage: f.quiver() is Q
            True
        """

        return self._quiver

    def base_ring(self):
        """
        Return the base ring of the representation in the codomain.

        OUTPUT:

        - ring, the base ring of the codomain

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: f = QuiverRepHom(P, P, {1: 1, 2: 1, 3: 1})
            sage: f.base_ring() is QQ
            True
        """

        return self._base_ring

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the homomorphism.         #
    #                                                                         #
    ###########################################################################

    def is_injective(self):
        """
        Tests whether the homomorphism is injective.

        OUTPUT:

        - bool, True if the homomorphism is injective, False otherwise

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: f = QuiverRepHom(P, P, {1: 1, 2: 1, 3: 1})
            sage: f.is_injective()
            True
            sage: g = QuiverRepHom(P, P)
            sage: g.is_injective()
            False
        """

        # The homomorphism is injective if and only if it is injective at every
        # vertex
        for v in self._quiver:
            if self.get_matrix(v).nullity() != 0:
                return False

        return True

    def is_surjective(self):
        """
        Tests whether the homomorphism is surjective.

        OUTPUT:

        - bool, True if the homomorphism is surjective, False otherwise

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: f = QuiverRepHom(P, P, {1: 1, 2: 1, 3: 1})
            sage: f.is_surjective()
            True
            sage: g = QuiverRepHom(P, P)
            sage: g.is_surjective()
            False
        """

        # The homomorphism is surjective if and only if it is surjective at
        # every vertex
        for v in self._quiver:
            m = self.get_matrix(v)
            if m.rank() != m.ncols():
                return False

        return True

    def is_isomorphism(self):
        """
        Tests whether the homomorphism is an isomorphism.

        OUTPUT:

        - bool, True if the homomorphism is bijective, False otherwise

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: f = QuiverRepHom(P, P, {1: 1, 2: 1, 3: 1})
            sage: f.is_isomorphism()
            True
            sage: g = QuiverRepHom(P, P)
            sage: g.is_isomorphism()
            False
        """

        # It's an iso if and only if it's an iso at every vertex
        for v in self._quiver:
            if not self.get_matrix(v).is_invertible():
                return False

        return True

    def is_zero(self):
        """
        Tests whether the homomorphism is the zero homomorphism.

        OUTPUT:

        - bool, True if the homomorphism is zero, False otherwise

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: f = QuiverRepHom(P, P, {1: 1, 2: 1, 3: 1})
            sage: f.is_zero()
            False
            sage: g = QuiverRepHom(P, P)
            sage: g.is_zero()
            True
        """

        # The homomorphism is zero if and only if it is zero at every vertex
        for v in self._quiver:
            if not self.get_matrix(v).is_zero():
                return False

        return True

    def is_endomorphism(self):
        """
        Tests whether the homomorphism is an endomorphism.

        OUTPUT:

        - bool, True if the domain equals the codomain, False otherwise

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: f = QuiverRepHom(P, P, {1: 1, 2: 1, 3: 1})
            sage: f.is_endomorphism()
            True
            sage: S = P/P.radical()
            sage: g = S.coerce_map_from(P)
            sage: g.is_endomorphism()
            False
        """

        return self._domain == self._codomain

    def rank(self):
        """
        Returns the rank.

        OUTPUT:

        - integer, the rank

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: S = P/P.radical()
            sage: f = S.coerce_map_from(P)
            sage: assert(f.rank() == 1)
        """

        # The rank is the sum of the ranks at each vertex
        r = 0
        for v in self._quiver:
            r += self.get_matrix(v).rank()

        return r

    ###########################################################################
    #                                                                         #
    # CONSTRUCTION FUNCTIONS                                                  #
    #    These functions create new homomorphisms, representations, and       #
    #    elements from the given homomorphism.                                #
    #                                                                         #
    ###########################################################################

    def kernel(self):
        """
        Returns the kernel of self.

        OUTPUT:

        - QuiverRep, the kernel

        .. NOTES::

            To get the inclusion map of the kernel, ``K``, into the domain, ``D``, use
            ``D.coerce_map_from(K)``.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^2, 3: QQ^1}
            sage: N = QuiverRep(QQ, Q, spaces2, {(2, 3, 'c'): [[1], [0]]})
            sage: maps2 = {2:[[1, 0], [0, 0]], 3:1}
            sage: g = QuiverRepHom(N, M, maps2)
            sage: g.kernel().dimension_vector()
            (0, 1, 0)
        """

        spaces = dict((v, self.get_map(v).kernel()) for v in self._quiver)
        return self._domain._submodule(spaces)

    def image(self):
        """
        Returns the image of self.

        OUTPUT:

        - QuiverRep, the image

        .. NOTES::

            To get the inclusion map of the image, ``I``, into the codomain, ``C``, use
            ``C.coerce_map_from(I)``.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^2, 3: QQ^1}
            sage: N = QuiverRep(QQ, Q, spaces2, {(2, 3, 'c'): [[1], [0]]})
            sage: maps2 = {2:[[1, 0], [0, 0]], 3:1}
            sage: g = QuiverRepHom(N, M, maps2)
            sage: g.image().dimension_vector()
            (0, 1, 1)
        """

        spaces = dict((v, self.get_map(v).image()) for v in self._quiver)
        return self._codomain._submodule(spaces)

    def cokernel(self):
        """
        Returns the cokernel of self.

        OUTPUT:

        - QuiverRep, the cokernel

        .. NOTES::

            To get the factor map of the codomain, ``D``, onto the cokernel, ``C``, use
            ``C.coerce_map_from(D)``.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepHom
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: spaces = {1: QQ^2, 2: QQ^2, 3:QQ^1}
            sage: maps = {(1, 2, 'a'): [[1, 0], [0, 0]], (1, 2, 'b'): [[0, 0], [0, 1]], (2, 3, 'c'): [[1], [1]]}
            sage: M = QuiverRep(QQ, Q, spaces, maps)
            sage: spaces2 = {2: QQ^2, 3: QQ^1}
            sage: N = QuiverRep(QQ, Q, spaces2, {(2, 3, 'c'): [[1], [0]]})
            sage: maps2 = {2:[[1, 0], [0, 0]], 3:1}
            sage: g = QuiverRepHom(N, M, maps2)
            sage: g.cokernel().dimension_vector()
            (2, 1, 0)
        """

        return self._codomain.quotient(self.image())

    def linear_dual(self):
        """
        Computes the linear dual Df:DN->DM of self = f:M->N where D(-) = Hom_k(-, k).

        OUTPUT:

        - QuiverRepHom, the map Df:DN->DM

        .. NOTES::

            If e is an edge of the quiver Q and g is an element of Hom_k(N, k) then we
            let (ga)(m) = g(ma).  This gives Hom_k(N, k) its structure as a module over
            the opposite quiver Q.reverse().  The map Hom_k(N, k) -> Hom_k(M, k)
            returned sends g to gf.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: S = P/P.radical()
            sage: f = S.coerce_map_from(P)

        The dual of a surjective map is injective and vice versa::

            sage: f.is_surjective()
            True
            sage: g = f.linear_dual()
            sage: g.is_injective()
            True

        The dual of a right module is a left module for the same quiver, Sage
        represents this as a right module for the opposite quiver::

            sage: g.quiver() == Q.reverse()
            True

        The double dual of a map is the original representation::

            sage: g.linear_dual() == f
            True
        """

        # The effect of the functor D is that it just transposes the matrix of
        # a hom
        maps = dict((v, self.get_matrix(v).transpose()) for v in self._quiver)
        return QuiverRepHom(self._codomain.linear_dual(), self._domain.linear_dual(), maps)

    def algebraic_dual(self):
        """
        Computes the algebraic dual f^t:N^t->M^t of self = f:M->N where (-)^t = Hom_Q(-, kQ).

        OUTPUT:

        - QuiverRepHom, the map f^t:N^t->M^t

        .. NOTES::

            If e is an edge of the quiver Q and g is an element of Hom_Q(N, kQ) then we
            let (ge)(m) = eg(m).  This gives Hom_Q(N, kQ) its structure as a module over
            the opposite quiver Q.reverse().  The map Hom_Q(N, kQ) -> Hom_Q(M, kQ)
            returned sends g to gf.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a'], 3:['b','c','d']}, 2:{4:['e','f']}, 3:{4:['g']}, 5:{2:['h','i']}})
            sage: P1 = Q.P(QQ, 4)
            sage: P1.algebraic_dual()
            Representation with dimension vector (5, 2, 1, 1, 4)

        The algebraic dual of an indecomposable projective is the indecomposable
        projective of the same vertex in the opposite quiver.

            sage: Q.reverse().P(QQ, 4)
            Representation with dimension vector (5, 2, 1, 1, 4)
        """

        # Get the domain, its basis, and the codomain
        domain, domain_gens = self._codomain.algebraic_dual(True)
        codomain, co_domain_gens = self._domain.algebraic_dual(True)

        # Find the images in the domain and create the module
        # H = QuiverHomSpace(self._domain, self._quiver.free_module(self._base_ring))
        im_gens = [QuiverRepElement(codomain, {v: (g*self)._vector})
                    for v in self._quiver for g in domain_gens[v]]
        return QuiverRepHom(domain, codomain, im_gens)

    def direct_sum(self, maps, return_maps=False, pinch=None):
        """
        Returns the direct sum of self with the maps in the list ``maps``.

        INPUT:

        - ``maps`` - QuiverRepHom or list of QuiverRepHoms

        - ``return_maps`` - bool (default: False), if False then the return value is a
          QuiverRepHom which is the direct sum of self with the QuiverRepHoms in ``maps``.
          If True then the return value is a tuple of length either 3 or 5.  The first
          entry of the tuple is the QuiverRepHom giving the direct sum.  If ``pinch`` is
          either None or 'codomain' then the next two entries in the tuple are lists
          giving respectively the inclusion and the projection maps for the factors of
          the direct sum.  Summands are ordered as given in maps with self as the
          zeroth summand.  If ``pinch`` is either None or 'domain' then the next two
          entries in the tuple are the inclusion and projection maps for the codomain.
          Thus if ``pinch`` is None then the tuple will have length 5.  If ``pinch`` is either
          'domain' or 'codomain' then the tuple will have length 3.

        - ``pinch`` - string or None (default: None), if equal to 'domain' then the domains
          of self and the given maps must be equal.  The direct sum of f: A -> B and
          g: A -> C returned is the map A -> B (+) C defined by sending x to
          (f(x), g(x)).  If ``pinch`` equals 'codomain' then the codomains of self and the
          given maps must be equal.  The direct sum of f: A -> C and g: B -> C returned
          is the map A (+) B -> C defined by sending (x, y) to f(x) + g(y).  Finally if
          ``pinch`` is anything other than 'domain' or 'codomain' then the direct sum of
          f: A -> B and g: C -> D returned is the map A (+) C -> B (+) D defined by
          sending (x, y) to (f(x), f(y)).

        OUTPUT:

            - QuiverRepHom or tuple

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a', 'b']}})
            sage: P1 = Q.P(GF(3), 1)
            sage: P2 = Q.P(GF(3), 2)
            sage: S1 = P1/P1.radical()
            sage: S2 = P2/P2.radical()
            sage: pi1 = S1.coerce_map_from(P1)
            sage: pi2 = S2.coerce_map_from(P2)
            sage: f = pi1.direct_sum(pi2)
            sage: f.domain().dimension_vector() == Q.free_module(GF(3)).dimension_vector()
            True
            sage: f.is_surjective()
            True
            sage: id = P1.Hom(P1).identity()
            sage: g = pi1.direct_sum(id, pinch='domain')
            sage: g.is_surjective()
            False
        """

        # Get the list of maps to be summed
        if isinstance(maps, QuiverRepHom):
            maplist = [self, maps]
        else:
            maplist = [self] + maps

        # Check that the quivers/base rings are the same.  If pinching also
        # check that the domain/codomains are correct
        for x in maplist:
            if not isinstance(x, QuiverRepHom):
                raise TypeError("maps must be a QuiverRepHom or list of QuiverRepHoms")
            if self._quiver is not x._quiver:
                raise ValueError("Cannot direct sum maps from different quivers")
            if self._base_ring is not x._base_ring:
                raise ValueError("Base rings must be identical")
            if pinch == 'domain' and self._domain is not x._domain:
                raise ValueError("Cannot pinch maps, domains do not agree")
            if pinch == 'codomain' and self._codomain is not x._codomain:
                raise ValueError("Cannot pinch maps, codomains do not agree")

        # Get the sums and their maps
        if pinch == 'domain':
            domain = self._domain
        else:
            domain, d_incl, d_proj = self._domain.direct_sum([x._domain for x in maplist[1:]], return_maps=True)
        if pinch == 'codomain':
            codomain = self._codomain
        else:
            codomain, c_incl, c_proj = self._codomain.direct_sum([x._codomain for x in maplist[1:]], return_maps=True)

        # Start with the zero map
        result = QuiverRepHom(domain, codomain)

        # Add each factor
        for i in range(0, len(maplist)):
            if pinch == 'domain':
                result += c_incl[i]*maplist[i]
            elif pinch == 'codomain':
                result += maplist[i]*d_proj[i]
            else:
                result += c_incl[i]*maplist[i]*d_proj[i]

        # Return the results
        if return_maps:
            if pinch == 'domain':
                return (result, c_incl, c_proj)
            elif pinch == 'codomain':
                return (result, d_incl, d_proj)
            else:
                return (result, d_incl, d_proj, c_incl, c_proj)
        else:
            return result

    def lift(self, x):
        """
        Given an element of the image, return an element of the codomain that maps onto
        it.

        INPUT:

        - ``x`` - QuiverRepElement

        OUTPUT:

        - QuiverRepElement

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['c','d']}})
            sage: P = Q.P(RR, 3)
            sage: S = P/P.radical()
            sage: proj = S.coerce_map_from(P)
            sage: x = S.an_element()
            sage: y = proj.lift(x)
            sage: proj(y) == x
            True
            sage: zero = S.hom(S, {})
            sage: zero.lift(x)
            Traceback (most recent call last):
            ...
            ValueError: element is not in the image
        """

        # Lift at each vertex
        elems = dict((v, self.get_map(v).lift(x._elems[v])) for v in self._quiver)
        return QuiverRepElement(self._domain, elems)

    ###########################################################################
    #                                                                         #
    # ADDITIONAL OPERATIONS                                                   #
    #    These functions operations that are not implemented via binary       #
    #    operators.                                                           #
    #                                                                         #
    ###########################################################################

    def scalar_mult(self, scalar):
        """
        Returns the result of the scalar multiplcation scalar*self.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}})
            sage: M = Q.P(QQ, 1)
            sage: f = M.Hom(M).an_element()
            sage: x = M.an_element()
            sage: g = f.scalar_mult(6)
            sage: g(x) == 6*f(x)
            True
        """

        return QuiverRepHom(self._domain, self._codomain, scalar*self._vector)

    def iscalar_mult(self, scalar):
        """
        Multiplies self by scalar in place.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a','b']}})
            sage: M = Q.P(QQ, 1)
            sage: f = M.Hom(M).an_element()
            sage: x = M.an_element()
            sage: y = f(x)
            sage: f.iscalar_mult(6)
            sage: f(x) == 6*y
            True
        """

        self._vector *= scalar

class QuiverRepElement(ModuleElement):
    """
    An element of a quiver representation is a choice of element from each of the
    spaces assigned to the vertices of the quiver.  Addition, subtraction, and
    scalar multiplication of these elements is done pointwise within these spaces.

    INPUT:

    - ``module`` - QuiverRep (default: None), the module to which the element belongs

    - ``elements`` - dict (default: empty), a dictionary associating to each vertex a
      vector or an object from which sage can create a vector.  Not all vertices must
      be specified, unspecified vertices will be assigned the zero vector of the space
      associated to that vertex in the given module.  Keys that do not correspond to a
      vertex are ignored.

    - ``name`` - string (default: None), the name of the element

    OUTPUT:

    - QuiverRepElement

    .. NOTES::

        The constructor needs to know the quiver in order to create an element of a
        representation over that quiver.  The default is to read this information
        from ``module`` as well as to fill in unspecified vectors with the zeros of the
        spaces in ``module``.  If ``module``=None then ``quiver`` MUST be a quiver and each
        vertex MUST be specified or an error will result.  If both ``module`` and
        ``quiver`` are given then ``quiver`` is ignored.

    EXAMPLES::

        sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
        sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
        sage: spaces = dict((v, GF(3)^2) for v in Q)
        sage: M = QuiverRep(GF(3), Q, spaces)
        sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
        sage: QuiverRepElement(M, elems)
        Element of quiver representation
        sage: v = QuiverRepElement(M, elems, 'v')
        sage: v
        v
        sage: (v + v + v).is_zero()
        True
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, module, elements={}, name=None):
        """
        Type QuiverRepElement? for more information.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: QuiverRepElement(M)
            Element of quiver representation
            sage: v = QuiverRepElement(M, elems, 'v')
            sage: v
            v
            sage: (v + v + v).is_zero()
            True
        """
        # The data describing an element is held in the following private
        # variables:
        #
        # * _elems
        #      A dictionary that assigns to each vertex of the quiver a choice
        #      of element from the space assigned to that vertex in the parent
        #      representation.
        # * _quiver
        #      The quiver of the representation.

        super(QuiverRepElement, self).__init__(module)

        self._elems = {}
        self._quiver = module._quiver
        for v in self._quiver:
            if v in elements:
                self._elems[v] = module._spaces[v](elements[v])
            else:
                self._elems[v] = module._spaces[v].zero()

        # Assign a name if supplied
        if name is not None:
            self.rename(name)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: Q.P(QQ, 3).an_element() # indirect doctest
            Element of quiver representation
        """
        return "Element of quiver representation"

    def _add_(left, right):
        """
        This overrides the + operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: (v + v + v).is_zero() # indirect doctest
            True
        """

        elements = {}
        for v in left._quiver:
            elements[v] = left._elems[v] + right._elems[v]

        return QuiverRepElement(left.parent(), elements)

    def _iadd_(self, right):
        """
        This overrides the += operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: w = QuiverRepElement(M, elems)
            sage: v += w # indirect doctest
            sage: v == w
            False
            sage: v.is_zero()
            False
        """

        for v in self._quiver:
            self._elems[v] += right._elems[v]

        return self

    def _sub_(left, right):
        """
        This overrides the - operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: (v - v).is_zero() # indirect doctest
            True
        """

        elements = {}
        for v in left._quiver:
            elements[v] = left._elems[v] - right._elems[v]

        return QuiverRepElement(left.parent(), elements)

    def _isub_(self, right):
        """
        This overrides the -= operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: v -= v # indirect doctest
            sage: v.is_zero()
            True
        """

        for v in self._quiver:
            self._elems[v] -= right._elems[v]

        return self

    def _neg_(self):
        """
        This overrides the unary - operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: v == -v # indirect doctest
            False
        """

        elements = {}
        for v in self._quiver:
            elements[v] = -self._elems[v]

        return QuiverRepElement(self.parent(), elements)

    def __mul__(self, other):
        """
        Implements * for right multiplication by Quiver Algebra elements.

        TESTS::

            sage: from sage.modules.quiver_module import *
            sage: Q = Quiver({1:{2:['a']}})
            sage: P = Q.P(QQ, 1)
            sage: A = Q.algebra(QQ)
            sage: m = P.an_element()
            sage: a = A((1, 2, 'a'))
            sage: e1 = A((1, 1))
            sage: m.support()
            [1, 2]
            sage: (m*a).support()
            [2]
            sage: (m*e1).support()
            [1]
            sage: (m*(a + e1)).support()
            [1, 2]
        """

        # Make sure the input is an element of the quiver algebra and get the
        # coefficients of the monomials in it
        parent = self.parent()
        mons = parent._base(other).monomial_coefficients()
        result = parent.zero()

        for path in mons:
            # Multiply by the scalar
            x = mons[path]*self._elems[path.initial_vertex()]

            # If the edge isn't trivial apply the corresponding maps
            if len(path) > 0:
                for edge in path:
                    x = parent.get_map(edge)(x)

            # Sum the results of acting by each monomial
            result._elems[path.terminal_vertex()] += x

        return result

    def __eq__(self, other):
        """
        This overrides the == operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: w = QuiverRepElement(M, elems)
            sage: v == w
            True
            sage: v += w
            sage: v == w
            False
        """

        # Return False if being compared to something other than a
        # QuiverRepElement or if comparing two elements from representations
        # with different quivers
        if not isinstance(other, QuiverRepElement) or self._quiver != other._quiver:
            return False

        # Return False if the elements differ at any vertex
        for v in self._quiver:
            if self._elems[v] != other._elems[v]:
                return False

        return True

    def __ne__(self, other):
        """
        This overrides the != operator.

        TESTS::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: w = QuiverRepElement(M, elems)
            sage: v != w
            False
            sage: v += w
            sage: v != w
            True
        """

        # Return True if being compared to something other than a
        # QuiverRepElement or if comparing two elements from representations
        # with different quivers
        if not isinstance(other, QuiverRepElement) or self._quiver != other._quiver:
            return True

        # Return True if the elements differ at any vertex
        for v in self._quiver:
            if self._elems[v] != other._elems[v]:
                return True

        return False

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view and modify the representation data. #
    #                                                                         #
    ###########################################################################

    def quiver(self):
        """
        Returns the quiver of the representation.

        OUTPUT:

        - Quiver, the quiver of the representation

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: P = Q.P(QQ, 1)
            sage: v = P.an_element()
            sage: v.quiver() is Q
            True
        """

        return self._quiver

    def get_element(self, vertex):
        """
        Returns the element at the given vertex.

        INPUT:

        - ``vertex`` - integer, a vertex of the quiver

        OUTPUT:

        - vector, the vaector assigned to the given vertex

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: v.get_element(1)
            (1, 0)
            sage: v.get_element(3)
            (2, 1)
        """

        return self._elems[vertex]

    def set_element(self, vector, vertex):
        """
        Sets the element at the given vertex.

        INPUT:

        - ``vector`` - a vector or an object from which the space associated to the given
          vertex in the parent can create a vector

        - ``vertex`` - integer, a vertex of the quiver

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: v.get_element(1)
            (1, 0)
            sage: v.set_element((1, 1), 1)
            sage: v.get_element(1)
            (1, 1)
          """

        self._elems[vertex] = self.parent()._spaces[vertex](vector)

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the homomorphism.         #
    #                                                                         #
    ###########################################################################

    def is_zero(self):
        """
        Tests whether self is zero.

        OUTPUT:

        - bool, True is the element is the zero element, False otherwise

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: v.is_zero()
            False
            sage: w = QuiverRepElement(M)
            sage: w.is_zero()
            True
        """

        for v in self._quiver:
            if not self._elems[v].is_zero():
                return False

        return True

    def support(self):
        """
        Returns the support of self as a list.

        The support is the set of vertices to which a nonzero vector is associated.

        OUTPUT

        - list, the support

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 0), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: v.support()
            [1, 3]
        """

        sup = []
        for v in self._quiver:
            if not self._elems[v].is_zero():
                sup.append(v)

        return sup

    ###########################################################################
    #                                                                         #
    # ADDITIONAL OPERATIONS                                                   #
    #    These functions operations that are not implemented via binary       #
    #    operators.                                                           #
    #                                                                         #
    ###########################################################################

    def copy(self):
        """
        Returns a copy of self.

        EXAMPLES::

            sage: from sage.modules.quiver_module import Quiver, QuiverRep, QuiverRepElement
            sage: Q = Quiver({1:{2:['a'], 3:['b']}, 2:{3:['c']}})
            sage: spaces = dict((v, GF(3)^2) for v in Q)
            sage: M = QuiverRep(GF(3), Q, spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = QuiverRepElement(M, elems)
            sage: w = v.copy()
            sage: w.set_element((0, 0), 1)
            sage: w.get_element(1)
            (0, 0)
            sage: v.get_element(1)
            (1, 0)
        """

        if hasattr(self, '__custom_name'):
            name = self.__custom_name
        else:
            name = None
        return QuiverRepElement(self.parent(), self._elems.copy(), name)
