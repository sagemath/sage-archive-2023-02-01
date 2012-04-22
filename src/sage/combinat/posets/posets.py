# -*- coding: utf-8 -*-
r"""
Posets
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
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
#*****************************************************************************

import random
import copy
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc import deprecated_function_alias
from sage.categories.category import Category
from sage.categories.sets_cat import Sets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.posets import Posets
from sage.categories.finite_posets import FinitePosets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.graphs.all import DiGraph
from sage.combinat.posets.hasse_diagram import HasseDiagram
from sage.combinat.posets.elements import PosetElement

def Poset(data=None, element_labels=None, cover_relations=False, linear_extension=False, category = None, facade = None, key = None):
    r"""
    Construct a poset from various forms of input data.

    INPUT:

    1. A two-element list or tuple (E, R), where E is a collection of
       elements of the poset and R is the set of relations.  Elements
       of R are two-element lists/tuples/iterables.  If
       cover_relations=True, then R is assumed to be the cover
       relations of the poset. If E is empty, then E is taken to be
       the set of elements appearing in the relations R.

    2. A two-element list or tuple (E, f), where E is the set of
       elements of the poset and f is a function such that f(x,y) is
       True if x <= y and False otherwise for all pairs of elements in
       E. If cover_relations=True, then f(x,y) should be True if and
       only if x is covered by y, and False otherwise.

    3. A dictionary, list or tuple of upper covers: data[x] is an
       list of the elements that cover the element x in the poset.

       .. note::

          If data is a list or tuple of length 2, then it is handled
          by the above cases.

    4. An acyclic, loop-free and multi-edge free DiGraph. If
       cover_relations is True, then the edges of the digraph
       correspond to cover relations in the poset. If cover_relations
       is False, then the cover relations are computed.

    5. A previously constructed poset (the poset itself is returned).


    - ``element_labels`` -- (default: None) an optional list or
      dictionary of objects that label the poset elements.

    - ``cover_relations`` -- a boolean (default: False); whether the
      data can be assumed to describe a directed acyclic graph whose
      arrows are cover relations; otherwise, the cover relations are
      first computed.

    - ``linear_extension`` -- a boolean (default: False); whether to
      use the provided list of elements as default linear extension
      for the poset; otherwise a linear extension is computed.

    OUTPUT:

        FinitePoset -- an instance of the FinitePoset class.

    If ``category`` is specified, then the poset is created in this
    category instead of :class:`FinitePosets`.

    .. seealso:: :class:`Posets`, :class:`~sage.categories.posets.Posets`, :class:`FinitePosets`

    EXAMPLES:

    1. Elements and cover relations::

          sage: elms = [1,2,3,4,5,6,7]
          sage: rels = [[1,2],[3,4],[4,5],[2,5]]
          sage: Poset((elms, rels), cover_relations = True)
          Finite poset containing 7 elements

       Elements and non-cover relations::

          sage: elms = [1,2,3,4]
          sage: rels = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
          sage: P = Poset( [elms,rels] ,cover_relations=False); P
          Finite poset containing 4 elements
          sage: P.cover_relations()
          [[1, 2], [2, 3], [3, 4]]

    2. Elements and function: the standard permutations of [1, 2, 3, 4]
       with the Bruhat order::

          sage: elms = Permutations(4)
          sage: fcn = lambda p,q : p.bruhat_lequal(q)
          sage: Poset((elms, fcn))
          Finite poset containing 24 elements

       With a function that identifies the cover relations: the set
       partitions of {1, 2, 3} ordered by refinement::

          sage: elms = SetPartitions(3)
          sage: def fcn(A, B):
          ...     if len(A) != len(B)+1:
          ...         return False
          ...     for a in A:
          ...         if not any(set(a).issubset(b) for b in B):
          ...             return False
          ...     return True
          sage: Poset((elms, fcn), cover_relations=True)
          Finite poset containing 5 elements

    3. A dictionary of upper covers::

          sage: Poset({'a':['b','c'], 'b':['d'], 'c':['d'], 'd':[]})
          Finite poset containing 4 elements

       A list of upper covers::

          sage: Poset([[1,2],[4],[3],[4],[]])
          Finite poset containing 5 elements

       A list of upper covers and a dictionary of labels::

          sage: elm_labs = {0:"a",1:"b",2:"c",3:"d",4:"e"}
          sage: P = Poset([[1,2],[4],[3],[4],[]],elm_labs)
          sage: P.list()
          [a, b, c, d, e]

       .. warning::

         The special case where the argument data is a list or tuple of
         length 2 is handled by the above cases. So you cannot use this
         method to input a 2-element poset.

    4. An acyclic DiGraph.

       ::

          sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
          sage: Poset(dag)
          Finite poset containing 6 elements

       Any directed acyclic graph without loops or multiple edges, as long
       as cover_relations=False::

          sage: dig = DiGraph({0:[2,3], 1:[3,4,5], 2:[5], 3:[5], 4:[5]})
          sage: dig.allows_multiple_edges()
          False
          sage: dig.allows_loops()
          False
          sage: dig.transitive_reduction() == dig
          False
          sage: Poset(dig, cover_relations=False)
          Finite poset containing 6 elements
          sage: Poset(dig, cover_relations=True)
          Traceback (most recent call last):
          ...
          ValueError: Hasse diagram is not transitively reduced.

    .. rubric:: Default Linear extension

    Every poset `P` obtained with ``Poset`` comes equipped with a
    default linear extension, which is also used for enumerating
    its elements. By default, this linear extension is computed,
    and has no particular significance::

        sage: P = Poset((divisors(12), attrcall("divides")))
        sage: P.list()
        [1, 2, 4, 3, 6, 12]
        sage: P.linear_extension()
        [1, 2, 4, 3, 6, 12]

    You may enforce a specific linear extension using the
    ``linear_extension`` option::

        sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
        sage: P.list()
        [1, 2, 3, 4, 6, 12]
        sage: P.linear_extension()
        [1, 2, 3, 4, 6, 12]

    Depending on popular request, ``Poset`` might eventually get
    modified to always use the provided list of elements as
    default linear extension, when it is one.

    .. seealso:: :meth:`FinitePoset.linear_extensions`

    .. rubric:: Facade posets

    By default, the elements of a poset are wrapped so as to make them
    aware that they belong to that poset::

        sage: P = Poset(DiGraph({'d':['c','b'],'c':['a'],'b':['a']}))
        sage: d,c,b,a = list(P)
        sage: a.parent() is P
        True

    This allows for comparing elements according to `P`::

        sage: c < a
        True

    As an experimental feature, one can construct instead facade posets::

        sage: P = Poset(DiGraph({'d':['c','b'],'c':['a'],'b':['a']}),
        ...             facade = True)

    In this example, the elements of the poset remain plain strings::

        sage: d,c,b,a = list(P)
        sage: type(a)
        <type 'str'>

    Of course, those strings are not aware of `P`. So to compare two
    such strings, one needs to query `P`::

        sage: a < b
        True
        sage: P.lt(a,b)
        False

    which models the usual mathematical notation `a <_P b`.

    Most operations seem to still work, but at this point there is no
    guarantee whatsoever::

        sage: P.list()
        ['d', 'b', 'c', 'a']
        sage: P.principal_order_ideal('a')
        ['d', 'b', 'c', 'a']
        sage: P.principal_order_ideal('b')
        ['d', 'b']
        sage: P.principal_order_ideal('d')
        ['d']
        sage: TestSuite(P).run()

    .. warning::

        :class:`DiGraph` is used to construct the poset, and the
        vertices of a :class:`DiGraph` are converted to plain Python
        :class:`int`'s if they are :class:`Integer`'s::

            sage: G = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: type(G.vertices()[0])
            <type 'int'>

        This is worked around by systematically converting back the
        vertices of a poset to :class:`Integer`'s if they are
        :class:`int`'s::

            sage: P = Poset((divisors(15), attrcall("divides")))
            sage: type(P.an_element().element)
            <type 'sage.rings.integer.Integer'>

            sage: P = Poset((divisors(15), attrcall("divides")), facade=True)
            sage: type(P.an_element())
            <type 'sage.rings.integer.Integer'>

        This may be abusive::

            sage: P = Poset((range(5), operator.le), facade = True)
            sage: P.an_element().parent()
            Integer Ring

    .. rubric:: Unique representation

    As most parents, :class:`Poset` have unique representation (see
    :class:`UniqueRepresentation`. Namely if two posets are created
    from two equal data, then they are not only equal but actually
    identical::

        sage: data1 = [[1,2],[3],[3]]
        sage: data2 = [[1,2],[3],[3]]
        sage: P1 = Poset(data1)
        sage: P2 = Poset(data2)
        sage: P1 == P2
        True
        sage: P1 is P2
        True

    In situations where this behaviour is not desired, one can use the
    ``key`` option::

        sage: P1 = Poset(data1, key = "foo")
        sage: P2 = Poset(data2, key = "bar")
        sage: P1 is P2
        False
        sage: P1 == P2
        False

    ``key`` can be any hashable value and is passed down to
    :class:`UniqueRepresentation`. It is otherwise ignored by the
    poset constructor.


    TESTS::

        sage: P = Poset([[1,2],[3],[3]])
        sage: type(hash(P))
        <type 'int'>
    """
    #Convert data to a DiGraph
    elements = None
    D = {}
    if isinstance(data, FinitePoset):
        if element_labels is None and category is None and facade is None:
            return data
        else:
            return FinitePoset(data, data._elements, category = category, facade = facade)
    elif data is None: # type 0
        D = DiGraph()
    elif isinstance(data, DiGraph): # type 4
        D = copy.deepcopy(data)
    elif isinstance(data, dict): # type 3: dictionary of upper covers
        D = DiGraph(data)
    elif isinstance(data,(list,tuple)): # types 1, 2, 3 (list/tuple)
        if len(data) == 2: # types 1 or 2
            if callable(data[1]): # type 2
                elements, function = data
                relations = []
                for x in elements:
                    for y in elements:
                        if function(x,y) is True:
                            relations.append([x,y])
            else: # type 1
                elements, relations = data
                # check that relations are relations
                for r in relations:
                    try:
                        u, v = r
                    except ValueError:
                        raise TypeError, "not a list of relations"
            D = DiGraph()
            D.add_vertices(elements)
            D.add_edges(relations)
        elif len(data) > 2:
            # type 3, list/tuple of upper covers
            D = DiGraph(dict([[Integer(i),data[i]] for i in range(len(data))]))
        else:
            raise ValueError, "not valid poset data."

    # DEBUG: At this point D should be a DiGraph.
    if not isinstance(D,DiGraph):
        raise TypeError, "BUG: D should be a digraph."

    # Determine cover relations, if necessary.
    if cover_relations is False:
        D = D.transitive_reduction()

    # Check that the digraph does not contain loops, multiple edges
    # and is transitively reduced.
    if D.has_loops():
        raise ValueError, "Hasse diagram contains loops."
    elif D.has_multiple_edges():
        raise ValueError, "Hasse diagram contains multiple edges."
    elif cover_relations is True and not D.is_transitively_reduced():
        raise ValueError, "Hasse diagram is not transitively reduced."

    if linear_extension and elements is not None:
        lin_ext = list(elements)
    else:
        # Compute a linear extension of the poset (a topological sort).
        try:
            lin_ext = D.topological_sort()
        except:
            raise ValueError, "Hasse diagram contains cycles."

    # Relabel using the linear_extension.
    # So range(len(D)) becomes a linear extension of the poset.
    rdict = dict([[lin_ext[i],i] for i in range(len(lin_ext))])
    D.relabel(rdict)

    # Set element labels.
    if element_labels is None:
        elements = lin_ext
        # Work around the fact that, currently, when a DiGraph is
        # created with Integer's as vertices, those vertices are
        # converted to plain int's. This is a bit abusive.
        elements = [ Integer(i) if isinstance(i,int) else i for i in elements ]
    else:
        elements = [element_labels[z] for z in lin_ext]

    return FinitePoset(D,elements, category = category, facade = facade, key = key)

class FinitePoset(UniqueRepresentation, Parent):
    r"""
    Constructs a (finite) `n`-element poset from a set of elements and a
    directed acyclic graph or poset.

    INPUT:

    - ``hasse_diagram`` -- an instance of this class (``FinitePoset``),
      or a digraph that is transitively-reduced, acyclic, loop-free,
      multiedge-free, and with vertices indexed by ``range(n)``. We also
      assume that ``range(n)`` is a linear extension of the poset.

    - ``elements`` - an optional list of elements, with ``element[i]``
      corresponding to vertex ``i``. If ``elements`` is ``None``, then it is
      set to be the vertex set of the digraph.

    - ``category`` -- :class:`FinitePosets`, or a subcategory thereof.

    - ``facade`` -- a boolean or ``None`` (the default).

    - ``key`` -- any hashable value (default: ``None``).

    EXAMPLES::

        sage: uc = [[2,3], [], [1], [1], [1], [3,4]]
        sage: from sage.combinat.posets.posets import FinitePoset
        sage: P = FinitePoset(DiGraph(dict([[i,uc[i]] for i in range(len(uc))]))); P
        Finite poset containing 6 elements
        sage: P.cover_relations()
        [[0, 2], [0, 3], [2, 1], [3, 1], [4, 1], [5, 3], [5, 4]]
        sage: TestSuite(P).run()
        sage: P.category()
        Category of finite posets
        sage: P.__class__
        <class 'sage.combinat.posets.posets.FinitePoset_with_category'>

        sage: Q = sage.combinat.posets.posets.FinitePoset(P); Q
        Finite poset containing 6 elements

        sage: Q is P
        True

    We keep the same underlying hasse diagram, but change the elements::

        sage: Q = sage.combinat.posets.posets.FinitePoset(P, elements=[1,2,3,4,5,6]); Q
        Finite poset containing 6 elements
        sage: Q.cover_relations()
        [[1, 3], [1, 4], [3, 2], [4, 2], [5, 2], [6, 4], [6, 5]]

    We test the facade argument::

        sage: P = Poset(DiGraph({'a':['b'],'b':['c'],'c':['d']}))
        sage: P.category()
        Category of finite posets
        sage: parent(P[0]) is P
        True

        sage: Q = Poset(DiGraph({'a':['b'],'b':['c'],'c':['d']}), facade = True)
        sage: Q.category()
        Category of facade finite posets
        sage: parent(Q[0]) is str
        True
        sage: TestSuite(Q).run(skip = ['_test_an_element']) # is_parent_of is not yet implemented

    Changing a non facade poset to a facade poset::

        sage: PQ = Poset(P, facade = True)
        sage: PQ.category()
        Category of facade finite posets
        sage: parent(PQ[0]) is str
        True
        sage: PQ is Q
        True

    Changing a facade poset to a non facade poset::

        sage: QP = Poset(Q, facade = False)
        sage: QP.category()
        Category of finite posets
        sage: parent(QP[0]) is QP
        True

    .. note::

       A class that inherits from this class needs to define
       ``Element``. This is the class of the elements that the inheriting
       class contains. For example, for this class, ``FinitePoset``,
       ``Element`` is ``PosetElement``.  It can also define ``_dual_class`` which
       is the class of dual posets of this
       class. E.g. ``FiniteMeetSemilattice._dual_class`` is
       ``FiniteJoinSemilattice``.

    TESTS:

    Equality is derived from :class:`UniqueRepresentation`. We check that this
    gives consistent results::

        sage: P = Poset([[1,2],[3],[3]])
        sage: P == P
        True
        sage: Q = Poset([[1,2],[],[1]])
        sage: Q == P
        False
        sage: p1, p2 = Posets(2).list()
        sage: p2 == p1, p1 != p2
        (False, True)
        sage: [[p1.__eq__(p2) for p1 in Posets(2)] for p2 in Posets(2)]
        [[True, False], [False, True]]
        sage: [[p2.__eq__(p1) for p1 in Posets(2)] for p2 in Posets(2)]
        [[True, False], [False, True]]
        sage: [[p2 == p1 for p1 in Posets(3)] for p2 in Posets(3)]
        [[True, False, False, False, False], [False, True, False, False, False], [False, False, True, False, False], [False, False, False, True, False], [False, False, False, False, True]]

        sage: [[p1.__ne__(p2) for p1 in Posets(2)] for p2 in Posets(2)]
        [[False, True], [True, False]]
        sage: P = Poset([[1,2,4],[3],[3]])
        sage: Q = Poset([[1,2],[],[1],[4]])
        sage: P != Q
        True
        sage: P != P
        False
        sage: Q != Q
        False
        sage: [[p1.__ne__(p2) for p1 in Posets(2)] for p2 in Posets(2)]
        [[False, True], [True, False]]
    """

    @staticmethod
    def __classcall__(cls, hasse_diagram, elements = None, category = None, facade = None, key = None):
        """
        Normalizes the arguments passed to the constructor

        TESTS::

            sage: P = sage.combinat.posets.posets.FinitePoset(DiGraph())
            sage: type(P)
            <class 'sage.combinat.posets.posets.FinitePoset_with_category'>
            sage: TestSuite(P).run()

        See also the extensive tests in the class documentation
        """
        assert isinstance(hasse_diagram, (FinitePoset, DiGraph))
        if isinstance(hasse_diagram, FinitePoset):
            if elements is None:
                elements = hasse_diagram._elements
            if category is None:
                category = hasse_diagram.category()
                if facade is False and category.is_subcategory(Sets().Facades()):
                    # We need to remove Sets().Facades() from the category
                    # This is fragile ...
                    from sage.categories.category import JoinCategory
                    assert isinstance(category, JoinCategory)
                    categories = list(category.super_categories())
                    categories.remove(Sets().Facades())
                    category = Category.join(categories)
            if facade is None:
                facade = hasse_diagram in Sets().Facades()
            hasse_diagram = hasse_diagram._hasse_diagram
        else:
            hasse_diagram = HasseDiagram(hasse_diagram)
            if elements is None:
                elements = hasse_diagram.vertices()
            if facade is None:
                facade = False
        elements = tuple(elements)
        category = Category.join([FinitePosets().or_subcategory(category), FiniteEnumeratedSets()])
        return super(FinitePoset, cls).__classcall__(cls, hasse_diagram = hasse_diagram, elements = elements,
                                                     category = category, facade = facade, key = key)

    def __init__(self, hasse_diagram, elements, category, facade, key):
        """
        EXAMPLES::

            sage: P = Poset(DiGraph({'a':['b'],'b':['c'],'c':['d']}))
            sage: type(P)
            <class 'sage.combinat.posets.posets.FinitePoset_with_category'>

        The internal data structure currently consists of:

        - the Hasse diagram of the poset, represented by a DiGraph
          with vertices labelled 0,...,n-1 according to a linear
          extension of the poset (that is if `i \mapsto j` is an edge
          then `i<j`), together with some extra methods (see
          :class:`sage.combinat.posets.hasse_diagram.HasseDiagram`)::

            sage: P._hasse_diagram
            Hasse diagram of a poset containing 4 elements
            sage: P._hasse_diagram.cover_relations()
            [(0, 1), (1, 2), (2, 3)]

        - a tuple of the original elements, not wrapped as elements of
          ``self`` (but see also ``P._list``)::

            sage: P._elements
            ('a', 'b', 'c', 'd')

          ``P._elements[i]`` gives the element of ``P`` corresponding
          to the vertex ``i``

        - a dictionary mapping back elements to vertices::

            sage: P._element_to_vertex_dict
            {'a': 0, 'c': 2, 'b': 1, 'd': 3}

        - and a boolean stating whether the poset is a facade poset::

            sage: P._is_facade
            False

        This internal data structure is subject to change at any
        point. Do not break encapsulation!

        TESTS::

            sage: TestSuite(P).run()

        See also the extensive tests in the class documentation.
        """
        Parent.__init__(self, category = category, facade = facade)
        self._hasse_diagram = hasse_diagram
        self._elements = elements
        self._element_to_vertex_dict = dict( (elements[i], i) for i in range(len(elements)) )
        self._is_facade = facade

    @lazy_attribute
    def _list(self):
        """
        The list of the elements of ``self``, each wrapped to have
        ``self`` as parent

        EXAMPLES::

            sage: P = Poset(DiGraph({'a':['b'],'b':['c'],'c':['d']}))
            sage: L = P._list; L
            (a, b, c, d)
            sage: type(L[0])
            <class 'sage.combinat.posets.elements.FinitePoset_with_category.element_class'>
            sage: L[0].parent() is P
            True

        Constructing them once for all makes future conversions
        between the vertex id and the element faster. This also
        ensures unique representation of the elements of this poset,
        which could be used to later speed up certain operations
        (equality test, ...)
        """
        if self._is_facade:
            return self._elements
        else:
            return tuple(self.element_class(self, self._elements[vertex], vertex)
                         for vertex in range(len(self._elements)))

    # This defines the type (class) of elements of poset.
    Element = PosetElement

    def _element_to_vertex(self, element):
        """
        Given an element of the poset (wrapped or not), returns the
        corresponding vertex of the Hasse diagram.

        EXAMPLES::

            sage: P = Poset((divisors(15), attrcall("divides")))
            sage: list(P)
            [1, 3, 5, 15]
            sage: x = P(5)
            sage: P._element_to_vertex(x)
            2

        The same with a non-wrapped element of `P`::

            sage: P._element_to_vertex(5)
            2

        TESTS::

            sage: P = Poset((divisors(15), attrcall("divides")))

        Testing for wrapped elements::

            sage: all(P._vertex_to_element(P._element_to_vertex(x)) is x for x in P)
            True

        Testing for non-wrapped elements::

            sage: all(P._vertex_to_element(P._element_to_vertex(x)) is P(x) for x in divisors(15))
            True

        Testing for non-wrapped elements for a facade poset::

            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)
            sage: all(P._vertex_to_element(P._element_to_vertex(x)) is x for x in P)
            True
        """
        if isinstance(element, self.element_class) and element.parent() is self:
            return element.vertex
        else:
            try:
                return self._element_to_vertex_dict[element]
            except KeyError:
                raise ValueError, "element (=%s) not in poset"%element

    def _vertex_to_element(self,vertex):
        """
        Returns the element corresponding to a vertex of the Hasse diagram.

        It is wrapped if ``self`` is not a facade poset.

        EXAMPLES::

            sage: P = Poset((divisors(15), attrcall("divides")))
            sage: x = P._vertex_to_element(2)
            sage: x
            5
            sage: x.parent() is P
            True

            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)
            sage: x = P._vertex_to_element(2)
            sage: x
            5
            sage: x.parent() is ZZ
            True
        """
        return self._list[vertex]

    def unwrap(self, element):
        """
        Unwraps an element of this poset

        INPUT:

        - ``element`` -- an element of ``self``

        EXAMPLES::

            sage: P = Poset((divisors(15), attrcall("divides")))
            sage: x = P.an_element(); x
            1
            sage: x.parent()
            Finite poset containing 4 elements
            sage: P.unwrap(x)
            1
            sage: P.unwrap(x).parent()
            Integer Ring

        For a non facade poset, this is equivalent to using the
        ``.element`` attribute::

            sage: P.unwrap(x) is x.element
            True

        For a facade poset, this does nothing::

            sage: P = Poset((divisors(15), attrcall("divides")), facade=True)
            sage: x = P.an_element()
            sage: P.unwrap(x) is x
            True

        This method is useful in code where we don't know if ``P`` is
        a facade poset or not.
        """
        if self._is_facade:
            return element
        else:
            return element.element

    def __contains__(self, x):
        r"""
        Returns True if x is an element of the poset.

        TESTS::

            sage: from sage.combinat.posets.posets import FinitePoset
            sage: P5 = FinitePoset(DiGraph({(5,):[(4,1),(3,2)], \
                    (4,1):[(3,1,1),(2,2,1)], \
                    (3,2):[(3,1,1),(2,2,1)], \
                    (3,1,1):[(2,1,1,1)], \
                    (2,2,1):[(2,1,1,1)], \
                    (2,1,1,1):[(1,1,1,1,1)], \
                    (1,1,1,1,1):[]}))
            sage: x = P5.list()[3]
            sage: x in P5
            True

        For the sake of speed, an element with the right class and
        parent is assumed to be in this parent. This can possibly be
        counterfeited by feeding garbage to the constructor::

            sage: x = P5.element_class(P5, "a", 5)
            sage: x in P5
            True
        """
        if isinstance(x, self.element_class):
            return x.parent() is self
        return x in self._element_to_vertex_dict

    is_parent_of = __contains__

    def _element_constructor_(self, element):
        """
        Constructs an element of ``self``

        EXAMPLES::

            sage: from sage.combinat.posets.posets import FinitePoset
            sage: P = FinitePoset(DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]}))
            sage: P(5)
            5
            sage: Q = FinitePoset(DiGraph({5:[2,3], 1:[3,4], 2:[0], 3:[0], 4:[0]}))
            sage: Q(5)
            5

        Accessing the n-th element of ``self`` as ``P(i)`` is deprecated::

            sage: P(5) == P(-1)
            doctest:...: DeprecationWarning: Accessing the i-th element of a poset as P(i) is deprecated. Please use P[i]
            True
            sage: Q(5) == Q(-1)
            True
            sage: R = FinitePoset(DiGraph({'a':['b','c'], 'b':['d'], 'c':['d'], 'd':[]}))
            sage: R(0)
            a
            sage: R('a') == R(0)
            True
            sage: R('d') == R(-1)
            True

        Please use instead ``P[i]``::

            sage: R('a') == R[0]
            True
            sage: R('d') == R[-1]
            True

        TESTS::

            sage: P = Poset(DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]}))
            sage: all(P(x) is x for x in P)
            True
            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)
            sage: all(P(x) is x for x in P)
            True
        """
        try:
            return self._list[self._element_to_vertex_dict[element]]
        except KeyError:
            if isinstance(element,Integer):
                import sage.misc.misc
                sage.misc.misc.deprecation("Accessing the i-th element of a poset as P(i) is deprecated. Please use P[i]")
                if element > -1:
                    return self.element_class(self, \
                        self._elements[element], element)
                else:
                    return self.element_class(self, \
                        self._elements[element], self.cardinality()+element)
            else:
                raise ValueError, "%s is not an element of this poset"%type(element)

    def __call__(self, element):
        """
        Creates elements of this poset

        This overrides the generic call method for all parents
        :meth:`Parent.__call__`, as a work around to allow for facade
        posets over plain Python objects (e.g. Python's
        int's). Indeed, the default __call__ method expects the result
        of :meth:`_element_constructor_` to be a Sage element (see
        :meth:`sage.structure.coerce_maps.DefaultConvertMap_unique._call_`)::

            sage: P = Poset(DiGraph({'d':['c','b'],'c':['a'],'b':['a']}),
            ...             facade = True)
            sage: P('a')              # indirect doctest
            'a'
            sage: TestSuite(P).run()
            sage: P = Poset(((False, True), operator.eq), facade = True)
            sage: P(True)
            1
        """
        if self._is_facade and element in self._element_to_vertex_dict:
            return element
        return super(FinitePoset, self).__call__(element)

    def hasse_diagram(self, wrapped = True):
        """
        Returns the Hasse diagram of ``self`` as a Sage :class:`DiGraph`.

        .. todo:: should the vertices of the diagram have the poset as
                  parent?

        EXAMPLES::

            sage: Q = Poset({5:[2,3], 1:[3,4], 2:[0], 3:[0], 4:[0]})
            sage: Q.hasse_diagram()
            Digraph on 6 vertices

            sage: P = Poset({'a':['b'],'b':['d'],'c':['d'],'d':['f'],'e':['f'],'f':[]})
            sage: H = P.hasse_diagram()
            sage: P.cover_relations()
            [[e, f], [c, d], [a, b], [b, d], [d, f]]
            sage: H.edges()
            [(a, b, None), (c, d, None), (b, d, None), (e, f, None), (d, f, None)]

            sage: P = Poset((divisors(15), attrcall("divides")))
            sage: H = P.hasse_diagram()
            sage: H.vertices()
            [1, 5, 3, 15]
            sage: H.edges()
            [(1, 3, None), (1, 5, None), (5, 15, None), (3, 15, None)]
            sage: H.set_latex_options(format = "dot2tex")   # optional
            sage: view(H, tight_page=True, pdflatex = True) # optional
        """
        return DiGraph(self._hasse_diagram).relabel(self._list, inplace = False)

    def _latex_(self):
        r"""
        Returns a latex method for the poset.

        EXAMPLES::

            sage: P = Poset(([1,2], [[1,2]]), cover_relations = True)
            sage: print P._latex_() #optional - dot2tex graphviz
            \begin{tikzpicture}
            %
            \useasboundingbox (0,0) rectangle (5.0cm,5.0cm);
            %
            \definecolor{cv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv0}{rgb}{1.0,1.0,1.0}
            \definecolor{clv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv1}{rgb}{1.0,1.0,1.0}
            \definecolor{clv1}{rgb}{0.0,0.0,0.0}
            \definecolor{cv0v1}{rgb}{0.0,0.0,0.0}
            %
            \Vertex[style={minimum size=1.0cm,draw=cv0,fill=cfv0,text=clv0,shape=circle},LabelOut=false,L=\hbox{$1$},x=0.0cm,y=0.0cm]{v0}
            \Vertex[style={minimum size=1.0cm,draw=cv1,fill=cfv1,text=clv1,shape=circle},LabelOut=false,L=\hbox{$2$},x=5.0cm,y=5.0cm]{v1}
            %
            \Edge[lw=0.1cm,style={post, bend right,color=cv0v1,},](v0)(v1)
            %
            \end{tikzpicture}
        """
        return self.hasse_diagram()._latex_()

    def _repr_(self):
        r"""
        Returns a string representation of the poset.

        TESTS::

            sage: partitions_of_five = {(5,):[(4,1),(3,2)], \
                    (4,1):[(3,1,1),(2,2,1)], \
                    (3,2):[(3,1,1),(2,2,1)], \
                    (3,1,1):[(2,1,1,1)], \
                    (2,2,1):[(2,1,1,1)], \
                    (2,1,1,1):[(1,1,1,1,1)], \
                    (1,1,1,1,1):[]}
            sage: P5 = Poset(partitions_of_five)
            sage: P5._repr_()
            'Finite poset containing 7 elements'
        """
        return "Finite poset containing %s elements"%self._hasse_diagram.order()

    def __iter__(self):
        """
        Iterates through the elements of a linear extension of the poset.

        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: sorted(D.__iter__())
            [0, 1, 2, 3, 4]
        """
        return iter(self._list)

    def linear_extension(self, linear_extension=None, check=True):
        """
        Returns a linear extension of this poset.

        INPUT:

        - ``linear_extension`` -- a list of the elements of ``self`` (default: ``None``)
        - ``check`` -- a boolean (default: True);
          whether to check that ``linear_extension`` is indeed a
          linear extension of ``self``.

        .. seealso:: :meth:`is_linear_extension`, :meth:`linear_extensions`

        EXAMPLES::

            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)

        Without optional argument, the default linear extension of the
        poset is returned, as a plain list::

            sage: P.linear_extension()
            [1, 3, 5, 15]

        Otherwise, a full-featured linear extension is constructed
        as an element of ``P.linear_extensions()``::

            sage: l = P.linear_extension([1,5,3,15]); l
            [1, 5, 3, 15]
            sage: type(l)
            <class 'sage.combinat.posets.linear_extensions.LinearExtensionsOfPoset_with_category.element_class'>
            sage: l.parent()
            The set of all linear extensions of Finite poset containing 4 elements

        By default, the linear extension is checked for correctness::

            sage: l = P.linear_extension([1,3,15,5])
            Traceback (most recent call last):
            ...
            ValueError: [1, 3, 15, 5] is not a linear extension of Finite poset containing 4 elements

        This can be disabled (at your own risks!) with::

            sage: P.linear_extension([1,3,15,5], check=False)
            [1, 3, 15, 5]

        .. todo::

            - Is it acceptable to have those two features for a single method?

            - In particular, we miss a short idiom to get the default
              linear extension
        """
        if linear_extension is not None:
            return self.linear_extensions()(linear_extension, check=check)
        else:
            # TODO: do we care whether this is a list or tuple?
            return list(self._list)

    @cached_method
    def linear_extensions(self, facade=False):
        """
        Returns the enumerated set of all the linear extensions of this poset

        INPUT:

        - ``facade`` -- a boolean (default: False);
          whether to return the linear extensions as plain lists

        .. seealso:: :meth:`linear_extension`, :meth:`is_linear_extension`

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: L = P.linear_extensions(); L
            The set of all linear extensions of Finite poset containing 6 elements
            sage: l = L.an_element(); l
            [1, 2, 3, 4, 6, 12]
            sage: L.cardinality()
            5
            sage: L.list()
            [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]

        Each element is aware that it is a linear extension of `P`::

            sage: type(l.parent())
            <class 'sage.combinat.posets.linear_extensions.LinearExtensionsOfPoset_with_category'>

        With ``facade=True``, the elements of ``L`` are plain lists instead::

            sage: L = P.linear_extensions(facade=True)
            sage: l = L.an_element()
            sage: type(l)
            <type 'list'>

        .. warning::

            In Sage <= 4.8, this function used to return a plain list
            of lists. To recover the previous functionality, please
            use::

                sage: L = list(P.linear_extensions(facade=True)); L
                [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]
                sage: type(L[0])
                <type 'list'>

        .. todo::

            The ``facade`` option is not yet fully functional::

                 sage: L = P.linear_extensions(facade=True); L
                 The set of all linear extensions of Finite poset containing 6 elements
                 sage: L([1, 2, 3, 4, 6, 12])
                 Traceback (most recent call last):
                 ...
                 TypeError: Cannot convert list to sage.structure.element.Element

        TESTS::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: list(D.linear_extensions())
            [[0, 1, 2, 3, 4], [0, 1, 2, 4, 3], [0, 2, 1, 3, 4], [0, 2, 1, 4, 3], [0, 2, 4, 1, 3]]

        """
        from linear_extensions import LinearExtensionsOfPoset
        return LinearExtensionsOfPoset(self, facade = facade)

    def is_linear_extension(self, l):
        """
        Returns whether ``l`` is a linear extension of ``self``

        INPUT:

        - ``l`` -- a list (or iterable) containing all of the elements of ``self`` exactly once

        .. seealso:: :meth:`linear_extension`, :meth:`linear_extensions`

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), facade=True, linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: P.is_linear_extension([1, 2, 4, 3, 6, 12])
            True
            sage: P.is_linear_extension([1, 2, 4, 6, 3, 12])
            False

            sage: [p for p in Permutations(list(P)) if P.is_linear_extension(p)]
            [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]
            sage: list(P.linear_extensions())
            [[1, 2, 3, 4, 6, 12], [1, 2, 3, 6, 4, 12], [1, 2, 4, 3, 6, 12], [1, 3, 2, 4, 6, 12], [1, 3, 2, 6, 4, 12]]

        .. note:: this is used and systematically tested in :class:`~sage.combinat.posets.linear_extensions.LinearExtensionsOfPosets`

        """
        index = { x:i for (i,x) in enumerate(l) }
        return all(index[i] < index[j] for (i,j) in self.cover_relations())

    def list(self):
        """
        List the elements of the poset. This just returns the result
        of :meth:`linear_extension`.

        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: D.list()
            [0, 1, 2, 3, 4]
            sage: type(D.list()[0])
            <class 'sage.combinat.posets.elements.FinitePoset_with_category.element_class'>
        """
        return list(self.linear_extension())

    def plot(self, label_elements=True, element_labels=None,
             label_font_size=12,label_font_color='black',
             vertex_size=300, vertex_colors=None,
             layout = 'acyclic',
             **kwds):
        """
        Returns a Graphic object corresponding the Hasse diagram of the
        poset. Optionally, it is labelled.

        INPUT:


        -  ``label_elements`` - whether to display element
           labels

        -  ``element_labels`` - a dictionary of element
           labels


        EXAMPLES::

            sage: D = Poset({ 1:[2,3], 2:[4], 3:[4,5] })
            sage: D.plot(label_elements=False)
            sage: D.plot()
            sage: type(D.plot())
            <class 'sage.plot.graphics.Graphics'>
            sage: elm_labs = {1:'a', 2:'b', 3:'c', 4:'d', 5:'e'}
            sage: D.plot(element_labels=elm_labs)

        ::

            sage: P = Poset({})
            sage: P.plot()

        ::

            sage: P = Poset(DiGraph('E@ACA@?'))
            sage: P.plot()

        TESTS:

        We check that ``label_elements`` and is honored::

            sage: def get_plot_labels(P): return sorted(t.string for t in P if isinstance(t, sage.plot.text.Text))
            sage: P1 = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: P2 = Poset({ 0:[1,2], 1:[3], 2:[3,4] }, facade=True)
            sage: get_plot_labels(P1.plot(label_elements=False))
            []
            sage: get_plot_labels(P1.plot(label_elements=True))
            ['0', '1', '2', '3', '4']
            sage: element_labels = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e'}
            sage: get_plot_labels(P1.plot(element_labels=element_labels))
            ['a', 'b', 'c', 'd', 'e']
            sage: get_plot_labels(P2.plot(element_labels=element_labels))
            ['a', 'b', 'c', 'd', 'e']

        """
        graph = self.hasse_diagram()
        if label_elements and element_labels is not None:
            graph = graph.relabel(dict((self(element),label) for (element,label) in element_labels.items()), inplace = False)
        return graph.plot(vertex_labels=label_elements,
                          label_font_size=label_font_size,
                          label_font_color=label_font_color,
                          vertex_size=vertex_size,
                          vertex_colors=vertex_colors,
                          layout = layout,
                          **kwds)

    def show(self, label_elements=True, element_labels=None,
            label_font_size=12,label_font_color='black',
            vertex_size=300, vertex_colors=None, layout='acyclic', **kwds):
        """
        Shows the Graphics object corresponding the Hasse diagram of the
        poset. Optionally, it is labelled.

        INPUT:


        -  ``label_elements`` - whether to display element
           labels

        -  ``element_labels`` - a dictionary of element
           labels


        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: D.plot(label_elements=False)
            sage: D.show()
            sage: elm_labs = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e'}
            sage: D.show(element_labels=elm_labs)
        """
        self.plot(label_elements=label_elements, element_labels=element_labels,
            label_font_size=label_font_size,label_font_color=label_font_color,
            vertex_size=vertex_size, vertex_colors=vertex_colors, layout=layout).show(**kwds)

    def level_sets(self):
        """
        Returns a list l such that l[i+1] is the set of minimal elements of
        the poset obtained by removing the elements in l[0], l[1], ...,
        l[i].

        EXAMPLES::

            sage: P = Poset({0:[1,2],1:[3],2:[3],3:[]})
            sage: [len(x) for x in P.level_sets()]
            [1, 2, 1]

        ::

            sage: Q = Poset({0:[1,2], 1:[3], 2:[4], 3:[4]})
            sage: [len(x) for x in Q.level_sets()]
            [1, 2, 1, 1]
        """
        return [map(self._vertex_to_element, level) for level in
                self._hasse_diagram.level_sets()]

    def cover_relations(self,element=None):
        """
        Returns the list of pairs [u,v] of elements of the poset such that
        u v is a cover relation (that is, u v and there does not exist z
        such that u z v).

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: Q.cover_relations()
            [[1, 2], [0, 2], [2, 3], [3, 4]]
        """
        return [c for c in self.cover_relations_iterator()]

    def cover_relations_iterator(self):
        """
        Returns an iterator for the cover relations of the poset.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.cover_relations_iterator())
            <type 'generator'>
            sage: [z for z in Q.cover_relations_iterator()]
            [[1, 2], [0, 2], [2, 3], [3, 4]]
        """
        for u,v,l in self._hasse_diagram.edge_iterator():
            yield map(self._vertex_to_element,(u,v))

    def relations(self):
        r"""
        Returns a list of all relations of the poset.

        OUTPUT:

        A list of pairs (each pair is a list), where the first element
        of the pair is less than or equal to the second element.

        Pairs are produced in a rough sort of lexicographic order,
        where earlier elements are from lower levels of the poset.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: Q.relations()
            [[1, 1], [1, 2], [1, 3], [1, 4], [0, 0], [0, 2], [0, 3], [0, 4], [2, 2], [2, 3], [2, 4], [3, 3], [3, 4], [4, 4]]

        AUTHOR:

        - Rob Beezer (2011-05-04)
        """
        return list(self.relations_iterator())

    def relations_iterator(self):
        r"""
        Returns an iterator for all the relations of the poset.

        OUTPUT:

        A generator that produces pairs (each pair is a list), where the
        first element of the pair is less than or equal to the second element.

        Pairs are produced in a rough sort of lexicographic order,
        where earlier elements are from lower levels of the poset.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.relations_iterator())
            <type 'generator'>
            sage: [z for z in Q.relations_iterator()]
            [[1, 1], [1, 2], [1, 3], [1, 4], [0, 0], [0, 2], [0, 3], [0, 4], [2, 2], [2, 3], [2, 4], [3, 3], [3, 4], [4, 4]]

        AUTHOR:

        - Rob Beezer (2011-05-04)
        """
        # Relies on vertices the fact that _elements correspond to the rows and
        # columns of the lequal matrix
        leq_mat = self.lequal_matrix()
        n = leq_mat.nrows()
        elements = self._elements
        for i in range(n):
            for j in range(i, n):
                if leq_mat[i,j]:
                    yield [elements[i], elements[j]]

    def is_lequal(self, x, y):
        """
        Returns True if x is less than or equal to y in the poset, and
        False otherwise.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = Q(0),Q(1),Q(4)
            sage: Q.is_lequal(x,y)
            False
            sage: Q.is_lequal(y,x)
            False
            sage: Q.is_lequal(x,z)
            True
            sage: Q.is_lequal(y,z)
            True
            sage: Q.is_lequal(z,z)
            True
        """
        r = self.compare_elements(x,y)
        return r == 0 or r == -1
    le = is_lequal

    def is_less_than(self, x, y):
        """
        Returns True if x is less than but not equal to y in the poset, and False
        otherwise.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = Q(0),Q(1),Q(4)
            sage: Q.is_less_than(x,y)
            False
            sage: Q.is_less_than(y,x)
            False
            sage: Q.is_less_than(x,z)
            True
            sage: Q.is_less_than(y,z)
            True
            sage: Q.is_less_than(z,z)
            False
        """
        return self.compare_elements(x,y) == -1
    lt = is_less_than

    def is_gequal(self, x, y):
        """
        Returns True if x is greater than or equal to y in the poset, and
        False otherwise.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = Q(0),Q(1),Q(4)
            sage: Q.is_gequal(x,y)
            False
            sage: Q.is_gequal(y,x)
            False
            sage: Q.is_gequal(x,z)
            False
            sage: Q.is_gequal(z,x)
            True
            sage: Q.is_gequal(z,y)
            True
            sage: Q.is_gequal(z,z)
            True
        """
        r = self.compare_elements(x,y)
        return r == 0 or r == 1
    ge = is_gequal

    def is_greater_than(self, x, y):
        """
        Returns True if ``x`` is greater than but not equal to ``y``
        in the poset, and False otherwise.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = Q(0),Q(1),Q(4)
            sage: Q.is_greater_than(x,y)
            False
            sage: Q.is_greater_than(y,x)
            False
            sage: Q.is_greater_than(x,z)
            False
            sage: Q.is_greater_than(z,x)
            True
            sage: Q.is_greater_than(z,y)
            True
            sage: Q.is_greater_than(z,z)
            False
        """
        return self.compare_elements(x,y) == 1
    gt = is_greater_than

    def compare_elements(self, x, y):
        r"""
        Compare ``x`` and ``y`` in the poset.

        If ``x`` = ``y``, then ``0`` is returned;
        if ``x`` < ``y``, then ``-1`` is returned;
        if ``x`` > ``y``, then ``1`` is returned;
        and if ``x`` and ``y`` are not comparable,
        then ``None`` is returned.

        EXAMPLES::

            sage: P = Poset([[1,2],[4],[3],[4],[]])
            sage: P.compare_elements(0,0)
            0
            sage: P.compare_elements(0,4)
            -1
            sage: P.compare_elements(4,0)
            1
            sage: P.compare_elements(1,2)

        """
        i, j = map(self._element_to_vertex,(x,y))
        if i == j:
            return 0
        elif self._hasse_diagram.is_less_than(i, j):
            return -1
        elif self._hasse_diagram.is_less_than(j, i):
            return  1
        else:
            return None

    def minimal_elements(self):
        """
        Returns a list of the minimal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P(0) in P.minimal_elements()
            True
            sage: P(1) in P.minimal_elements()
            True
            sage: P(2) in P.minimal_elements()
            True
        """
        return map(self._vertex_to_element, self._hasse_diagram.minimal_elements())

    def maximal_elements(self):
        """
        Returns a list of the maximal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.maximal_elements()
            [4]
        """
        return map(self._vertex_to_element, self._hasse_diagram.maximal_elements())

    def bottom(self):
        """
        Returns the bottom element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.bottom() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.bottom()
            0
        """
        hasse_bot = self._hasse_diagram.bottom()
        if hasse_bot is None:
            return None
        else:
            return self._vertex_to_element(hasse_bot)

    def has_bottom(self):
        """
        Returns True if the poset has a unique minimal element.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.has_bottom()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_bottom()
            True
        """
        return self._hasse_diagram.has_bottom()

    def top(self):
        """
        Returns the top element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.top() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.top()
            1

        TESTS::

            sage: R = Poset([[0],[]])
            sage: R.list()
            [0]
            sage: R.top() #Trac #10776
            0

        """
        hasse_top = self._hasse_diagram.top()
        if hasse_top is None:
            return None
        else:
            return self._vertex_to_element(hasse_top)

    def has_top(self):
        """
        Returns True if the poset contains a unique maximal element, and
        False otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.has_top()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_top()
            True
        """
        return self._hasse_diagram.has_top()

    def is_bounded(self):
        """
        Returns True if the poset contains a unique maximal element and a
        unique minimal element, and False otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.is_bounded()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.is_bounded()
            True
        """
        return self._hasse_diagram.is_bounded()

    def is_chain(self):
        """
        Returns True if the poset is totally ordered, and False otherwise.

        EXAMPLES::

            sage: L = Poset({0:[1],1:[2],2:[3],3:[4]})
            sage: L.is_chain()
            True
            sage: V = Poset({0:[1,2]})
            sage: V.is_chain()
            False
        """
        return self._hasse_diagram.is_chain()

    def rank_function(self):
        r"""
        Returns a rank function of the poset, if it exists.

        A *rank function* of a poset `P` is a function `r`
        from that maps elements of `P` to integers and satisfies:
        `r(x) = r(y) + 1` if `x` covers
        `y`.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: P.rank_function() is not None
            True
            sage: r = P.rank_function()
            sage: for u,v in P.cover_relations_iterator():
            ...    if r(v) != r(u) + 1:
            ...        print "Bug in rank_function!"

        ::

            sage: Q = Poset([[1,2],[4],[3],[4],[]])
            sage: Q.rank_function() is None
            True
        """
        hasse_rf = self._hasse_diagram.rank_function()
        if hasse_rf is None:
            return None
        else:
            return lambda z: hasse_rf(self._element_to_vertex(z))

    def rank(self,element=None):
        r"""
        Returns the rank of an element, or the rank of the poset if element
        is None. (The rank of a poset is the length of the longest chain of
        elements of the poset.)

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: P.rank(5)
            2
            sage: P.rank()
            3
            sage: Q = Poset([[1,2],[3],[],[]])

            sage: P = Posets.SymmetricGroupBruhatOrderPoset(4)

            sage: [(v,P.rank(v)) for v in P]
            [(1234, 0),
             (1324, 1),
            ...
             (4231, 5),
             (4321, 6)]
        """
        if element is None:
            return len(self.level_sets())-1
        elif self.is_ranked():
            return self.rank_function()(element)
        else:
            raise ValueError, "Poset is not ranked."

    def is_ranked(self):
        r"""
        Returns True if the poset is ranked, and False otherwise.

        A poset is {ranked} if it admits a rank function.

        EXAMPLES::

            sage: P = Poset([[1],[2],[3],[4],[]])
            sage: P.is_ranked()
            True
            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.is_ranked()
            False
        """
        return bool(self.rank_function())

    def is_graded(self):
        r"""
        Returns True if the poset is graded, and False otherwise.

        A poset is *graded* if it admits a rank function.

        EXAMPLES::

            sage: P = Poset([[1],[2],[3],[4],[]])
            sage: P.is_graded()
            True
            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.is_graded()
            False
        """
        return self.is_ranked()

    def covers(self,x,y):
        """
        Returns True if y covers x and False otherwise.

        EXAMPLES::

            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.covers(Q(1),Q(6))
            True
            sage: Q.covers(Q(1),Q(4))
            False
        """
        return self._hasse_diagram.has_edge(*map(self._element_to_vertex,(x,y)))

    def upper_covers_iterator(self,y):
        """
        Returns an iterator for the upper covers of the element y. An upper
        cover of y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.upper_covers_iterator(0))
            <type 'generator'>
        """
        for x in self._hasse_diagram.neighbor_out_iterator(self._element_to_vertex(y)):
            yield self._vertex_to_element(x)

    def upper_covers(self,y):
        """
        Returns a list of upper covers of the element y. An upper cover of
        y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: map(Q.upper_covers,Q.list())
            [[2], [2], [3], [4], []]
        """
        return [x for x in self.upper_covers_iterator(y)]

    def lower_covers_iterator(self,y):
        """
        Returns an iterator for the lower covers of the element y. An lower
        cover of y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.lower_covers_iterator(0))
            <type 'generator'>
        """
        for x in self._hasse_diagram.neighbor_in_iterator(self._element_to_vertex(y)):
            yield self._vertex_to_element(x)

    def lower_covers(self,y):
        """
        Returns a list of lower covers of the element y. An lower cover of
        y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: map(Q.lower_covers,Q.list())
            [[], [], [1, 0], [2], [3]]
        """
        return [x for x in self.lower_covers_iterator(y)]

    def cardinality(self):
        """
        Returns the number of elements in the poset.

        EXAMPLES::

            sage: Poset([[1,2,3],[4],[4],[4],[]]).cardinality()
            5
        """
        return Integer(self._hasse_diagram.order())

    size = deprecated_function_alias(cardinality, 'Sage Version 4.4 (2010-05)')

    def mobius_function(self,x,y):
        r"""
        Returns the value of the Mobius function of the poset on the
        elements x and y.

        EXAMPLES::

            sage: P = Poset([[1,2,3],[4],[4],[4],[]])
            sage: P.mobius_function(P(0),P(4))
            2
            sage: sum([P.mobius_function(P(0),v) for v in P])
            0
            sage: sum([abs(P.mobius_function(P(0),v)) \
            ...        for v in P])
            6
            sage: for u,v in P.cover_relations_iterator():
            ...    if P.mobius_function(u,v) != -1:
            ...        print "Bug in mobius_function!"

        ::

            sage: Q = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: Q.mobius_function(Q(0),Q(7))
            0
            sage: Q.mobius_function(Q(0),Q(5))
            0
            sage: Q.mobius_function(Q(2),Q(7))
            2
            sage: Q.mobius_function(Q(3),Q(3))
            1
            sage: sum([Q.mobius_function(Q(0),v) for v in Q])
            0
        """
        i,j = map(self._element_to_vertex,(x,y))
        return self._hasse_diagram.mobius_function(i,j)

    def mobius_function_matrix(self, ring = ZZ, sparse = False):
        r"""
        Returns a matrix whose ``(i,j)`` entry is the value of the Mobius
        function evaluated at ``self.linear_extension()[i]`` and
        ``self.linear_extension()[j]``.

        INPUT:

        - ``ring`` -- the ring of coefficients (default: ``ZZ``)
        - ``sparse`` -- whether the returned matrix is sparse or not
          (default: ``True``)

        EXAMPLES::

            sage: P = Poset([[4,2,3],[],[1],[1],[1]])
            sage: x,y = (P.linear_extension()[0],P.linear_extension()[1])
            sage: P.mobius_function(x,y)
            -1
            sage: M = P.mobius_function_matrix(); M
            [ 1 -1 -1 -1  2]
            [ 0  1  0  0 -1]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
            [ 0  0  0  0  1]
            sage: M[0,4]
            2
            sage: M[0,1]
            -1

        We now demonstrate the usage of the optional parameters::

            sage: P.mobius_function_matrix(ring=QQ, sparse=False).parent()
            Full MatrixSpace of 5 by 5 dense matrices over Rational Field
        """
        M = self._hasse_diagram.mobius_function_matrix()
        if ring is not ZZ:
            M = M.change_ring(ring)
        if not sparse:
            M = M.dense_matrix()
        return M


    def lequal_matrix(self, ring = ZZ, sparse = False):
        """
        Computes the matrix whose ``(i,j)`` entry is 1 if
        ``self.linear_extension()[i] < self.linear_extension()[j]`` and 0
        otherwise.

        INPUT:

        - ``ring`` -- the ring of coefficients (default: ``ZZ``)
        - ``sparse`` -- whether the returned matrix is sparse or not
          (default: ``True``)

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: LEQM = P.lequal_matrix(); LEQM
            [1 1 1 1 1 1 1 1]
            [0 1 0 1 0 0 0 1]
            [0 0 1 1 1 0 1 1]
            [0 0 0 1 0 0 0 1]
            [0 0 0 0 1 0 0 1]
            [0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 1 1]
            [0 0 0 0 0 0 0 1]
            sage: LEQM[1,3]
            1
            sage: P.linear_extension()[1] < P.linear_extension()[3]
            True
            sage: LEQM[2,5]
            0
            sage: P.linear_extension()[2] < P.linear_extension()[5]
            False

        We now demonstrate the usage of the optional parameters::

            sage: P.lequal_matrix(ring=QQ, sparse=False).parent()
            Full MatrixSpace of 8 by 8 dense matrices over Rational Field
        """
        M = self._hasse_diagram.lequal_matrix()
        if ring is not ZZ:
            M = M.change_ring(ring)
        if not sparse:
            M = M.dense_matrix()
        return M

    def coxeter_transformation(self):
        r"""
        Returns the matrix of the Auslander-Reiten translation acting
        on the Grothendieck group of the derived category of modules
        on the poset ``self``, in the basis of simple modules. This matrix is
        usually called the Coxeter transformation.

        EXAMPLES::

            sage: Posets.PentagonPoset().coxeter_transformation()
            [ 0  0  0  0 -1]
            [ 0  0  0  1 -1]
            [ 0  1  0  0 -1]
            [-1  1  1  0 -1]
            [-1  1  0  1 -1]

        TESTS::

            sage: M = Posets.PentagonPoset().coxeter_transformation()
            sage: M**8 == 1
            True
        """
        return self._hasse_diagram.coxeter_transformation()

    def meet_matrix(self):
        """
        Returns a matrix whose ``(i,j)`` entry is ``k``, where
        ``self.linear_extension()[k]`` is the meet (greatest lower bound) of
        ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: M = P.meet_matrix(); M
            [0 0 0 0 0 0 0 0]
            [0 1 0 1 0 0 0 1]
            [0 0 2 2 2 0 2 2]
            [0 1 2 3 2 0 2 3]
            [0 0 2 2 4 0 2 4]
            [0 0 0 0 0 5 5 5]
            [0 0 2 2 2 5 6 6]
            [0 1 2 3 4 5 6 7]
            sage: M[P(4).vertex,P(3).vertex] == P(0).vertex
            True
            sage: M[P(5).vertex,P(2).vertex] == P(2).vertex
            True
            sage: M[P(5).vertex,P(2).vertex] == P(5).vertex
            False
        """
        return self._hasse_diagram.meet_matrix()

    def is_meet_semilattice(self):
        r"""
        Returns True if self has a meet operation, and False otherwise.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: P.is_meet_semilattice()
            True

            sage: P = Poset([[1,2],[3],[3],[]])
            sage: P.is_meet_semilattice()
            True

            sage: P = Poset({0:[2,3],1:[2,3]})
            sage: P.is_meet_semilattice()
            False
        """
        return self._hasse_diagram.is_meet_semilattice()

    def join_matrix(self):
        """
        Returns a matrix whose ``(i,j)`` entry is ``k``, where
        ``self.linear_extension()[k]`` is the join (least upper bound) of
        ``self.linear_extension()[i]`` and ``self.linear_extension()[j]``.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: J = P.join_matrix(); J
            [0 1 2 3 4 5 6 7]
            [1 1 3 3 7 7 7 7]
            [2 3 2 3 4 6 6 7]
            [3 3 3 3 7 7 7 7]
            [4 7 4 7 4 7 7 7]
            [5 7 6 7 7 5 6 7]
            [6 7 6 7 7 6 6 7]
            [7 7 7 7 7 7 7 7]
            sage: J[P(4).vertex,P(3).vertex] == P(7).vertex
            True
            sage: J[P(5).vertex,P(2).vertex] == P(5).vertex
            True
            sage: J[P(5).vertex,P(2).vertex] == P(2).vertex
            False
        """
        return self._hasse_diagram.join_matrix()

    def is_join_semilattice(self):
        """
        Returns True is the poset has a join operation, and False
        otherwise.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: P.is_join_semilattice()
            True

            sage: P = Poset([[1,2],[3],[3],[]])
            sage: P.is_join_semilattice()
            True

            sage: P = Poset({0:[2,3],1:[2,3]})
            sage: P.is_join_semilattice()
            False
        """
        return self._hasse_diagram.is_join_semilattice()

    def is_isomorphic(self,other):
        """
        Returns True if both posets are isomorphic.

        EXAMPLES::

            sage: P = Poset(([1,2,3],[[1,3],[2,3]]))
            sage: Q = Poset(([4,5,6],[[4,6],[5,6]]))
            sage: P.is_isomorphic( Q )
            True
        """
        if hasattr(other,'hasse_diagram'):
            return self.hasse_diagram().is_isomorphic( other.hasse_diagram() )
        else:
            raise ValueError, 'The input is not a finite poset.'

    import __builtin__ # Caveat: list is overridden by the method list above!!!
    def antichains(self, element_constructor = __builtin__.list):
        """
        Returns the antichains of the poset.

        INPUT:

         - ``element_constructor`` -- a function taking an iterable as
           argument (default: list)

        OUTPUT: an enumerated set

        An *antichain* of a poset is a collection of elements of the
        poset that are pairwise incomparable.

        EXAMPLES::

            sage: A = Posets.PentagonPoset().antichains(); A
            Set of antichains of Finite lattice containing 5 elements
            sage: list(A)
            [[], [0], [1], [1, 2], [1, 3], [2], [3], [4]]
            sage: A.cardinality()
            8
            sage: A[3]
            [1, 2]
            sage: list(Posets.AntichainPoset(3).antichains())
            [[], [2], [2, 1], [2, 1, 0], [2, 0], [1], [1, 0], [0]]
            sage: list(Posets.ChainPoset(3).antichains())
            [[], [0], [1], [2]]

        To get the antichains of a given size one can currently use::

            sage: list(A.elements_of_depth_iterator(2))
            [[1, 2], [1, 3]]

        Eventually the following syntax will be accepted::

            sage: A.subset(size = 2) # todo: not implemented

        To get the antichains as, say, sets, one may use the
        ``element_constructor`` option::

            sage: list(Posets.ChainPoset(3).antichains(element_constructor = set))
            [set([]), set([0]), set([1]), set([2])]

        .. note:: this function used to return a list; this change is
            slightly backward incompatible; e.g. ``len(A)`` does not work.

        .. note:: Internally, this uses
            :class:`sage.combinat.subsets_pairwise.PairwiseCompatibleSubsets`
            and :class:`SearchForest`. At this point, iterating
            through this set is about twice slower than using
            :meth:`antichains_iterator` (tested on
            ``posets.AntichainPoset(15)``). The algorithm is the same
            (depth first search through the tree), but
            :meth:`antichains_iterator` manually inlines things which
            apparently avoids some infrastructure overhead.

            On the other hand, this returns a full featured enumerated
            set, with containment testing, etc.

        """
        vertex_to_element = self._vertex_to_element
        def f(antichain):
            return element_constructor(vertex_to_element(x) for x in antichain)
        result = self._hasse_diagram.antichains(element_class = f)
        result.rename("Set of antichains of %s"%self)
        return result

    def antichains_iterator(self):
        """
        Returns an iterator over the antichains of the poset.

        EXAMPLES::

            sage: Posets.PentagonPoset().antichains_iterator()
            <generator object antichains_iterator at ...>

        .. seealso:: :meth:`antichains`
        """
        vertex_to_element = self._vertex_to_element
        for antichain in self._hasse_diagram.antichains_iterator():
            yield map(vertex_to_element, antichain)

    def chains(self, element_constructor = __builtin__.list):
        """
        Returns all the chains of ``self``

        INPUT:

         - ``element_constructor`` -- a function taking an iterable as
           argument (default: list)

        OUTPUT: an enumerated set

        A *chain* of a poset is a collection of elements of the poset
        that are pairwise comparable.

        EXAMPLES::

            sage: A = Posets.PentagonPoset().chains(); A
            Set of chains of Finite lattice containing 5 elements
            sage: list(A)
            [[], [0], [0, 1], [0, 1, 4], [0, 2], [0, 2, 3], [0, 2, 3, 4], [0, 2, 4], [0, 3], [0, 3, 4], [0, 4], [1], [1, 4], [2], [2, 3], [2, 3, 4], [2, 4], [3], [3, 4], [4]]

        To get the chains of a given size one can currently use::

            sage: list(A.elements_of_depth_iterator(2))
            [[0, 1], [0, 2], [0, 3], [0, 4], [1, 4], [2, 3], [2, 4], [3, 4]]

        Eventually the following syntax will be accepted::

            sage: A.subset(size = 2) # todo: not implemented

        .. seealso:: :meth:`maximal_chains`, :meth:`antichains`
        """
        vertex_to_element = self._vertex_to_element
        def f(chain):
            return element_constructor(vertex_to_element(x) for x in chain)
        result = self._hasse_diagram.chains(element_class = f)
        result.rename("Set of chains of %s"%self)
        return result

    def dual(self):
        """
        Returns the dual poset of the given poset.

        EXAMPLE::

            sage: P = Poset(([1,2,3],[[1,2],[1,3]]))
            sage: P.cover_relations()
            [[1, 2], [1, 3]]
            sage: Q = P.dual()
            sage: Q.cover_relations()
            [[3, 1], [2, 1]]

            sage: P = LatticePoset([[1,2],[3],[3]], facade = True)
            sage: P.cover_relations()
            [[0, 1], [0, 2], [1, 3], [2, 3]]
            sage: Q = P.dual()
            sage: Q.cover_relations()
            [[3, 2], [3, 1], [2, 0], [1, 0]]
            sage: Q.category()
            Category of facade finite lattice posets
            sage: Q.__class__
            <class 'sage.combinat.posets.lattices.FiniteLatticePoset_with_category'>

            sage: P = MeetSemilattice([[1,2],[3],[3]])
            sage: P.dual().__class__
            <class 'sage.combinat.posets.lattices.FiniteJoinSemilattice_with_category'>
            sage: P = JoinSemilattice([[1,2],[3],[3]])
            sage: P.dual().__class__
            <class 'sage.combinat.posets.lattices.FiniteMeetSemilattice_with_category'>
        """
        n = self.cardinality()
        dual_hasse_digraph = DiGraph(self._hasse_diagram).reverse()
        dual_hasse_digraph.relabel(lambda i: n-i-1)

        return self._dual_class(dual_hasse_digraph,
                                elements = reversed(self._elements),
                                category = self.category(),
                                facade = self._is_facade)

    def relabel(self, relabelling):
        r"""
        Returns a copy of this poset with its elements relabelled

        INPUT:

        - ``relabelling`` -- a function or dictionnary

          This function should map each (non-wrapped) element of
          ``self`` to some distinct object.

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: P.cover_relations()
            [[1, 2], [1, 3], [2, 4], [2, 6], [3, 6], [4, 12], [6, 12]]
            sage: Q = P.relabel(lambda x: 12/x)
            sage: Q.list()
            [12, 6, 4, 3, 2, 1]
            sage: Q.cover_relations()
            [[12, 6], [12, 4], [6, 3], [6, 2], [4, 2], [3, 1], [2, 1]]

        Here we relabel the elements of a poset by {0,1,2, ...}, using
        a dictionary::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
            sage: relabelling = {c.element:i for (i,c) in enumerate(P)}; relabelling
            {1: 0, 2: 1, 3: 2, 4: 3, 6: 4, 12: 5}
            sage: Q = P.relabel(relabelling)
            sage: Q.list()
            [0, 1, 2, 3, 4, 5]
            sage: Q.cover_relations()
            [[0, 1], [0, 2], [1, 3], [1, 4], [2, 4], [3, 5], [4, 5]]

        Mind the ``c.element``; this is because the relabelling is
        applied to the elements of the poset without the wrapping.
        Thanks to this convention, the same relabelling function can
        be used both for facade or non facade posets::

            sage: P = Poset((divisors(12), attrcall("divides")), facade = True, linear_extension=True)
            sage: P.list()
            [1, 2, 3, 4, 6, 12]
            sage: Q = P.relabel(lambda x: 12/x)
            sage: Q.list()
            [12, 6, 4, 3, 2, 1]
            sage: Q.cover_relations()
            [[12, 6], [12, 4], [6, 3], [6, 2], [4, 2], [3, 1], [2, 1]]

        .. note::

            As can be seen in the above examples, the default linear
            extension of ``Q`` is that of ``P`` after relabelling. In
            particular, ``P`` and ``Q`` share the same internal Hasse
            diagram.
        """
        assert not isinstance(relabelling, (tuple, list)), "relabelling by tuple or list not yet defined"
        if isinstance(relabelling, dict):
            relabelling = relabelling.__getitem__
        elements = tuple(relabelling(x) for x in self._elements)
        return FinitePoset(self._hasse_diagram,
                           elements = elements,
                           category=self.category(),
                           facade=self._is_facade)

    def with_linear_extension(self, linear_extension):
        """
        Returns a copy of ``self`` with a different default linear extension

        EXAMPLES::

            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
            sage: P.cover_relations()
            [[1, 2], [1, 3], [2, 4], [2, 6], [3, 6], [4, 12], [6, 12]]
            sage: list(P)
            [1, 2, 3, 4, 6, 12]
            sage: Q = P.with_linear_extension([1,3,6,2,4,12])
            sage: list(Q)
            [1, 3, 6, 2, 4, 12]
            sage: Q.cover_relations()
            [[1, 3], [1, 2], [3, 6], [6, 12], [2, 6], [2, 4], [4, 12]]

        TESTS:

        We check that we can pass in a list of elements of P instead::

            sage: Q = P.with_linear_extension(map(P, [1,3,6,2,4,12]))
            sage: list(Q)
            [1, 3, 6, 2, 4, 12]
            sage: Q.cover_relations()
            [[1, 3], [1, 2], [3, 6], [6, 12], [2, 6], [2, 4], [4, 12]]

        We check that this works for facade posets too::

            sage: P = Poset((divisors(12), attrcall("divides")), facade=True)
            sage: Q = P.with_linear_extension([1,3,6,2,4,12])
            sage: list(Q)
            [1, 3, 6, 2, 4, 12]
            sage: Q.cover_relations()
            [[1, 3], [1, 2], [3, 6], [6, 12], [2, 6], [2, 4], [4, 12]]
            sage: sorted(Q.cover_relations()) == sorted(P.cover_relations())
            True

        .. note::

            With the current implementation, this requires relabelling
            the internal Dynkin diagram which is `O(n+m)`, where `n`
            is the number of elements and `m` the number of cover
            relations.

        """
        new_vertices = [ self._element_to_vertex(element) for element in linear_extension]
        new_elements = [ self._elements[i] for i in new_vertices ]
        vertex_relabelling = dict(zip(new_vertices, range(len(new_vertices))))
        return FinitePoset(self._hasse_diagram.relabel(vertex_relabelling, inplace=False),
                           elements = new_elements,
                           category=self.category(),
                           facade=self._is_facade)

    def graphviz_string(self,graph_string="graph",edge_string="--"):
        r"""
        Returns a representation in the DOT language, ready to render in
        graphviz.

        REFERENCES:

        - http://www.graphviz.org/doc/info/lang.html

        EXAMPLES::

            sage: P = Poset({'a':['b'],'b':['d'],'c':['d'],'d':['f'],'e':['f'],'f':[]})
            sage: print P.graphviz_string()
            graph {
            "f";"d";"b";"a";"c";"e";
            "f"--"e";"d"--"c";"b"--"a";"d"--"b";"f"--"d";
            }
        """
        s = '%s {\n' % graph_string
        for v in reversed(self.list()):
            s+= '"%s";'%v
        s+= '\n'
        for u, v in self.cover_relations_iterator():
            s+= '"%s"%s"%s";' % (v, edge_string, u)
        s+= "\n}"
        return s

    def subposet(self, elements):
        """
        Returns the poset containing elements with partial order induced by
        that of self.

        EXAMPLES::

            sage: P = Poset({"a":["c","d"], "b":["d","e"], "c":["f"], "d":["f"], "e":["f"]})
            sage: Q = P.subposet(["a","b","f"]); Q
            Finite poset containing 3 elements
            sage: Q.cover_relations()
            [[b, f], [a, f]]

        A subposet of a facade poset is again a facade poset::

            sage: P = Poset({"a":["c","d"], "b":["d","e"], "c":["f"], "d":["f"], "e":["f"]}, facade=True)
            sage: Q = P.subposet(["a","b","f"]); Q
            Finite poset containing 3 elements
            sage: Q.cover_relations()
            [['b', 'f'], ['a', 'f']]

        One may specified wrapped elements or not::

            sage: P = Poset({"a":["c","d"], "b":["d","e"], "c":["f"], "d":["f"], "e":["f"]})
            sage: Q = P.subposet([P("a"),P("b"),P("f")]); Q
            Finite poset containing 3 elements
            sage: Q.cover_relations()
            [[b, f], [a, f]]

            sage: B = posets.BooleanLattice(2)
            sage: above = B.principal_order_filter(0)
            sage: Q = B.subposet(above)
            sage: above_new = Q.principal_order_filter(Q.list()[0])
            sage: Q.subposet(above_new)
            Finite poset containing 4 elements

        TESTS::

            sage: P.subposet(("a","b","f"))
            Finite poset containing 3 elements
            sage: P.subposet(["a","b","x"])
            Traceback (most recent call last):
            ...
            ValueError: <type 'str'> is not an element of this poset
            sage: P.subposet(3)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
        """
        # Type checking is performed by the following line:
        elements  = [self(e) for e in elements]
        relations = []
        for u in elements:
            for v in elements:
                if self.is_less_than(u,v):
                    relations.append([u,v])
        if not self._is_facade:
            elements = [e.element for e in elements]
            relations = [[u.element,v.element] for u,v in relations]
        return Poset((elements, relations), cover_relations=False, facade=self._is_facade)

    def random_subposet(self, p):
        """
        Returns a random subposet that contains each element with
        probability p.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: Q = P.random_subposet(.25)
        """
        elements = []
        p = float(p)
        for v in self:
            if random.random() <= p:
                elements.append(v)
        return self.subposet(elements)

    def order_filter(self,elements):
        """
        Returns the order filter generated by a list of elements.

        `I` is an order filter if, for any `x` in `I` and `y` such that
        `y \ge x`, then `y` is in `I`.

        EXAMPLES::

            sage: B = Posets.BooleanLattice(4)
            sage: B.order_filter([3,8])
            [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        vertices = sorted(map(self._element_to_vertex,elements))
        of = self._hasse_diagram.order_filter(vertices)
        return map(self._vertex_to_element,of)

    def order_ideal(self,elements):
        """
        Returns the order ideal generated by a list of elements.

        `I` is an order ideal if, for any `x` in `I` and `y` such that
        `y \le x`, then `y` is in `I`.

        EXAMPLES::

            sage: B = Posets.BooleanLattice(4)
            sage: B.order_ideal([7,10])
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
        """
        vertices = map(self._element_to_vertex,elements)
        oi = self._hasse_diagram.order_ideal(vertices)
        return map(self._vertex_to_element,oi)

    def interval(self, x, y):
        """
        Returns a list of the elements `z` such that `x \le z \le y`.
        The order is that induced by the ordering in
        ``self.linear_extension()``.

        INPUT:


        -  ``x`` - any element of the poset

        -  ``y`` - any element of the poset


        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: P = Poset(dag)
            sage: I = set(map(P,[2,5,6,4,7]))
            sage: I == set(P.interval(2,7))
            True

        ::

            sage: dg = DiGraph({"a":["b","c"], "b":["d"], "c":["d"]})
            sage: P = Poset(dg)
            sage: P.interval("a","d")
            [a, c, b, d]
        """
        return map(self._vertex_to_element,self._hasse_diagram.interval(
                self._element_to_vertex(x),self._element_to_vertex(y)))

    def closed_interval(self, x, y):
        """
        Returns a list of the elements `z` such that `x \le z \le y`.
        The order is that induced by the ordering in
        ``self.linear_extension()``.

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: P = Poset(dag)
            sage: I = set(map(P,[2,5,6,4,7]))
            sage: I == set(P.closed_interval(2,7))
            True
        """
        return self.interval(x,y)

    def open_interval(self, x, y):
        """
        Returns a list of the elements `z` such that `x < z < y`. The
        order is that induced by the ordering in
        ``self.linear_extension()``.

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: P = Poset(dag)
            sage: I = set(map(P,[5,6,4]))
            sage: I == set(P.open_interval(2,7))
            True

        ::

            sage: dg = DiGraph({"a":["b","c"], "b":["d"], "c":["d"]})
            sage: P = Poset(dg)
            sage: P.open_interval("a","d")
            [c, b]
        """
        return map(self._vertex_to_element,self._hasse_diagram.open_interval(
                self._element_to_vertex(x),self._element_to_vertex(y)))

    def maximal_chains(self, partial=None):
        """
        Returns all maximal chains of this poset.  Each chain
        is listed in increasing order.

        INPUT:


        -  ``partial`` - list (optional).  If present, find all maximal
           chains starting with the elements in partial.

        Returns list of the maximal chains of this poset.

        This is used in constructing the order complex for the poset.

        EXAMPLES::

            sage: P = Posets.BooleanLattice(3)
            sage: P.maximal_chains()
            [[0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7], [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 6, 7]]
            sage: P.maximal_chains(partial=[0,2])
            [[0, 2, 3, 7], [0, 2, 6, 7]]
            sage: Q = Posets.ChainPoset(6)
            sage: Q.maximal_chains()
            [[0, 1, 2, 3, 4, 5]]
        """
        if partial is None or len(partial) == 0:
            start = self.minimal_elements()
            partial = []
        else:
            start = self.upper_covers(partial[-1])
        if len(start) == 0:
            return [partial]
        if len(start) == 1:
            return self.maximal_chains(partial=partial + start)
        parts = [partial + [x] for x in start]
        answer = []
        for new in parts:
            answer += self.maximal_chains(partial=new)
        return answer

    def order_complex(self, on_ints=False):
        """
        Returns the order complex associated to this poset.

        The order complex is the simplicial complex with vertices equal
        to the elements of the poset, and faces given by the chains.

        INPUT:

        - ``on_ints`` -- a boolean (default: False)

        EXAMPLES::

            sage: P = Posets.BooleanLattice(3)
            sage: S = P.order_complex(); S
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6, 7) and 6 facets
            sage: S.f_vector()
            [1, 8, 19, 18, 6]
            sage: S.homology()      # S is contractible
            {0: 0, 1: 0, 2: 0, 3: 0}
            sage: Q = P.subposet([1,2,3,4,5,6])
            sage: Q.order_complex().homology()    # a circle
            {0: 0, 1: Z}

            sage: P = Poset((divisors(15), attrcall("divides")), facade = True)
            sage: P.order_complex()
            Simplicial complex with vertex set (1, 3, 5, 15) and facets {(1, 3, 15), (1, 5, 15)}

        If ``on_ints``, then the elements of the poset are labelled
        `0,1,\dots` in the chain complex::

            sage: P.order_complex(on_ints=True)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 3)}
        """
        from sage.homology.simplicial_complex import SimplicialComplex
        L = self.list()
        if on_ints:
            vertices = range(len(L))
            iso = dict( [ (L[i],i) for i in vertices ] )
        else:
            vertices = self._elements

        facets = []
        for f in self.maximal_chains():
            # TODO: factor out the logic for on_ints / facade / ...
            # We will want to do similar things elsewhere
            if on_ints:
                facets.append([iso[a]    for a in f])
            elif self._is_facade:
                facets.append([a         for a in f])
            else:
                facets.append([a.element for a in f])

        return SimplicialComplex(vertices, facets)

    def promotion(self, i=1):
        r"""
        Computes the (extended) promotion on the linear extension of the poset ``self``

        INPUT:

        - ``i`` -- an integer between `1` and `n` (default: `1`)

        OUTPUT:

        - an isomorphic poset, with the same default linear extension

        The extended promotion is defined on a poset ``self`` of size
        `n` by applying the promotion operator `\tau_i \tau_{i+1}
        \cdots \tau_{n-1}` to the default linear extension `\pi` of ``self``
        (see :meth:`~sage.combinat.posets.linear_extensions.LinearExtensionOfPoset.promotion`),
        and relabelling ``self`` accordingly. For more details see [St2009]_.

        When the vertices of the poset ``self`` are labelled by
        `\{1,2,\ldots,n\}`, the linear extension is the identity, and
        `i=1`, the above algorithm corresponds to the promotion
        operator on posets defined by Schützenberger as
        follows. Remove `1` from ``self`` and replace it by the
        minimum `j` of all labels covering `1` in the poset. Then,
        remove `j` and replace it by the minimum of all labels
        covering `j`, and so on.  This process ends when a label is a
        local maximum. Place the label `n+1` at this vertex.  Finally,
        decrease all labels by `1`.

        REFERENCES:

            .. [St2009] Richard Stanley,
               *Promotion and evacuation*,
               Electron. J. Combin. 16 (2009), no. 2, Special volume in honor of Anders Björner,
               Research Paper 9, 24 pp.

        EXAMPLES::

            sage: P = Poset(([1,2], [[1,2]]))
            sage: P.promotion()
            Finite poset containing 2 elements
            sage: P == P.promotion()
            True

            sage: P = Poset(([1,2,3,4,5,6,7], [[1,2],[1,4],[2,3],[2,5],[3,6],[4,7],[5,6]]))
            sage: P.list()
            [1, 2, 3, 5, 6, 4, 7]
            sage: Q = P.promotion(4); Q
            Finite poset containing 7 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 6], [2, 3], [2, 5], [3, 7], [5, 7], [6, 4]]

        Note that if one wants to obtain the promotion defined by
        Schützenberger's algorithm directly on the poset, one needs
        to make sure the linear extension is the identity::

            sage: P = P.with_linear_extension([1,2,3,4,5,6,7])
            sage: P.list()
            [1, 2, 3, 4, 5, 6, 7]
            sage: Q = P.promotion(4); Q
            Finite poset containing 7 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 6], [2, 3], [2, 4], [3, 5], [4, 5], [6, 7]]
            sage: Q = P.promotion()
            sage: Q.cover_relations()
            [[1, 2], [1, 3], [2, 4], [2, 5], [3, 6], [4, 7], [5, 7]]

        Here is an example for a poset not labelled by `\{1,2,\ldots,n\}`::

            sage: P = Poset((divisors(30), attrcall("divides")), linear_extension = True)
            sage: P.list()
            [1, 2, 3, 5, 6, 10, 15, 30]
            sage: P.cover_relations()
            [[1, 2], [1, 3], [1, 5], [2, 6], [2, 10], [3, 6], [3, 15], [5, 10], [5, 15], [6, 30], [10, 30], [15, 30]]
            sage: Q = P.promotion(4); Q
            Finite poset containing 8 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 3], [1, 6], [2, 5], [2, 15], [3, 5], [3, 10], [5, 30], [6, 10], [6, 15], [10, 30], [15, 30]]

        .. seealso::

            - :meth:`linear_extension`
            - :meth:`with_linear_extension` and the ``linear_extension`` option of :func:`Poset`
            - :meth:`~sage.combinat.posets.linear_extensions.LinearExtensionOfPoset.promotion`
            - :meth:`evacuation`

        AUTHOR:

        - Anne Schilling (2012-02-18)
        """
        return self.linear_extension(self.linear_extension()).promotion(i).to_poset()

    def evacuation(self):
        r"""
        Computes evacuation on the linear extension associated to the poset ``self``.

        OUTPUT:

        - an isomorphic poset, with the same default linear extension

        Evacuation is defined on a poset ``self`` of size `n` by
        applying the evacuation operator
        `(\tau_1 \cdots \tau_{n-1}) (\tau_1 \cdots \tau_{n-2}) \cdots (\tau_1)`,
        to the default linear extension `\pi` of ``self``
        (see :meth:`~sage.combinat.posets.linear_extensions.LinearExtensionOfPoset.evacuation`),
        and relabelling ``self`` accordingly. For more details see [Stan2009]_.

        .. seealso::

            - :meth:`linear_extension`
            - :meth:`with_linear_extension` and the ``linear_extension`` option of :func:`Poset`
            - :meth:`~sage.combinat.posets.linear_extensions.LinearExtensionOfPoset.evacuation`
            - :meth:`promotion`

        REFERENCES:

            .. [Stan2009] Richard Stanley,
               *Promotion and evacuation*,
               Electron. J. Combin. 16 (2009), no. 2, Special volume in honor of Anders Björner,
               Research Paper 9, 24 pp.

        EXAMPLES::

            sage: P = Poset(([1,2], [[1,2]]))
            sage: P.evacuation()
            Finite poset containing 2 elements
            sage: P.evacuation() == P
            True

            sage: P = Poset(([1,2,3,4,5,6,7], [[1,2],[1,4],[2,3],[2,5],[3,6],[4,7],[5,6]]), linear_extension = True)
            sage: P.list()
            [1, 2, 3, 4, 5, 6, 7]
            sage: Q = P.evacuation(); Q
            Finite poset containing 7 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 3], [2, 5], [3, 4], [3, 6], [4, 7], [6, 7]]

        Note that the results depend on the linear extension associated to the poset::

            sage: P = Poset(([1,2,3,4,5,6,7], [[1,2],[1,4],[2,3],[2,5],[3,6],[4,7],[5,6]]))
            sage: P.list()
            [1, 2, 3, 5, 6, 4, 7]
            sage: Q = P.evacuation(); Q
            Finite poset containing 7 elements
            sage: Q.cover_relations()
            [[1, 2], [1, 5], [2, 3], [5, 6], [5, 4], [6, 7], [4, 7]]

        Here is an example of a poset where the vertices are not labelled by `\{1,2,\ldots,n\}`::

            sage: P = Poset((divisors(15), attrcall("divides")), linear_extension = True)
            sage: P.list()
            [1, 3, 5, 15]
            sage: Q = P.evacuation(); Q
            Finite poset containing 4 elements
            sage: Q.cover_relations()
            [[1, 3], [1, 5], [3, 15], [5, 15]]

        AUTHOR:

        - Anne Schilling (2012-02-18)
        """
        return self.linear_extension(self.linear_extension()).evacuation().to_poset()


FinitePoset._dual_class = FinitePoset

##### Posets #####

class FinitePosets_n(UniqueRepresentation, Parent):
    r"""
    The finite enumerated set of all posets on `n` vertices, up to an isomorphism.

    EXAMPLES::

        sage: P = Posets(3)
        sage: P.cardinality()
        5
        sage: for p in P: print p.cover_relations()
        []
        [[1, 2]]
        [[0, 1], [0, 2]]
        [[0, 1], [1, 2]]
        [[0, 2], [1, 2]]
    """

    def __init__(self, n):
        r"""
        EXAMPLES::

            sage: P = Posets(3); P
            Posets containing 3 vertices
            sage: P.category()
            Category of finite enumerated sets
            sage: P.__class__
            <class 'sage.combinat.posets.posets.FinitePosets_n_with_category'>
            sage: TestSuite(P).run()
        """
        Parent.__init__(self, category = FiniteEnumeratedSets())
        self._n = n

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: P = Posets(3)
            sage: P._repr_()
            'Posets containing 3 vertices'
        """
        return "Posets containing %s vertices" % self._n

    def __contains__(self, P):
        """
        EXAMPLES::

            sage: posets.PentagonPoset() in Posets(5)
            True
            sage: posets.PentagonPoset() in Posets(3)
            False
            sage: 1 in Posets(3)
            False
        """
        return P in FinitePosets() and P.cardinality() == self._n

    def __iter__(self):
        """
        Returns an iterator of representatives of the isomorphism classes
        of finite posets of a given size.

        .. note::

           This uses the DiGraph iterator as a backend to construct
           transitively-reduced, acyclic digraphs.

        EXAMPLES::
            sage: P = Posets(2)
            sage: list(P)
            [Finite poset containing 2 elements, Finite poset containing 2 elements]
        """
        from sage.graphs.digraph_generators import DiGraphGenerators
        for dig in DiGraphGenerators()(self._n, is_poset):
            # We need to relabel the digraph since range(self._n) must be a linear
            # extension. Too bad we need to compute this again. TODO: Fix this.
            label_dict = dict(zip(dig.topological_sort(),range(dig.order())))
            yield FinitePoset(dig.relabel(label_dict,inplace=False))

    def cardinality(self, from_iterator=False):
        r"""
        Return the cardinality of this object.

        .. note::

            By default, this returns pre-computed values obtained from
            the On-Line Encyclopedia of Integer Sequences (A000112).
            To override this, pass the argument ``from_iterator=True``.

        EXAMPLES::

            sage: P = Posets(3)
            sage: P.cardinality()
            5
            sage: P.cardinality(from_iterator=True)
            5
        """
        # Obtained from The On-Line Encyclopedia of Integer Sequences;
        # this is sequence number A000112.
        known_values = [1, 1, 2, 5, 16, 63, 318, 2045, 16999, 183231,
                2567284, 46749427, 1104891746, 33823827452, 1338193159771,
                68275077901156, 4483130665195087]
        if not from_iterator and self._n < len(known_values):
            return Integer(known_values[self._n])
        else:
            return super(FinitePosets_n, self).cardinality()

# For backward compatibility of pickles of the former Posets()
Posets_all = Posets

##### Miscellaneous functions #####

def is_poset(dig):
    r"""
    Tests whether a directed graph is acyclic and transitively
    reduced.

    EXAMPLES::

        sage: from sage.combinat.posets.posets import is_poset
        sage: dig = DiGraph({0:[2,3], 1:[3,4,5], 2:[5], 3:[5], 4:[5]})
        sage: is_poset(dig)
        False
        sage: is_poset(dig.transitive_reduction())
        True
    """
    return dig.is_directed_acyclic() and dig.is_transitively_reduced()
