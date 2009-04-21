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

import __builtin__
from copy import copy, deepcopy
from sage.structure.parent_base import ParentWithBase
from sage.rings.integer import Integer
from sage.graphs.graph import DiGraph
from sage.combinat.posets.hasse_diagram import HasseDiagram
from sage.combinat.posets.elements import *
from random import random

def Poset(data, element_labels=None, cover_relations=True):
    r"""
    Construct a poset from various forms of input data.

    INPUT:


    -  ``data`` - can be any of the following.

       1. a dictionary of upper covers: data[x] is a list of the
          elements that cover x.

       2. a list of upper covers,

       3. a tuple of elements and relations. Relations are two-element
          lists or tuples. (If these are not cover relations, then set
          cover_relations=False)

       4. a suitable digraph: transitively-reduced, directed, acyclic,
          no loops nor multiple edges.

       5. a previously constructed poset.

       6. In combination with cover_relations=False, data can be any
          directed acyclic graph without loops or multiple edges; the
          transitive reduction will be computed.

    -  ``element_labels`` - an optional list or dictionary
       of objects that label the poset elements.

    -  ``cover_relations`` - If True, then the data is
       assumed to describe a directed acyclic graph whose arrows are cover
       relations. If False, then Poset first determines the cover
       relations.


    OUTPUT:


    -  ``FinitePoset`` - an instance of the FinitePoset
       class.


    EXAMPLES: 1. A dictionary of upper covers.

    ::

        sage: Poset({'a':['b','c'], 'b':['d'], 'c':['d'], 'd':[]})
        Finite poset containing 4 elements

    2. A list of upper covers.

    ::

        sage: Poset([[1,2],[4],[3],[4],[]])
        Finite poset containing 5 elements

    A list of upper covers and a dictionary of labels.

    ::

        sage: elm_labs = {0:"a",1:"b",2:"c",3:"d",4:"e"}
        sage: P = Poset([[1,2],[4],[3],[4],[]],elm_labs)
        sage: P.list()
        [a, b, c, d, e]

    3. Elements and cover relations.

    ::

        sage: Poset(([1,2,3,4,5,6,7],[[1,2],[3,4],[4,5],[2,5]]))
        Finite poset containing 7 elements

    Elements and non-cover relations.

    ::

        sage: elms = [1,2,3,4]
        sage: rels = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
        sage: P = Poset( (elms,rels) ,cover_relations=False); P
        Finite poset containing 4 elements
        sage: P.cover_relations()
        [[1, 2], [2, 3], [3, 4]]

    4. An acyclic DiGraph.

    ::

        sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
        sage: Poset(dag)
        Finite poset containing 6 elements

    5. A previously constructed poset (returns the poset).

    ::

        sage: P = Poset([[1,2],[4],[3],[4],[]])
        sage: Poset(P)
        Finite poset containing 5 elements
        sage: Poset(P) == P
        True

    6. Any directed acyclic graph without loops or multiple edges as
    long as cover_relations=False.

    ::

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
    """
    D = {}
    if data is None: # type 0: empty poset.
        D = DiGraph()
    elif isinstance(data,dict): # type 1: dict of upper covers.
        D = DiGraph(data)
    elif isinstance(data,list): # type 2: list of upper covers.
        if [] in data:
            D = DiGraph(dict([[i,data[i]] for i in range(len(data))]))
    elif isinstance(data,tuple): # type 3: elements and relations.
        elements,relations = data
        for r in relations:
            if len(r) != 2:
                raise ValueError, \
                  "data[1](=%s) is not a list of relations"%relations
        D = DiGraph()
        D.add_vertices(elements)
        D.add_edges(relations)
    elif isinstance(data,DiGraph): # type 4: DiGraph
        D = deepcopy(data)
    elif isinstance(data,FinitePoset): # type 5: FinitePoset
        return data
    else:
        raise ValueError, "not valid poset data."

    # At this point D should be a DiGraph.
    if not isinstance(D,DiGraph):
        raise TypeError, "BUG: D should be a digraph."

    # Determine cover relations, if necessary.
    if cover_relations is False: # type 6: non-cover relations.
        D = D.transitive_reduction()

    # Check that the digraph does not contain loops, multiple edges
    # and is transitively reduced.
    if D.has_loops():
        raise ValueError, "Hasse diagram contains loops."
    elif D.has_multiple_edges():
        raise ValueError, "Hasse diagram contains multiple edges."
    elif cover_relations is True \
            and not D.transitive_reduction() == D:
        raise ValueError, "Hasse diagram is not transitively reduced."

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
    else:
        elements = [element_labels[z] for z in lin_ext]

    return FinitePoset(D,elements)

class FinitePoset(ParentWithBase):
    def __init__(self,digraph,elements=None):
        """
        Constructs a (finite) n-element poset from a set of elements and a
        directed acyclic graph.

        INPUT:


        -  ``digraph`` - a transitively-reduced, directed,
           acyclic graph without loops or multiple edges, and with vertices
           indexed by range(n). We assume that range(n) is a linear extension
           of the poset. (See DiGraph.relabel.)

        -  ``elements`` - an optional list of elements, with
           element[i] corresponding to vertex i. If elements==None, then it is
           set to be the vertex set of the digraph.


        EXAMPLES::

            sage: uc = [[2,3], [], [1], [1], [1], [3,4]]
            sage: elms = ["a","b","c","d","e","f"]
            sage: dag = DiGraph(dict([[i,uc[i]] for i in range(len(uc))]))
            sage: dag.relabel(elms)
            sage: P = Poset(dag); P
            Finite poset containing 6 elements
        """
        # A FinitePoset is a parent object; initialize it.
        ParentWithBase.__init__(self,None)

        # Construct the Hasse diagram.
        if isinstance(digraph,FinitePoset):
            self._hasse_diagram = digraph.hasse_diagram()
        else:
            self._hasse_diagram = HasseDiagram(digraph)

        # Store the elements of the poset.
        if elements is None:
            elements = range(self._hasse_diagram.order())
        self._elements = elements

        # Set the element type
        classname = self.__class__.__name__
        if classname.startswith('Finite'):
            element_type = classname[6:] + "Element"
        else:
            element_type = classname + "Element"
        self._element_type = eval(element_type)

    def _element_to_vertex(self,x):
        """
        Given an element of the poset, returns the corresponding vertex in
        the Hasse diagram.

        EXAMPLES::

            sage: P = Poset([[2,3], [], [1], [1], [1], [3,4]])
            sage: x = P.list()[3]
            sage: P._element_to_vertex(x)
            3
            sage: P._vertex_to_element(P._element_to_vertex(x)) == x
            True
        """
        try:
            if isinstance(x,self._element_type):
                return x.vertex
            else:
                return self._elements.index(x)
        except ValueError:
            raise ValueError, "element (=%s) not in poset"%x

    def _vertex_to_element(self,vertex):
        """
        Given a vertex of the Hasse diagram, returns the corresponding
        PosetElement.

        EXAMPLES::

            sage: P = Poset([[2,3], [], [1], [1], [1], [3,4]])
            sage: P._element_to_vertex(P._vertex_to_element(3)) == 3
            True
        """
        return self._element_type(self, self._elements[vertex], vertex)

    def __contains__(self,x):
        r"""
        Returns True if x is an element of the poset.

        TESTS::

            sage: P5 = Poset({(5,):[(4,1),(3,2)], \
                    (4,1):[(3,1,1),(2,2,1)], \
                    (3,2):[(3,1,1),(2,2,1)], \
                    (3,1,1):[(2,1,1,1)], \
                    (2,2,1):[(2,1,1,1)], \
                    (2,1,1,1):[(1,1,1,1,1)], \
                    (1,1,1,1,1):[]})
            sage: x = P5.list()[3]
            sage: x in P5
            True
        """
        if isinstance(x,self._element_type):
            if x.parent() is not self:
                return False
            y = x.element
        else:
            y = x
        return y in self._elements

    def __call__(self,element):
        """
        Returns the element labelled by n, if n is a label, or the n-th
        element of self.linear_extension.

        EXAMPLES::

            sage: P = Poset({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P(5)
            5
            sage: P(5) in P.minimal_elements()
            False
            sage: P(5) == P(-1)
            True
            sage: Q = Poset({5:[2,3], 1:[3,4], 2:[0], 3:[0], 4:[0]})
            sage: Q(5)
            5
            sage: Q(5) in Q.minimal_elements()
            True
            sage: Q(5) == Q(-6)
            True
            sage: R = Poset({'a':['b','c'], 'b':['d'], 'c':['d'], 'd':[]})
            sage: R(0)
            a
            sage: R('a') == R(0)
            True
            sage: R('d') == R(-1)
            True
        """
        if isinstance(element, self._element_type) and element.parent() == self:
            return element
        elif element in self._elements:
            return self._element_type(self, element, self._elements.index(element))
        elif isinstance(element,Integer):
            if element > -1:
                return self._element_type(self, \
                        self._elements[element], element)
            else:
                return self._element_type(self, \
                        self._elements[element], self.size()+element)
        else:
            raise ValueError, "__call__ accepts a poset element or an integer; you passed %s"%type(element)

    def hasse_diagram(self):
        """
        Returns the hasse_diagram of the poset.

        EXAMPLES::

            sage: Q = Poset({5:[2,3], 1:[3,4], 2:[0], 3:[0], 4:[0]})
            sage: Q.hasse_diagram()
            Hasse diagram of a poset containing 6 elements
        """
        return self._hasse_diagram

    def __repr__(self):
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
            sage: P5.__repr__()
            'Finite poset containing 7 elements'
        """
        return "Finite poset containing %s elements"%self.hasse_diagram().order()

    def __iter__(self):
        """
        Iterates through the elements of a linear extension of the poset.
        """
        return self.linear_extension().__iter__()

    def linear_extension(self):
        """
        Returns a linear extension of the poset.

        EXAMPLES::

            sage: from sage.combinat.posets.poset_examples import BooleanLattice
            sage: B = BooleanLattice(3)
            sage: B.linear_extension()
            [0, 1, 2, 3, 4, 5, 6, 7]
        """
        return map(self._vertex_to_element,range(self.size()))

    def linear_extensions(self):
        """
        Returns a list of all the linear extensions of the poset.

        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: D.linear_extensions()
            [[0, 1, 2, 3, 4], [0, 1, 2, 4, 3], [0, 2, 1, 3, 4], [0, 2, 1, 4, 3], [0, 2, 4, 1, 3]]
        """
        return [map(self._vertex_to_element,lin_ext) for lin_ext in
                self.hasse_diagram().linear_extensions()]

    def list(self):
        """
        List the elements of the poset. This just returns
        linear_extension().

        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: D.list()
            [0, 1, 2, 3, 4]
            sage: type(D.list()[0])
            <class 'sage.combinat.posets.elements.PosetElement'>
        """
        return self.linear_extension()

    def plot(self, label_elements=True, element_labels=None,
            label_font_size=12,label_font_color='black',
            vertex_size=300, vertex_colors=None,**kwds):
        """
        Returns a Graphic object corresponding the Hasse diagram of the
        poset. Optionally, it is labelled.

        INPUT:


        -  ``label_elements`` - whether to display element
           labels

        -  ``element_labels`` - a dictionary of element
           labels


        EXAMPLES::

            sage: D = Poset({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: D.plot(label_elements=False)
            sage: D.plot()
            sage: type(D.plot())
            <class 'sage.plot.plot.Graphics'>
            sage: elm_labs = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e'}
            sage: D.plot(element_labels=elm_labs)

        ::

            sage: P = Poset({})
            sage: P.plot()

        ::

            sage: P = Poset(DiGraph('E@ACA@?'))
            sage: P.plot()
        """
        if label_elements and element_labels is None:
            element_labels = self._elements
        return self.hasse_diagram().plot(label_elements=label_elements,
                        element_labels=element_labels,
                        label_font_size=label_font_size,
                        label_font_color=label_font_color,
                        vertex_size=vertex_size,vertex_colors=vertex_colors,**kwds)

    def show(self, label_elements=True, element_labels=None,
            label_font_size=12,label_font_color='black',
            vertex_size=300, vertex_colors=None,**kwds):
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
            vertex_size=vertex_size, vertex_colors=vertex_colors).show(**kwds)

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
                self.hasse_diagram().level_sets()]

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
        for u,v,l in self.hasse_diagram().edge_iterator():
            yield map(self._vertex_to_element,(u,v))

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
        return self.hasse_diagram().is_lequal(*map(self._element_to_vertex,(x,y)))

    def is_less_than(self, x, y):
        if x == y:
            return False
        else:
            return self.is_lequal(x,y)

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
        return self.hasse_diagram().is_lequal(*map(self._element_to_vertex,(y,x)))

    def is_greater_than(self, x, y):
        if x == y:
            return False
        else:
            return self.is_gequal(x,y)

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
        return map(self._vertex_to_element, self.hasse_diagram().minimal_elements())

    def maximal_elements(self):
        """
        Returns a list of the maximal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.maximal_elements()
            [4]
        """
        return map(self._vertex_to_element, self.hasse_diagram().maximal_elements())

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
        hasse_bot = self.hasse_diagram().bottom()
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
        return self.hasse_diagram().has_bottom()

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
        """
        hasse_top = self.hasse_diagram().top()
        if hasse_top:
            return self._vertex_to_element(hasse_top)
        else:
            return None

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
        return self.hasse_diagram().has_top()

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
        return self.hasse_diagram().is_bounded()

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
        return self.hasse_diagram().is_chain()

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
        hasse_rf = self.hasse_diagram().rank_function()
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

        ::

            sage: P = SymmetricGroupBruhatOrderPoset(4)
            sage: [(v,P.rank(v)) for v in P]
            [(1234, 0),
             (2134, 1),
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
        return self.hasse_diagram().has_edge(*map(self._element_to_vertex,(x,y)))

    def upper_covers_iterator(self,y):
        """
        Returns an iterator for the upper covers of the element y. An upper
        cover of y is an element x such that y x is a cover relation.

        EXAMPLES::

            sage: Q = Poset({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: type(Q.upper_covers_iterator(0))
            <type 'generator'>
        """
        for x in self.hasse_diagram().successor_iterator(self._element_to_vertex(y)):
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
        for x in self.hasse_diagram().predecessor_iterator(self._element_to_vertex(y)):
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

    def size(self):
        """
        Returns the number of elements in the poset.

        EXAMPLES::

            sage: Poset([[1,2,3],[4],[4],[4],[]]).size()
            5
        """
        return self.hasse_diagram().order()

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
            sage: Q.mobius_function(Q(0),Q(-1))
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
        return self.hasse_diagram().mobius_function(i,j)

    def mobius_function_matrix(self):
        r"""
        Returns a matrix whose (i,j) entry is the value of the Mobius
        function evaluated at self.linear_extension()[i] and
        self.linear_extension()[j].

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
        """
        return self.hasse_diagram().mobius_function_matrix()

    def lequal_matrix(self,**kwds):
        """
        Computes the matrix whose [i,j] entry is 1 if
        self.linear_extension()[i] self.linear_extension()[j] 0
        otherwise.

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
        """
        return self.hasse_diagram().lequal_matrix(**kwds)

    def meet_matrix(self):
        """
        Returns a matrix whose (i,j) entry is k, where
        self.linear_extension()[k] is the meet (greatest lower bound) of
        self.linear_extension()[i] and self.linear_extension()[j].

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
        return self.hasse_diagram().meet_matrix()

    def is_meet_semilattice(self):
        """
        Returns True is the poset has a meet operation, and False
        otherwise.
        """
        return self.hasse_diagram().is_meet_semilattice()

    def join_matrix(self):
        """
        Returns a matrix whose (i,j) entry is k, where
        self.linear_extension()[k] is the join (least upper bound) of
        self.linear_extension()[i] and self.linear_extension()[j].

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
        return self.hasse_diagram().join_matrix()

    def is_join_semilattice(self):
        """
        Returns True is the poset has a join operation, and False
        otherwise.
        """
        return self.hasse_diagram().is_join_semilattice()

    def antichains(self):
        """
        Returns a list of all antichains of the poset.

        An antichain of a poset is a collection of elements of the poset
        that are pairwise incomparable.

        EXAMPLES::

            sage: PentagonPoset().antichains()
            [[], [0], [1], [2], [3], [4], [1, 2], [1, 3]]
            sage: AntichainPoset(3).antichains()
            [[], [2], [1], [0], [1, 0], [2, 1], [2, 0], [2, 1, 0]]
            sage: ChainPoset(3).antichains()
            [[], [0], [1], [2]]
        """
        return [map(self._vertex_to_element,antichain) for
                antichain in self.hasse_diagram().antichains()]

    def dual(self):
        """
        Returns the dual poset of the given poset.

        EXAMPLE::

            sage: P = Poset([[1,2],[4],[3],[4],[]])
            sage: P.dual()
            Finite poset containing 5 elements
        """
        dual_graph = self.hasse_diagram().reverse()
        dual_graph.relabel(self._elements)
        return FinitePoset(dual_graph)

    def graphviz_string(self,graph_string="graph",edge_string="--"):
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
        EXAMPLES::

            sage: P = Poset({"a":["c","d"], "b":["d","e"], "c":["f"], "d":["f"], "e":["f"]})
            sage: P.subposet(["a","b","f"])
            Finite poset containing 3 elements
        """
        if not isinstance(elements,list):
            raise ValueError, "not a list."
        for element in elements:
            if element not in self:
                raise ValueError, "element not in self"
        relations = []
        for u in elements:
            for v in elements:
                if self(u) < self(v): relations.append([u,v])
        return Poset((elements, relations), cover_relations=False)

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
            if random() <= p:
                elements.append(v)
        return self.subposet(elements)

    def order_filter(self,elements):
        """
        Returns the order filter generated by a list of elements.

        I is an order filter if it satisfies: x in I and y = x implies y in
        I

        EXAMPLES::

            sage: B = BooleanLattice(4)
            sage: B.order_filter([3,8])
            [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        vertices = map(self._element_to_vertex,elements)
        of = self.hasse_diagram().order_filter(vertices)
        return map(self._vertex_to_element,of)

    def principal_order_filter(self, x):
        """
        Returns the order filter generated by an element x.

        EXAMPLES::

            sage: B = BooleanLattice(4)
            sage: B.principal_order_filter(2)
            [2, 3, 6, 7, 10, 11, 14, 15]
        """
        return self.order_filter([x])

    def order_ideal(self,elements):
        """
        Returns the order ideal generated by a list of elements.

        I is an order ideal if it satisfies: x in I and y = x implies y in
        I

        EXAMPLES::

            sage: B = BooleanLattice(4)
            sage: B.order_ideal([7,10])
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
        """
        vertices = map(self._element_to_vertex,elements)
        oi = self.hasse_diagram().order_ideal(vertices)
        return map(self._vertex_to_element,oi)

    def principal_order_ideal(self, x):
        """
        Returns the order ideal generated by an element x.

        EXAMPLES::

            sage: B = BooleanLattice(4)
            sage: B.principal_order_ideal(6)
            [0, 2, 4, 6]
        """
        return self.order_ideal([x])

    def interval(self, x, y):
        """
        Returns a list of the elements z such that x = z = y. The order is
        that induced by the ordering in self.linear_extension().

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
        """
        return map(self._vertex_to_element,self.hasse_diagram().interval(x,y))

    def closed_interval(self, x, y):
        """
        Returns a list of the elements z such that x = z = y. The order is
        that induced by the ordering in self.linear_extension().

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: P = Poset(dag)
            sage: I = set(map(P,[2,5,6,4,7]))
            sage: I == set(P.interval(2,7))
            True
        """
        return self.interval(x,y)

    def open_interval(self, x, y):
        """
        Returns a list of the elements z such that x z y. The order is that
        induced by the ordering in self.linear_extension().

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: P = Poset(dag)
            sage: I = set(map(P,[5,6,4]))
            sage: I == set(P.open_interval(2,7))
            True
        """
        return map(self._vertex_to_element,self.hasse_diagram().open_interval(x,y))

    def _all_maximal_chains(self, partial=None):
        """
        Returns list of the maximal chains of this poset.  Each chain
        is listed in increasing order.

        INPUT:


        -  ``partial`` - list (optional).  If present, find all maximal
           chains starting with the elements in partial.

        Returns list of the maximal chains of this poset.

        This is used in constructing the order complex for the poset.

        EXAMPLES::

            sage: P = BooleanLattice(3)
            sage: P._all_maximal_chains()
            [[0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7], [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 6, 7]]
            sage: P._all_maximal_chains(partial=[0,2])
            [[0, 2, 3, 7], [0, 2, 6, 7]]
            sage: Q = ChainPoset(6)
            sage: Q._all_maximal_chains()
            [[0, 1, 2, 3, 4, 5]]
        """
        chains = []
        if partial is None or len(partial) == 0:
            start = self.minimal_elements()
            partial = []
        else:
            start = self.upper_covers(partial[-1])
        if len(start) == 0:
            return [partial]
        if len(start) == 1:
            return self._all_maximal_chains(partial=partial + start)
        parts = [partial + [x] for x in start]
        answer = []
        for new in parts:
            answer += self._all_maximal_chains(partial=new)
        return answer

    def order_complex(self):
        """
        Returns the order complex associated to this poset.

        The order complex is the simplical complex with vertices equal
        to the elements of the poset, and faces given by the chains.

        EXAMPLES::

            sage: P = BooleanLattice(3)
            sage: S = P.order_complex(); S
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6, 7) and 6 facets
            sage: S.f_vector()
            [1, 8, 19, 18, 6]
            sage: S.homology()      # S is contractible
            {0: 0, 1: 0, 2: 0, 3: 0}
            sage: Q = P.subposet([1,2,3,4,5,6])
            sage: Q.order_complex().homology()    # a circle
            {0: 0, 1: Z}
        """
        from sage.homology.simplicial_complex import SimplicialComplex
        vertices = [a.element for a in self.list()]
        facets = []
        for f in self._all_maximal_chains():
            facets.append([a.element for a in f])
        return SimplicialComplex(vertices, facets)
