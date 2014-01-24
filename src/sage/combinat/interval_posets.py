r"""
Tamari Interval-posets

This module implements the combinatorial object Tamari interval-poset which 
represents an interval of the Tamari order. It has been introduced in [PCH]_
and allows for many combinatorial operations on Tamari intervals. In particular,
it is linked to :class:`DyckWords` and :class:`BinaryTrees`.

REFERENCES:

    .. [PCH] Counting smaller trees in the Tamari order, G. Chatel, V. Pons, 2013

**AUTHORS:**

- Viviane Pons 2014: initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2013 Viviane Pons <viviane.pons@univie.ac.at>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.finite_posets import FinitePosets
from sage.categories.finite_posets import Posets
from sage.combinat.binary_tree import BinaryTree
from sage.combinat.binary_tree import BinaryTrees
from sage.combinat.binary_tree import LabelledBinaryTrees
from sage.combinat.dyck_word import DyckWords
from sage.combinat.permutation import Permutation
from sage.combinat.posets.posets import Poset
from sage.combinat.posets.posets import FinitePoset 
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.rings.integer import Integer
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.structure.element import Element
from sage.structure.global_options import GlobalOptions
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

TamariIntervalPosetOptions = GlobalOptions(name="Tamari Interval-posets",
    doc=r"""
    Set and display the global options for Tamari interval-posets. If no parameters
    are set, then the function returns a copy of the options dictionary.

    The ``options`` to Tamari interval-posets can be accessed as the method
    :obj:`TamariIntervalPosets.global_options` of :class:`TamariIntervalPosets` and
    related parent classes.
    """,
    end_doc=r"""
    EXAMPLES::

        sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
        sage: ip.latex_options()["color_decreasing"]
        'red'
        sage: TamariIntervalPosets.global_options(latex_color_decreasing='green')
        sage: ip.latex_options()["color_decreasing"]
        'green'
        sage: TamariIntervalPosets.global_options.reset()
        sage: ip.latex_options()["color_decreasing"]
        'red'

    """,
    latex_tikz_scale=dict(default=1,
                          description='The default value for the tikz scale when latexed',
                          checker=lambda x: True), # More trouble than it's worth to check
    latex_line_width_scalar=dict(default=0.5,
                                 description='The default value for the line width as a'
                                             'multiple of the tikz scale when latexed',
                                 checker=lambda x: True), # More trouble than it's worth to check
    latex_color_decreasing = dict(default="red",
                                  description='The default color of decreasing relations when latexed',
                                  checker=lambda x: True), # More trouble than it's worth to check
    latex_color_increasing = dict(default="blue",
                                  description='The default color of increasing relations when latexed',
                                  checker=lambda x: True), # More trouble than it's worth to check
    latex_hspace = dict(default=1,
                        description='The default difference between horizontal coordinates of vertices when latexed',
                        checker=lambda x: True), # More trouble than it's worth to check
    latex_vspace = dict(default=1,
                        description='The default difference between vertical coordinates of vertices when latexed',
                        checker=lambda x: True), # More trouble than it's worth to check
)

class TamariIntervalPoset(Element):
    r"""
    The class of Tamari Interval-posets.
    
    An interval-poset is a labelled poset of size n, with labelled `1,\dots,n`
    satisfying the following conditions:

    - if `a < c` (as a number) and `a` precedes `c` in the poset, then, 
      for all `b` such that `a<b<c`, `b` precedes `c`,
    
    - if `a<c` (as a number) and `c` precedes `a` in the poset, then,
      for all `b` such that `a<b<c`, `b` precedes `a`.

    They are in bijection with intervals of the Tamari lattice.

    INPUT:

    - ``size`` -- an integer, the size of the interval-posets (number of 
      vertices)

    - ``relations`` -- an iterable of couples (a,b) (list, tuple or iterable) 
      representing a relation 'a precedes b' in the poset. 
      
    - ``check`` -- (default: True) whether to check the interval-poset 
      condition or not.

    EXAMPLES::

        sage: TamariIntervalPoset(0,[])
        The tamari interval of size 0 induced by relations []
        sage: TamariIntervalPoset(3,[])
        The tamari interval of size 3 induced by relations []
        sage: TamariIntervalPoset(3,[(1,2)])
        The tamari interval of size 3 induced by relations [(1, 2)]
        sage: TamariIntervalPoset(3,[(1,2),(2,3)])
        The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[(1,2),(2,3),(1,3)])
        The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[(1,2),(3,2)])
        The tamari interval of size 3 induced by relations [(1, 2), (3, 2)]

        sage: TamariIntervalPoset(3,[(3,4)])
        Traceback (most recent call last):
        ...
        ValueError: The relations do not correspond to the size of the poset.

        sage: TamariIntervalPoset(2,[(2,1),(1,2)])
        Traceback (most recent call last):
        ...
        ValueError: Hasse diagram contains cycles.
        
        sage: TamariIntervalPoset(3,[(1,3)])
        Traceback (most recent call last):
        ...
        ValueError: This does not satisfy the Tamari interval-poset condition.

    It is also possible to transform a poset directly into an interval-poset::
        
        sage: TIP = TamariIntervalPosets()
        sage: p = Poset( ([1,2,3], [(1,2)]))
        sage: TIP(p)
        The tamari interval of size 3 induced by relations [(1, 2)]

    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        Ensure that interval-posets created by the enumerated sets and directly
        are the same and that they are instances of :class:`TamariIntervalPoset`

        TESTS::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.parent()
            Interval-posets
            sage: type(ip)
            <class 'sage.combinat.interval_posets.TamariIntervalPosets_all_with_category.element_class'>

            sage: ip2 = TamariIntervalPosets()(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip2.parent() is ip.parent()
            True
            sage: type(ip) is type(ip2)
            True

            sage: ip3 = TamariIntervalPosets(4)([(2,4),(3,4),(2,1),(3,1)])
            sage: ip3.parent() is ip.parent()
            True
            sage: type(ip3) is type(ip)
            True
        """
        return TamariIntervalPosets_all().element_class(*args, **opts)

    def __init__(self, size, relations, check=True):
        r"""
        TESTS::

            sage: TamariIntervalPoset(3,[(1,2),(3,2)]).parent()
            Interval-posets

        """
        parent = TamariIntervalPosets()
        self._size = size
        self._poset = Poset( ([i for i in xrange(1,size+1)], relations) )
        if(self._poset.cardinality()!=size):
            raise ValueError, "The relations do not correspond to the size of the poset."%()

        if check and not TamariIntervalPosets.check_poset(self._poset):
            raise ValueError, "This does not satisfy the Tamari interval-poset condition."%()
            
        Element.__init__(self, parent)

        self._cover_relations = tuple(self._poset.cover_relations())

        self._latex_options = dict()

    def set_latex_options(self, D):
        r"""
        Set the latex options for use in the ``_latex_`` function.  The
        default values are set in the ``__init__`` function.

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package.

        - ``line_width`` -- (default: 1*``tikz_scale``) value representing the
          line width.

        - ``color_decreasing`` -- (default: red) the color for decreasing relations.

        - ``color_increasing`` -- (default: blue) the color for increasing relations.

        - ``hspace`` -- (default: 1) the difference between horizontal coordinates of adjacent vertices
        
        - ``vpsace`` -- (default: 1) the difference between vertical coordinates of adjacent vertices

        INPUT:

        - ``D`` -- a dictionary with a list of latex parameters to change

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.latex_options()["color_decreasing"]
            'red'
            sage: ip.set_latex_options({"color_decreasing":'green'})
            sage: ip.latex_options()["color_decreasing"]
            'green'
            
        To change the options for all interval-posets, use the parent's
        global latex options::
        
            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.latex_options()["color_decreasing"]
            'red'
            sage: ip2.latex_options()["color_decreasing"]
            'red'
            sage: TamariIntervalPosets.global_options(latex_color_decreasing='green')
            sage: ip.latex_options()["color_decreasing"]
            'green'
            sage: ip2.latex_options()["color_decreasing"]
            'green'
            sage: TamariIntervalPosets.global_options.reset()

        """
        for opt in D:
            self._latex_options[opt] = D[opt]

    def latex_options(self):
        r"""
        Return the latex options for use in the ``_latex_`` function as a
        dictionary. The default values are set using the global options.

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package.

        - ``line_width`` -- (default: 1*``tikz_scale``) value representing the
          line width.

        - ``color_decreasing`` -- (default: red) the color for decreasing relations.

        - ``color_increasing`` -- (default: blue) the color for increasing relations.

        - ``hspace`` -- (default: 1) the difference between horizontal coordinates of adjacent vertices
        
        - ``vpsace`` -- (default: 1) the difference between vertical coordinates of adjacent vertices

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.latex_options()['color_decreasing']
            'red'
            sage: ip.latex_options()['hspace']
            1

        """
        d = self._latex_options.copy()
        if "tikz_scale" not in d:
            d["tikz_scale"] = self.parent().global_options["latex_tikz_scale"]
        if "line_width" not in d:
            d["line_width"] = self.parent().global_options["latex_line_width_scalar"]*d["tikz_scale"]
        if "color_decreasing" not in d:
            d["color_decreasing"] = self.parent().global_options["latex_color_decreasing"]
        if "color_increasinf" not in d:
            d["color_increasing"] = self.parent().global_options["latex_color_increasing"]
        if "hspace" not in d:
            d["hspace"] = self.parent().global_options["latex_hspace"]
        if "vspace" not in d:
            d["vspace"] = self.parent().global_options["latex_vspace"]
        return d
        
    def _latex_(self):
        r"""
        A latex representation of ``self`` using the tikzpicture package.
        
        If `x` precedes `y` than, `y` will always be placed on top of `x` and / or on the right of `x`.
        Decreasing relations are drawn vertically and increasing relations horizontally.
        The algorithm tries to avoid superposition but on big interval-posets, it might happen.
        
        You can use ``self.set_latex_options()`` to change default latex options. 
        Or you can use the parent's global options.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: print ip._latex_()
            \begin{tikzpicture}[scale=1]
            \node(T2) at (0,-1) {2};
            \node(T3) at (1,-2) {3};
            \node(T1) at (1,0) {1};
            \node(T4) at (2,-1) {4};
            \draw[line width = 0.5, color=blue] (T2) -- (T4);
            \draw[line width = 0.5, color=red] (T2) -- (T1);
            \draw[line width = 0.5, color=blue] (T3) -- (T4);
            \draw[line width = 0.5, color=red] (T3) -- (T1);
            \end{tikzpicture}
        """
        latex.add_package_to_preamble_if_available("tikz")
        latex_options = self.latex_options()
        start = "\\begin{tikzpicture}[scale="+str(latex_options['tikz_scale'])+"]\n"
        end = "\\end{tikzpicture}"
        def draw_node(j,x,y):
            r"""
            Internal method to draw vertices
            """
            return "\\node(T" + str(j) + ") at (" + str(x) + "," + str(y) + ") {" + str(j) + "};\n"
        def draw_increasing(i,j):
            r"""
            Internal method to draw increasing relations
            """
            return "\\draw[line width = " + str(latex_options["line_width"]) + ", color=" + latex_options["color_increasing"] + "] (T" + str(i) + ") -- (T" + str(j) +");\n"
        def draw_decreasing(i,j):
            r"""
            Internal method to draw decreasing relations
            """
            return "\\draw[line width = " + str(latex_options["line_width"]) + ", color=" + latex_options["color_decreasing"] + "] (T" + str(i) + ") -- (T" + str(j) +");\n"    
        if self.size() == 0:
            nodes = "\\node(T0) at (0,0){$\emptyset$};"
            relations = ""
        else:
            nodes = "" # latex for node decraltions
            relations = "" # latex for drawing relations
            to_draw = []
            to_draw.append((1,0)) # a pilo of nodes to be drawn with their y position

            current_parent = [self.increasing_parent(1)] # a pilo for the current increasing parents
            parenty = [0] # a pilo for the current parent y positions
            if current_parent[-1]!= None:
                relations+= draw_increasing(1,current_parent[-1])
            vspace = latex_options["vspace"]
            hspace = latex_options["hspace"]
            x = 0
            y = 0
            
            # the idea is that we draw the nodes from left to right and save their y position
            for i in xrange(2,self.size()+1): 
                # at each, we draw all possible nodes and add the current node to the to_draw pilo
                decreasing_parent = self.decreasing_parent(i)
                increasing_parent = self.increasing_parent(i)
                while len(to_draw)>0 and (decreasing_parent is None or decreasing_parent < to_draw[-1][0]):
                    # we draw all the nodes which can be placed at x
                    # we know these nodes won't have any more decreasing children (so their horizontal position is fixed)
                    n = to_draw.pop()
                    nodes+= draw_node(n[0],x,n[1])
                if i != current_parent[-1]:
                    #i is not the current increasing parent
                    if(not self.le(i,i-1) and decreasing_parent != None):
                        # there is no decreasing relation between i and i-1
                        #they share a decreasing parent and are placed alongside horizontally
                        x+=hspace
                        if current_parent[-1] != None:
                            y-= vspace
                    else:
                        #otherwise, they are placed alongside vertically
                        y-=vspace
                    if increasing_parent != current_parent[-1]:
                        current_parent.append(increasing_parent)
                        parenty.append(y)
                    nodey = y
                else:
                    # i is the current increasing parent so it takes the current vertical position
                    current_parent.pop()
                    x+=hspace 
                    nodey = parenty.pop()
                    if len(current_parent)==0 or increasing_parent != current_parent[-1]:
                        current_parent.append(increasing_parent)
                        parenty.append(nodey)
                to_draw.append((i,nodey))
                if increasing_parent != None:
                    relations+=draw_increasing(i,increasing_parent)
                if decreasing_parent != None:
                    relations+=draw_decreasing(i,decreasing_parent)
            for n in to_draw:
                # we draw all remaining nodes
                nodes+= draw_node(n[0],x,n[1])
        return start + nodes + relations + end


    @cached_method
    def increasing_cover_relations(self):
        r"""
        Return the cover relations of the increasing poset of ``self`` (the poset
        formed by keeping only relations `a` precedes `b` with `a<b`)

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)]).increasing_cover_relations()
            [(1, 2), (2, 4), (3, 4)]
            sage: TamariIntervalPoset(3,[(1,2),(1,3),(2,3)]).increasing_cover_relations()
            [(1, 2), (2, 3)]

        """
        relations = []
        for i in xrange(1,self.size()):
            for j in xrange(i+1, self.size()+1):
                if self.le(i,j):
                    relations.append((i,j))
                    break
        return relations

    def increasing_roots(self):
        r"""
        Return the root vertices of the initial forest of ``self``, 
        i.e., the vertices `a` such that there is no `b>a` with `a` precedes `b`.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_roots()
            [6, 5, 2]
            sage: ip.initial_forest().increasing_roots()
            [6, 5, 2]
        """
        if self.size()==0:
            return []
        roots = [self.size()]
        root = self.size()
        for i in xrange(self.size()-1,0,-1):
            if not self.le(i,root):
                roots.append(i)
                root = i
        return roots

    def increasing_children(self, v):
        r"""
        Return the children of ``v`` in the initial forest of ``self``.

        INPUT:
        - ``v`` -- an integer representing a vertex of ``self`` (between 1 and ``size``)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_children(2)
            [1]
            sage: ip.increasing_children(5)
            [4, 3]
            sage: ip.increasing_children(1)
            []
        """
        children = []
        root = None
        for i in xrange(v-1,0,-1):
            if not self.le(i,v):
                break
            if root is None or not self.le(i,root):
                children.append(i)
                root = i
        return children
    
    def increasing_parent(self, v):
        r"""
        Return the vertex parent of ``v`` in the initial forest of ``self``.
        This is the minimal `b>v` such that ``v`` precedes ``b``. If there
        is no such vertex (``v`` is an increasing root), then ``None`` is 
        returned.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self`` (between 1 and ``size``)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_parent(1)
            2
            sage: ip.increasing_parent(3)
            5
            sage: ip.increasing_parent(4)
            5
            sage: ip.increasing_parent(5) is None
            True

        """
        parent = None
        for i in xrange(self.size(),v,-1):
            if self.le(v,i):
                parent = i
        return parent

    @cached_method
    def decreasing_cover_relations(self):
        r"""
        Return the cover relations of the increasing poset of ``self`` (the poset
        formed by keeping only relations a preceded b with a<b)

        EXAMPLES::
        
            sage: TamariIntervalPoset(4,[(2,1),(3,2),(3,4),(4,2)]).decreasing_cover_relations()
            [(4, 2), (3, 2), (2, 1)]
            sage: TamariIntervalPoset(4,[(2,1),(4,3),(2,3)]).decreasing_cover_relations()
            [(4, 3), (2, 1)]
            sage: TamariIntervalPoset(3,[(2,1),(3,1),(3,2)]).decreasing_cover_relations()
            [(3, 2), (2, 1)]
        """
        relations = []
        for i in xrange(self.size(),1,-1):
            for j in xrange(i-1,0,-1):
                if self.le(i,j):
                    relations.append((i,j))
                    break
        return relations
        
    def decreasing_roots(self):
        r"""
        Return the root vertices of the final forest of ``self`` , 
        i.e., the vertices `b` such that there is no `a<b` with `b` precedes `a`.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_roots()
            [1, 2]
            sage: ip.final_forest().decreasing_roots()
            [1, 2]
        """
        if self.size()==0:
            return []
        roots = [1]
        root = 1
        for i in xrange(2,self.size()+1):
            if not self.le(i,root):
                roots.append(i)
                root = i
        return roots

    def decreasing_children(self, v):
        r"""
        Return the children of ``v`` in the final forest of ``self``.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self`` (between 1 and ``size``)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_children(2)
            [3, 5]
            sage: ip.decreasing_children(3)
            [4]
            sage: ip.decreasing_children(1)
            []
        """
        children = []
        root = None
        for i in xrange(v+1,self.size()+1):
            if not self.le(i,v):
                break
            if root is None or not self.le(i,root):
                children.append(i)
                root = i
        return children
    
    def decreasing_parent(self, v):
        r"""
        Return the vertex parent of ``v`` in the final forest of ``self``.
        This is the maximal `a < v` such that ``v`` precedes ``a``. If there
        is no such vertex (``v`` is a decreasing root), then ``None`` is 
        returned.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self`` (between 1 and ``size``)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_parent(4)
            3
            sage: ip.decreasing_parent(3)
            2
            sage: ip.decreasing_parent(5)
            2
            sage: ip.decreasing_parent(2) is None
            True

        """
        parent = None
        for i in xrange(1,v):
            if self.le(v,i):
                parent = i
        return parent

    def le(self, e1, e2):
        r"""
        Return whether ``e1`` precedes or equals ``e2`` in ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.le(1,2)
            True
            sage: ip.le(1,3)
            True
            sage: ip.le(2,3)
            True
            sage: ip.le(3,4)
            False
            sage: ip.le(1,1)
            True

        """
        return self._poset.le(e1,e2)
        
    def lt(self, e1, e2):
        r"""
        Return whether ``e1`` strictly precedes ``e2`` in ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.lt(1,2)
            True
            sage: ip.lt(1,3)
            True
            sage: ip.lt(2,3)
            True
            sage: ip.lt(3,4)
            False
            sage: ip.lt(1,1)
            False

        """
        return self._poset.lt(e1,e2)
        
    def ge(self, e1, e2):
        r"""
        Return whether ``e2`` precedes or equals ``e1`` in ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.ge(2,1)
            True
            sage: ip.ge(3,1)
            True
            sage: ip.ge(3,2)
            True
            sage: ip.ge(4,3)
            False
            sage: ip.ge(1,1)
            True
        """
        return self._poset.ge(e1, e2)
        
    def gt(self, e1, e2):
        r"""
        Return whether ``e2`` strictly precedes ``e1`` in ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.gt(2,1)
            True
            sage: ip.gt(3,1)
            True
            sage: ip.gt(3,2)
            True
            sage: ip.gt(4,3)
            False
            sage: ip.gt(1,1)
            False
        """
        return self._poset.gt(e1,e2)
        
    def size(self):
        r"""
        Return the size (number of vertices) of the interval-poset.

        EXAMPLES::

            sage: TamariIntervalPoset(3,[(2,1),(3,1)]).size()
            3
        """
        return self._size
    
    def _repr_(self):
        r"""
        TESTS::

            sage: TamariIntervalPoset(3,[(2,1),(3,1)])
            The tamari interval of size 3 induced by relations [(3, 1), (2, 1)]
            sage: TamariIntervalPoset(3,[(3,1),(2,1)])
            The tamari interval of size 3 induced by relations [(3, 1), (2, 1)]
            sage: TamariIntervalPoset(3,[(2,3),(2,1)])
            The tamari interval of size 3 induced by relations [(2, 3), (2, 1)]
        """
        return "The tamari interval of size %s induced by relations %s"%(self.size(),str(self.increasing_cover_relations() + self.decreasing_cover_relations()))

    def __eq__(self, other):
        r"""
        TESTS::

            sage: TamariIntervalPoset(0,[]) == TamariIntervalPoset(0,[])
            True
            sage: TamariIntervalPoset(1,[]) == TamariIntervalPoset(0,[])
            False
            sage: TamariIntervalPoset(3,[(1,2),(3,2)]) == TamariIntervalPoset(3,[(3,2),(1,2)])
            True
            sage: TamariIntervalPoset(3,[(1,2),(3,2)]) == TamariIntervalPoset(3,[(1,2)])
            False

        """
        if(not isinstance(other, TamariIntervalPoset)):
            return False
        return self.size() == other.size() and self._cover_relations == other._cover_relations
        
    def __le__(self, el2):
        r"""
        TESTS::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 <= ip2
            True
            sage: ip1 <= ip1
            True
            sage: ip2 <= ip1
            False
        """
        return self.parent().le(self, el2)
    
    def __lt__(self, el2):
        r"""
        TESTS::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 < ip2
            True
            sage: ip1 < ip1
            False
            sage: ip2 < ip1
            False

        """
        return self.parent().lt(self, el2)
    
    def __ge__(self, el2):
        r"""
        TESTS::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 >= ip2
            False
            sage: ip1 >= ip1
            True
            sage: ip2 >= ip1
            True
        """
        return self.parent().ge(self,el2)
        
    def __gt__(self, el2):
        r"""
        TESTS::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 > ip2
            False
            sage: ip1 > ip1
            False
            sage: ip2 > ip1
            True

        """
        return self.parent().gt(self, el2)
        
    def contains_interval(self, other):
        r"""
        Return whether the interval represented by ``other`` is contained
        into ``self`` as an interval of the Tamari lattice. 
        
        In terms of interval-posets, it means that all relations of ``self``
        are relations of ``other``.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(2,3)])
            sage: ip2.contains_interval(ip1)
            True
            sage: ip3 = TamariIntervalPoset(4,[(2,1)])
            sage: ip2.contains_interval(ip3)
            False
            sage: ip4 = TamariIntervalPoset(3,[(2,3)])
            sage: ip2.contains_interval(ip4)
            False
        """
        if other.size() != self.size():
            return False
        for (i,j) in self._cover_relations:
            if not other.le(i,j):
                return False
        return True
    
    def lower_contains_interval(self,other):
        r"""
        Return whether the interval represented by ``other`` is contained
        into ``self`` as an interval of the Tamari lattice and if they share
        the same lower bound. 
        
        As interval-posets, it means that ``other`` contains the relations
        of ``self`` plus some extra increasing relations.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)]);
            sage: ip2 = TamariIntervalPoset(4,[(4,3)])
            sage: ip2.lower_contains_interval(ip1)
            True
            sage: ip2.contains_interval(ip1) and ip2.lower_binary_tree() == ip1.lower_binary_tree()
            True
            sage: ip3 = TamariIntervalPoset(4,[(4,3),(2,1)])
            sage: ip2.contains_interval(ip3)
            True
            sage: ip2.lower_binary_tree() == ip3.lower_binary_tree()
            False
            sage: ip2.lower_contains_interval(ip3)
            False
        """
        if not self.contains_interval(other):
            return False
        for (i,j) in other.decreasing_cover_relations():
            if not self.le(i,j):
                return False
        return True

    def upper_contains_interval(self,other):
        r"""
        Return whether the interval represented by ``other`` is contained
        into ``self`` as an interval of the Tamari lattice and if they share
        the same upper bound. 
        
        As interval-posets, it means that ``other`` contains the relations
        of ``self`` plus some extra decreasing relations.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip2.upper_contains_interval(ip1)
            True
            sage: ip2.contains_interval(ip1) and ip2.upper_binary_tree() == ip1.upper_binary_tree()
            True
            sage: ip3 = TamariIntervalPoset(4,[(1,2),(2,3),(3,4)])
            sage: ip2.upper_contains_interval(ip3)
            False
            sage: ip2.contains_interval(ip3)
            True
            sage: ip2.upper_binary_tree() == ip3.upper_binary_tree()
            False
        """
        if not self.contains_interval(other):
            return False
        for (i,j) in other.increasing_cover_relations():
            if not self.le(i,j):
                return False
        return True

    def is_linear_extension(self, perm):
        r"""
        Return whether ``perm`` is a linear extension of ``self``.

        INPUT:

        - ``perm`` -- a pemutation

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.is_linear_extension([1,4,2,3])
            True
            sage: ip.is_linear_extension(Permutation([1,4,2,3]))
            True
            sage: ip.is_linear_extension(Permutation([1,4,3,2]))
            False
        """
        return self._poset.is_linear_extension(perm)
    
    def contains_binary_tree(self, binary_tree):
        r"""
        Return whether the interval represented by ``self`` contains 
        the binary tree ``binary_tree``.

        INPUT:

        - ``binary_tree`` -- a binary tree

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.contains_binary_tree(BinaryTree([[None,[None,[]]],None]))
            True
            sage: ip.contains_binary_tree(BinaryTree([None,[[[],None],None]]))
            True
            sage: ip.contains_binary_tree(BinaryTree([[],[[],None]]))
            False
            sage: ip.contains_binary_tree(ip.lower_binary_tree())
            True
            sage: ip.contains_binary_tree(ip.upper_binary_tree())
            True
            sage: all([ip.contains_binary_tree(bt) for bt in ip.binary_trees()])
            True

        """
        return self.is_linear_extension(binary_tree.to_132_avoiding_permutation())
    
    def contains_dyck_word(self, dyck_word):
        r"""
        Return whether the interval represented by ``self`` contains 
        the Dyck word ``dyck_word``.

        INPUT:

        - ``dyck_word`` -- a Dyck word

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.contains_dyck_word(DyckWord([1,1,1,0,0,0,1,0]))
            True
            sage: ip.contains_dyck_word(DyckWord([1,1,0,1,0,1,0,0]))
            True
            sage: ip.contains_dyck_word(DyckWord([1,0,1,1,0,1,0,0]))
            False        
            sage: ip.contains_dyck_word(ip.lower_dyck_word())
            True
            sage: ip.contains_dyck_word(ip.upper_dyck_word())
            True
            sage: all([ip.contains_dyck_word(bt) for bt in ip.dyck_words()])
            True
        """
        return self.contains_binary_tree(dyck_word.to_binary_tree_tamari())

    def intersection(self, other):
        r"""
        Return the interval-poset formed of relations from both ``self``
        and ``other``. It corresponds to the intersection of the two 
        corresponding intervals of Tamari.
        
        INPUT:

        - ``other`` -- an interval-poset of the same size as ``self``

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip2 = TamariIntervalPoset(4,[(4,3)])
            sage: ip1.intersection(ip2)
            The tamari interval of size 4 induced by relations [(1, 2), (2, 3), (4, 3)]
            sage: ip3 = TamariIntervalPoset(4,[(2,1)])
            sage: ip1.intersection(ip3)
            Traceback (most recent call last):
            ...
            ValueError: This intersection is empty, it does not correspond to an interval-poset.
            sage: ip4 = TamariIntervalPoset(3,[(2,3)])
            sage: ip2.intersection(ip4)
            Traceback (most recent call last):
            ...
            ValueError: Intersections are only possible on interval-posets of the same size.
        """
        if other.size() != self.size():
            raise ValueError, "Intersections are only possible on interval-posets of the same size."
        try:
            return TamariIntervalPoset(self.size(),self._cover_relations + other._cover_relations)
        except ValueError:
            raise ValueError, "This intersection is empty, it does not correspond to an interval-poset."%()

    def initial_forest(self):
        r"""
        Return the initial forest of ``self``, i.e., the interval-poset
        formed with only the increasng relations of ``self``.

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)]).initial_forest()
            The tamari interval of size 4 induced by relations [(1, 2), (2, 4), (3, 4)]
            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.initial_forest() == ip
            True
        """
        return TamariIntervalPoset(self.size(),self.increasing_cover_relations())
        
    def final_forest(self):
        r"""
        Return the final forest of ``self``, i.e., the interval-poset
        formed with only the decreasing relations of ``self``.

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(2,1),(3,2),(3,4),(4,2)]).final_forest()
            The tamari interval of size 4 induced by relations [(4, 2), (3, 2), (2, 1)]
            sage: ip = TamariIntervalPoset(3,[(2,1),(3,1)])
            sage: ip.final_forest() == ip
            True
        """
        return TamariIntervalPoset(self.size(),self.decreasing_cover_relations())

    def lower_binary_tree(self):
        r"""
        Return the lower binary tree of the interval represented by ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.lower_binary_tree()
            [[., .], [[., [., .]], [., .]]]
            sage: TamariIntervalPosets.final_forest(ip.lower_binary_tree()) == ip.final_forest()
            True
            sage: ip == TamariIntervalPosets.from_binary_trees(ip.lower_binary_tree(),ip.upper_binary_tree())
            True
        """
        return self.min_linear_extension().binary_search_tree_shape(left_to_right=False)

    def lower_dyck_word(self):
        r"""
        Return the lower Dyck word of the interval represented by ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.lower_dyck_word()
            [1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0]
            sage: TamariIntervalPosets.final_forest(ip.lower_dyck_word()) == ip.final_forest()
            True
            sage: ip == TamariIntervalPosets.from_dyck_words(ip.lower_dyck_word(),ip.upper_dyck_word())
            True
        """
        return self.lower_binary_tree().to_dyck_word_tamari()

    def upper_binary_tree(self):
        r"""
        Return the upper binary tree of the interval represented by ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.upper_binary_tree()
            [[., .], [., [[., .], [., .]]]]
            sage: TamariIntervalPosets.initial_forest(ip.upper_binary_tree()) == ip.initial_forest()
            True
            sage: ip == TamariIntervalPosets.from_binary_trees(ip.lower_binary_tree(),ip.upper_binary_tree())
            True
        """
        return self.max_linear_extension().binary_search_tree_shape(left_to_right=False)
        
    def upper_dyck_word(self):
        r"""
        Return the upper Dyck word of the interval represented by ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.upper_dyck_word()
            [1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0]
            sage: TamariIntervalPosets.initial_forest(ip.upper_dyck_word()) == ip.initial_forest()
            True
            sage: ip == TamariIntervalPosets.from_dyck_words(ip.lower_dyck_word(),ip.upper_dyck_word())
            True
        """
        return self.upper_binary_tree().to_dyck_word_tamari()

    def sub_poset(self, start, end):
        r"""
        Return the remormalized sub-poset of ``self`` from ``start`` (inclusive)
        to ``end`` (not inclusive).
        
        INPUT:

        - ``start`` -- an integer, the starting vertex (inclusive)
        - ``end`` -- an integer, the ending vertex (not inclusive)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.sub_poset(1,3)
            The tamari interval of size 2 induced by relations [(1, 2)]
            sage: ip.sub_poset(1,4)
            The tamari interval of size 3 induced by relations [(1, 2), (3, 2)]
            sage: ip.sub_poset(1,5)
            The tamari interval of size 4 induced by relations [(1, 2), (4, 3), (3, 2)]
            sage: ip.sub_poset(1,7) == ip
            True
            sage: ip.sub_poset(1,1)
            The tamari interval of size 0 induced by relations []
        """
        if start<1 or start>end or end>self.size()+1:
            raise ValueError, "Invalid starting or ending value, accepted: 1 <= start <= end <= size+1"%()
        if start==end:
            return TamariIntervalPoset(0,[])
        relations = [(i-start+1,j-start+1) for (i,j) in self.increasing_cover_relations() if i>=start and j<end]
        relations.extend([(j-start+1,i-start+1) for (j,i) in self.decreasing_cover_relations() if i>=start and j<end])
        return TamariIntervalPoset(end-start,relations)

    def min_linear_extension(self):
        r"""
        Return the minimal permutation for the right weak order which is 
        a linear extension of ``self``. It corresponds to the minimal 
        permutation of the sylvester class of ``self.lower_binary_tree()``
        and is a 312-avoiding permutation.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.min_linear_extension()
            [1, 2, 4, 3]
            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)])
            sage: ip.min_linear_extension()
            [1, 4, 3, 6, 5, 2]
            sage: ip = TamariIntervalPoset(0,[])
            sage: ip.min_linear_extension()
            []

        """
        def add(perm,i):
            r"""
            Internal recursive method to compute the min linear extension.
            """
            for j in self.decreasing_children(i):
                add(perm,j)
            perm.append(i)
        perm = []
        for i in self.decreasing_roots():
            add(perm,i)
        return Permutation(perm)
        
    def max_linear_extension(self):
        r"""
        Return the maximal permutation for the right weak order which is 
        a linear extension of ``self``. It corresponds to the maximal 
        permutation of the sylvester class of ``self.upper_binary_tree()``
        and is a 132-avoiding permutation.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.max_linear_extension()
            [4, 1, 2, 3]
            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.max_linear_extension()
            [6, 4, 5, 3, 1, 2]
            sage: ip = TamariIntervalPoset(0,[]); ip
            The tamari interval of size 0 induced by relations []
            sage: ip.max_linear_extension()
            []
        """
        def add(perm,i):
            r"""
            Internal recursive method to compute the max linear extension.
            """
            for j in self.increasing_children(i):
                add(perm,j)
            perm.append(i)
        perm = []
        for i in self.increasing_roots():
            add(perm,i)
        return Permutation(perm)


    def linear_extensions(self):
        r"""
        Return an iterator on the permutations which are linear extensions of ``self``.
        They correspond to an interval of the right weak order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(3,[(1,2),(3,2)])
            sage: list(ip.linear_extensions())
            [[3, 1, 2], [1, 3, 2]]
            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: list(ip.linear_extensions())
            [[4, 1, 2, 3], [1, 4, 2, 3], [1, 2, 4, 3]]
        """
        for ext in self._poset.linear_extensions():
            yield(Permutation(ext))
            
    def lower_contained_intervals(self):
        r"""
        If ``self`` represents the interval `[t_1, t_2]`, 
        return an interator on all intervals `[t_1,t]` with `t \leq t2` for the Tamari lattice.
        
        In terms of interval-posets, it corresponds to add all possible relations
        `n` precedes `m` with `n<m`.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.lower_contained_intervals())
            [The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(1, 4), (2, 4), (3, 4), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(2, 3), (3, 4), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(1, 4), (2, 3), (3, 4), (3, 1), (2, 1)]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: len(list(ip.lower_contained_intervals()))
            14

        """
        size = self._size
        yield(self)        
        r""" 
        we try to add links recursively in this order :
        1 -> 2
        2 -> 3
        1 -> 3
        3 -> 4
        2 -> 4
        1 -> 4
        ...
        """
        def add_relations(poset,n, m):
            r"""
            Internal recursive method to generate all possible intervals
            we always have n < m
            """
            if(n<=0):
                #if n<= 0, then we go to the next m
                n=m
                m+=1
            if(m>size):
                #if m>size, it's finished
                return 
            
            if(poset.le(n,m)):
                #there is a already a link n->m, so we go to the next n
                for pos in add_relations(poset, n-1, m):
                    yield pos 
            elif(poset.le(m,n)):
                #there is an inverse link m->n, we know we won't be able 
                #to create a link i-> m with i<n, so we go to the next m
                for pos in add_relations(poset,m,m+1):
                    yield pos
            else:
                # there is no link n->m
                # first option : we don't create the link and go to the next m
                for pos in add_relations(poset,m,m+1):
                    yield pos 
                #second option : we create the link
                poset = TamariIntervalPoset(poset.size(), poset._cover_relations + ((n,m),))
                yield poset 
                #and then, we go to the next n
                for pos in add_relations(poset, n-1,m):
                    yield pos
                
        for inter in add_relations(self, 1, 2):
            yield inter

    def interval_size(self):
        r"""
        Return the size of the interval, i.e., the number of elements 
        (binary trees or Dyck words) in the interval represented by ``self``.
        
        Not to be confused with ``self.size()`` which is the number of 
        vertices.

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).interval_size()
            4
            sage: TamariIntervalPoset(4,[]).interval_size()
            14
            sage: TamariIntervalPoset(4,[(1,2),(2,3),(3,4)]).interval_size()
            1

        """
        return len(list(self.lower_contained_intervals()))
        
    def binary_trees(self):
        r"""
        Return an iterator on all the binary trees of the interval 
        represented by ``self``.

        EXAMPLES::

            sage: list(TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).binary_trees())
            [[., [[., [., .]], .]],
             [[., [., [., .]]], .],
             [., [[[., .], .], .]],
             [[., [[., .], .]], .]]
            sage: set(TamariIntervalPoset(4,[]).binary_trees()) == set(BinaryTrees(4))
            True
        """
        for ip in self.lower_contained_intervals():
            yield ip.upper_binary_tree()
            
    def dyck_words(self):
        r"""
        Return an iterator on all the Dyck words of the interval 
        represented by ``self``.

        EXAMPLES::

            sage: list(TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).dyck_words())
            [[1, 1, 1, 0, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 1, 0, 0],
             [1, 1, 0, 1, 0, 0, 1, 0]]
            sage: set(TamariIntervalPoset(4,[]).dyck_words()) == set(DyckWords(4))
            True
        """
        for ip in self.lower_contained_intervals():
            yield ip.upper_dyck_word()
            
    def maximal_chain_intervals(self):
        r"""
        Return an iterator on the upper contained inverals of one maximal chain of ``self``
        
        If ``self`` represents the interval `[T_1,T_2]`, it returns intervals 
        `[T',T_2]` with `T'` following one maximal chain between `T_1` and `T_2`.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_intervals())
            [The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (4, 1), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (4, 1), (3, 2), (2, 1)]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_intervals())
            [The tamari interval of size 4 induced by relations [],
             The tamari interval of size 4 induced by relations [(2, 1)],
             The tamari interval of size 4 induced by relations [(3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(4, 1), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(4, 1), (3, 2), (2, 1)],
             The tamari interval of size 4 induced by relations [(4, 2), (3, 2), (2, 1)],
             The tamari interval of size 4 induced by relations [(4, 3), (3, 2), (2, 1)]]
        """

        yield(self)
        rel = list(self._cover_relations)
        ti = self
        # we add relations in this order
        # 2 -> 1
        # 3 -> 1
        # 4 -> 1
        # 3 -> 2
        # 4 -> 2
        # 4 -> 3
        for i in xrange(1,self.size()):
            for j in xrange(i+1,self.size()+1):
                if ti.le(j,i):
                    #the relation j->i is already there, we go to the next j
                    continue
                if ti.le(i,j):
                    #there is a relation i->j which forbid any (>j)->i
                    # we go to the next i
                    break 
                # there is no j->i or i->j, so we add j->i
                rel.append((j,i))
                ti = TamariIntervalPoset(self.size(),rel)
                yield ti

    def maximal_chain_binary_trees(self):
        r"""
        Return an iterator on the binary trees forming a maximal chain of ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_binary_trees())
            [[[., [[., .], .]], .], [., [[[., .], .], .]], [., [[., [., .]], .]]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_binary_trees())
            [[[[[., .], .], .], .],
             [[[., [., .]], .], .],
             [[., [[., .], .]], .],
             [., [[[., .], .], .]],
             [., [[., [., .]], .]],
             [., [., [[., .], .]]],
             [., [., [., [., .]]]]]

        """
        for it in self.maximal_chain_intervals():
            yield it.lower_binary_tree()
            
    def maximal_chain_dyck_words(self):
        r"""
        Return an iterator on the Dyck words forming a maximal chain of ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_dyck_words())
            [[1, 1, 0, 1, 0, 0, 1, 0], [1, 1, 0, 1, 0, 1, 0, 0], [1, 1, 1, 0, 0, 1, 0, 0]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_dyck_words())
            [[1, 0, 1, 0, 1, 0, 1, 0],
             [1, 1, 0, 0, 1, 0, 1, 0],
             [1, 1, 0, 1, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 1, 0, 0],
             [1, 1, 1, 0, 1, 0, 0, 0],
             [1, 1, 1, 1, 0, 0, 0, 0]]
        """
        for it in self.maximal_chain_intervals():
            yield it.lower_dyck_word()
    
    def length_of_maximal_chain(self):
        r"""
        Return the length of a maximal chain of ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.length_of_maximal_chain()
            3
            sage: ip = TamariIntervalPoset(4,[])
            sage: ip.length_of_maximal_chain()
            7
            sage: ip = TamariIntervalPoset(3,[])
            sage: ip.length_of_maximal_chain()
            4
        """
        return len(list(self.maximal_chain_intervals()))

# Abstract class to serve as a Factory no instance are created.
class TamariIntervalPosets(UniqueRepresentation, Parent):
    r"""
    Factory for interval-posets.

    INPUT:

    - ``size`` -- (optional) an integer

    OUPUT:

    - the set of all interval-posets (of the given ``size`` if specified)

    EXAMPLES::

        sage: TamariIntervalPosets()
        Interval-posets

        sage: TamariIntervalPosets(2)
        Interval-posets of size 2

    .. NOTE:: this in a factory class whose constructor returns instances of
              subclasses.
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        r"""
        TESTS::

            sage: from sage.combinat.interval_posets import TamariIntervalPosets_all, TamariIntervalPosets_size
            sage: isinstance(TamariIntervalPosets(2), TamariIntervalPosets)
            True
            sage: isinstance(TamariIntervalPosets(), TamariIntervalPosets)
            True
            sage: TamariIntervalPosets(2) is TamariIntervalPosets_size(2)
            True
            sage: TamariIntervalPosets() is TamariIntervalPosets_all()
            True
        """
        if n is None:
            return TamariIntervalPosets_all()
        else:
            if not (isinstance(n, (Integer, int)) and n >= 0):
                raise ValueError("n must be a non negative integer")
            return TamariIntervalPosets_size(Integer(n))

    @staticmethod
    def check_poset(poset):
        r"""
        Checks if the given poset is a interval-poset, that is

        - if its labels are exactly `1,\dots,n`
        - if a<c (as numbers) and a precedes c, then b precedes c for all
          b such that a<b<c
        - if a<c (as numbers) and c precedes a, then b precedes a for all
          b such that a<b<c

        INPUT:

        - ``poset`` -- a finite labeled poset
        
        EXAMPLES::

            sage: p = Poset(([1,2,3],[(1,2),(3,2)]))
            sage: TamariIntervalPosets.check_poset(p)
            True
            sage: p = Poset(([2,3],[(3,2)]))
            sage: TamariIntervalPosets.check_poset(p)
            False
            sage: p = Poset(([1,2,3],[(3,1)]))
            sage: TamariIntervalPosets.check_poset(p)
            False
            sage: p = Poset(([1,2,3],[(1,3)]))
            sage: TamariIntervalPosets.check_poset(p)
            False

        """
        if not set(poset._elements) == set([i for i in xrange(1,poset.cardinality()+1)]):
            return False

        for i in xrange(1,poset.cardinality() +1):
            stop = False
            for j in xrange(i-1,0,-1):
                if(not poset.le(j,i)):
                    stop = True # j does not precede i so no j'<j should
                elif stop:
                    return False
            stop = False
            for j in xrange(i+1, poset.cardinality()+1):
                if(not poset.le(j,i)):
                    stop = True # j does not precede i so no j'>j should
                elif stop:
                    return False
        return True
        
    @staticmethod
    def final_forest(element):
        r"""
        Return the final forest of a binary tree, an interval-poset or a Dyck word. 
        An final forest is an interval-poset corresponding to a final
        interval of the Tamari lattice, i.e., containing only decreasing 
        relations.
        
        It can be constructed from a binary tree by its binary
        search tree labeling  with the rule: `b` precedes
        `a` in the final forest iif `b` is in the right subtree of `a`.

        INPUT:

        - ``element`` -- a binary tree, a Dyck word or an interval-poset

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: TamariIntervalPosets.final_forest(ip)
            The tamari interval of size 4 induced by relations [(1, 2), (2, 3)]

        from binary trees::

            sage: bt = BinaryTree(); bt
            .
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 0 induced by relations []
            sage: bt = BinaryTree([]); bt
            [., .]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 1 induced by relations []
            sage: bt = BinaryTree([[],None]); bt
            [[., .], .]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 2 induced by relations []
            sage: bt = BinaryTree([None,[]]); bt
            [., [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 2 induced by relations [(2, 1)]
            sage: bt = BinaryTree([[],[]]); bt
            [[., .], [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 3 induced by relations [(3, 2)]
            sage: bt = BinaryTree([[None,[[],None]],[]]); bt
            [[., [[., .], .]], [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 5 induced by relations [(5, 4), (3, 1), (2, 1)]

        from Dyck words::

            sage: dw = DyckWord([1,0])
            sage: TamariIntervalPosets.final_forest(dw)
            The tamari interval of size 1 induced by relations []
            sage: dw = DyckWord([1,1,0,1,0,0,1,1,0,0])
            sage: TamariIntervalPosets.final_forest(dw)
            The tamari interval of size 5 induced by relations [(5, 4), (3, 1), (2, 1)]
        """

        if isinstance(element,TamariIntervalPoset):
            return element.initial_forest()
        elif element in DyckWords():
            binary_tree = element.to_binary_tree_tamari()
        elif element in BinaryTrees() or element in LabelledBinaryTrees():
            binary_tree = element 
        else:
            raise ValueError, "Do not know how to construct the initial forest of %s"%(str(element))

        def get_relations(bt, start = 1):
            r"""
            Recursive method to get the binary tree final forest relations
            with only one recursive read of the tree.

            OUTPUT:

            - the indexes of the nodes on the left border of the tree
            - the relations of the final forest (as a list of tuples)
            - the next available index for a node (size of tree +1)
            """
            if not bt:
                return [],[], start # leaf
            roots, relations, index = get_relations(bt[0],start=start)
            rroots, rrelations, rindex = get_relations(bt[1],start=index+1)
            roots.append(index)
            relations.extend(rrelations)
            relations.extend([(j,index) for j in rroots])
            return roots,relations,rindex
        
        roots, relations, index = get_relations(binary_tree)
        return TamariIntervalPoset(index-1,relations)

    @staticmethod
    def initial_forest(element):
        r"""
        Return the inital forest of a binary tree, an interval-poset or a Dyck word. 
        An initial forest is an interval-poset corresponding to an initial
        interval of the Tamari lattice, i.e., containing only increasing 
        relations.
        
        It can be constructed from a binary tree by its binary
        search tree labeling  with the rule: `a` precedes
        `b` in the initial forest iif `a` is in the left subtree of `b`.

        INPUT:

        - ``element`` -- a binary tree, a Dyck word or an interval-poset

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: TamariIntervalPosets.initial_forest(ip)
            The tamari interval of size 4 induced by relations [(1, 2), (2, 3)]

        with binary trees::

            sage: bt = BinaryTree(); bt
            .
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 0 induced by relations []
            sage: bt = BinaryTree([]); bt
            [., .]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 1 induced by relations []
            sage: bt = BinaryTree([[],None]); bt
            [[., .], .]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 2 induced by relations [(1, 2)]
            sage: bt = BinaryTree([None,[]]); bt
            [., [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 2 induced by relations []
            sage: bt = BinaryTree([[],[]]); bt
            [[., .], [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: bt = BinaryTree([[None,[[],None]],[]]); bt
            [[., [[., .], .]], [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 5 induced by relations [(1, 4), (2, 3), (3, 4)]
            
        from Dyck words::

            sage: dw = DyckWord([1,0])
            sage: TamariIntervalPosets.initial_forest(dw)
            The tamari interval of size 1 induced by relations []
            sage: dw = DyckWord([1,1,0,1,0,0,1,1,0,0])
            sage: TamariIntervalPosets.initial_forest(dw)
            The tamari interval of size 5 induced by relations [(1, 4), (2, 3), (3, 4)]
        """
        if isinstance(element,TamariIntervalPoset):
            return element.initial_forest()
        elif element in DyckWords():
            binary_tree = element.to_binary_tree_tamari()
        elif element in BinaryTrees() or element in LabelledBinaryTrees():
            binary_tree = element 
        else:
            raise ValueError, "Do not know how to construct the initial forest of %s"%(str(element))
        
        def get_relations(bt, start = 1):
            r"""
            Recursive method to get the binary tree final forest relations
            with only one recursive read of the tree.

            OUTPUT:

            - the indexes of the nodes on the right border of the tree
            - the relations of the initial forest (as a list of tuples)
            - the next available index for a node (size of tree +1)
            """
            if not bt:
                return [],[], start # leaf
            lroots, lrelations, index = get_relations(bt[0],start=start)
            roots, relations, rindex = get_relations(bt[1],start=index+1)
            roots.append(index)
            relations.extend(lrelations)
            relations.extend([(j,index) for j in lroots])
            return roots,relations,rindex

        roots, relations, index = get_relations(binary_tree)
        return TamariIntervalPoset(index-1,relations)

    @staticmethod
    def from_binary_trees(tree1,tree2):
        r"""
        Return the interval-poset corresponding to the interval [``tree1``,``tree2``]
        in the Tamari lattice. Raise an exception if the two trees are not comparable.

        INPUT:

        - ``tree1`` -- a binary tree
        - ``tree2`` -- a binary tree biger or equal than ``tree1`` for the Tamari lattice.
        
        EXAMPLES::

            sage: tree1 = BinaryTree([[],None])
            sage: tree2 = BinaryTree([None,[]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree2)
            The tamari interval of size 2 induced by relations []
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree1)
            The tamari interval of size 2 induced by relations [(1, 2)]
            sage: TamariIntervalPosets.from_binary_trees(tree2,tree2)
            The tamari interval of size 2 induced by relations [(2, 1)]

            sage: tree1 = BinaryTree([[],[[None,[]],[]]])
            sage: tree2 = BinaryTree([None,[None,[None,[[],[]]]]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree2)
            The tamari interval of size 6 induced by relations [(4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]

            sage: tree3 = BinaryTree([None,[None,[[],[None,[]]]]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree3)
            Traceback (most recent call last):
            ...
            ValueError: The two binary trees are not comparable on the Tamari lattice.
            sage: TamariIntervalPosets.from_binary_trees(tree1,BinaryTree())
            Traceback (most recent call last):
            ...
            ValueError: The two binary trees are not comparable on the Tamari lattice.

        """
        initial_forest = TamariIntervalPosets.initial_forest(tree2)
        final_forest = TamariIntervalPosets.final_forest(tree1)
        try:
            return initial_forest.intersection(final_forest)
        except:
            raise ValueError, "The two binary trees are not comparable on the Tamari lattice."%()

    @staticmethod
    def from_dyck_words(dw1, dw2):
        r"""
        Return the interval-poset corresponding to the interval [``dw1``,``dw2``]
        in the Tamari lattice. Raise an exception if the two Dyck words are not comparable.

        INPUT:

        - ``dw1`` -- a Dyck word
        - ``dw2`` -- a Dyck word biger or equal than ``dw1`` for the Tamari lattice.

        EXAMPLES::

            sage: dw1 = DyckWord([1,0,1,0])
            sage: dw2 = DyckWord([1,1,0,0]) 
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw2)
            The tamari interval of size 2 induced by relations []
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw1)
            The tamari interval of size 2 induced by relations [(1, 2)]
            sage: TamariIntervalPosets.from_dyck_words(dw2,dw2)
            The tamari interval of size 2 induced by relations [(2, 1)]

            sage: dw1 = DyckWord([1,0,1,1,1,0,0,1,1,0,0,0])
            sage: dw2 = DyckWord([1,1,1,1,0,1,1,0,0,0,0,0])
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw2)
            The tamari interval of size 6 induced by relations [(4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]

            sage: dw3 = DyckWord([1,1,1,0,1,1,1,0,0,0,0,0])
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw3)
            Traceback (most recent call last):
            ...
            ValueError: The two Dyck words are not comparable on the Tamari lattice.
            sage: TamariIntervalPosets.from_dyck_words(dw1,DyckWord([1,0]))
            Traceback (most recent call last):
            ...
            ValueError: The two Dyck words are not comparable on the Tamari lattice.

        """
        tree1 = dw1.to_binary_tree_tamari()
        tree2 = dw2.to_binary_tree_tamari()
        try:
            return TamariIntervalPosets.from_binary_trees(tree1,tree2)
        except:
            raise ValueError, "The two Dyck words are not comparable on the Tamari lattice."%()

    def __call__(self, *args, **keywords):
        r"""
        Allows for a poset to be directly transformed into an interval-poset
        It is some kind of coercion but cannot be made throught the cohercion
        system because posets do not have parents.

        EXAMPLES::

            sage: TIP = TamariIntervalPosets()
            sage: p = Poset( ([1,2,3], [(1,2)]))
            sage: TIP(p)
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP(3,[(1,2)])
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: p = Poset(([1,2,3],[(1,3)]))
            sage: TIP(p)
            Traceback (most recent call last):
            ...
            ValueError: This does not satisfy the Tamari interval-poset condition.

        """
        if isinstance(args[0],TamariIntervalPoset):
            return args[0]
        if len(args)==1 and isinstance(args[0],FinitePoset):
            return self.element_class(args[0].cardinality(), args[0].cover_relations())

        return super(TamariIntervalPosets, self).__call__(*args, **keywords)
    
    def le(self,el1,el2):
        r"""
        Poset stucture on interval-poset through interval containment.
        Return whether the interval represented by ``el1`` is contained in
        the interval represented by ``el2``.

        INPUT:

        - ``el1`` -- an interval-poset
        - ``el2`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: TamariIntervalPosets().le(ip1,ip2)
            True
            sage: TamariIntervalPosets().le(ip2,ip1)
            False
        """
        return el2.contains_interval(el1)

    global_options = TamariIntervalPosetOptions


#################################################################
# Enumerated set of all Tamari Interval-posets
#################################################################
class TamariIntervalPosets_all(DisjointUnionEnumeratedSets, TamariIntervalPosets):

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.combinat.interval_posets import TamariIntervalPosets_all
            sage: S = TamariIntervalPosets_all()
            sage: S.cardinality()
            +Infinity

            sage: it = iter(S)
            sage: [it.next() for i in xrange(5)]
            [The tamari interval of size 0 induced by relations [],
             The tamari interval of size 1 induced by relations [],
             The tamari interval of size 2 induced by relations [],
             The tamari interval of size 2 induced by relations [(2, 1)],
             The tamari interval of size 2 induced by relations [(1, 2)]]
            sage: it.next().parent()
            Interval-posets
            sage: S(0,[])
            The tamari interval of size 0 induced by relations []

            sage: S is TamariIntervalPosets_all()
            True
            sage: TestSuite(S).run()
            """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), TamariIntervalPosets_size),
            facade=True, keepkey = False, category=(Posets(),EnumeratedSets()))

    def _repr_(self):
        r"""
        TEST::

            sage: TamariIntervalPosets() # indirect doctest
            Interval-posets

        """
        return "Interval-posets"
        
    def _element_constructor_(self, size, relations):
        r"""
            EXAMPLES::

                sage: TIP = TamariIntervalPosets()
                sage: TIP(3,[(1,2)]) # indirect doctest
                The tamari interval of size 3 induced by relations [(1, 2)]
        """
        return self.element_class(size,relations)
        
    def __contains__(self, x):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets()
            sage: 1 in S
            False
            sage: S(0,[]) in S
            True
        """
        return isinstance(x, self.element_class)

    Element = TamariIntervalPoset

#################################################################
# Enumerated set of Tamari interval-posets of a given size
#################################################################
class TamariIntervalPosets_size(TamariIntervalPosets):
    r"""
    The enumerated sets of interval-posets of a given size

    TESTS::

        sage: from sage.combinat.interval_posets import TamariIntervalPosets_size
        sage: for i in xrange(6): TestSuite(TamariIntervalPosets_size(i)).run()
    """
    def __init__(self, size):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: S == loads(dumps(S))
            True

            sage: S is TamariIntervalPosets(3)
            True
        """
        # there is a natural order on interval-posets throught inclusions
        # that is why we use the FinitePosets category
        super(TamariIntervalPosets_size, self).__init__(category = (FinitePosets(),FiniteEnumeratedSets()))

        self._size = size
        
    def _repr_(self):
        r"""
        TESTS::

            sage: TamariIntervalPosets(3)
            Interval-posets of size 3
        """
        return "Interval-posets of size %s"%(self._size)

    def __contains__(self, x):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: 1 in S
            False
            sage: S([]) in S
            True
        """
        return isinstance(x, self.element_class) and x.size() == self._size
        
    def cardinality(self):
        r"""
        The cardinality of ``self``

        The formula was given in [CHA]_ `\frac{2(4n+1)!}{(n+1)!(3n+2)!}`

        REFERENCES:

        .. [CHA] Sur le nombre d'intervalles dans les treillis de Tamari, F. Chapoton, 2008

        EXAMPLES::

            sage: TamariIntervalPosets(2).cardinality()
            3
            sage: TamariIntervalPosets(3).cardinality()
            13
            sage: TamariIntervalPosets(4).cardinality()
            68
            sage: TamariIntervalPosets(5).cardinality()
            399
        """
        from sage.functions.other import factorial
        n = self._size
        return Integer(2*factorial(4*n+1)/(factorial(n+1)*factorial(3*n+2)))

    def __iter__(self):
        r"""
        Recursive generation: we iterate through all interval-posets of
        `size -1` and add all possible relations to the last vertex.

        TESTS::

            sage: TIP1 = TamariIntervalPosets(1)
            sage: list(TIP1)
            [The tamari interval of size 1 induced by relations []]
            sage: TIP2 = TamariIntervalPosets(2)
            sage: list(TIP2)
            [The tamari interval of size 2 induced by relations [],
             The tamari interval of size 2 induced by relations [(2, 1)],
             The tamari interval of size 2 induced by relations [(1, 2)]]
            sage: TIP3 = TamariIntervalPosets(3)
            sage: list(TIP3)
            [The tamari interval of size 3 induced by relations [],
             The tamari interval of size 3 induced by relations [(3, 2)],
             The tamari interval of size 3 induced by relations [(2, 3)],
             The tamari interval of size 3 induced by relations [(1, 3), (2, 3)],
             The tamari interval of size 3 induced by relations [(2, 1)],
             The tamari interval of size 3 induced by relations [(3, 2), (2, 1)],
             The tamari interval of size 3 induced by relations [(3, 1), (2, 1)],
             The tamari interval of size 3 induced by relations [(2, 3), (2, 1)],
             The tamari interval of size 3 induced by relations [(2, 3), (3, 1), (2, 1)],
             The tamari interval of size 3 induced by relations [(1, 3), (2, 3), (2, 1)],
             The tamari interval of size 3 induced by relations [(1, 2)],
             The tamari interval of size 3 induced by relations [(1, 2), (3, 2)],
             The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]]
            sage: all([len(list(TamariIntervalPosets(i)))==TamariIntervalPosets(i).cardinality() for i in xrange(6)])
            True

        """
        n = self._size
        if(n==0):
            yield TamariIntervalPoset(0,[])
        elif(n==1):
            yield TamariIntervalPoset(1,[])
        else:
            for tip in TamariIntervalPosets(n-1):
                new_tip = TamariIntervalPoset(n,tip._cover_relations)
                yield new_tip # we have added an extra vertex but no relations
                
                # adding a decreasing relation n>>m2 with m2<n
                for m2 in xrange(n-1,0,-1):
                    if new_tip.le(n-1,m2):
                        yield TamariIntervalPoset(n,new_tip._cover_relations + ((n,m2),))

                for m in xrange(n-1,0,-1):

                    # addind an increasing relation m>>n
                    if not new_tip.le(m,n):
                        new_tip = TamariIntervalPoset(n, new_tip._cover_relations + ((m,n),))
                        yield new_tip
                    else:
                        continue

                    # adding a decreasing relation n>>m2 with m2<m
                    for m2 in xrange(m-1,0,-1):
                        if new_tip.le(n-1,m2):
                            yield TamariIntervalPoset(n,new_tip._cover_relations + ((n,m2),))

    @lazy_attribute
    def _parent_for(self):
        r"""
        The parent of the element generated by ``self``

        TESTS::

            sage: TIP3 = TamariIntervalPosets(3)
            sage: TIP3._parent_for
            Interval-posets
        """
        return TamariIntervalPosets_all()

    @lazy_attribute
    def element_class(self):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: S.element_class
            <class 'sage.combinat.interval_posets.TamariIntervalPosets_all_with_category.element_class'>
            sage: S.first().__class__ == TamariIntervalPosets().first().__class__
            True

        """
        return self._parent_for.element_class

    def _element_constructor_(self, relations):
        r"""
        EXAMPLES::

            sage: TIP3 = TamariIntervalPosets(3)
            sage: TIP3([(1,2)]) # indirect doctest
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP3([(3,4)])
            Traceback (most recent call last):
            ...
            ValueError: The relations do not correspond to the size of the poset.

        """
        return self.element_class(self._size,relations) 

