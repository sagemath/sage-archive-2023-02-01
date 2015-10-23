r"""
The Parallelogram Polyominoes
=============================

The goal of this module is to give some tools to manipulate the
parallelogram polyominoes.
"""
#******************************************************************************
#  Copyright (C) 2014 Adrien Boussicault (boussica@labri.fr),
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.list_clone import ClonableList
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.set_factories import (
    SetFactory, ParentWithSetFactory, TopMostParentPolicy
)
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.global_options import GlobalOptions
from sage.sets.set import Set
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.disjoint_union_enumerated_sets \
    import DisjointUnionEnumeratedSets
from sage.rings.integer import Integer
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from copy import deepcopy
from sage.matrix.constructor import matrix
from sage.combinat.combinat import catalan_number
from sage.combinat.combinatorial_map import combinatorial_map
from sage.functions.trig import cos, sin
from sage.functions.other import sqrt

default_tikz_options = dict(
    scale=1, line_size=1, point_size=3.5,
    color_line='black', color_point='black',
    color_bounce_0='red', color_bounce_1='blue',
    translation=[0, 0], rotation=0,
    mirror=None,
)

ParallelogramPolyominoesOptions = GlobalOptions(
    name='Parallelogram Polyominoes',
    doc=r"""
    """,
    end_doc=r"""
    """,
    tikz_options=dict(
        default=default_tikz_options,
        description='the tikz options',
        checker=lambda x: Set(x.keys()).issubset(
            Set(
                [
                    'scale', 'line_size', 'point_size',
                    'color_line', 'color_point', 'translation', 'mirror',
                    'rotation', 'color_bounce_0', 'color_bounce_1',
                ]
            )
        )
    ),
    drawing_components=dict(
        default=dict(diagram=True),
        description='Different tree-like tableaux components to draw',
        checker=lambda x: Set(x.keys()).issubset(
            Set(['diagram', 'tree', 'bounce_0', 'bounce_1', ])
        )
    ),
    display=dict(
        default="list",
        values=dict(
            list='displayed as list',
            drawing='as a drawing',
        )
    ),
    latex=dict(
        default="drawing",
        values=dict(
            list='displayed as list',
            drawing='as a drawing',
        )
    )
)


class _drawing_tool:
    def __init__(self, options, XY=lambda v: v):
        self._XY = lambda v: XY([float(v[0]), float(v[1])])
        self._translation = options['translation']
        self._mirror = options['mirror']
        self._rotation = options['rotation']
        self._color_line = options['color_line']
        self._line_size = options['line_size']
        self._point_size = options['point_size']
        self._color_point = options['color_point']

    def XY(self, v):
        """
        TODO
        """
        def translate(pos, v):
            return [pos[0]+v[0], pos[1]+v[1]]

        def rotate(pos, angle):
            [x, y] = pos
            return [x*cos(angle) - y*sin(angle), x*sin(angle) + y*cos(angle)]

        def mirror(pos,  axe):
            if axe is None:
                return pos
            n = float(sqrt(axe[0]**2 + axe[1]**2))
            axe[0] = float(axe[0]/n)
            axe[1] = float(axe[1]/n)
            sp = (pos[0]*axe[0] + pos[1]*axe[1])
            sn = (- pos[0]*axe[1] + pos[1]*axe[0])
            return [
                sp*axe[0] + sn*axe[1],
                sp*axe[1] - sn*axe[0]
            ]
        return rotate(
            mirror(
                translate(self._XY(v), self._translation),
                self._mirror
            ), self._rotation
        )

    def draw_line(self, v1, v2, color=None, size=None):
        """
        TODO
        """
        if color is None:
            color = self._color_line
        if size is None:
            size = self._line_size
        [x1, y1] = self.XY(v1)
        [x2, y2] = self.XY(v2)
        return "\n  \\draw[color=%s, line width=%s] (%f, %f) -- (%f, %f);" % (
            color, size, float(x1), float(y1), float(x2), float(y2)
        )

    def draw_polyline(self, list_of_vertices, color=None, size=None):
        """
        TODO
        """
        res = ""
        for i in range(len(list_of_vertices)-1):
            res += self.draw_line(
                list_of_vertices[i], list_of_vertices[i+1], color, size
            )
        return res

    def draw_point(self, p1, color=None, size=None):
        """
        TODO
        """
        if color is None:
            color = self._color_point
        if size is None:
            size = self._point_size
        [x1, y1] = self.XY(p1)
        return "\n  \\filldraw[color=%s] (%f, %f) circle (%spt);" % (
            color, float(x1), float(y1), size
        )


class ParallelogramPolyomino(ClonableList):
    r"""
    The class of Parallelogram Polyominoes
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        """
        return cls._auto_parent._element_constructor_(*args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        r"""
        """
        return ParallelogramPolyominoes()

    def check(self):
        r"""
        This method raises an error if the internal data of the class does not
        represent a parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp = ParallelogramPolyomino([[1], [1]])

            sage: pp = ParallelogramPolyomino([[1, 0], [0, 1]])
            Traceback (most recent call last):
            ...
            ValueError: Lower and upper path are crossing.

            sage: pp = ParallelogramPolyomino([[1], [0, 1]])
            Traceback (most recent call last):
            ...
            ValueError: Lower upper paht have different size (2 != 1).

            sage: pp = ParallelogramPolyomino([[1], [0]])
            Traceback (most recent call last):
            ...
            ValueError: The two paths don't join together at the end.

            sage: pp = ParallelogramPolyomino([[0], [1]])
            Traceback (most recent call last):
            ...
            ValueError: The two paths don't join together at the end.

            sage: pp = ParallelogramPolyomino([[0], [0]])
            Traceback (most recent call last):
            ...
            ValueError: A Parallelogam Polyomino can have the path [[0], [0]].

            sage: pp = ParallelogramPolyomino([[], [0]])
            Traceback (most recent call last):
            ...
            ValueError: A Parallelogam Polyomino can have lower or upper path equals to [].

            sage: pp = ParallelogramPolyomino([[0], []])
            Traceback (most recent call last):
            ...
            ValueError: A Parallelogam Polyomino can have lower or upper path equals to [].

            sage: pp = ParallelogramPolyomino([[], []])
            Traceback (most recent call last):
            ...
            ValueError: A Parallelogam Polyomino can have lower or upper path equals to [].
        """
        lower_path = self.lower_path()
        upper_path = self.upper_path()
        if lower_path == [0] and upper_path == [0]:
            raise ValueError(
                "A Parallelogam Polyomino can have the path [[0], [0]]."
            )
        if lower_path == [] or upper_path == []:
            raise ValueError(
                "A Parallelogam Polyomino can have lower or upper path equals "
                "to []."
            )
        if len(upper_path) != len(lower_path):
            raise ValueError(
                "Lower upper paht have different size (%s != %s)." % (
                    len(upper_path), len(lower_path)
                )
            )
        p_up = [0, 0]
        p_down = [0, 0]
        for i in range(len(upper_path)-1):
            p_up[1-upper_path[i]] += 1
            p_down[1-lower_path[i]] += 1
            if(p_up[0] <= p_down[0] or p_down[1] <= p_up[1]):
                raise ValueError("Lower and upper path are crossing.")
        p_up[1-upper_path[-1]] += 1
        p_down[1-lower_path[-1]] += 1
        if(p_up[0] != p_down[0] or p_up[1] != p_down[1]):
            raise ValueError("The two paths don't join together at the end.")

    def __hash__(self):
        r"""
        Return the hash code of the parallelogram polyomino

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: hash(pp) == hash((
            ....:     (0, 0, 0, 1, 0, 1, 0, 1, 1), (1, 0, 1, 1, 0, 0, 1, 0, 0)
            ....: ))
            True

            sage: PPS = ParallelogramPolyominoes(7)
            sage: D = { PPS[0]: True, PPS[1]: True }
            sage: D[PPS[0]] = False
            sage: import pprint
            sage: pp = pprint.PrettyPrinter()
            sage: pp.pprint(D)
            {[[0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0]]: False,
            [[0, 0, 0, 0, 0, 0, 1, 1], [1, 0, 0, 0, 0, 0, 1, 0]]: True}
        """
        return hash(tuple(map(tuple, list(self))))

    def __copy__(self):
        r"""
        Copy a parallelogram Polyomino

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp1 = copy(pp)
            sage: pp1 is pp
            False
            sage: pp1 == pp
            True
            sage: pp1
            [[0, 0, 0, 1, 0, 1, 0, 1, 1], [1, 0, 1, 1, 0, 0, 1, 0, 0]]
        """
        return ParallelogramPolyomino([self.lower_path(), self.upper_path()])

    def __init__(self, parent, value, check=True):
        r"""
        Construct a parallelogram polyomino.

        The input is a pair of lower path and upper path.

        The lower and upper paths of the empty parallelogram polyomino are
        [1] and [1].

        EXAMPLES::

            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino([lower_path, upper_path])
            sage: pp
            [[0, 0, 1, 0, 1, 1], [1, 1, 0, 1, 0, 0]]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp
            [[0, 1], [1, 0]]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp
            [[1], [1]]
        """
        ClonableList.__init__(self, parent, value)
        if check:
            if not isinstance(value, (list, tuple)):
                raise ValueError(
                    "Value %s must be a list or a tuple." % (value)
                )
            self.check()
        self._options = None

    def _to_dyck_delest_viennot(self):
        r"""
        Convert to a Dyck word using the Delest-Viennot bijection.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]
            ....: )
            sage: pp._to_dyck_delest_viennot()
            [1, 1, 0, 1, 1, 0, 1, 0, 0, 0]
        """
        from sage.combinat.dyck_word import DyckWord
        size = self.size()
        dyck = []
        upper_path = self.upper_path()
        lower_path = self.lower_path()
        dyck.append(1 - lower_path[0])
        for i in range(1, size):
            dyck.append(upper_path[i])
            dyck.append(1 - lower_path[i])
        dyck.append(upper_path[size])
        return DyckWord(dyck)

    @combinatorial_map(name="To Dyck word")
    def to_dyck_word(self, bijection=None):
        r"""
        Convert to a Dyck word.

        TODO: verifier:
        ref. Delest, M.-P. and Viennot, G.
        "Algebraic Languages and Polyominoes [sic] Enumeration."
        Theoret. Comput. Sci. 34, 169-206, 1984.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]
            ....: )
            sage: pp.to_dyck_word()
            [1, 1, 0, 1, 1, 0, 1, 0, 0, 0]
            sage: pp.to_dyck_word(bijection='Delest-Viennot')
            [1, 1, 0, 1, 1, 0, 1, 0, 0, 0]
        """
        if bijection is None or bijection == 'Delest-Viennot':
            return self._to_dyck_delest_viennot()

    @staticmethod
    def _from_dyck_word_delest_viennot(dyck):
        r"""
        Convert dyck word to parallelogram polyomino using the Delest Viennot
        bijection.

        TODO: verifier:
        ref. Delest, M.-P. and Viennot, G.
        "Algebraic Languages and Polyominoes [sic] Enumeration."
        Theoret. Comput. Sci. 34, 169-206, 1984.

        INPUTS:

        - `dyck` -- a dyck word

        OUTPUT:

        A parallelogram polyomino.

        EXAMPLES::

            sage: dyck = DyckWord( [1, 1, 0, 1, 1, 0, 1, 0, 0, 0] )
            sage: pp = ParallelogramPolyomino._from_dyck_word_delest_viennot(
            ....:     dyck
            ....: )
            sage: pp
            [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]
        """
        l = [1] + list(dyck) + [0]
        word_up = []
        word_down = []
        for i in range(0, len(l), 2):
            word_up.append(l[i])
            word_down.append(1 - l[i+1])
        return ParallelogramPolyomino([word_down, word_up])

    @staticmethod
    def from_dyck_word(dyck, bijection=None):
        r"""
        Convert dyck word to parallelogram polyomino.

        INPUTS:

        - `dyck` -- a dyck word
        - `bijection` -- (default: 'Delest-Viennot') the bijection to use.
          Expected values : 'Delest-Viennot'.

        OUTPUT:

        A parallelogram polyomino.

        EXAMPLES::

            sage: dyck = DyckWord( [1, 1, 0, 1, 1, 0, 1, 0, 0, 0] )
            sage: pp = ParallelogramPolyomino.from_dyck_word( dyck )
            sage: pp
            [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]
            sage: pp = ParallelogramPolyomino.from_dyck_word(
            ....:     dyck, bijection='Delest-Viennot'
            ....: )
            sage: pp
            [[0, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0]]
        """
        if bijection is None or bijection == 'Delest-Viennot':
            return ParallelogramPolyomino._from_dyck_word_delest_viennot(dyck)

    def _to_binary_tree_Aval_Boussicault(self, position=[0, 0]):
        r"""
        Convert to a binary tree using the Aval-Boussicault algorithm.

        Ref.:
        J.C Aval, A. Boussicault, M. Bouvel, M. Silimbani,
        "Combinatorics of non-ambiguous trees",
        :arxiv:`1305.3716`

        INPUT:

        - ``bijection`` -- ``None`` (default) The name of bijection to use for
          the convertion. The possible value are, 'Aval-Boussicault'.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp._to_binary_tree_Aval_Boussicault()
            [[., [[., .], [[., [., .]], .]]], [[., .], .]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp._to_binary_tree_Aval_Boussicault()
            [., .]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp._to_binary_tree_Aval_Boussicault()
            .
        """
        from sage.combinat.binary_tree import BinaryTree
        if self.size() == 0:
            return BinaryTree()
        result = [BinaryTree(), BinaryTree()]
        right_son = list(position)
        left_son = list(position)
        w = right_son[1] + 1
        h = left_son[0] + 1
        while w < self.width():
            if self[right_son[0]][w] == 1:
                if self[right_son[0]-1][w] == 0:
                    right_son[1] = w
                    result[1] = self._to_binary_tree_Aval_Boussicault(
                        right_son
                    )
                    break
            w += 1
        while h < self.height():
            if self[h][left_son[1]] == 1:
                if self[h][left_son[1]-1] == 0:
                    left_son[0] = h
                    result[0] = self._to_binary_tree_Aval_Boussicault(left_son)
                    break
            h += 1
        return BinaryTree(result)

    def get_sub_binary_tree(self, position):
        r"""
        TODO
        """
        from sage.combinat.binary_tree import BinaryTree
        if self[position[0]][position[1]] == 0:
            return BinaryTree()
        return self._to_tree_Aval_Boussicault(position)

    @combinatorial_map(name="To binary tree")
    def to_binary_tree(self, bijection=None):
        r"""
        Convert to a binary tree

        INPUT:

        - ``bijection`` -- ``None`` (default) The name of bijection to use for
          the convertion. The possible value are, 'Aval-Boussicault'.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp.to_binary_tree()
            [[., [[., .], [[., [., .]], .]]], [[., .], .]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.to_binary_tree()
            [., .]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.to_binary_tree()
            .
        """
        if bijection is None or bijection == 'Aval-Boussicault':
            return self._to_binary_tree_Aval_Boussicault([0, 0])

    def _to_ordered_tree_via_dyck(self):
        return self._to_dyck_delest_viennot().to_ordered_tree()

    def _to_ordered_tree_Bou_Socci(self):
        r"""
        Return the ordered tree using the Boussicault-Socci bijection.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp.to_ordered_tree()
            [[[[[]], [[[]]]]], [[]]]
            sage: pp.to_ordered_tree(bijection='Boussicault-Socci')
            [[[[[]], [[[]]]]], [[]]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.to_ordered_tree()
            [[]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.to_ordered_tree()
            []
        """
        from sage.combinat.ordered_tree import OrderedTree
        from sage.combinat.binary_tree import BinaryTree

        def make_tree(b_tree, d):
            if b_tree == BinaryTree():
                return OrderedTree([])
            res = []
            res.append(make_tree(b_tree[1-d], 1-d))
            res += make_tree(b_tree[d], d)
            return OrderedTree(res)
        return make_tree(
            self.to_binary_tree(bijection='Aval-Boussicault'), 1
        )

    @combinatorial_map(name="To ordered tree")
    def to_ordered_tree(self, bijection=None):
        r"""
        Return an ordered tree from the parallelogram polyomino

        INPUT:

        - ``bijection`` -- ``None`` (default) The name of bijection to use for
          the convertion. The possible value are, 'Boussicault-Socci',
          'via dyck', 'via dyck and Delest-Viennot'. The default bijection is
          'Boussicault-Socci'.

          TODO: Faire de la biblio.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 1, 0, 1, 1, 0, 0, 0, 1, 0]
            ....:     ]
            ....: )
            sage: pp.to_ordered_tree()
            [[[[[]], [[[]]]]], [[]]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.to_ordered_tree()
            [[]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.to_ordered_tree()
            []
        """
        if bijection is None or bijection == 'Boussicault-Socci':
            return self._to_ordered_tree_Bou_Socci()
        if bijection == 'via dyck and Delest-Viennot':
            return self._to_ordered_tree_via_dyck()

    def get_options(self):
        r"""
        Return all the options of the object.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.get_options()
            options for Parallelogram Polyominoes
        """
        if self._options is None:
            return self.parent().get_options()
        return self._options

    def set_options(self, *get_value, **set_value):
        r"""
        Set new options to the object.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
            ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: pp
            [[0, 0, 0, 0, 1, 0, 1, 0, 1], [1, 0, 0, 0, 1, 1, 0, 0, 0]]
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 0 0]
            [1 0 0]
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 0 1]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: view(PP) # not tested
            sage: pp.set_options(
            ....:     drawing_components=dict(
            ....:         diagram = True,
            ....:         bounce_0 = True,
            ....:         bounce_1 = True,
            ....:     )
            ....: )
            sage: view(PP) # not tested
        """
        if self._options is None:
            self._options = deepcopy(self.get_options())
        self._options(*get_value, **set_value)

    def upper_path(self):
        r"""
        Get the upper path of the parallelogram polyomino.

        EXAMPLES::

            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino([lower_path, upper_path])
            sage: pp.upper_path()
            [1, 1, 0, 1, 0, 0]
        """
        return list(ClonableList.__getitem__(self, 1))

    def lower_path(self):
        r"""
        Get the lower path of the parallelogram polyomino.

        EXAMPLES::

            sage: lower_path = [0, 0, 1, 0, 1, 1]
            sage: upper_path = [1, 1, 0, 1, 0, 0]
            sage: pp = ParallelogramPolyomino([lower_path, upper_path])
            sage: pp.lower_path()
            [0, 0, 1, 0, 1, 1]
        """
        return list(ClonableList.__getitem__(self, 0))

    def _path_heights(self, word, up):
        res = []
        h = 0
        for e in word:
            if e == up:
                h += 1
            else:
                res.append(h)
        return res

    def upper_heights(self):
        return self._path_heights(self.upper_path(), 0)

    def lower_heights(self):
        return self._path_heights(self.lower_path(), 0)

    def upper_widths(self):
        return self._path_heights(self.upper_path(), 1)

    def lower_widths(self):
        return self._path_heights(self.lower_path(), 1)

    def widths(self):
        r"""
        This method return a list of the widths of the parallelogram
        polyomino: the parallelogram polyomino is splitted row by row and the
        method returns the list containing the sizes of the rows.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.widths()
            [1, 3, 3, 3, 2]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.widths()
            [1]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.widths()
            []
        """
        widths = []
        uw = self.upper_widths()
        lw = self.lower_widths()
        for i in range(len(lw)):
            widths.append(uw[i] - lw[i])
        return widths

    def degree_convexity(self):
        r"""
        Return the degree convexity of a parallelogram polyomino.

        A convex polyomino is said to be k-convex if every pair of its cells
        can be connected by a monotone path (path with south and east steps)
        with at most k changes of direction.
        The degree of convexity of a convex polyomino P is the smallest integer
        k such that P is k-convex.

        If the parallelogram polyomino is empty, the function return -1.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.degree_convexity()
            3

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.degree_convexity()
            0

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.degree_convexity()
            -1
        """
        l0 = len(self.bounce_path(direction=0))
        l1 = len(self.bounce_path(direction=1))
        return min(l0, l1) - 1

    def is_flat(self):
        r"""
        Return true if the two bounce path join together in the rightmost cell
        of the bottom row of P.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.is_flat()
            False

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.is_flat()
            True

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.is_flat()
            True
        """

        l0 = len(self.bounce_path(direction=0))
        l1 = len(self.bounce_path(direction=1))
        return l0 == l1

    def is_k_directed(self, k):
        r"""
        Return true if the Polyomino Parallelogram is k-directed.

        A convex polyomino is said to be k-convex if every pair of its cells
        can be connected by a monotone path (path with south and east steps)
        with at most k changes of direction.
        The degree of convexity of a convex polyomino P is the smallest integer
        k such that P is k-convex.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.is_k_directed(3)
            True
            sage: pp.is_k_directed(4)
            True
            sage: pp.is_k_directed(5)
            True
            sage: pp.is_k_directed(0)
            False
            sage: pp.is_k_directed(1)
            False
            sage: pp.is_k_directed(2)
            False

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.is_k_directed(0)
            True
            sage: pp.is_k_directed(1)
            True

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.is_k_directed(0)
            True
            sage: pp.is_k_directed(1)
            True
        """
        return self.degree_convexity() <= k

    def heights(self):
        r"""
        This method return a list of heights of the parallelogram
        polyomino: the parallelogram polyomino is splitted column by column and
        the method returns the list containing the sizes of the columns.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:         [1, 0, 1, 1, 0, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.heights()
            [3, 3, 4, 2]

            sage: pp = ParallelogramPolyomino([[0, 1], [1, 0]])
            sage: pp.heights()
            [1]

            sage: pp = ParallelogramPolyomino([[1], [1]])
            sage: pp.heights()
            [0]
        """
        heights = []
        uh = self.upper_heights()
        lh = self.lower_heights()
        for i in range(len(uh)):
            heights.append(lh[i] - uh[i])
        return heights

    def width(self):
        r"""
        Return the width of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 1, 0, 0, 1, 1, 0, 1, 1, 1],
            ....:         [1, 1, 1, 0, 1, 0, 0, 1, 1, 0]
            ....:     ]
            ....: )
            sage: pp.width()
            6

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.width()
            1

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.width()
            1
        """
        if self.size() == 0:
            return 1
        return self.upper_widths()[-1]

    def height(self):
        r"""
        Return the height of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 1, 0, 0, 1, 1, 0, 1, 1, 1],
            ....:         [1, 1, 1, 0, 1, 0, 0, 1, 1, 0]
            ....:     ]
            ....: )
            sage: pp.height()
            4

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.height()
            1

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.height()
            0
        """
        if self.size() == 0:
            return 0
        return self.lower_heights()[-1]

    @cached_method
    def get_array(self):
        r"""
        Return an array of 0s and 1s such that the 1s represent the boxes of
        the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 0, 0, 1, 0, 1, 0, 1],
            ....:         [1, 0, 0, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: matrix(pp.get_array())
            [1 0 0]
            [1 0 0]
            [1 0 0]
            [1 1 1]
            [0 1 1]
            [0 0 1]

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.get_array()
            [[1]]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.get_array()
            []
        """
        width = self.width()
        height = self.height()
        lower_widths = self.lower_widths()
        widths = self.widths()

        def val(w, h):
            if w >= len(widths) or w < 0:
                return 0
            if lower_widths[w] <= h and h < lower_widths[w]+widths[w]:
                return 1
            return 0
        return [[val(h, w) for w in range(width)] for h in range(height)]

    class _polyomino_row:
        def __init__(self, polyomino, row):
            self.polyomino = polyomino
            self.row = row

        def __getitem__(self, column):
            if(
                self.is_inside()
                and 0 <= column and column < self.polyomino.width()
            ):
                return self.polyomino.get_array()[self.row][column]
            return 0

        def is_inside(self):
            return 0 <= self.row and self.row < self.polyomino.height()

        def is_outside(self):
            return not self.is_inside()

        def __repr__(self):
            if self.is_outside():
                return "The (outside) row %s of the parallelogram" % (self.row)
            else:
                return "The row %s of the parallelogram polyomino" % (self.row)

    def __getitem__(self, row):
        r"""
        Return the row of the parallelogram polyomino.

        The index of the row can be out of range of the height of
        the parallelogram polyomino.
        In that case, the row returned is outside the parallelogram
        polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp[0].is_inside()
            True
            sage: pp[3].is_outside()
            True
            sage: pp[-1].is_outside()
            True
            sage: pp[0][0]
            1
            sage: pp[[0, 1]]
            1
            sage: pp[2][0]
            0
            sage: pp[-1][0]
            0
            sage: pp[[4, 4]]
            0
        """
        if isinstance(row, list):
            return self._polyomino_row(self, row[0])[row[1]]
        return self._polyomino_row(self, row)

    def bounce_path(self, direction=1):
        r"""
        Return the bounce path of the parallelogram polyomino.

        The bounce path is a path with two steps (1, 0) and (0, 1).

        If 'direction' is 1 (resp. 0), the bounce path is the path
        starting at position position (h=1, w=0) (resp. (h=0, w=1)) with
        initial direction, the vector (0, 1) (resp. (1, 0)), and turning
        each time the path crosses the perimeter of the parallelogram
        polyomino.

        The path is coded by a list of integers. Each integer represents
        the size of the path between two turnings.

        You can visualize the two bounce paths by using the following
        commands:

        EXAMPLES::

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: PP.bounce_path(direction=1)
            [2, 2, 1]
            sage: PP.bounce_path(direction=0)
            [2, 1, 1, 1]

            sage: PP = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 1, 1, 0, 0, 1, 1],
            ....:         [1, 1, 1, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: PP.bounce_path(direction=1)
            [3, 1, 2, 2]
            sage: PP.bounce_path(direction=0)
            [2, 4, 2]

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: PP.set_options(
            ....:     drawing_components=dict(
            ....:         diagram = True
            ....:         , bounce_0 = True
            ....:         , bounce_1 = True
            ....:     )
            ....: )
            sage: view(PP) # not tested

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: PP.bounce_path(direction=1)
            [1]
            sage: PP.bounce_path(direction=0)
            [1]

            sage: PP = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: PP.bounce_path(direction=1)
            []
            sage: PP.bounce_path(direction=0)
            []

        TESTS::

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: PP.bounce_path(direction=1)
            [2, 2, 1]
            sage: PP.bounce_path(direction=0)
            [2, 1, 1, 1]
        """
        result = []
        pos = [0, 0]
        pos[direction] -= 1
        old = list(pos)
        ne = list(pos)
        ne[direction] += 1
        while self[ne] == 1:
            pos[direction] += 1
            while self[pos] == 1:
                pos[direction] += 1
            pos[direction] -= 1
            result.append(pos[direction]-old[direction])
            direction = 1 - direction
            old[0], old[1] = pos
            ne[0], ne[1] = pos
            ne[direction] += 1
        return result

    def bounce(self, direction=1):
        r"""
        Return the bounce of the parallelogram polyomino.

        Les p be the bounce path of the parallelogram
        polyomino. (p=self.bounce_path())
        The bounce is defined by:
        sum([(1+ floor(i/2))*p[i] for i in range(len(p))])

        EXAMPLES::

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: PP.bounce(direction=1)
            6
            sage: PP.bounce(direction=0)
            7

            sage: PP = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 1, 1, 0, 0, 1, 1],
            ....:         [1, 1, 1, 0, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: PP.bounce(direction=1)
            12
            sage: PP.bounce(direction=0)
            10

            sage: PP = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: PP.bounce(direction=1)
            1
            sage: PP.bounce(direction=0)
            1

            sage: PP = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: PP.bounce(direction=1)
            0
            sage: PP.bounce(direction=0)
            0
        """
        result = 0
        path = self.bounce_path(direction)
        for i in range(len(path)):
            result += (1+(int(i/2)))*path[i]
        return result

    def area(self):
        r"""
        Returns the area of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1],
            ....:         [1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0]
            ....:     ]
            ....: )
            sage: pp.area()
            13

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.area()
            1

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.area()
            0
        """
        res = 0
        for h in self.heights():
            res += h
        return res

    def _repr_(self):
        r"""
        Return a string representation of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp
            [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            sage: pp.set_options(display='drawing')
            sage: pp
            [1 1 0]
            [1 1 0]
            [0 1 1]
        """
        return self.get_options().dispatch(self, '_repr_', 'display')

    def _repr_list(self):
        r"""
        Return a string representation with list style.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp._repr_list()
            '[[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]'
        """
        return ClonableList._repr_(self)

    def _repr_drawing(self):
        r"""
        Return a string representing a drawing of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 0]]
            ....: )
            sage: pp._repr_drawing()
            '[1 1 0]\n[1 1 0]\n[0 1 1]'
        """
        return str(matrix(self.get_array()))

    def get_tikz_options(self):
        return self.get_options()['tikz_options']

    def _to_tikz_diagram(self):
        tikz_options = self.get_tikz_options()
        grid_width = self.width() + 1
        grid_height = self.height() + 1
        drawing_tool = _drawing_tool(
            tikz_options,
            XY=lambda v: [v[0], grid_height-1-v[1]]
        )
        res = ""
        if self.size() == 0:
            res += drawing_tool.draw_line([0, 0], [1, 0])
            return res
        res += drawing_tool.draw_line([0, 0], [0, self.lower_heights()[0]])
        res += drawing_tool.draw_line(
            [grid_width-1, self.upper_heights()[grid_width-2]],
            [grid_width-1, self.lower_heights()[grid_width-2]]
        )
        res += drawing_tool.draw_line([0, 0], [self.upper_widths()[0], 0])
        res += drawing_tool.draw_line(
            [self.lower_widths()[grid_height-2], grid_height-1],
            [self.upper_widths()[grid_height-2], grid_height-1]
        )
        for w in xrange(1, grid_width-1):
            h1 = self.upper_heights()[w-1]
            h2 = self.lower_heights()[w]
            res += drawing_tool.draw_line([w, h1], [w, h2])
        for h in xrange(1, grid_height-1):
            w1 = self.lower_widths()[h-1]
            w2 = self.upper_widths()[h]
            res += drawing_tool.draw_line([w1, h], [w2, h])
        return res

    def _to_tikz_tree_with_bounce(self, directions=[0, 1]):
        res = ""
        tikz_options = self.get_tikz_options()
        if self.size() == 0:
            return res
        grid_width = self.width() + 1
        grid_height = self.height() + 1
        drawing_tool = _drawing_tool(
            tikz_options,
            XY=lambda v: [v[0] + .5, grid_height-1-v[1] - .5]
        )
        if 0 in directions:
            for node in self.get_right_nodes():
                res += drawing_tool.draw_point(
                    [node[1], node[0]], tikz_options['color_bounce_0']
                )
        if 1 in directions:
            for node in self.get_left_nodes():
                res += drawing_tool.draw_point(
                    [node[1], node[0]], tikz_options['color_bounce_1']
                )
        res += drawing_tool.draw_point([0, 0])
        return res

    def _to_tikz_bounce(self, directions=[0, 1]):
        res = ""
        tikz_options = self.get_tikz_options()
        grid_width = self.width() + 1
        grid_height = self.height() + 1
        drawing_tool = _drawing_tool(
            tikz_options,
            XY=lambda v: [v[0], grid_height-1-v[1]]
        )

        def draw_bounce(direction, color):
            if(
                len(self.bounce_path(direction))
                > len(self.bounce_path(1-direction))
            ):
                increase_size_line = 1
            else:
                increase_size_line = 0
            res = ""
            bp = self.bounce_path(direction)
            pos = [0, 0]
            pos[1-direction] += 1
            old = list(pos)
            for e in bp:
                pos[direction] += e
                res += drawing_tool.draw_line(
                    [old[1], old[0]], [pos[1], pos[0]],
                    color=color,
                    size=2*tikz_options['line_size'] + increase_size_line,
                )
                old[0], old[1] = pos
                direction = 1-direction
            return res
        if len(self.bounce_path(0)) > len(self.bounce_path(1)):
            if 0 in directions:
                res += draw_bounce(0, tikz_options['color_bounce_0'])
            if 1 in directions:
                res += draw_bounce(1, tikz_options['color_bounce_1'])
        else:
            if 1 in directions:
                res += draw_bounce(1, tikz_options['color_bounce_1'])
            if 0 in directions:
                res += draw_bounce(0, tikz_options['color_bounce_0'])
        return res

    def _to_tikz_tree(self):
        res = ""
        tikz_options = self.get_tikz_options()
        if self.size() == 0:
            return res
        grid_width = self.width() + 1
        grid_height = self.height() + 1
        drawing_tool = _drawing_tool(
            tikz_options,
            XY=lambda v: [v[0] + .5, grid_height-1-v[1] - .5]
        )
        for node in self.get_nodes():
            res += drawing_tool.draw_point([node[1], node[0]])
        res += drawing_tool.draw_point([0, 0])
        return res

    def get_node_position_at_row(self, row):
        h = row
        for w in range(self.width()):
            if self[h][w] == 1:
                return [h, w]
        return None

    def get_node_position_at_line(self, line):
        w = line
        for h in range(self.height()):
            if self[h][w] == 1:
                return [h, w]
        return None

    def get_node_position_from_box(
        self, box_position, direction, nb_crossed_nodes=[0]
    ):
        pos = list(box_position)
        if self[pos[0]][pos[1]] == 0:
            return None
        while self[pos[0]][pos[1]] != 0:
            pos[direction] -= 1
            if self.box_is_node(pos):
                nb_crossed_nodes[0] += 1
        pos[direction] += 1
        return pos

    def box_is_node(self, pos):
        if self[pos[0]][pos[1]] == 0:
            return False
        if self[pos[0]-1][pos[1]] == 0:
            return True
        if self[pos[0]][pos[1]-1] == 0:
            return True

    def box_is_root(self, box):
        return box[0] == 0 and box[1] == 0

    def get_path_in_pair_of_tree_from_box(self, box, direction):
        path = []
        while not self.box_is_root(box):
            nb_sons = [0]
            box = self.get_node_position_from_box(box, direction, nb_sons)
            direction = 1 - direction
            path.append(nb_sons[0]-1)
        path.reverse()
        return path

    def get_path_in_pair_of_tree_from_row(self, line):
        pos = self.get_node_position_at_row(line)
        return self.get_path_in_pair_of_tree_from_box(pos, 0)

    def get_path_in_pair_of_tree_from_line(self, line):
        pos = self.get_node_position_at_line(line)
        return self.get_path_in_pair_of_tree_from_box(pos, 1)

    def get_nodes(self):
        result = []
        for h in range(self.height()):
            result.append(self.get_node_position_at_row(h))
        for w in range(self.width()):
            result.append(self.get_node_position_at_line(w))
        return result

    def get_right_nodes(self):
        result = []
        for h in range(self.height()):
            path2 = self.get_path_in_pair_of_tree_from_row(h)
            if len(path2) % 2 == 1:
                result.append(self.get_node_position_at_row(h))
        for w in range(self.width()):
            path2 = self.get_path_in_pair_of_tree_from_line(w)
            if len(path2) % 2 == 0:
                result.append(self.get_node_position_at_line(w))
        return result

    def get_left_nodes(self):
        result = []
        for h in range(self.height()):
            path2 = self.get_path_in_pair_of_tree_from_row(h)
            if len(path2) % 2 == 0:
                result.append(self.get_node_position_at_row(h))
        for w in range(self.width()):
            path2 = self.get_path_in_pair_of_tree_from_line(w)
            if len(path2) % 2 == 1:
                result.append(self.get_node_position_at_line(w))
        return result

    def to_tikz(self):
        r"""
        Return the tikz code of the parallelogram polyomino.

        This code is the code present inside a tikz latex environment.
        """
        res = ""
        drawing_components = self.get_options()['drawing_components']
        if 'diagram' in drawing_components:
            res += self._to_tikz_diagram()
        directions = []
        if 'bounce_0' in drawing_components:
            directions.append(0)
        if 'bounce_1' in drawing_components:
            directions.append(1)
        if len(directions) != 0:
            res += self._to_tikz_bounce(directions)
        if 'tree' in drawing_components:
            res += self._to_tikz_tree()
            if len(directions) != 0:
                res += self._to_tikz_tree_with_bounce(directions)
        return res

    def geometry(self):
        r"""
        Return a pair [h, w] containing the height and the width of the
        parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1, 1, 1, 1], [1, 1, 1, 1, 0]]
            ....: )
            sage: pp.geometry()
            [1, 4]

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.geometry()
            [1, 1]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.geometry()
            [0, 1]
        """
        return [self.height(), self.width()]

    def size(self):
        r"""
        Return the size of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 0, 0, 0, 1, 0, 1, 1], [1, 0, 0, 0, 1, 1, 0, 0]]
            ....: )
            sage: pp.size()
            7

            sage: pp = ParallelogramPolyomino(
            ....:     [[0, 1], [1, 0]]
            ....: )
            sage: pp.size()
            1

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: pp.size()
            0
        """
        return len(self.upper_path())-1

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        For more on the latex options, see
        :meth:`ParallelogramPolyominoes.global_options`.
        """
        return self.get_options().dispatch(self, '_latex_', 'latex')

    def _latex_drawing(self):
        r"""
        Return a LaTeX version of ``self`` in a drawing style.
        """
        latex.add_package_to_preamble_if_available("tikz")
        tikz_options = self.get_tikz_options()
        res = "\n\\begin{tikzpicture}[scale=%s]" % (tikz_options['scale'])
        res += self.to_tikz()
        res += "\n\\end{tikzpicture}"
        return res

    def _latex_list(self):
        r"""
        Return a LaTeX version of ``self`` in a list style.
        """
        return "\\[%s\\]" % self._repr_list()
        NotImplemented


class ParallelogramPolyominoesFactory(SetFactory):
    r"""
    The parallelogram polyominoes factory.
    """
    def __call__(self, size=None, policy=None):
        r"""
        Return a family of parallelogram polyominoes enumerated with the 
        parameter constraints.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes(size=3)
            sage: PPS
            parallelogram polyominoes of size 3
            sage: sorted( list(PPS) )
            [[[0, 0, 0, 1], [1, 0, 0, 0]],
             [[0, 0, 1, 1], [1, 0, 1, 0]],
             [[0, 0, 1, 1], [1, 1, 0, 0]],
             [[0, 1, 0, 1], [1, 1, 0, 0]],
             [[0, 1, 1, 1], [1, 1, 1, 0]]]
        """
        if policy is None:
            policy = self._default_policy

        if isinstance(size, (Integer, int)):
            return ParallelogramPolyominoes_size(size, policy)
        if size is None:
            return ParallelogramPolyominoes_all(policy)
        raise ValueError("Invalid argument for Parallelogram Polyominoes "
                         "Factory.")

    def add_constraints(self, cons, (args, opts)):
        r"""
        This function permit to add some enumeration constraint to the 
        factory. The factory make a family using the given constraints.
        """
        return cons + args

    @lazy_attribute
    def _default_policy(self):
        return TopMostParentPolicy(self, (), ParallelogramPolyomino)

    def _repr_(self):
        r"""
        Return the string representation of the parallelogram polyominoes 
        factory.

        EXAMPLES::

            sage: ParallelogramPolyominoes
            Factory for parallelogram polyominoes
        """
        return "Factory for parallelogram polyominoes"

ParallelogramPolyominoes = ParallelogramPolyominoesFactory()
ParallelogramPolyominoes.__doc__ = \
    ParallelogramPolyominoesFactory.__call__.__doc__


class ParallelogramPolyominoes_size(
    ParentWithSetFactory, UniqueRepresentation
):
    r"""
    The parallelogram polyominoes of size `n`.

    EXAMPLES::

        sage: PPS = ParallelogramPolyominoes(3)
        sage: PPS
        parallelogram polyominoes of size 3
        sage: sorted( list(PPS) )
        [[[0, 0, 0, 1], [1, 0, 0, 0]],
         [[0, 0, 1, 1], [1, 0, 1, 0]],
         [[0, 0, 1, 1], [1, 1, 0, 0]],
         [[0, 1, 0, 1], [1, 1, 0, 0]],
         [[0, 1, 1, 1], [1, 1, 1, 0]]]
    """
    def __init__(self, size, policy):
        r"""
        Construct a set of Parallelogram Polyominoes of a given size.

        EXAMPLES::

            sage: ParallelogramPolyominoes(3)
            parallelogram polyominoes of size 3
        """
        self._size = size
        ParentWithSetFactory.__init__(
            self, (size, ), policy, category=FiniteEnumeratedSets()
        )

    def _repr_(self):
        r"""
        Return the string representation of the set of
        parallelogram polyominoes

        EXAMPLES::

            sage: ParallelogramPolyominoes(3)
            parallelogram polyominoes of size 3
        """
        return "parallelogram polyominoes of size %s" % (self._size)

    def an_element(self):
        r"""
        Return an element of a parallelogram polyomino of a given size.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes(3)
            sage: PPS.an_element() in PPS
            True
        """
        return next(self.__iter__())

    def check_element(self, el, check):
        r"""
        Check is a given element `el` is in the set of parallelogram
        polyominoes of a fixed size.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes(2)
            sage: ParallelogramPolyomino([[0, 1, 1], [1, 1, 0]])  in PPS
            True
        """
        if el.size() != self.size():
            raise ValueError(
                "The parallelogram polyomino have a Wrong size: %s" % (
                    el.size()
                )
            )

    def cardinality(self):
        r"""
        Return the number of parallelogram polyominoes.

        EXAMPLES::

            sage: ParallelogramPolyominoes(1).cardinality()
            1
            sage: ParallelogramPolyominoes(2).cardinality()
            2
            sage: ParallelogramPolyominoes(3).cardinality()
            5

            sage: all([
            ....:     ParallelogramPolyominoes(i).cardinality()
            ....:     == catalan_number(i)
            ....:     for i in range(6)
            ....: ])
            True

            sage: all([
            ....:     ParallelogramPolyominoes(i).cardinality()
            ....:     == len(list(ParallelogramPolyominoes(i)))
            ....:     for i in range(6)
            ....: ])
            True
        """
        return catalan_number(self.size())

    def __iter__(self):
        r"""
        Return a parallelogram polyomino generator.

        EXAMPLES::

            sage: len(list(ParallelogramPolyominoes(3))) == 5
            True
            sage: all([
            ....:     pp in ParallelogramPolyominoes()
            ....:     for pp in ParallelogramPolyominoes(3)
            ....: ])
            True
        """
        from sage.combinat.dyck_word import DyckWords
        for dyck in DyckWords(self.size()):
            yield ParallelogramPolyomino.from_dyck_word(dyck)

    def get_options(self):
        return self.global_options

    def size(self):
        r"""
        Return the size of the parallelogram polyominoes generated by this
        parent.

        EXAMPLES::

            sage: ParallelogramPolyominoes(0).size()
            0
            sage: ParallelogramPolyominoes(1).size()
            1
            sage: ParallelogramPolyominoes(5).size()
            5
        """
        return self._size

    def set_options(self, *get_value, **set_value):
        r"""
        Set new options to the object.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes(3)
            sage: PPS.set_options(
            ....:     drawing_components=dict(
            ....:         diagram = True,
            ....:         bounce_0 = True,
            ....:         bounce_1 = True,
            ....:     )
            ....: )
            sage: pp = PPS[0]
            sage: view(pp) # not tested
        """
        self.global_options(*get_value, **set_value)

    global_options = ParallelogramPolyominoesOptions


class ParallelogramPolyominoes_all(
    ParentWithSetFactory, DisjointUnionEnumeratedSets
):
    r"""
    This class enumerates all the parallelogram polyominoes.

    EXAMPLES::

        sage: PPS = ParallelogramPolyominoes()
        sage: PPS
        parallelogram polyominoes
    """
    def __init__(self, policy):
        r"""
        Construct the set of all parallelogram polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: PPS
            parallelogram polyominoes

            sage: ParallelogramPolyomino([[0, 1, 1], [1, 1, 0]])  in PPS
            True

            sage: PPS = ParallelogramPolyominoes()
            sage: next(PPS.__iter__()) in PPS
            True
        """
        ParentWithSetFactory.__init__(
            self, (), policy, category=FiniteEnumeratedSets()
        )
        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                NonNegativeIntegers(),
                lambda  n: ParallelogramPolyominoes_size(
                    n, policy=self.facade_policy()
                )
            ),
            facade=True, keepkey=False, category=self.category()
        )

    def _repr_(self):
        r"""
        Return a string representation of the set of parallelogram polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: PPS
            parallelogram polyominoes
        """
        return "parallelogram polyominoes"

    def check_element(self, el, check):
        r"""
        Check is a given element `el` is in the set of parallelogram
        polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: ParallelogramPolyomino([[0, 1, 1], [1, 1, 0]])  in PPS
            True
        """
        pass

    def get_options(self):
        r"""
        Return all the options associated with the set of
        parallelogram polyominoes.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: options = PPS.get_options()
            sage: options
            options for Parallelogram Polyominoes
            sage: options()
            Current options for Parallelogram Polyominoes
              - display:            list
            ...
        """
        return self.global_options

    def set_options(self, *get_value, **set_value):
        r"""
        Set new options to the object.

        EXAMPLES::

            sage: PPS = ParallelogramPolyominoes()
            sage: PPS.set_options(
            ....:     drawing_components=dict(
            ....:         diagram = True,
            ....:         bounce_0 = True,
            ....:         bounce_1 = True,
            ....:     )
            ....: )
            sage: pp = next(iter(PPS))
            sage: view(pp) # not tested
        """
        self.global_options(*get_value, **set_value)

    global_options = ParallelogramPolyominoesOptions
