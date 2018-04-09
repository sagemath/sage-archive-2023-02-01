r"""
Super Partitions

AUTHORS:

- Mike Zabrocki

A super partition of size ``n`` and fermionic sector ``m`` is a
pair consisting of a strict partition of some integer `r` of 
length ``m``(that may end in a `0`) and an integer partition of
``n-r``.

This module provides tools for manipulating super partitions.

Super partitions are the indexing set for symmetric functions in
super space.

Super partitions may be input in two different formats: one as a pair
consisiting of fermionic (strict partition) and a bosonic (partition) part
and as a list of integer values where the negative entries come first
and are listed in strict order followed by the positive values in weak
order.

A super partition is displayed as two partitions separated by a semicolon
as a default.  Super partitions may also be displayed as a weakly increasing
sequence of integers that are strict if the numbers are not positive.

These combinatorial objects index the space of symmetric polynomials in
two sets of variables, one commuting and one anti-commuting, and they
are known as symmetric functions in super space (hence the origin of the
name super partitions).

EXAMPLES::

    sage: SuperPartitions()
    SuperPartitions
    sage: SuperPartitions(2)
    SuperPartitions of 2
    sage: SuperPartitions(2).cardinality()
    8
    sage: SuperPartitions(4,2)
    SuperPartitions of 4 and of fermionic sector 2
    sage: [[2,0],[1,1]] in SuperPartitions(4,2)
    True
    sage: [[1,0],[1,1]] in SuperPartitions(4,2)
    False
    sage: [[1,0],[2,1]] in SuperPartitions(4)
    True
    sage: [[1,0],[2,2,1]] in SuperPartitions(4)
    False
    sage: [[1,0],[2,1]] in SuperPartitions()
    True
    sage: [[1,1],[2,1]] in SuperPartitions()
    False
    sage: [-2, 0, 1, 1] in SuperPartitions(4,2)
    True
    sage: [-1, 0, 1, 1] in SuperPartitions(4,2)
    False
    sage: [-2, -2, 2, 1] in SuperPartitions(7,2)
    False
"""
from sage.combinat.combinat import CombinatorialElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from six.moves import builtins
from sage.combinat.partition import Partition, Partitions
from sage.combinat.permutation import Permutations
from sage.combinat.composition import Composition
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.rings.integer import Integer
from sage.structure.global_options import GlobalOptions
from sage.rings.all import ZZ
from sage.misc.all import uniq

class SuperPartition(CombinatorialElement):
    r"""
    The element class of super partitions

    A super partition of size ``n`` and fermionic sector ``m`` is a
    pair consisting of a strict partition of some integer `r` of 
    length ``m`` (that may end in a `0`) and an integer partition of
    ``n-r``.

    Super Partitions are a list of two lists.  The degree (or bosonic
    degree) is equal to the sum of the entries in the lists.  The
    fermionic degree (or fermionic sector) is equal to the length of
    the first list.
    
    EXAMPLES::

        sage: sp = SuperPartition([[1,0],[2,2,1]]); sp
        [1, 0; 2, 2, 1]
        sage: sp[0]
        [1, 0]
        sage: sp[1]
        [2, 2, 1]
        sage: sp.fermionic_degree()
        2
        sage: sp.bosonic_degree()
        6
        sage: sp.length()
        5
        sage: sp.conjugate()
        [4, 2; ]
    """
    @staticmethod
    def __classcall_private__(cls, lst):
        r"""

        EXAMPLES::

            sage: SuperPartition([[1],[1]]).parent()
            SuperPartitions
            sage: SuperPartition([[1],[1]])
            [1; 1]
            sage: SuperPartition([-1, 1])
            [1; 1]
            sage: SuperPartition([[1,1],[1]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1], [1]] must be an element of SuperPartitions
            sage: SuperPartition([-1,1])
            [1; 1]
            sage: SuperPartition([])
            [; ]
        """
        if lst in SuperPartitions():
            if len(lst)==0:
                return SuperPartitions()([[],[]])
            elif isinstance(lst[0], builtins.list):
                return SuperPartitions()([[Integer(a) for a in lst[0]],\
                    [Integer(a) for a in lst[1]]])
            else:
                return SuperPartitions()([[-a for a in lst if a<=0], \
                    [a for a in lst if a>0]])
        else:
            raise ValueError("%s must be an element of SuperPartitions"%lst)

    def _repr_(self):
        r"""
        Represention of a super partition

        A superpartition is represented by the antisymmetric and symmetric
        parts separated by a semicolon.

        EXAMPLES::

            sage: SuperPartition([[1],[1]])
            [1; 1]
            sage: SuperPartition([[],[1]])
            [; 1]
            sage: SuperPartition([])
            [; ]
            sage: SuperPartitions.options.display = "list"
            sage: SuperPartition([[1],[1]])
            [-1, 1]
            sage: SuperPartition([[],[1]])
            [1]
            sage: SuperPartition([-2,-1,0,2,1])
            [-2, -1, 0, 2, 1]
            sage: SuperPartitions.options.display = "pair"
            sage: SuperPartition([[1],[1]])
            [[1], [1]]
            sage: SuperPartition([[],[1]])
            [[], [1]]
            sage: SuperPartition([-2,-1,0,2,1])
            [[2, 1, 0], [2, 1]]
            sage: SuperPartitions.options.display = "default"
        """
        display = self.parent().options.display
        if display=="default":
            return '['+', '.join(str(a) for a in self.antisymmetric_part())+\
                '; '+', '.join(str(a) for a in self.symmetric_part())+']'
        elif display=="pair":
            return self._repr_pair_()
        elif display=="list":
            return self._repr_list_()

    def _repr_pair_(self):
        r"""
        Represention of a super partition as a pair

        A superpartition is represented by a list consisting of the
        antisymmetric and symmetric parts.

        EXAMPLES::

            sage: SuperPartition([[1],[1]])._repr_pair_()
            '[[1], [1]]'
            sage: SuperPartition([[],[1]])._repr_pair_()
            '[[], [1]]'
            sage: SuperPartition([[],[]])._repr_pair_()
            '[[], []]'
        """
        return str(self._list)

    def _repr_list_(self):
        r"""
        Represention of a super partition as a list

        A superpartition is represented by a list consisting of the
        negative values for the antisymmetric part listed first followed by
        positive values for the symmetric part

        EXAMPLES::

            sage: SuperPartition([[1],[1]])._repr_list_()
            '[-1, 1]'
            sage: SuperPartition([[],[1]])._repr_list_()
            '[1]'
            sage: SuperPartition([[],[]])._repr_list_()
            '[]'
        """
        return str([-a for a in self._list[0]]+self._list[1])

    def _latex_(self):
        r"""
        Latex a super partition

        A superpartition is represented by the antisymmetric and symmetric
        parts separated by a semicolon.

        EXAMPLES::

            sage: latex(SuperPartition([[1],[1]]))
            (1; 1)
            sage: latex(SuperPartition([[],[1]]))
            (; 1)
        """
        return '('+','.join(str(a) for a in self.antisymmetric_part())+\
            '; '+', '.join(str(a) for a in self.symmetric_part())+')'

    def to_list(self):
        r"""
        The list of two lists with the antisymmetric and symmetric parts

        EXAMPLES::

            sage: SuperPartition([[1],[1]]).to_list()
            [[1], [1]]
            sage: SuperPartition([[],[1]]).to_list()
            [[], [1]]
        """
        return self._list

    def to_composition(self):
        r"""
        Concatenate the antisymmetric and symmetric parts to a composition

        OUTPUT:

        - a (possibly weak) composition

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).to_composition()
            [3, 1, 2, 2, 1]
            sage: SuperPartition([[2,1,0],[3,3]]).to_composition()
            [2, 1, 0, 3, 3]
            sage: SuperPartition([[2,1,0],[3,3]]).to_composition().parent()
            Compositions of non-negative integers
        """
        return Composition(self._list[0]+self._list[1])

    def to_partition(self):
        r"""
        Concatenate and sort the antisymmetric and symmetric parts to a partition

        OUTPUT:

        - a partition

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).to_partition()
            [3, 2, 2, 1, 1]
            sage: SuperPartition([[2,1,0],[3,3]]).to_partition()
            [3, 3, 2, 1]
            sage: SuperPartition([[2,1,0],[3,3]]).to_partition().parent()
            Partitions
        """
        return Partition(sorted(self._list[0]+self._list[1],reverse=True))

    def antisymmetric_part(self):
        r"""
        The antisymmetric part as a list of strictly decreasing integers

        OUTPUT:

        - a list

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).antisymmetric_part()
            [3, 1]
            sage: SuperPartition([[2,1,0],[3,3]]).antisymmetric_part()
            [2, 1, 0]
        """
        return self._list[0]
    a_part = antisymmetric_part

    def symmetric_part(self):
        r"""
        The symmetric part as a list of weakly decreasing integers

        OUTPUT:

        - a list

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).symmetric_part()
            [2, 2, 1]
            sage: SuperPartition([[2,1,0],[3,3]]).symmetric_part()
            [3, 3]
        """
        return self._list[1]
    s_part = symmetric_part

    def bosonic_degree(self):
        r"""
        The sum of the sizes of the antisymmetric and symmetric parts

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).bosonic_degree()
            9
            sage: SuperPartition([[2,1,0],[3,3]]).bosonic_degree()
            9
        """
        return sum(self.antisymmetric_part()+self.symmetric_part())
    degree = bosonic_degree
    
    def fermionic_degree(self):
        r"""
        The length of the antisymmetric part

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).fermionic_degree()
            2
            sage: SuperPartition([[2,1,0],[3,3]]).fermionic_degree()
            3
        """
        return len(self.antisymmetric_part())

    fermionic_sector = fermionic_degree

    def bi_degree(self):
        r"""
        A pair consisting of the bosonic and fermionic degree

        OUTPUT:

        - a tuple of two integers

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).bi_degree()
            (9, 2)
            sage: SuperPartition([[2,1,0],[3,3]]).bi_degree()
            (9, 3)
        """
        return (self.bosonic_degree(), self.fermionic_degree())

    def length(self):
        r"""
        The sum of the lengths of the antisymmetric and symmetric part

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).length()
            5
            sage: SuperPartition([[2,1,0],[3,3]]).length()
            5
        """
        return self.fermionic_degree()+len(self.symmetric_part())

    def bosonic_length(self):
        r"""
        The length of the partition of the symmetric part

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).bosonic_length()
            3
            sage: SuperPartition([[2,1,0],[3,3]]).bosonic_length()
            2
        """
        return len(self.symmetric_part())

    def shape_circled_diagram(self):
        r"""
        A concatenated partition with an extra cell for each antisymmetric part

        OUTPUT:

        - a partition

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).shape_circled_diagram()
            [4, 2, 2, 2, 1]
            sage: SuperPartition([[2,1,0],[3,3]]).shape_circled_diagram()
            [3, 3, 3, 2, 1]
        """
        return Partition(sorted([a+1 for a in self.antisymmetric_part()] \
                                +self.symmetric_part(),reverse=True))

    @staticmethod
    def from_circled_diagram(shape, corners):
        r"""
        Construct a SuperPartition from a circled diagram

        A circled diagram consists of a partition of the concatenation of
        the antisymmetric and symmetric parts and a list of addable cells
        of the partition which indicate the location of the circled cells

        INPUT:
        
        - ``shape`` -- a partition or list of integers
        - ``corners`` -- a list of removable cells of ``shape``

        OUTPUT:

        - a SuperPartition

        EXAMPLES::

            sage: SuperPartition.from_circled_diagram([3, 2, 2, 1, 1], [(0, 3), (3, 1)])
            [3, 1; 2, 2, 1]
            sage: SuperPartition.from_circled_diagram([3, 3, 2, 1], [(2, 2), (3, 1), (4, 0)])
            [2, 1, 0; 3, 3]
            sage: from_cd = SuperPartition.from_circled_diagram
            sage: all(sp==from_cd(*sp.to_circled_diagram()) for sp in SuperPartitions(4))
            True
        """
        return SuperPartition([sorted([c[1] for c in corners], reverse=True),
                 [shape[i] for i in range(len(shape)) if
                    i not in [c[0] for c in corners]]])

    def to_circled_diagram(self):
        r"""
        The shape of the circled diagram and a list of addable cells

        A circled diagram consists of a partition for the outer shape
        and a list of removable cells of the partition indicating the
        location of the circled cells

        OUTPUT:

        - a list consisting of a partition and a list of pairs of integers

        EXAMPLES::

            sage: SuperPartition([[3,1],[2,2,1]]).to_circled_diagram()
            [[3, 2, 2, 1, 1], [(0, 3), (3, 1)]]
            sage: SuperPartition([[2,1,0],[3,3]]).to_circled_diagram()
            [[3, 3, 2, 1], [(2, 2), (3, 1), (4, 0)]]
            sage: from_cd = SuperPartition.from_circled_diagram
            sage: all(sp==from_cd(*sp.to_circled_diagram()) for sp in SuperPartitions(4))
            True
        """
        shape = self.to_partition()
        corners = [c for c in shape.addable_cells() if c[1] in self.antisymmetric_part()]
        return [shape, corners]
    
    def conjugate(self):
        r"""
        Conjugate of a super-partition

        The conjugate of a super-partition is defined by conjugating the circled
        diagram.

        OUPUT:

        - a SuperPartition        

        EXAMPLES::

            sage: SuperPartition([[3, 1, 0], [4, 3, 2, 1]]).conjugate()
            [6, 4, 1; 3]
            sage: all(sp==sp.conjugate().conjugate() for sp in SuperPartitions(4))
            True
            sage: all(sp.conjugate() in SuperPartitions(3,2) for sp in SuperPartitions(3,2))
            True
        """
        sd = self.to_circled_diagram()
        return SuperPartition.from_circled_diagram(sd[0].conjugate(), [(j,i) for (i,j) in sd[1]])

    def zee(self):
        r"""
        The centralizer size of a permutation of cycle type symmetric part of ``self``

        OUTPUT:

        - a positive integer

        EXAMPLES::

            sage: SuperPartition([[1,0],[3,1,1]]).zee()
            6
            sage: SuperPartition([[1],[2,2,1]]).zee()
            8
            sage: sum(1/sp.zee() for sp in SuperPartitions(6,0))
            1
        """
        return Partition(self.symmetric_part()).centralizer_size()

    def sign(self):
        r"""
        The sign of a permutation of cycle type the symmetric part of ``self``

        OUTPUT:

        - either 1 or -1

        EXAMPLES::

            sage: SuperPartition([[1,0],[3,1,1]]).sign()
            -1
            sage: SuperPartition([[1,0],[3,2,1]]).sign()
            1
            sage: sum(sp.sign()/sp.zee() for sp in SuperPartitions(6,0))
            0
        """
        return (-1)**(self.degree()-len(self.symmetric_part()))

    def Bruhat_less(self):
        return NotImplemented
    def S_order_less(self):        
        return NotImplemented
    def T_order_less(self):        
        return NotImplemented
    def dominates(self, other):
        r"""
        Return ``True`` if and only if ``self`` dominates ``other``.

        If the symmetric and anti-symmetric parts of ``self`` and ``other``
        are not the same size then the result is ``False``.

        EXAMPLES::

            sage: LA = SuperPartition([[2,1],[2,1,1]])
            sage: LA.dominates([[2,1],[3,1]])
            False
            sage: LA.dominates([[2,1],[1,1,1,1]])
            True
            sage: LA.dominates([[3],[2,1,1]])
            False
            sage: LA.dominates([[1],[1]*6])
            False
        """
        return self.degree()==sum(other[0])+sum(other[1]) and \
            Partition(self.antisymmetric_part()).dominates(other[0]) and \
            Partition(self.symmetric_part()).dominates(other[1])

    def add_horizontal_border_strip_star(self, h):
        r"""
        a list of super partitions that differ from ``self`` by a horizontal strip

        INPUT:

        - ``h`` -- number of cells in the horizontal strip

        OUPUT:

        - a list of super partitions

        EXAMPLES::

            sage: SuperPartition([[4,1],[3]]).add_horizontal_border_strip_star(3)
            [[4, 1; 3, 3],
             [4, 1; 4, 2],
             [3, 1; 5, 2],
             [4, 1; 5, 1],
             [3, 1; 6, 1],
             [4, 0; 4, 3],
             [3, 0; 5, 3],
             [4, 0; 5, 2],
             [3, 0; 6, 2],
             [4, 1; 6],
             [3, 1; 7]]
            sage: SuperPartition([[2,1],[3]]).add_horizontal_border_strip_star(2)
            [[2, 1; 3, 2], [2, 1; 4, 1], [2, 0; 3, 3], [2, 0; 4, 2], [2, 1; 5]]
        """
        (sp1, circ_list) = SuperPartition(self).to_circled_diagram()
        nsp=[list(la)+[0] for la in sp1.add_horizontal_border_strip(h)]
        sp1= sp1+[0]
        out=[]
        rows_with_circle=[a[0] for a in circ_list]
        for elt in nsp:
            row_changed = [row1-row2 for row1,row2 in zip(elt,sp1)]
            new_sp = [elt,[(i[0]+1,elt[i[0]+1]) for i in circ_list \
                if row_changed[i[0]]!=0]+[(i) for i in circ_list \
                if row_changed[i[0]]==0]]
            if len(uniq([k for (j,k) in new_sp[1]]))==len(new_sp[1]): 
                out+=[SuperPartition.from_circled_diagram(*new_sp)]
        return out

    def add_horizontal_border_strip_bar(self, h):
        r"""
        a list of super partitions that differ from ``self`` by a horizontal strip

        INPUT:

        - ``h`` -- number of cells in the horizontal strip

        OUPUT:

        - a list of super partitions

        EXAMPLES::

            sage: SuperPartition([[4,1],[5,4]]).add_horizontal_border_strip_bar(3)
            [[4, 3; 5, 4, 1],
             [4, 1; 5, 4, 3],
             [4, 2; 5, 5, 1],
             [4, 1; 5, 5, 2],
             [4, 2; 6, 4, 1],
             [4, 1; 6, 4, 2],
             [4, 1; 6, 5, 1],
             [4, 1; 7, 4, 1],
             [4, 3; 5, 5],
             [4, 3; 6, 4],
             [4, 2; 6, 5],
             [4, 2; 7, 4],
             [4, 1; 7, 5],
             [4, 1; 8, 4]]
            sage: SuperPartition([[3,1],[5]]).add_horizontal_border_strip_bar(2)
            [[3, 2; 5, 1],
             [3, 1; 5, 2],
             [4, 1; 5, 1],
             [3, 1; 6, 1],
             [4, 2; 5],
             [3, 2; 6],
             [4, 1; 6],
             [3, 1; 7]]
        """
        (sp1, circ_list) = SuperPartition(self).to_circled_diagram()
        nsp=[list(la)+[0] for la in sp1.add_horizontal_border_strip(h)]
        sp1= sp1+[0]
        out=[]
        for asp in nsp:
            asp = asp+[0]
            change_in_rows = [asp[i]-sp1[i] for i in range(len(sp1))]
            #map(operator.sub,asp[:len(sp1)],sp1)
            moved_circ_list= [[]for i in range(len(circ_list))]
            for i in range(len(circ_list)):
                if change_in_rows[circ_list[i][0]]==0:
                    moved_circ_list[i].append(circ_list[i])
                else:
                    if circ_list[i][0]==0:
                        moved_circ_list[i].append((0,circ_list[i][1]+change_in_rows[0]))
                        if circ_list[i][1]==asp[1]:
                            moved_circ_list[i].append((1,asp[1]))
                    else:
                        if circ_list[i][1]+change_in_rows[circ_list[i][0]]<sp1[circ_list[i][0]-1]:
                            moved_circ_list[i].append((circ_list[i][0],circ_list[i][1]+change_in_rows[circ_list[i][0]]))
                        if asp[circ_list[i][0]+1]==sp1[circ_list[i][0]]:
                            moved_circ_list[i].append((circ_list[i][0]+1,circ_list[i][1]))
            out+=[[moved_circ_list,asp]]
        result = []
        for i in out:
            if not i[0]:
                result += [[i[1],i[0]]]
            else:
                x = reduce(lambda a,b: [item_a+item_b for item_a in a for item_b in b], i[0])
                for j in x:
                    result += [[i[1],zip(j,j[1:])[::2]]]
        return [SuperPartition.from_circled_diagram(*i) for i in result if len(i[1])==len(self[0])]

class SuperPartitions(UniqueRepresentation, Parent):
    r"""
    The class of super partitions
    
    A super partition of size ``n`` and fermionic sector ``m`` is a
    pair consisting of a strict partition of some integer `r` of 
    length ``m``(that may end in a `0`) and an integer partition of
    ``n-r``.
    
    INPUT:
    
    - ``n`` -- an integer (optional: default ``None``)

    - ``m`` -- if ``n`` is specified, an integer (optional: default ``None``)

    Super partitions are the indexing set for symmetric functions in
    super space.

    EXAMPLES::

        sage: SuperPartitions()
        SuperPartitions
        sage: SuperPartitions(2)
        SuperPartitions of 2
        sage: SuperPartitions(2).cardinality()
        8
        sage: SuperPartitions(4,2)
        SuperPartitions of 4 and of fermionic sector 2
        sage: [[2,0],[1,1]] in SuperPartitions(4,2)
        True
        sage: [[1,0],[1,1]] in SuperPartitions(4,2)
        False
        sage: [[1,0],[2,1]] in SuperPartitions(4)
        True
        sage: [[1,0],[2,2,1]] in SuperPartitions(4)
        False
        sage: [[1,0],[2,1]] in SuperPartitions()
        True
        sage: [[1,1],[2,1]] in SuperPartitions()
        False
    """
    @staticmethod
    def __classcall_private__(self, n=None, m=None, **kwargs):
        r"""
        The classes of super partitions of all, ``n`` and ``n|m``
        """
        if n is None:
            return SuperPartitions_all()
        elif isinstance(n, (int,Integer)):
            if m is None:
                return SuperPartitions_n(n)
            elif isinstance(m, (int,Integer)):
                return SuperPartitions_n_m(n,m)
            raise ValueError("m must be an integer")
        raise ValueError("n must be an integer")

    def __init__(self, is_infinite=False):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: SP = SuperPartitions()
            sage: TestSuite(SP).run()
        """
        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())

    Element = SuperPartition

    options=GlobalOptions('SuperPartitions',
        #module='sage.combinat.superpartition',
        doc=r"""
        Set the global options for elements of the SuperPartition class. The
        defaults are for SuperPartitions to be displayed in a list notation
        with the fermionic part and the bosonic part separated by a semicolon.
        There is a slight disadvantage to this notation because a list
        containing a semicolon can not be used as input for a super partition.

        """,
        end_doc="""
        EXAMPLES::

            sage: sp = SuperPartition([[1,0],[2,2,1]])
            sage: SuperPartitions.options.display
            'default'
            sage: sp
            [1, 0; 2, 2, 1]
            sage: Permutations.options.display='list'
            sage: sp
            [[1, 0], [2, 2, 1]]
            sage: Permutations.options.display='default'
        """,
        display=dict(default="default",
                     description="Specifies how the super partitions should "
                     "be printed",
                     values=dict(list="the super partitions are displayed in "
                                  "a list of two lists",
                                  pair="the super partition is displayed as a "
                                  "list of integers",
                                 default="the super partition is displayed in "
                                      "a form [fermionic part; bosonic part]"),
                     case_sensitive=False),
    )

    def _element_constructor_(self, lst):
        """
        Construct an element with ``self`` as parent.

        EXAMPLES::

            sage: SP = SuperPartitions()
            sage: SP([[],[3,3,1]])
            [; 3, 3, 1]
            sage: SP([[],[3,3,1]]) in SP
            True
            sage: SP([[],[3,3,1]]).parent()
            SuperPartitions
            sage: SuperPartitions(7)([[],[3,3,1]])
            [; 3, 3, 1]
            sage: SuperPartitions(7,0)([[],[3,3,1]])
            [; 3, 3, 1]
            sage: SuperPartitions(7,1)([[],[3,3,1]])
            Traceback (most recent call last):
            ...
            ValueError: [[], [3, 3, 1]] not in SuperPartitions of 7 and of fermionic sector 1
        """
        if isinstance(lst, SuperPartition):
            lst = lst._list
        if lst not in self:
            raise ValueError("%s not in %s"%(lst, self))
        if len(lst)==0:
            return self.element_class(self, [[], []])
        elif isinstance(lst[0], builtins.list):
            return self.element_class(self, [lst[0],[a for a in lst[1] if a>0]])
        else:
            return self.element_class(self, [[-a for a in lst if a<=0], [a for a in lst if a>0]])

    def __contains__(self, x):
        """
        TESTS::

            sage: [[1],[2,1]] in SuperPartitions()
            True
            sage: [[],[]] in SuperPartitions()
            True
            sage: [[0],[]] in SuperPartitions()
            True
            sage: [[],[0]] in SuperPartitions()
            True
            sage: [-1, 2, 1] in SuperPartitions()
            True
            sage: [2, -1, 1, 0] in SuperPartitions()
            True
            sage: [2, 0, 1, -1] in SuperPartitions()
            False
            sage: [] in SuperPartitions()
            True
            sage: [0] in SuperPartitions()
            True
        """
        if isinstance(x, SuperPartition):
            return True
        elif isinstance(x, builtins.list) and all(isinstance(i, (int, Integer))
            or i in ZZ for i in x):
            sp = [a for a in x if a<=0]
            return all(sp[i]>sp[i-1] for i in range(1,len(sp))) \
                and [a for a in x if a>0] in Partitions()
        elif isinstance(x, builtins.list) and len(x)==2 and \
          isinstance(x[0], builtins.list) and isinstance(x[1], builtins.list):
            for i in x[0]+x[1]:
                if (not isinstance(i, (int, Integer))) and i not in ZZ:
                    return False
                if i < 0:
                    return False
            return all(x[0][i]>x[0][i+1] for i in range(len(x[0])-1)) \
                and all(x[1][i]>=x[1][i+1] for i in range(len(x[1])-1)) \
                and (len(x[0])==0 or x[0][-1]>=0) and (len(x[1])==0 or x[1][-1]>=0)
        else:
            return False

class SuperPartitions_n_m(SuperPartitions):
    @staticmethod
    def __classcall_private__(cls, n, m):
        """
        Standardize input to ensure a unique representation.

        EXAMPLES::

            sage: SP = SuperPartitions(5,2)
            sage: SP2 = SuperPartitions(int(5),int(2))
            sage: SP3 = SuperPartitions(ZZ(5),int(2))
            sage: SP is SP2
            True
            sage: SP is SP3
            True
        """
        return super(SuperPartitions_n_m, cls).__classcall__(cls, \
                     Integer(n), Integer(m))

    def __init__(self, n, m):
        """
        Initialize ``self``.

        TESTS::

            sage: SP = SuperPartitions(3,2)
            sage: TestSuite(SP).run()
        """
        self.n = n
        self.m = m
        SuperPartitions.__init__(self, False)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(SuperPartitions(3,2))
            'SuperPartitions of 3 and of fermionic sector 2'
        """
        return "SuperPartitions of %s and of fermionic sector %s"%(self.n, self.m)

    def __contains__(self, x):
        """
        TESTS::

            sage: [[3,2,1,0],[2]] in SuperPartitions(8,4)
            True
            sage: [[3,2,1,0],[]] in SuperPartitions(6,3)
            False
            sage: [[],[]] in SuperPartitions(0,0)
            True
            sage: [[0],[]] in SuperPartitions(0,1)
            True
            sage: [[],[]] in SuperPartitions(0,1)
            False
            sage: [-3,-2,-1,0,2] in SuperPartitions(8,4)
            True
            sage: [0] in SuperPartitions(0,0)
            False
            sage: [] in SuperPartitions(0,0)
            True
            sage: [0] in SuperPartitions(0,1)
            True
        """
        if x in SuperPartitions():
            if len(x)==0:
                return self.n==0 and self.m==0
            if isinstance(x[0], builtins.list):
                n=sum(x[0]+x[1])
                m=len(x[0])
            else:
                n=sum(abs(a) for a in x)
                m=len([a for a in x if a<=0])
            return n==self.n and m==self.m
        else:
            return False
    
    def __iter__(self):
        r"""
        An iterator for SuperPartitions of degree ``n`` and sector ``m``

        EXAMPLES::

            sage: SuperPartitions(6,2).cardinality()
            28
            sage: SuperPartitions(6,4).first()
            [3, 2, 1, 0; ]
        """
        for r in range(self.n+1):
            for p1 in Partitions(r):
                for p0 in Partitions(self.n-r, max_slope=-1, length=self.m):
                    yield self.element_class(self, [list(p0), list(p1)])
                for p0 in Partitions(self.n-r, max_slope=-1, length=self.m-1):
                    yield self.element_class(self, [list(p0)+[0], list(p1)])

class SuperPartitions_n(SuperPartitions):
    @staticmethod
    def __classcall_private__(cls, n):
        """
        Standardize input to ensure a unique representation.

        EXAMPLES::

            sage: SP = SuperPartitions(5)
            sage: SP2 = SuperPartitions(int(5))
            sage: SP3 = SuperPartitions(ZZ(5))
            sage: SP is SP2
            True
            sage: SP is SP3
            True
        """
        return super(SuperPartitions_n, cls).__classcall__(cls, Integer(n))

    def __init__(self, n):
        """
        Initialize ``self``.

        TESTS::

            sage: SP = SuperPartitions(3)
            sage: TestSuite(SP).run()
        """
        self.n = n
        SuperPartitions.__init__(self, False)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(SuperPartitions(3))
            'SuperPartitions of 3'
        """
        return "SuperPartitions of %s"%self.n

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SuperPartitions(7)([[],[3,3,1]]) in SuperPartitions()
            True
            sage: SuperPartitions()([[],[3,3,1]]) in SuperPartitions(7)
            True
            sage: [[],[]] in SuperPartitions(0)
            True
            sage: [[0],[]] in SuperPartitions(0)
            True
            sage: [0] in SuperPartitions(0)
            True
            sage: [] in SuperPartitions(0)
            True
            sage: [1] in SuperPartitions(0)
            False
        """
        if x in SuperPartitions():
            if len(x)==0:
                return self.n==0
            if isinstance(x[0], builtins.list):
                n = sum(x[0]+x[1])
            else:
                n = sum(abs(a) for a in x)
            return n==self.n
        else:
            return False

    def __iter__(self):
        r"""
        An iterator for SuperPartitions of degree ``n``

        EXAMPLES::

            sage: SuperPartitions(1).list()
            [[; 1], [1; ], [0; 1], [1, 0; ]]
            sage: SuperPartitions(6).cardinality()
            80
        """
        m = 0
        while self.n >= m*(m-1)/2:
            for LA in SuperPartitions(self.n, m):
                yield self.element_class(self, LA)
            m+=1

class SuperPartitions_all(SuperPartitions):
    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: SP = SuperPartitions()
            sage: TestSuite(SP).run()
        """
        SuperPartitions.__init__(self, True)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(SuperPartitions())
            'SuperPartitions'
        """
        return "SuperPartitions"

    def __iter__(self):
        """
        Iterate over all SuperPartitions.

        EXAMPLES::

            sage: SP = SuperPartitions()
            sage: it = SP.__iter__()
            sage: [next(it) for i in range(6)]
            [[; ], [0; ], [; 1], [1; ], [0; 1], [1, 0; ]]
        """
        n = 0
        while True:
            for sp in SuperPartitions(n):
                yield self.element_class(self, list(sp))
            n += 1
