r"""
Ribbon Shaped Tableaux
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.skew_tableau import SkewTableau, StandardSkewTableaux, from_expr
from sage.combinat.tableau import TableauOptions
from sage.combinat.permutation import Permutation, descents_composition_first, descents_composition_list, descents_composition_last
from sage.combinat.skew_partition import SkewPartition
from sage.rings.integer import Integer
from combinat import CombinatorialObject
from sage.combinat.words.words import Words
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

class RibbonShapedTableau(SkewTableau):
    r"""
    A ribbon.

    A ribbon is a skew tableau that does not contain a `2 \times 2` box. A
    ribbon is given by a list of the rows from top to bottom.

    EXAMPLES::

        sage: x = RibbonShapedTableau([[None, None, None, 2, 3], [None, 1, 4, 5], [3, 2]]); x
        [[None, None, None, 2, 3], [None, 1, 4, 5], [3, 2]]
        sage: x.pp()
          .  .  .  2  3
          .  1  4  5
          3  2
        sage: x.shape()
        [5, 4, 2] / [3, 1]

    The entries labeled by ``None`` correspond to the inner partition.
    Using ``None`` is optional, the entries will be shifted accordingly.  ::

        sage: x = RibbonShapedTableau([[2,3],[1,4,5],[3,2]]); x.pp()
          .  .  .  2  3
          .  1  4  5
          3  2
    """
    @staticmethod
    def __classcall_private__(cls, r):
        r"""
        Return a ribbon tableau object.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]])
            [[None, None, 2, 3], [1, 4, 5]]
        """
        if isinstance(r, list):
            if len(r) == 0:
                return StandardRibbonShapedTableaux()(r)
            if all(isinstance(i, list) for i in r):
                if all(all(j is None or (isinstance(j, (int, Integer)) and j>0) for j in i) for i in r):
                    return StandardRibbonShapedTableaux()(r)
        raise TypeError("r must be a list of positive integers")

    def __init__(self, parent, t):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R = RibbonShapedTableau([[2,3],[1,4,5]])
            sage: TestSuite(R).run()
        """
        if not isinstance(t, SkewTableau):
            #scrubing None
            t = [ [i for i in row if i is not None] for row in t]

            st = []
            space_count = 0
            for row in reversed(t):
                st.append( [None]*space_count + row )
                space_count += len(row) - 1
            st.reverse()
            t = st
        else:
            t = list(t)

        SkewTableau.__init__(self, parent, t)

    def __setstate__(self, state):
        r"""
        In order to maintain backwards compatibility and be able to unpickle
        a old pickle from ``Ribbon_class`` we have to override the
        default ``__setstate__``.

        EXAMPLES::

            sage: loads( 'x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+\xcaLJ\xca\xcf\xe3\n\x02S\xf1\xc99\x89\xc5\xc5\\\x85\x8c\x9a\x8d\x85L\xb5\x85\xcc\x1a\xa1\xac\xf1\x19\x89\xc5\x19\x85,~@VNfqI!kl!\x9bFl!\xbb\x06\xc4\x9c\xa2\xcc\xbc\xf4b\xbd\xcc\xbc\x92\xd4\xf4\xd4"\xae\xdc\xc4\xec\xd4x\x18\xa7\x90#\x94\xd1\xb05\xa8\x903\x03\xc80\x022\xb8Rc\x0b\xb95@<c \x8f\x07\xc40\x012xSSK\x93\xf4\x00l\x811\x17')
            [[1, 2], [3, 4]]
            sage: loads(dumps( RibbonShapedTableau([[3,2,1], [1,1]]) ))  # indirect doctest
            [[None, 3, 2, 1], [1, 1]]
        """
        if isinstance(state, dict):   # for old pickles from Ribbon_class
            self._set_parent(StandardRibbonShapedTableaux())
            self.__dict__ = state
        else:
            self._set_parent(state[0])
            self.__dict__ = state[1]

    def height(self):
        """
        Return the height of ``self``.

        The height is given by the number of rows in the outer partition.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]]).height()
            2
        """
        return len(self.outer_shape())

    def spin(self):
        """
        Return the spin of ``self``.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]]).spin()
            1/2
        """
        return Integer(self.height()-1)/2

    def width(self):
        """
        Return the width of the ribbon.

        This is given by the length of the longest row in the outer partition.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]]).width()
            4
            sage: RibbonShapedTableau([]).width()
            0
        """
        #return 1+sum([len(r)-1 for r in self])
        return len(self[0]) if len(self) > 0 else 0

    def to_skew_tableau(self):
        """
        This is deprecated in :trac:`14101` since :class:`RibbonShapedTableau`
        inherits from :class:`SkewTableau`.

        EXAMPLES::

            sage: RibbonShapedTableau([[2,3],[1,4,5]]).to_skew_tableau()
            doctest:...: DeprecationWarning: this method is deprecated since ribbon shaped tableaux are skew partitions
            See http://trac.sagemath.org/14101 for details.
            [[None, None, 2, 3], [1, 4, 5]]
        """
        from sage.misc.superseded import deprecation
        deprecation(14101,'this method is deprecated since ribbon shaped tableaux are skew partitions')
        return self

    def to_permutation(self):
        """
        Returns the permutation corresponding to the ribbon tableau.

        EXAMPLES::

            sage: r = RibbonShapedTableau([[1], [2,3], [4, 5, 6]])
            sage: r.to_permutation()
            [4, 5, 6, 2, 3, 1]
        """
        return Permutation(self.to_word())

    def evaluation(self):
        """
        Return the evaluation of the word from ribbon.

        EXAMPLES::

            sage: RibbonShapedTableau([[1,2],[3,4]]).evaluation()
            [1, 1, 1, 1]
        """
        ed = self.to_word().evaluation_dict()
        entries = ed.keys()
        m = max(entries) + 1 if entries else -1
        return [ed.get(k,0) for k in range(1,m)]


class StandardRibbonShapedTableaux(StandardSkewTableaux):
    """
    The set of all standard ribbon shaped tableaux.

    INPUT:

    - ``shape`` -- (optional) the composition shape of the rows
    """
    @staticmethod
    def __classcall_private__(cls, shape=None, **kwds):
        """
        Normalize input to ensure a unique representation and pick the correct
        class based on input.

        EXAMPLES::

            sage: S1 = StandardRibbonShapedTableaux([4, 2, 2, 1])
            sage: S2 = StandardRibbonShapedTableaux((4, 2, 2, 1))
            sage: S1 is S2
            True
        """
        if shape is not None:
            from sage.combinat.partition import Partition
            return StandardRibbonShapedTableaux_shape(Partition(shape))

        # Otherwise arg0 takes the place of the category in pickling
        return super(StandardRibbonShapedTableaux, cls).__classcall__(cls, **kwds)

    def __init__(self, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = StandardRibbonShapedTableaux()
            sage: TestSuite(S).run()
        """
        if category is None:
            category = InfiniteEnumeratedSets()

        StandardSkewTableaux.__init__(self, category=category)

    def __iter__(self):
        """
        Iterate through ``self``.

        EXAMPLES::

            sage: it = StandardRibbonShapedTableaux().__iter__()
            sage: [it.next() for x in range(10)]
            [[],
             [[1]],
             [[1, 2]],
             [[1], [2]],
             [[1, 2, 3]],
             [[None, 2], [1, 3]],
             [[None, 1], [2, 3]],
             [[1], [2], [3]],
             [[1, 2, 3, 4]],
             [[None, None, 3],
             [1, 2, 4]]]
        """
        from sage.combinat.partition import _Partitions
        for p in _Partitions:
            for r in StandardRibbonShapedTableaux_shape(p):
                yield self.element_class(self, r)

    Element = RibbonShapedTableau
    global_options = TableauOptions

    def from_shape_and_word(self, shape, word):
        """
        Return the ribbon corresponding to the given ribbon shape and word.

        EXAMPLES::

            sage: StandardRibbonShapedTableaux().from_shape_and_word([2,3],[1,2,3,4,5])
            [[None, None, 1, 2], [3, 4, 5]]
        """
        pos = 0
        r = []
        for l in shape:
            r.append(word[pos:pos+l])
            pos += l
        return self.element_class(self, r)

    def from_permutation(self, p):
        """
        Return a standard ribbon of size ``len(p)`` from a permutation ``p``. The
        lengths of each row are given by the distance between the descents
        of the permutation ``p``.

        EXAMPLES::

            sage: import sage.combinat.ribbon_shaped_tableau as rst
            sage: [StandardRibbonShapedTableaux().from_permutation(p) for p in Permutations(3)]
            [[[1, 2, 3]],
             [[None, 2], [1, 3]],
             [[1, 3], [2]],
             [[None, 1], [2, 3]],
             [[1, 2], [3]],
             [[1], [2], [3]]]
        """
        if p == []:
            return self.element_class(self, [])

        comp = p.descents()

        if comp == []:
            return self.element_class(self, [p[:]])

        #[p[j]$j=compo[i]+1..compo[i+1]] $i=1..nops(compo)-1, [p[j]$j=compo[nops(compo)]+1..nops(p)]
        r = []
        r.append([p[j] for j in range(comp[0]+1)])
        for i in range(len(comp)-1):
            r.append([ p[j] for j in range(comp[i]+1,comp[i+1]+1) ])
        r.append( [ p[j] for j in range(comp[-1]+1, len(p))] )
        r.reverse()
        return self.element_class(self, r)

class StandardRibbonShapedTableaux_shape(StandardRibbonShapedTableaux):
    """
    Class of standard ribbon tableaux of ribbon shape ``shape``.

    EXAMPLES::

        sage: StandardRibbonShapedTableaux([2,2])
        Standard ribbon tableaux of shape [2, 2]
        sage: StandardRibbonShapedTableaux([2,2]).first()
        [[None, 2, 4], [1, 3]]
        sage: StandardRibbonShapedTableaux([2,2]).last()
        [[None, 1, 2], [3, 4]]
        sage: StandardRibbonShapedTableaux([2,2]).cardinality()
        5
        sage: StandardRibbonShapedTableaux([2,2]).list()
        [[[None, 2, 4], [1, 3]],
         [[None, 2, 3], [1, 4]],
         [[None, 1, 4], [2, 3]],
         [[None, 1, 3], [2, 4]],
         [[None, 1, 2], [3, 4]]]
        sage: StandardRibbonShapedTableaux([3,2,2]).cardinality()
        155
    """
    @staticmethod
    def __classcall_private__(cls, shape):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = StandardRibbonShapedTableaux([2,2])
            sage: S2 = StandardRibbonShapedTableaux((2,2))
            sage: S is S2
            True
        """
        return super(StandardRibbonShapedTableaux, cls).__classcall__(cls, tuple(shape))

    def __init__(self, shape):
        """
        TESTS::

            sage: S = StandardRibbonShapedTableaux([2,2])
            sage: TestSuite(S).run()
        """
        self.shape = shape
        StandardRibbonShapedTableaux.__init__(self, FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: StandardRibbonShapedTableaux([2,2])
            Standard ribbon tableaux of shape [2, 2]
        """
        return "Standard ribbon tableaux of shape %s"%list(self.shape)

    def first(self):
        """
        Return the first standard ribbon of ``self``.

        EXAMPLES::

            sage: StandardRibbonShapedTableaux([2,2]).first()
            [[None, 2, 4], [1, 3]]
        """
        return self.from_permutation(descents_composition_first(self.shape))

    def last(self):
        """
        Return the last standard ribbon of ``self``.

        EXAMPLES::

            sage: StandardRibbonShapedTableaux([2,2]).last()
            [[None, 1, 2], [3, 4]]
        """
        return self.from_permutation(descents_composition_last(self.shape))

    def __iter__(self):
        """
        An iterator for the standard ribbon of ``self``.

        EXAMPLES::

            sage: [t for t in StandardRibbonShapedTableaux([2,2])]
            [[[None, 2, 4], [1, 3]],
             [[None, 2, 3], [1, 4]],
             [[None, 1, 4], [2, 3]],
             [[None, 1, 3], [2, 4]],
             [[None, 1, 2], [3, 4]]]
        """
        for p in descents_composition_list(self.shape):
            yield self.from_permutation(p)

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.ribbon', 'Ribbon_class', RibbonShapedTableau)

