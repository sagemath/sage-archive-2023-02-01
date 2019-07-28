r"""
Super Tableaux

AUTHORS:

- Chaman Agrawal (2019-07-23): Initial version
"""

from __future__ import print_function, absolute_import

from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.structure.parent import Parent
from sage.arith.all import factorial
from sage.rings.integer import Integer
from sage.misc.all import prod
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.shifted_primed_tableau import PrimedEntry
from sage.combinat.tableau import Tableau, Tableaux, SemistandardTableaux


class SemistandardSuperTableau(Tableau):
    """
    A class to model a semistandard super tableau.

    INPUT:

    - ``t`` -- a tableau, a list of iterables, or an empty list

    OUTPUT:

    - A SemistandardSuperTableau object constructed from ``t``.

    A semistandard super tableau is a tableau whose entries are positive 
    integers, with even and odd parity and are weakly increasing in rows 
    and strictly increasing down columns.

    EXAMPLES::

        sage: from sage.combinat.super_tableau import SemistandardSuperTableau
        sage: t = SemistandardSuperTableau([['1p',2,"3'"],[2,3]]); t
        [[1', 2, 3'], [2, 3]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty printing
        1' 2 3'
        2 3
        sage: t = Tableau([["1p",2],[2]])
        sage: s = SemistandardSuperTableau(t); s
        [[1', 2], [2]]
        sage: SemistandardSuperTableau([]) # The empty tableau
        []

    TESTS::

        sage: from sage.combinat.super_tableau import SemistandardSuperTableaux
        sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
        sage: from sage.combinat.tableau import Tableaux
        sage: t = Tableaux()([[1,1],[2]])
        sage: s = SemistandardSuperTableaux()([[PrimedEntry(1),PrimedEntry(1)],
        ....:                                   [PrimedEntry(2)]])
        sage: s == t
        True
        sage: s.parent()
        Semistandard super tableaux
        sage: r = SemistandardSuperTableaux()(t); r.parent()
        Semistandard super tableaux
        sage: isinstance(r, Tableau)
        True
        sage: s2 = SemistandardSuperTableaux()([[PrimedEntry(1),PrimedEntry(1)],[
        ....:                                   PrimedEntry(2)]])
        sage: s2 == s
        True
        sage: s2.parent()
        Semistandard super tableaux
    """
    @staticmethod
    def __classcall_private__(cls, t):
        r"""
        This ensures that a SemistandardSuperTableau is only ever constructed 
        as an element_class call of an appropriate parent.

        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableau, SemistandardSuperTableaux
            sage: t = SemistandardSuperTableau([[1,1],[2]])
            sage: TestSuite(t).run()
            sage: t.parent()
            Semistandard super tableaux
            sage: t.category()
            Category of elements of Semistandard super tableaux
            sage: type(t)
            <class 'sage.combinat.super_tableau.SemistandardSuperTableaux_all_with_category.element_class'>
        """
        if isinstance(t, cls):
            return t

        # We must verify ``t`` is a list of iterables, and also
        # normalize it to be a list of tuples.
        try:
            t = [tuple(_) for _ in t]
        except TypeError:
            raise ValueError("A tableau must be a list of iterables.")
        return SemistandardSuperTableaux_all().element_class(SemistandardSuperTableaux_all(), t)

    def __init__(self, parent, t, check=True, preprocessed=False):
        r"""
        Initialize a semistandard super tableau

        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableau, SemistandardSuperTableaux
            sage: s = SemistandardSuperTableau([[1,"2'","3'",3], [2,"3'"]])
            sage: t = SemistandardSuperTableaux()([[1,"2p","3p",3], [2,"3p"]])
            sage: s == t
            True
            sage: t.parent()
            Semistandard super tableaux
            sage: s.parent()
            Semistandard super tableaux
            sage: r = SemistandardSuperTableaux()(s); r.parent()
            Semistandard super tableaux
            sage: s is t  # identical Semistandard super tableaux are distinct objects
            False
        """
        if not preprocessed:
            t = self._preprocess(t)
        super(SemistandardSuperTableau, self).__init__(parent, t, check=check)

    @staticmethod
    def _preprocess(t):
        """
        Preprocessing list ``T`` to initialize the tableau.
        The output is a list of rows as tuples, with entries being
        ``PrimedEntry``s.

        Trailing empty rows are removed.

        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableau
            sage: SemistandardSuperTableau._preprocess([["2'", "3p", 3.5]])
            [[2', 3', 4']]
            sage: SemistandardSuperTableau._preprocess([[None]])
            []
            sage: SemistandardSuperTableau._preprocess([])
            []
        """
        if isinstance(t, SemistandardSuperTableau):
            return t
        # Preprocessing list t for primes and other symbols
        t = [[PrimedEntry(entry) for entry in row if entry is not None]
             for row in t]
        while len(t) > 0 and len(t[-1]) == 0:
            t = t[:-1]
        return t
    
    def check(self):
        """
        Check that ``self`` is a valid semistandard super tableau.

        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableau
            sage: SemistandardSuperTableau([[1,2,3],[1]])
            Traceback (most recent call last):
            ...
            ValueError: the unprimed entries of each column must 
            be strictly increasing

            sage: SemistandardSuperTableau([[1,2,1]])
            Traceback (most recent call last):
            ...
            ValueError: the entries in each row of a semistandard super 
            tableau must be weakly increasing

            sage: SemistandardSuperTableau([[0,1]])
            Traceback (most recent call last):
            ...
            ValueError: the entries of a semistandard super tableau must be 
            non-negative primed integers
        """
        super(SemistandardSuperTableau, self).check()

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows

        for row in self:
            if not all(isinstance(c, PrimedEntry) and c > 0 for c in row):
                raise ValueError("""the entries of a semistandard super tableau
                                     must be non-negative primed integers""")
            if any(row[c] > row[c+1] for c in range(len(row)-1)):
                raise ValueError("""the entries in each row of a semistandard 
                                    super tableau must be weakly increasing""")

        if self:
            for row, next in zip(self, self[1:]):
                # Check that letters are weakly increasing down columns
                if any(row[c] > next[c] for c in range(len(next))):
                    raise ValueError("""the entries of each column of a 
                       semistandard super tableau must be weakly increasing""")
                # Check that unprimed letters are column strict
                if not all(row[c] < next[c] 
                        for c in range(len(next)) 
                        if (row[c].is_unprimed() or next[c].is_unprimed())):
                    raise ValueError("""the unprimed entries of each column 
                                            must be strictly increasing""")

            # Check that primed letters are row strict
            for row in self:
                if not all(row[c] < row[c+1] 
                        for c in range(len(row)-1) 
                        if (row[c].is_primed() or row[c+1].is_primed())):
                    raise ValueError("""the primed entries in each row must be 
                                        strictly increasing""")


class StandardSuperTableau(SemistandardSuperTableau):
    r"""
        A class to model a standard super tableau.

    INPUT:

      - ``t`` -- a Tableau, a list of iterables, or an empty list

    A standard super tableau is a semistandard super tableau whose entries 
    are in bijection with positive primed integers `1', 1, 2' \ldots n`.

    EXAMPLES::

        sage: from sage.combinat.super_tableau import StandardSuperTableau
        sage: t = StandardSuperTableau([["1'",1,"2'",2,"3'"],[3,"4'"]]); t
        [[1', 1, 2', 2, 3'], [3, 4']]
        sage: t.shape()
        [5, 2]
        sage: t.pp() # pretty printing
        1' 1 2' 2 3'
        3 4'
        sage: t.is_standard()
        True
        sage: StandardSuperTableau([]) # The empty tableau
        []

    When using code that will generate a lot of tableaux, it is more
    efficient to construct a StandardSuperTableau from the appropriate
    :class:`Parent` object::

        sage: from sage.combinat.super_tableau import StandardSuperTableaux
        sage: ST = StandardSuperTableaux()
        sage: ST([["1'",1,"2'",2,"3'"],[3,"4'"]])
        [[1', 1, 2', 2, 3'], [3, 4']]

    TESTS::

        sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
        sage: t = Tableaux()([[PrimedEntry('1p'), PrimedEntry('2p')], 
        ....:                [PrimedEntry(1)]])
        sage: s = StandardSuperTableaux()([['1p','2p'],[1]])
        sage: s == t
        True
        sage: s.parent()
        Standard super tableaux
        sage: r = StandardSuperTableaux()([]); r.parent()
        Standard super tableaux
        sage: isinstance(r, Tableau)
        True
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a :class:`StandardSuperTableau` is only ever 
        constructed as an ``element_class`` call of an appropriate parent.

        TESTS::

            sage: from sage.combinat.super_tableau import StandardSuperTableau
            sage: t = StandardSuperTableau([['1p','2p'],[1]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Standard super tableaux
            sage: type(t)
            <class 'sage.combinat.super_tableau.StandardSuperTableaux_all_with_category.element_class'>
        """
        if isinstance(t, StandardSuperTableau):
            return t

        return StandardSuperTableaux_all().element_class(StandardSuperTableaux_all(), t)

    def check(self):
        r"""
        Check that ``self`` is a standard tableau.

        TESTS::

            sage: from sage.combinat.super_tableau import StandardSuperTableau
            sage: StandardSuperTableau([[1,2,3],[4,5]])
            Traceback (most recent call last):
            ...
            ValueError: the entries in a standard tableau must be in 
            bijection with 1',1,2',2,...,n

            sage: StandardSuperTableau([[1,3,2]])
            Traceback (most recent call last):
            ...
            ValueError: the entries in each row of a semistandard super 
            tableau must be weakly increasing
        """
        super(StandardSuperTableau, self).check()
        # t is semistandard so we only need to check
        # that its entries are in bijection with {1',1,2', 2, ..., n}
        flattened_list = [i for row in self for i in row]
        a = PrimedEntry('1p')
        primed_list = []
        for i in range(len(flattened_list)):
            primed_list.append(a)
            a = a.increase_half()

        if sorted(flattened_list) != primed_list:
            raise ValueError("""the entries in a standard tableau must be in 
                                bijection with 1',1,2',2,...,n""")

    def is_standard(self):
        """
        Return ``True`` since ``self`` is a standard super tableau.

        EXAMPLES::

            sage: from sage.combinat.super_tableau import StandardSuperTableau
            sage: StandardSuperTableau([['1p', 1], ['2p', 2]]).is_standard()
            True
        """
        return True


################################
# Semi-standard Super tableaux #
################################
class SemistandardSuperTableaux(SemistandardTableaux):
    r"""
    The set of semistandard super tableaux.

    OUTPUT:

    - the class of all semistandard super tableaux

    A semistandard super tableau is a tableau whose entries are primed positive 
    integers, which are weakly increasing in rows and down columns as a whole 
    while the primed letters are row strict and the unprimed letters are column 
    strict. Note that Sage uses the English convention for partitions and 
    tableaux; the longer rows are displayed on top.

    EXAMPLES::

        sage: from sage.combinat.super_tableau import SemistandardSuperTableaux
        sage: SST = SemistandardSuperTableaux(); SST
        Semistandard super tableaux
    """
    @staticmethod
    def __classcall_private__(cls):
        r"""
        Normalize and process input to return the correct parent and
        ensure a unique representation.

        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableaux
            sage: SemistandardSuperTableaux()
            Semistandard super tableaux
        """
        return SemistandardSuperTableaux_all()

    Element = SemistandardSuperTableau

    def __contains__(self, x):
        """
        Return ``True`` if ``t`` can be interpreted as a
        :class:`SemistandardSuperTableau`.

        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableau
            sage: from sage.combinat.super_tableau import StandardSuperTableau
            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: T = sage.combinat.super_tableau.SemistandardSuperTableaux_all()
            sage: [[1,2],[2]] in T
            True
            sage: [[PrimedEntry('1p'),PrimedEntry(2)],[PrimedEntry(2)]] in T
            True
            sage: [] in T
            True
            sage: Tableau([[1]]) in T
            True
            sage: StandardSuperTableau([[1]]) in T
            Traceback (most recent call last):
            ...
            ValueError: the entries in a standard tableau must be in bijection 
            with 1',1,2',2,...,n

            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            True
            sage: [[1,3,2]] in T
            False
        """
        if isinstance(x, SemistandardSuperTableau):
            return True
        elif Tableaux.__contains__(self, x):
            x = SemistandardSuperTableau._preprocess(x)
            for row in x:
                if any(row[c] > row[c+1] for c in range(len(row)-1)):
                    return False
                if not all(row[c] < row[c+1] 
                        for c in range(len(row)-1) 
                        if (row[c].is_primed() or row[c+1].is_primed())):
                    return False
            for row, next in zip(x, x[1:]):
                if any(row[c] > next[c] for c in range(len(next))):
                    return False
                if not all(row[c] < next[c] 
                        for c in range(len(next)) 
                        if (row[c].is_unprimed() or next[c].is_unprimed())):
                    return False
            return True
        else:
            return False


class SemistandardSuperTableaux_all(SemistandardSuperTableaux):
    """
    All semistandard super tableaux.
    """
    def __init__(self):
        r"""
        Initializes the class of all semistandard super tableaux.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardSuperTableaux` 
            to ensure the options are properly parsed.

        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableaux_all
            sage: SST = SemistandardSuperTableaux_all(); SST
            Semistandard super tableaux
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())
        SemistandardSuperTableaux.__init__(self)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.super_tableau import SemistandardSuperTableaux
            sage: repr(SemistandardSuperTableaux())
            'Semistandard super tableaux'
        """
        return "Semistandard super tableaux"


################################
# Standard Super tableaux #
################################
class StandardSuperTableaux(SemistandardSuperTableaux, Parent):
    r"""
    The set of standard super tableaux.

    INPUT:

    - Either a non-negative integer (possibly specified with the keyword ``n``)
      or a partition.
    
    OUTPUT:

    - The class of all standard super tableaux

    - With a non-negative integer argument, ``n``, the class of all standard
      super tableaux of size ``n``

    - With a partition argument, the class of all standard super tableaux of 
      that shape.

    A standard super tableau is a tableau whose entries are primed positive 
    integers, which are strictly increasing in rows and down columns and 
    contains each letters from 1',1,2'...n exactly once.

    EXAMPLES::
        sage: from sage.combinat.super_tableau import StandardSuperTableaux
        sage: SST = StandardSuperTableaux(); SST
        Standard super tableaux
        sage: SST = StandardSuperTableaux(3); SST
        Standard super tableaux of size 3
        sage: SST.first()
        [[1', 1, 2']]
        sage: SST.last()
        [[1'], [1], [2']]
        sage: SST.cardinality()
        4
        sage: SST.list()
        [[[1', 1, 2']], [[1', 2'], [1]], [[1', 1], [2']], [[1'], [1], [2']]]

    TESTS::

        sage: StandardSuperTableaux()([])
        []
        sage: SST = StandardSuperTableaux([3,2]); SST
        Standard super tableaux of shape [3, 2]
        sage: SST.first()
        [[1', 2', 3'], [1, 2]]
        sage: SST.last()
        [[1', 1, 2'], [2, 3']]
        sage: SST.cardinality()
        5
        sage: SST.list()
        [[[1', 2', 3'], [1, 2]],
         [[1', 1, 3'], [2', 2]],
         [[1', 2', 2], [1, 3']],
         [[1', 1, 2], [2', 3']],
         [[1', 1, 2'], [2, 3']]]
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This class returns the appropriate parent based on arguments. 
        See the documentation for :class:`StandardTableaux` for more 
        information.

        TESTS::
            sage: from sage.combinat.super_tableau import StandardSuperTableaux
            sage: SST = StandardSuperTableaux(); SST
            Standard super tableaux
            sage: StandardSuperTableaux(3)
            Standard super tableaux of size 3
            sage: StandardSuperTableaux([2,2])
            Standard super tableaux of shape [2, 2]
            sage: StandardSuperTableaux(-1)
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a non-negative integer or a partition
            sage: StandardSuperTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a non-negative integer or a partition
        """
        from sage.combinat.partition import _Partitions
        from sage.combinat.skew_partition import SkewPartitions

        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs['n']
        else:
            n = None

        if n is None:
            return StandardSuperTableaux_all()

        elif n in _Partitions:
            return StandardSuperTableaux_shape(_Partitions(n))

        elif n in SkewPartitions():
            raise NotImplementedError("""Standard super tableau for skew partitions 
                                        is not implemented yet""")

        if not isinstance(n, (int, Integer)) or n < 0:
            raise ValueError("""the argument must be a non-negative integer 
                                or a partition""")

        return StandardSuperTableaux_size(n)

    Element = StandardSuperTableau

    def __contains__(self, x):
        """
        Return ``True`` if ``t`` can be interpreted as a
        :class:`StandardSuperTableau`.

        TESTS::

            sage: from sage.combinat.super_tableau import StandardSuperTableau
            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: T = sage.combinat.super_tableau.StandardSuperTableaux_all()
            sage: [[0.5,1],[1.5]] in T
            True
            sage: [[PrimedEntry('1p'),PrimedEntry('2p')],[PrimedEntry(1)]] in T
            True
            sage: [] in T
            True
            sage: Tableau([['1p']]) in T
            True
            sage: StandardSuperTableau([['1p']]) in T
            True

            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            False
            sage: [[1,3,2]] in T
            False
        """
        if isinstance(x, StandardSuperTableau):
            return True
        elif Tableaux.__contains__(self, x):
            x = SemistandardSuperTableau._preprocess(x)
            flattened_list = [i for row in x for i in row]
            a = PrimedEntry('1p')
            primed_list = []
            for i in range(len(flattened_list)):
                primed_list.append(a)
                a = a.increase_half()
            # return True
            return sorted(flattened_list) == primed_list and (len(x) == 0 or
                     (all(row[i]<row[i+1] for row in x for i in range(len(row)-1)) and
                       all(x[r][c]<x[r+1][c] for r in range(len(x)-1)
                                              for c in range(len(x[r+1])) )
                     ))
        else:
            return False


class StandardSuperTableaux_all(StandardSuperTableaux):
    """
    All standard super tableaux.
    """
    def __init__(self):
        r"""
        Initializes the class of all standard super tableaux.

        .. WARNING::

            Input is not checked; please use :class:`StandardSuperTableaux` 
            to ensure the options are properly parsed.

        TESTS::

            sage: from sage.combinat.super_tableau import StandardSuperTableaux_all
            sage: SST = StandardSuperTableaux_all(); SST
            Standard super tableaux
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())
        StandardSuperTableaux.__init__(self)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.super_tableau import StandardSuperTableaux
            sage: repr(StandardSuperTableaux())
            'Standard super tableaux'
        """
        return "Standard super tableaux"


class StandardSuperTableaux_size(StandardSuperTableaux, DisjointUnionEnumeratedSets):
    """
    Standard super tableaux of fixed size `n`.

    EXAMPLES::

        sage: [ t for t in StandardSuperTableaux(1) ]
        [[[1']]]
        sage: [ t for t in StandardSuperTableaux(2) ]
        [[[1', 1]], [[1'], [1]]]
        sage: [ t for t in StandardSuperTableaux(3) ]
        [[[1', 1, 2']], [[1', 2'], [1]], [[1', 1], [2']], [[1'], [1], [2']]]
        sage: StandardSuperTableaux(4)[:]
        [[[1', 1, 2', 2]],
         [[1', 2', 2], [1]],
         [[1', 1, 2], [2']],
         [[1', 1, 2'], [2]],
         [[1', 2'], [1, 2]],
         [[1', 1], [2', 2]],
         [[1', 2], [1], [2']],
         [[1', 2'], [1], [2]],
         [[1', 1], [2'], [2]],
         [[1'], [1], [2'], [2]]]
    """
    def __init__(self, n):
        r"""
        Initializes the class of all standard super tableaux of size ``n``.

        .. WARNING::

            Input is not checked; please use :class:`StandardSuperTableaux` to
            ensure the options are properly parsed.
        """
        StandardSuperTableaux.__init__(self)
        from sage.combinat.partition import Partitions_n
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(Partitions_n(n), StandardSuperTableaux_shape),
                                             category=FiniteEnumeratedSets(),
                                             facade=True, keepkey=False)
        self.size = Integer(n)

    def _repr_(self):
        """
        TESTS::

            sage: StandardSuperTableaux(3)
            Standard super tableaux of size 3
        """
        return "Standard super tableaux of size %s" % self.size

    def __contains__(self, x):
        """
        TESTS::

            sage: SST3 = StandardSuperTableaux(3)
            sage: all(st in SST3 for st in SST3)
            True
            sage: SST4 = StandardSuperTableaux(4)
            sage: [x for x in SST4 if x in SST3]
            []
            sage: 1 in StandardSuperTableaux(4)
            False
        """
        return StandardSuperTableaux.__contains__(self, x) and sum(map(len, x)) == self.size

    def cardinality(self):
        r"""
        Return the number of all standard super tableaux of size ``n``.

        The number of standard super tableaux of size `n` is equal to the
        number of standard tableaux of size `n` which is equal to the number 
        of involutions in the symmetric group `S_n`.

        ALGORITHM:

        The algorithm uses the fact that standard super tableaux are in 
        bijection with standard tableaux of size ``n`` which themselves are 
        in bijection with the involutions of size ``n``, (see page 41 in 
        section 4.1 of [Ful1997]_).  For each number of fixed points, you 
        count the number of ways to choose those fixed points multiplied 
        by the number of perfect matchings on the remaining values.

        EXAMPLES::

            sage: StandardSuperTableaux(3).cardinality()
            4
            sage: ns = [1,2,3,4,5,6]
            sage: sts = [StandardSuperTableaux(n) for n in ns]
            sage: all(st.cardinality() == len(st.list()) for st in sts)
            True
            sage: StandardSuperTableaux(50).cardinality()  # long time
            27886995605342342839104615869259776

        TESTS::

            sage: def cardinality_using_hook_formula(n):
            ....:     c = 0
            ....:     for p in Partitions(n):
            ....:         c += StandardSuperTableaux(p).cardinality()
            ....:     return c
            sage: all(cardinality_using_hook_formula(i) == 
            ....:       StandardSuperTableaux(i).cardinality() 
            ....:       for i in range(10))
            True
        """
        tableaux_number = self.size % 2  # identity involution
        fixed_point_numbers = list(range(tableaux_number, self.size + 1 - tableaux_number, 2))

        # number of involutions of size "size" (number of ways to
        # choose "fixed_point_number" out of "size" elements *
        # number of involutions without fixed point of size
        # "size" - "fixed_point_number")
        for fixed_point_number in fixed_point_numbers:
            tableaux_number += (self.size.binomial(fixed_point_number) *
                                prod(range(1, self.size - fixed_point_number, 2)))

        return tableaux_number


class StandardSuperTableaux_shape(StandardSuperTableaux):
    """
    Standard super tableaux of a fixed shape `p`.
    """
    def __init__(self, p):
        r"""
        Initializes the class of all standard super tableaux of a given shape.

        .. WARNING::

            Input is not checked; please use :class:`StandardSuperTableaux` to
            ensure the options are properly parsed.
        """
        super(StandardSuperTableaux_shape, self).__init__(category=FiniteEnumeratedSets())
        StandardSuperTableaux.__init__(self)
        self.shape = p

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: ST = StandardSuperTableaux([2,1,1])
            sage: all(st in ST for st in ST)
            True
            sage: len([x for x in StandardSuperTableaux(4) if x in ST])
            3
            sage: ST.cardinality()
            3
            sage: 1 in StandardSuperTableaux([2,1,1])
            False
        """
        return StandardSuperTableaux.__contains__(self, x) and [len(_) for _ in x] == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardSuperTableaux([2,1,1]))
            'Standard super tableaux of shape [2, 1, 1]'
        """
        return "Standard super tableaux of shape %s"%str(self.shape)

    def cardinality(self):
        r"""
        Return the number of standard super tableaux of given shape.

        This method uses the so-called *hook length formula*, a formula
        for the number of Young tableaux associated with a given
        partition. The formula says the following: Let `\lambda` be a
        partition. For each cell `c` of the Young diagram of `\lambda`,
        let the *hook length* of `c` be defined as `1` plus the number of
        cells horizontally to the right of `c` plus the number of cells
        vertically below `c`. The number of standard Young tableaux of
        shape `\lambda` is then `n!` divided by the product of the hook
        lengths of the shape of `\lambda`, where `n = |\lambda|`.

        For example, consider the partition ``[3,2,1]`` of ``6`` with
        Ferrers diagram::

            # # #
            # #
            #

        When we fill in the cells with their respective hook lengths, we
        obtain::

            5 3 1
            3 1
            1

        The hook length formula returns

        .. MATH::

            \frac{6!}{5 \cdot 3 \cdot 1 \cdot 3 \cdot 1 \cdot 1} = 16.

        EXAMPLES::

            sage: StandardSuperTableaux([3,2,1]).cardinality()
            16
            sage: StandardSuperTableaux([2,2]).cardinality()
            2
            sage: StandardSuperTableaux([5]).cardinality()
            1
            sage: StandardSuperTableaux([6,5,5,3]).cardinality()
            6651216
            sage: StandardSuperTableaux([]).cardinality()
            1

        REFERENCES:

        - http://mathworld.wolfram.com/HookLengthFormula.html
        """
        pi = self.shape

        number = factorial(sum(pi))
        hook = pi.hook_lengths()

        for row in hook:
            for col in row:
                #Divide the hook length by the entry
                number /= col

        return Integer(number)

    def __iter__(self):
        r"""
        An iterator for the standard super tableaux associated to the
        shape `p` of ``self``.

        EXAMPLES::

            sage: [t for t in StandardSuperTableaux([2,2])]
            [[[1', 2'], [1, 2]], [[1', 1], [2', 2]]]
            sage: [t for t in StandardSuperTableaux([3,2])]
            [[[1', 2', 3'], [1, 2]],
             [[1', 1, 3'], [2', 2]],
             [[1', 2', 2], [1, 3']],
             [[1', 1, 2], [2', 3']],
             [[1', 1, 2'], [2, 3']]]
            sage: st = StandardSuperTableaux([2,1])
            sage: st[0].parent() is st
            True
        """
        import copy
        pi = self.shape
        #Set the initial tableau by filling it in going down the columns
        tableau = [[None]*n for n in pi]
        size = sum(pi)
        row = 0
        col = 0
        for i in range(size):
            tableau[row][col] = i+1
            #If we can move down, then do it;
            #otherwise, move to the next column over
            if ( row + 1 < len(pi) and col < pi[row+1]):
                row += 1
            else:
                row = 0
                col += 1
        
        primedTableau = copy.deepcopy(tableau)
        for i, row in enumerate(primedTableau):
            for j, val in enumerate(row):
                primedTableau[i][j] = PrimedEntry(float(val)/2)
        yield self.element_class(self, primedTableau)

        # iterate until we reach the last tableau which is
        # filled with the row indices.
        last_tableau = sum([[row]*l for (row,l) in enumerate(pi)], [])

        #Convert the tableau to "vector format"
        #tableau_vector[i] is the row that number i
        #is in
        tableau_vector = [None]*size
        for row in range(len(pi)):
            for col in range(pi[row]):
                tableau_vector[tableau[row][col]-1] = row

        while tableau_vector!=last_tableau:
            #Locate the smallest integer j such that j is not
            #in the lowest corner of the subtableau T_j formed by
            #1,...,j.  This happens to be first j such that
            #tableau_vector[j]<tableau_vector[j-1].
            #l will correspond to the shape of T_j
            l = [0]*size
            l[0] = 1
            j = 0
            for i in range(1,size):
                l[tableau_vector[i]] += 1
                if ( tableau_vector[i] < tableau_vector[i-1] ):
                    j = i
                    break

            #Find the last nonzero row of l and store it in k
            i = size - 1
            while ( l[i] == 0 ):
                i -= 1
            k = i

            #Find a new row for the letter j (next lowest corner)
            t = l[ 1 + tableau_vector[j] ]
            i = k
            while ( l[i] != t ):
                i -= 1

            #Move the letter j to row i
            tableau_vector[j] = i
            l[i] -= 1

            #Fill in the columns of T_j using 1,...,j-1 in increasing order
            m = 0
            while ( m < j ):
                r = 0
                while ( l[r] != 0 ):
                    tableau_vector[m] = r
                    l[r] -= 1
                    m += 1
                    r += 1

            #Convert the tableau vector back to the regular tableau
            #format
            row_count= [0]*len(pi)
            tableau = [[None]*n for n in pi]

            for i in range(size):
                tableau[tableau_vector[i]][row_count[tableau_vector[i]]] = i+1
                row_count[tableau_vector[i]] += 1

            primedTableau = copy.deepcopy(tableau)
            for i, row in enumerate(primedTableau):
                for j, val in enumerate(row):
                    primedTableau[i][j] = PrimedEntry(float(val)/2)
            yield self.element_class(self, primedTableau)

        return

    def list(self):
        r"""
        Return a list of the standard super tableaux of the specified shape.

        EXAMPLES::

            sage: StandardSuperTableaux([2,2]).list()
            [[[1', 2'], [1, 2]], [[1', 1], [2', 2]]]
            sage: StandardSuperTableaux([5]).list()
            [[[1', 1, 2', 2, 3']]]
            sage: StandardSuperTableaux([3,2,1]).list()
            [[[1', 2, 3], [1, 3'], [2']],
             [[1', 2', 3], [1, 3'], [2]],
             [[1', 1, 3], [2', 3'], [2]],
             [[1', 2', 3], [1, 2], [3']],
             [[1', 1, 3], [2', 2], [3']],
             [[1', 2, 3'], [1, 3], [2']],
             [[1', 2', 3'], [1, 3], [2]],
             [[1', 1, 3'], [2', 3], [2]],
             [[1', 2', 2], [1, 3], [3']],
             [[1', 1, 2], [2', 3], [3']],
             [[1', 1, 2'], [2, 3], [3']],
             [[1', 2', 3'], [1, 2], [3]],
             [[1', 1, 3'], [2', 2], [3]],
             [[1', 2', 2], [1, 3'], [3]],
             [[1', 1, 2], [2', 3'], [3]],
             [[1', 1, 2'], [2, 3'], [3]]]
        """
        return [y for y in self]
