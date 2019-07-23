r"""
Super Tableaux

AUTHORS:

- Chaman Agrawal (2019-07-23): Initial version
"""

from __future__ import print_function, absolute_import
from six import add_metaclass

from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.structure.parent import Parent
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
            if not all(isinstance(c, (PrimedEntry)) and c > 0 for c in row):
                raise ValueError("the entries of a semistandard super tableau must be non-negative primed integers")
            if any(row[c] > row[c+1] for c in range(len(row)-1)):
                raise ValueError("the entries in each row of a semistandard super tableau must be weakly increasing")

        if self:
            for row, next in zip(self, self[1:]):
                # Check that letters are weakly increasing down columns
                if any(row[c] > next[c] for c in range(len(next))):
                    raise ValueError("the entries of each column of a semistandard super tableau must be weakly increasing")
                # Check that unprimed letters are column strict
                if not all(row[c] < next[c] for c in range(len(next)) if (row[c].is_unprimed() or next[c].is_unprimed())):
                    raise ValueError("the unprimed entries of each column must be strictly increasing")

            # Check that primed letters are row strict
            for row in self:
                if not all(row[c] < row[c+1] for c in range(len(row)-1) if (row[c].is_primed() or row[c+1].is_primed())):
                    raise ValueError("the primed entries in each row must be strictly increasing")


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
        This ensures that a :class:`StandardSuperTableau` is only ever constructed
        as an ``element_class`` call of an appropriate parent.

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
            ValueError: the entries in each row of a semistandard super tableau must be weakly increasing
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
            raise ValueError("the entries in a standard tableau must be in bijection with 1',1,2',2,...,n")

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
    A factory class for the various classes of semistandard super tableaux.

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
                if not all(row[c] < row[c+1] for c in range(len(row)-1) if (row[c].is_primed() or row[c+1].is_primed())):
                    return False
            for row, next in zip(x, x[1:]):
                if any(row[c] > next[c] for c in range(len(next))):
                    return False
                if not all(row[c] < next[c] for c in range(len(next)) if (row[c].is_unprimed() or next[c].is_unprimed())):
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
class StandardSuperTableaux(SemistandardSuperTableaux):
    r"""
    A factory class for the various classes of standard super tableaux.

    OUTPUT:

    - The class of all standard super tableaux

    A standard super tableau is a tableau whose entries are primed positive 
    integers, which are strictly increasing in rows and down columns and 
    contains each letters from 1',1,2'...n exactly once.

    EXAMPLES::
        sage: from sage.combinat.super_tableau import StandardSuperTableaux
        sage: SST = StandardSuperTableaux(); SST
        Standard super tableaux
    """
    @staticmethod
    def __classcall_private__(cls):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`StandardTableaux` for
        more information.

        TESTS::
            sage: from sage.combinat.super_tableau import StandardSuperTableaux
            sage: SST = StandardSuperTableaux(); SST
            Standard super tableaux
        """
        return StandardSuperTableaux_all()

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