from __future__ import print_function, absolute_import
from six import add_metaclass
from six.moves import range, zip, map

import numpy as np

from sage.combinat.partition import Partition
from sage.combinat.permutation import Permutation
from sage.combinat.posets.posets import Poset
from sage.combinat.tableau import Tableaux

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.functions.other import factorial
from sage.misc.cachefunc import cached_method
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.misc_c import prod
from sage.misc.prandom import randrange
from sage.rings.integer import Integer
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


@add_metaclass(InheritComparisonClasscallMetaclass)
class ShiftedPrimedTableau(ClonableArray):
    """
    A shifted primed tableau.

    EXAMPLES::

        sage: T = ShiftedPrimedTableaux([4,2])
        sage: T([[1,"2'","3'",3],[2,"3'"]])[1]
        [2, 3']
        sage: t = ShiftedPrimedTableau([[1,1.5,2.5,3],[0,2,2.5]])
        sage: t[1]
        [2, 3']
        sage: t = ShiftedPrimedTableau([[1,2,2.5,3],[0,2,2.5]])
        ValueError: not a primed tableau
    """
    
    @staticmethod
    def __classcall_private__(cls, t):
        
        """
        This ensures that a shifted tableau is only ever constructed
        as an ``element_class`` call of an appropriate parent.
        EXAMPLES::

            sage: data = [[1,2',3',3],[2,3']]
            sage: t = ShiftedPrimedTableau(data)
            sage: T = ShiftedPrimedTableaux([4,2])
            sage: t == T(data)
            True
        """
    
        if isinstance(t, cls):
            return t

        #Accounting for zeros at the beginning and at the end of a row
        i=0
        while i < len(t):
            row = t[i]
            try:
                while row[0]==0 or row[0]=='0':
                    row.pop(0)
                while row[-1]==0 or row[-1]=='0':
                    row.pop(-1)
            except IndexError:
                t.remove(t[i])
                continue
            t[i] = row
            i += 1

        shape = [len(row) for row in t]
        return ShiftedPrimedTableaux(shape = shape)(t)
    

    def __init__(self,parent, t):
        """
        Preprocessing of list t and initializing a shifted tableau.

        TESTS::

            sage: s = ShiftedPrimedTableau([[1,"2'","3'",3],[2,"3'"]])
            sage: t = ShiftedPrimedTableaux([4,2])([[1,"2p","3p",3],[0, 2,"3p"]])
            sage: s==t
            True
            sage: t.parent()
            Shifted primed tableaux of shape [4, 2]
            sage: s.parent()
            Shifted primed tableaux of shape [4, 2]
            sage: r = ShiftedPrimedTableaux([4, 2])(s); r.parent()
            Shifted primed tableaux of shape [4, 2]
            sage: s is t # identical shifted tableaux are distinct objects
            False

        A shifted primed tableau is shallowly immutable. The entries
        themselves may be mutable objects, though in that case the
        resulting ShiftedPrimedTableau should be unhashable.

            sage: t = ShiftedPrimedTableau([[1,"2p","3p",3],[2,"3p"]])
            sage: t0 = t[0]
            sage: t0[1] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
            sage: t[0][1] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
        """
        if isinstance(t, np.ndarray):
            t = t.tolist()
        
        # Preprocessing list t for primes and other symbols
        for i in range(len(t)):
            row = t[i]
            for j in range(len(row)):
                string = str(row[j])
                if string[-1] == "'" and string[:-1].isdigit() == True:
                    row[j] = float(float(string[:-1]) - .5)
                    continue
                if string[-1] == "p" and string[:-1].isdigit() == True:
                    row[j] = float(float(string[:-1]) - .5)
                    continue
                try:
                    row[j] = float(string)
                except ValueError:
                    row.pop(string)
            t[i] = row
                
        #Accounting for zeros at the beginning and at the end of a row
        i=0
        while i < len(t):
            row = t[i]
            try:
                while row[0]==0:
                    row.pop(0)
                while row[-1]==0:
                    row.pop(-1)
            except IndexError:
                t.remove(t[i])
                continue
            t[i] = row
            i += 1

        # Normalize t to be a list of tuples.
        t = [tuple(_) for _ in t]

        if isinstance(t, ShiftedPrimedTableau):
            # Since we are (supposed to be) immutable, we can share the underlying data
            ClonableArray.__init__(self, parent, t)
            return

        ClonableArray.__init__(self, parent, t)
        # This dispatches the input verification to the :meth:`check`
        # method.

    
    def to_matrix(self):
        array = []
        m = len(self[0])
        for i in range(len(self)):
            array.append([0]*i + self[i] + [0]*(m-i-len(self[i])))
        array = np.array(array)
        return array
    
    def check(self):
        """
        Check that ``self`` is a valid primed tableaux.

        EXAMPLES::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: t = T([[1,'2p',2,2],[2,'3p']])
            sage: t.check()
        """
        if [len(_) for _ in self] not in StrictPartitions():
            raise ValueError('shape must be a strict partition')
        for i,row in enumerate(self):
            if i > 0:
                if not all(val > self[i-1][j+1] for j,val in enumerate(row) if val in ZZ):
                    raise ValueError('column is not strictly increasing in non-primes')
                if not all(val >= self[i-1][j+1] for j,val in enumerate(row) if val not in ZZ):
                    raise ValueError('column is not weakly increasing in primes')
            if not all(row[j] <= row[j+1] for j in range(len(row)-1) if row[j] in ZZ):
                raise ValueError('row is not weakly increasing in non-primes')
            if not all(row[j] < row[j+1] for j in range(len(row)-1) if row[j] not in ZZ):
                raise ValueError('row is not strictly increasing in primes')
        if not all(row[0] in ZZ for row in self):
            raise ValueError('diagonal elements must be non-primes')

    
    def _repr_(self):
        
        string_list = ''
        for i,row in enumerate(self):
            string_list += '  '*i
            for val in row:
                if val in ZZ:
                    string_list += str(int(val))+' '
                else:
                    string_list += str(int(val+.5))+"'"
            string_list += '\n'
        
        return (string_list)
    
    def shape(self):
        return
    def weight(self):
        return
    def entry(self,position):
        return
    def reading_word(self):
        return


class ShiftedPrimedTableaux(UniqueRepresentation, Parent):
    r"""
    A factory for the various classes of shifted standard tableaux.

    INPUT:

    - a weight and/or a partition

    OUTPUT:

    - with no argument, the class of all primed tableaux

    - with a partition argument, the class of all primed tableaux of that
        shape (infinite set if we don't specify the weight)
        
    - with a weight argument, the class of all primed tableaux of that
        weight (finite set)
    

    A primed tableau is a shifted tableau on the alphabet 
        X' = {1' < 1 < 2' < 2 <...< n' < n} such that 
        1). the entries are weakly increasing along rows and columns
        2). a row can't have two repeated primed elements, and a column
            can't have two repeated non-primed elements
        3). there are only non-primed elements along the main diagonal
        
    The weight of a tableau is defined to be the vector with i-th component
        equal to the number of letters i and i' in the tableau.
    The sum of the entries in the weight vector must be equal to the number 
        of boxes in the partition.


    TESTS::

        sage: ShiftedPrimedTableaux([])
        []
        sage: SPT = ShiftedPrimedTableaux((1,2,2),[3,2]); SPT
        Shifted primed tableaux of shape [3, 2] and weight (1,2,2)
        sage: SPT.cardinality()
        2
        sage: SPT.list()
        [[[1, 2, 3], [4, 5]], [[1, 2, 4], [3, 5]]]
    """
    # use Tableaux options
    options = Tableaux.options

    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`ShiftedPrimedTableaux` for
        more information.

        TESTS::

            sage: ShiftedPrimedTableaux()
            Shifted primed tableaux
            sage: ShiftedPrimedTableaux(3)
            Shifted tableaux of size 3
            sage: ShiftedTableaux([2,1])
            Shifted tableaux of shape [2, 1]
            sage: ShiftedTableaux(0)
            Shifted tableaux of size 0

            sage: ShiftedTableaux(-1)
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a non-negative integer or a partition
            sage: ShiftedTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a non-negative integer or a partition
        """
        weight = None
        shape = None
        
        if 'size' in kwargs and isinstance(kwargs['size'],(list,tuple,Partition)):
            shape = Partition(kwargs['size'])
        
        if 'shape' in kwargs:
            shape = Partition(kwargs['shape'])

        if 'weight' in kwargs:
            weight = tuple(kwargs['weight'])

        if args:
            if isinstance(args[0],tuple) and weight==None:
                weight = args[0]
                if len(args)>1:
                    if (isinstance(args[1],(list,Partition)) or args[1]==None) and shape==None:
                        shape = args[1]
                    else:
                        raise ValueError('weight argument must be a tuple and shape argument must be a strictly increasing partition')
            
            elif isinstance(args[0],(list,Partition)) and shape==None:
                shape = args[0]
                if len(args)>1:
                    if (isinstance(args[1],tuple) or args[1]==None) and weight==None:
                        weight = args[1]
                    else:
                        raise ValueError('weight argument must be a tuple and shape argument must be a strictly increasing partition')
            else:
                raise ValueError('weight argument must be a tuple and shape argument must be a strictly increasing partition')
        
        if shape is not None:
            try:
                shape = Partition(shape)
            except ValueError:
                raise ValueError('{} is not a strict partition'.format(shape))
    
        if shape is None and weight is None:
            return ShiftedPrimedTableaux_all()
        

        elif weight is None and shape in StrictPartitions():
            return ShiftedPrimedTableaux_shape(Partition(shape))
        
        elif shape is None:
            return ShiftedPrimedTableaux_weight(weight)

        if shape not in StrictPartitions():
            raise ValueError("{} is not a strict partition".format(shape))
            
        if sum(shape) != sum(weight):
            raise ValueError("the sum of weights should match the sum of parts of the shape")

        return ShiftedPrimedTableaux_weight_shape(weight,shape)

    def __contains__(self, t):
        """
        EXAMPLES::

        """
        if isinstance(t, ShiftedPrimedTableau) or t == []:
            return True

        if isinstance(t, np.ndarray):
            t = t.tolist()
        
        # Preprocessing list t for primes and other symbols
        for i in range(len(t)):
            row = t[i]
            for j in range(len(row)):
                string = str(row[j])
                if string[-1] == "'" and string[:-1].isdigit() == True:
                    row[j] = float(float(string[:-1]) - .5)
                    continue
                if string[-1] == "p" and string[:-1].isdigit() == True:
                    row[j] = float(float(string[:-1]) - .5)
                    continue
                try:
                    row[j] = float(string)
                except ValueError:
                    row.pop(string)
            t[i] = row
                
        #Accounting for zeros at the beginning and at the end of a row
        i=0
        while i < len(t):
            row = t[i]
            try:
                while row[0]==0:
                    row.pop(0)
                while row[-1]==0:
                    row.pop(-1)
            except IndexError:
                t.remove(t[i])
                continue
            t[i] = row
            i += 1

        # Normalize t to be a list of tuples.
        t = [tuple(_) for _ in t]

        for i,row in enumerate(t):
            if i > 0:
                if not all(val > t[i-1][j+1] for j,val in enumerate(row) if val in ZZ):
                    return False
                if not all(val >= t[i-1][j+1] for j,val in enumerate(row) if val not in ZZ):
                    return False
            if not all(row[j] <= row[j+1] for j in range(len(row)-1) if row[j] in ZZ):
                return False
            if not all(row[j] < row[j+1] for j in range(len(row)-1) if row[j] not in ZZ):
                return False
        if not all(row[0] in ZZ for row in t):
            return False

        return True

    _is_a = __contains__

    def an_element(self):
        r"""
        Return a particular shifted tableaux in the class.

        TESTS::

            sage: ShiftedTableaux().an_element()
            []
            sage: ShiftedTableaux(4).an_element()
            [[1, 2, 3, 4]]
            sage: ShiftedTableaux([3,1]).an_element()
            [[1, 2, 3], [4]]
        """
        return self[0]

class ShiftedPrimedTableaux_shape(ShiftedPrimedTableaux):
    """
    Shifted tableaux of a fixed shape.
    """
    Element = ShiftedPrimedTableau

    def __init__(self, shape):
        r"""
        Initializes the class of semistandard tableaux of shape ``p`` and no
        maximum entry.

        TESTS::

            sage: TestSuite( ShiftedTableaux([3,2,1]) ).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._shape = shape

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedTableaux([3,2,1])    # indirect doctest
            Shifted tableaux of shape [3, 2, 1]
        """
        return "Shifted tableaux of shape {}".format(self._shape)

    def __contains__(self, x):
        """
        TESTS::

            sage: ST421 = ShiftedTableaux([4,2,1])
            sage: all([st in ST421 for st in ST421])
            True
            sage: ST42 = ShiftedTableaux([4,2])
            sage: filter(lambda x: x in ST42, ST421)
            []
        """
        if isinstance(x, ShiftedPrimedTableau):
            return [len(row) for row in x] == self._shape

        return (ShiftedPrimedTableaux.__contains__(self, x)
                and [len(row) for row in x] == self._shape)

    def _element_constructor_(self, t):
        r"""
        Constructs an object from ``t`` as an element of ``self``, if
        possible.

        INPUT:

        - ``t`` -- data which can be interpreted as a tableau

        OUTPUT:

        - the corresponding tableau object

        TESTS::

            sage: ShiftedTableaux([3])([[1,2,3]]).parent() is ShiftedTableaux([3])
            True
            sage: ShiftedTableaux([3])([[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 2]] is not an element of Shifted tableaux of shape [3]
        """
        if not t in self:
            raise ValueError("{} is not an element of {}".format(t, self))
        
        return self.element_class(self, t)

    def shape(self):
        """
        Return the shape of the shifted tableaux ``self``.

        EXAMPLES::

            sage: ShiftedTableaux([6,4,3,1]).shape()
            [6, 4, 3, 1]
        """
        return self._shape

    @cached_method
    def size(self):
        """
        Return the shape of the shifted tableaux ``self``.

        EXAMPLES::

            sage: ShiftedTableaux([6,4,3,1]).size()
            14
        """
        return self._shape.size()

######################
# Strict Partitions #
######################
class StrictPartitions(Partitions):
    r"""
    The class of **strict partitions**.

    A strict partition of an integer `n` is a :class:`Partition` with
    distinct parts.

    INPUT:

    - ``n`` -- a non-negative integer, the size of the strict partitions
    """
    @staticmethod
    def __classcall_private__(cls, size=None):
        """
        Normalize the input to ensure a unique representation.

        TESTS::

            sage: from sage.combinat.partition import StrictPartitions
            sage: P = StrictPartitions()
            sage: P12 = StrictPartitions(12)
            sage: P is P12
            False
        """
        if size is not None:
            return StrictPartitions_size( ZZ(size) )

        return StrictPartitions_all()

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is contained in ``self`` and ``False``
        otherwise.

        EXAMPLES::

            sage: from sage.combinat.partition import StrictPartitions
            sage: [5,2] in StrictPartitions()
            True
            sage: [5,2] in StrictPartitions(7)
            True
            sage: [5,2] in StrictPartitions(5)
            False
            sage: [5,2,2] in StrictPartitions(5)
            False
            sage: [5,2,0] in StrictPartitions(7)
            True
            sage: Partition([5,2]) in StrictPartitions()
            True
            sage: Partition([5,2]) in StrictPartitions(7)
            True
            sage: Partition([5,2]) in StrictPartitions(3)
            False
        """
        if not Partitions.__contains__(self, x):
            return False
        return len(x) == 0 or (x[-1] in NN and all(x[i]>x[i+1] for i in range(len(x)-1) if x[i]>0))

class StrictPartitions_all(StrictPartitions, DisjointUnionEnumeratedSets):
    """
    Class of all strict partitions.

    EXAMPLES::

        sage: from sage.combinat.partition import StrictPartitions
        sage: StrictPartitions()
        Strict Partitions
    """
    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.partition import StrictPartitions
            sage: TestSuite(StrictPartitions()).run()
        """
        # I'd rather use super() here but Partitions() complains
        DisjointUnionEnumeratedSets.__init__(self,
                family=Family(NonNegativeIntegers(), StrictPartitions_size),
                facade=True, keepkey=False
        )
        Partitions.__init__(self, is_infinite=True)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.combinat.partition import StrictPartitions
            sage: StrictPartitions(5)
            Strict Partitions of the integer 5
        """
        return "Strict Partitions"

class StrictPartitions_size(StrictPartitions):
    """
    Strict Partitions of the integer ``size``.

    TESTS::

        sage: TestSuite( sage.combinat.partition.StrictPartitions_size(0) ).run()
        sage: TestSuite( sage.combinat.partition.StrictPartitions_size(10) ).run()
    """

    def __init__(self, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.partition import StrictPartitions
            sage: StrictPartitions(3)
            Strict Partitions of the integer 3

        TESTS::

            sage: from sage.combinat.partition import StrictPartitions
            sage: TestSuite(StrictPartitions(9)).run()
        """
        Partitions.__init__(self, is_infinite=False)
        self.n = n # would be better to call this size, but for consistency...
        self.size = n

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.combinat.partition import StrictPartitions
            sage: StrictPartitions(5)
            Strict Partitions of the integer 5
        """
        return "Strict Partitions of the integer {}".format(self.n)

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is contained in ``self`` and ``False``
        otherwise.

        Examples::

            sage: from sage.combinat.partition import StrictPartitions
            sage: [5,2] in StrictPartitions(7)
            True
            sage: [5,2] in StrictPartitions(5)
            False
            sage: [5,2,2] in StrictPartitions(5)
            False
            sage: [5,2,0,0] in StrictPartitions(7)
            True
            sage: Partition([5,2]) in StrictPartitions(7)
            True
            sage: Partition([5,2,1]) in StrictPartitions(7)
            False
        """
        return StrictPartitions.__contains__(self, x) and sum(x) == self.size

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: from sage.combinat.partition import StrictPartitions
            sage: StrictPartitions(10)[:] #indirect doct test
            [[10],
             [9, 1],
             [8, 2],
             [7, 3],
             [7, 2, 1],
             [6, 4],
             [6, 3, 1],
             [5, 4, 1],
             [5, 3, 2],
             [4, 3, 2, 1]]
        """
        for p in self._fast_iterator(self.n, self.n+1):
            yield self.element_class(self, p)

    def _fast_iterator(self, size, max):
        """
        A fast (recursive) iterator which returns a list.

        This method is not intended to be called directy.

        INPUT:

        - ``size`` -- a positive integer, giving the size of the partitions

        - ``max`` -- a positive integer giving the maximu size of the parts of
          the partitions to be returned

        OUTPUT:

        - an iterator  of the strict partitions of size ``size``

        EXAMPLES::

            sage: from sage.combinat.partition import StrictPartitions
            sage: StrictPartitions(7)[:]  # indirect doc test
            [[7], [6, 1], [5, 2], [4, 3], [4, 2, 1]]
        """
        if size < max:
            yield [size]

        for m in reversed(range(1, min(size, max))):
            for mu in self._fast_iterator(size-m, m):
                yield [m] + mu
        return

    def an_element(self):
        """
        Returns a partition in ``self``.

        EXAMPLES::

            sage: from sage.combinat.partition import StrictPartitions
            sage: StrictPartitions(4).an_element()  # indirect doctest
            [3, 1]
            sage: StrictPartitions(0).an_element()
            []
            sage: StrictPartitions(1).an_element()
            [1]
        """
        if self.n == 0:
            elt = []
        elif self.n == 1:
            elt = [1]
        else:
            elt = [self.n-1, 1]
        return self.element_class(self, elt)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: from sage.combinat.partition import StrictPartitions
            sage: StrictPartitions(7).cardinality()
            5
        """
        return ZZ(len( [1 for p in self._fast_iterator(self.n, self.n+1)] ))

    def random_element(self, measure = 'uniform'):
        r"""
        Return a random strict partition.

        EXAMPLE:

            sage: from sage.combinat.partition import StrictPartitions
            sage: StrictPartitions(7).cardinality()  # random
        """
        from sage.misc.prandom import randrange
        partitions = self.list()
        return partitions[randrange(len(partitions))]
