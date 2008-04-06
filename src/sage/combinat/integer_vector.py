"""
Integer vectors
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

from combinat import CombinatorialClass
from __builtin__ import list as builtinlist
from sage.rings.integer import Integer
from sage.rings.arith import binomial
import misc
from sage.rings.infinity import PlusInfinity
import integer_list
import cartesian_product
import functools


def _default_function(l, default, i):
    """
    EXAMPLES:
        sage: from sage.combinat.integer_vector import _default_function
        sage: import functools
        sage: f = functools.partial(_default_function, [1,2,3], 99)
        sage: f(0)
        99
        sage: f(1)
        1
        sage: f(2)
        2
        sage: f(3)
        3
        sage: f(4)
        99
    """
    try:
        if i <= 0:
            return default
        return l[i-1]
    except IndexError:
        return default

infinity = PlusInfinity()
def list2func(l, default=None):
    """
    Given a list l, return a function that takes in a value
    i and return l[i-1].  If default is not None, then the function
    will return the default value for out of range i's.

    EXAMPLES:
        sage: f = sage.combinat.integer_vector.list2func([1,2,3])
        sage: f(1)
        1
        sage: f(2)
        2
        sage: f(3)
        3
        sage: f(4)
        Traceback (most recent call last):
        ...
        IndexError: list index out of range

        sage: f = sage.combinat.integer_vector.list2func([1,2,3], 0)
        sage: f(3)
        3
        sage: f(4)
        0
    """
    if default is None:
        return lambda i: l[i-1]
    else:
        return functools.partial(_default_function, l, default)


def constant_func(i):
    """
    Returns the constant function i.

    EXAMPLES:
        sage: f = sage.combinat.integer_vector.constant_func(3)
        sage: f(-1)
        3
        sage: f('asf')
        3
    """
    return lambda x: i

def IntegerVectors(n=None, k=None, **kwargs):
    """
    Returns the combinatorial class of integer vectors.

    EXAMPLES:
      If n is not specified, it returns the class of all
      integer vectors.

        sage: IntegerVectors()
        Integer vectors
        sage: [] in IntegerVectors()
        True
        sage: [1,2,1] in IntegerVectors()
        True
        sage: [1, 0, 0] in IntegerVectors()
        True

      If n is specified, then it returns the class of all
      integer vectors which sum to n.

        sage: IV3 = IntegerVectors(3); IV3
        Integer vectors that sum to 3

      Note that trailing zeros are ignored so that [3, 0]
      does not show up in the following list (since [3] does)

        sage: IntegerVectors(3, max_length=2).list()
        [[3], [2, 1], [1, 2], [0, 3]]

      If n and k are both specified, then it returns the class
      of integer vectors that sum to n and are of length k.

        sage: IV53 = IntegerVectors(5,3); IV53
        Integer vectors of length 3 that sum to 5
        sage: IV53.count()
        21
        sage: IV53.first()
        [5, 0, 0]
        sage: IV53.last()
        [0, 0, 5]
        sage: IV53.random() #random
        [0, 1, 4]

    """
    if n is None:
        return IntegerVectors_all()
    elif k is None:
        return IntegerVectors_nconstraints(n,kwargs)
    else:
        if isinstance(k, builtinlist):
            return IntegerVectors_nnondescents(n,k)
        else:
            if len(kwargs) == 0:
                return IntegerVectors_nk(n,k)
            else:
                return IntegerVectors_nkconstraints(n,k,kwargs)


class IntegerVectors_all(CombinatorialClass):
    def __repr__(self):
        """
        EXAMPLES:
            sage: IntegerVectors()
            Integer vectors
        """
        return "Integer vectors"

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: [] in IntegerVectors()
            True
            sage: [3,2,2,1] in IntegerVectors()
            True
        """
        if not isinstance(x, builtinlist):
            return False
        for i in x:
            if not isinstance(i, (int, Integer)):
                return False
            if i < 0:
                return False

        return True

    def list(self):
        """
        EXAMPLES:
            sage: IntegerVectors().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def count(self):
        """
        EXAMPLES:
            sage: IntegerVectors().count()
            +Infinity
        """
        return infinity

class IntegerVectors_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS:
            sage: IV = IntegerVectors(2,3)
            sage: IV == loads(dumps(IV))
            True

        AUTHORS:
            --Martin Albrecht
            --Mike Hansen
        """
        self.n = n
        self.k = k


    def _list_rec(self, n, k):
        """
        Return a list of a exponent tuples of length $size$ such that the
        degree of the associated monomial is $D$.

        INPUT:
            n -- degree (must be > 0)
            k -- length of exponent tuples (must be > 0)

        EXAMPLES:
            sage: IV = IntegerVectors(2,3)
            sage: IV._list_rec(2,3)
            [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
        """
        res = []

        if k == 1:
            return [ (n, ) ]

        for nbar in range(n+1):
            n_diff = n-nbar
            for rest in self._list_rec( nbar , k-1):
                res.append((n_diff,)+rest)
        return res

    def list(self):
        """
        EXAMPLE:
            sage: IV = IntegerVectors(2,3)
            sage: IV.list()
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
            sage: IntegerVectors(3, 0).list()
            []
            sage: IntegerVectors(3, 1).list()
            [[3]]
            sage: IntegerVectors(0, 1).list()
            [[0]]
            sage: IntegerVectors(0, 2).list()
            [[0, 0]]
            sage: IntegerVectors(2, 2).list()
            [[2, 0], [1, 1], [0, 2]]
        """
        if self.n < 0:
            return []

        if self.k == 0:
            if self.n == 0:
                return [[]]
            else:
                return []
        elif self.k == 1:
            return [[self.n]]

        res = self._list_rec(self.n, self.k)
        return map(list, res)


    def iterator(self):
        """
        EXAMPLE:
            sage: IV = IntegerVectors(2,3)
            sage: list(IV)
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
            sage: list(IntegerVectors(3, 0))
            []
            sage: list(IntegerVectors(3, 1))
            [[3]]
            sage: list(IntegerVectors(0, 1))
            [[0]]
            sage: list(IntegerVectors(0, 2))
            [[0, 0]]
            sage: list(IntegerVectors(2, 2))
            [[2, 0], [1, 1], [0, 2]]
            sage: IntegerVectors(0,0).list()
            [[]]
            sage: IntegerVectors(1,0).list()
            []
            sage: IntegerVectors(0,1).list()
            [[0]]
            sage: IntegerVectors(2,2).list()
            [[2, 0], [1, 1], [0, 2]]
            sage: IntegerVectors(-1,0).list()
            []
            sage: IntegerVectors(-1,2).list()
            []

        """
        if self.n < 0:
            return

        if self.k == 0:
            if self.n == 0:
                yield []
            return
        elif self.k == 1:
            yield [self.n]
            return

        for nbar in range(self.n+1):
            n = self.n-nbar
            for rest in IntegerVectors_nk(nbar , self.k-1):
                yield [n] + rest

    def __repr__(self):
        """
        TESTS:
            sage: IV = IntegerVectors(2,3)
            sage: repr(IV)
            'Integer vectors of length 3 that sum to 2'
        """
        return "Integer vectors of length %s that sum to %s"%(self.k, self.n)

    def __contains__(self, x):
        """
        TESTS:
            sage: IV = IntegerVectors(2,3)
            sage: all([i in IV for i in IV])
            True
            sage: [0,1,2] in IV
            False
            sage: [2.0, 0, 0] in IV
            False
            sage: [0,1,0,1] in IV
            False
            sage: [0,1,1] in IV
            True
            sage: [-1,2,1] in IV
            False
        """
        if x not in IntegerVectors():
            return False

        if sum(x) != self.n:
            return False

        if len(x) != self.k:
            return False

        if len(x) > 0 and min(x) < 0:
            return False

        return True

class IntegerVectors_nkconstraints(CombinatorialClass):
    def __init__(self, n, k, constraints):
        """
        EXAMPLES:
            sage: IV = IntegerVectors(2,3,min_slope=0)
            sage: IV == loads(dumps(IV))
            True
        """
        self.n = n
        self.k = k
        self.constraints = constraints

    def __repr__(self):
        """
        EXAMPLES:
            sage: IntegerVectors(2,3,min_slope=0).__repr__()
            'Integer vectors of length 3 that sum to 2 with constraints: min_slope=0'
        """
        return "Integer vectors of length %s that sum to %s with constraints: %s"%(self.k, self.n, ", ".join( ["%s=%s"%(key, self.constraints[key]) for key in sorted(self.constraints.keys())] ))


    def __contains__(self, x):
        """
        TESTS:
            sage: [0] in IntegerVectors(0)
            True
            sage: [0] in IntegerVectors(0, 1)
            True
            sage: [] in IntegerVectors(0, 0)
            True
            sage: [] in IntegerVectors(0, 1)
            False
            sage: [] in IntegerVectors(1, 0)
            False
            sage: [3] in IntegerVectors(3)
            True
            sage: [3] in IntegerVectors(2,1)
            False
            sage: [3] in IntegerVectors(2)
            False
            sage: [3] in IntegerVectors(3,1)
            True
            sage: [3,2,2,1] in IntegerVectors(9)
            False
            sage: [3,2,2,1] in IntegerVectors(9,5)
            False
            sage: [3,2,2,1] in IntegerVectors(8)
            True
            sage: [3,2,2,1] in IntegerVectors(8,5)
            False
            sage: [3,2,2,1] in IntegerVectors(8,4)
            True
            sage: [3,2,2,1] in IntegerVectors(8,4, min_part = 1)
            True
            sage: [3,2,2,1] in IntegerVectors(8,4, min_part = 2)
            False
        """
        if x not in IntegerVectors():
            return False

        if sum(x) != self.n:
            return False

        if len(x) != self.k:
            return False

        if self.constraints:
            if not misc.check_integer_list_constraints(x, singleton=True, **self.constraints):
                return False

        return True

    def count(self):
        """
        EXAMPLES:
            sage: IntegerVectors(3,3, min_part=1).count()
            1
            sage: IntegerVectors(5,3, min_part=1).count()
            6
            sage: IntegerVectors(13, 4, min_part=2, max_part=4).count()
            16
        """
        if not self.constraints:
            if self.n >= 0:
                return binomial(self.n+self.k-1,self.n)
            else:
                return 0
        else:
            if len(self.constraints) == 1 and 'max_part' in self.constraints and self.constraints['max_part'] != infinity:
                m = self.constraints['max_part']
                if m >= self.n:
                    return binomial(self.n+self.k-1,self.n)
                else: #do by inclusion / exclusion on the number
                      #i of parts greater than m
                    return sum( [(-1)**i * binomial(self.n+self.k-1-i*(m+1), self.k-1)*binomial(self.k,i) for i in range(0, self.n/(m+1)+1)])
            else:
                return len(self.list())


    def _parameters(self):
        """
        Returns a tuple (min_length, max_length, floor, ceiling, min_slope, max_slope)
        for the parameters of self.

        EXAMPLES:
            sage: IV = IntegerVectors(2,3,min_slope=0)
            sage: min_length, max_length, floor, ceiling, min_slope, max_slope = IV._parameters()
            sage: min_length
            3
            sage: max_length
            3
            sage: [floor(i) for i in range(1,10)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: [ceiling(i) for i in range(1,5)]
            [+Infinity, +Infinity, +Infinity, +Infinity]
            sage: min_slope
            0
            sage: max_slope
            +Infinity

        """
        constraints = self.constraints
        #n, min_length, max_length, floor, ceiling, min_slope, max_slope
        if self.k == -1:
            min_length = constraints.get('min_length', 0)
            max_length = constraints.get('max_length', infinity)
        else:
            min_length = self.k
            max_length = self.k

        min_part = constraints.get('min_part', 0)
        max_part = constraints.get('max_part', infinity)
        min_slope = constraints.get('min_slope', -infinity)
        max_slope = constraints.get('max_slope', infinity)
        if 'outer' in self.constraints:
            ceiling = list2func( map(lambda i: min(max_part, i), self.constraints['outer']), default=max_part )
        else:
            ceiling = constant_func(max_part)

        if 'inner' in self.constraints:
            floor = list2func( map(lambda i: max(min_part, i), self.constraints['outer']), default=min_part )
        else:
            floor = constant_func(min_part)

        return (min_length, max_length, floor, ceiling, min_slope, max_slope)


    def first(self):
        """
        EXAMPLES:
            sage: IntegerVectors(2,3,min_slope=0).first()
            [0, 1, 1]
        """
        return integer_list.first(self.n, *self._parameters())

    def next(self, x):
        """
        EXAMPLES:
            sage: IntegerVectors(2,3,min_slope=0).last()
            [0, 0, 2]
        """
        return integer_list.next(x, *self._parameters())

    def iterator(self):
        """
        EXAMPLES:
            sage: IntegerVectors(-1, 0, min_part = 1).list()
            []
            sage: IntegerVectors(-1, 2, min_part = 1).list()
            []
            sage: IntegerVectors(0, 0, min_part=1).list()
            [[]]
            sage: IntegerVectors(3, 0, min_part=1).list()
            []
            sage: IntegerVectors(0, 1, min_part=1).list()
            []
            sage: IntegerVectors(2, 2, min_part=1).list()
            [[1, 1]]
            sage: IntegerVectors(2, 3, min_part=1).list()
            []
            sage: IntegerVectors(4, 2, min_part=1).list()
            [[3, 1], [2, 2], [1, 3]]

            sage: IntegerVectors(0, 3, outer=[0,0,0]).list()
            [[0, 0, 0]]
            sage: IntegerVectors(1, 3, outer=[0,0,0]).list()
            []
            sage: IntegerVectors(2, 3, outer=[0,2,0]).list()
            [[0, 2, 0]]
            sage: IntegerVectors(2, 3, outer=[1,2,1]).list()
            [[1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1]]
            sage: IntegerVectors(2, 3, outer=[1,1,1]).list()
            [[1, 1, 0], [1, 0, 1], [0, 1, 1]]
            sage: IntegerVectors(2, 5, outer=[1,1,1,1,1]).list()
            [[1, 1, 0, 0, 0],
             [1, 0, 1, 0, 0],
             [1, 0, 0, 1, 0],
             [1, 0, 0, 0, 1],
             [0, 1, 1, 0, 0],
             [0, 1, 0, 1, 0],
             [0, 1, 0, 0, 1],
             [0, 0, 1, 1, 0],
             [0, 0, 1, 0, 1],
             [0, 0, 0, 1, 1]]

            sage: iv = [ IntegerVectors(n,k) for n in range(-2, 7) for k in range(7) ]
            sage: all(map(lambda x: x.count() == len(x.list()), iv))
            True
            sage: essai = [[1,1,1], [2,5,6], [6,5,2]]
            sage: iv = [ IntegerVectors(x[0], x[1], max_part = x[2]-1) for x in essai ]
            sage: all(map(lambda x: x.count() == len(x.list()), iv))
            True

        """
        return integer_list.iterator(self.n, *self._parameters())

class IntegerVectors_nconstraints(IntegerVectors_nkconstraints):
    def __init__(self, n, constraints):
        """
        TESTS:
            sage: IV = IntegerVectors(3, max_length=2)
            sage: IV == loads(dumps(IV))
            True
        """
        IntegerVectors_nkconstraints.__init__(self, n, -1, constraints)

    def __repr__(self):
        """
        EXAMPLES:
            sage: repr(IntegerVectors(3))
            'Integer vectors that sum to 3'
            sage: repr(IntegerVectors(3, max_length=2))
            'Integer vectors that sum to 3 with constraints: max_length=2'
        """
        if self.constraints:
            return "Integer vectors that sum to %s with constraints: %s"%(self.n,", ".join( ["%s=%s"%(key, self.constraints[key]) for key in sorted(self.constraints.keys())] ))
        else:
            return "Integer vectors that sum to %s"%(self.n,)

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: [0,3,0,1,2] in IntegerVectors(6)
            True
            sage: [0,3,0,1,2] in IntegerVectors(6, max_length=3)
            False
        """
        if self.constraints:
            return x in IntegerVectors_all() and misc.check_integer_list_constraints(x, singleton=True, **self.constraints)
        else:
            return x in IntegerVectors_all() and sum(x) == self.n

    def count(self):
        """
        EXAMPLES:
            sage: IntegerVectors(3).count()
            +Infinity
        """
        return infinity

    def list(self):
        """
        EXAMPLES:
            sage: IntegerVectors(3, max_length=2).list()
            [[3], [2, 1], [1, 2], [0, 3]]
            sage: IntegerVectors(3).list()
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list
        """
        if 'max_length' not in self.constraints:
            raise NotImplementedError, "infinite list"
        else:
            return list(self.iterator())


class IntegerVectors_nnondescents(CombinatorialClass):
    r"""
    The combinatorial class of integer vectors v graded by two parameters:
     - n: the sum of the parts of v
     - comp: the non descents composition of v

    In other words: the length of v equals c[1]+...+c[k], and v is
    descreasing in the consecutive blocs of length c[1], ..., c[k]

    Those are the integer vectors of sum n which are lexicographically
    maximal (for the natural left->right reading) in their orbit by the
    young subgroup S_{c_1} x \dots x S_{c_k}.  In particular, they form
    a set of orbit representative of integer vectors w.r.t. this young
    subgroup.
    """
    def __init__(self, n, comp):
        """
        EXAMPLES:
            sage: IV = IntegerVectors(4, [2])
            sage: IV == loads(dumps(IV))
            True
        """
        self.n = n
        self.comp = comp

    def __repr__(self):
        """
        EXAMPLES:
            sage: IntegerVectors(4, [2]).__repr__()
            'Integer vectors of 4 with non-descents composition [2]'
        """
        return "Integer vectors of %s with non-descents composition %s"%(self.n, self.comp)

    def iterator(self):
        """
        TESTS:
            sage: IntegerVectors(0, []).list()
            [[]]
            sage: IntegerVectors(5, []).list()
            []
            sage: IntegerVectors(0, [1]).list()
            [[0]]
            sage: IntegerVectors(4, [1]).list()
            [[4]]
            sage: IntegerVectors(4, [2]).list()
            [[4, 0], [3, 1], [2, 2]]
            sage: IntegerVectors(4, [2,2]).list()
             [[4, 0, 0, 0],
             [3, 0, 1, 0],
             [2, 0, 2, 0],
             [2, 0, 1, 1],
             [1, 0, 3, 0],
             [1, 0, 2, 1],
             [0, 0, 4, 0],
             [0, 0, 3, 1],
             [0, 0, 2, 2]]
            sage: IntegerVectors(5, [1,1,1]).list()
            [[5, 0, 0],
             [4, 1, 0],
             [4, 0, 1],
             [3, 2, 0],
             [3, 1, 1],
             [3, 0, 2],
             [2, 3, 0],
             [2, 2, 1],
             [2, 1, 2],
             [2, 0, 3],
             [1, 4, 0],
             [1, 3, 1],
             [1, 2, 2],
             [1, 1, 3],
             [1, 0, 4],
             [0, 5, 0],
             [0, 4, 1],
             [0, 3, 2],
             [0, 2, 3],
             [0, 1, 4],
             [0, 0, 5]]
            sage: IntegerVectors(0, [2,3]).list()
            [[0, 0, 0, 0, 0]]

         """
        for iv in IntegerVectors(self.n, len(self.comp)):
            blocks = [ IntegerVectors(iv[i], self.comp[i], max_slope=0).iterator() for i in range(len(self.comp))]
            for parts in cartesian_product.CartesianProduct(*blocks):
                res = []
                for part in parts:
                    res += part
                yield res


