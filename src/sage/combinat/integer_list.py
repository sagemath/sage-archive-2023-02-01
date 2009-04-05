r"""
Tools for generating lists of integers in lexicographic order.

IMPORTANT NOTE (2009/02):
The internal functions in this file will be deprecated soon.
Please only use them through IntegerListsLex.

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

import generator
from sage.calculus.calculus import floor as flr
from sage.rings.arith import binomial
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.misc.lazy_attribute import lazy_attribute
import __builtin__

def first(n, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Returns the lexicographically smallest valid composition of n
    satisfying the conditions.

    .. warning::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    Preconditions:

    - minslope < maxslope

    - floor and ceiling need to satisfy the slope constraints,
      e.g. be obtained fromcomp2floor or comp2ceil

    - floor must be below ceiling to ensure
      the existence a valid composition

    TESTS::

        sage: import sage.combinat.integer_list as integer_list
        sage: f = lambda l: lambda i: l[i-1]
        sage: f([0,1,2,3,4,5])(1)
        0
        sage: integer_list.first(12, 4, 4, f([0,0,0,0]), f([4,4,4,4]), -1, 1)
        [4, 3, 3, 2]
        sage: integer_list.first(36, 9, 9, f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -1, 1)
        [7, 6, 5, 5, 4, 3, 3, 2, 1]
        sage: integer_list.first(25, 9, 9, f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        [7, 6, 5, 4, 2, 1, 0, 0, 0]
        sage: integer_list.first(36, 9, 9, f([3,3,3,2,1,4,2,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        [7, 6, 5, 5, 5, 4, 3, 1, 0]
    """

    #Increase minl until n <= sum([ceiling(i) for i in range(min_length)])
    #This may run forever!
    N = sum([ceiling(i) for i in range(1,min_length+1)])
    while N < n:
        min_length += 1
        if min_length > max_length:
            return None

        if ceiling(min_length) == 0 and max_slope <= 0:
            return None

        N += ceiling(min_length)


    #This is the place where it's required that floor(i)
    #respects the floor conditions.
    n -= sum([floor(i) for i in range(1, min_length+1)])
    if n < 0:
        return None

    #Now we know that we can build the composition inside
    #the "tube" [1 ... min_length] * [floor, ceiling]

    if min_slope == -infinity:
        #Easy case: min_slope == -infinity
        result = []
        i = min_length
        for i in range(1,min_length+1):
            if n <= ceiling(i) - floor(i):
                result.append(floor(i) + n)
                break
            else:
                result.append(ceiling(i))
                n -= ceiling(i) - floor(i)

        result += [floor(j) for j in range(i+1,min_length+1)]
        return result

    else:
        if n == 0 and min_length == 0:
            return []

        low_x = 1
        low_y = floor(1)
        high_x = 1
        high_y = floor(1)

        while n > 0:
            #invariant after each iteration of the loop:
            #[low_x, low_y] is the coordinate of the rightmost point of the
            #current diagonal s.t. floor(low_x) < low_y
            low_y += 1
            while low_x < min_length and low_y+min_slope > floor(low_x + 1):
                low_x += 1
                low_y += min_slope

            high_y += 1
            while high_y > ceiling(high_x):
                high_x += 1
                high_y += min_slope

            n -= low_x - high_x + 1

        #print "lx, ly, hw, hy, n", low_x, low_y, high_x, high_y, n
        #print (high_x-1 - 1 + 1)  + (n + 1 - 0 + 1) + ( low_x-high_x - n + 1) + (min_length - (low_x + 1) +1)
        result = []
        result += [ ceiling(j) for j in range(1,high_x)]
        result += [ high_y + min_slope*i - 1 for i in range(0, -n) ]
        result += [ high_y + min_slope*i for i in range(-n, low_x-high_x+1)]
        result += [ floor(j) for j in range(low_x+1,min_length+1) ]

        return result


def lower_regular(comp, min_slope, max_slope):
    """
    Returns the uppest regular composition below comp

    TESTS::

        sage: import sage.combinat.integer_list as integer_list
        sage: integer_list.lower_regular([4,2,6], -1, 1)
        [3, 2, 3]
        sage: integer_list.lower_regular([4,2,6], -1, infinity)
        [3, 2, 6]
        sage: integer_list.lower_regular([1,4,2], -1, 1)
        [1, 2, 2]
        sage: integer_list.lower_regular([4,2,6,3,7], -2, 1)
        [4, 2, 3, 3, 4]
        sage: integer_list.lower_regular([4,2,infinity,3,7], -2, 1)
        [4, 2, 3, 3, 4]
        sage: integer_list.lower_regular([1, infinity, 2], -1, 1)
        [1, 2, 2]
        sage: integer_list.lower_regular([infinity, 4, 2], -1, 1)
        [4, 3, 2]
    """

    new_comp = comp[:]
    for i in range(1, len(new_comp)):
        new_comp[i] = min(new_comp[i], new_comp[i-1] + max_slope)

    for i in reversed(range(len(new_comp)-1)):
        new_comp[i] = min( new_comp[i], new_comp[i+1] - min_slope)

    return new_comp

def rightmost_pivot(comp, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    TESTS::

        sage: import sage.combinat.integer_list as integer_list
        sage: f = lambda l: lambda i: l[i-1]
        sage: integer_list.rightmost_pivot([7,6,5,5,4,3,3,2,1], 9, 9, f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -1, 0)
        [7, 2]
        sage: integer_list.rightmost_pivot([7,6,5,5,4,3,3,2,1], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 0)
        [7, 1]
        sage: integer_list.rightmost_pivot([7,6,5,5,4,3,3,2,1], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 4)
        [8, 1]
        sage: integer_list.rightmost_pivot([7,6,5,5,4,3,3,2,1], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        [8, 1]
        sage: integer_list.rightmost_pivot([7,6,5,5,5,5,5,4,4], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        sage: integer_list.rightmost_pivot([3,3,3,2,1,1,0,0,0], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        sage: g = lambda x: lambda i: x
        sage: integer_list.rightmost_pivot([1],1,1,g(0),g(2),-10, 10)
        sage: integer_list.rightmost_pivot([1,2],2,2,g(0),g(2),-10, 10)
        sage: integer_list.rightmost_pivot([1,2],2,2,g(1),g(2), -10, 10)
        sage: integer_list.rightmost_pivot([1,2],2,3,g(1),g(2), -10, 10)
        [2, 1]
        sage: integer_list.rightmost_pivot([2,2],2,3,g(2),g(2),-10, 10)
        sage: integer_list.rightmost_pivot([2,3],2,3,g(2),g(2),-10,+10)
        sage: integer_list.rightmost_pivot([3,2],2,3,g(2),g(2),-10,+10)
        sage: integer_list.rightmost_pivot([3,3],2,3,g(2),g(2),-10,+10)
        [1, 2]
        sage: integer_list.rightmost_pivot([6],1,3,g(0),g(6),-1,0)
        [1, 0]
        sage: integer_list.rightmost_pivot([6],1,3,g(0),g(6),-2,0)
        [1, 0]
        sage: integer_list.rightmost_pivot([7,9,8,7],1,5,g(0),g(10),-1,10)
        [2, 6]
        sage: integer_list.rightmost_pivot([7,9,8,7],1,5,g(5),g(10),-10,10)
        [3, 5]
        sage: integer_list.rightmost_pivot([7,9,8,7],1,5,g(5),g(10),-1,10)
        [2, 6]
        sage: integer_list.rightmost_pivot([7,9,8,7],1,5,g(4),g(10),-2,10)
        [3, 7]
        sage: integer_list.rightmost_pivot([9,8,7],1,4,g(4),g(10),-2,0)
        [1, 4]
        sage: integer_list.rightmost_pivot([1,3],1,5,lambda i: i,g(10),-10,10)
        sage: integer_list.rightmost_pivot([1,4],1,5,lambda i: i,g(10),-10,10)
        sage: integer_list.rightmost_pivot([2,4],1,5,lambda i: i,g(10),-10,10)
        [1, 1]
    """

    y = len(comp) + 1
    while y <= max_length:
        if ceiling(y) > 0:
            break
        if max_slope <= 0:
            y = max_length + 1
            break
        y += 1

    x = len(comp)
    if x == 0:
        return None

    ceilingx_x = comp[x-1]-1
    floorx_x = floor(x)
    if x > 1:
        floorx_x = max(floorx_x, comp[x-2]+min_slope)

    F = comp[x-1] - floorx_x
    G = ceilingx_x - comp[x-1] #this is -1

    highX = x
    lowX  = x

    while not (ceilingx_x >= floorx_x and
               (G >= 0 or
               ( y < max_length +1 and
                 F - max(floor(y), floorx_x + (y-x)*min_slope) >= 0 and
                 G + min(ceiling(y), ceilingx_x + (y-x)*max_slope) >= 0 ))):

        if x == 1:
            return None

        x -= 1

        oldfloorx_x    = floorx_x
        ceilingx_x = comp[x-1] - 1
        floorx_x   = floor(x)
        if x > 1:
            floorx_x = max(floorx_x, comp[x-2]+min_slope)

        min_slope_lowX  = min_slope*(lowX - x)
        max_slope_highX = max_slope*(highX - x)


        #Update G
        if max_slope == infinity:
            #In this case, we have
            #  -- ceiling_x(i) = ceiling(i) for i > x
            #  --G >= 0 or G = -1
            G += ceiling(x+1)-comp[x]
        else:
            G += (highX - x)*( (comp[x-1]+max_slope) - comp[x]) - 1
            temp = (ceilingx_x + max_slope_highX) - ceiling(highX)
            while highX > x and ( temp >= 0 ):
                G  -= temp
                highX -= 1
                max_slope_highX = max_slope*(highX-x)
                temp = (ceilingx_x + max_slope_highX) - ceiling(highX)

        if G >= 0 and comp[x-1] > floorx_x:
            #By case 1, x is at the rightmost pivot position
            break

        #Update F
        if y < max_length+1:
            F += comp[x-1] - floorx_x
            if min_slope != -infinity:
                F += (lowX - x) * (oldfloorx_x - (floorx_x + min_slope))
                temp = floor(lowX) - (floorx_x + min_slope_lowX)
                while lowX > x and temp >= 0:
                    F -= temp
                    lowX -= 1
                    min_slope_lowX = min_slope*(lowX-x)
                    temp = floor(lowX) - (floorx_x + min_slope_lowX)

    return [x, floorx_x]

def next(comp, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Returns the next integer list after comp that satisfies the
    constraints.

    .. warning::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    EXAMPLES::

        sage: from sage.combinat.integer_list import next
        sage: IV = IntegerVectors(2,3,min_slope=0)
        sage: params = IV._parameters()
        sage: next([0,1,1], *params)
        [0, 0, 2]
    """
    x = rightmost_pivot( comp, min_length, max_length, floor, ceiling, min_slope, max_slope)
    if x == None:
        return None
    [x, low] = x
    high = comp[x-1]-1

##     // Build wrappers around floor and ceiling to take into
##     // account the new constraints on the value of compo[x].
##     //
##     // Efficiency note: they are not wrapped more than once, since
##     // the method Next calls first, but not the converse.

    if min_slope == -infinity:
        new_floor = lambda i: floor(x+(i-1))
    else:
        new_floor = lambda i: max(floor(x+(i-1)), low+(i-1)*min_slope)

    if max_slope == infinity:
        new_ceiling = lambda i: comp[x-1] - 1 if i == 1 else ceiling(x+(i-1))
    else:
        new_ceiling = lambda i: min(ceiling(x+(i-1)), high+(i-1)*max_slope)


    res = []
    res += comp[:x-1]
    res += first(sum(comp[x-1:]), max(min_length-x+1, 0), max_length-x+1,
                 new_floor, new_ceiling, min_slope, max_slope)
    return res



def iterator(n, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    .. warning::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    EXAMPLES::

        sage: from sage.combinat.integer_list import iterator
        sage: IV = IntegerVectors(2,3,min_slope=0)
        sage: params = IV._parameters()
        sage: list(iterator(2,*params))
        [[0, 1, 1], [0, 0, 2]]
    """
    #from sage.misc.misc import deprecation
    #deprecation("sage.combinat.integer_list.iterator is deprecated. Please use IntegerListsLex(...)")
    succ = lambda x: next(x, min_length, max_length, floor, ceiling, min_slope, max_slope)

    #Handle the case where n is a list of integers
    if isinstance(n, __builtin__.list):
        iterators = [iterator(i, min_length, max_length, floor, ceiling, min_slope, max_slope) for i in range(n[0], min(n[1]+1,upper_bound(min_length, max_length, floor, ceiling, min_slope, max_slope)))]

        return generator.concat(iterators)
    else:
        f =  first(n, min_length, max_length, floor, ceiling, min_slope, max_slope)
        if f == None:
            return generator.element(None, 0)
        return generator.successor(f, succ)

def list(n, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    .. warning::

    THIS FUNCTION IS DEPRECATED!

    Please use IntegersListsLex(...) instead

    EXAMPLES::

        sage: import sage.combinat.integer_list as integer_list
        sage: g = lambda x: lambda i: x
        sage: integer_list.list(0,0,infinity,g(1),g(infinity),0,infinity)
        [[]]
        sage: integer_list.list(0,0,infinity,g(1),g(infinity),0,0)
        [[]]
        sage: integer_list.list(0, 0, 0, g(1), g(infinity), 0, 0)
        [[]]
        sage: integer_list.list(0, 0, 0, g(0), g(infinity), 0, 0)
        [[]]
        sage: integer_list.list(0, 0, infinity, g(1), g(infinity), 0, infinity)
        [[]]
        sage: integer_list.list(1, 0, infinity, g(1), g(infinity), 0, infinity)
        [[1]]
        sage: integer_list.list(0, 1, infinity, g(1), g(infinity), 0, infinity)
        []
        sage: integer_list.list(0, 1, infinity, g(0), g(infinity), 0, infinity)
        [[0]]
        sage: integer_list.list(3, 0, 2, g(0), g(infinity), -infinity, infinity)
        [[3], [2, 1], [1, 2], [0, 3]]
        sage: partitions = (0, infinity, g(0), g(infinity), -infinity, 0)
        sage: partitions_min_2 = (0, infinity, g(2), g(infinity), -infinity, 0)
        sage: compositions = (0, infinity, g(1), g(infinity), -infinity, infinity)
        sage: integer_vectors = lambda l: (l, l, g(0), g(infinity), -infinity, infinity)
        sage: lower_monomials = lambda c: (len(c), len(c), g(0), lambda i: c[i], -infinity, infinity)
        sage: upper_monomials = lambda c: (len(c), len(c), g(0), lambda i: c[i], -infinity, infinity)
        sage: constraints = (0, infinity, g(1), g(infinity), -1, 0)
        sage: integer_list.list(6, *partitions)
        [[6],
         [5, 1],
         [4, 2],
         [4, 1, 1],
         [3, 3],
         [3, 2, 1],
         [3, 1, 1, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: integer_list.list(6, *constraints)
        [[6],
         [3, 3],
         [3, 2, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: integer_list.list(1, *partitions_min_2)
        []
        sage: integer_list.list(2, *partitions_min_2)
        [[2]]
        sage: integer_list.list(3, *partitions_min_2)
        [[3]]
        sage: integer_list.list(4, *partitions_min_2)
        [[4], [2, 2]]
        sage: integer_list.list(5, *partitions_min_2)
        [[5], [3, 2]]
        sage: integer_list.list(6, *partitions_min_2)
        [[6], [4, 2], [3, 3], [2, 2, 2]]
        sage: integer_list.list(7, *partitions_min_2)
        [[7], [5, 2], [4, 3], [3, 2, 2]]
        sage: integer_list.list(9, *partitions_min_2)
        [[9], [7, 2], [6, 3], [5, 4], [5, 2, 2], [4, 3, 2], [3, 3, 3], [3, 2, 2, 2]]
        sage: integer_list.list(10, *partitions_min_2)
        [[10],
         [8, 2],
         [7, 3],
         [6, 4],
         [6, 2, 2],
         [5, 5],
         [5, 3, 2],
         [4, 4, 2],
         [4, 3, 3],
         [4, 2, 2, 2],
         [3, 3, 2, 2],
         [2, 2, 2, 2, 2]]
        sage: integer_list.list(4, *compositions)
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
    """
    #deprecation("sage.combinat.integer_list.list(...) is deprecated. Please use list(IntegerListsLex(...))")
    return __builtin__.list(iterator(n, min_length, max_length, floor, ceiling, min_slope, max_slope))

def upper_regular(comp, min_slope, max_slope):
    """
    Returns the uppest regular composition above comp.

    TESTS::

        sage: import sage.combinat.integer_list as integer_list
        sage: integer_list.upper_regular([4,2,6],-1,1)
        [4, 5, 6]
        sage: integer_list.upper_regular([4,2,6],-2, 1)
        [4, 5, 6]
        sage: integer_list.upper_regular([4,2,6,3,7],-2, 1)
        [4, 5, 6, 6, 7]
        sage: integer_list.upper_regular([4,2,6,1], -2, 1)
        [4, 5, 6, 4]
    """

    new_comp = comp[:]
    for i in range(1, len(new_comp)):
        new_comp[i] = max(new_comp[i], new_comp[i-1] + min_slope)

    for i in reversed(range(len(new_comp)-1)):
        new_comp[i] = max( new_comp[i], new_comp[i+1] - max_slope)

    return new_comp

def comp2floor(f, min_slope, max_slope):
    """
    Given a composition, returns the lowest regular function N->N above
    it.

    EXAMPLES::

        sage: from sage.combinat.integer_list import comp2floor
        sage: f = comp2floor([2, 1, 1],-1,0)
        sage: [f(i) for i in range(10)]
        [2, 1, 1, 1, 2, 3, 4, 5, 6, 7]
    """
    if len(f) == 0: return lambda i: 0
    floor = upper_regular(f, min_slope, max_slope)
    return lambda i: floor[i] if i < len(floor) else max(0, floor[-1]-(i-len(floor))*min_slope)


def comp2ceil(c, min_slope, max_slope):
    """
    Given a composition, returns the lowest regular function N->N below
    it.

    EXAMPLES::

        sage: from sage.combinat.integer_list import comp2ceil
        sage: f = comp2ceil([2, 1, 1],-1,0)
        sage: [f(i) for i in range(10)]
        [2, 1, 1, 1, 2, 3, 4, 5, 6, 7]
    """
    if len(c) == 0: return lambda i: 0
    ceil = lower_regular(c, min_slope, max_slope)
    return lambda i: ceil[i] if i < len(ceil) else max(0, ceil[-1]-(i-len(ceil))*min_slope)


def upper_bound(min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Compute a coarse upper bound on the size of a vector satisfying the
    constraints.

    TESTS::

        sage: import sage.combinat.integer_list as integer_list
        sage: f = lambda x: lambda i: x
        sage: integer_list.upper_bound(0,4,f(0), f(1),-infinity,infinity)
        4
        sage: integer_list.upper_bound(0, infinity, f(0), f(1), -infinity, infinity)
        +Infinity
        sage: integer_list.upper_bound(0, infinity, f(0), f(1), -infinity, -1)
        1
        sage: integer_list.upper_bound(0, infinity, f(0), f(5), -infinity, -1)
        15
        sage: integer_list.upper_bound(0, infinity, f(0), f(5), -infinity, -2)
        9
    """

    if max_length < infinity:
        return sum( [ ceiling(j) for j in range(max_length)] )
    elif max_slope < 0 and ceiling(1) < infinity:
        maxl = flr(-ceiling(1)/max_slope)
        return ceiling(1)*(maxl+1) + binomial(maxl+1,2)*max_slope
    #FIXME: only checking the first 10000 values, but that should generally
    #be enough
    elif [ceiling(j) for j in range(10000)] == [0]*10000:
        return 0
    else:
        return infinity



def is_a(comp, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Returns True if comp meets the constraints imposed by the
    arguments.

    .. warning::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    EXAMPLES::

        sage: from sage.combinat.integer_list import is_a
        sage: IV = IntegerVectors(2,3,min_slope=0)
        sage: params = IV._parameters()
        sage: all([is_a(iv, *params) for iv in IV])
        True
    """
    if len(comp) < min_length or len(comp) > max_length:
        return False
    for i in range(len(comp)):
        if comp[i] < floor(i+1):
            return False
        if comp[i] > ceiling(i+1):
            return False
    for i in range(len(comp)-1):
        slope = comp[i+1] - comp[i]
        if slope < min_slope or slope > max_slope:
            return False
    return True

from combinat import CombinatorialClass

class IntegerListsLex(CombinatorialClass):
    r"""
    A combinatorial class `C` for integer lists satisfying certain
    sum, length, upper/lower bound and regularity constraints. The
    purpose of this tool is mostly to provide a Constant Amortized
    Time iterator through those lists, in lexicographic order.

    INPUT:

    - ``n`` -  a non negative integer
    - ``min_length`` -  a non negative integer
    - ``max_length`` -  an integer or `\infty`
    - ``length`` -  an integer; overrides min_length and max_length if specified
    - ``floor`` -  a function `f` (or list);    defaults to ``lambda i: 0``
    - ``ceiling`` -  a function `f` (or list);  defaults to ``lambda i: infinity``
    - ``min_slope`` -  an integer or `-\infty`; defaults to `-\infty`
    - ``max_slope`` -  an integer or `+\infty`; defaults to `+\infty`

    An *integer list* is a list `l` of nonnegative integers, its
    *parts*. The *length* of `l` is the number of its parts;
    the *sum* of `l` is the sum of its parts.

    .. note::

       Two valid integer lists are considered equivalent if they only
       differ by trailing zeroes. In this case, only the list with the
       least number of trailing zeroes will be produced.

    The constraints on the lists are as follow:

    - Sum: `sum(l) == n`

    - Length: `min\_length \leq len(l) \leq max\_length`

    - Lower and upper bounds: `floor(i) \leq l[i] \leq ceiling(i)`, for `i=0,\dots,len(l)`

    - Regularity condition: `minSlope \leq l[i+1]-l[i] \leq maxSlope`, for `i=0,\dots,len(l)-1`

    This is a generic low level tool. The interface has been designed
    with efficiency in mind. It is subject to incompatible changes in
    the future. More user friendly interfaces are provided by high
    level tools like Partitions or Compositions.

    EXAMPLES:

    We create the combinatorial class of lists of length 3 and sum 2::

        sage: C = IntegerListsLex(2, length=3)
        sage: C
        Integer lists of sum 2 satisfying certain constraints
        sage: C.cardinality()
        6
        sage: [p for p in C]
        [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

        sage: [2, 0, 0] in C
        True
        sage: [2, 0, 1] in C
        False
        sage: "a" in C
        False
        sage: ["a"] in C
        False

        sage: C.first()
        [2, 0, 0]

    One can specify lower and upper bound on each part::

        sage: list(IntegerListsLex(5, length = 3, floor = [1,2,0], ceiling = [3,2,3]))
        [[3, 2, 0], [2, 2, 1], [1, 2, 2]]

    Using the slope condition, one can generate integer partitions
    (but see :mod:`sage.combinat.partition.Partitions`)::

        sage: list(IntegerListsLex(4, max_slope=0))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

    This is the list of all partitions of `7` with parts at least `2`::

        sage: list(IntegerListsLex(7, max_slope = 0, min_part = 2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]

    This is the list of all partitions of `5` and length at most 3
    which are bounded below by [2,1,1]::

        sage: list(IntegerListsLex(5, max_slope = 0, max_length = 3, floor = [2,1,1]))
        [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1]]

    Note that [5] is considered valid, because the lower bound
    constraint only apply to existing positions in the list. To
    obtain instead the partitions containing [2,1,1], one need to
    use min_length::

        sage: list(IntegerListsLex(5, max_slope = 0, min_length = 3, max_length = 3, floor = [2,1,1]))
        [[3, 1, 1], [2, 2, 1]]

    This is the list of all partitions of `5` which are contained in [3,2,2]::

        sage: list(IntegerListsLex(5, max_slope = 0, max_length = 3, ceiling = [3,2,2]))
        [[3, 2], [3, 1, 1], [2, 2, 1]]


    This is the list of all compositions of `4` (but see Compositions)::

        sage: list(IntegerListsLex(4, min_part = 1))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]

    This is the list of all integer vectors of sum `4` and length `3`::

        sage: list(IntegerListsLex(4, length = 3))
        [[4, 0, 0], [3, 1, 0], [3, 0, 1], [2, 2, 0], [2, 1, 1], [2, 0, 2], [1, 3, 0], [1, 2, 1], [1, 1, 2], [1, 0, 3], [0, 4, 0], [0, 3, 1], [0, 2, 2], [0, 1, 3], [0, 0, 4]]


    There are all the lists of sum 4 and length 4 such that l[i] <= i::

        sage: list(IntegerListsLex(4, length=4, ceiling=lambda i: i))
        [[0, 1, 2, 1], [0, 1, 1, 2], [0, 1, 0, 3], [0, 0, 2, 2], [0, 0, 1, 3]]

    This is the list of all monomials of degree `4` which divide the
    monomial `x^3y^1z^2` (a monomial being identified with its
    exponent vector)::

        sage: R.<x,y,z> = QQ[]
        sage: m = [3,1,2]
        sage: def term(exponents):
        ...       return x^exponents[0] * y^exponents[1] * z^exponents[2]
        ...
        sage: list( IntegerListsLex(4, length = len(m), ceiling = m, element_constructor = term) )
        [x^3*y, x^3*z, x^2*y*z, x^2*z^2, x*y*z^2]

    Note the use of the element_constructor feature.

    In general, the complexity of the iteration algorithm is constant
    time amortized for each integer list produced.  There is one
    degenerate case though where the algorithm may run forever without
    producing anything. If max_length is `+\infty` and max_slope `>0`,
    testing whether there exists a valid integer list of sum `n` may
    be only semi-decidable. In the following example, the algorithm
    will enter an infinite loop, because it needs to decide whether
    `ceiling(i)` is nonzero for some `i`::

        sage: list( IntegerListsLex(1, ceiling = lambda i: 0) ) # todo: not implemented


    .. note::

       Caveat: counting is done by brute force generation. In some
       special cases, it would be possible to do better by counting
       techniques for integral point in polytopes.

    .. note::

       Caveat: with the current implementation, the constraints should
       satisfy the following conditions:

       - The upper and lower bounds themselves should satisfy the
         slope constraints.

       - The maximal and minimal slopes values should not be equal.

       - The maximal and minimal part values should not be equal.

    Those conditions are not checked by the algorithm, and the
    result may be completely incorrect if they are not satisfied:

    In the following example, the slope condition is not satisfied
    by the upper bound on the parts, and [3,3] is erroneously
    included in the result::

        sage: list(IntegerListsLex(6, max_part=3,max_slope=-1))
        [[3, 3], [3, 2, 1]]


    With some work, this could be fixed withoug affecting the overall
    complexity and efficiency. Also, the generation algorithm could be
    extended to deal with non-constant slope constraints and with
    negative parts, as well as to accept a range parameter instead of
    a single integer for the sum `n` of the lists (the later was
    readilly implemented in MuPAD-Combinat). Encouragements,
    suggestions, and help are welcome.

    TODO: integrate all remaining tests from http://mupad-combinat.svn.sourceforge.net/viewvc/mupad-combinat/trunk/MuPAD-Combinat/lib/COMBINAT/TEST/MachineIntegerListsLex.tst

    TESTS:
        sage: g = lambda x: lambda i: x
        sage: list(IntegerListsLex(0, floor = g(1), min_slope = 0))
        [[]]
        sage: list(IntegerListsLex(0, floor = g(1), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(IntegerListsLex(0, max_length=0, floor = g(1), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(IntegerListsLex(0, max_length=0, floor = g(0), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(IntegerListsLex(0, min_part = 1, min_slope = 0))
        [[]]
        sage: list(IntegerListsLex(1, min_part = 1, min_slope = 0))
        [[1]]
        sage: list(IntegerListsLex(0, min_length = 1, min_part = 1, min_slope = 0))
        []
        sage: list(IntegerListsLex(0, min_length = 1, min_slope = 0))
        [[0]]
        sage: list(IntegerListsLex(3, max_length=2, ))
        [[3], [2, 1], [1, 2], [0, 3]]
        sage: partitions = {"min_part": 1, "max_slope": 0}
        sage: partitions_min_2 = {"floor": g(2), "max_slope": 0}
        sage: compositions = {"min_part": 1}
        sage: integer_vectors = lambda l: {"length": l}
        sage: lower_monomials = lambda c: {"length": c, "floor": lambda i: c[i]}
        sage: upper_monomials = lambda c: {"length": c, "ceiling": lambda i: c[i]}
        sage: constraints = { "min_part":1, "min_slope": -1, "max_slope": 0}
        sage: list(IntegerListsLex(6, **partitions))
        [[6],
         [5, 1],
         [4, 2],
         [4, 1, 1],
         [3, 3],
         [3, 2, 1],
         [3, 1, 1, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(IntegerListsLex(6, **constraints))
        [[6],
         [3, 3],
         [3, 2, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(IntegerListsLex(1, **partitions_min_2))
        []
        sage: list(IntegerListsLex(2, **partitions_min_2))
        [[2]]
        sage: list(IntegerListsLex(3, **partitions_min_2))
        [[3]]
        sage: list(IntegerListsLex(4, **partitions_min_2))
        [[4], [2, 2]]
        sage: list(IntegerListsLex(5, **partitions_min_2))
        [[5], [3, 2]]
        sage: list(IntegerListsLex(6, **partitions_min_2))
        [[6], [4, 2], [3, 3], [2, 2, 2]]
        sage: list(IntegerListsLex(7, **partitions_min_2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]
        sage: list(IntegerListsLex(9, **partitions_min_2))
        [[9], [7, 2], [6, 3], [5, 4], [5, 2, 2], [4, 3, 2], [3, 3, 3], [3, 2, 2, 2]]
        sage: list(IntegerListsLex(10, **partitions_min_2))
        [[10],
         [8, 2],
         [7, 3],
         [6, 4],
         [6, 2, 2],
         [5, 5],
         [5, 3, 2],
         [4, 4, 2],
         [4, 3, 3],
         [4, 2, 2, 2],
         [3, 3, 2, 2],
         [2, 2, 2, 2, 2]]
        sage: list(IntegerListsLex(4, **compositions))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]

    """
    def __init__(self,
                 n,
                 length = None, min_length=0, max_length=infinity,
                 floor=None, ceiling = None,
                 min_part = 0, max_part = infinity,
                 min_slope=-infinity, max_slope=infinity,
                 name = None,
                 element_constructor = None):
        """
        TESTS::

            sage: C = IntegerListsLex(2, length=3)
            sage: C == loads(dumps(C))
            True
            sage: C == loads(dumps(C)) # this did fail at some point, really!
            True
            sage: C is loads(dumps(C)) # todo: not implemented
            True
            sage: C.cardinality().parent() is ZZ
            True
        """
        if floor is None:
            self.floor_list = []
        elif type(floor) is type([]): # FIXME: how to refer to type list rather than the function list above?
            self.floor_list = floor
        else:
            self.floor = floor
        if ceiling is None:
            self.ceiling_list = []
        elif type(ceiling) is type([]):
            self.ceiling_list = ceiling
        else:
            self.ceiling = ceiling
        if length is not None:
            min_length = length
            max_length = length
        if name is not None:
            self._name = name
        if n in ZZ:
            self.n = n
            self.n_range = [n]
        else:
            self.n_range = n
        self.min_length = min_length
        self.max_length = max_length
        self.min_part = min_part
        self.max_part = max_part
        # FIXME: the internal functions currently assume that floor and ceiling start at 1
        # this is a workaround
        self.max_slope = max_slope
        self.min_slope = min_slope
        if element_constructor is not None:
            self._element_constructor_ = element_constructor

    _element_constructor_ = type([])

    @lazy_attribute
    def _name(self):
        """
        Returns the name of this combinatorial class.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: C._name
            'Integer lists of sum 2 satisfying certain constraints'

            sage: C = IntegerListsLex([1,2,4], length=3)
            sage: C._name
            'Integer lists of sum in [1, 2, 4] satisfying certain constraints'
        """
        if hasattr(self, "n"):
            return "Integer lists of sum %s satisfying certain constraints"%self.n
        else:
            return "Integer lists of sum in %s satisfying certain constraints"%self.n_range

    def floor(self, i):
        """
        Returns the minimum part that can appear in the $i^{th}$ position in any
        list produced.

        EXAMPLES::

            sage: C = IntegerListsLex(4, length=2, min_part=1)
            sage: C.floor(0)
            1
            sage: C = IntegerListsLex(4, length=2, floor=[1,2])
            sage: C.floor(0)
            1
            sage: C.floor(1)
            2

        """
        return self.floor_list[i]   if i < len(self.floor_list  ) else self.min_part

    def ceiling(self, i):
        """
        Returns the maximum part that can appear in the $i^{th}
        position in any list produced.

        EXAMPLES::

            sage: C = IntegerListsLex(4, length=2, max_part=3)
            sage: C.ceiling(0)
            3
            sage: C = IntegerListsLex(4, length=2, ceiling=[3,2])
            sage: C.ceiling(0)
            3
            sage: C.ceiling(1)
            2

        """
        return self.ceiling_list[i] if i < len(self.ceiling_list) else self.max_part


    # Temporary adapter to use the preexisting list/iterator/is_a function above.
    # FIXME: fix their specs so that floor and ceiling start from 0 instead of 1...
    # FIXME: integrate them as methods of this class
    def build_args(self):
        """
        Returns a list of arguments that can be passed into the prexisting
        first,next,is_a, ... functions in this module.

        n is currently not included in this list.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: C.build_args()
            [3,
             3,
             <function <lambda> at 0x...>,
             <function <lambda> at 0x...>,
             -Infinity,
             +Infinity]

        """
        return [self.min_length, self.max_length,
                lambda i: self.floor(i-1), lambda i: self.ceiling(i-1),
                self.min_slope, self.max_slope]

    def first(self):
        """
        Returns the lexicographically maximal element in this
        combinatorial class.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: C.first()
            [2, 0, 0]
        """
        return self._element_constructor_(first(self.n_range[0], *(self.build_args())))

    def __iter__(self):
        """
        Returns an iterator for the elements of this combinatorial class.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: list(C) #indirect doctest
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

        """
        args = self.build_args()
        for n in self.n_range:
            l =  first(n, *args)
            while l is not None:
                yield self._element_constructor_(l)
                l = next(l, *args)

    def count(self):
        """
        Default brute force implementation of count by iteration
        through all the objects.

        Note that this skips the call to _element_constructor, unlike
        the default implementation from CombinatorialClass

        TODO: do the iteration in place to save on copying time
        """
        args = self.build_args()
        c = ZZ(0)
        for n in self.n_range:
            l =  first(n, *args)
            while l is not None:
                c += 1
                l = next(l, *args)
        return c

    def __contains__(self, v):
        """
        Returns True if and only if v is in this combinatorial class.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: [2, 0, 0] in C
            True
            sage: [2, 0] in C
            False
            sage: [3, 0, 0] in C
            False
            sage: all(v in C for v in C)
            True
        """
        return type(v) is type([]) and is_a(v, *(self.build_args())) and sum(v) in self.n_range
