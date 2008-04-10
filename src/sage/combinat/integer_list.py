r"""
Tools for generating lists of integers in lexicographic order.
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
from sage.rings.infinity import PlusInfinity
import __builtin__

infinity = PlusInfinity()

def first(n, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Returns the lexicographically smallest valid composition of n
    satisfying the conditions.

    Preconditions:
    //  - minslope < maxslope
    //  - floor and ceiling need to satisfy the slope constraints
    //    e.g. be obtained from comp2floor or comp2ceil
    //  - floor must be below ceiling to ensure the existence a
    //    valid composition

    TESTS:
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

    TESTS:
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

    TESTS:
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
    Returns the next integer list after comp that satisfies the contraints.

    EXAMPLES:
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
    EXAMPLES:
        sage: from sage.combinat.integer_list import iterator
        sage: IV = IntegerVectors(2,3,min_slope=0)
        sage: params = IV._parameters()
        sage: list(iterator(2,*params))
        [[0, 1, 1], [0, 0, 2]]
    """
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
    EXAMPLES:
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
    return __builtin__.list(iterator(n, min_length, max_length, floor, ceiling, min_slope, max_slope))

def upper_regular(comp, min_slope, max_slope):
    """
    Returns the uppest regular composition above comp.

    TESTS:
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

    EXAMPLES:
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

    EXAMPLES:
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
    Compute a coarse upper bound on the size of a vector satisfying
    the constraints.

    TESTS:
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
    Returns True if comp meets the constraints imposed by the arguments.

    EXAMPLES:
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
