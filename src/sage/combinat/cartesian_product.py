r"""
Cartesian Products
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

import sage.misc.prandom as rnd
import __builtin__
from combinat import CombinatorialClass


def CartesianProduct(*iters):
    """
    Returns the combinatorial class of the cartesian
    product of *iters.

    EXAMPLES:
        sage: cp = CartesianProduct([1,2], [3,4]); cp
        Cartesian product of [1, 2], [3, 4]
        sage: cp.list()
        [[1, 3], [1, 4], [2, 3], [2, 4]]

      Note that if you have a generator-type object that is returned
      by a function, then you should use IterableFunctionCall class
      defined in sage.combinat.misc.

        sage: def a(): yield 1; yield 2
        sage: def b(): yield 'a'; yield 'b'
        sage: CartesianProduct(a(), b()).list()
        [[1, 'a'], [1, 'b']]
        sage: from sage.combinat.misc import IterableFunctionCall
        sage: CartesianProduct(IterableFunctionCall(a), IterableFunctionCall(b)).list()
        [[1, 'a'], [1, 'b'], [2, 'a'], [2, 'b']]
    """
    return CartesianProduct_iters(*iters)

class CartesianProduct_iters(CombinatorialClass):
    def __init__(self, *iters):
        """
        TESTS:
            sage: import sage.combinat.cartesian_product as cartesian_product
            sage: cp = cartesian_product.CartesianProduct_iters([1,2],[3,4]); cp
            Cartesian product of [1, 2], [3, 4]
            sage: loads(dumps(cp)) == cp
            True
        """
        self.iters = iters

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: cp = CartesianProduct([1,2],[3,4])
            sage: [1,3] in cp
            True
            sage: [1,2] in cp
            False
            sage: [1, 3, 1] in cp
            False
        """
        try:
            return len(x) == len(self.iters) and all(x[i] in self.iters[i] for i in range(len(self.iters)))
        except (TypeError, IndexError):
            return False

    def __repr__(self):
        """
        EXAMPLES:
            sage: CartesianProduct(range(2), range(3))
            Cartesian product of [0, 1], [0, 1, 2]
        """
        return "Cartesian product of " + ", ".join(map(str, self.iters))

    def count(self):
        """
        Returns the number of elements in the cartesian product of
        everything in *iters.

        EXAMPLES:
            sage: CartesianProduct(range(2), range(3)).count()
            6
            sage: CartesianProduct(range(2), xrange(3)).count()
            6
            sage: CartesianProduct(range(2), xrange(3), xrange(4)).count()
            24
        """
        total = 1
        for it in self.iters:
            total *= len(it)
        return total


    def list(self):
        """
        Returns

        EXAMPLES:
            sage: CartesianProduct(range(3), range(3)).list()
            [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
            sage: CartesianProduct('dog', 'cat').list()
            [['d', 'c'],
             ['d', 'a'],
             ['d', 't'],
             ['o', 'c'],
             ['o', 'a'],
             ['o', 't'],
             ['g', 'c'],
             ['g', 'a'],
             ['g', 't']]
        """
        return [e for e in self.iterator()]


    def iterator(self):
        """
        An iterator for the elements in the cartesian product
        of the iterables *iters.

        From Recipe 19.9 in the Python Cookbook by Alex Martelli
        and David Ascher.

        EXAMPLES:
            sage: [e for e in CartesianProduct(range(3), range(3))]
            [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
            sage: [e for e in CartesianProduct('dog', 'cat')]
            [['d', 'c'],
             ['d', 'a'],
             ['d', 't'],
             ['o', 'c'],
             ['o', 'a'],
             ['o', 't'],
             ['g', 'c'],
             ['g', 'a'],
             ['g', 't']]
        """

        # visualize an odometer, with "wheels" displaying "digits"...:
        wheels = map(iter, self.iters)
        digits = [it.next() for it in wheels]
        while True:
            yield __builtin__.list(digits)
            for i in range(len(digits)-1, -1, -1):
                try:
                    digits[i] = wheels[i].next()
                    break
                except StopIteration:
                    wheels[i] = iter(self.iters[i])
                    digits[i] = wheels[i].next()
            else:
                break

    def random(self):
        """
        Returns a random element from the cartesian product
        of *iters.

        EXAMPLES:
            sage: CartesianProduct('dog', 'cat').random()
            ['d', 'a']
        """
        return list(map(rnd.choice, self.iters))
