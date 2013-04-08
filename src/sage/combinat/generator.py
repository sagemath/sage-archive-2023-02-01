"""
Generators
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

def concat(gens):
    r"""
    Returns a generator that is the concatenation of the generators in
    the list.

    EXAMPLES::

        sage: list(sage.combinat.generator.concat([[1,2,3],[4,5,6]]))
        [1, 2, 3, 4, 5, 6]
    """

    for gen in gens:
        for element in gen:
            yield element


def map(f, gen):
    """
    Returns a generator that returns f(g) for g in gen.

    EXAMPLES::

        sage: f = lambda x: x*2
        sage: list(sage.combinat.generator.map(f,[4,5,6]))
        [8, 10, 12]
    """
    for element in gen:
        yield f(element)

def element(element, n = 1):
    """
    Returns a generator that yield a single element n times.

    EXAMPLES::

        sage: list(sage.combinat.generator.element(1))
        [1]
        sage: list(sage.combinat.generator.element(1, n=3))
        [1, 1, 1]
    """
    for i in range(n):
        yield element

def select(f, gen):
    """
    Returns a generator for all the elements g of gen such that f(g) is
    True.

    EXAMPLES::

        sage: f = lambda x: x % 2 == 0
        sage: list(sage.combinat.generator.select(f,range(7)))
        [0, 2, 4, 6]
    """
    for element in gen:
        if f(element):
            yield element


def successor(initial, succ):
    """
    Given an initial value and a successor function, yield the initial
    value and each following successor. The generator will continue to
    generate values until the successor function yields None.

    EXAMPLES::

        sage: f = lambda x: x+1 if x < 10 else None
        sage: list(sage.combinat.generator.successor(0,f))
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    """
    yield initial

    s = succ(initial)
    while s:
        yield s
        s = succ(s)


