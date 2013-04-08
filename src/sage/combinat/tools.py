r"""
Tools
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

def transitive_ideal(f, x):
    """
    Given an initial value x and a successor function f, return a list
    containing x and all of its successors. The successor function
    should return a list of all the successors of f.

    Note that if x has an infinite number of successors,
    transitive_ideal won't return.

    EXAMPLES::

        sage: f = lambda x: [x-1] if x > 0 else []
        sage: sage.combinat.tools.transitive_ideal(f, 10)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    """

    #succ_rec = lambda g: map(succ_rec, f(g))
    #succ_rec(x)
    result = [x]
    known  = [x]
    todo   = [x]

    while todo != []:
        y = todo[0]
        todo.remove(y)
        for z in f(y):
            if z == [] or z in known:
                continue
            result.append(z)
            known.append(z)
            todo.append(z)


    result.sort()
    return result
