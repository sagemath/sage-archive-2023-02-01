"""
Backtracking
"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
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
class GenericBacktracker(object):
    def __init__(self, initial_data, initial_state):
        """
        EXAMPLES::

            sage: from sage.combinat.backtrack import GenericBacktracker
            sage: p = GenericBacktracker([], 1)
            sage: loads(dumps(p))
            <sage.combinat.backtrack.GenericBacktracker object at 0x...>
        """
        self._initial_data = initial_data
        self._initial_state = initial_state

    def iterator(self):
        """
        EXAMPLES::

            sage: from sage.combinat.permutation import PatternAvoider
            sage: p = PatternAvoider(4, [[1,3,2]])
            sage: len(list(p.iterator()))
            14
        """
        #Initialize the stack of generators with the initial data.
        #The generator in stack[i] is a generator for the i^th level
        #of the search tree.
        stack = []
        stack.append(self._rec(self._initial_data, self._initial_state))

        done = False
        while not done:
            #Try to get the next object in this level
            try:
                obj, state, yld = stack[-1].next()
            except StopIteration:
                #If there are no more, go back up the tree
                #We also need to check if we've exhausted all
                #possibilities
                stack.pop()
                done = len(stack) == 0
                continue

            #If the return state is None, then obj is a leaf
            #of the search tree.  If yld is True, then obj
            #should be yielded.
            if yld is True:
                yield obj
            if state is not None:
                stack.append( self._rec(obj, state) )

    __iter__ = iterator


