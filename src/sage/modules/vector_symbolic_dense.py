"""
Vectors over the symbolic ring.

Implements vectors over the symbolic ring.  Currently, this class only
provides methods for the simplification of symbolic vectors, as this
functionality was needed during the development of Trac #10132.  In the
long run, this class could be extended along the lines of
``sage.matrix.matrix_symbolic_dense``.


AUTHOR:

    -- Joris Vankerschaver (2011-05-15)

EXAMPLES::

    sage: x, y = var('x, y')
    sage: u = vector([sin(x)^2 + cos(x)^2, log(2*y) + log(3*y)]); u
    (sin(x)^2 + cos(x)^2, log(2*y) + log(3*y))
    sage: type(u)
    <class 'sage.modules.vector_symbolic_dense.Vector_symbolic_dense'>
    sage: u.simplify_full()
    (1, log(6) + 2*log(y))

TESTS::

    sage: u = vector(SR, [sin(x^2)])
    sage: loads(dumps(u)) == u
    True

"""

#*****************************************************************************
#       Copyright (C) 2011 Joris Vankerschaver (jv@caltech.edu)
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

import free_module_element
from sage.symbolic.ring import SR


class Vector_symbolic_dense(free_module_element.FreeModuleElement_generic_dense):

    def simplify_full(self):
        """
        Applies :meth:`simplify_full` to the entries of self.

        EXAMPLES::

            sage: u = vector([sin(x)^2 + cos(x)^2, 1])
            sage: u.simplify_full()
            (1, 1)
            sage: v = vector([log(exp(x))])
            sage: v.simplify_full()
            (x)

        """
        return (SR**len(self))([fun.simplify_full() for fun in self])
