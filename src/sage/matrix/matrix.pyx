"""
Abstract base class for matrices

For design documentation see matrix/docs.py.
"""

################################################################################
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL).
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
################################################################################

from sage.misc.superseded import deprecation
deprecation(24096, "the module sage.modules.module_element is deprecated, import from sage.structure.element instead")

def is_Matrix(x):
    """
    EXAMPLES::

        sage: from sage.structure.element import is_Matrix
        sage: is_Matrix(0)
        False
        sage: is_Matrix(matrix([[1,2],[3,4]]))
        True
    """
    return isinstance(x, Matrix)
