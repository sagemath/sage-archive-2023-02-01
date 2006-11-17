r"""
Base class for parent objects
"""

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cimport sage_object

include '../ext/python_object.pxi'
include '../ext/python_bool.pxi'

def is_Parent(x):
    """
    Return True if x is a parent object, i.e., derives from
    sage.structure.parent.Parent and False otherwise.

    EXAMPLES:
        sage: is_Parent(2/3)
        False
        sage: is_Parent(ZZ)
        True
        sage: is_Parent(Primes())
        True
    """
    return PyBool_FromLong(PyObject_TypeCheck(x, Parent))

cdef class Parent(sage_object.SageObject):
    """
    Parents are the SAGE/mathematical analogues of container objects in computer science.
    """
    #############################################################################
    # Containment testing
    #############################################################################

    def __contains__(self, x):
        r"""
        True if any coercion of $x$ into self is possible (this is not
        necessarily canonical coercion).

        EXAMPLES:
            sage: 2 in Integers(7)
            True
            sage: 2 in ZZ
            True
            sage: Integers(7)(3) in ZZ
            True
            sage: 3/1 in ZZ
            True
            sage: 5 in QQ
            True
            sage: I in RR
            False
        """
        try:
            self(x)
        except TypeError:
            return False
        return True




