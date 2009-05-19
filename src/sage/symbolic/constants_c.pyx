"""
Wrapper around Pynac's constants
"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#                     2009 Mike Hansen <mhansen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################
from sage.symbolic.expression cimport new_Expression_from_GEx
from sage.symbolic.ring import SR

include "../libs/ginac/decl.pxi"

cdef extern from "pynac/constant.h":
    pass

cdef class PynacConstant:
    def __cinit__(self, name, texname, domain_string):
        """
        Creates a constant in Pynac.

        EXAMPLES:

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f
            foo

        Note that this just creates a 'constant' object and not an
        expression.  If you want to work with this constant, you'll
        have to use the result of :meth:`expression`::

            sage: foo = f.expression()
            sage: foo + 2
            foo + 2
        """
        cdef unsigned domain
        if domain_string == 'complex':
            domain = domain_complex
        elif domain_string == 'real':
            domain = domain_real
        elif domain_string == 'positive':
            domain = domain_positive
        else:
            raise ValueError

        self._name = name

        #For the constants explicitly defined in constant.cpp in the Pynac
        #library, we use those symbols.
        if self._name == "pi":
            self.pointer = <GConstant *>&g_Pi
        elif self._name == "catalan":
            self.pointer = <GConstant *>&g_Catalan
        elif self._name == "euler_gamma":
            self.pointer = <GConstant *>&g_Euler
        else:
            GConstant_construct(&self.object, name, texname, domain)
            self.pointer = &self.object

    def __dealloc__(self):
        if self.pointer == & self.object:
            GConstant_destruct(&self.object)

    def serial(self):
        """
        Returns the underlying Pynac serial for this constant.

        EXAMPLES::

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f.serial()  #random
            15
        """
        return int(self.pointer.get_serial())

    def name(self):
        """
        Returns the name of this constant.

        EXAMPLES::

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f.name()
            'foo'
        """
        return self._name

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real'); f
            foo
        """
        return self.name()

    def expression(self):
        """
        Returns this constant as an Expression.

        EXAMPLES::

            sage: from sage.symbolic.constants_c import PynacConstant
            sage: f = PynacConstant('foo', 'foo', 'real')
            sage: f + 2
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': '<type 'sage.symbolic.constants_c.PynacConstant'>' and 'Integer Ring'

            sage: foo = f.expression(); foo
            foo
            sage: foo + 2
            foo + 2
        """
        return new_Expression_from_GEx(SR, <GEx>(self.pointer[0]))
