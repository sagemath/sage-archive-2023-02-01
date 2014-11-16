"""
The display hook for plain Python

This is not used directly in interactive Sage (where we use the
IPython system for display hooks).  This class provides a way to use
the Sage "plain text" display formatting when not using interactive
Sage, for example when running doctests.
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import __builtin__


from sage.repl.display.formatter import SageDoctestTextFormatter


class DoctestDisplayHook(object):

    def __init__(self):
        """
        Python constructor

        EXAMPLES::

            sage: from sage.repl.display.python_hook import DoctestDisplayHook
            sage: d = DoctestDisplayHook()
            sage: d(set([1, 2, 3]))       # Sage commandline output
            {1, 2, 3}
            sage: print(set([1, 2, 3]))   # Plain Python output
            set([1, 2, 3])
        """
        self.formatter = SageDoctestTextFormatter()

    def __call__(self, obj):
        """
        Format the object using Sage's formatting, or format it using the old
        display hook if Sage does not want to handle the object.

        EXAMPLES::

            sage: from sage.repl.display.python_hook import DoctestDisplayHook
            sage: d = DoctestDisplayHook()
            sage: d((identity_matrix(3), identity_matrix(3)))
            (
            [1 0 0]  [1 0 0]
            [0 1 0]  [0 1 0]
            [0 0 1], [0 0 1]
            )
        """
        if obj is None:
            return
        print(self.formatter(obj))
        __builtin__._ = obj
