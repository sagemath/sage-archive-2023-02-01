r"""
User interface errors

This module provides subclasses of ``RuntimeError`` to indicate error
conditions in the user interface.

AUTHORS:

- Julian Rueth: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


class OperationCancelledError(RuntimeError):
    r"""
    Indicates that the user cancelled an interactive action, e.g., the user was
    asked to edit a file but left the editor with a non-zero exit code.

    EXAMPLES::

        sage: from sage.dev.user_interface_error import OperationCancelledError
        sage: raise OperationCancelledError("cancelled")
        Traceback (most recent call last):
        ...
        OperationCancelledError: cancelled
    """
    def __init__(self, reason):
        r"""
        Initialization.

        EXAMPLES::

            sage: from sage.dev.user_interface_error import OperationCancelledError
            sage: type(OperationCancelledError("cancelled"))
            <class 'sage.dev.user_interface_error.OperationCancelledError'>
        """
        RuntimeError.__init__(self, reason)
