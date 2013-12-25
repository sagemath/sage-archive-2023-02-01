r"""
Errors related to trac access

This module provides exception classes to report errors in trac sessions via
:class:`trac_interface.TracInterface`.

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

class TracError(RuntimeError):
    r"""
    Abstract base class for an error which ocurred during a trac session.

    EXAMPLES::

        sage: from sage.dev.trac_error import TracError
        sage: TracError()
        TracError()
    """
    pass

class TracConnectionError(TracError):
    r"""
    Trac connection error.

    EXAMPLES::

        sage: from sage.dev.trac_error import TracConnectionError
        sage: TracConnectionError()
        TracConnectionError('Connection to trac server failed.',)
    """
    def __init__(self):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.trac_error import TracConnectionError
            sage: type(TracConnectionError())
            <class 'sage.dev.trac_error.TracConnectionError'>

        """
        TracError.__init__(self, "Connection to trac server failed.")

class TracAuthenticationError(TracError):
    r"""
    Trac authentication error.

    EXAMPLES::

        sage: from sage.dev.trac_error import TracAuthenticationError
        sage: TracAuthenticationError()
        TracAuthenticationError('Authentication with trac server failed.',)
    """
    def __init__(self):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.trac_error import TracAuthenticationError
            sage: type(TracAuthenticationError())
            <class 'sage.dev.trac_error.TracAuthenticationError'>

        """
        TracError.__init__(self, "Authentication with trac server failed.")

class TracInternalError(TracError):
    r"""
    Error to indicate that the XML-RPC interface of trac returned an error.

    EXAMPLES::

        sage: from sage.dev.trac_error import TracInternalError
        sage: import xmlrpclib
        sage: raise TracInternalError(xmlrpclib.Fault(403, 
        ....:     "TICKET_CREATE privileges are required to perform this operation."
        ....:     " You don't have the required permissions."))
        Traceback (most recent call last):
        ...
        TracInternalError: <Fault 403: "TICKET_CREATE privileges are required to 
        perform this operation. You don't have the required permissions.">
    """
    def __init__(self, fault):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.trac_error import TracInternalError
            sage: import xmlrpclib
            sage: type(TracInternalError(xmlrpclib.Fault(403, 
            ....:      "TICKET_CREATE privileges are required to perform this operation."
            ....:      " You don't have the required permissions.")))
            <class 'sage.dev.trac_error.TracInternalError'>
        """
        self._fault = fault
        self.faultCode = fault.faultCode
        TracError.__init__(self, self._fault)
