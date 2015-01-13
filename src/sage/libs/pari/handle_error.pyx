"""
Functions for handling PARI errors

AUTHORS:

- Peter Bruin (September 2013): initial version (:trac:`9640`)

- Jeroen Demeyer (January 2015): use ``cb_pari_err_handle`` (:trac:`14894`)

"""

#*****************************************************************************
#       Copyright (C) 2013 Peter Bruin
#       Copyright (C) 2015 Jeroen Demeyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"
include "decl.pxi"

from cpython cimport PyErr_Occurred
from pari_instance cimport PariInstance


cdef void _pari_init_error_handling():
    """
    Set up our code for handling PARI errors.

    TESTS::

        sage: try:
        ....:     p = pari.polcyclo(-1)
        ....: except PariError as e:
        ....:     print e.errtext()
        domain error in polcyclo: index <= 0

    Warnings still work just like in GP::

        sage: pari('warning("test")')
          ***   user warning: test
    """
    global cb_pari_err_handle
    global cb_pari_err_recover
    cb_pari_err_handle = _pari_err_handle
    cb_pari_err_recover = _pari_err_recover


cdef int _pari_err_handle(GEN E) except 0:
    """
    Convert a PARI error into a Sage exception, unless the error was
    a stack overflow, in which case we enlarge the stack.

    This function is a callback from the PARI error handler.

    EXAMPLES::

        sage: pari('error("test")')
        Traceback (most recent call last):
        ...
        PariError: error: user error: test
        sage: pari(1)/pari(0)
        Traceback (most recent call last):
        ...
        PariError: impossible inverse in gdiv: 0

    """
    from sage.libs.pari.all import PariError, pari
    cdef PariInstance P = pari

    cdef long errnum = E[1]
    if errnum == e_STACK:
        # PARI is out of memory.  We double the size of the PARI stack
        # and retry the computation.
        P.allocatemem(silent=True)
        return 0

    sig_block()
    cdef char* s = pari_err2str(E)
    try:
        pari_error_string = s.decode('ascii')
    finally:
        pari_free(s)
        sig_unblock()

    s = closure_func_err()
    if s is not NULL:
        pari_error_string = s.decode('ascii') + ": " + pari_error_string

    raise PariError(errnum, pari_error_string, P.new_gen_noclear(E))


cdef void _pari_err_recover(long errnum):
    """
    Reset the error string and jump back to ``sig_on()``, either to
    retry the code (in case of no error) or to make the already-raised
    exception known to Python.

    TEST:

    Perform a computation that requires doubling the default stack
    several times::

        sage: pari.allocatemem(2^12)
        PARI stack size set to 4096 bytes
        sage: x = pari('2^(2^26)')
        sage: x == 2^(2^26)
        True

    """
    if not PyErr_Occurred():
        # No exception was raised => retry the computation starting
        # from sig_on(). This can happen if we successfully enlarged the
        # PARI stack in _pari_handle_exception().
        sig_retry()
    else:
        # An exception was raised.  Jump to the signal-handling code
        # which will cause sig_on() to see the exception.
        sig_error()
