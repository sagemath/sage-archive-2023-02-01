r"""
Interface between Python and c_lib.

This allows Python code to access a few parts of c_lib.  This is not
needed for Cython code, since such code can access c_lib directly.


AUTHORS:

- Jeroen Demeyer (2010-10-13): initial version

"""
#*****************************************************************************
#       Copyright (C) 2011 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include '../ext/stdsage.pxi'
include '../ext/interrupt.pxi'

def _init_csage():
    """
    Call init_csage() and enable interrupts.

    This is normally done exactly once during Sage startup from
    sage/all.py
    """
    # Set the Python-level interrupt handler. When a SIGINT occurs,
    # this will not be called directly. Instead, a SIGINT is caught by
    # the libcsage (c_lib) interrupt handler. If it happens during pure
    # Python code (not within sig_on()/sig_off()), the handler will set
    # Python's interrupt flag. Python regularly checks this and will
    # call its interrupt handler (which is the one we set now). This
    # handler issues a sig_check() which finally raises the
    # KeyboardInterrupt exception.
    import signal
    signal.signal(signal.SIGINT, sage_python_check_interrupt)

    init_csage()


def sig_on_count():
    """
    Return _signals.sig_on_count.

    This is normally zero, but increased by every sig_on() or sig_str()
    and decreased by every sig_off(). In pure Python code, this should
    always be zero.

    EXAMPLES::

        sage: sig_on_count()
        0
    """
    return _signals.sig_on_count


def sage_python_check_interrupt(sig, frame):
    """
    Python-level interrupt handler for interrupts raised in Python
    code. This simply delegates to the interrupt handling code in
    libcsage (c_lib).
    """
    sig_check()
