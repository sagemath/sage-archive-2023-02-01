"""
Grab signal handlers back from PARI or other C libraries
"""

import signal
import sage.ext.sig

def my_sigint(x, n):
    sage.ext.sig.get_bad_sigs()
    raise KeyboardInterrupt

def my_sigfpe(x, n):
    sage.ext.sig.get_bad_sigs()
    raise RuntimeError, "A floating point exception occured."

def get_sigs():
    signal.signal(signal.SIGINT, my_sigint)
    signal.signal(signal.SIGABRT, my_sigint)
    signal.signal(signal.SIGFPE, my_sigfpe)
    signal.signal(signal.SIGALRM, my_sigint)


