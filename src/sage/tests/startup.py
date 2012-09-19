r"""

Test that on sage.math (a fast computer), the Sage startup time (from
a warm cache) has not got out of hand.::

    sage: if os.uname()[1] == 'sage.math.washington.edu':
    ...     startup = float(os.popen("sage --startuptime>/dev/null; sage --startuptime|grep 'Total time'").readlines()[0].split()[-1].rstrip('ms'))/1000
    ...     if startup < 2:
    ...         print True
    ...     else:
    ...         print "'sage --startuptime' took %s seconds." % startup
    ... else: print True   # nothing when not on sage.math
    True

Ensure that certain modules are not loaded on startup::

    sage: sage0("'numpy' in sys.modules")
    False
"""
