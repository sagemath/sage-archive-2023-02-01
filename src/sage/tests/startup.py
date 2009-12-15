r"""

Test that on sage.math (a fast computer), the Sage startup time (from
a warm cache) has not got out of hand.::

    sage: if os.uname()[1] == 'sage.math.washington.edu':
    ...     print float(os.popen("sage -startuptime>/dev/null; sage -startuptime|grep sage.all").readlines()[0].split()[1]) < 1.5
    ... else: print True   # nothing when not on sage.math
    True

"""
