"""
IPython's timeit magic functionality ported to Sage.

This is an implemention of nice timeit functionality, like the
%timeit magic command in IPython.  To use it, use the timeit
command.

AUTHOR:
    -- William Stein, based on code copied from Fernand Perez's Ipython
"""

import timeit as timeit_, time, math, preparser, interpreter

def sage_timeit(stmt, globals, preparse=None,
                   number = 0, repeat = 3, precision = 3):
    """
    INPUT:
        stmt -- a text string which may
        globals -- evaluate stmt in the context of the globals dictionary
        preparse -- (default: use global preparser default) if True preparse
                    stmt using the Sage preparser.
    EXAMPLES:
        sage: from sage.misc.sage_timeit import sage_timeit
        sage: sage_timeit('3^100000', globals(), preparse=True, number=50)      # random output
        '50 loops, best of 3: 1.97 ms per loop'
        sage: sage_timeit('3^100000', globals(), preparse=False, number=50)     # random output
        '50 loops, best of 3: 67.1 ns per loop'
        sage: a = 10
        sage: sage_timeit('a^2', globals(), number=50)                            # random output
        '50 loops, best of 3: 4.26 us per loop'

    It's usually better to use the timeit object, usually:
        sage: timeit('10^2', number=50)
        50 loops, best of 3: ... per loop

    The input expression can contain newlines:
        sage: from sage.misc.sage_timeit import sage_timeit
        sage: sage_timeit("a = 2\nb=131\nfactor(a^b-1)", globals(), number=10)
        '10 loops, best of 3: ... per loop'

    Test to make sure that timeit behaves well with output:
        sage: timeit("print 'Hi'", number=50)
        50 loops, best of 3: ... per loop


    Make sure that garbage collection is renabled after an exception
    occurs in timeit.

    TESTS:
        sage: def f(): raise ValueError
        sage: import gc
        sage: gc.isenabled()
        True
        sage: timeit("f()")
        Traceback (most recent call last):
        ...
        ValueError
        sage: gc.isenabled()
        True

    """
    number=int(number)
    repeat=int(repeat)
    precision=int(precision)
    if preparse is None:
        preparse = interpreter.do_preparse
    if preparse:
        stmt = preparser.preparse(stmt)
    if stmt == "":
        return ''

    units = ["s", "ms", "\xc2\xb5s", "ns"]
    scaling = [1, 1e3, 1e6, 1e9]

    timer = timeit_.Timer()

    # this code has tight coupling to the inner workings of timeit.Timer,
    # but is there a better way to achieve that the code stmt has access
    # to the shell namespace?

    src = timeit_.template % {'stmt': timeit_.reindent(stmt, 8),
                             'setup': "pass"}
    code = compile(src, "<magic-timeit>", "exec")
    ns = {}
    exec code in globals, ns
    timer.inner = ns["inner"]


    try:
        import sys
        f = sys.stdout
        sys.stdout = open('/dev/null', 'w')

        if number == 0:
            # determine number so that 0.2 <= total time < 2.0
            number = 1
            for i in range(1, 5):
                number *= 5
                if timer.timeit(number) >= 0.2:
                    break

        best = min(timer.repeat(repeat, number)) / number

    finally:
        sys.stdout.close()
        sys.stdout = f
        import gc
        gc.enable()

    if best > 0.0:
        order = min(-int(math.floor(math.log10(best)) // 3), 3)
    else:
        order = 3
    return "%d loops, best of %d: %.*g %s per loop" % (number, repeat,
                                                      precision,
                                                      best * scaling[order],
                                                      units[order])


