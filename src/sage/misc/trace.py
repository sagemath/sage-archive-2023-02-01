"""
Interactively tracing execution of a command
"""

def trace(code, preparse=True):
    r"""
    Evaluate \sage code using the interactive tracer and return the
    result.  The string \var{code} must be a valid expression enclosed in
    quotes (no assignments -- the result of the expression is returned).

    INPUT:
        code -- str
        preparse -- bool (default: True); if True, run expression
                    through the Sage preparser.

    REMARKS: This function is extremely powerful!  For example, if you
    want to step through each line of execution of, e.g., \samp{factor(100)},
    type

        sage: trace("factor(100)")             # not tested

    then at the (Pdb) prompt type \kbd{s} (or \kbd{step}), then press return
    over and over to step through every line of Python that is called
    in the course of the above computation.   Type \kbd{?} at any time for
    help on how to use the debugger (e.g., \kbd{l} lists 11 lines around
    the current line; \kbd{bt} gives a back trace, etc.).

    Setting a break point: If you have some code in a file and would
    like to drop into the debugger at a given point, put the following
    code at that point in the file:

             \code{import pdb; pdb.set_trace()}

    For an article on how to use the Python debuger, see
       \url{http://www.onlamp.com/pub/a/python/2005/09/01/debugger.html}
    """

    import IPython.Debugger
    pdb = IPython.Debugger.Pdb()

    import IPython.ipapi
    _ip = IPython.ipapi.get()

    import preparser

    code = preparser.preparse(code)
    return pdb.runeval(code, _ip.user_ns)

    # this could also be useful; it drops
    # us into a debugger in an except block:
    #     import pdb; pdb.post_mortem(sys.exc_info()[2])
