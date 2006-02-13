###############################################################
# Interactive debugger

import pdb

import preparser


def trace(code, preparse=True):
    """
    Evaluate SAGE code using the interactive tracer and return the
    result.  The string code must be a valid expression enclosed in
    quotes (no assignments -- the result of the expression is returned).

    INPUT:
        code -- str
        preparse -- bool (default: True); if True, run expression
                    through the SAGE preparser.

    REMARKS: This function is extremely powerful!  For example, if you
    want to step through each line of execution of, e.g.,"factor(100)",
    type

        sage.: trace("factor(100)")

    then at the (Pdb) prompt type "s" (or "step"), then press return
    over and over to step through every line of Python that is called
    in the course of the above computation.   Type ? at any time for
    help on how to use the debugger (e.g., "l" lists 11 lines around
    the current line; "bt" gives a back trace, etc.).

    SETTING A BREAK POINT: If you have some code in a file and would
    like to drop into the debugger at a given point, put the following
    code at that point in the file: \code{import pdb; pdb.set_trace()}
    For an article on how to use the Python debuger, see
       http://www.onlamp.com/pub/a/python/2005/09/01/debugger.html
    """
    code = preparser.preparse(code)
    return pdb.runeval(code)

    # this could also be useful; it drops
    # us into a debugger in an except block:
    #     import pdb; pdb.post_mortem(sys.exc_info()[2])

