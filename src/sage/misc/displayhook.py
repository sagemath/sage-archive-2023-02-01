"""
Implements a displayhook for Sage. The main improvement over the default
displayhook is a new facility for displaying lists of matrices in an easier to
read format.

AUTHORS:

- Bill Cauchois (2009): initial version
"""

import IPython, sys, __builtin__
from sage.matrix.matrix import is_Matrix

# This is used to wrap lines when printing "tall" lists.
MAX_COLUMN = 70

def _check_tall_list_and_print(out_stream, the_list):
    """
    First check whether a list is "tall" -- whether the reprs of the elements of
    the list will span multiple lines and cause the list to be printed awkwardly.
    If not, this function returns False and does nothing; you should revert back
    to the normal method for printing an object (its repr). If so, return True and
    print the list in the special format. Note that the special format isn't just
    for matrices. Any object with a multiline repr will be formatted.

    INPUT:

    - ``out_stream`` - The output stream to use.

    - ``the_list`` - The list (or a tuple).

    TESTS::

        sage: import sage.misc.displayhook, sys

    We test _check_tall_list_and_print() indirectly by calling print_obj() on
    a list of matrices::

        sage: sage.misc.displayhook.print_obj(sys.stdout, \
                [matrix([[1, 2, 3, 4], [5, 6, 7, 8]]) for i in xrange(7)])
        [
        [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]
        [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8],
        <BLANKLINE>
        [1 2 3 4]
        [5 6 7 8]
        ]
    """
    # For every object to be printed, split its repr on newlines and store the
    # result in this list.
    split_reprs = []
    tall = False
    for elem in the_list:
        split_reprs.append(`elem`.split('\n'))
        if len(split_reprs[-1]) > 1:
            # Meanwhile, check to make sure the list is actually "tall".
            tall = True
    if not tall:
        return False
    # Figure out which type of parenthesis to use, based on the type of the_list.
    if isinstance(the_list, tuple):
        parens = '()'
    elif isinstance(the_list, list):
        parens = '[]'
    else:
        raise TypeError, 'expected list or tuple'

    # running_lines is a list of lines, which are stored as lists of strings
    # to be joined later. For each split repr, we add its lines to the
    # running_lines array. When current_column exceeds MAX_COLUMN, process
    # and output running_lines using _print_tall_list_row.
    running_lines = [[]]
    current_column = 0
    print >>out_stream, parens[0]
    for split_repr in split_reprs:
        width = max(len(x) for x in split_repr)
        if current_column + width > MAX_COLUMN and not (width > MAX_COLUMN):
            _print_tall_list_row(out_stream, running_lines)
            running_lines = [[]]
            current_column = 0
        current_column += width + 2
        # Add the lines from split_repr to the running_lines array. It may
        # be necessary to add or remove lines from either one so that the
        # number of lines matches up.
        for i in xrange(len(running_lines), len(split_repr)):
            running_lines.insert(0, [' ' * len(x) for x in running_lines[-1]])
        line_diff = len(running_lines) - len(split_repr)
        for i, x in enumerate(split_repr):
            running_lines[i + line_diff].append(x.ljust(width))
        for i in xrange(line_diff):
            running_lines[i].append(' ' * width)
    # Output any remaining entries.
    if len(running_lines[0]) > 0:
        _print_tall_list_row(out_stream, running_lines, True)
    print >>out_stream, parens[1]
    return True

# This helper function for _print_tall_list processes and outputs the
# contents of the running_lines array.
def _print_tall_list_row(out_stream, running_lines, last_row=False):
    for i, line in enumerate(running_lines):
        if i + 1 != len(running_lines):
            sep, tail = '  ', ''
        else:
            # The commas go on the bottom line of this row.
            sep, tail = ', ', '' if last_row else ','
        print >>out_stream, sep.join(line) + tail
    # Separate rows with a newline to make them stand out.
    if not last_row:
        print >>out_stream

def print_obj(out_stream, obj):
    """
    Print an object. This function is used internally by the displayhook.

    EXAMPLES::

        sage: import sage.misc.displayhook, sys

    For most objects, printing is done simply using their repr::

        sage: sage.misc.displayhook.print_obj(sys.stdout, 'Hello, world!')
        'Hello, world!'
        sage: sage.misc.displayhook.print_obj(sys.stdout, (1, 2, 3, 4))
        (1, 2, 3, 4)

    We demonstrate the special format for lists of matrices::

        sage: sage.misc.displayhook.print_obj(sys.stdout, \
                [matrix([[1], [2]]), matrix([[3], [4]])])
        [
        [1]  [3]
        [2], [4]
        ]
    """
    # We only apply the special formatting to lists (or tuples) where the first
    # element is a matrix. This should cover most cases.
    if isinstance(obj, (tuple, list)):
        if len(obj) > 0 and is_Matrix(obj[0]):
            if _check_tall_list_and_print(out_stream, obj):
                return
    print >>out_stream, `obj`

def result_display(ip_self, obj):
    """
    This function implements the ``result_display`` hook for IPython.
    """
    # IPython's default result_display() uses the IPython.genutils.Term.cout stream.
    # See also local/lib/python2.6/site-packages/IPython/hooks.py.
    print_obj(IPython.genutils.Term.cout, obj)

def displayhook(obj):
    """
    This function adheres to the displayhook protocol described in `PEP 217`_.
    In order to mimic the behavior of the default displayhook, we update the
    variable ``__builtin__._`` every time an object is printed.

    .. _`PEP 217`: http://www.python.org/dev/peps/pep-0217/

    TESTS::

        sage: import sage.misc.displayhook, sys
        sage: sage.misc.displayhook.displayhook(1)
        1
        sage: print _
        1
        sage: sage.misc.displayhook.displayhook(None)
        sage: print _
        1
    """
    if obj is None:
        return
    __builtin__._ = None
    print_obj(sys.stdout, obj)
    __builtin__._ = obj

def install():
    """
    Install the new displayhook, so that subsequent output from the interpreter
    will be preprocessed by the mechanisms in this module.
    """
    # First, try to install the hook using the IPython hook API.
    ipapi = IPython.ipapi.get()
    if ipapi:
        ipapi.set_hook('result_display', result_display)
    else:
        # In certain modes where IPython is not in use, it is necessary to fall
        # back to setting Python's sys.displayhook.
        sys.displayhook = displayhook
