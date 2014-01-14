# -*- coding: utf-8 -*-
r"""
Implements a displayhook for Sage.

This displayhook has two new facilities, by default the displayhook contains a
new facility for displaying lists of matrices in an easier to read format::

    sage: [identity_matrix(i) for i in range(2,5)]
    [
                    [1 0 0 0]
           [1 0 0]  [0 1 0 0]
    [1 0]  [0 1 0]  [0 0 1 0]
    [0 1], [0 0 1], [0 0 0 1]
    ]

This facility uses :meth:`_repr_` (and a simple string) to try do a nice read
format (see :meth:`sage.structure.parent._repr_option` for details).

With this displayhook there exists an other way for displaying object and more
generally, all sage expression as an ASCII art object::

    sage: from sage.misc.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell('%display ascii_art')
    sage: shell.run_cell('integral(x^2/pi^x, x)')
     / 2    2                      \  -x*log(pi)
    -\x *log (pi) + 2*x*log(pi) + 2/*e
    ---------------------------------------------
                         3
                      log (pi)
    sage: shell.run_cell("i = var('i')")
    sage: shell.run_cell('sum(i*x^i, i, 0, 10)')
        10      9      8      7      6      5      4      3      2
    10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x
    sage: shell.run_cell('StandardTableaux(4).list()')
    [
    [                                                                  1  4    1  3
    [                 1  3  4    1  2  4    1  2  3    1  3    1  2    2       2
    [   1  2  3  4,   2      ,   3      ,   4      ,   2  4,   3  4,   3   ,   4
    <BLANKLINE>
                1 ]
        1  2    2 ]
        3       3 ]
    ,   4   ,   4 ]
    sage: shell.run_cell('%display simple')

This other facility uses a simple `AsciiArt` object
(see :class:`sage.misc.ascii_art.AsciiArt` and
:meth:`sage.structure.parent._ascii_art_`).

AUTHORS:

- Bill Cauchois (2009): initial version
- Jean-Baptiste Priez <jbp@kerios.fr> (2013): ASCII art
- Volker Braun (2013): refactored into DisplayHookBase
"""

import sys, __builtin__

class ListFormatter(object):

    # This is used to wrap lines when printing "tall" lists.
    MAX_COLUMN = 70

    def _check_tall_list_and_format(self, the_list):
        """
        First check whether a list is "tall" -- whether the reprs of the
        elements of the list will span multiple lines and cause the list
        to be printed awkwardly.  If not, this function returns ``None`` and
        does nothing; you should revert back to the normal method for
        printing an object (its repr). If so, return the string in the
        special format. Note that the special format isn't just for
        matrices. Any object with a multiline repr will be formatted.

        INPUT:

        - ``the_list`` - The list (or a tuple).

        TESTS::

            sage: from sage.misc.displayhook import DisplayHookBase
            sage: dhb = DisplayHookBase()

        We test :meth:`_check_tall_list_and_format` indirectly by
        calling :meth:`simple_format_obj` on a list of matrices::

            sage: print dhb.simple_format_obj(
            ....:        [matrix([[1, 2, 3, 4], [5, 6, 7, 8]]) for i in xrange(7)])
            [
            [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]  [1 2 3 4]
            [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8], [5 6 7 8],
            <BLANKLINE>
            [1 2 3 4]
            [5 6 7 8]
            ]

        We return ``None`` if we don't have anything special to do::

            sage: dhb.simple_format_obj('one-line string')
            sage: dhb.simple_format_obj(matrix([[1,2,3]]))
        """
        # For every object to be printed, split its repr on newlines and store the
        # result in this list.
        split_reprs = []
        tall = False
        for elem in the_list:
            split_reprs.append(repr(elem).split('\n'))
            if len(split_reprs[-1]) > 1:
                # Meanwhile, check to make sure the list is actually "tall".
                tall = True
        if not tall:
            return None
        # Figure out which type of parenthesis to use, based on the type of the_list.
        if isinstance(the_list, tuple):
            parens = '()'
        elif isinstance(the_list, list):
            parens = '[]'
        else:
            raise TypeError('expected list or tuple')

        # running_lines is a list of lines, which are stored as lists of strings
        # to be joined later. For each split repr, we add its lines to the
        # running_lines array. When current_column exceeds MAX_COLUMN, process
        # and output running_lines using _print_tall_list_row.
        running_lines = [[]]
        current_column = 0
        s = [parens[0]]
        for split_repr in split_reprs:
            width = max(len(x) for x in split_repr)
            if current_column + width > self.MAX_COLUMN and not (width > self.MAX_COLUMN):
                s.extend(self._tall_list_row(running_lines))
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
            s.extend(self._tall_list_row(running_lines, True))
        s.append(parens[1])
        return "\n".join(s)

    def _tall_list_row(self, running_lines, last_row=False):
        """
        Helper for :meth:`_check_tall_list_and_format`

        This helper function processes and outputs the contents of the
        running_lines array.

        TESTS::

            sage: from sage.misc.displayhook import format_list
            sage: format_list._tall_list_row(['a   b', 'b  c', 'c'])
            ['a           b', 'b        c', 'c,', '']
        """
        s=[]
        for i, line in enumerate(running_lines):
            if i + 1 != len(running_lines):
                sep, tail = '  ', ''
            else:
                # The commas go on the bottom line of this row.
                sep, tail = ', ', '' if last_row else ','
            s.append(sep.join(line) + tail)
        # Separate rows with a newline to make them stand out.
        if not last_row:
            s.append("")
        return s

    def try_format_list(self, obj):
        """
        Format list/tuple.

        OUTPUT:

        A string representation or ``None``. The latter means that no
        Sage-specific formatting is defined and the default should be
        used.

        EXAMPLES::

            sage: from sage.misc.displayhook import format_list
            sage: format_list.try_format_list('Hello, world!')
            sage: format_list('Hello, world!')
            "'Hello, world!'"
        """
        ascii_art_repr = False
        for o in obj:
            try:
                ascii_art_repr = ascii_art_repr or o.parent()._repr_option('element_ascii_art')
            except (AttributeError, TypeError):
                pass
        if ascii_art_repr:
            return self._check_tall_list_and_format(obj)
        else:
            return None

    def __call__(self, obj):
        """
        Return a string formatting.

        This method is like :meth:`try_format_list` except that it
        will always return a string.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.misc.displayhook import format_list
            sage: format_list.try_format_list('Hello, world!')
            sage: format_list('Hello, world!')
            "'Hello, world!'"
        """
        s = self.try_format_list(obj)
        if s is None:
            return repr(obj)
        else:
            return s

format_list = ListFormatter()


class DisplayHookBase(object):

    def simple_format_obj(self, obj):
        """This function is used internally by the displayhook.

        We attempt to keep ascii art of list/tuple members intact as we
        print them. See :meth:`sage.structure.parent._repr_option` for
        details.

        OUTPUT:

        Return a string if we want to print it in a special way;
        otherwise, return ``None``.

        EXAMPLES::

            sage: from sage.misc.displayhook import DisplayHookBase
            sage: dhb = DisplayHookBase()

        For most objects, nothing is done (``None`` is returned):

            sage: dhb.simple_format_obj('Hello, world!')
            sage: dhb.simple_format_obj((1, 2, 3, 4))

        We demonstrate the special format for lists of matrices::

            sage: dhb.simple_format_obj(
            ....:     [matrix([[1], [2]]), matrix([[3], [4]])])
            '[\n[1]  [3]\n[2], [4]\n]'

        TESTS:

        In :trac:`14466` we override IPython's special printing of
        ``type`` objects and revert it to Python's standard string
        representation::

            sage: shell=sage.misc.interpreter.get_test_shell()
            sage: shell.displayhook(type)
            <type 'type'>
        """
        if isinstance(obj, type):
            return repr(obj)
        if isinstance(obj, (tuple, list)) and len(obj) > 0:
            return format_list.try_format_list(obj)
        else:
            return None

    def set_display(self, mode):
        r"""
        Select the text formatting method.

        INPUT:

        - ``mode`` -- string. One of ``simple``, ``ascii_art``, or ``typeset``.

        See :func:`simple_format_obj` or :func:`sage.misc.ascii_art.ascii_art`.

        TESTS::

            sage: [identity_matrix(i) for i in range(3,7)]
            [
                                             [1 0 0 0 0 0]
                                [1 0 0 0 0]  [0 1 0 0 0 0]
                     [1 0 0 0]  [0 1 0 0 0]  [0 0 1 0 0 0]
            [1 0 0]  [0 1 0 0]  [0 0 1 0 0]  [0 0 0 1 0 0]
            [0 1 0]  [0 0 1 0]  [0 0 0 1 0]  [0 0 0 0 1 0]
            [0 0 1], [0 0 0 1], [0 0 0 0 1], [0 0 0 0 0 1]
            ]
            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%display ascii_art')   # indirect doctest
            sage: shell.run_cell("i = var('i')")
            sage: shell.run_cell('sum(i*x^i, i, 0, 10)')
                10      9      8      7      6      5      4      3      2
            10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x
            sage: shell.run_cell('%display simple')
        """
        if mode not in ['simple', 'ascii_art', 'typeset']:
            raise ValueError('invalid mode set')
        self.mode = mode

    mode = 'simple'

    @property
    def simple(self):
        """
        Whether the mode is the "simple" (default) display.

        EXAMPLES::

            sage: sys.displayhook.simple
            True
            sage: sys.displayhook.ascii_art
            False
        """
        return self.mode == 'simple'

    @property
    def ascii_art(self):
        """
        Whether the mode is the ascii art display.

        EXAMPLES::

            sage: sys.displayhook.simple
            True
            sage: sys.displayhook.ascii_art
            False
        """
        return self.mode == 'ascii_art'

    @property
    def typeset(self):
        """
        Whether the mode is the notebook "Typeset" display.

        EXAMPLES::

            sage: sys.displayhook.simple
            True
            sage: sys.displayhook.typeset
            False
        """
        return self.mode == 'typeset'

    def try_format_obj(self, obj):
        """
        Format non-graphics object.

        OUTPUT:

        A string representation or ``None``. The latter means that no
        Sage-specific formatting is defined and the default should be
        used.

        TESTS::

            sage: from sage.misc.displayhook import DisplayHookBase
            sage: dhb = DisplayHookBase()
            sage: dhb.try_format_obj('Hello, world!')
        """
        if self.simple:
            return self.simple_format_obj(obj)
        if self.ascii_art:
            from sage.misc.ascii_art import ascii_art
            return ascii_art(obj)
        if self.typeset:
            from sage.misc.latex import pretty_print
            pretty_print(obj)
            return ''
        assert(False)

    def try_format_graphics(self, obj):
        """
        Format graphics object.

        OUTPUT:

        Boolean. Whether the object is graphics and was successfully
        displayed.

        TESTS::

            sage: from sage.misc.displayhook import DisplayHookBase
            sage: dhb = DisplayHookBase()
            sage: dhb.try_format_graphics('Hello, world!')
            False
        """
        from sage.structure.sage_object import SageObject
        if isinstance(obj, SageObject) and hasattr(obj, '_graphics_'):
            return obj._graphics_()
        return False


class DisplayHook(DisplayHookBase):
    """
    Display hook for Sage.

    This is not used directly in interactive Sage (where we use the
    IPython system for display hooks).  This class provides a way to
    use the Sage display formatting when not using interactive Sage.
    """
    def __init__(self, oldhook=sys.__displayhook__):
        """
        Set the old display hook (default to repr)

        EXAMPLES::

            sage: from sage.misc.displayhook import DisplayHook
            sage: def f(o): print repr(o)[:5], "..."
            sage: d = DisplayHook(f)
            sage: d(range(10))
            [0, 1 ...
        """
        self.oldhook = oldhook

    def __call__(self, obj):
        """
        Format the object using Sage's formatting, or format it using the old
        display hook if Sage does not want to handle the object.

        EXAMPLES::

            sage: from sage.misc.displayhook import DisplayHook
            sage: d = DisplayHook()
            sage: d((identity_matrix(3), identity_matrix(3)))
            (
            [1 0 0]  [1 0 0]
            [0 1 0]  [0 1 0]
            [0 0 1], [0 0 1]
            )
        """
        if self.try_format_graphics(obj):
            return
        s = self.try_format_obj(obj)
        if s is not None:
            print s
            __builtin__._ = obj
        else:
            self.oldhook(obj)


from IPython.core.formatters import PlainTextFormatter
class SagePlainTextFormatter(DisplayHookBase, PlainTextFormatter):
    r"""
    A replacement for the plain text formatter which can use two facilities:

    - correctly print lists of matrices or other objects (see
      :meth:`sage.structure.parent._repr_option`),
    - print ASCII art objects (like expressions) (see
      :meth:`sage.structure.parent._ascii_art_`).

    EXAMPLES::

        sage: from sage.misc.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.display_formatter.formatters['text/plain']
        <...displayhook.SagePlainTextFormatter object at 0x...>
        sage: shell.run_cell('a = identity_matrix(ZZ, 2); [a,a]')
        [
        [1 0]  [1 0]
        [0 1], [0 1]
        ]
    """
    def __call__(self, obj):
        r"""
        Computes the format data of ``result``.  If the
        :func:`sage.misc.displayhook.simple_format_obj` writes a string, then
        we override IPython's :class:`DisplayHook` formatting.

        EXAMPLES::

            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: fmt = shell.display_formatter.formatters['text/plain']
            sage: fmt
            <...displayhook.SagePlainTextFormatter object at 0x...>
            sage: shell.displayhook.compute_format_data(2)
            {u'text/plain': '2'}

            sage: a = identity_matrix(ZZ, 2)
            sage: shell.displayhook.compute_format_data([a,a])
            {u'text/plain': '[\n[1 0]  [1 0]\n[0 1], [0 1]\n]'}
            sage: fmt.set_display('ascii_art')
            sage: shell.displayhook.compute_format_data([a,a])
            {u'text/plain': [ [1 0]  [1 0] ]
                            [ [0 1], [0 1] ]}

            sage: i = var('i')
            sage: shell.displayhook.compute_format_data(sum(i*x^i, i, 0, 10))
            {u'text/plain':     10      9      8      7      6      5      4      3      2
                            10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x}
            sage: fmt.set_display('simple')
        """
        if self.try_format_graphics(obj):
            return ''
        s = self.try_format_obj(obj)
        if s is None:
            s = super(SagePlainTextFormatter, self).__call__(obj)
        return s


SPTextFormatter = None
