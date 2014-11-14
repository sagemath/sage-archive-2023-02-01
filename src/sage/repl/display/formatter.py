# -*- coding: utf-8 -*-
r"""
IPython Displayhook Formatters

The classes in this module can be used as IPython displayhook
formatters. It has two main features, by default the displayhook
contains a new facility for displaying lists of matrices in an easier
to read format::

    sage: [identity_matrix(i) for i in range(2,5)]
    [
                    [1 0 0 0]
           [1 0 0]  [0 1 0 0]
    [1 0]  [0 1 0]  [0 0 1 0]
    [0 1], [0 0 1], [0 0 0 1]
    ]

This facility uses :meth:`_repr_` (and a simple string) to try do a nice read
format (see :meth:`sage.structure.parent.Parent._repr_option` for details).

With this displayhook there exists an other way for displaying object and more
generally, all sage expression as an ASCII art object::

    sage: from sage.repl.interpreter import get_test_shell
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
    sage: shell.run_cell('%display default')
    sage: shell.quit()

This other facility uses a simple
:class:`~sage.misc.ascii_art.AsciiArt` object (see and
:meth:`sage.structure.sage_object.SageObject._ascii_art_`).  """

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from IPython.core.formatters import DisplayFormatter, PlainTextFormatter, warn_format_error
from IPython.utils.py3compat import str_to_unicode, unicode_to_str

from sage.structure.sage_object import SageObject
from sage.repl.display.pretty_print import SagePrettyPrinter



class SageDisplayFormatter(DisplayFormatter):

    def __init__(self, *args, **kwds):
        """
        This is where the Sage rich objects are translated to IPython
        """
        super(SageDisplayFormatter, self).__init__(*args, **kwds)
        from sage.repl.rich_output.display_manager import get_display_manager
        self.dm = get_display_manager()
        from sage.repl.rich_output.backend_ipython import BackendIPython
        self.dm.check_backend_class(BackendIPython)

    def format(self, obj, include=None, exclude=None):
        """
        Use the Sage rich output instead of IPython

        EXAMPLES::

            sage: [identity_matrix(i) for i in range(3,7)]
            [
                                             [1 0 0 0 0 0]
                                [1 0 0 0 0]  [0 1 0 0 0 0]
                     [1 0 0 0]  [0 1 0 0 0]  [0 0 1 0 0 0]
            [1 0 0]  [0 1 0 0]  [0 0 1 0 0]  [0 0 0 1 0 0]
            [0 1 0]  [0 0 1 0]  [0 0 0 1 0]  [0 0 0 0 1 0]
            [0 0 1], [0 0 0 1], [0 0 0 0 1], [0 0 0 0 0 1]
            ]
            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%display ascii_art')   # indirect doctest
            sage: shell.run_cell("i = var('i')")
            sage: shell.run_cell('sum(i*x^i, i, 0, 10)')
                10      9      8      7      6      5      4      3      2
            10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x
            sage: shell.run_cell('%display default')
            sage: shell.quit()
        """
        return self.dm.displayhook(obj)
        


class SagePlainTextFormatter(PlainTextFormatter):

    def __init__(self, *args, **kwds):
        r"""
        Improved plain text formatter.

        .. warning::

            This formatter is usually NOT called. It is part of the
            IPython output formatting system, but Sage uses its
            own. Only when :meth:`SageDisplayFormatter.formatter`
            falls back to IPython's formatter will the
            :class:`SagePlainTextFormatter` be used.

        In particular, it correctly print lists of matrices or other
        objects (see
        :meth:`sage.structure.parent.Parent._repr_option`).

        EXAMPLES::
    
            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.display_formatter.formatters['text/plain']
            <sage.repl.display.formatter.SagePlainTextFormatter object at 0x...>
            sage: shell.run_cell('a = identity_matrix(ZZ, 2); [a,a]')
            [
            [1 0]  [1 0]
            [0 1], [0 1]
            ]
            sage: shell.quit()
        """
        super(SagePlainTextFormatter, self).__init__(*args, **kwds)

    @warn_format_error
    def __call__(self, obj):
        """
        Compute the pretty representation of the object.

        Adapted from ``IPython.core.formatters.PlainTextPrettyPrint``.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        String. The plain text representation.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: display_manager = get_display_manager()
            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: fmt = shell.display_formatter.formatters['text/plain']
            sage: fmt
            <sage.repl.display.formatter.SagePlainTextFormatter object at 0x...>
            sage: shell.displayhook.compute_format_data(2)
            ({u'text/plain': '2'}, {})
            sage: a = identity_matrix(ZZ, 2)
            sage: shell.displayhook.compute_format_data([a,a])
            ({u'text/plain': '[\n[1 0]  [1 0]\n[0 1], [0 1]\n]'}, {})
            sage: display_manager.preferences.text = 'ascii_art'
            sage: shell.displayhook.compute_format_data([a,a])
            ({u'text/plain': '[ [1 0]  [1 0] ]\n[ [0 1], [0 1] ]'}, {})
            sage: i = var('i')
            sage: shell.displayhook.compute_format_data(sum(i*x^i, i, 0, 10))
            ({u'text/plain': '    10      9      8      7      6      5      4      3
              2    \n10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x'},
             {})
            sage: del display_manager.preferences.text   # reset to default
            sage: shell.quit()
        """
        print('---- calling ipython formatter ----')
        import StringIO
        stream = StringIO.StringIO()
        printer = SagePrettyPrinter(
            stream, self.max_width, unicode_to_str(self.newline))
        printer.pretty(obj)
        printer.flush()
        return stream.getvalue()

