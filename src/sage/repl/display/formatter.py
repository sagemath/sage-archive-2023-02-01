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
    sage: shell.run_cell('%display simple')

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


from IPython.core.formatters import PlainTextFormatter, warn_format_error
from IPython.utils.py3compat import str_to_unicode, unicode_to_str

from sage.repl.display.pretty_print import (
    SagePrettyPrinter, AsciiArtPrettyPrinter, TypesetPrettyPrinter
)


class SagePlainTextFormatter(PlainTextFormatter):

    def __init__(self, *args, **kwds):
        r"""
        Improved plain text formatter.

        In particular, it has the following two features:
    
        - correctly print lists of matrices or other objects (see
          :meth:`sage.structure.parent.Parent._repr_option`),

        - print ASCII art objects (like expressions) (see
          :meth:`sage.structure.sage_object.SageObject._ascii_art_`).
    
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
        """
        super(SagePlainTextFormatter, self).__init__(*args, **kwds)
        self.set_display('simple')

    def set_display(self, mode):
        r"""
        Select the text formatting method.

        INPUT:

        - ``mode`` -- string. One of ``simple``, ``ascii_art``, or ``typeset``.

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
            sage: shell.run_cell('%display simple')
        """
        if mode not in ['simple', 'ascii_art', 'typeset']:
            raise ValueError('invalid mode set')
        self._mode = mode
        self._pretty_printer_class = dict(
            simple=SagePrettyPrinter, 
            ascii_art=AsciiArtPrettyPrinter,
            typeset=TypesetPrettyPrinter,
        )[mode]

    _mode = 'simple'

    @property
    def simple(self):
        """
        Whether the mode is the "simple" (default) display.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: sys.displayhook.formatter.simple
            True
            sage: sys.displayhook.formatter.ascii_art
            False
        """
        return self._mode == 'simple'

    @property
    def ascii_art(self):
        """
        Whether the mode is the ascii art display.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: sys.displayhook.formatter.simple
            True
            sage: sys.displayhook.formatter.ascii_art
            False
        """
        return self._mode == 'ascii_art'

    @property
    def typeset(self):
        """
        Whether the mode is the notebook "Typeset" display.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: sys.displayhook.formatter.simple
            True
            sage: sys.displayhook.formatter.typeset
            False
        """
        return self._mode == 'typeset'

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
            sage: fmt.set_display('ascii_art')
            sage: shell.displayhook.compute_format_data([a,a])
            ({u'text/plain': '[ [1 0]  [1 0] ]\n[ [0 1], [0 1] ]'}, {})
            sage: i = var('i')
            sage: shell.displayhook.compute_format_data(sum(i*x^i, i, 0, 10))
            ({u'text/plain': '    10      9      8      7      6      5      4      3
              2    \n10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x'},
             {})
            sage: fmt.set_display('simple')
        """
        import StringIO
        stream = StringIO.StringIO()
        printer = self._pretty_printer_class(
            stream, self.max_width, unicode_to_str(self.newline))
        printer.pretty(obj)
        printer.flush()
        return stream.getvalue()


class SageDoctestTextFormatter(SagePlainTextFormatter):

    @warn_format_error
    def __call__(self, obj):
        """
        Display ``obj``.

        For doctests, we both 

        * Print the textual representation. This makes it clear in the
          documentation that the command returns graphics.

        * Run ``show()`` on graphics objects, to test that it
          correctly generates graphics. Note that, in
          ``DOCTEST_MODE``, the ``show()`` method will save graphics
          to temporary files but not launch a viewer.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        String. The plain text representation.

        EXAMPLES::

            sage: class FooGraphics(SageObject):
            ....:     def _graphics_(self, **kwds):
            ....:         print('showing graphics') 
            ....:         from sage.structure.graphics_file import GraphicsFile
            ....:         return GraphicsFile('/nonexistent.png', 'image/png')
            ....:     def _repr_(self):
            ....:         return 'Textual representation'
            sage: from sage.repl.display.formatter import SageDoctestTextFormatter
            sage: fmt = SageDoctestTextFormatter()
            sage: fmt(FooGraphics())
            showing graphics
            'Textual representation'
        """
        from sage.structure.sage_object import SageObject
        if isinstance(obj, SageObject) and hasattr(obj, '_graphics_'):
            obj._graphics_()      # ignore whether there actually is graphics
        return super(SageDoctestTextFormatter, self).__call__(obj)


class SageNBTextFormatter(SagePlainTextFormatter):

    @warn_format_error
    def __call__(self, obj):
        """
        Display ``obj``.

        .. warning::

            This is mostly a workaround for the old Notebook. Do not
            use it as a model for your own code.

        This is the default formatter for the old Sage Notebook
        (SageNB).

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        String. The plain text representation. Graphics output is a
        side effect.

        EXAMPLES::

            sage: class FooGraphics(SageObject):
            ....:     def _graphics_(self, **kwds):
            ....:         print('showing graphics') 
            ....:         from sage.structure.graphics_file import GraphicsFile
            ....:         return GraphicsFile('/nonexistent.png', 'image/png')
            ....:     def _repr_(self):
            ....:         return 'Textual representation'
            sage: from sage.repl.display.formatter import SageNBTextFormatter
            sage: fmt = SageNBTextFormatter()
            sage: fmt(FooGraphics())
            showing graphics
            ''
        """
        from sage.plot.plot3d.base import Graphics3d
        if isinstance(obj, Graphics3d):
            # TODO: Clean up Graphics3d and remove all the hardcoded SageNB stuff
            obj.show()
            return ''
        from sage.structure.sage_object import SageObject
        if isinstance(obj, SageObject) and hasattr(obj, '_graphics_'):
            gfx = obj._graphics_()
            if gfx:
                gfx.sagenb_embedding()
                return ''
        return super(SageNBTextFormatter, self).__call__(obj)


class SageConsoleTextFormatter(SagePlainTextFormatter):

    @warn_format_error
    def __call__(self, obj):
        """
        Display ``obj``.

        This is the default formatter for the Sage command line
        interface.

        If the object is graphics, the empty string is written to the
        output. An external viewer is started in a separate process as
        a side effect.

        Otherwise, the usual textual representation is generated.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        String. The plain text representation. Third-party viewers for
        media files can be launched as a side effect.

        EXAMPLES::

            sage: class FooGraphics(SageObject):
            ....:     def _graphics_(self, **kwds):
            ....:         print('showing graphics') 
            ....:         from sage.structure.graphics_file import GraphicsFile
            ....:         return GraphicsFile('/nonexistent.png', 'image/png')
            ....:     def _repr_(self):
            ....:         return 'Textual representation'
            sage: from sage.repl.display.formatter import SageConsoleTextFormatter
            sage: fmt = SageConsoleTextFormatter()
            sage: fmt(FooGraphics())
            showing graphics
            ''
        """
        from sage.structure.sage_object import SageObject
        if isinstance(obj, SageObject) and hasattr(obj, '_graphics_'):
            gfx = obj._graphics_()
            if gfx: 
                gfx.launch_viewer()
                return ''
        return super(SageConsoleTextFormatter, self).__call__(obj)
