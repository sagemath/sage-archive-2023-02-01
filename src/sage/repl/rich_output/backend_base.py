# -*- encoding: utf-8 -*-
r"""
Base class for Backends

The display backends are the commandline, the SageNB notebook, the
ipython notebook, the Emacs sage mode, the Sage doctester, .... All of
these have different capabilities for what they can display.

To implement a new display backend, you need to subclass
:class:`BackendBase`. All backend-specific handlig of rich output
should be in :meth:`~BackendBase.displayhook` and
:meth:`~BackendBase.display_immediately`.

"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject


class BackendBase(SageObject):

    def _repr_(self):
        """
        Return string representation of the backend

        Every concrete backend must implement this method.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend._repr_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def get_display_manager(self):
        """
        Return the display manager singleton

        This is a convenience method to access the display manager
        singleton.

        OUTPUT:

        The unique
        :class:`~sage.repl.rich_output.display_manager.DisplayManager`
        instance.

        EXAMPLES::
        
            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.get_display_manager()
            The Sage display manager using the doctest backend
        """
        from sage.repl.rich_output import get_display_manager
        return get_display_manager()

    def install(self, **kwds):
        """
        Hook that will be called once before the backend is used for the
        first time.

        INPUT:

        - ``kwds`` -- optional keyword arguments that are passed
          through by the
          :meth:`~sage.repl.rich_output.display_manager.DisplayManager.switch_backend`
          method.
        """
        pass

    def uninstall(self):
        """
        Hook that will be called once before the backend is removed.
        """
        pass

    def default_preferences(self):
        """
        Return the backend's display preferences

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.preferences.DisplayPreferences`.
        """
        from sage.repl.rich_output.preferences import DisplayPreferences
        return DisplayPreferences()

    def supported_output(self):
        """
        Return the outputs that are supported by the backend.

        OUTPUT:

        Iterable of output container classes, that is, subclass of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`).
        The order is ignored.
        """
        raise NotImplementedError

    def max_width(self):
        """
        Return the number of characters that fit into one output line

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.max_width()
            79
        """
        return 79

    def newline(self):
        r"""
        Return the newline string.

        OUTPUT:

        String for starting a new line of output.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.newline()
            '\n'
        """
        return '\n'

    def _apply_pretty_printer(self, pretty_printer_class, obj):
        """
        Helper method to format ``obj`` as text

        INPUT:

        - ``pretty_printer_class`` -- subclass of
          :class:`sage.repl.display.pretty_print.SagePrettyPrinter`.

        - ``obj`` -- anything.

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: from sage.repl.display.pretty_print import SagePrettyPrinter
            sage: backend._apply_pretty_printer(SagePrettyPrinter, 1/2)
            '1/2'
        """
        import StringIO
        stream = StringIO.StringIO()
        printer = pretty_printer_class(
            stream, self.max_width(), self.newline())
        printer.pretty(obj)
        printer.flush()
        return stream.getvalue()

    def plain_text_formatter(self, obj):
        r"""
        Hook to override how plain text is being formatted.

        If the object does not have a ``_rich_repr_`` method, or if it
        does not return a rich output object
        (:class:`~sage.repl.rich_output.output_basic.OutputBase`),
        then this method is used to generate plain text output.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.output_basic.OutputPlainText`
        containing the string representation of the object.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: out = backend.plain_text_formatter(range(30))
            sage: out
            OutputPlainText container
            sage: out.text
            buffer containing 139 bytes
            sage: out.text.get()
            '[0,\n 1,\n 2,\n 3,\n 4,\n 5,\n 6,\n 7,\n 8,\n 9,\n 
            10,\n 11,\n 12,\n 13,\n 14,\n 15,\n 16,\n 17,\n 18,\n 
            19,\n 20,\n 21,\n 22,\n 23,\n 24,\n 25,\n 26,\n 27,\n 
            28,\n 29]'
        """
        from sage.repl.display.pretty_print import SagePrettyPrinter
        plain_text = self._apply_pretty_printer(SagePrettyPrinter, obj)
        from sage.repl.rich_output.output_basic import OutputPlainText
        return OutputPlainText(plain_text)

    def ascii_art_formatter(self, obj):
        """
        Hook to override how ascii art is being formatted.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.output_basic.OutputAsciiArt`
        containing the ascii art string representation of the object.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: out = backend.ascii_art_formatter(range(30))
            sage: out
            OutputAsciiArt container
            sage: out.ascii_art
            buffer containing 228 bytes
            sage: print(out.ascii_art.get())
            [                                                                              
            [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
            <BLANKLINE>
                                            ]
             22, 23, 24, 25, 26, 27, 28, 29 ]
        """
        from sage.repl.display.pretty_print import AsciiArtPrettyPrinter
        ascii_art = self._apply_pretty_printer(AsciiArtPrettyPrinter, obj)
        from sage.repl.rich_output.output_basic import OutputAsciiArt
        return OutputAsciiArt(ascii_art)

    def mathjax_formatter(self, obj):
        r"""
        Hook to override how MathJax (math/tex) is being formatted.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.output_basic.OutputMathJax`
        containing the math/tex string representation of the object.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: out = backend.mathjax_formatter(1/2)
            sage: out
            OutputMathJax container
            sage: out.math_tex
            buffer containing 91 bytes
            sage: out.math_tex.get()
            '<html><script type="math/tex">\\newcommand{\\Bold}[1]{\\mathbf{#1}}\\frac{1}{2}</script></html>'
        """
        from sage.misc.latex import MathJax
        mathjax = MathJax().eval(obj, mode='inline', combine_all=True)
        from sage.repl.rich_output.output_basic import OutputMathJax
        return OutputMathJax(str(mathjax))

    def set_underscore_variable(self, obj):
        """
        Set the ``_`` builtin variable.

        INPUT:

        - ``obj`` -- result of the most recent evaluation.

        Typically this sets the special ``_`` variable. Backends that
        organize the history differently (e.g. IPython) can override
        this method.

        EXAMPLES::

            sage: 'foo'
            'foo'
            sage: _
            'foo'
        """
        import __builtin__
        __builtin__._ = obj
    
    def displayhook(self, plain_text, rich_output):
        """
        Backend implementation of the displayhook
        
        The value of the last statement on a REPL input line or
        notebook cell are usually handed to the Python displayhook and
        shown on screen.  By overriding this method you define how
        your backend handles output. The difference to the usual
        displayhook is that Sage already converted the value to the
        most suitable rich output container.

        INPUT:

        - ``plain_text`` -- instance of
          :class:`~sage.repl.rich_output.output_basic.OutputPlainText`. The
          plain text version of the output.

        - ``rich_output`` -- instance of an output container class
          (subclass of
          :class:`~sage.repl.rich_output.output_basic.OutputBase`). Guaranteed
          to be one of the output containers returned from
          :meth:`supported_output`, possibly the same as
          ``plain_text``.

        OUTPUT:

        This method may return something, which is then returned from
        the display manager's
        :meth:`~sage.repl.rich_output.display_manager.DisplayManager.displayhook`
        method.
        """
        return self.display_immediately(plain_text, rich_output)

    def display_immediately(self, plain_text, rich_output):
        """
        Show output without going back to the command line prompt.

        This method is similar to the rich output :meth:`displayhook`,
        except that it can be invoked at any time. Typically, it ends
        up being called by :meth:`sage.plot.graphics.Graphics.show`.
        
        INPUT:

        Same as :meth:`displayhook`.
        """
        raise NotImplementedError


class BackendSimple(BackendBase):
    """
    Simple Backend

    This backend only supports plain text.
    
    EXAMPLES::

        sage: from sage.repl.rich_output.backend_base import BackendSimple
        sage: BackendSimple()
        simple
    """
    
    def _repr_(self):
        r"""
        Return string representation of the backend

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendSimple
            sage: backend = BackendSimple()
            sage: backend._repr_()
            'simple'
        """
        return 'simple'

    def supported_output(self):
        from sage.repl.rich_output.output_basic import OutputPlainText
        return set([OutputPlainText])

    def display_immediately(self, plain_text, rich_output):
        print(rich_output.text.get())

