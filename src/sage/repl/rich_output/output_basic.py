# -*- encoding: utf-8 -*-
r"""
Basic Output Types

The Sage rich representation system requires a special container class
to hold the data for each type of rich output. They all inherit from
:class:`OutputBase`, though a more typical example is
:class:`OutputPlainText`. Some output classes consist of more than one
data buffer, for example jmol or certain animation formats. The output
class is independent of user preferences and of the display
backend. 

The display backends can form derive classes to attach
backend-specific display functionality to, for example how to launch a
viewer. But they must not change how the output container is
created. To enforce this, the Sage ``_rich_repr_`` magic method will
only ever see the output class defined here. The display manager will
promote it to a backend-specific subclass if necessary prior to
displaying it.

To create new types of output, you must create your own subclass of
:class:`OutputBase` and register it in the
:class:`~sage.repl.rich_output.display_manager.DisplayManager`.

.. warning::

    All rich output data in sublasses of :class:`OutputBase` must be
    contained in :class:`~sage.repl.rich_output.buffer.OutputBuffer`
    instances. You must never reference any files on the local file
    system, as there is no guarantee that the notebook server and the
    worker process are on the same computer. Or even share a common
    file system.  
"""

from sage.structure.sage_object import SageObject
from sage.repl.rich_output.buffer import OutputBuffer


class OutputBase(SageObject):
    """
    Base class for all rich output containers.
    """

    def _repr_(self):
        """
        Return a string representation.
        """
        return '{0} container'.format(self.__class__.__name__)

    @classmethod
    def example(cls):
        """
        Construct a sample instance

        This static method is meant for doctests, so they can easily
        construt an example.

        OUTPUT:

        An instance of the :class:`OutputBase` subclass.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputBase
            sage: OutputBase.example()
            Traceback (most recent call last):
            ...
            NotImplementedError: derived classes must implement this class method
        """
        raise NotImplementedError('derived classes must implement this class method')


class OutputPlainText(OutputBase):

    def __init__(self, plain_text):
        """
        Plain Text Output
        
        This should always be exactly the same as the (non-rich)
        output from the ``_repr_`` method. Every backend object must
        support plain text output as fallback.
        """
        self.text = OutputBuffer(plain_text)        

    @classmethod
    def example(cls):
        """
        Construct a sample plain text output container

        This static method is meant for doctests, so they can easily
        construt an example.

        OUTPUT:

        An instance of :class:`OutputPlainText`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputPlainText
            sage: OutputPlainText.example()
            OutputPlainText container
            sage: OutputPlainText.example().text.get()
            'Example plain text output'
        """
        return cls('Example plain text output')

    def print_to_stdout(self):
        """
        Write the data to stdout.
        """
        print(self.text.get())


class OutputAsciiArt(OutputBase):

    def __init__(self, ascii_art):
        """
        ASCII Art Output
        
        INPUT:

        - ``ascii_art`` -- ascii art rendered into a string.
        """
        self.ascii_art = OutputBuffer(ascii_art)        

    @classmethod
    def example(cls):
        r"""
        Construct a sample ascii art output container

        This static method is meant for doctests, so they can easily
        construt an example.

        OUTPUT:

        An instance of :class:`OutputAsciiArt`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputAsciiArt
            sage: OutputAsciiArt.example()
            OutputAsciiArt container
            sage: OutputAsciiArt.example().ascii_art.get()
            '[                        *   *   *    * ]\n[      **   **   *    *  *   *  *    *  ]\n[ ***, * , *  , **, ** , *, * , * , *   ]'
        """
        return cls('[                        *   *   *    * ]\n'
                   '[      **   **   *    *  *   *  *    *  ]\n'
                   '[ ***, * , *  , **, ** , *, * , * , *   ]')

    def print_to_stdout(self):
        """
        Write the data to stdout.
        """
        print(self.ascii_art.get())


class OutputMathJax(OutputBase):

    def __init__(self, math_tex):
        """
        MathJax Output
        
        INPUT:

        - ``math_tex`` -- string containing the math/tex
          code. Includes the surrounding ``<html>`` tag.
        """
        self.math_tex = OutputBuffer(math_tex)        

    @classmethod
    def example(cls):
        r"""
        Construct a sample MathJax output container

        This static method is meant for doctests, so they can easily
        construt an example.

        OUTPUT:

        An instance of :class:`OutputMathJax`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputMathJax
            sage: OutputMathJax.example()
            OutputMathJax container
            sage: OutputMathJax.example().math_tex.get()
            '<html><script type="math/tex; mode=display">\newcommand{\\Bold}[1]{\\mathbf{#1}}\\int \\sin\\left(x\right)\\,{d x}</script></html>'
        """
        return cls('<html><script type="math/tex; mode=display">'
                   '\newcommand{\Bold}[1]{\mathbf{#1}}'
                   '\int \sin\left(x\right)\,{d x}'
                   '</script></html>')

    def print_to_stdout(self):
        """
        Write the data to stdout.
        """
        print(self.math_tex.get())
