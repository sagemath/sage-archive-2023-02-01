"""
The Sage ZMQ Kernel

Version of the IPython kernel when running Sage inside the IPython
notebook.
"""

from IPython.kernel.zmq.ipkernel import Kernel
from IPython.kernel.zmq.zmqshell import ZMQInteractiveShell
from IPython.utils.traitlets import Type
from IPython.core.formatters import DisplayFormatter

from sage.structure.graphics_file import Mime
from sage.repl.interpreter import SageInteractiveShell
from sage.repl.display.formatter import SagePlainTextFormatter
from sage.misc.temporary_file import tmp_filename
from sage.structure.sage_object import SageObject


class SageZMQDisplayFormatter(DisplayFormatter):

    def __init__(self, *args, **kwds):
        shell = kwds['parent']
        self.plain_text = SagePlainTextFormatter(config=shell.config)

    _format_types = frozenset([
        Mime.TEXT,
        Mime.HTML,
        Mime.LATEX,
        Mime.JSON,
        Mime.JAVASCRIPT,
        Mime.PDF,
        Mime.PNG,
        Mime.JPG,
        Mime.SVG,
    ])

    @property
    def format_types(self):
        """
        Return the enabled format types (MIME types)

        OUTPUT:

        Set of mime types (as strings).

        EXAMPLES::

            sage: from sage.repl.zmq_kernel import SageZMQDisplayFormatter
            sage: from sage.repl.interpreter import get_test_shell
            sage: fmt = SageZMQDisplayFormatter(parent=get_test_shell())
            sage: fmt.format_types
            frozenset({u'application/javascript',
                       u'application/json',
                       u'application/pdf',
                       u'image/jpeg',
                       u'image/png',
                       u'image/svg+xml',
                       u'text/html',
                       u'text/latex',
                       u'text/plain'})
        """
        return self._format_types

    # TODO: setter for format_types

    def format(self, obj, include=None, exclude=None):
        """
        Return a format data dict for an object
        """
        output = dict()
        if isinstance(obj, SageObject) and hasattr(obj, '_graphics_'):
            gfx = obj._graphics_(mime_types=self.format_types)
            if gfx is not None: 
                output[gfx.mime()] = gfx.data()
        if Mime.TEXT not in output:
            output[Mime.TEXT] = self.plain_text(obj)
        return (output, {})


class SageZMQInteractiveShell(SageInteractiveShell, ZMQInteractiveShell):

    def init_display_formatter(self):
        self.display_formatter = SageZMQDisplayFormatter(parent=self)
        self.configurables.append(self.display_formatter)


class SageKernel(Kernel):    
    shell_class = Type(SageZMQInteractiveShell)

