"""
The Sage ZMQ Kernel

Version of the IPython kernel when running Sage inside the IPython
notebook or remote IPython sessions.
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sys
from IPython.kernel.zmq.ipkernel import IPythonKernel
from IPython.kernel.zmq.zmqshell import ZMQInteractiveShell
from IPython.utils.traitlets import Type

from sage.env import SAGE_VERSION, SAGE_EXTCODE, SAGE_DOC
from sage.repl.interpreter import SageNotebookInteractiveShell
from sage.repl.ipython_extension import SageCustomizations

class SageZMQInteractiveShell(SageNotebookInteractiveShell, ZMQInteractiveShell):
    pass


class SageKernel(IPythonKernel):    
    implementation = 'sage'
    implementation_version = SAGE_VERSION

    shell_class = Type(SageZMQInteractiveShell)

    def __init__(self, **kwds):
        super(SageKernel, self).__init__(**kwds)
        SageCustomizations(self.shell)

    @property
    def banner(self):
        from sage.misc.banner import banner_text
        return banner_text()

    help_links = [
        {
            'text': 'Sage Tutorial',
            'url': 'http://www.sagemath.org/doc/tutorial/index.html',
        },
        {
            'text': 'Thematic Tutorials',
            'url': 'http://www.sagemath.org/doc/thematic_tutorials/index.html',
        },
        {
            'text': 'FAQs',
            'url': 'http://www.sagemath.org/doc/faq/index.html',
        },
        {
            'text': 'PREP Tutorials',
            'url': 'http://www.sagemath.org/doc/prep/index.html',
        },
        {
            'text': 'Sage Reference',
            'url': 'http://www.sagemath.org/doc/reference/index.html',
        },
        {
            'text': 'Developers Guide',
            'url': 'http://www.sagemath.org/doc/developer/index.html',
        },
        {
            'text': "Python",
            'url': "http://docs.python.org/%i.%i" % sys.version_info[:2],
        },
        {
            'text': "IPython",
            'url': "http://ipython.org/documentation.html",
        },
        {
            'text': 'Singular',
            'url': 'http://www.singular.uni-kl.de/Manual/latest/index.htm',
        },
        {
            'text': 'GAP',
            'url': 'http://gap-system.org/Manuals/doc/ref/chap0.html',
        },
        {
            'text': "NumPy",
            'url': "http://docs.scipy.org/doc/numpy/reference/",
        },
        {
            'text': "SciPy",
            'url': "http://docs.scipy.org/doc/scipy/reference/",
        },
        {
            'text': "SymPy",
            'url': 'http://docs.sympy.org/latest/index.html',
        },
        {
            'text': "Matplotlib",
            'url': "http://matplotlib.org/contents.html",
        },
        {
            'text': "Markdown",
            'url': "http://help.github.com/articles/github-flavored-markdown",
        },
    ]
