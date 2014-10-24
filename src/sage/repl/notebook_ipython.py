"""
Sage-Enhanced IPython Notebook

.. note::

    The customized Jinja2 templates for the IPython notebook are in
    ``SAGE_LOCAL/share/sage/ext/ipython-notebook``. You may have to
    update them as well when updating to a new IPython version.
"""

import os
import copy

from IPython.html.notebookapp import NotebookApp
from IPython import Config

from sage.env import DOT_SAGE, SAGE_EXTCODE, SAGE_DOC
from sage.repl.interpreter import (
    SageCrashHandler, DEFAULT_SAGE_CONFIG, SageInteractiveShell,
)


# The notebook Jinja2 templates and static files
TEMPLATE_PATH = os.path.join(SAGE_EXTCODE, 'notebook-ipython', 'templates')
STATIC_PATH = os.path.join(SAGE_EXTCODE, 'notebook-ipython', 'static')
DOCS_PATH = os.path.join(SAGE_DOC, 'output', 'html', 'en')


# Note: sage.repl.interpreter.DEFAULT_SAGE_CONFIG will be applied, too
DEFAULT_SAGE_NOTEBOOK_CONFIG = Config(
    SageNotebookApp = Config(
        # log_level = 'DEBUG',       # if you want more logs
        # open_browser = False,      # if you want to avoid browser restart
        webapp_settings = Config(
            template_path = TEMPLATE_PATH,
        ),
        extra_static_paths = [STATIC_PATH, DOCS_PATH],
    ),
)


class SageNotebookApp(NotebookApp):
    name = u'sage-notebook-ipython'
    crash_handler_class = SageCrashHandler

    def load_config_file(self, *args, **kwds):
        r"""
        Merges a config file with the default sage notebook config.

        EXAMPLES::

            sage: from sage.misc.temporary_file import tmp_dir
            sage: from sage.repl.notebook_ipython import SageNotebookApp
            sage: d = tmp_dir()
            sage: IPYTHONDIR = os.environ['IPYTHONDIR']
            sage: os.environ['IPYTHONDIR'] = d
            sage: app = SageNotebookApp()
            sage: app.load_config_file()    # random output
            2014-09-16 23:57:35.6 [SageNotebookApp] Created profile dir: 
            u'/home/vbraun/.sage/temp/desktop.localdomain/1490/dir_ZQupP5/profile_default'
            sage: app.notebook_dir          # random output
            u'/home/vbraun/'
            sage: os.environ['IPYTHONDIR'] = IPYTHONDIR
        """
        super(SageNotebookApp, self).load_config_file(*args, **kwds)
        newconfig = copy.deepcopy(DEFAULT_SAGE_CONFIG)
        newconfig.merge(DEFAULT_SAGE_NOTEBOOK_CONFIG)
        newconfig.merge(self.config)
        self.config = newconfig

    def init_kernel_argv(self):
        """
        Construct the kernel arguments
        
        The kernel is running in a separate process, so it does not
        see our config dictionary. Any options need to be passed
        through the command line interface.

        EXAMPLES::

            sage: from sage.repl.notebook_ipython import SageNotebookApp
            sage: d = tmp_dir()
            sage: IPYTHONDIR = os.environ['IPYTHONDIR']
            sage: os.environ['IPYTHONDIR'] = d
            sage: app = SageNotebookApp()
            sage: app.kernel_argv
            []
            sage: app.init_kernel_argv()    # random output
            2014-09-16 23:57:35.6 [SageNotebookApp] Created profile dir: 
            u'/home/vbraun/.sage/temp/desktop.localdomain/1490/dir_ZQupP5/profile_default'
            sage: app.kernel_argv
            [u"--IPKernelApp.parent_appname='sage-notebook-ipython'",
             '--profile-dir',
             u'/.../profile_default',
             u'--IPKernelApp.kernel_class=sage.repl.zmq_kernel.SageKernel',
             u'--IPKernelApp.extra_extension=sage']
            sage: os.environ['IPYTHONDIR'] = IPYTHONDIR
        """
        super(SageNotebookApp, self).init_kernel_argv()
        self.kernel_argv.append(
            u'--IPKernelApp.kernel_class=sage.repl.zmq_kernel.SageKernel')
        self.kernel_argv.append(
            u'--IPKernelApp.extra_extension=sage')
