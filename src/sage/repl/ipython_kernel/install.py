"""
Sage-Enhanced IPython Notebook

.. note::

    The customized Jinja2 templates for the IPython notebook are in
    ``SAGE_LOCAL/share/sage/ext/ipython-notebook``. You may have to
    update them as well when updating to a new IPython version.
"""

import os
import copy

from IPython.kernel.kernelspec import (
    get_kernel_spec, install_kernel_spec, NoSuchKernel)
    

from sage.env import (
    SAGE_ROOT, DOT_SAGE, SAGE_LOCAL,
    SAGE_VERSION
)
from sage.misc.temporary_file import tmp_dir



# # The notebook Jinja2 templates and static files
# TEMPLATE_PATH = os.path.join(SAGE_EXTCODE, 'notebook-ipython', 'templates')
# STATIC_PATH = os.path.join(SAGE_EXTCODE, 'notebook-ipython', 'static')
# DOCS_PATH = os.path.join(SAGE_DOC, 'output', 'html', 'en')


# # Note: sage.repl.interpreter.DEFAULT_SAGE_CONFIG will be applied, too
# DEFAULT_SAGE_NOTEBOOK_CONFIG = Config(
#     SageNotebookApp = Config(
#         # log_level = 'DEBUG',       # if you want more logs
#         # open_browser = False,      # if you want to avoid browser restart
#         webapp_settings = Config(
#             template_path = TEMPLATE_PATH,
#         ),
#         extra_static_paths = [STATIC_PATH, DOCS_PATH],
#     ),
# )


class SageKernelSpec(object):

    def __init__(self):
        """
        Utility to manage Sage kernels

        EXAMPLES::

            sage: from sage.repl.notebook_ipython import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec._display_name
            sage: spec._identifier
        """
        self._display_name = 'Sage {0}'.format(SAGE_VERSION)
        self._identifier = self._display_name.lower().replace(' ', '_').replace('.', '_')

    def use_local_mathjax(self):
        ipython_dir = os.environ['IPYTHONDIR']
        src = os.path.join(SAGE_LOCAL, 'share', 'mathjax')
        dst = os.path.join(ipython_dir, 'nbextensions', 'mathjax')
        print('symlink', src, dst)
        if not os.path.exists(dst):
            os.symlink(src, dst)
        
    def _kernel_cmd(self):
        return [
            os.path.join(SAGE_ROOT, 'sage'),
            '-python',
            '-m', 'sage.repl.ipython_kernel',
            '-f', '{connection_file}',
        ]
        
    def kernel_spec(self):
        """
        Return the kernel spec as Python dictionary
        """
        return dict(
            argv=self._kernel_cmd(),
            display_name=self._display_name,
        )
    
    def install(self):
        """
        Install the Sage IPython kernel
        
        It is safe to call this method multiple times, only one Sage
        kernel spec is ever installed.
        """
        import json
        temp = tmp_dir()
        kernel_spec = os.path.join(temp, 'kernel.json')
        with open(kernel_spec, 'w') as f:
            json.dump(self.kernel_spec(), f)
        install_kernel_spec(temp, self._identifier, user=True, replace=True)

    @classmethod
    def update(cls):
        instance = cls()
        instance.use_local_mathjax()
        instance.install()

        

def have_prerequisites():
    """
    Check that we have all prerequisites to run the IPython notebook.

    In particular, the IPython notebook requires OpenSSL whether or
    not you are using https. See trac:`17318`.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: from sage.repl.notebook_ipython import have_prerequisites
        sage: have_prerequisites() in [True, False]
        True
    """
    try:
        from IPython.html.notebookapp import NotebookApp
        return True
    except ImportError as err:
        return False




