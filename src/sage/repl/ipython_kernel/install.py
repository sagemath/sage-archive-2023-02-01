"""
Sage-Enhanced IPython Notebook

.. note::

    The customized Jinja2 templates for the IPython notebook are in
    ``SAGE_LOCAL/share/sage/ext/ipython-notebook``. You may have to
    update them as well when updating to a new IPython version.
"""

import os
import errno
import copy

from IPython.kernel.kernelspec import (
    get_kernel_spec, install_kernel_spec, NoSuchKernel)
    

from sage.env import (
    SAGE_ROOT, SAGE_DOC, SAGE_LOCAL, SAGE_EXTCODE,
    SAGE_VERSION
)
from sage.misc.temporary_file import tmp_dir


class SageKernelSpec(object):

    def __init__(self):
        """
        Utility to manage Sage kernels

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: spec._display_name    # not tested
            'Sage 6.6.beta2'
            sage: spec._identifier      # not tested
            'sage_6_6_beta2'
        """
        self._display_name = 'Sage {0}'.format(SAGE_VERSION)

    @classmethod
    def identifier(self):
        return 'Sage {0}'.format(SAGE_VERSION).lower().replace(' ', '_').replace('.', '_')
        
    def symlink(self, src, dst):
        """
        Symlink ``src`` to ``dst``

        This is not an atomic operation.

        Already-existing symlinks will be deleted, already existing
        non-empty directories will be kept.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.install import SageKernelSpec
            sage: spec = SageKernelSpec()
            sage: path = tmp_dir()
            sage: spec.symlink(os.path.join(path, 'a'), os.path.join(path, 'b'))
            sage: os.listdir(path)
            ['b']
        """
        try:
            os.remove(dst)
        except OSError as err:
            if err.errno == errno.EEXIST:
                return
        os.symlink(src, dst)
        
    def use_local_mathjax(self):
        ipython_dir = os.environ['IPYTHONDIR']
        src = os.path.join(SAGE_LOCAL, 'share', 'mathjax')
        dst = os.path.join(ipython_dir, 'nbextensions', 'mathjax')
        self.symlink(src, dst)

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
    
    def install_spec(self):
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
        identifier = self.identifier()
        install_kernel_spec(temp, identifier, user=True, replace=True)
        self._spec = get_kernel_spec(identifier)

    def symlink_resources(self):
        spec_dir = self._spec.resource_dir
        path = os.path.join(SAGE_EXTCODE, 'notebook-ipython')
        for filename in os.listdir(path):
            self.symlink(
                os.path.join(path, filename),
                os.path.join(spec_dir, filename)
            )
        self.symlink(
            os.path.join(SAGE_DOC, 'output', 'html', 'en'),
            os.path.join(spec_dir, 'doc')
        )
        
    @classmethod
    def update(cls):
        instance = cls()
        instance.use_local_mathjax()
        instance.install_spec()
        instance.symlink_resources()
        

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




