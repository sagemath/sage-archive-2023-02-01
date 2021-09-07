__all__ = ['all']

# Set sage.__version__ to the current version number. This is analogous
# to many other Python packages.
from sage.version import version as __version__

import sys
# Make sure that the correct zlib library is loaded. This is needed
# to prevent the system zlib to be loaded instead of the Sage one.
# See https://trac.sagemath.org/ticket/23122
import zlib

# IPython calls this when starting up
def load_ipython_extension(*args):
    import sage.repl.ipython_extension
    sage.repl.ipython_extension.load_ipython_extension(*args)


# Deprecated leftover of monkey-patching inspect.isfunction() to support Cython functions.
# We cannot use lazy_import for the deprecation here.
def isfunction(obj):
    """
    Check whether something is a function.

    This is a variant of ``inspect.isfunction``:
    We assume that anything which has a genuine ``__code__``
    attribute (not using ``__getattr__`` overrides) is a function.
    This is meant to support Cython functions.

    This function is deprecated.  Most uses of ``isfunction``
    can be replaced by ``callable``.

    EXAMPLES::

        sage: from sage import isfunction
        sage: def f(): pass
        sage: isfunction(f)
        doctest:warning...
        DeprecationWarning: sage.isfunction is deprecated; use callable or sage.misc.sageinspect.is_function_or_cython_function instead
        See https://trac.sagemath.org/32479 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(32479, "sage.isfunction is deprecated; use callable or sage.misc.sageinspect.is_function_or_cython_function instead")
    from sage.misc.sageinspect import is_function_or_cython_function
    return is_function_or_cython_function(obj)


# Monkey-patch ExtensionFileLoader to allow IPython to find the sources
# of Cython files. See https://trac.sagemath.org/ticket/24681
try:
    from importlib.machinery import ExtensionFileLoader
except ImportError:
    pass  # Python 2
else:
    del ExtensionFileLoader.get_source


# Work around a Cygwin-specific bug caused by sqlite3; see
# https://trac.sagemath.org/ticket/30157 and the docstring for
# fix_for_ticket_30157
# Here we monkey-patch the sqlite3 module to ensure the fix is
# applied the very first time a connection is made to a sqlite3
# database
if sys.platform == 'cygwin':
    def patch_sqlite3():
        try:
            from sage.misc.sage_ostools import fix_for_ticket_30157
        except ImportError:
            # The module might not have been re-built yet; don't worry about it
            # then
            return

        import sqlite3
        import functools
        orig_sqlite3_connect = sqlite3.connect

        @functools.wraps(orig_sqlite3_connect)
        def connect(*args, **kwargs):
            if fix_for_ticket_30157():
                raise RuntimeError(
                    'patch for Trac ticket #30157 failed; please report this '
                    'bug to https://trac.sagemath.org')

            # Undo the monkey-patch
            try:
                return orig_sqlite3_connect(*args, **kwargs)
            finally:
                sqlite3.connect = orig_sqlite3_connect

        sqlite3.connect = connect

    patch_sqlite3()
    del patch_sqlite3
