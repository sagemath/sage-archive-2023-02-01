"""
Support for Notebook Introspection and Setup

AUTHORS:

- William Stein (much of this code is from IPython).

- Nick Alexander
"""

EMBEDDED_MODE = False

from sage.misc.lazy_import import lazy_import
lazy_import('sagenb.misc.support',
        ('sage_globals', 'globals_at_init',
        'global_names_at_init', 'init', 'setup_systems', 'help',
        'get_rightmost_identifier', 'completions', 'docstring',
        'source_code', 'save_session', 'load_session', 'variables',
        'syseval'),
    deprecation=2891)


######################################################################
# Cython
######################################################################
import sage.misc.cython
import sys
import __builtin__

def cython_import(filename, verbose=False, compile_message=False,
                 use_cache=False, create_local_c_file=True, **kwds):
    """
    Compile a file containing Cython code, then import and return the
    module.  Raises an ``ImportError`` if anything goes wrong.

    INPUT:

    - ``filename`` - a string; name of a file that contains Cython
      code

    See the function :func:`sage.misc.cython.cython` for documentation
    for the other inputs.

    OUTPUT:

    - the module that contains the compiled Cython code.
    """
    name, build_dir = sage.misc.cython.cython(filename, verbose=verbose,
                                            compile_message=compile_message,
                                            use_cache=use_cache,
                                            create_local_c_file=create_local_c_file,
                                            **kwds)
    sys.path.append(build_dir)
    return __builtin__.__import__(name)


def cython_import_all(filename, globals, verbose=False, compile_message=False,
                     use_cache=False, create_local_c_file=True):
    """
    Imports all non-private (i.e., not beginning with an underscore)
    attributes of the specified Cython module into the given context.
    This is similar to::

        from module import *

    Raises an ``ImportError`` exception if anything goes wrong.

    INPUT:

    - ``filename`` - a string; name of a file that contains Cython
      code
    """
    m = cython_import(filename, verbose=verbose, compile_message=compile_message,
                     use_cache=use_cache,
                     create_local_c_file=create_local_c_file)
    for k, x in m.__dict__.iteritems():
        if k[0] != '_':
            globals[k] = x

