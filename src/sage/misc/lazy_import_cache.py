"""
This is a pure Python file with no dependancies so it can be used in setup.py.
"""

import os

def get_cache_file():
    """
    Returns a per-branch file for caching names of lazily imported modules.

    EXAMPLES::

        sage: from sage.misc.lazy_import_cache import get_cache_file
        sage: get_cache_file()
        '...-lazy_import_cache.pickle'
        sage: get_cache_file().startswith(DOT_SAGE)
        True
        sage: from sage.misc.misc import branch_current_hg
        sage: branch_current_hg() in get_cache_file()
        True

    It shouldn't matter whether DOT_SAGE ends with a slash::

        sage: OLD = DOT_SAGE
        sage: os.environ['DOT_SAGE'] = '/tmp'
        sage: get_cache_file().startswith('/tmp/')
        True
        sage: os.environ['DOT_SAGE'] = OLD
    """
    mangled = os.path.realpath(os.path.join(os.environ['SAGE_ROOT'], 'devel', 'sage')).replace(os.sep, '_')
    return os.path.join(os.environ['DOT_SAGE'],
                        "%s-lazy_import_cache.pickle" % mangled)
