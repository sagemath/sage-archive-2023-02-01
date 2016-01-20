"""
List of assigned names in GAP

EXAMPLES::

    sage: from sage.libs.gap.assigned_names import KEYWORDS, GLOBALS, FUNCTIONS
    sage: 'fi' in KEYWORDS
    True
    sage: 'ZassenhausIntersection' in GLOBALS
    True
    sage: 'SubdirectProduct' in FUNCTIONS
    True
"""

###############################################################################
#       Copyright (C) 2016, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################


import cPickle
import string
from sage.libs.gap.libgap import libgap
from sage.libs.gap.saved_workspace import workspace


def load_or_compute(name, function):
    """
    Helper to load a cached value or compute it

    INPUT:

    - ``name`` -- string. Part of the cache filename

    - ``function`` -- function. To compute the value if not cached.

    OUTPUT:

    The value of ``function``, possibly cached.

    EXAMPLES::

        sage: from sage.libs.gap.assigned_names import GLOBALS
        sage: len(GLOBALS) > 1000    # indirect doctest
        True
        sage: from sage.libs.gap.saved_workspace import workspace
        sage: workspace(name='globals')
        ('...', True)
    """
    filename, up_to_date = workspace(name=name)
    if up_to_date:
        with open(filename, 'rb') as f:
            return cPickle.load(f)
    else:
        value = function()
        from sage.misc.temporary_file import atomic_write
        with atomic_write(filename) as f:
            cPickle.dump(value, f)
        return value


def list_keywords():
    """
    Return the GAP reserved keywords

    OUTPUT:

    Tuple of strings.

    EXAMPLES:

        sage: from sage.libs.gap.assigned_names import KEYWORDS
        sage: 'fi' in KEYWORDS
        True
    """
    keywords = libgap.get_global('GAPInfo')['Keywords'].sage()
    return tuple(sorted(keywords))


KEYWORDS = list_keywords()


def list_globals():
    """
    Return the GAP reserved keywords

    OUTPUT:

    Tuple of strings.

    EXAMPLES:

        sage: from sage.libs.gap.assigned_names import GLOBALS
        sage: 'ZassenhausIntersection' in GLOBALS
        True
    """
    NamesGVars = libgap.function_factory('NamesGVars')
    IsBoundGlobal = libgap.function_factory('IsBoundGlobal')
    gvars = set(
        name.sage() for name in NamesGVars()
        if IsBoundGlobal(name)
    )
    gvars.difference_update(KEYWORDS)
    return tuple(sorted(gvars))


GLOBALS = load_or_compute('globals', list_globals)


def list_functions():
    """
    Return the GAP documented global functions

    OUTPUT:

    Tuple of strings.

    EXAMPLES:

        sage: from sage.libs.gap.assigned_names import FUNCTIONS
        sage: 'IsBound' in FUNCTIONS    # is a keyword
        False
        sage: 'SubdirectProduct' in FUNCTIONS
        True
    """
    Filtered =libgap.function_factory('Filtered')
    IsDocumentedWord = libgap.function_factory('IsDocumentedWord')
    fnames = libgap.eval('GLOBAL_FUNCTION_NAMES')
    documented = Filtered(fnames, IsDocumentedWord)
    return tuple(sorted(documented.sage()))


FUNCTIONS = load_or_compute('functions', list_functions)

    




