r"""
Lazy imports

This module allows one to lazily import callable objects into the
global namespace, where the actual import is delayed until the object
is actually called or inspected. This is useful for modules that are
expensive to import or may cause circular references, though there is
some overhead in its use.

EXAMPLES::

    sage: from sage.misc.lazy_import import lazy_import
    sage: lazy_import('sage.rings.all', 'ZZ')
    sage: type(ZZ)
    <class 'sage.misc.lazy_import.LazyImport'>
    sage: ZZ(4.0)
    4

AUTHOR:

 - Robert Bradshaw
"""

#*****************************************************************************
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import inspect
import sageinspect

class LazyImport(object):
    """
    EXAMPLES::

        sage: from sage.misc.lazy_import import LazyImport
        sage: my_integer = LazyImport('sage.rings.all', 'Integer')
        sage: my_integer(4)
        4
        sage: my_integer('101', base=2)
        5
        sage: my_integer(3/2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    def __init__(self, module, name):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: my_isprime(5)
            True
            sage: my_isprime(55)
            False
        """
        self._module = module
        self._name = name
        self._object = None

    def _get_object(self):
        """
        Return the wrapped object, importing it if necessary.

        EXAMPLE::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_integer_ring = LazyImport('sage.rings.all', 'ZZ')
            sage: my_integer_ring._object is None
            True
            sage: my_integer_ring._get_object()
            Integer Ring
            sage: my_integer_ring._object is None
            False
        """
        if self._object is None:
            self._object = getattr(__import__(self._module, {}, {}, [self._name]), self._name)
        return self._object

    def _sage_doc_(self):
        """
        Return the docstring of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: my_isprime._sage_doc_() is is_prime.__doc__
            True
        """
        return sageinspect._sage_getdoc_unformatted(self._get_object())

    def _sage_src_(self):
        """
        Returns the source of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: 'def is_prime(' in my_isprime._sage_src_()
            True
        """
        return sageinspect.sage_getsource(self._get_object())

    def _sage_argspec_(self):
        """
        Returns the argspec of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: rm = LazyImport('sage.all', 'random_matrix')
            sage: rm._sage_argspec_()
            (['ring', 'nrows', 'ncols', 'algorithm'], 'args', 'kwds', (None, 'randomize'))
        """
        return sageinspect.sage_getargspec(self._get_object())

    def __getattr__(self, attr):
        """
        Attribute lookup on self defers to attribute lookup on the
        wrapped object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_integer = LazyImport('sage.rings.all', 'Integer')
            sage: my_integer.sqrt is Integer.sqrt
            True
        """
        return getattr(self._get_object(), attr)

    def __dir__(self):
        """
        Tab completion on self defers to completion on the wrapped
        object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: dir(my_ZZ) == dir(ZZ)
            True
        """
        return dir(self._get_object())

    def __call__(self, *args, **kwds):
        """
        Calling self calls the wrapped object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: my_isprime(12)
            False
            sage: my_isprime(13)
            True
        """
        return self._get_object()(*args, **kwds)

def lazy_import(module, names, _as=None):
    """
    Create a lazy import object and inject it into the caller's global
    namespace. For the purposes of introspection and calling, this
    like performing a lazy "from module import name" where the import
    is delayed until the object actually is used or inspected.

    EXAMPLES::

        sage: from sage.misc.lazy_import import lazy_import
        sage: lazy_import('sage.rings.all', 'ZZ')
        sage: type(ZZ)
        <class 'sage.misc.lazy_import.LazyImport'>
        sage: ZZ(4.0)
        4
        sage: lazy_import('sage.rings.all', 'RDF', 'my_RDF')
        sage: my_RDF(1/2)
        0.5
        sage: my_RDF._get_object() is RDF
        True

        sage: lazy_import('sage.all', ['QQ', 'RR'], ['my_QQ', 'my_RR'])
        sage: my_QQ._get_object() is QQ
        True
        sage: my_RR._get_object() is RR
        True
    """
    if _as is None:
        _as = names
    if isinstance(names, str):
        names = [names]
        _as = [_as]
    calling_globals = inspect.currentframe().f_back.f_globals
    for name, alias in zip(names, _as):
        if alias is None:
            alias = name
        calling_globals[alias] = LazyImport(module, name)
