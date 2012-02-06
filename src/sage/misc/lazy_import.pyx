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
    <type 'sage.misc.lazy_import.LazyImport'>
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

cdef extern from *:
    cdef int Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE

import os, shutil, tempfile, cPickle as pickle, operator
import inspect
import sageinspect

from lazy_import_cache import get_cache_file

cdef binop(op, left, right):
    if isinstance(left, LazyImport):
        left = (<LazyImport>left)._get_object()
    if isinstance(right, LazyImport):
        right = (<LazyImport>right)._get_object()
    return op(left, right)

cdef class LazyImport(object):
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

    cdef readonly _object
    cdef _module
    cdef _name
    cdef _as_name
    cdef _namespace

    def __init__(self, module, name, as_name=None, namespace=None):
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
        self._as_name = as_name
        self._namespace = namespace

    cpdef _get_object(self, owner=None):
        """
        Return the wrapped object, importing it if necessary.

        INPUT:

        - ``owner`` -- ``None`` or the class (or subclass thereof)
          which contains this :class:`LazyImport` object in its
          ``__dict__``.

        OUTPUT:

        - the wrapped object

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_integer_ring = LazyImport('sage.rings.all', 'ZZ')
            sage: my_integer_ring._object is None
            True
            sage: my_integer_ring._get_object()
            Integer Ring
            sage: my_integer_ring._object is None
            False

        .. note::

           For a :class:`LazyImport` object that appears in a class
           namespace, we need to do something special. Indeed, the
           class namespace dictionary at the time of the class
           definition is not the one that actually gets used. Thus,
           when this function is called, :meth:`__get__`, ``owner``
           should be set to the ``owner`` class passed into
           ``__get__``::

               sage: class Foo(object):
               ...       lazy_import('sage.all', 'plot')
               sage: class Bar(Foo):
               ...       pass
               sage: type(Foo.__dict__['plot'])
               <type 'sage.misc.lazy_import.LazyImport'>

           Here is how :meth:`_get_object` is called internally upon
           ``Bar.plot``::

               sage: Foo.__dict__['plot']._get_object(Bar)
               <function plot at ...>

           Now ``Bar`` has been replaced in the dictionary of ``Foo``::

               sage: type(Foo.__dict__['plot'])
               <type 'function'>
        """
        if self._object is None:
            self._object = getattr(__import__(self._module, {}, {}, [self._name]), self._name)
            alias = self._as_name or self._name
            if owner is None:
                if self._namespace and self._namespace[alias] is self:
                    self._namespace[alias] = self._object
            else:
                from inspect import getmro
                for cls in getmro(owner):
                    if cls.__dict__.get(alias, None) is self:
                        setattr(cls, alias, self._object)
                        break
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
            ArgSpec(args=['ring', 'nrows', 'ncols', 'algorithm'], varargs='args', keywords='kwds', defaults=(None, 'randomize'))
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

    # We need to wrap all the slot methods, as they are not forwarded
    # via getattr.

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

    def __repr__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: type(lazy_ZZ)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: lazy_ZZ
            Integer Ring
            sage: repr(lazy_ZZ)
            'Integer Ring'
        """
        return repr(self._get_object())

    def __str__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: str(lazy_ZZ)
            'Integer Ring'
        """
        return str(self._get_object())

    def __unicode__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: unicode(lazy_ZZ)
            u'Integer Ring'
        """
        return unicode(self._get_object())

    def __nonzero__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: not lazy_ZZ
            True
        """
        return not self._get_object()

    def __hash__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: hash(lazy_ZZ) == hash(1.parent())
            True
        """
        return hash(self._get_object())

    def __cmp__(left, right):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: cmp(lazy_ZZ, ZZ)
            0
            sage: cmp(lazy_ZZ, QQ)
            -1
        """
        return binop(cmp, left, right)

    def __richcmp__(left, right, int op):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: lazy_ZZ == RR
            False
            sage: lazy_ZZ == 1.parent()
            True
        """
        if isinstance(left, LazyImport):
            left = (<LazyImport>left)._get_object()
        if isinstance(right, LazyImport):
            right = (<LazyImport>right)._get_object()
        if op == Py_LT:
            return left <  right
        elif op == Py_LE:
            return left <= right
        elif op == Py_EQ:
            return left == right
        elif op == Py_NE:
            return left != right
        elif op == Py_GT:
            return left >  right
        elif op == Py_GE:
            return left >= right

    def __len__(self):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: len(version_info)
            5
        """
        return len(self._get_object())

    def __get__(self, instance, owner):
        """
        EXAMPLES:

        Here we show how to take a function in a module, and lazy
        import it as a method of a class. For the sake of this
        example, we add manually a function in sage.all::

            sage: def my_method(self): return self
            sage: sage.all.my_method = my_method

        Now we lazy import it as a method of a new class ``Foo``::

            sage: from sage.misc.lazy_import import LazyImport
            sage: class Foo:
            ...       my_method = LazyImport('sage.all', 'my_method')

        Now we can use it as a usual method::

            sage: Foo().my_method()
            <__main__.Foo instance at ...>
            sage: Foo.my_method
            <unbound method Foo.my_method>
            sage: Foo().my_method
            <bound method Foo.my_method of <__main__.Foo instance at ...>>

        When a :class:`LazyImport` method is a method (or attribute)
        of a class, then extra work must be done to replace this
        :class:`LazyImport` object with the actual object. See the
        documentation of :meth:`_get_object` for an explanation of
        this.
        """
        obj = self._get_object(owner=owner)
        if hasattr(obj, "__get__"):
            return obj.__get__(instance, owner)
        return obj

    def __getitem__(self, key):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: version_info[0]
            2
        """
        return self._get_object()[key]

    def __setitem__(self, key, value):
        """
        TESTS::

            sage: sage.all.foo = range(10)
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo[1] = 100
            sage: print foo
            [0, 100, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        self._get_object()[key] = value

    def __delitem__(self, key):
        """
        TESTS::

            sage: sage.all.foo = range(10)
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: del foo[1]
            sage: print foo
            [0, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        del self._get_object()[key]

    def __iter__(self):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: iter(version_info)
            <iterator object at ...>
        """
        return iter(self._get_object())

    def __contains__(self, item):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: 2 in version_info
            True

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: 2000 not in version_info
            True
        """
        return item in self._get_object()

    def __add__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo + 1
            11
        """
        return binop(operator.add, left, right)

    def __sub__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo - 1
            9
        """
        return binop(operator.sub, left, right)

    def __mul__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo * 2
            20
        """
        return binop(operator.mul, left, right)

    def __div__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo / 2
            5
        """
        return binop(operator.div, left, right)

    def __floordiv__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo  // 3
            3
        """
        return binop(operator.floordiv, left, right)

    def __truediv__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: operator.truediv(foo, 3)
            10/3
        """
        return binop(operator.truediv, left, right)

    def __pow__(left, right, mod):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo ** 2
            100
        """
        if isinstance(left, LazyImport):
            left = (<LazyImport>left)._get_object()
        if isinstance(right, LazyImport):
            right = (<LazyImport>right)._get_object()
        if mod is None:
            return left ** right
        else:
            return left.__pow__(right, mod)

    def __mod__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo % 7
            3
        """
        return binop(operator.mod, left, right)

    def __lshift__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo << 3
            80
        """
        return binop(operator.lshift, left, right)

    def __rshift__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo >> 2
            2
        """
        return binop(operator.rshift, left, right)

    def __and__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo & 7
            2
        """
        return binop(operator.and_, left, right)

    def __or__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo | 7
            15
        """
        return binop(operator.or_, left, right)

    def __xor__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo ^^ 7
            13
        """
        return binop(operator.xor, left, right)

    def __neg__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: -foo
            -10
        """
        return -self._get_object()

    def __pos__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: +foo
            10
        """
        return +self._get_object()

    def __abs__(self):
        """
        TESTS::

            sage: sage.all.foo = -1000
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: abs(foo)
            1000
        """
        return abs(self._get_object())

    def __invert__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: ~foo
            1/10
        """
        return ~self._get_object()

    def __complex__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: complex(foo)
            (10+0j)
        """
        return complex(self._get_object())

    def __int__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: int(foo)
            10
        """
        return int(self._get_object())

    def __long__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: long(foo)
            10L
        """
        return long(self._get_object())

    def __float__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: float(foo)
            10.0
        """
        return float(self._get_object())

    def __oct__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: oct(foo)
            '12'
        """
        return oct(self._get_object())

    def __hex__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: hex(foo)
            'a'
        """
        return hex(self._get_object())

    def __index__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: range(100)[foo]
            10
        """
        return operator.index(self._get_object())


def lazy_import(module, names, _as=None, namespace=None, bint overwrite=True):
    """
    Create a lazy import object and inject it into the caller's global
    namespace. For the purposes of introspection and calling, this
    like performing a lazy "from module import name" where the import
    is delayed until the object actually is used or inspected.

    INPUT:

    - module: a string representing the module to import.
    - names: a string or list of strings representing the names to import from
      module.
    - _as: (optional) a string or list of strings representing the aliases of the names
      imported.
    - namespace: The namespace where importing the names. By default, import
      the names to current namespace.
    - overwrite: (default: True) If set to True and a name is already in the
      namespace, overwrite it with the lazy_import-ed name.

    EXAMPLES::

        sage: from sage.misc.lazy_import import lazy_import
        sage: lazy_import('sage.rings.all', 'ZZ')
        sage: type(ZZ)
        <type 'sage.misc.lazy_import.LazyImport'>
        sage: ZZ(4.0)
        4
        sage: lazy_import('sage.rings.all', 'RDF', 'my_RDF')
        sage: my_RDF._get_object() is RDF
        True
        sage: my_RDF(1/2)
        0.5

        sage: lazy_import('sage.all', ['QQ', 'RR'], ['my_QQ', 'my_RR'])
        sage: my_QQ._get_object() is QQ
        True
        sage: my_RR._get_object() is RR
        True

    Upon the first use, the object is injected directly into
    the calling namespace::

        sage: lazy_import('sage.all', 'ZZ', 'my_ZZ')
        sage: my_ZZ is ZZ
        False
        sage: my_ZZ(37)
        37
        sage: my_ZZ is ZZ
        True

    We check that :func:`lazy_import` also works for methods::

        sage: class Foo(object):
        ...       lazy_import('sage.all', 'plot')
        sage: class Bar(Foo):
        ...       pass
        sage: type(Foo.__dict__['plot'])
        <type 'sage.misc.lazy_import.LazyImport'>
        sage: 'EXAMPLES' in Bar.plot.__doc__
        True
        sage: type(Foo.__dict__['plot'])
        <type 'function'>
    """
    if _as is None:
        _as = names
    if isinstance(names, str):
        names = [names]
        _as = [_as]
    if namespace is None:
        namespace = inspect.currentframe().f_locals
    if "*" in names:
        ix = names.index("*")
        names[ix:ix+1] = get_star_imports(module)
        _as[ix:ix+1] = [None] * (len(names) - len(_as) + 1)
    for name, alias in zip(names, _as):
        if not overwrite and (alias or name) in namespace:
            continue
        namespace[alias or name] = LazyImport(module, name, alias, namespace)


star_imports = None

def save_cache_file():
    """
    Used to save the cached * import names.

    TESTS::

        sage: import sage.misc.lazy_import
        sage: sage.misc.lazy_import.save_cache_file()
    """
    global star_imports
    if star_imports is None:
        star_imports = {}
    _, tmp_file = tempfile.mkstemp()
    pickle.dump(star_imports, open(tmp_file, "w"))
    shutil.move(tmp_file, get_cache_file())

def get_star_imports(module_name):
    """
    Lookup the list of names in a module that would be imported with "import \\*"
    either via a cache or actually importing.

    EXAMPLES::

        sage: from sage.misc.lazy_import import get_star_imports
        sage: 'get_star_imports' in get_star_imports('sage.misc.lazy_import')
        True
        sage: 'EllipticCurve' in get_star_imports('sage.schemes.all')
        True
    """
    global star_imports
    if star_imports is None:
        cache_file = get_cache_file()
        if os.path.exists(cache_file):
            star_imports = pickle.load(open(cache_file))
        else:
            star_imports = {}
    try:
        return star_imports[module_name]
    except KeyError:
        module = __import__(module_name, {}, {}, ["*"])
        if hasattr(module, "__all__"):
            all = module.__all__
        else:
            all = [key for key in dir(module) if key[0] != "_"]
        star_imports[module_name] = all
        return all
