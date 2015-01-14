"""
Abstract base class for modules

Two classes for modules are defined here: :class:`.Module_old` and
:class:`.Module`. The former is a legacy class which should not be used for new
code anymore as it does not conform to the coercion model. Eventually all
occurences shall be ported to :class:`.Module`.

AUTHORS:

- William Stein: initial version

- Julian Rueth (2014-05-10): category parameter for Module, doc cleanup

EXAMPLES:

A minimal example of a module::

    sage: class MyElement(sage.structure.element.ModuleElement):
    ....:     def __init__(self, parent, x):
    ....:         self.x = x
    ....:         sage.structure.element.ModuleElement.__init__(self, parent=parent)
    ....:     def _rmul_(self, c):
    ....:         return self.parent()(c*self.x)
    ....:     def _add_(self, other):
    ....:         return self.parent()(self.x + other.x)
    ....:     def __cmp__(self, other):
    ....:         return cmp(self.x, other.x)
    ....:     def __hash__(self):
    ....:         return hash(self.x)
    ....:     def _repr_(self):
    ....:         return repr(self.x)

    sage: class MyModule(sage.modules.module.Module):
    ....:     Element = MyElement
    ....:     def _element_constructor_(self, x):
    ....:         if isinstance(x, MyElement): x = x.x
    ....:         return self.element_class(self, self.base_ring()(x))
    ....:     def __cmp__(self, other):
    ....:         if not isinstance(other, MyModule): return cmp(type(other),MyModule)
    ....:         return cmp(self.base_ring(),other.base_ring())

    sage: M = MyModule(QQ)
    sage: M(1)
    1

    sage: import __main__
    sage: __main__.MyModule = MyModule
    sage: __main__.MyElement = MyElement
    sage: TestSuite(M).run()

"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
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

from sage.structure.parent_gens cimport ParentWithAdditiveAbelianGens
from sage.structure.parent cimport Parent

cdef class Module_old(sage.structure.parent_gens.ParentWithAdditiveAbelianGens):
    """
    Generic module class.
    """
    def category(self):
        """
        Return the category to which this module belongs.
        """
        # Defining a category method is deprecated for parents.
        # Instead, the category should be specified in the constructor.
        # See: http://sagetrac.org/sage_trac/wiki/CategoriesRoadMap
        if self._is_category_initialized():
            return Parent.category(self)
        from sage.categories.modules import Modules
        return Modules(self.base_ring())

    def endomorphism_ring(self):
        """
        Return the endomorphism ring of this module in its category.
        """
        try:
            return self.__endomorphism_ring
        except AttributeError:
            import sage.categories.all
            E = sage.categories.all.End(self)
            # the following is invalid code, you can't add attributes
            # to Cython classes like this. I guess nobody calls this
            # method.
            self.__endomorphism_ring = E
            return E

cdef class Module(sage.structure.parent.Parent):
    """
    Generic module class.

    INPUT:

    - ``base`` -- a ring. The base ring of the module.

    - ``category`` -- a category (default: ``None``), the category for this
      module. If ``None``, then this is set to the category of modules/vector
      spaces over ``base``.

    EXAMPLES::

        sage: from sage.modules.module import Module
        sage: M = Module(ZZ)
        sage: M.base_ring()
        Integer Ring
        sage: M.category()
        Category of modules over Integer Ring

    Normally the category is set to the category of modules over ``base``. If
    ``base`` is a field, then the category is the category of vector spaces
    over ``base``::

        sage: M_QQ = Module(QQ)
        sage: M_QQ.category()
        Category of vector spaces over Rational Field

    The ``category`` parameter can be used to set a more specific category::

        sage: N = Module(ZZ, category=FiniteDimensionalModulesWithBasis(ZZ))
        sage: N.category()
        Category of finite dimensional modules with basis over Integer Ring

     TESTS:

     We check that :trac:`8119` has been resolved::

        sage: M = ZZ^3
        sage: h = M.__hash__()
        sage: M.rename('toto')
        sage: h == M.__hash__()
        True

    """
    def __init__(self, base, category=None):
        """
        Initialization.

        TESTS::

            sage: from sage.modules.module import Module
            sage: M = Module(ZZ)
            sage: type(M)
            <type 'sage.modules.module.Module'>

        """
        from sage.categories.modules import Modules
        if category is None:
            category = Modules(base)
        Parent.__init__(self, base=base, category=category)

    def endomorphism_ring(self):
        """
        Return the endomorphism ring of this module in its category.

        EXAMPLES::

            sage: from sage.modules.module import Module
            sage: M = Module(ZZ)
            sage: M.endomorphism_ring()
            Set of Morphisms from <type 'sage.modules.module.Module'> to <type 'sage.modules.module.Module'> in Category of modules over Integer Ring

        """
        from sage.categories.all import End
        return End(self)

def is_Module(x):
    """
    Return ``True`` if ``x`` is a module, ``False`` otherwise.

    INPUT:

    - ``x`` -- anything.

    EXAMPLES::

        sage: from sage.modules.module import is_Module
        sage: M = FreeModule(RationalField(),30)
        sage: is_Module(M)
        True
        sage: is_Module(10)
        False

    """
    return isinstance(x, Module) or isinstance(x, Module_old)

def is_VectorSpace(x):
    """
    Return ``True`` if ``x`` is a vector space, ``False`` otherwise.

    INPUT:

    - ``x`` -- anything.

    EXAMPLES::

        sage: from sage.modules.module import is_Module, is_VectorSpace
        sage: M = FreeModule(RationalField(),30)
        sage: is_VectorSpace(M)
        True
        sage: M = FreeModule(IntegerRing(),30)
        sage: is_Module(M)
        True
        sage: is_VectorSpace(M)
        False

    """
    try:
        return is_Module(x) and x.base_ring().is_field()
    except AttributeError:
        return False


