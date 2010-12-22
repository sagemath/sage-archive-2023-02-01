"""
Abstract base class for modules
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import random
from sage.structure.parent_gens cimport ParentWithAdditiveAbelianGens
from sage.structure.parent cimport Parent



# This is the old Module class which does not conform to the new
# coercion model. Eventually all occurences shall be ported to the new
# Module class, which is defined below.

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


    def is_atomic_repr(self):
        """
        True if the elements have atomic string representations, in the
        sense that they print if they print at s, then -s means the
        negative of s. For example, integers are atomic but polynomials are
        not.
        """
        return False

    def __hash__(self):
        return hash(self.__repr__())





### The new Module class that should be the base of all Modules
### The derived Module class must implement the element
### constructor:
#
# class MyModule(sage.modules.module.Module):
#     Element = MyElement
#     def _element_constructor_(self, x):
#         return self.Element(x)
#
### The Element should also implement _rmul_ (or _lmul_)
#
# class MyElement(sage.structure.element.ModuleElement):
#     def _rmul_(self, c):
#         ...

cdef class Module(sage.structure.parent.Parent):
    """
    Generic module class.

    EXAMPLES::

        sage: from sage.modules.module import Module
        sage: M = Module(ZZ)
        sage: M.category()
        Category of modules over Integer Ring
        sage: M.category().required_methods()
        {'parent': {'required': ['__contains__'],
                    'optional': []},
         'element': {'required': [],
                     'optional': ['_add_']}}
        sage: M_QQ = Module(QQ)
        sage: M_QQ.category()
        Category of vector spaces over Rational Field
    """

    def __init__(self, base):
        """
        Construct a module and set the category.

        INPUT:

        - ``base`` -- a ring. The base ring of the module.

        EXAMPLES::

            sage: from sage.modules.module import Module
            sage: M = Module(ZZ); M
            <type 'sage.modules.module.Module'>
            sage: M.base_ring()
            Integer Ring
        """
        from sage.categories.modules import Modules
        Parent.__init__(self, base=base, category=Modules(base))


    def endomorphism_ring(self):
        """
        Return the endomorphism ring of this module in its category.

        EXAMPLES::

            sage: from sage.modules.module import Module
            sage: M = Module(ZZ); M
            <type 'sage.modules.module.Module'>
            sage: M.endomorphism_ring()
            Set of Morphisms from <type 'sage.modules.module.Module'> to
            <type 'sage.modules.module.Module'> in Category of
            modules over Integer Ring
        """
        from sage.categories.all import End
        return End(self)


    def is_atomic_repr(self):
        """
        Whether the elements have atomic string representations.

        OUTPUT:

        Boolean. ``True`` if the elements have atomic string
        representations, in the sense that they print if they print
        ``s``, then ``-s`` means the negative of ``s``. For example,
        integers are atomic but polynomials are not.

        EXAMPLES::

            sage: from sage.modules.module import Module
            sage: M = Module(ZZ)
            sage: M.is_atomic_repr()
            False
            sage: ZZ.is_atomic_repr()
            True
        """
        return False


    def __hash__(self):
        """
        Return a hash value for the :class:`Module` instance.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: from sage.modules.module import Module
            sage: M = Module(ZZ)
            sage: M.__hash__()
            6190647798068218210
        """
        return hash(self.__repr__())



def is_Module(x):
    """
    Return True if x is a module.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

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
    Return True if x is a vector space.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

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


