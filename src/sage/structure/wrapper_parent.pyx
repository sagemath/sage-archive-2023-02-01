"""
This module provides a general class that wraps other parents and element types.

The purpose of this functionality is to be able to add morphisms and actions to the coercion system that interact with objects that are otherwise unique.  For example, one would use wrappers to implement:
0. Residue Fields where one has a fast implementation for the arithmetic of the quotient, for example $O_K / p O_K$ for $K$ a number field, $O_K$ its maximal order and $p$ a prime of $O_K$.  This residue field includes a coercion map from $O_K$.
1. Completions at a places of number fields.  There is a morphism added that maps to the completed ring, which is represented as a finite extension of $\QQ_p$
2. G-sets for groups acting on an object that would otherwise be unique.

There are three options for implementing such a wrapper:
Option 0: WrapperParent_model0
  With this option, one inherits from the class WrapperParent_model0.  Elements are not wrapped, but instead merely have their parents changed to the wrapping parent.  This will work fine for most parent and element types, but can break when parents and elements are either highly optimized or implemented strangely.  In particularly, the following are failure modes for this option:
    a) If the parent caches elements (as in IntegerMod_int), then the elements in that cache will point to the wrong parent.  Thus the sum of two elements of the wrapper could be an element of the wrapped object, which is undesirable.
    b) If the parent is written in Cython, more work must be done, which is not necessarily possible.  See FiniteField_givaro for an example where option 0 just fails (because various methods in that class assume that the parent is of type FiniteField_givaro).
  This option should generally be fast because there are no wrappers around elements.  It should also not be all that difficult to implement.
  When inheriting from this class, you should do the following:
    a) Check to ensure that your element's _make_new_with_parent_c correctly returns an element with parent set to the given argument (default behavior is to set parent and return self).  This would not work in IntegerMod_int for example, because self.__modulus.table would still be wrong.
    b) Make sure that your elements don't call cdefed functions or access cdefed attributes of the parent.  If they do, you need to duplicate that behavior in your subclass of the wrapper.
    c) Override __call__
    d) Pass in the correct morphisms/actions to the __init__ method of WrapperParent_model0.  This may involve writing a custom Morphism extension.
    e) Write any functions that you want in addition to the functions of the wrapped parent or to replace the wrapped parents method.
    f) Update this list to reflect any lessons you've learned.

Option 1: WrapperParent_model1
  With this option, one wraps both the parent and the element.  On the downside, there will be a speed penalty to arithmetic.  On the plus side, this method should work for classes that fail option 0.  It should also be easier to work with.
  When using this option, you should do the following:
    a) Write any morphisms that you need in order to make the desired coercion work.  Construct a WrapperParent_model1 with the appropriate arguments.
    b) Update this list to reflect any lessons you've learned.

Option 2: Write your own custom.
  When you want speed but can't use option 1.  See sage/rings/finite_rings/residue_field.pyx for examples (ResidueFiniteField_prime_modn, ResidueFiniteField_givaro and ResidueFiniteField_pari_ffelt)

AUTHORS:
  -- David Roe (2007-10-3)
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/cdefs.pxi"
include "sage/ext/stdsage.pxi"

import copy
import inspect
from sage.categories.all import Objects
from sage.categories.homset import Hom
from sage.rings.ring import Ring
from sage.categories.morphism import CallMorphism

cdef class StealMorphism(Morphism):
    """
    Represents a morphism S -> R that just sets the parent of x to R.
    """
    def _repr_type(self):
        return "Steal"

    cpdef Element _call_(self, x):
        """
        Applies this morphism to the element x.

        NB: make sure the function _make_new_with_parent_c changes all of the state of the element to make it behave correctly with self.codomain() as a parent.
        """
        y = copy.copy(x)
        return (<Element>y)._make_new_with_parent_c(self.codomain())

cdef class WrapperParent_model0(Parent):
    """
    This class is designed to wrap around unique parents to provide facility for different coercions, actions and embeddings.

    NB: See the comments at the top of the module for advantages and disadvantages of this class.
    """
    def __init__(self, R, cat = None, coerce_from = [], actions = [], embeddings = []):
        self.R = R
        if cat is None:
            cat = Objects()
        stealS = StealMorphism(Hom(R, self, cat))
        stealR = StealMorphism(Hom(self, R, cat))
        # need to modify actions
        Parent.__init__(self, R.base(), [stealS * f for f in coerce_from], actions, [f * stealR for f in embeddings])

    def __call__(self, x):
        """
        This call method needs to basically create a copy and then change appropriate fields such as parent so that the element behaves correctly as an element in self.
        """
        raise NotImplementedError

    def trait_names(self):
        """
        This function is here to provide tab completion for this parent.
        """
        return [e[0] for e in inspect.getmembers(self.R)]

    def __getattr__(self, name):
        """
        Transparently passes calls to methods and attributes of self down to the wrapped ring.

        Note that this doesn't work for cdefed methods and attributes, so if you parent is Cythoned you need to deal with such specially.
        """
        try:
            return self.R.__getattribute__(name)
        except AttributeError:
            return Parent.__getattribute__(name)

cdef class WrapperParent_model1(Parent):
    """
    """
    def __init__(self, R, cat = None, coerce_from = [], actions = [], embeddings = []):
        self.R = R
        if cat is None:
            cat = Objects()
        stealS = CallMorphism(Hom(R, self, cat))
        stealR = CallMorphism(Hom(self, R, cat))
        # need to modify actions
        Parent.__init__(self, R.base(), [stealS * f for f in coerce_from], actions, [f * stealR for f in embeddings])

    def __call__(self, x):
        """
        Creates an element of self from x.

        NB: You need to override this if you subclass the WrapperElement class.
        """
        return WrapperElement(self, x)

    def trait_names(self):
        """
        This function is here to provide tab completion for this parent.
        """
        return [e[0] for e in inspect.getmembers(self.R)]

    def __getattr__(self, name):
        """
        Transparently passes calls to methods and attributes of self down to the wrapped ring.

        Note that this doesn't work for cdefed methods and attributes, so if you parent is Cythoned you need to deal with such specially.
        """
        try:
            return self.R.__getattribute__(name)
        except AttributeError:
            return Parent.__getattribute__(name)

cdef class WrapperElement(AlgebraElement):
    """
    Should work for Elements, ModuleElements, RingElements and AlgebraElements
    """
    def __init__(self, parent, x):
        self.val = x
        Element.__init__(parent, x)

    def trait_names(self):
        """
        This function is here to provide tab completion for this parent.
        """
        return [e[0] for e in inspect.getmembers(self.val)]

    def __getattr__(self, name):
        """
        Transparently passes calls to methods and attributes of self down to the wrapped ring.

        Note that this doesn't work for cdefed methods and attributes, so if you parent is Cythoned you need to deal with such specially.
        """
        try:
            return self.val.__getattribute__(name)
        except AttributeError:
            return AlgebraElement.__getattribute__(name)

    cdef _richcmp_c_impl(left, Element right, int op):
        return (<Element>left.val)._richcmp_c(right, op)
    cdef int _cmp_c_impl(left, Element right) except -2:
        return (<Element>left.val)._cmp_c(right)
    cdef base_extend_c_impl(self, Parent R):
        return (<Element>self.val).base_extend_c(R)
    cpdef ModuleElement _add_(self, ModuleElement right):
        return self.val + (<WrapperElement>right).val
    cpdef ModuleElement _sub_(self, ModuleElement right):
        return self.val - (<WrapperElement>right).val
    cpdef ModuleElement _neg_(self):
        return -self.val
    cpdef ModuleElement _lmul_(self, RingElement right):
        return self.val * right
    cpdef ModuleElement _rmul_(self, RingElement left):
        return left * self.val
    cdef RingElement coerce_to_base_ring(self, x):
        if PY_TYPE_CHECK(self.val, RingElement):
            return (<RingElement>self.val).coerce_to_base_ring(x)
        else:
            raise TypeError, "self is not actually a RingElement"
    cpdef RingElement _mul_(self, RingElement right):
        return self.val * (<WrapperElement>right).val
    cpdef RingElement _div_(self, RingElement right):
        return self.val / (<WrapperElement>right).val

