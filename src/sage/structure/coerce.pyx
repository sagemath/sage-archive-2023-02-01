#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "../ext/stdsage.pxi"
include "coerce.pxi"

import operator

cdef class CoercionModel_simple(CoercionModel):
    """
    This is the original coercion model, as of SAGE 2.6 (2007-06-02)
    """

    cdef canonical_coercion_c(self, x, y):
        cdef int i
        xp = parent_c(x)
        yp = parent_c(y)
        if xp is yp:
            return x, y

        if PY_IS_NUMERIC(x):
            try:
                x = yp(x)
            except TypeError:
                y = x.__class__(y)
                return x, y
            # Calling this every time incurs overhead -- however, if a mistake
            # gets through then one can get infinite loops in C code hence core
            # dumps.  And users define _coerce_ and __call__ for rings, which
            # can easily have bugs in it, i.e., not really make the element
            # have the correct parent.  Thus this check is *crucial*.
            return _verify_canonical_coercion_c(x,y)

        elif PY_IS_NUMERIC(y):
            try:
                y = xp(y)
            except TypeError:
                x = y.__class__(x)
                return x, y
            return _verify_canonical_coercion_c(x,y)

        try:
            if xp.has_coerce_map_from(yp):
                y = (<Parent>xp)._coerce_c(y)
                return _verify_canonical_coercion_c(x,y)
        except AttributeError:
            pass
        try:
            if yp.has_coerce_map_from(xp):
                x = (<Parent>yp)._coerce_c(x)
                return _verify_canonical_coercion_c(x,y)
        except AttributeError:
            pass
        raise TypeError, "no common canonical parent for objects with parents: '%s' and '%s'"%(xp, yp)

    cdef canonical_base_coercion_c(self, Element x, Element y):
        if not have_same_base(x, y):
            if (<Parent> x._parent._base).has_coerce_map_from_c(y._parent._base):
                # coerce all elements of y to the base ring of x
                y = y.base_extend_c(x._parent._base)
            elif (<Parent> y._parent._base).has_coerce_map_from_c(x._parent._base):
                # coerce x to have elements in the base ring of y
                x = x.base_extend_c(y._parent._base)
        return x,y

    def canonical_base_coercion(self, x, y):
        try:
            xb = x.base_ring()
        except AttributeError:
            #raise TypeError, "unable to find base ring for %s (parent: %s)"%(x,x.parent())
            raise TypeError, "unable to find base ring"
        try:
            yb = y.base_ring()
        except AttributeError:
            raise TypeError, "unable to find base ring"
            #raise TypeError, "unable to find base ring for %s (parent: %s)"%(y,y.parent())
        try:
            b = self.canonical_coercion_c(xb(0),yb(0))[0].parent()
        except TypeError:
            raise TypeError, "unable to find base ring"
            #raise TypeError, "unable to find a common base ring for %s (base ring: %s) and %s (base ring %s)"%(x,xb,y,yb)
        return x.change_ring(b), y.change_ring(b)


    cdef bin_op_c(self, x, y, op):
        """
        Compute x op y, where coercion of x and y works according to
        SAGE's coercion rules.
        """
        # Try canonical element coercion.
        try:
            x1, y1 = self.canonical_coercion_c(x, y)
            return op(x1,y1)
        except TypeError, msg:
            #print msg  # this can be useful for debugging.
            if not op is operator.mul:
                raise TypeError, arith_error_message(x,y,op)

        # If the op is multiplication, then some other algebra multiplications
        # may be defined

        # 2. Try scalar multiplication.
        # No way to multiply x and y using the ``coerce into a canonical
        # parent'' rule.
        # The next rule to try is scalar multiplication by coercing
        # into the base ring.
        cdef bint x_is_modelt, y_is_modelt

        y_is_modelt = PY_TYPE_CHECK(y, ModuleElement)
        if y_is_modelt:
            # First try to coerce x into the base ring of y if y is an element.
            try:
                R = (<ModuleElement> y)._parent._base
                if R is None:
                    raise RuntimeError, "base of '%s' must be set to a ring (but it is None)!"%((<ModuleElement> y)._parent)
                x = (<Parent>R)._coerce_c(x)
                return (<ModuleElement> y)._rmul_c(x)     # the product x * y
            except TypeError, msg:
                pass

        x_is_modelt = PY_TYPE_CHECK(x, ModuleElement)
        if x_is_modelt:
            # That did not work.  Try to coerce y into the base ring of x.
            try:
                R = (<ModuleElement> x)._parent._base
                if R is None:
                    raise RuntimeError, "base of '%s' must be set to a ring (but it is None)!"%((<ModuleElement> x)._parent)
                y = (<Parent> R)._coerce_c(y)
                return (<ModuleElement> x)._lmul_c(y)    # the product x * y
            except TypeError:
                pass

        if y_is_modelt and x_is_modelt:
            # 3. Both canonical coercion failed, but both are module elements.
            # Try base extending the right object by the parent of the left

            ## TODO -- WORRY -- only unambiguous if one succeeds!
            if  PY_TYPE_CHECK(x, RingElement):
                try:
                    return x * y.base_extend((<RingElement>x)._parent)
                except (TypeError, AttributeError), msg:
                    pass
            # Also try to base extending the left object by the parent of the right
            if  PY_TYPE_CHECK(y, RingElement):
                try:
                    return y * x.base_extend((<Element>y)._parent)
                except (TypeError, AttributeError), msg:
                    pass

        # 4. Try _l_action or _r_action.
        # Test to see if an _r_action or _l_action is
        # defined on either side.
        try:
            return x._l_action(y)
        except (AttributeError, TypeError):
            pass
        try:
            return y._r_action(x)
        except (AttributeError, TypeError):
            pass

        raise TypeError, arith_error_message(x,y,op)

cdef class CoercionModel_cache_maps(CoercionModel_simple):

    cdef coercion_maps_c(self, R, S):
        try:
            return self._coercion_maps[R,S]
        except KeyError:
            homs = self.discover_coercion_c(R, S)
            self._coercion_maps[R,S] = homs
            self._coercion_maps[S,R] = (homs[1], homs[0])
            return homs

    cdef action_maps_c(self, R, S, op):
        try:
            return self._action_maps[R,S,op]
        except KeyError:
            action = self.discover_action_c(R, S, op)
            self._action_maps[R,S,op] = action
            return action

    cdef discover_coercion_c(self, R, S):
        from sage.categories.homset import Hom
        if R is S:
            return None, None

        # See if there is a natural coercion from R to S
        try:
            return Hom(R, S).natural_map(), None
        except TypeError:
            pass

        # See if there is a natural coercion from S to R
        try:
            return None, Hom(S, R).natural_map()
        except TypeError:
            pass


    cdef discover_action_c(self, R, S, op):
        pass
