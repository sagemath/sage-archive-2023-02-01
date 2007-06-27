#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "../ext/stdsage.pxi"
include "../ext/python_tuple.pxi"
include "coerce.pxi"

import operator
import sage.categories.morphism

cdef class CoercionModel_original(CoercionModel):
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
            # print msg  # this can be useful for debugging.
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


# just so we can detect these fast to avoid action-searching
cdef op_add, op_sub
from operator import add as op_add, sub as op_sub

cdef class CoercionModel_cache_maps(CoercionModel_original):

    def __init__(self):
        # This MUST be a mapping of tuples, where each
        # tuple contains at least two elements that are either
        # None or of type Morphism.
        self._coercion_maps = {}
        # This MUST be a mapping of actions.
        self._action_maps = {}

    cdef bin_op_c(self, x, y, op):

        if (op is not op_add) and (op is not op_sub):
            # Actions take preference over common-parent coercions.
            xp = parent_c(x)
            yp = parent_c(y)
            if xp is yp:
                return op(x,y)
            action = self.action_maps_c(xp, yp, op)
            if action is not None:
                return (<Action>action)._call_c(x, y)

        try:
            xy = self.canonical_coercion_c(x,y)
            return op(<object>PyTuple_GET_ITEM(xy, 0), <object>PyTuple_GET_ITEM(xy, 1))

        except TypeError:
            pass

        if op is operator.mul:
            if isinstance(x, (str, list, tuple)):
                # __mul__ overridden with special meaning for sequences
                if y == int(y):
                    return x * int(y)

            # TODO: roll into actions (here now so doctests pass)
            try:
                return x._l_action(y)
            except (AttributeError, TypeError):
                pass
            try:
                return y._r_action(x)
            except (AttributeError, TypeError):
                pass

        raise TypeError, arith_error_message(x,y,op)



    cdef canonical_coercion_c(self, x, y):
        xp = parent_c(x)
        yp = parent_c(y)
        if xp is yp:
            return x,y

        cdef Element x_elt, y_elt
        coercions = self.coercion_maps_c(xp, yp)
        if coercions is not None:
            # Unsafe access for speed
            x_map = <object>PyTuple_GET_ITEM(coercions, 0)
            y_map = <object>PyTuple_GET_ITEM(coercions, 1)
            if x_map is not None:
                x_elt = (<Morphism>x_map)._call_c(x)
            else:
                x_elt = x
            if y_map is not None:
                y_elt = (<Morphism>y_map)._call_c(y)
            else:
                y_elt = y
            if x_elt._parent is y_elt._parent:
                # We must verify this as otherwise we are prone to
                # getting into an infinite loop in c, and the above
                # morphisms may be written by (imperfect) users.
                return x_elt,y_elt
            elif x_elt._parent == y_elt._parent:
                # TODO: Non-uniqueness of parents strikes again!
                # print parent_c(x_elt), " is not ", parent_c(y_elt)
                y_elt = parent_c(x_elt)(y_elt)
                if x_elt._parent is y_elt._parent:
                    return x_elt,y_elt
            self._coercion_error(x, x_map, x_elt, y, y_map, y_elt)

        # Now handle the native python + sage object cases
        # that were not taken care of above.
        elif PY_IS_NUMERIC(x):
            try:
                x = yp(x)
                if PY_TYPE_CHECK(yp, type): return x,y
            except TypeError:
                y = x.__class__(y)
                return x, y
            return _verify_canonical_coercion_c(x,y)

        elif PY_IS_NUMERIC(y):
            try:
                y = xp(y)
                if PY_TYPE_CHECK(xp, type): return x,y
            except TypeError:
                x = y.__class__(x)
                return x, y
            return _verify_canonical_coercion_c(x,y)

        raise TypeError, "no common canonical parent for objects with parents: '%s' and '%s'"%(xp, yp)


    def _coercion_error(self, x, x_map, x_elt, y, y_map, y_elt):
        raise RuntimeError, """There is a bug in the coercion code in SAGE.
Both x (=%r) and y (=%r) are supposed to have identical parents but they don't.
In fact, x has parent '%s'
whereas y has parent '%s'

Original elements %r (parent %s) and %r (parent %s) and morphisms
%s %r
%s %r"""%( x_elt, y_elt, parent_c(x_elt), parent_c(y_elt),
            x, parent_c(x), y, parent_c(y),
            type(x_map), x_map, type(y_map), y_map)

    def coercion_maps(self, R, S):
        return self.coercion_maps_c(R, S)

    cdef coercion_maps_c(self, R, S):
        try:
            return self._coercion_maps[R,S]
        except KeyError:
            homs = self.discover_coercion_c(R, S)
            if homs is not None:
                self._coercion_maps[R,S] = homs
                self._coercion_maps[S,R] = (homs[1], homs[0])
            else:
                self._coercion_maps[R,S] = self._coercion_maps[S,R] = None
            return homs

    def action_maps(self, R, S, op):
        return self.action_maps_c(R, S, op)

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
        if PY_TYPE_CHECK(R, Parent):
            mor = (<Parent>R).coerce_map_from_c(S)
            if mor is not None:
                return None, mor

        # See if there is a natural coercion from S to R
        if PY_TYPE_CHECK(S, Parent):
            mor = (<Parent>S).coerce_map_from_c(R)
            if mor is not None:
                return mor, None

        # Try base extending to left and right
        # TODO: This is simple and ambiguous, add sophistication
        if PY_TYPE_CHECK(R, ParentWithBase) and PY_TYPE_CHECK(S, Parent):
            if (<Parent>S).coerce_map_from_c((<ParentWithBase>R)._base) is not None:
                Z = R.base_extend(S) # should there be a base-extension morphism?
            elif (<Parent>R).coerce_map_from_c((<ParentWithBase>S)._base) is not None:
                Z = S.base_extend(R)
            else:
                Z = None
            if Z is not None:
                from sage.categories.homset import Hom
                # Can I trust always __call__() to do the right thing in this case?
                return sage.categories.morphism.CallMorphism(Hom(R, Z)), sage.categories.morphism.CallMorphism(Hom(S, Z))

        return None

    cdef discover_action_c(self, R, S, op):
        return None # for now
        if op is operator.mul:
            try:
                return LeftModuleAction(R, S)
            except TypeError:
                pass

            try:
                return RightModuleAction(S, R)
            except TypeError:
                pass

        elif op is operator.div:
            try:
                return ~RightModuleAction(S, R)
            except TypeError:
                pass



from sage.structure.element cimport Element # workaround SageX bug


cdef class LeftModuleAction(Action):
    def __init__(self, G, S):
        # TODO: detect this better
        cdef ModuleElement a
        cdef RingElement g
        # if this is bad it will raise a type error in the subsequent lines, which we propagate
        g = rand_elt(G)
        a = rand_elt(S)
        res = a._rmul_c(g) # g * a
        print g, a, res
        if parent_c(res) is not S:
            raise TypeError
        Action.__init__(self, G, S, 1)

    cdef Element _call_c_impl(self, Element g, Element a):
        return (<ModuleElement>a)._rmul_c(g)

cdef class RightModuleAction(Action):
    def __init__(self, G, S):
        # TODO: detect this better
        cdef ModuleElement a
        cdef RingElement g
        # if this is bad it will raise a type error in the subsequent lines, which we propagate
        g = rand_elt(G)
        a = rand_elt(S)
        res = a._lmul_c(g) # a * g
        if parent_c(res) is not S:
            raise TypeError
        Action.__init__(self, G, S, 0)

    cdef Element _call_c_impl(self, Element a, Element g):
        return (<ModuleElement>a)._lmul_c(g)


def rand_elt(R):
    """
    Want to return something sufficiently generic that __call__()
    won't work in general (as _lmul_ and _rmul_ are to liberal).
    """
    try:
        z = R.random_element()
        if z not in ZZ:
            return z
    except:
        pass
    try:
        from sage.functions.constants import pi
        return R(pi)
    except:
        pass
    try:
        from sage.rings.integer_ring import ZZ
        return R.gen()*ZZ.random_element()
    except:
        pass
    try:
        from sage.rings.integer_ring import ZZ
        return R(ZZ.random_element())
    except:
        pass
    raise TypeError