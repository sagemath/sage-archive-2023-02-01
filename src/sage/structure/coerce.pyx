"""
Coercion helper functions
"""

#*****************************************************************************
#       Copyright (C) 2004-2005 William Stein <wstein@gmail.com>
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

import __builtin__
import operator

import  element

cimport sage.modules.module
import  sage.modules.module

def lcm(x,y):
    from sage.rings.arith import lcm
    return lcm(x,y)

def gcd(x,y):
    from sage.rings.arith import gcd
    return gcd(x,y)

def xgcd(x,y):
    from sage.rings.arith import xgcd
    return xgcd(x,y)

def parent(x):
    try:
        return x.parent()
    except AttributeError:
        return type(x)

def coerce(p, x):
    try:
        return p._coerce_(x)
    except AttributeError:
        return p(x)


def canonical_coercion(x, y):
    cdef int i
    xp = parent(x)
    yp = parent(y)
    if xp is yp:
        return x, y

    try:
        if x.__class__ in [int, long, float, complex]:
            try:
                x = yp(x)
            except TypeError:
                y = x.__class__(y)
                return x, y
        elif y.__class__ in [int, long, float, complex]:
            try:
                y = xp(y)
            except TypeError:
                x = y.__class__(x)
                return x, y
        else:
            i = 0
            try:
                x0 = x
                x = coerce(yp, x)
            except TypeError, msg:
                i = i + 1
            try:
                y = coerce(xp, y)
            except TypeError, msg:
                i = i + 1
            if i == 0:
                # Both succeed.  But we must be careful to take x before
                # it was coerced, or we end up *switching* to the parents,
                # which is no good.
                return x0, y
            elif i == 2:
                import  sage.rings.ring
                if isinstance(x, sage.rings.ring.Ring) or isinstance(y, sage.rings.ring.Ring):
                    raise TypeError, "you cannot +,*,/ a ring with a number."
                raise TypeError, "unable to find a common parent for %s (parent: %s) and %s (parent: %s)"%(x,xp, y, yp)
        return x, y
    except AttributeError:
        raise TypeError, "unable to find a common canonical parent"

def canonical_base_coercion(x, y):
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
        b = canonical_coercion(xb(0),yb(0))[0].parent()
    except TypeError:
        raise TypeError, "unable to find base ring"
        #raise TypeError, "unable to find a common base ring for %s (base ring: %s) and %s (base ring %s)"%(x,xb,y,yb)
    return x.change_ring(b), y.change_ring(b)

def bin_op(x, y, op):
    """
    Compute x op y, where coercion of x and y works according to SAGE's coercion
    rules.
    """
    #print "bin_op(%s,%s,%s)"%(x,y,op)   # debug
    if isinstance(y, element.InfinityElement):
        return op(y,x)
    if op == operator.mul and \
           isinstance(y, (\
                      element.ModuleElement,
                      element.AlgebraElement,
                      sage.modules.module.Module)) and \
           isinstance(x, (element.RingElement, int, long, float)):
        return op(y,x)
    try:
        #print 1, x, y, x.parent(), y.parent()
        x, y = canonical_coercion(x, y)
        #print 2, x, y, x.parent(), y.parent()
    except TypeError, mesg:
        try:
            return y._r_action(x)
        except AttributeError:
            raise TypeError, '%s: x parent: %s, y parent: %s'%(mesg, parent(x), parent(y))
        except TypeError:
            raise TypeError, "No right action defined"
        try:
            return x._l_action(y)
        except AttributeError:
            raise TypeError, mesg
        except TypeError:
            raise TypeError, "No left action defined"
            #raise TypeError, "No left action of %s on %s defined"%(x,y)

    return op(x,y)


N = type(None)

cdef class Coerce:
    cdef cmp_c(self, x, y):
        tx = type(x); ty = type(y)
        if (tx == N and ty != N) or (tx != N and ty == N):
            return -1
        elif isinstance(y, element.InfinityElement):
            return -y.__cmp__(x)

        xp = parent(x)
        yp = parent(y)
        if xp is yp:
            return __builtin__.cmp(x,y)

        cdef int fails
        fails = 0
        if isinstance(x, (int, long)):

            return __builtin__.cmp(yp(x), y)

        elif isinstance(y, (int, long)):

            return __builtin__.cmp(x, xp(y))

        else:

            fails = 0
            try:
                x0 = x
                x = coerce(yp, x)
            except (TypeError, ValueError):
            #    print "coercing %s to %s failed"%(x0,yp)
                fails = fails + 1
            #else:
            #    print "coercion %s to %s suceeded with %s (parent=%s)"%(x0,yp,x,parent(x))

            try:
                y0 = y
                y = coerce(xp, y)
            except (TypeError, ValueError):
            #    print "coercing %s to %s failed"%(y0,xp)
                fails = fails + 1
            #else:
            #    print "coercion %s to %s suceeded with %s (parent=%s)"%(y0,xp,y,parent(y))


            if fails == 0:
                assert (parent(x0) is parent(y))  # debug
                return __builtin__.cmp(x0,y)

            elif fails == 2:

                return -1

            else:
                if not (parent(x) is parent(y)):
                    raise RuntimeError, "There is a bug in coercion: x=%s (parent=%s), y=%s (parent=%s)"%(x, parent(x), y, parent(y))
                return __builtin__.cmp(x,y)


# TODO -- don't use a classic, use "extern public".

cdef Coerce functions
functions = Coerce()

def cmp(x,y):  # external interface to cmp_cdef
    return functions.cmp_c(x,y)

