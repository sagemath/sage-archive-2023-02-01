"""
Coercion helper functions
"""

#*****************************************************************************
#       Copyright (C) 2004-2005 William Stein <wstein@ucsd.edu>
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

cimport element
import  element

cimport module
import  module

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

    if isinstance(x, (int, long, float, complex)):
        x = yp(x)
    elif isinstance(y, (int, long, float, complex)):
        y = xp(y)
    else:
        i = 0
        try:
            x = coerce(yp, x)
        except TypeError, ValueError:
            i = i + 1
        try:
            y = coerce(xp, y)
        except TypeError, ValueError:
            i = i + 1
        if i == 0:
            raise TypeError, "unable to find an unambiguous parent for %s (parent: %s) and %s (parent: %s)"%(x,xp, y, yp)
        elif i == 2:
            import  sage.ext.ring
            if isinstance(x, sage.ext.ring.Ring) or isinstance(y, sage.ext.ring.Ring):
                raise TypeError, "you cannot combine ring (=%s) with a number or another ring (=%s)!"%(x, y)
            raise TypeError, "unable to find a common parent for %s (parent: %s) and %s (parent: %s)"%(x,xp, y, yp)
    return x, y


def bin_op(x, y, op):
    if isinstance(y, element.InfinityElement):
        return op(y,x)
    if op == operator.mul and \
           isinstance(y, (\
                      element.ModuleElement,
                      element.AlgebraElement,
                      module.Module)) and \
           isinstance(x, (element.RingElement, int, long, float)):
        return op(y,x)
    try:
        x, y = canonical_coercion(x, y)
    except TypeError, mesg:
        try:
            return y._r_action(x)
        except AttributeError:
            raise TypeError, mesg
        except TypeError:
            raise TypeError, "No right action of %s on %s defined"%(y,x)
        try:
            return x._l_action(y)
        except AttributeError:
            raise TypeError, mesg
        except TypeError:
            raise TypeError, "No left action of %s on %s defined"%(x,y)

    return op(x,y)



N = type(None)
def cmp(x, y):
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
            fails = fails + 1

        try:
            y0 = y
            y = coerce(xp, y)
        except (TypeError, ValueError):
            fails = fails + 1

        if fails == 0:
            c0 = __builtin__.cmp(x0,y)
            c1 = __builtin__.cmp(x,y0)
            if c0 == c1:
                return c0
            else:
                return -1

        elif fails == 2:

            return -1

        else:
            if not (parent(x) is parent(y)):
                raise RuntimeError, "There is a bug in coercion: x=%s (parent=%s), y=%s (parent=%s)"%(x, parent(x), y, parent(y))
            return __builtin__.cmp(x,y)
