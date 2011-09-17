"""
Coercion helper functions
(This is the pure-python version, which is not used.)
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

import operator

import sage.structure.element
import sage.modules.module
import sage.rings.ring

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
    xp = parent(x)
    yp = parent(y)
    if isinstance(x, (int, long)):
        x = yp(x)
        #x = coerce(yp, x)
    elif isinstance(y, (int, long)):
        y = xp(y)
        #y = coerce(xp, y)
    else:
        i = 0
        try:
            x = coerce(yp, x)
        except TypeError, ValueError:
            i += 1
        try:
            y = coerce(xp, y)
        except TypeError, ValueError:
            i += 1
        if i == 0:
            raise TypeError, "unable to find an unambiguous parent for %s (parent: %s) and %s (parent: %s)"%(x,xp, y, yp)
        elif i == 2:
            if isinstance(x, sage.rings.ring.Ring) or isinstance(y, sage.rings.ring.Ring):
                raise TypeError, "you cannot add ring to a number or to another ring!"
            raise TypeError, "unable to find any common parent for %s (parent: %s) and %s (parent: %s)"%(x,xp, y, yp)
    return x, y

def bin_op(x, y, op):
    if op == operator.mul and \
           isinstance(y, (\
                      sage.structure.element.ModuleElement,
                      sage.structure.element.AlgebraElement,
                      sage.modules.module.Module,
                      sage.structure.element.InfinityElement)) and \
       isinstance(x, (sage.structure.element.RingElement, int, long, float)):
        return op(y,x)
    x, y = canonical_coercion(x, y)
    return op(x,y)


def cmp(x, y):
    if type(x) == type(None) and type(y) != type(None):
        return -1
    elif type(x) != type(None) and type(y) == type(None):
        return -1
    elif isinstance(y, sage.structure.element.InfinityElement):
        return -y.__cmp__(x)
    x, y = canonical_coercion(x, y)
    return x.__cmp__(y)
