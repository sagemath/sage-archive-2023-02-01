from __future__ import print_function, absolute_import, division

import re
from copy import copy
from sage.rings.infinity import infinity
from sage.structure.factorization import Factorization
from sage.structure.element import Element, RingElement, AlgebraElement, ModuleElement
from sage.structure.parent import Parent
from sage.structure.parent_gens import ParentWithGens
from sage.misc.abstract_method import abstract_method
from sage.categories.homset import Hom
from sage.categories.fields import Fields
from sage.rings.integer import Integer
#from cpython.object import PyObject_RichCompare
from sage.categories.map import Map
from sage.rings.morphism import RingHomomorphism
from sage.rings.polynomial.polynomial_element import _dict_to_list
from sage.structure.element import coerce_binop
from sage.misc.superseded import experimental
from sage.functions.other import binomial

class OrePolynomial(AlgebraElement):
    def __init__(self, parent, x = None, check = 1, construct = False, **kwds):
        AlgebraElement.__init__(self, parent)
        self._parent = parent
        if x is None:
            self._coeffs = []
            return
        R = parent.base_ring()
        if isinstance(x, list):
            if check:
                self._coeffs = [R(t) for t in x]
                self.__normalize()
            else:
                self._coeffs = x
            return
        if isinstance(x, OrePolynomial):
            if x._parent is self._parent:
                x = list(x.list())
            elif R.has_coerce_map_from(x._parent):
                try:
                    if x.is_zero():
                        self._coeffs = []
                        return
                except (AttributeError, TypeError):
                    pass
                x = [x]
            else:
                self._coeffs = [R(a, **kwds) for a in x.list()]
                if check:
                    self.__normalize()
                return
        elif isinstance(x, int) and x == 0:
            self._coeffs = []
            return
        elif isinstance(x, dict):
            x = _dict_to_list(x, R.zero())
        elif not isinstance(x, list):
            x = [x]
        if check:
            self._coeffs = [R(z, **kwds) for z in x]
            self.__normalize()
        else:
            self._coeffs = x


    def _repr_(self, name = None):
        s = " "
        m = self.degree() + 1
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        coeffs = self.list()
        for n in reversed(range(m)):
            x = coeffs[n]
            if x:
                if n != m-1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        s = s.replace(" + -", " - ")
        s = re.sub(r' 1(\.0+)?\*',' ', s)
        s = re.sub(r' -1(\.0+)?\*',' -', s)
        if s == " ":
            return "0"
        return s[1:]

    
    def degree(self):
        return len(self._coeffs) -1


    def list(self, copy = True):
        if copy:
            return list(self._coeffs)
        else:
            return self._coeffs
    
    def __normalize(self):
        x = self._coeffs
        n = len(x) -1
        while n >= 0 and not x[n]:
            del x[n]
            n-= 1

    def _add_(self, other):
        x = self._coeffs
        y = other._coeffs
        dx = len(x)
        dy = len(y)
        if dx > dy:
            r = [x[i] + y[i] for i in range(dy)] + x[dy:]
        elif dy > dx:
            r = [x[i] + y[i] for i in range(dx)] + y[dx:]
        else:
            r = [x[i] + y[i] for i in range(dx)]
        return OrePolynomial(self._parent, r)

    def _sub_(self, other):
        x = self._coeffs
        y = other._coeffs
        dx = len(x)
        dy = len(y)
        if dx > dy:
            r = [x[i] - y[i] for i in range(dy)] + x[dy:]
        elif dy > dx:
            temp = [-t for t in y[dx:]]
            r = [x[i] - y[i] for i in range(dx)] + temp
        else:
            r = [x[i] - y[i] for i in range(dx)]
        return OrePolynomial(self._parent, r)

    def px_n(self, n):
        return OrePolynomial(self._parent, [0 for _ in range(n)] + self._coeffs) 
        
    def cx_n_a(self, n, a, c=1):
        if n == 0:
            return OrePolynomial(self._parent, [c*a])
        res = OrePolynomial(self._parent, [c*a]).px_n(n)
        for i in range(1, n+1):
            a = self._parent._derivation(a)
            temp = OrePolynomial(self._parent, [c*binomial(n,i)*a]).px_n(n-i)
            res = res + temp
        return res
        
    def cx_n_p(self, n, p, c = 1):
        coeffs = p._coeffs
        res = self._parent.zero()
        for i in range(p.degree()+1):
            temp = p.cx_n_a(n+i, coeffs[i], c)
            res = res + temp
        return res


    def is_zero(self):
        if self._coeffs == [] or self.degree() == -1:
            return True
        else: return False


        
    def _mul_(self, other):
        if self.degree() == 0:
            return OrePolynomial(self._parent, [self._coeffs[0]*x for x in other._coeffs])
        x, y = self._coeffs, other._coeffs
        dx, dy = self.degree(), other.degree()
        parent = self._parent
        if dx == -1 or dy == -1:
            return parent.zero()
        if self._parent._derivation.is_zero():
            coeffs = []
            for k in range(dx+dy+1):
                start = 0 if k <= dy else k - dy
                end = k if k <= dx else dx
                sum = x[start] * parent.twist_map(start)(y[k-start])
                for i in range(start, end+1):
                    sum += x[i]*parent.twist_map(i)(y[k-i])
                coeffs.append(sum)
            return OrePolynomial(parent, coeffs)
        elif self._parent._map.is_identity():
            res = self.cx_n_p(dx, other, x[dx])
            for i in range(dx-1, -1, -1):
                temp = self.cx_n_p(i, other, x[i])
                res = res + temp
            return res
        else:
            raise NotImplementedError()


    def div_euc(self, other):
        
        x = self._coeffs
        y = other._coeffs
        dx = self.degree()
        dy = other.degree()
        if dy > dx:
            return self._parent.zero(), self
        c = y[dy]
        res = self._parent.zero()
        i = 0
        while i < dx - dy + 1:
            if not c.divides(x[dx-i]):
                raise ValueError()
            factor = x[dx-i]/c
            temp = OrePolynomial(self._parent, [0 for _ in range(dx - dy)] + [factor])
            res = res + temp
            temp = temp*other
            reste = self - temp
            if reste.degree() < 1:
                break
            x = reste._coeffs
            i = len(x)
            
        return res, reste


    def rgcd(self, other):
        R = [self, other]
        i = 1
        while not R[i].is_zero():
            Q, temp = R[i-1].div_euc(R[i])
            R += [temp]
            i += 1
        temp = R[i-1].list()
        n = len(temp)
        temp = [x/temp[n] for x in temp]
        return OrePolynomial(self._parent, temp)
        



class OrePolynomial1(AlgebraElement):        
    def __init__(self, parent, x = None, check = 1, construct = False, **kwds):
        AlgebraElement.__init__(self, parent)
        self._parent = parent
        if x is None:
            self._coeffs = []
            return
        R = parent.base_ring()
        if isinstance(x, list):
            if check:
                self._coeffs = [R(t) for t in x]
                self.__normalize()
            else:
                self._coeffs = x
            return
        if isinstance(x, OrePolynomial):
            if x._parent is self._parent:
                x = list(x.list())
            elif R.has_coerce_map_from(x._parent):
                try:
                    if x.is_zero():
                        self._coeffs = []
                        return
                except (AttributeError, TypeError):
                    pass
                x = [x]
            else:
                self._coeffs = [R(a, **kwds) for a in x.list()]
                if check:
                    self.__normalize()
                return
        elif isinstance(x, int) and x == 0:
            self._coeffs = []
            return
        elif isinstance(x, dict):
            x = _dict_to_list(x, R.zero())
        elif not isinstance(x, list):
            x = [x]
        if check:
            self._coeffs = [R(z, **kwds) for z in x]
            self.__normalize()
        else:
            self._coeffs = x

    def _repr_(self, name = None):
        s = " "
        m = self.degree() + 1
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        coeffs = self.list()
        for n in reversed(range(m)):
            x = coeffs[n]
            if x:
                if n != m-1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        s = s.replace(" + -", " - ")
        s = re.sub(r' 1(\.0+)?\*',' ', s)
        s = re.sub(r' -1(\.0+)?\*',' -', s)
        if s == " ":
            return "0"
        return s[1:]

    def __call__(self, x):
        if self.is_zero():
            return 0
        res = x*self._coeffs[0]
        temp = x
        n = self.degree()+1
        for i in range(1, n):
            temp = self._parent._derivation(temp)
            res += temp*self._coeffs[i]
        return res
    
    def degree(self):
        return len(self._coeffs) -1

    def list(self, copy = True):
        if copy:
            return list(self._coeffs)
        else:
            return self._coeffs

    def __normalize(self):
        x = self._coeffs
        n = len(x) -1
        while n >= 0 and not x[n]:
            del x[n]
            n-= 1

    def _add_(self, other):
        x = self._coeffs
        y = other._coeffs
        dx = len(x)
        dy = len(y)
        if dx > dy:
            r = [x[i] + y[i] for i in range(dy)] + x[dy:]
        elif dy > dx:
            r = [x[i] + y[i] for i in range(dx)] + y[dx:]
        else:
            r = [x[i] + y[i] for i in range(dx)]
        return OrePolynomial1(self._parent, r)

    def _sub_(self, other):
        x = self._coeffs
        y = other._coeffs
        dx = len(x)
        dy = len(y)
        if dx > dy:
            r = [x[i] - y[i] for i in range(dy)] + x[dy:]
        elif dy > dx:
            temp = [-t for t in y[dx:]]
            r = [x[i] - y[i] for i in range(dx)] + temp
        else:
            r = [x[i] - y[i] for i in range(dx)]
        return OrePolynomial1(self._parent, r)

    def is_zero(self):
        if self._coeffs == [] or self.degree() == -1:
            return True
        else: return False

    def dom(self):
        temp = self.list()
        return temp[self.degree()]
    
    def px_n(self, n):
        return OrePolynomial1(self._parent, [0 for _ in range(n)] + self._coeffs) 

    def cx_n_a(self, n, a, c=1):
        if n == 0:
            return OrePolynomial1(self._parent, [c*a])
        res = OrePolynomial1(self._parent, [c*a]).px_n(n)
        for i in range(1, n+1):
            a = self._parent._derivation(a)
            temp = OrePolynomial1(self._parent, [c*binomial(n,i)*a]).px_n(n-i)
            res = res + temp
        return res
        
    def cx_n_p(self, n, p, c = 1):
        coeffs = p._coeffs
        res = OrePolynomial1(self._parent, [0])
        for i in range(len(coeffs)):
            temp = p.cx_n_a(n, coeffs[i], c)
            temp = temp.px_n(i)
            res = res + temp
        return res
    
    def _mul_(self, other):
        if self.is_zero() or other.is_zero():
            return OrePolynomial1(self._parent, [0])
        if self.degree() == 0:
            return OrePolynomial1(self._parent, [self._coeffs[0]*x for x in other._coeffs])
        x, y = self._coeffs, other._coeffs
        dx, dy = self.degree(), other.degree()
        res = OrePolynomial1(self._parent, [0])
        for i in range(dx+1):
            temp = self.cx_n_p(i, other, x[i])
            res += temp
        return res

    def div_euc(self, other):
        q = OrePolynomial1(self._parent, [0])
        r = self
        if other.is_zero():
            raise ZeroDivisionError
        while not r.is_zero() and r.degree() >= other.degree():
            tempr = r.list()
            tempother = other.list()
            dr = r.degree()
            dother = other.degree()
            tempt = tempr[dr]/tempother[dother]
            t = OrePolynomial1(self._parent, [0 for _ in range(dr - dother)] + [tempt])
            q = q+t
            r = r-t*other
        return q, r
    


    def rgcd(self, other):
        R = [self, other]
        i = 1
        while not R[i].is_zero():
            Q, temp = R[i-1].div_euc(R[i])
            R += [temp]
            i += 1
        temp = R[i-1].list()
        n = len(temp)
        temp = [x/temp[n-1] for x in temp]
        return OrePolynomial1(self._parent, temp)

    def llcm(self, other):
        R = [self, other]
        U = [OrePolynomial1(self._parent, [1]), OrePolynomial1(self._parent, [0])]
        V = [OrePolynomial1(self._parent, [0]), OrePolynomial1(self._parent, [1])]
        i = 1
        while not R[i].is_zero():
            Q, temp = R[i-1].div_euc(R[i])
            R += [temp]
            U += [U[i-1] - Q*U[i]]
            V += [V[i-1] - Q*V[i]]
            i += 1
        tempU = (U[i]*self).list()
        n = tempU[len(tempU)-1]
        U = [x/n for x in tempU]
        U = OrePolynomial1(self._parent, U)
        return U
        

    def bezout(self, other):
        R = [self, other]
        U = [OrePolynomial1(self._parent, [1]), OrePolynomial1(self._parent, [0])]
        V = [OrePolynomial1(self._parent, [0]), OrePolynomial1(self._parent, [1])]
        i = 1
        while not R[i].is_zero():
            Q, temp = R[i-1].div_euc(R[i])
            R += [temp]
            U += [U[i-1] - Q*U[i]]
            V += [V[i-1] - Q*V[i]]
            i += 1
        tempR = R[i-1].list()
        n = len(tempR)
        R = [x/tempR[n-1] for x in tempR]
        tempU = U[i-1].list()
        U = [x/tempR[n-1] for x in tempU]
        tempV = V[i-1].list()
        V = [x/tempR[n-1] for x in tempV]

        return OrePolynomial1(self._parent, U), OrePolynomial1(self._parent, V), OrePolynomial1(self._parent, R)
