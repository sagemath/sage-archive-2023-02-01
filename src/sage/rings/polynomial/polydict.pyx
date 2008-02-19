"""
PolyDict functionality -- an implementation of the underlying
arithmetic for multi-variate polynomial rings using Python dicts.

This class is not meant for end users, but instead for implementing
multivariate polynomial rings over a completely general base.  It does
not do strong type checking or have parents, etc.  For more special
bases, etc., one would implement something similar in Pyrex for speed.

AUTHOR: William Stein, David Joyner and Martin Albrecht (ETuple)
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein (wstein@ucsd.edu)
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


import copy
from sage.structure.element import generic_power
from sage.misc.misc import cputime
from sage.misc.latex import latex

"""
SUMMARY:
    The functions in this file use the 'dictionary representation'
    of multivariate polynomials

        {(e1,...,er):c1,...} <-> c1*x1^e1*...*xr^er+...,

    which we call a polydict. The exponent tuple (e1,...,er) in this
    representation is an instance of the class ETuple. This class
    behaves like a normal Python tuple but also offers advanced access
    methods for sparse monomials like positions of non-zero exponents
    etc.

    This file implements arithmetic functionality for polydicts.



AUTHORS:  William Stein and Martin Albrecht (ETuples)
"""


import sage.rings.ring_element as ring_element

cdef class PolyDict:
    cdef object __repn
    cdef object __zero

    def __init__(PolyDict self, pdict, zero=0, remove_zero=False, force_int_exponents=True, force_etuples=True):
        """
        INPUT:
            pdict -- list, which represents a multi-variable polynomial
                     with the distribute representation (a copy is not made)

            zero --  (optional) zero in the base ring

            force_int_exponents -- bool (optional) arithmetic with int exponents is much
                      faster than some of the alternatives, so this is True by default.

            force_etuples -- bool (optional) enforce that the exponent tuples are instances
                             of ETuple class

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PolyDict({(2,3):2, (1,2):3, (2,1):4})
            PolyDict with representation  {(1, 2): 3, (2, 1): 4, (2, 3): 2}

            sage: PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1):4}, force_int_exponents=False)
            PolyDict with representation {(2, 1): 4, (1, 2, 1): 3, (2/3, 3, 5): 2}

            sage: PolyDict({(2,3):0, (1,2):3, (2,1):4}, remove_zero=True)
            PolyDict with representation {(1, 2): 3, (2, 1): 4}
        """
        if not isinstance(pdict, dict):
            if isinstance(pdict, list):
                v = {}
                for w in pdict:
                    if w[0] != 0:
                        v[ETuple(w[1])] = w[0]
                remove_zero = False
                pdict = v
            else:
                raise TypeError, "pdict must be a list."

        if isinstance(pdict,dict) and force_etuples==True:
            pdict2 = []
            for k,v in pdict.iteritems():
                pdict2.append((ETuple(k),v))

            pdict = dict(pdict2)

        if force_int_exponents:
            new_pdict = {}
            if remove_zero:
                for k, c in pdict.iteritems():
                    if c != zero:
                        new_pdict[ETuple(map(int,k))] = c
            else:
                for k, c in pdict.iteritems():
                    new_pdict[ETuple(map(int,k))] = c
            pdict = new_pdict
        else:
            if remove_zero:
                for k in pdict.keys():
                    if pdict[k] == zero:
                        del pdict[k]
        self.__repn  = pdict
        self.__zero = zero

    def __cmp__(self, PolyDict right):
        return cmp(self.__repn, right.__repn)

    def compare(PolyDict self, PolyDict other, fn=None):
        if fn is None:
            return cmp(self.__repn, other.__repn)

        left  = iter(sorted( self.__repn,fn,reverse=True)) #start with biggest
        right = iter(sorted(other.__repn,fn,reverse=True))

        for m in left:
            try:
                n = right.next()
            except StopIteration:
                return 1 # left has terms, right doesn't
            ret =  fn(m,n)
            if ret!=0:
                return ret # we have a difference
            ret = cmp(self.__repn[m],other.__repn[n]) #compare coefficents
            if ret!=0:
                return ret #if they differ use it
            #try next pair

        try:
            n = right.next()
        except StopIteration:
            return 0 # both have no terms

        return -1 # right has terms, left doesn't

    def list(PolyDict self):
        """
        Return a list that defines self. It is safe to change this.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.list()
            [[3, [1, 2]], [4, [2, 1]], [2, [2, 3]]]
        """
        ret = []
        for e,c in self.__repn.iteritems():
            ret.append([c,list(e)])
        return ret

    def dict(PolyDict self):
        """
        Return a copy of the dict that defines self.  It is
        safe to change this.  For a reference, use dictref.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.dict()
            {(1, 2): 3, (2, 1): 4, (2, 3): 2}
        """
        return self.__repn.copy()

    def coefficients(PolyDict self):
        """
        Return the coefficients of self.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.coefficients()
            [3, 4, 2]
        """
        return self.__repn.values()

    def exponents(PolyDict self):
        """
        Return the exponents of self.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.exponents()
            [(1, 2), (2, 1), (2, 3)]
        """
        return self.__repn.keys()

    def __len__(PolyDict self):
        """
        Return the number of terms of the polynomial.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: len(f)
            3
        """
        return len(self.__repn)

    def __getitem__(PolyDict self, e):
        """
        Return a coefficient of the polynomial.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f[1,2]
            3
            sage: f[(2,1)]
            4
        """
        return self.__repn[ETuple(e)]

    def __repr__(PolyDict self):
        return 'PolyDict with representation %s'%self.__repn

    def degree(PolyDict self, PolyDict x=None):
        if x is None:
            return self.total_degree()
        L = x.__repn.keys()
        if len(L) != 1:
            raise TypeError, "x must be one of the generators of the parent."
        L = L[0]
        nonzero_positions = L.nonzero_positions()
        if len(nonzero_positions) != 1:
            raise TypeError, "x must be one of the generators of the parent."
        i = nonzero_positions[0]
        if L[i] != 1:
            raise TypeError, "x must be one of the generators of the parent."
        _max = []
        for v in self.__repn.keys():
            _max.append(v[i])
        return max(_max)

    def total_degree(PolyDict self):
        return max([-1] + map(sum,self.__repn.keys()))

    def monomial_coefficient(PolyDict self, mon):
        K = mon.keys()[0]
        if not self.__repn.has_key(K):
            return 0
        return self.__repn[K]

    def coefficient(PolyDict self, mon):
        """
        Return a polydict that defines a polynomial in 1 less number
        of variables that gives the coefficient of mon in this
        polynomial.

        The coefficient is defined as follows.  If f is this
        polynomial, then the coefficient is the sum T/mon where the
        sum is over terms T in f that are exactly divisible by mon.
        """
        K = mon.keys()[0]
        nz = K.nonzero_positions() #set([i for i in range(len(K)) if K[i] != 0])
        ans = {}
        for S in self.__repn.keys():
            exactly_divides = True
            for j in nz:
                if S[j] != K[j]:
                    exactly_divides = False
                    break
            if exactly_divides:
                t = list(S)
                for m in nz:
                    t[m] = 0
                ans[ETuple(t)] = self.__repn[S]
        return PolyDict(ans, force_etuples=False)



    def is_homogeneous(PolyDict self):
        K = self.__repn.keys()
        if len(K) == 0:
            return True
        # A polynomial is homogeneous if the number of different
        # exponent sums is at most 1.
        return len(set(map(sum,K))) <= 1

    def homogenize(PolyDict self):
        R = self.__repn
        H = {}
        deg = self.degree()
        for e, val in R.iteritems():
            i = deg - sum(e)
            f = tuple(e) + (i,)
            H[ETuple(f)] = val
        return PolyDict(H, zero=self.__zero, force_etuples=False)

    def latex(PolyDict self, vars, atomic_exponents=True, atomic_coefficients=True, cmpfn = None):
        """
        Return a nice polynomial latex representation of this PolyDict, where
        the vars are substituted in.

        INPUT:
            vars -- list
            atomic_exponents -- bool (default: True)
            atomic_coefficients -- bool (default: True)

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.latex(['a','WW'])
            '2 a^{2}WW^{3} + 4 a^{2}WW + 3 aWW^{2}'

        When atomic_exponents is False, the exponents are surrounded
        in parenthesis, since ^ has such high precedence.
            sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            sage: f.latex(['a','b','c'], atomic_exponents=False)
            '4 a^{2}bc + 3 ab^{2}c + 2 a^{2/3}b^{3}c^{5}'
        """
        n = len(vars)
        poly = ""
        E = self.__repn.keys()
        if cmpfn:
            E.sort(cmp = cmpfn, reverse=True)
        else:
            E.sort(reverse=True)
        for e in E:
            c = self.__repn[e]
            if c != 0:
                sign_switch = False
                # First determine the multinomial:
                multi = ""
                for j in e.nonzero_positions(sort=True):
                    multi = multi + vars[j]
                    if e[j] != 1:
                        multi =  multi +"^{%s}"%e[j]
                # Next determine coefficient of multinomial
                if len(multi) == 0:
                    multi = latex(c)
                elif c != 1:
                    if not atomic_coefficients:
                        c = latex(c)
                        if c.find("+") != -1 or c.find("-") != -1 or c.find(" ") != -1:
                            c = "(%s)"%c
                    if len(poly) > 0 and c == -1:
                        sign_switch = True
                    else:
                        multi = "%s %s"%(latex(c),multi)

                # Now add on coefficiented multinomials
                if len(poly) > 0:
                    if sign_switch:
                        poly = poly + " - "
                    else:
                        poly = poly + " + "
                poly = poly + multi
        poly = poly.lstrip().rstrip()
        poly = poly.replace("+ -","- ")
        if len(poly) == 0:
            return "0"
        return poly


    def poly_repr(PolyDict self, vars, atomic_exponents=True, atomic_coefficients=True, cmpfn = None):
        """
        Return a nice polynomial string representation of this PolyDict, where
        the vars are substituted in.

        INPUT:
            vars -- list
            atomic_exponents -- bool (default: True)
            atomic_coefficients -- bool (default: True)

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.poly_repr(['a','WW'])
            '2*a^2*WW^3 + 4*a^2*WW + 3*a*WW^2'

        When atomic_exponents is False, the exponents are surrounded
        in parenthesis, since ^ has such high precedence.
            sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            sage: f.poly_repr(['a','b','c'], atomic_exponents=False)
            '4*a^(2)*b*c + 3*a*b^(2)*c + 2*a^(2/3)*b^(3)*c^(5)'
        """
        n = len(vars)
        poly = ""
        E = self.__repn.keys()
        if cmpfn:
            E.sort(cmp = cmpfn, reverse=True)
        else:
            E.sort(reverse=True)
        for e in E:
            c = self.__repn[e]
            if c != 0:
                sign_switch = False
                # First determine the multinomial:
                multi = ""
                for j in e.nonzero_positions(sort=True):
                    if len(multi) > 0:
                        multi = multi + "*"
                    multi = multi + vars[j]
                    if e[j] != 1:
                        if atomic_exponents:
                            multi = multi + "^%s"%e[j]
                        else:
                            multi = multi + "^(%s)"%e[j]
                # Next determine coefficient of multinomial
                if len(multi) == 0:
                    multi = str(c)
                elif c != 1:
                    if not atomic_coefficients:
                        c = str(c)
                        if c.find("+") != -1 or c.find("-") != -1 or c.find(" ") != -1:
                            c = "(%s)"%c
                    if len(poly) > 0 and c == -1:
                        sign_switch = True
                    else:
                        multi = "%s*%s"%(c,multi)

                # Now add on coefficiented multinomials
                if len(poly) > 0:
                    if sign_switch:
                        poly = poly + " - "
                    else:
                        poly = poly + " + "
                poly = poly + multi
        poly = poly.lstrip().rstrip()
        poly = poly.replace("+ -","- ")
        if len(poly) == 0:
            return "0"
        return poly


    def __add__(PolyDict self,PolyDict other):
        """
        Add two PolyDict's in the same number of variables.

        EXAMPLES:
        We add two polynomials in 2 variables:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(1,5):-3, (2,3):-2, (1,1):3})
            sage: f+g
            PolyDict with representation  {(1, 1): 3, (1, 2): 3, (2, 1): 4, (1, 5): -3}

        Next we add two polynomials with fractional exponents in 3 variables:
            sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            sage: g = PolyDict({(2/3,3,5):3}, force_int_exponents=False)
            sage: f+g
            PolyDict with representation {(1, 2, 1): 3, (2/3, 3, 5): 5, (2, 1, 1): 4}
        """
        zero = self.__zero
        #D = copy.copy(self.__repn)
        D = self.__repn.copy() #faster!
        R = other.__repn
        for e,c in R.iteritems():
            try:
                D[e] += c
            except KeyError:
                D[e] = c
        F = PolyDict(D, zero=zero, remove_zero=True, force_int_exponents=False,  force_etuples=False)
        return F

    def __mul__(PolyDict self,PolyDict  right):
        """
        Multiply two PolyDict's in the same number of variables.

        EXAMPLES:
        We multiply two polynomials in 2 variables:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(1,5):-3, (2,3):-2, (1,1):3})
            sage: f*g
            PolyDict with representation {(3, 4): 6, (3, 5): -6, (4, 4): -8, (3, 2): 12, (4, 6): -4, (2, 3): 9, (2, 7): -9, (3, 8): -6, (3, 6): -12}

        Next we multiply two polynomials with fractional exponents in 3 variables:
            sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            sage: g = PolyDict({(2/3,3,5):3}, force_int_exponents=False)
            sage: f*g
            PolyDict with representation {(8/3, 4, 6): 12, (5/3, 5, 6): 9, (4/3, 6, 10): 6}

        Finally we print the result in a nice format.
            sage: (f*g).poly_repr(['a','b','c'], atomic_exponents = False)
            '12*a^(8/3)*b^(4)*c^(6) + 9*a^(5/3)*b^(5)*c^(6) + 6*a^(4/3)*b^(6)*c^(10)'
        """
        newpoly = {}
        k = self.__repn.keys()
        if len(k) == 0:   # product is zero anyways
            return self
        n = len(self.__repn.keys()[0])
        #r = range(n)
        for e0, c0 in self.__repn.iteritems():
            for e1, c1 in right.__repn.iteritems():
                #r = e1.nonzero_positions().union(e2.nonzero_positions())
                #e = etuple([(i,e0[i] + e1[i]) for i in r], n )
                e = e0.eadd(e1)
                c = c0*c1
                if e in newpoly:
                    newpoly[e] = newpoly[e] + c
                else:
                    newpoly[e] = c
        F = PolyDict(newpoly, self.__zero, force_int_exponents=False, remove_zero=True, force_etuples=False)
        return F

    def scalar_rmult(PolyDict self, s):
        """
        Right Scalar Multiplication

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: x,y=FreeMonoid(2,'x,y').gens()  # a strange object to live in a polydict, but non-commutative!
            sage: f = PolyDict({(2,3):x})
            sage: f.scalar_rmult(y)
            PolyDict with representation {(2, 3): x*y}
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.scalar_rmult(-2)
            PolyDict with representation {(1, 2): -6, (2, 1): -8, (2, 3): -4}
        """
        v = {}
        # if s is 0, then all the products will be zero
        if s != self.__zero:
            for e, c in self.__repn.iteritems():
                v[e] = c*s
        return PolyDict(v, self.__zero, force_int_exponents=False, force_etuples=False)

    def scalar_lmult(PolyDict self, s):
        """
        Left Scalar Multiplication

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: x,y=FreeMonoid(2,'x,y').gens()  # a strange object to live in a polydict, but non-commutative!
            sage: f = PolyDict({(2,3):x})
            sage: f.scalar_lmult(y)
            PolyDict with representation {(2, 3): y*x}
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.scalar_lmult(-2)
            PolyDict with representation {(1, 2): -6, (2, 1): -8, (2, 3): -4}
        """
        v = {}
        # if s is 0, then all the products will be zero
        if s != self.__zero:
            for e, c in self.__repn.iteritems():
                v[e] = s*c
        return PolyDict(v, self.__zero, force_int_exponents=False, force_etuples=False)

    def __sub__(PolyDict self,PolyDict  other):
        """
        Subtract two PolyDict's.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(2,3):2, (1,1):-10})
            sage: f - g
            PolyDict with representation {(1, 1): 10, (1, 2): 3, (2, 1): 4}
            sage: g - f
            PolyDict with representation {(1, 1): -10, (1, 2): -3, (2, 1): -4}
        """

        # TOOD: should refactor add, make abstract operator, so can do both +/-; or copy code.
        return self + other.scalar_lmult(-1)

    def __one(PolyDict self):
        one = self.__zero + 1
        if len(self.__repn) == 0:
            v = {(0):one}
        else:
            v = {ETuple({},len(self.__repn.keys()[0])):one}
        return PolyDict(v, self.__zero, force_int_exponents=False, force_etuples=False)

    def __pow__(PolyDict self, n, ignored):
        """
        Return the n-th nonnegative power of this PolyDict.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f**2
            PolyDict with representation {(3, 3): 24, (3, 5): 12, (4, 4): 16, (4, 2): 16, (4, 6): 4, (2, 4): 9}
            sage: f**0
            PolyDict with representation {(0, 0): 1}
            sage: (f-f)**0
            Traceback (most recent call last):
            ...
            ArithmeticError: 0^0 is undefined.
        """
        return generic_power(self, n, self.__one())

    def lcmt(PolyDict self, greater_etuple):
        """
        Provides functionality of lc, lm, and lt by calling the tuple
        compare function on the provided term order T.

        INPUT:
            greater_etuple -- a term order
        """
        try:
            return ETuple(reduce(greater_etuple,self.__repn.keys()))
        except KeyError:
            raise ArithmeticError, "%s not supported",greater_etuple

    def __reduce__(PolyDict self):
        """

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: loads(dumps(f)) == f
            True
        """
        return make_PolyDict,(self.__repn,)


cdef class ETupleIter:
    cdef int _i
    cdef int _length
    cdef object _data

    def __init__(self, data, length):
        self._data = data
        self._length = length
        self._i = -1

    def __next__(self):
        self._i = self._i + 1

        if self._i == self._length:
            self._i = -1
            raise StopIteration

        return self._data.get(self._i,0)

cdef class ETuple:
    """
    Representation of the exponents of a polydict monomial. If
    (0,0,3,0,5) is the exponent tuple of x_2^3*x_4^5 then this class
    only stores {2:3,4:5} instead of the full tuple. This sparse
    information may be optained by provided methods.
    """
    cdef object _data
    cdef object _length

    def __init__(ETuple self, data, length=None):
        """
        ETuple() -> an empty ETuple
        ETuple(sequence) -> ETuple initialized from sequence's items

        If the argument is an ETuple, the return value is the same object.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1,1,0])
            (1, 1, 0)
            sage: ETuple({int(1):int(2)},int(3))
            (0, 2, 0)
        """
        if isinstance(data,ETuple):
            self._data = (<ETuple>data)._data
            self._length = (<ETuple>data)._length
        elif isinstance(data,dict) and isinstance(length,int):
            self._data = data
            self._length = length
        elif isinstance(data,(tuple,list)):
            tpl = zip(range(len(data)),data)
            tmp = []
            for (i,v) in tpl:
                if v:
                    tmp.append((i,v))
            self._data = dict( tmp )
            self._length = len(data)
        else:
            raise TypeError

    # methods to simulate tuple

    def __add__(ETuple self,ETuple other):
        """
        x.__add__(n) <==> x+n

        concatenates two ETuples

        EXAMPLE:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1,1,0]) + ETuple({int(1):int(2)},int(3))
            (1, 1, 0, 0, 2, 0)

        """
        data = self._data.copy()

        for i,v in other._data.iteritems():
            data[i+self._length]=v

        return ETuple(data,self._length+other._length)

    def __rmul__(ETuple self, factor):
        """
        x.__rmul__(n) <==> n*x
        """
        if factor <= 0:
            return ETuple({},0)
        d = {}
        for k,v in self._data.iteritems():
            for i in range(0,factor):
                d[(i*self._length)+k]=v
        return ETuple(d,int(self._length*factor))

    def __mul__(ETuple self,factor):
        """
        x.__mul__(n) <==> x*n

        EXAMPLE:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1,2,3])*2
            (1, 2, 3, 1, 2, 3)
        """
        if factor <= 0:
            return ETuple({},0)
        d = {}
        for k,v in self._data.iteritems():
            for i in range(0,factor):
                d[(i*self._length)+k]=v
        return ETuple(d,int(self._length*factor))

    def __getitem__(ETuple self,int i):
        """
        x.__getitem__(i) <==> x[i]
        """

        return self._data.get(i%self._length,0)

    def __getslice__(ETuple self, Py_ssize_t i, Py_ssize_t j):
        """
        x.__getslice(i,j) <==> x[i:j]

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e=ETuple([1,2,3])
            sage: e[1:]
            (2, 3)
            sage: e[:1]
            (1,)

        """
        if i<0:
            i = i % self._length
        elif i>self._length:
            i = self._length

        if j<0:
            j = j % self._length
        elif j>self._length:
            j = self._length

        d = {}
        for k,v in self._data.iteritems():
            if i<=k and k<j:
                d[k-i]=v
        return ETuple(d,j-i)

    def __hash__(ETuple self):
        """
        x.__hash__() <==> hash(x)
        """
        return hash((tuple(sorted(self._data.iteritems())),self._length))

    def __len__(ETuple self):
        """
        x.__len__() <==> len(x)
        """
        return self._length

    def __contains__(ETuple self, elem):
        """
        x.__contains__(n) <==> n in x


        EXAMPLE:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple({int(1):int(2)},int(3)); e
            (0, 2, 0)
            sage: 1 in e
            False
            sage: 2 in e
            True
        """
        if elem!=0:
            return elem in self._data.values()
        else:
            return len(self._data)!=self._length

    def __richcmp__(ETuple self, ETuple other, op):
        if op == 0:  #<

            for k in sorted(set(other._data.iterkeys()).union(self._data.iterkeys())):
                sk = self._data.get(k,0)
                ok = other._data.get(k,0)
                if sk > ok:
                    return False
                elif sk < ok:
                    return True

            if self._data==other._data:
                return False
            else:
                return True

        if op == 2: #==

            return self._data == other._data and self._length == other._length

        if op == 4: #>

            for k in sorted(set(other._data.iterkeys()).union(self._data.iterkeys())):
                sk = self._data.get(k,0)
                ok = other._data.get(k,0)
                if sk<ok:
                    return False
                elif sk>ok:
                    return True

            if self._data==other._data:
                return False
            else:
                return True


        if op == 1: #<=

            if self._data==other._data:
                return True

            for k in sorted(set(other._data.iterkeys()).union(self._data.iterkeys())):
                sk = self._data.get(k,0)
                ok = other._data.get(k,0)
                if sk > ok:
                    return False
                elif sk < ok:
                    return True

            return True # should never get here


        if op == 3: #!=

            return self._data != other._data

        if op == 5: #>=

            for k in sorted(set(other._data.iterkeys()).union(self._data.iterkeys())):
                sk = self._data.get(k,0)
                ok = other._data.get(k,0)
                if sk<ok:
                    return False
                elif sk>ok:
                    return True
            return True

    def __iter__(ETuple self):
        """
        x.__iter__() <==> iter(x)
        """
        return ETupleIter(self._data,self._length)

    def __str__(ETuple self):
        return self.__repr__()

    def __repr__(ETuple self):
        res = [0,]*self._length
        for i,v in self._data.iteritems():
            res[i]=v
        return str(tuple(res))

    def __reduce__(ETuple self):
        """
        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,1,0])
            sage: bool(e == loads(dumps(e)))
            True
        """

        return make_ETuple,(self._data,self._length)

    # additional methods

    def eadd(ETuple self,ETuple other):
        """
        Vector addition of self with other.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: f = ETuple([0,1,1])
            sage: e.eadd(f)
            (1, 1, 3)

        """
        if self._length!=other._length:
            raise ArithmeticError

        d = self._data.copy()
        for k,v in other._data.iteritems():
            if d.has_key(k):
                d[k] = d[k] + v
            else:
                d[k] = v
        return ETuple(d,self._length)

    def esub(ETuple self,ETuple other):
        """
        Vector subtraction of self with other.

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: f = ETuple([0,1,1])
            sage: e.esub(f)
            (1, -1, 1)
        """
        if self._length!=other._length:
            raise ArithmeticError

        d = self._data.copy()
        for k,v in other._data.iteritems():
            if d.has_key(k):
                if d[k] - v != 0:
                    d[k] = d[k] - v
                else:
                    del d[k]
            else:
                d[k] = -v
        return ETuple(d,self._length)

    def emul(ETuple self,int factor):
        """
        Scalar Vector multiplication of self.

        EXAMPLE:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: e.emul(2)
            (2, 0, 4)
        """
        d = {}
        for k,v in self._data.iteritems():
            d[k]=int(factor*v)
        return ETuple(d,self._length)

    def nonzero_positions(ETuple self,sort=False):
        """
        Returns the positions of non-zero exponents in the tuple.

        INPUT:
            sort -- if True a sorted list is returned. If False a an
                    unsorted list is returned. (default: False)

        EXAMPLES:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: e.nonzero_positions()
            [0, 2]
        """
        if sort:
            return sorted(self._data.iterkeys())
        else:
            return self._data.keys()

    def common_nonzero_positions(ETuple self, ETuple other,sort=False):
        """
        Returns an optionally sorted list of non zero positions either
        in self or other, i.e. the only positions that need to be
        considered for any vector operation.

        EXAMPLE:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: f = ETuple([0,0,1])
            sage: e.common_nonzero_positions(f)
            set([0, 2])
            sage: e.common_nonzero_positions(f,sort=True)
            [0, 2]
        """


        res = set(self._data.iterkeys()).union(other._data.iterkeys())
        if sort:
            return sorted(res)
        else:
            return res

    def nonzero_values(ETuple self, sort=True):
        """
        Returns the non-zero values of the tuple.

        INPUT:
            sort -- if True the values are sorted by their indices. Otherwise a
                    the values are returned unsorted. (default: True)

        EXAMPLE:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([2,0,1])
            sage: e.nonzero_values()
            [2, 1]
            sage: f = ETuple([0,-1,1])
            sage: f.nonzero_values(sort=True)
            [-1, 1]
        """
        if sort:
            tmp = []
            for e in self.nonzero_positions(True):
                tmp.append(self._data[e])
            return tmp
        else:
            return self._data.values()

    def reversed(ETuple self):
        """
        Returns the reversed ETuple of self.

        EXAMPLE:
            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,2,3])
            sage: e.reversed()
            (3, 2, 1)
        """
        data = {}
        length = self._length-1
        for k,v in self._data.iteritems():
            data[length-k]=v
        return ETuple(data,length+1)

    def sparse_iter(ETuple self):
        """
        Iterator over the elements of self where the elements are
        returned as $(i,e)$ where $i$ is the position of $e$ in the
        tuple.
        """
        return self._data.iteritems()

def make_PolyDict(data):
    return PolyDict(data,remove_zero=False, force_int_exponents=False, force_etuples=False)

def make_ETuple(data,length):
    return ETuple(data,length)
