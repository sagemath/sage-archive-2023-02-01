"""
PolyDict engine for generic multivariate polynomial rings

This module provides an implementation of the underlying arithmetic for
multi-variate polynomial rings using Python dicts.

This class is not meant for end users, but instead for implementing
multivariate polynomial rings over a completely general base.  It does
not do strong type checking or have parents, etc. For speed, it has been
implemented in Cython.

The functions in this file use the 'dictionary representation' of multivariate
polynomials

``{(e1,...,er):c1,...} <-> c1*x1^e1*...*xr^er+...``,

which we call a polydict. The exponent tuple ``(e1,...,er)`` in this
representation is an instance of the class :class:`ETuple`. This class behaves
like a normal Python tuple but also offers advanced access methods for sparse
monomials like positions of non-zero exponents etc.

AUTHORS:

- William Stein
- David Joyner
- Martin Albrecht (ETuple)
- Joel B. Mohler (2008-03-17) -- ETuple rewrite as sparse C array
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

include "sage/ext/stdsage.pxi"
include 'sage/ext/cdefs.pxi'
from cpython.dict cimport *

import copy
from functools import reduce
from sage.structure.element import generic_power
from sage.misc.misc import cputime
from sage.misc.latex import latex

import sage.rings.ring_element as ring_element

cdef class PolyDict:
    def __init__(PolyDict self, pdict, zero=0, remove_zero=False, force_int_exponents=True, force_etuples=True):
        """
        INPUT:

        - ``pdict`` -- list, which represents a multi-variable polynomial with
          the distribute representation (a copy is not made)

        - ``zero`` --  (optional) zero in the base ring

        - ``force_int_exponents`` -- bool (optional) arithmetic with int
          exponents is much faster than some of the alternatives, so this is
          True by default.

        - ``force_etuples`` -- bool (optional) enforce that the exponent tuples
          are instances of ETuple class

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: PolyDict({(2,3):2, (1,2):3, (2,1):4})
            PolyDict with representation {(1, 2): 3, (2, 3): 2, (2, 1): 4}

            # I've removed fractional exponent support in ETuple when moving to a sparse C integer array
            #sage: PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1):4}, force_int_exponents=False)
            #PolyDict with representation {(2, 1): 4, (1, 2, 1): 3, (2/3, 3, 5): 2}

            sage: PolyDict({(2,3):0, (1,2):3, (2,1):4}, remove_zero=True)
            PolyDict with representation {(1, 2): 3, (2, 1): 4}

            sage: PolyDict({(0,0):RIF(-1,1)}, remove_zero=True)
            PolyDict with representation {(0, 0): 0.?}
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
                raise TypeError("pdict must be a list.")

        if isinstance(pdict,dict) and force_etuples==True:
            pdict2 = []
            for k,v in pdict.iteritems():
                pdict2.append((ETuple(k),v))

            pdict = dict(pdict2)

        if force_int_exponents:
            new_pdict = {}
            if remove_zero:
                for k, c in pdict.iteritems():
                    if not c == zero:
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
            ret = cmp(self.__repn[m],other.__repn[n]) #compare coefficients
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

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.list()
            [[3, [1, 2]], [2, [2, 3]], [4, [2, 1]]]
        """
        ret = []
        for e,c in self.__repn.iteritems():
            ret.append([c,list(e)])
        return ret

    def dict(PolyDict self):
        """
        Return a copy of the dict that defines self.  It is
        safe to change this.  For a reference, use dictref.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.dict()
            {(1, 2): 3, (2, 3): 2, (2, 1): 4}
        """
        return self.__repn.copy()

    def coefficients(PolyDict self):
        """
        Return the coefficients of self.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.coefficients()
            [3, 2, 4]
        """
        return self.__repn.values()

    def exponents(PolyDict self):
        """
        Return the exponents of self.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.exponents()
            [(1, 2), (2, 3), (2, 1)]
        """
        return self.__repn.keys()

    def __len__(PolyDict self):
        """
        Return the number of terms of the polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: len(f)
            3
        """
        return len(self.__repn)

    def __getitem__(PolyDict self, e):
        """
        Return a coefficient of the polynomial.

        EXAMPLES::

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
            raise TypeError("x must be one of the generators of the parent.")
        L = L[0]
        nonzero_positions = L.nonzero_positions()
        if len(nonzero_positions) != 1:
            raise TypeError("x must be one of the generators of the parent.")
        i = nonzero_positions[0]
        if L[i] != 1:
            raise TypeError("x must be one of the generators of the parent.")
        _max = []
        for v in self.__repn.keys():
            _max.append(v[i])
        return max(_max)

    def valuation(PolyDict self, PolyDict x = None):
        L = x.__repn.keys()
        if x is None:
            _min = []
            negative = False
            for k in self.__repn.keys():
                _sum = 0
                for m in self.__repn[k].nonzero_values(sort=False):
                    if m < 0:
                        negative = True
                        break
                    _sum += m
                if negative:
                    break
                _min.append(_sum)
            else:
                return min(_min)
            for k in self.__repn.keys():
                _min.append(sum([m for m in self.__repn[k].nonzero_values(sort=False) if m < 0]))
            return min(_min)
        L = x.__repn.keys()
        if len(L) != 1:
            raise TypeError("x must be one of the generators of the parent.")
        L = L[0]
        nonzero_positions = L.nonzero_positions()
        if len(nonzero_positions) != 1:
            raise TypeError("x must be one of the generators of the parent.")
        i = nonzero_positions[0]
        if L[i] != 1:
            raise TypeError("x must be one of the generators of the parent.")
        _min = []
        for v in self.__repn.keys():
            _min.append(v[i])
        return min(_min)

    def total_degree(PolyDict self):
        return max([-1] + map(sum,self.__repn.keys()))

    def monomial_coefficient(PolyDict self, mon):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.monomial_coefficient(PolyDict({(2,1):1}).dict())
            4
        """
        K = mon.keys()[0]
        if K not in self.__repn:
            return 0
        return self.__repn[K]

    def polynomial_coefficient(PolyDict self, degrees):
        """
        Return a polydict that defines the coefficient in the current
        polynomial viewed as a tower of polynomial extensions.

        INPUT:

        - ``degrees`` -- a list of degree restrictions; list elements are None
          if the variable in that position should be unrestricted

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.polynomial_coefficient([2,None])
            PolyDict with representation {(0, 3): 2, (0, 1): 4}
            sage: f = PolyDict({(0,3):2, (0,2):3, (2,1):4})
            sage: f.polynomial_coefficient([0,None])
            PolyDict with representation {(0, 3): 2, (0, 2): 3}
        """
        nz = []
        cdef int i
        for i from 0<=i<len(degrees):
            if degrees[i] is not None:
                nz.append(i)
        ans = {}
        for S in self.__repn.keys():
            exactly_divides = True
            for j in nz:
                if S[j] != degrees[j]:
                    exactly_divides = False
                    break
            if exactly_divides:
                t = list(S)
                for m in nz:
                    t[m] = 0
                ans[ETuple(t)] = self.__repn[S]
        return PolyDict(ans, force_etuples=False)

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

    def homogenize(PolyDict self, var):
        R = self.__repn
        H = {}
        deg = self.degree()
        for e, val in R.iteritems():
            i = deg - sum(e)
            f = list(e)
            f[var] += i
            H[ETuple(f)] = val
        return PolyDict(H, zero=self.__zero, force_etuples=False)

    def latex(PolyDict self, vars, atomic_exponents=True, atomic_coefficients=True, cmpfn = None):
        r"""
        Return a nice polynomial latex representation of this PolyDict, where
        the vars are substituted in.

        INPUT:

        - ``vars`` -- list
        - ``atomic_exponents`` -- bool (default: True)
        - ``atomic_coefficients`` -- bool (default: True)

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.latex(['a','WW'])
            '2 a^{2} WW^{3} + 4 a^{2} WW + 3 a WW^{2}'

        When ``atomic_exponents`` is False, the exponents are surrounded in
        parenthesis, since ^ has such high precedence::

            # I've removed fractional exponent support in ETuple when moving to a sparse C integer array
            #sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            #sage: f.latex(['a','b','c'], atomic_exponents=False)
            #'4 a^{2}bc + 3 ab^{2}c + 2 a^{2/3}b^{3}c^{5}'

        TESTS:

        We check that the issue on Trac 9478 is resolved::

            sage: R2.<a> = QQ[]
            sage: R3.<xi, x> = R2[]
            sage: print latex(xi*x)
            \xi x
        """
        n = len(vars)
        poly = ""
        E = self.__repn.keys()
        if cmpfn:
            E.sort(cmp = cmpfn, reverse=True)
        else:
            E.sort(reverse=True)
        try:
            pos_one = self.__zero.parent()(1)
            neg_one = -pos_one
        except AttributeError:
            # probably self.__zero is not a ring element
            pos_one = 1
            neg_one = -1
        for e in E:
            c = self.__repn[e]
            if not c == self.__zero:
                sign_switch = False
                # First determine the multinomial:
                multi = " ".join([vars[j] +
                                  ("^{%s}" % e[j] if e[j] != 1 else "")
                                 for j in e.nonzero_positions(sort=True)])
                # Next determine coefficient of multinomial
                if len(multi) == 0:
                    multi = latex(c)
                elif c == neg_one:
                    # handle -1 specially because it's a pain
                    if len(poly) > 0:
                        sign_switch = True
                    else:
                        multi = "-%s"%(multi)
                elif c != pos_one:
                    if not atomic_coefficients:
                        c = latex(c)
                        if c.find("+") != -1 or c.find("-") != -1 or c.find(" ") != -1:
                            c = "(%s)"%c
                    multi = "%s %s"%(c,multi)

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

        - ``vars`` -- list
        - ``atomic_exponents`` -- bool (default: True)
        - ``atomic_coefficients`` -- bool (default: True)

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.poly_repr(['a','WW'])
            '2*a^2*WW^3 + 4*a^2*WW + 3*a*WW^2'

        When atomic_exponents is False, the exponents are surrounded
        in parenthesis, since ^ has such high precedence. ::

            # I've removed fractional exponent support in ETuple when moving to a sparse C integer array
            #sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            #sage: f.poly_repr(['a','b','c'], atomic_exponents=False)
            #'4*a^(2)*b*c + 3*a*b^(2)*c + 2*a^(2/3)*b^(3)*c^(5)'

        We check to make sure that when we are in characteristic two, we
        don't put negative signs on the generators. ::

            sage: Integers(2)['x,y'].gens()
            (x, y)

        We make sure that intervals are correctly represented. ::

            sage: f = PolyDict({(2,3):RIF(1/2,3/2), (1,2):RIF(-1,1)})
            sage: f.poly_repr(['x','y'])
            '1.?*x^2*y^3 + 0.?*x*y^2'
        """
        n = len(vars)
        poly = ""
        E = self.__repn.keys()
        if cmpfn:
            E.sort(cmp = cmpfn, reverse=True)
        else:
            E.sort(reverse=True)
        try:
            pos_one = self.__zero.parent()(1)
            neg_one = -pos_one
        except AttributeError:
            # probably self.__zero is not a ring element
            pos_one = 1
            neg_one = -1

        is_characteristic_2 = bool(pos_one == neg_one)

        for e in E:
            c = self.__repn[e]
            if not c == self.__zero:
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
                elif c == neg_one and not is_characteristic_2:
                    # handle -1 specially because it's a pain
                    if len(poly) > 0:
                        sign_switch = True
                    else:
                        multi = "-%s"%(multi)
                elif not c == pos_one:
                    if not atomic_coefficients:
                        c = str(c)
                        if c.find("+") != -1 or c.find("-") != -1 or c.find(" ") != -1:
                            c = "(%s)"%c
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
        We add two polynomials in 2 variables::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(1,5):-3, (2,3):-2, (1,1):3})
            sage: f+g
            PolyDict with representation {(1, 2): 3, (2, 1): 4, (1, 1): 3, (1, 5): -3}

        Next we add two polynomials with fractional exponents in 3 variables::

            # I've removed fractional exponent support in ETuple when moving to a sparse C integer array
            #sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            #sage: g = PolyDict({(2/3,3,5):3}, force_int_exponents=False)
            #sage: f+g
            #PolyDict with representation {(1, 2, 1): 3, (2/3, 3, 5): 5, (2, 1, 1): 4}
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
        We multiply two polynomials in 2 variables::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(1,5):-3, (2,3):-2, (1,1):3})
            sage: f*g
            PolyDict with representation {(2, 3): 9, (3, 2): 12, (2, 7): -9, (3, 5): -6, (4, 6): -4, (3, 4): 6, (3, 6): -12, (4, 4): -8, (3, 8): -6}

        Next we multiply two polynomials with fractional exponents in 3 variables::

            # I've removed fractional exponent support in ETuple when moving to a sparse C integer array
            #sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            #sage: g = PolyDict({(2/3,3,5):3}, force_int_exponents=False)
            #sage: f*g
            #PolyDict with representation {(8/3, 4, 6): 12, (5/3, 5, 6): 9, (4/3, 6, 10): 6}

        Finally we print the result in a nice format. ::

            # I've removed fractional exponent support in ETuple when moving to a sparse C integer array
            #sage: (f*g).poly_repr(['a','b','c'], atomic_exponents = False)
            #'12*a^(8/3)*b^(4)*c^(6) + 9*a^(5/3)*b^(5)*c^(6) + 6*a^(4/3)*b^(6)*c^(10)'
        """
        cdef PyObject *cc
        newpoly = {}
        if len(self.__repn) == 0:   # product is zero anyways
            return self
        for e0, c0 in self.__repn.iteritems():
            for e1, c1 in right.__repn.iteritems():
                e = (<ETuple>e0).eadd(<ETuple>e1)
                c = c0*c1
                cc = PyDict_GetItem(newpoly,e)
                if cc == <PyObject*>0:
                    PyDict_SetItem(newpoly,e,c)
                else:
                    PyDict_SetItem(newpoly,e,<object>cc+c)
        F = PolyDict(newpoly, self.__zero, force_int_exponents=False, remove_zero=True, force_etuples=False)
        return F

    def scalar_rmult(PolyDict self, s):
        """
        Right Scalar Multiplication

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: x,y=FreeMonoid(2,'x,y').gens()  # a strange object to live in a polydict, but non-commutative!
            sage: f = PolyDict({(2,3):x})
            sage: f.scalar_rmult(y)
            PolyDict with representation {(2, 3): x*y}
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.scalar_rmult(-2)
            PolyDict with representation {(1, 2): -6, (2, 3): -4, (2, 1): -8}
            sage: f.scalar_rmult(RIF(-1,1))
            PolyDict with representation {(1, 2): 0.?e1, (2, 3): 0.?e1, (2, 1): 0.?e1}
        """
        v = {}
        # if s is 0, then all the products will be zero
        if not s == self.__zero:
            for e, c in self.__repn.iteritems():
                v[e] = c*s
        return PolyDict(v, self.__zero, force_int_exponents=False, force_etuples=False)

    def scalar_lmult(PolyDict self, s):
        """
        Left Scalar Multiplication

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: x,y=FreeMonoid(2,'x,y').gens()  # a strange object to live in a polydict, but non-commutative!
            sage: f = PolyDict({(2,3):x})
            sage: f.scalar_lmult(y)
            PolyDict with representation {(2, 3): y*x}
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.scalar_lmult(-2)
            PolyDict with representation {(1, 2): -6, (2, 3): -4, (2, 1): -8}
            sage: f.scalar_lmult(RIF(-1,1))
            PolyDict with representation {(1, 2): 0.?e1, (2, 3): 0.?e1, (2, 1): 0.?e1}
        """
        v = {}
        # if s is 0, then all the products will be zero
        if not s == self.__zero:
            for e, c in self.__repn.iteritems():
                v[e] = s*c
        return PolyDict(v, self.__zero, force_int_exponents=False, force_etuples=False)

    def __sub__(PolyDict self,PolyDict  other):
        """
        Subtract two PolyDict's.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(2,3):2, (1,1):-10})
            sage: f - g
            PolyDict with representation {(1, 2): 3, (2, 1): 4, (1, 1): 10}
            sage: g - f
            PolyDict with representation {(1, 2): -3, (1, 1): -10, (2, 1): -4}
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

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f**2
            PolyDict with representation {(4, 2): 16, (3, 3): 24, (3, 5): 12, (4, 6): 4, (2, 4): 9, (4, 4): 16}
            sage: f**0
            PolyDict with representation {(0, 0): 1}
            sage: (f-f)**0
            PolyDict with representation {0: 1}
        """
        return generic_power(self, n, self.__one())

    def lcmt(PolyDict self, greater_etuple):
        """
        Provides functionality of lc, lm, and lt by calling the tuple
        compare function on the provided term order T.

        INPUT:

        - ``greater_etuple`` -- a term order
        """
        try:
            return ETuple(reduce(greater_etuple,self.__repn.keys()))
        except KeyError:
            raise ArithmeticError("%s not supported",greater_etuple)

    def __reduce__(PolyDict self):
        """

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: loads(dumps(f)) == f
            True
        """
        return make_PolyDict,(self.__repn,)

    def min_exp(self):
        """
        Returns an ETuple containing the minimum exponents appearing.  If
        there are no terms at all in the PolyDict, it returns None.

        The nvars parameter is necessary because a PolyDict doesn't know it
        from the data it has (and an empty PolyDict offers no clues).

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.min_exp()
            (1, 1)
            sage: PolyDict({}).min_exp() # returns None
        """
        cdef ETuple r
        ETuples = self.__repn.keys()
        if len(ETuples)>0:
            r = <ETuple>ETuples[0]
            for e in ETuples:
                r = r.emin(e)
            return r
        else:
            return None

    def max_exp(self):
        """
        Returns an ETuple containing the maximum exponents appearing.  If
        there are no terms at all in the PolyDict, it returns None.

        The nvars parameter is necessary because a PolyDict doesn't know it
        from the data it has (and an empty PolyDict offers no clues).

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.max_exp()
            (2, 3)
            sage: PolyDict({}).max_exp() # returns None
        """
        cdef ETuple r
        ETuples = self.__repn.keys()
        if len(ETuples)>0:
            r = <ETuple>ETuples[0]
            for e in ETuples:
                r = r.emax(e)
            return r
        else:
            return None

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

cdef inline bint dual_etuple_iter(ETuple self, ETuple other, size_t *ind1, size_t *ind2, size_t *index, int *exp1, int *exp2):
    """
    This function is a crucial helper function for a number of methods of
    the ETuple class.

    This is a rather fragile function.  Perhaps some Cython guru could make
    it appear a little less stilted -- a principal difficulty is passing
    C types by reference.  In any case, the complicated features of looping
    through two ETuple _data members is all packaged up right here and
    shouldn't be allowed to spread.
    """
    if ind1[0] >= self._nonzero and ind2[0] >= other._nonzero:
        return 0
    if ind1[0] < self._nonzero and ind2[0] < other._nonzero:
        if self._data[2*ind1[0]] == other._data[2*ind2[0]]:
            exp1[0] = self._data[2*ind1[0]+1]
            exp2[0] = other._data[2*ind2[0]+1]
            index[0] = self._data[2*ind1[0]]
            ind1[0] += 1
            ind2[0] += 1
        elif self._data[2*ind1[0]] > other._data[2*ind2[0]]:
            exp1[0] = 0
            exp2[0] = other._data[2*ind2[0]+1]
            index[0] = other._data[2*ind2[0]]
            ind2[0] += 1
        else:
            exp1[0] = self._data[2*ind1[0]+1]
            exp2[0] = 0
            index[0] = self._data[2*ind1[0]]
            ind1[0] += 1
    else:
        if ind2[0] >= other._nonzero:
            exp1[0] = self._data[2*ind1[0]+1]
            exp2[0] = 0
            index[0] = self._data[2*ind1[0]]
            ind1[0] += 1
        elif ind1[0] >= self._nonzero:
            exp1[0] = 0
            exp2[0] = other._data[2*ind2[0]+1]
            index[0] = other._data[2*ind2[0]]
            ind2[0] += 1
    return 1

cdef class ETuple:
    """
    Representation of the exponents of a polydict monomial. If
    (0,0,3,0,5) is the exponent tuple of x_2^3*x_4^5 then this class
    only stores {2:3,4:5} instead of the full tuple. This sparse
    information may be obtained by provided methods.

    The index/value data is all stored in the _data C int array member
    variable.  For the example above, the C array would contain
    2,3,4,5.  The indices are interlaced with the values.

    This data structure is very nice to work with for some functions
    implemented in this class, but tricky for others.  One reason that
    I really like the format is that it requires a single memory
    allocation for all of the values.  A hash table would require more
    allocations and presumably be slower.  I didn't benchmark this
    question (although, there is no question that this is much faster
    than the prior use of python dicts).
    """
    cdef ETuple _new(ETuple self):
        """
        Quickly creates a new initialized ETuple with the
        same length as self.
        """
        cdef ETuple x
        x = <ETuple>PY_NEW_SAME_TYPE(self)
        x._length = self._length
        return x

    def __init__(ETuple self, data=None, length=None):
        """
        - ``ETuple()`` -> an empty ETuple
        - ``ETuple(sequence)`` -> ETuple initialized from sequence's items

        If the argument is an ETuple, the return value is the same object.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1,1,0])
            (1, 1, 0)
            sage: ETuple({int(1):int(2)},int(3))
            (0, 2, 0)
        """
        if data is None:
            return
        cdef size_t ind
        cdef int v
        if PY_TYPE_CHECK(data,ETuple):
            self._length = (<ETuple>data)._length
            self._nonzero = (<ETuple>data)._nonzero
            self._data = <int*>sage_malloc(sizeof(int)*self._nonzero*2)
            memcpy(self._data,(<ETuple>data)._data,sizeof(int)*self._nonzero*2)
        elif PY_TYPE_CHECK(data,dict) and PY_TYPE_CHECK(length,int):
            self._length = length
            self._nonzero = len(data)
            self._data = <int*>sage_malloc(sizeof(int)*self._nonzero*2)
            nz_elts = sorted(data.iteritems())
            ind = 0
            for index,exp in nz_elts:
                self._data[2*ind] = index
                self._data[2*ind+1] = exp
                ind += 1
        elif PY_TYPE_CHECK(data,list) or PY_TYPE_CHECK(data,tuple):
            self._length = len(data)
            self._nonzero = 0
            for v in data:
                if v != 0:
                    self._nonzero += 1
            ind = 0
            self._data = <int*>sage_malloc(sizeof(int)*self._nonzero*2)
            for i from 0 <= i < self._length:
                v = data[i]
                if v != 0:
                    self._data[ind] = i
                    self._data[ind+1] = v
                    ind += 2
        else:
            raise TypeError

    def __cinit__(ETuple self):
        self._data = <int*>0

    def __dealloc__(self):
        if self._data != <int*>0:
            sage_free(self._data)

    # methods to simulate tuple

    def __add__(ETuple self,ETuple other):
        """
        x.__add__(n) <==> x+n

        concatenates two ETuples

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1,1,0]) + ETuple({int(1):int(2)},int(3))
            (1, 1, 0, 0, 2, 0)
        """
        cdef size_t index = 0
        cdef ETuple result = <ETuple>PY_NEW(ETuple)
        result._length = self._length+other._length
        result._nonzero = self._nonzero+other._nonzero
        result._data = <int*>sage_malloc(sizeof(int)*result._nonzero*2)
        for index from 0 <= index < self._nonzero:
            result._data[2*index] = self._data[2*index]
            result._data[2*index+1] = self._data[2*index+1]
        for index from 0 <= index < other._nonzero:
            result._data[2*(index+self._nonzero)] = other._data[2*index]+self._length # offset the second tuple (append to end!)
            result._data[2*(index+self._nonzero)+1] = other._data[2*index+1]
        return result

    def __mul__(ETuple self,factor):
        """
        x.__mul__(n) <==> x*n

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1,2,3])*2
            (1, 2, 3, 1, 2, 3)
        """
        cdef int _factor = factor
        cdef ETuple result = <ETuple>PY_NEW(ETuple)
        if factor <= 0:
            result._length = 0
            result._nonzero = 0
            return result
        cdef size_t index
        cdef size_t f
        result._length = self._length*factor
        result._nonzero = self._nonzero*factor
        result._data = <int*>sage_malloc(sizeof(int)*result._nonzero*2)
        for index from 0 <= index < self._nonzero:
            for f from 0 <= f < factor:
                result._data[2*(f*self._nonzero+index)] = self._data[2*index]+f*self._length
                result._data[2*(f*self._nonzero+index)+1] = self._data[2*index+1]
        return result

    def __getitem__(ETuple self, i):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: m = ETuple([1,2,0,3])
            sage: m[2]
            0
            sage: m[1]
            2
            sage: e = ETuple([1,2,3])
            sage: e[1:]
            (2, 3)
            sage: e[:1]
            (1,)
        """
        cdef size_t ind
        if isinstance(i, slice):
            start, stop = i.start, i.stop
            if start is None:
                start = 0
            elif start < 0:
                start = start % self._length
            elif start > self._length:
                start = self._length

            if stop > self._length or stop is None:
                stop = self._length
            elif stop < 0:
                stop = stop % self._length

            # this is not particularly fast, but I doubt many people care
            # if you do, feel free to tweak!
            d = [self[ind] for ind from start <= ind < stop]
            return ETuple(d)
        else:
            ind = 0
            for ind from 0 <= ind < self._nonzero:
                if self._data[2*ind] == i:
                    return self._data[2*ind+1]
                elif self._data[2*ind] > i:
                    # the indices are sorted in _data, we are beyond, so quit
                    return 0
            return 0

    def __hash__(ETuple self):
        """
        x.__hash__() <==> hash(x)
        """
        cdef int i
        cdef int result = 0
        for i from 0 <= i < self._nonzero:
            result += (1000003 * result) ^ self._data[2*i]
            result += (1000003 * result) ^ self._data[2*i+1]
        result = (1000003 * result) ^ self._length
        return result

    def __len__(ETuple self):
        """
        x.__len__() <==> len(x)

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e=ETuple([1,0,2,0,3])
            sage: len(e)
            5
        """
        return self._length

    def __contains__(ETuple self, elem):
        """
        x.__contains__(n) <==> n in x


        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple({int(1):int(2)},int(3)); e
            (0, 2, 0)
            sage: 1 in e
            False
            sage: 2 in e
            True
        """
        if elem==0:
            return self._length > self._nonzero

        cdef size_t ind = 0
        for ind from 0 <= ind < self._nonzero:
            if elem == self._data[2*ind+1]:
                return True
        return False

    def __richcmp__(ETuple self, ETuple other, op):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: ETuple([1,1,0])<ETuple([1,1,0])
            False
            sage: ETuple([1,1,0])<ETuple([1,0,0])
            False
            sage: ETuple([1,1,0])<ETuple([1,2,0])
            True
            sage: ETuple([1,1,0])<ETuple([1,-1,0])
            False
            sage: ETuple([0,-2,0])<ETuple([1,-1,0])
            True
            sage: ETuple([1,1,0])>ETuple([1,1,0])
            False
            sage: ETuple([1,1,0])>ETuple([1,0,0])
            True
            sage: ETuple([1,1,0])>ETuple([1,2,0])
            False
            sage: ETuple([1,1,0])>ETuple([1,-1,0])
            True
            sage: ETuple([0,-2,0])>ETuple([1,-1,0])
            False
        """
        cdef size_t ind = 0
        if op == 2: #==
            if self._nonzero != other._nonzero:
                return False
            for ind from 0 <= ind < self._nonzero:
                if self._data[2*ind] != other._data[2*ind]:
                    return False
                if self._data[2*ind+1] != other._data[2*ind+1]:
                    return False
            return self._length == other._length

        if op == 0:  #<
            while ind < self._nonzero and ind < other._nonzero:
                if self._data[2*ind] < other._data[2*ind]:
                    return self._data[2*ind+1] < 0
                if self._data[2*ind] > other._data[2*ind]:
                    return other._data[2*ind+1] > 0
                if self._data[2*ind] == other._data[2*ind] and self._data[2*ind+1] != other._data[2*ind+1]:
                    return self._data[2*ind+1] < other._data[2*ind+1]
                ind += 1
            if ind < self._nonzero and ind == other._nonzero:
                return self._data[2*ind+1] < 0
            if ind < other._nonzero and ind == self._nonzero:
                return other._data[2*ind+1] > 0
            return self._length < other._length

        if op == 4: #>
            while ind < self._nonzero and ind < other._nonzero:
                if self._data[2*ind] < other._data[2*ind]:
                    return self._data[2*ind+1] > 0
                if self._data[2*ind] > other._data[2*ind]:
                    return other._data[2*ind+1] < 0
                if self._data[2*ind] == other._data[2*ind] and self._data[2*ind+1] != other._data[2*ind+1]:
                    return self._data[2*ind+1] > other._data[2*ind+1]
                ind += 1
            if ind < self._nonzero and ind == other._nonzero:
                return self._data[2*ind+1] > 0
            if ind < other._nonzero and ind == self._nonzero:
                return other._data[2*ind+1] < 0
            return self._length < other._length

        # the rest of these are not particularly fast

        if op == 1: #<=
            return tuple(self) <= tuple(other)

        if op == 3: #!=
            return tuple(self) != tuple(other)

        if op == 5: #>=
            return tuple(self) >= tuple(other)

    def __iter__(ETuple self):
        """
        x.__iter__() <==> iter(x)
        """
        cdef size_t ind
        # this is not particularly fast, but I doubt many people care
        # if you do, feel free to tweak!
        d = dict([(self._data[2*ind],self._data[2*ind+1]) for ind from 0<=ind<self._nonzero])
        return ETupleIter(d,self._length)

    def __str__(ETuple self):
        return self.__repr__()

    def __repr__(ETuple self):
        res = [0,]*self._length
        cdef size_t ind = 0
        for ind from 0 <= ind < self._nonzero:
            res[self._data[2*ind]] = self._data[2*ind+1]
        return str(tuple(res))

    def __reduce__(ETuple self):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,1,0])
            sage: bool(e == loads(dumps(e)))
            True
        """
        cdef size_t ind
        # this is not particularly fast, but I doubt many people care
        # if you do, feel free to tweak!
        d = dict([(self._data[2*ind],self._data[2*ind+1]) for ind from 0<=ind<self._nonzero])
        return make_ETuple,(d,int(self._length))

    # additional methods

    cpdef ETuple eadd(ETuple self,ETuple other):
        """
        Vector addition of self with other.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: f = ETuple([0,1,1])
            sage: e.eadd(f)
            (1, 1, 3)

        Verify that trac 6428 has been addressed::

            sage: R.<y,z> = Frac(QQ['x'])[]
            sage: type(y)
            <class 'sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict'>
            sage: y^(2^32)
            Traceback (most recent call last):
            ...
            OverflowError: Exponent overflow (2147483648).
        """
        if self._length!=other._length:
            raise ArithmeticError

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef int s  # sum
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sage_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self,other,&ind1,&ind2,&index,&exp1,&exp2):
            s = exp1 + exp2
            # Check for overflow and underflow
            if (exp2 > 0 and s < exp1) or (exp2 < 0 and s > exp1):
                raise OverflowError("Exponent overflow (%s)."%(int(exp1)+int(exp2)))
            if s != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = s
                result._nonzero += 1
        return result

    cpdef ETuple eadd_p(ETuple self, int other, int pos):
        """
        Adds other to self at position pos.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: e.eadd_p(5, 1)
            (1, 5, 2)
            sage: e = ETuple([0]*7)
            sage: e.eadd_p(5,4)
            (0, 0, 0, 0, 5, 0, 0)

            sage: ETuple([0,1]).eadd_p(1, 0) == ETuple([1,1])
            True

        """
        cdef size_t index = 0
        cdef size_t rindex = 0
        cdef int new_value
        cdef int need_to_add = 1
        if pos < 0 or pos >= self._length:
            raise ValueError("pos must be between 0 and %s"%self._length)

        cdef size_t alloc_len = self._nonzero + 1

        cdef ETuple result = <ETuple>self._new()
        result._nonzero = self._nonzero
        result._data = <int*>sage_malloc(sizeof(int)*alloc_len*2)

        for index from 0 <= index < self._nonzero:
            if self._data[2*index] == pos:
                new_value = self._data[2*index+1] + other
                if new_value != 0:
                    result._data[2*rindex] = pos
                    result._data[2*rindex+1] = new_value
                else:
                    result._nonzero -= 1
                    rindex -= 1
                need_to_add = 0
            else:
                result._data[2*rindex] = self._data[2*index]
                result._data[2*rindex+1] = self._data[2*index+1]

            rindex += 1

        rindex = 0
        if need_to_add:
            for index from 0 <= index < self._nonzero:
                if self._data[2*index] > pos:
                    result._data[2*rindex] = pos
                    result._data[2*rindex+1] = other
                    rindex += 1
                    result._nonzero += 1
                result._data[2*rindex] = self._data[2*index]
                result._data[2*rindex+1] = self._data[2*index+1]
                rindex += 1

            if rindex == index and other:
                result._data[2*rindex] = pos
                result._data[2*rindex+1] = other
                result._nonzero += 1

        return result

    cpdef ETuple esub(ETuple self,ETuple other):
        """
        Vector subtraction of self with other.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: f = ETuple([0,1,1])
            sage: e.esub(f)
            (1, -1, 1)
        """
        if self._length!=other._length:
            raise ArithmeticError

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef int d  # difference
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sage_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self,other,&ind1,&ind2,&index,&exp1,&exp2):
            # Check for overflow and underflow
            d = exp1 - exp2
            if (exp2 > 0 and d > exp1) or (exp2 < 0 and d < exp1):
                raise OverflowError("Exponent overflow (%s)."%(int(exp1)-int(exp2)))
            if d != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = d
                result._nonzero += 1
        return result

    cpdef ETuple emul(ETuple self,int factor):
        """
        Scalar Vector multiplication of self.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: e.emul(2)
            (2, 0, 4)
        """
        cdef size_t ind
        cdef ETuple result = <ETuple>self._new()
        if factor == 0:
            result._nonzero = 0  # all zero, no non-zero entries!
            result._data = <int*>sage_malloc(sizeof(int)*result._nonzero*2)
        else:
            result._nonzero = self._nonzero
            result._data = <int*>sage_malloc(sizeof(int)*result._nonzero*2)
            for ind from 0 <= ind < self._nonzero:
                result._data[2*ind] = self._data[2*ind]
                result._data[2*ind+1] = self._data[2*ind+1]*factor
        return result

    cpdef ETuple emax(ETuple self,ETuple other):
        """
        Vector of maximum of components of self and other.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: f = ETuple([0,1,1])
            sage: e.emax(f)
            (1, 1, 2)
            sage: e=ETuple((1,2,3,4))
            sage: f=ETuple((4,0,2,1))
            sage: f.emax(e)
            (4, 2, 3, 4)
            sage: e=ETuple((1,-2,-2,4))
            sage: f=ETuple((4,0,0,0))
            sage: f.emax(e)
            (4, 0, 0, 4)
            sage: f.emax(e).nonzero_positions()
            [0, 3]
        """
        if self._length!=other._length:
            raise ArithmeticError

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sage_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self,other,&ind1,&ind2,&index,&exp1,&exp2):
            if exp1 >= exp2 and exp1 != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = exp1
                result._nonzero += 1
            elif exp2 >= exp1 and exp2 != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = exp2
                result._nonzero += 1
        return result

    cpdef ETuple emin(ETuple self,ETuple other):
        """
        Vector of minimum of components of self and other.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: f = ETuple([0,1,1])
            sage: e.emin(f)
            (0, 0, 1)
            sage: e = ETuple([1,0,-1])
            sage: f = ETuple([0,-2,1])
            sage: e.emin(f)
            (0, -2, -1)
        """
        if self._length!=other._length:
            raise ArithmeticError

        cdef size_t ind1 = 0
        cdef size_t ind2 = 0
        cdef size_t index
        cdef int exp1
        cdef int exp2
        cdef size_t alloc_len = self._nonzero + other._nonzero  # we simply guesstimate the length -- there might be double the correct amount allocated -- who cares?
        if alloc_len > self._length:
            alloc_len = self._length
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = 0  # we don't know the correct length quite yet
        result._data = <int*>sage_malloc(sizeof(int)*alloc_len*2)
        while dual_etuple_iter(self,other,&ind1,&ind2,&index,&exp1,&exp2):
            if exp1 <= exp2 and exp1 != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = exp1
                result._nonzero += 1
            elif exp2 <= exp1 and exp2 != 0:
                result._data[2*result._nonzero] = index
                result._data[2*result._nonzero+1] = exp2
                result._nonzero += 1
        return result

    def nonzero_positions(ETuple self,sort=False):
        """
        Returns the positions of non-zero exponents in the tuple.

        INPUT:

        - ``sort`` -- if True a sorted list is returned. If False an unsorted
          list is returned. (default: False)

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: e.nonzero_positions()
            [0, 2]
        """
        cdef size_t ind
        return [self._data[2*ind] for ind from 0 <= ind < self._nonzero]

    def common_nonzero_positions(ETuple self, ETuple other,sort=False):
        """
        Returns an optionally sorted list of non zero positions either
        in self or other, i.e. the only positions that need to be
        considered for any vector operation.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2])
            sage: f = ETuple([0,0,1])
            sage: e.common_nonzero_positions(f)
            set([0, 2])
            sage: e.common_nonzero_positions(f,sort=True)
            [0, 2]
        """
        # TODO:  we should probably make a fast version of this!
        res = set(self.nonzero_positions()).union(other.nonzero_positions())
        if sort:
            return sorted(res)
        else:
            return res

    def nonzero_values(ETuple self, sort=True):
        """
        Returns the non-zero values of the tuple.

        INPUT:

        - ``sort`` -- if True the values are sorted by their indices. Otherwise
          the values are returned unsorted. (default: True)

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([2,0,1])
            sage: e.nonzero_values()
            [2, 1]
            sage: f = ETuple([0,-1,1])
            sage: f.nonzero_values(sort=True)
            [-1, 1]
        """
        cdef size_t ind
        return [self._data[2*ind+1] for ind from 0 <= ind < self._nonzero]

    def reversed(ETuple self):
        """
        Returns the reversed ETuple of self.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,2,3])
            sage: e.reversed()
            (3, 2, 1)
        """
        cdef size_t ind
        cdef ETuple result = <ETuple>self._new()
        result._nonzero = self._nonzero
        result._data = <int*>sage_malloc(sizeof(int)*result._nonzero*2)
        for ind from 0 <= ind < self._nonzero:
            result._data[2*(result._nonzero-ind-1)] = self._length-self._data[2*ind]-1
            result._data[2*(result._nonzero-ind-1)+1] = self._data[2*ind+1]
        return result

    def sparse_iter(ETuple self):
        """
        Iterator over the elements of self where the elements are returned as
        ``(i,e)`` where ``i`` is the position of ``e`` in the tuple.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([1,0,2,0,3])
            sage: list(e.sparse_iter())
            [(0, 1), (2, 2), (4, 3)]
        """
        cdef size_t ind
        d = dict([(self._data[2*ind],self._data[2*ind+1]) for ind from 0<=ind<self._nonzero])
        return d.iteritems()

    def combine_to_positives(ETuple self, ETuple other):
        """
        Given a pair of ETuples (self, other), returns a triple of ETuples (a, b, c) so that self = a + b,other = a + c and b and c have all positive entries.

        EXAMPLES::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: e = ETuple([-2,1,-5, 3, 1,0])
            sage: f = ETuple([1,-3,-3,4,0,2])
            sage: e.combine_to_positives(f)
            ((-2, -3, -5, 3, 0, 0), (0, 4, 0, 0, 1, 0), (3, 0, 2, 1, 0, 2))
        """
        m = self.emin(other)
        return m, self.esub(m), other.esub(m)

def make_PolyDict(data):
    return PolyDict(data,remove_zero=False, force_int_exponents=False, force_etuples=False)

def make_ETuple(data,length):
    return ETuple(data,length)
