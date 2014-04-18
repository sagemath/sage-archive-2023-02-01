"""
Free algebra elements

AUTHORS:

- David Kohel (2005-09)

TESTS::

    sage: R.<x,y> = FreeAlgebra(QQ,2)
    sage: x == loads(dumps(x))
    True
    sage: x*y
    x*y
    sage: (x*y)^0
    1
    sage: (x*y)^3
    x*y*x*y*x*y
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.misc import repr_lincomb
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.algebras.algebra_element import AlgebraElement

class FreeAlgebraElement(AlgebraElement):
    """
    A free algebra element.
    """
    def __init__(self, A, x):
        """
        Create the element ``x`` of the FreeAlgebra ``A``.

        TESTS::

            sage: F.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: elt = x^3 * y - z^2*x
            sage: TestSuite(elt).run()
        """
        if isinstance(x, FreeAlgebraElement):
            x = x.__monomial_coefficients
        AlgebraElement.__init__(self, A)
        R = A.base_ring()
        if isinstance(x, AlgebraElement): #and x.parent() == A.base_ring():
            self.__monomial_coefficients = { A.monoid()(1):R(x) }
        elif isinstance(x, FreeMonoidElement):
            self.__monomial_coefficients = { x:R(1) }
        elif True:
            self.__monomial_coefficients = dict([ (A.monoid()(e1),R(e2)) for e1,e2 in x.items()])
        else:
            raise TypeError("Argument x (= %s) is of the wrong type."%x)

    def __iter__(self):
        """
        Returns an iterator which yields tuples of coefficient and monomial.

        EXAMPLES::

            sage: a = FreeAlgebra(QQ, 5, 'a').gens()
            sage: list(3*a[0]*a[1]*a[4]**3*a[0]+1)
            [(1, 1), (3, a0*a1*a4^3*a0)]
        """
        for (key,val) in self.__monomial_coefficients.iteritems():
            yield (val,key)

    def _repr_(self):
        """
        Return string representation of self.

        EXAMPLES::

            sage: A.<x,y,z>=FreeAlgebra(ZZ,3)
            sage: repr(-x+3*y*z)    # indirect doctest
            '-x + 3*y*z'

        Trac ticket #11068 enables the use of local variable names::

            sage: from sage.structure.parent_gens import localvars
            sage: with localvars(A, ['a','b','c']):
            ...    print -x+3*y*z
            ...
            -a + 3*b*c

        """
        v = sorted(self.__monomial_coefficients.items())
        P = self.parent()
        M = P.monoid()
        from sage.structure.parent_gens import localvars
        with localvars(M, P.variable_names(), normalize=False):
            x = repr_lincomb(v, strip_one=True)
        return x

    def _latex_(self):
        r"""
        Return latex representation of self.

        EXAMPLES::

            sage: A.<x,y,z>=FreeAlgebra(ZZ,3)
            sage: latex(-x+3*y^20*z)   # indirect doctest
            -x + 3y^{20}z
            sage: alpha,beta,gamma=FreeAlgebra(ZZ,3,'alpha,beta,gamma').gens()
            sage: latex(alpha-beta)
            \alpha - \beta
        """
        v = sorted(self.__monomial_coefficients.items())
        return repr_lincomb(v, strip_one=True, is_latex=True)

    def __call__(self, *x, **kwds):
        """
        EXAMPLES::

            sage: A.<x,y,z>=FreeAlgebra(ZZ,3)
            sage: (x+3*y).subs(x=1,y=2,z=14)
            7
            sage: (2*x+y).subs({x:1,y:z})
            2 + z
            sage: f=x+3*y+z
            sage: f(1,2,1/2)
            15/2
            sage: f(1,2)
            Traceback (most recent call last):
            ...
            ValueError: must specify as many values as generators in parent

        AUTHORS:

        - Joel B. Mohler (2007-10-27)
        """
        if len(kwds)>0 and len(x)>0:
            raise ValueError("must not specify both a keyword and positional argument")

        if len(kwds)>0:
            p = self.parent()
            def extract_from(kwds,g):
                for x in g:
                    try:
                        return kwds[x]
                    except KeyError:
                        pass
                return None

            x = [extract_from(kwds,(p.gen(i),p.variable_name(i))) for i in range(p.ngens())]
        elif isinstance(x[0],tuple):
            x = x[0]

        if len(x) != self.parent().ngens():
            raise ValueError("must specify as many values as generators in parent")

        # I don't start with 0, because I don't want to preclude evaluation with
        #arbitrary objects (e.g. matrices) because of funny coercion.
        result = None
        for m, c in self.__monomial_coefficients.iteritems():
            if result is None:
                result = c*m(x)
            else:
                result += c*m(x)

        if result is None:
            return self.parent()(0)
        return result

    def __cmp__(left, right):
        """
        Compare two free algebra elements with the same parents.

        The ordering is the one on the underlying sorted list of
        (monomial,coefficients) pairs.

        EXAMPLES::

            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: x < y
            True
            sage: x * y < y * x
            True
            sage: y * x < x * y
            False
        """
        v = sorted(left.__monomial_coefficients.items())
        w = sorted(right.__monomial_coefficients.items())
        return cmp(v, w)

    def _add_(self, y):
        """
        Return sum of self and y (another free algebra element with the
        same parents)

        EXAMPLES::

            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: x + y    # indirect doctest
            x + y
        """
        A = self.parent()
##         if isinstance(y, (int, long, Integer)):
##             z_elt = dict(self.__monomial_coefficients)
##             e = A.monoid()(1)
##             if e in z_elt:
##                 z_elt[e] += A.base_ring()(y)
##             else:
##                 z_elt[e] = A.base_ring()(y)
##             z = A(0)
##             z.__monomial_coefficients = z_elt
##             return z
##         if not isinstance(y, FreeAlgebraElement) or not A == y.parent():
##             raise TypeError, "Argument y (= %s) is of the wrong type."%y
        z_elt = dict(self.__monomial_coefficients)
        for m, c in y.__monomial_coefficients.iteritems():
            if m in z_elt:
                cm = z_elt[m] + c
                if cm == 0:
                    del z_elt[m]
                else:
                    z_elt[m] = cm
            else:
                z_elt[m] = c
        z = A(0)
        z.__monomial_coefficients = z_elt
        return z

    def _neg_(self):
        """
        Return negation of self

        EXAMPLES::

            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: -(x+y)    # indirect doctest
            -x - y
        """
        y = self.parent()(0)
        y_elt = {}
        for m, c in self.__monomial_coefficients.iteritems():
            y_elt[m] = -c
        y.__monomial_coefficients = y_elt
        return y

    def _sub_(self, y):
        """
        Return self minus y (another free algebra element with the same
        parents)

        EXAMPLES::

            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: x - y    # indirect doctest
            x - y
        """
        A = self.parent()
##         if isinstance(y, (int, long, Integer)):
##             z_elt = dict(self.__monomial_coefficients)
##             e = A.monoid()(1)
##             if e in z_elt:
##                 z_elt[e] += A.base_ring()(-y)
##             else:
##                 z_elt[e] = A.base_ring()(-y)
##             z = A(0)
##             z.__monomial_coefficients = z_elt
##             return z
##         if not isinstance(y, FreeAlgebraElement) or not A == y.parent():
##             raise TypeError, "Argument y (= %s) is of the wrong type."%y
        z_elt = dict(self.__monomial_coefficients)
        for m, c in y.__monomial_coefficients.iteritems():
            if m in z_elt:
                cm = z_elt[m] - c
                if cm == 0:
                    del z_elt[m]
                else:
                    z_elt[m] = cm
            else:
                z_elt[m] = -c
        z = A(0)
        z.__monomial_coefficients = z_elt
        return z

    def _mul_(self, y):
        """
        Return product of self and y (another free algebra element with the
        same parents)

        EXAMPLES::

            sage: A.<x,y,z>=FreeAlgebra(ZZ,3)
            sage: (x+y+x*y)*(x+y+1)    # indirect doctest
            x + y + x^2 + 2*x*y + y*x + y^2 + x*y*x + x*y^2
        """
        A = self.parent()
        z_elt = {}
        for mx, cx in self.__monomial_coefficients.iteritems():
            for my, cy in y.__monomial_coefficients.iteritems():
                key = mx*my
                if key in z_elt:
                    z_elt[key] += cx*cy
                else:
                    z_elt[key] = cx*cy
        z = A(0)
        z.__monomial_coefficients = z_elt
        return z

    def to_pbw_basis(self):
        """
        Return ``self`` in the Poincare-Birkhoff-Witt (PBW) basis.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(ZZ, 3)
            sage: p = x^2*y + 3*y*x + 2
            sage: p.to_pbw_basis()
            2*PBW[1] + 3*PBW[y]*PBW[x] + PBW[x^2*y] + PBW[x*y]*PBW[x] + PBW[y]*PBW[x]^2
        """
        return self.parent().pbw_element(self)

