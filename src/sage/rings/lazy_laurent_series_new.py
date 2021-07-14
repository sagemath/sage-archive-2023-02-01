r"""
Lazy Laurent Series

A lazy Laurent series is a Laurent series whose coefficients are computed as
demanded or needed. Unlike the usual Laurent series in Sage, lazy Laurent
series do not have precisions because a lazy Laurent series knows (can be
computed, lazily) all its coefficients.

EXAMPLES:

Generating functions are Laurent series over the integer ring::

    sage: L.<z> = LLSRing(ZZ)

This defines the generating function of Fibonacci sequence::

    sage: def coeff(s, i):
    ....:     if i in [0, 1]:
    ....:         return 1
    ....:     else:
    ....:         return s.coefficient(i - 1) + s.coefficient(i - 2)
    sage: f = L.series(coeff, True, approximate_valuation=0); f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
    sage: f = L.series(coeff, False, approximate_valuation=0); f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...

The 100th element of Fibonacci sequence can be obtained from the generating
function::

    sage: f.coefficient(100)
    573147844013817084101

Coefficients are computed and cached only when necessary::

    sage: f._aux._cache[100]
    573147844013817084101
    sage: f._aux._cache[101]
    Traceback (most recent call last):
    ...
    KeyError: 101

You can do arithmetic with lazy power series::

    sage: f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
    sage: f^-1
    1 - z - z^2 + ...
    sage: f + f^-1
    2 + z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
    sage: g = (f + f^-1)*(f - f^-1); g
    4*z + 6*z^2 + 8*z^3 + 19*z^4 + 38*z^5 + 71*z^6 + ...

You may need to change the base ring::

    sage: h = g.change_ring(QQ)
    sage: h.parent()
    Lazy Laurent Series Ring in z over Rational Field
    sage: h
    4*z + 6*z^2 + 8*z^3 + 19*z^4 + 38*z^5 + 71*z^6 + ...
    sage: h^-1
    1/4*z^-1 - 3/8 + 1/16*z - 17/32*z^2 + 5/64*z^3 - 29/128*z^4 + 165/256*z^5 + ...
    sage: _.valuation()
    -1

AUTHORS:

- Kwankyu Lee (2019-02-24): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from .infinity import infinity
from sage.structure.element import ModuleElement
from .integer_ring import ZZ
from sage.structure.richcmp import op_EQ, op_NE
from sage.arith.power import generic_power


class LLS(ModuleElement):
    r"""
    Return a lazy Laurent series.

    INPUT:

    - ``coefficient_function`` -- Python function that computes coefficients

    - ``issparse`` -- Boolean that determines whether the implementation is sparse or dense

    - ``approximate_valuation`` -- integer; approximate valuation of the series

    - ``constant`` -- either ``None`` or pair of an element of the base ring and an integer

    Let the coefficient of index `i` mean the coefficient of the term of the
    series with exponent `i`.

    Python function ``coefficient`` returns the value of the coefficient of
    index `i` from input.

    Let ``approximate_valuation`` be `n`. All coefficients of index below `n` are zero.  If
    ``constant`` is ``None``, then the ``coefficient`` function is responsible
    to compute the values of all coefficients of index `\ge n`. If ``constant``
    is a pair `(c,m)`, then the ``coefficient`` function is responsible to
    compute the values of all coefficients of index `\ge n` and `< m` and all
    the coefficients of index `\ge m` is the constant `c`.

    EXAMPLES::

        sage: L = LLSRing(ZZ, 'z')
        sage: L.series(lambda i: i, is_sparse=True, approximate_valuation=-3, constant=(-1,3))
        -3*z^-3 - 2*z^-2 - z^-1 + z + 2*z^2 - z^3 - z^4 - z^5 + ...
        sage: L.series(lambda i: i, is_sparse=False, approximate_valuation=-3, constant=(-1,3))
        -3*z^-3 - 2*z^-2 - z^-1 + z + 2*z^2 - z^3 - z^4 - z^5 + ...

    ::

        sage: def coeff(s, i):
        ....:     if i in [0, 1]:
        ....:         return 1
        ....:     else:
        ....:         return s.coefficient(i - 1) + s.coefficient(i - 2)
        sage: f = L.series(coeff, True, approximate_valuation=0); f
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: f.coefficient(100)
        573147844013817084101

    Lazy Laurent series is picklable::

        sage: z = L.gen()
        sage: f = 1/(1 - z - z^2)
        sage: f
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: g = loads(dumps(f))
        sage: g
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: g == f
        True
    """
    def __init__(self, parent, aux):
        """
        Initialize the series.

        TESTS::

            sage: L = LLSRing(GF(2), 'z')
            sage: z = L.gen()
            sage: TestSuite(z).run()
        """
        ModuleElement.__init__(self, parent)
        self._aux = aux

    def __getitem__(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = z/(1 - 2*z^3)
            sage: [f[n] for n in range(20)]
            [0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 16, 0, 0, 32, 0, 0, 64]
            sage: M = L.series(lambda n: n, True, 0)                                                      
            sage: [M[n] for n in range(20)]                                                               
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: M = L.series(lambda n: n, False, 0)                                                      
            sage: [M[n] for n in range(20)]                                                               
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

        """
        return self.base_ring()(self._aux[n])

    def _mul_(self, other):
        """
        Return the product of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
            sage: M = L.series(lambda n: n, True, 0)                                                      
            sage: M                                                                                       
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = M * (1 - M)                                                                         
            sage: N                                                                                       
            z + z^2 - z^3 - 6*z^4 - 15*z^5 - 29*z^6 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L.series(lambda n: 1, False, 0); N                                                  
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M * N; P                                                                            
            z + 3*z^2 + 6*z^3 + 10*z^4 + 15*z^5 + 21*z^6 + ...
        """
        P = self.parent()
        return P.element_class(P, LLS_mul(self._aux, other._aux))
    
    def _add_(self, other):
        """
        Return the sum of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
            sage: z + z                                                                                   
            2*z
            sage: z^2 + 3*z^2                                                                             
            4*z^2
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L.series(lambda n: 1, True, 0); N                                                   
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M + N; P                                                                            
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L.series(lambda n: 1, False, 0); N                                                  
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M + N; P                                                                            
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + ...
        """
        P = self.parent()
        return P.element_class(P, LLS_add(self._aux, other._aux))

    def _div_(self, other):
        """
        Return ``self`` divided by ``other``.

        INPUT:

        - ``other`` -- nonzero series

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: z/(1 - z)
            z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L.series(lambda n: 1, False, 0); N                                                  
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M / N; P                                                                            
            z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + ...
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L.series(lambda n: 1, True, 0); N                                                   
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M / N; P                                                                            
            z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + ...
        """
        P = self.parent()
        return P.element_class(P, LLS_div(self._aux, other._aux))
    
    def _rmul_(self, scalar):
        """
        Return the scalar multiplication of this series by ``scalar``.

        INPUT:

        - ``scalar`` -- an element of the base ring

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: 2*z
            2*z
            sage: -1*z
            -z
            sage: 0*z
            0
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...     
            sage: M._rmul_(3)                                                                                    
            3*z + 6*z^2 + 9*z^3 + 12*z^4 + 15*z^5 + 18*z^6 + ...        
            sage: N = L.series(lambda n: 1, False, 0); N                                                         
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...       
            sage: N._rmul_(4)                                                                                    
            4 + 4*z + 4*z^2 + 4*z^3 + 4*z^4 + 4*z^5 + 4*z^6 + ...       
        """
        P = self.parent()
        return P.element_class(P, LLS_rmul(self._aux, scalar))
    
    def _sub_(self, other):
        """
        Return the series of this series minus ``other`` series.

        INPUT:

        - ``other`` -- other series

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: z - z
            0
            sage: 3*z - 2*z                                                                               
            z
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L.series(lambda n: 1, False, 0); N                                                  
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M - N; P                                                                            
            -1 + z^2 + 2*z^3 + 3*z^4 + 4*z^5 + 5*z^6 + ...
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L.series(lambda n: 1, True, 0); N                                                   
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M - N; P                                                                            
            -1 + z^2 + 2*z^3 + 3*z^4 + 4*z^5 + 5*z^6 + ...
        """
        P = self.parent()
        return P.element_class(P, LLS_sub(self._aux, other._aux))
    
    def _neg_(self):
        """
        Return the negative of this series.

        TESTS::

            sage: L = LLSRing(ZZ, 'z')
            sage: z = L.gen()
            sage: -(1 - z)
            -1 + z
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: P = -M; P                                                                               
            -z - 2*z^2 - 3*z^3 - 4*z^4 - 5*z^5 - 6*z^6 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: P = -M; P                                                                               
            -z - 2*z^2 - 3*z^3 - 4*z^4 - 5*z^5 - 6*z^6 + ...
            sage: -(z^2 + 3*z - 4*z^3)                                                                    
            -3*z - z^2 + 4*z^3
        """
        P = self.parent()
        return P.element_class(P, LLS_neg(self._aux))
    
    def __invert__(self):
        """
        Return the multiplicative inverse of the element.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: ~(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: P = ~M; P                                                                               
            z^-1 - 2 + z + ...
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: P = ~M; P                                                                               
            z^-1 - 2 + z + ...
        """
        P = self.parent()
        return P.element_class(P, LLS_inv(self._aux))
    
    def coefficient(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: L = LLSRing(ZZ, 'z')
            sage: def g(s, i):
            ....:     if i == 0:
            ....:         return 1
            ....:     else:
            ....:         return sum(s.coefficient(j)*s.coefficient(i - 1 -j) for j in [0..i-1])
            sage: e = L.series(g, True, approximate_valuation=0)
            sage: e.coefficient(10)
            16796
            sage: e
            1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...

        TESTS::

            sage: def g(s, i):
            ....:     if i == 0:
            ....:         return 1
            ....:     else:
            ....:         return sum(s.coefficient(j)*s.coefficient(i - 1 - j) for j in [0..i-1])
            ....:
            sage: L = LLSRing(ZZ, 'z')
            sage: e = L.series(g, False, approximate_valuation = 0)
            sage: e
            1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...
            sage: e._aux._cache
            [1, 1, 2, 5, 14, 42, 132]
            sage: e.coefficient(10)
            16796
            sage: e._aux._cache
            [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]
            sage: M = L.series(lambda n: n^2, False, 0); M                                                
            z + 4*z^2 + 9*z^3 + 16*z^4 + 25*z^5 + 36*z^6 + ...
            sage: M._aux._cache                                                                           
            [0, 1, 4, 9, 16, 25, 36]
            sage: M.coefficient(9)                                                                        
            81
            sage: M._aux._cache                                                                           
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
            sage: M = L.series(lambda n: n^2, True, 0); M                                                 
            z + 4*z^2 + 9*z^3 + 16*z^4 + 25*z^5 + 36*z^6 + ...
            sage: M._aux._cache                                                                           
            {0: 0, 1: 1, 2: 4, 3: 9, 4: 16, 5: 25, 6: 36}
            sage: M.coefficient(10)                                                                       
            100
            sage: M._aux._cache                                                                           
            {0: 0, 1: 1, 2: 4, 3: 9, 4: 16, 5: 25, 6: 36, 10: 100}
        """
        return self.__getitem__(n)
    
    def apply_to_coefficients(self, newfunction, ring):
        """
        Return the series with ``newfunction`` applied to each coefficient of this series.

        INPUT:

        - ``newfunction`` -- Python function

        - ``ring`` -- Base Ring of the series

        Python function ``newfunction`` returns a new coefficient for input coefficient.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: s = z/(1 - 2*z)
            sage: t = s.apply_to_coefficients(lambda c: c + 1, L)
            sage: s
            z + 2*z^2 + 4*z^3 + 8*z^4 + 16*z^5 + 32*z^6 + 64*z^7 + ...
            sage: t
            2*z + 3*z^2 + 5*z^3 + 9*z^4 + 17*z^5 + 33*z^6 + 65*z^7 + ...
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = M.apply_to_coefficients(lambda c: c + 1, L); N                                         
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = M.apply_to_coefficients(lambda c: c + 1, L); N                                         
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + ...
        """
        P = self.parent()
        return P.element_class(P, LLS_apply_coeff(self._aux, newfunction, ring))
    
    def change_ring(self, ring):
        """
        Return this series with coefficients converted to elements of ``ring``.

        INPUT:

        - ``ring`` -- a ring

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: s = 2 + z
            sage: t = s.change_ring(QQ)
            sage: t^-1
            1/2 - 1/4*z + 1/8*z^2 - 1/16*z^3 + 1/32*z^4 - 1/64*z^5 + 1/128*z^6 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = M.change_ring(QQ)                                                                   
            sage: N.parent()                                                                              
            Lazy Laurent Series Ring in z over Rational Field
            sage: M.parent()                                                                              
            Lazy Laurent Series Ring in z over Integer Ring
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M.parent()                                                                              
            Lazy Laurent Series Ring in z over Integer Ring
            sage: N = M.change_ring(QQ)                                                                   
            sage: N.parent()                                                                              
            Lazy Laurent Series Ring in z over Rational Field
            sage: M ^-1                                                                                   
            z^-1 - 2 + z + ...
        """
        from .lazy_laurent_series_ring_new import LLSRing
        Q = LLSRing(ring, names=self.parent().variable_name())
        return Q.element_class(Q, self._aux)

    def truncate(self, d):
        """
        Return this series with its terms of degree >= ``d`` truncated.

        INPUT:

        - ``d`` -- integer

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: alpha = 1/(1-z)
            sage: alpha
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: beta = alpha.truncate(5)
            sage: beta
            1 + z + z^2 + z^3 + z^4
            sage: alpha - beta
            z^5 + z^6 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M.truncate(4)                                                                           
            z + 2*z^2 + 3*z^3
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M.truncate(4)                                                                           
            z + 2*z^2 + 3*z^3
        """
        P = self.parent()
        return P.element_class(P, LLS_trunc(self._aux, d))

    def __pow__(self, n):
        """
        Return the `n`-th power of the series.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: (1 - z)^-1
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: (1 - z)^0
            1
            sage: (1 - z)^3
            1 - 3*z + 3*z^2 - z^3
            sage: (1 - z)^-3
            1 + 3*z + 6*z^2 + 10*z^3 + 15*z^4 + 21*z^5 + 28*z^6 + ...
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M ^ 2                                                                                   
            z^2 + 4*z^3 + 10*z^4 + 20*z^5 + 35*z^6 + ...
            sage: M = L.series(lambda n: n, False, 0); M                                                  
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M ^ 2                                                                                   
            z^2 + 4*z^3 + 10*z^4 + 20*z^5 + 35*z^6 + ...
        """
        if n == 0:
            return self.parent().one()

        return generic_power(self, n)      
    
    def approximate_series(self, prec, name=None):
        """
        Return the Laurent series with absolute precision ``prec`` approximated
        from this series.

        INPUT:

        - ``prec`` -- an integer

        - ``name`` -- name of the variable; if it is ``None``, the name of the variable
          of the series is used

        OUTPUT: a Laurent series with absolute precision ``prec``

        TESTS::

            sage: L = LLSRing(ZZ, 'z')
            sage: z = L.gen()
            sage: f = (z - 2*z^3)^5/(1 - 2*z)
            sage: f
            z^5 + 2*z^6 - 6*z^7 - 12*z^8 + 16*z^9 + 32*z^10 - 16*z^11 + ...
            sage: g = f.approximate_series(10)
            sage: g
            z^5 + 2*z^6 - 6*z^7 - 12*z^8 + 16*z^9 + O(z^10)
            sage: g.parent()
            Power Series Ring in z over Integer Ring
            sage: h = (f^-1).approximate_series(3)
            sage: h
            z^-5 - 2*z^-4 + 10*z^-3 - 20*z^-2 + 60*z^-1 - 120 + 280*z - 560*z^2 + O(z^3)
            sage: h.parent()
            Laurent Series Ring in z over Integer Ring
        """
        S = self.parent()

        if name is None:
            name = S.variable_name()

        if self.valuation() < 0:
            from sage.rings.all import LaurentSeriesRing
            R = LaurentSeriesRing(S.base_ring(), name=name)
            n = self.valuation()
            return R([self[i] for i in range(n, prec)], n).add_bigoh(prec)
        else:
            from sage.rings.all import PowerSeriesRing
            R = PowerSeriesRing(S.base_ring(), name=name)
            return R([self[i] for i in range(prec)]).add_bigoh(prec)

    def prec(self):
        """
        Return the precision of the series, which is infinity.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z)
            sage: f.prec()
            +Infinity
        """
        return infinity
    
    def polynomial(self, degree=None, name=None):
        """
        Return the polynomial or Laurent polynomial if the series is actually so.

        INPUT:

        - ``degree`` -- ``None`` or an integer

        - ``name`` -- name of the variable; if it is ``None``, the name of the variable
          of the series is used

        OUTPUT: a Laurent polynomial if the valuation of the series is negative or
        a polynomial otherwise.

        If ``degree`` is not ``None``, the terms of the series of degree
        greater than ``degree`` are truncated first. If ``degree`` is ``None``
        and the series is not a polynomial or a Laurent polynomial, a
        ``ValueError`` is raised.

        EXAMPLES::

            sage: L = LLSRing(ZZ, 'z')
            sage: f = L.series([1,0,0,2,0,0,0,3], True, 5); f
            z^5 + 2*z^8 + 3*z^12
            sage: f.polynomial()
            3*z^12 + 2*z^8 + z^5

        TESTS::

            sage: g = L.series([1,0,0,2,0,0,0,3], True, -5); g
            z^-5 + 2*z^-2 + 3*z^2
            sage: g.polynomial()
            z^-5 + 2*z^-2 + 3*z^2
            sage: z = L.gen()
            sage: f = (1 + z)/(z^3 - z^5)
            sage: f
            z^-3 + z^-2 + z^-1 + 1 + z + z^2 + z^3 + ...
            sage: f.polynomial(5)
            z^-3 + z^-2 + z^-1 + 1 + z + z^2 + z^3 + z^4 + z^5
            sage: f.polynomial(0)
            z^-3 + z^-2 + z^-1 + 1
            sage: f.polynomial(-5)
            0
            sage: M = L.series(lambda n: n^2, False, 0)                                                   
            sage: M.polynomial(3)                                                                         
            9*z^3 + 4*z^2 + z
            sage: M = L.series(lambda n: n^2, True, 0)                                                    
            sage: M.polynomial(5)                                                                         
            25*z^5 + 16*z^4 + 9*z^3 + 4*z^2 + z
        """
        if degree is None:
            if self._aux._constant is None or not self._aux._constant[0].is_zero():
                raise ValueError("not a polynomial")
            m = self._aux._constant[1]
        else:
            m = degree + 1

        S = self.parent()

        if name is None:
            name = S.variable_name()

        if self.valuation() < 0:
            from sage.rings.all import LaurentPolynomialRing
            R = LaurentPolynomialRing(S.base_ring(), name=name)
            n = self.valuation()
            return R([self[i] for i in range(n, m)]).shift(n)
        else:
            from sage.rings.all import PolynomialRing
            R = PolynomialRing(S.base_ring(), name=name)
            return R([self[i] for i in range(m)])
    
    def valuation(self):
        """
        Return the valuation of the series.

        This method determines the valuation of the series by looking for a
        nonzero coefficient. Hence if the series happens to be zero, then it
        may run forever.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: s = 1/(1 - z) - 1/(1 - 2*z)
            sage: s.valuation()
            1
            sage: t = z - z
            sage: t.valuation()
            +Infinity
            sage: M = L.series(lambda n: n^2, True, 0)                                                    
            sage: M.valuation()                                                                           
            1
            sage: M = L.series(lambda n: n^2, False, 0)                                                   
            sage: M.valuation()                                                                           
            1
        """
        return self._aux.valuation()

    def _repr_(self):
        """
        Return the string representation of this Laurent series.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: -1/(1 + 2*z)
            -1 + 2*z - 4*z^2 + 8*z^3 - 16*z^4 + 32*z^5 - 64*z^6 + ...
        """
        atomic_repr = self.base_ring()._repr_option('element_is_atomic')
        X = self.parent().variable_name()
        n = self._aux._approximate_valuation

        if self._aux._constant is None:
            m = n + 7  # long enough
        elif self._aux._constant[0] != 0:
            m = self._aux._constant[1] + 3
        else:
            m = self._aux._constant[1]

        s = ' '
        first = True
        while n < m:
            x = repr(self._aux[n])
            if x != '0':
                if not first:
                    s += ' + '
                if not atomic_repr and n > 0 and (x[1:].find('+') != -1 or x[1:].find('-') != -1):
                    x = '({})'.format(x)
                if n > 1 or n < 0:
                    var = '*{}^{}'.format(X, n)
                elif n == 1:
                    var = '*{}'.format(X)
                else:  # n == 0
                    var = ''
                s += '{}{}'.format(x, var)
                first = False
            n += 1

        s = s.replace(" + -", " - ").replace(" 1*", " ").replace(" -1*", " -")[1:]

        if not s:  # zero series
            s = '0'

        if self._aux._constant is None or self._aux._constant[1] > m or self._aux._constant[0] != 0:
            s += ' + {}'.format('...')

        return s
    
    def _richcmp_(self, other, op):
        """
        Compare ``self` with ``other`` with respect to the comparison operator ``op``.

        Equality is verified if the corresponding coefficients of both series
        can be checked for equality without computing coefficients
        indefinitely.  Otherwise an exception is raised to declare that
        equality is not decidable.

        Inequality is not defined for lazy Laurent series.

        TESTS::

            sage: L.<z> = LLSRing(QQ)
            sage: z + z^2 == z^2 + z
            True
            sage: z + z^2 != z^2 + z
            False
            sage: z + z^2 > z^2 + z
            False
            sage: z + z^2 < z^2 + z
            False
        """
        if op is op_EQ:
            if self._aux._constant is None:
                if other._aux._constant is None:
                    # Implement the checking of the caches here.
                    n = min(self._aux._approximate_valuation, other._aux._approximate_valuation)
                    m = max(self._aux._approximate_valuation, other._aux._approximate_valuation)
                    for i in range(n, m):
                        if self[i] != other[i]:
                            return False
                    if self._aux == other._aux:
                        return True
                    raise ValueError("undecidable as lazy Laurent series")
                else:
                    raise ValueError("undecidable as lazy Laurent series")
            elif other._aux._constant is None:
                raise ValueError("undecidable as lazy Laurent series")

            sc, sm = self._aux._constant
            oc, om = other._aux._constant

            if sc != oc:
                return False

            n = self._aux._approximate_valuation
            m = max(sm, om)

            for i in range(n, m):
                if self[i] != other[i]:
                    return False

            return True

        if op is op_NE:
            return not (self == other)

        return False
    
    def __hash__(self):
        """
        Return the hash of ``self``

        TESTS::

            sage: L = LLSRing(ZZ, 'z')
            sage: f = L.series([1,2,3,4], True, -5)
            sage: g = (1 + f)/(1 - f)^2
            sage: {g: 1}
            {z^5 - 2*z^6 + z^7 + 5*z^9 - 11*z^10 + z^11 + ...: 1}
        """
        return hash((type(self), self._aux._approximate_valuation, self._aux._constant))
    
    def __bool__(self):
        """
        Test whether ``self`` is not zero.

        TESTS::

            sage: L.<z> = LLSRing(GF(2))
            sage: (z-z).is_zero()
            True
            sage: f = 1/(1 - z)
            sage: f.is_zero()
            False
            sage: M = L.series(lambda n: n, True, 0); M                                                   
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M.is_zero()                                                                             
            False
        """
        if self._aux._constant is None:
            for a in self._aux._cache:
                if a:
                    return True
            if self[self._aux._approximate_valuation]:
                return True
            raise ValueError("undecidable as lazy Laurent series")

        sc, sm = self._aux._constant

        if sc:
            return True

        for i in range(self._aux._approximate_valuation, sm):
            if self[i]:
                return True

        return False


class LLS_aux():
    def __init__(self, is_sparse, approximate_valuation, constant=None):
        self._approximate_valuation = approximate_valuation
        self._constant = constant
        self._is_sparse = is_sparse
        
        if self._is_sparse:
            self._cache = dict()  # cache of known coefficients
        else:
            self._cache = list()
            self._offset = approximate_valuation
            self._iter = self.iterate_coefficients()

    def __getitem__(self, n):
        if self._approximate_valuation is infinity:
            return ZZ.zero()
        elif n < self._approximate_valuation:
            return ZZ.zero()
        elif self._constant is not None and n >= self._constant[1]:
            return self._constant[0]

        if self._is_sparse:
            try:
                c = self._cache[n]
            except KeyError:
                c = self.get_coefficient(n)
                self._cache[n] = c
        else:
            i = n - self._offset
            if i >= len(self._cache):
                a = len(self._cache) + self._offset
                # it is important to extend by generator:
                # self._coefficient_function might recurse, and
                # thereby extend the cache itself, too
                self._cache.extend(next(self._iter) for _ in range(a, n+1))
            c = self._cache[i]

        return c
    
    def valuation(self):
        if self._constant is not None:
            n = self._approximate_valuation
            m = self._constant[1]
            while n <= m:
                if self[n] != 0:
                    self._approximate_valuation = n
                    return n
                n += 1
            return infinity

        if self._is_sparse:
            n = self._approximate_valuation
            cache = self._cache
            while True:
                if n in cache:
                    if cache[n]:
                        self._approximate_valuation = n
                        return n
                    n += 1
                else:
                    if self[n] != 0:
                        self._approximate_valuation = n
                        return n
                    n += 1
        else:
            n = self._approximate_valuation
            cache = self._cache
            while True:
                if n - self._offset < len(cache):
                    if cache[n - self._offset]:
                        self._approximate_valuation = n
                        return n
                    n += 1
                else:
                    if self[n] != 0:
                        self._approximate_valuation = n
                        return n
                    n += 1
    
    # def __eq__(self, other):
    # Implement for sparse and dense variations separately.
    #     n = min(self._approximate_valuation, other._approximate_valuation)
    #     m = max(self._approximate_valuation, other._approximate_valuation)
    #     for i in range(n, m):
    #         if self[i] != other[i]:
    #             return False



class LLS_unary(LLS_aux):
    """
    Abstract base class for unary operators.

    INPUT:

    - ``series`` -- series upon which the operator operates

    """
    def __init__(self, series, *args, **kwargs):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = -1/(1 - z)
            sage: f
            -1 - z - z^2 - z^3 - z^4 - z^5 - z^6 + ...
            sage: loads(dumps(f)) == f
            True
        """
        self._series = series
        super().__init__(*args, **kwargs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = ~(1 - z)
            sage: {f: 1}
            {1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...: 1}
        """
        return hash((type(self), self._series))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: g = 1/(1 - z) + 1/(1 + z)
            sage: f == g
            True
            sage: f = ~(1 - z)
            sage: g = ~(1 - z)
            sage: f == g
            True
        """
        return isinstance(other, type(self)) and self._series == other._series


class LLS_binary(LLS_aux):
    
    def __init__(self, left, right, *args, **kwargs):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
            sage: f = 1/(1 - z) - 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
        """
        self._left = left
        self._right = right
        super().__init__(*args, **kwargs)
    
    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: {f: 1}
            {2 + 2*z^2 + 2*z^4 + 2*z^6 + ...: 1}
        """
        return hash((type(self), self._left, self._right))
    
    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: g = 1/(1 - z) + 1/(1 + z)
            sage: f == g
            True
        """
        if not isinstance(other, type(self)):
            return False
        return self._left == other._left and self._right == other._right


class LLS_binary_commutative(LLS_binary):

    def __hash__(self):
        return hash((type(self), frozenset([self._left, self._right])))

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        if self._left == other._left and self._right == other._right:
            return True
        if self._left == other._right and self._right == other._left:
            return True
        return False

class LLS_coefficient_function(LLS_aux):
    def __init__(self, coefficient_function, is_sparse, approximate_valuation, constant=None):
        self._coefficient_function = coefficient_function
        super().__init__(is_sparse, approximate_valuation, constant)

    def get_coefficient(self, n):
        return self._coefficient_function(n)

    def iterate_coefficients(self):
        n = self._offset
        while True:
            yield self._coefficient_function(n)
            n += 1
    
           
class LLS_mul(LLS_binary):
    """
    Operator for multiplication.

    We are assuming commutativity of the coefficient ring here. 
    """
    def __init__(self, left, right):
        a = left._approximate_valuation + right._approximate_valuation
        c = None
        if left._constant is not None and right._constant is not None:
            if left._constant[0] == 0 and right._constant[0] == 0:
                c = (left._constant[0], left._constant[1] + right._constant[1] - 1)

        if left._is_sparse != right._is_sparse:
            raise NotImplementedError
        super().__init__(left, right, left._is_sparse, a, c)

    def get_coefficient(self, n):
        c = ZZ.zero()
        for k in range(self._left._approximate_valuation,
                       n - self._right._approximate_valuation + 1):
            l = self._left[k]
            if l:
                c += l * self._right[n-k]
        return c

    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = ZZ.zero()
            for k in range(self._left._approximate_valuation,
                           n - self._right._approximate_valuation + 1):
                l = self._left[k]
                if l:
                    c += l * self._right[n-k]
            yield c
            n += 1


class LLS_add(LLS_binary):
    """
    Operator for addition.
    """
    def __init__(self, left, right):
        a = min(left._approximate_valuation, right._approximate_valuation)

        if left._constant is not None and right._constant is not None:
            c = (left._constant[0] + right._constant[0],
                 max(left._constant[1], right._constant[1]))
        else:
            c = None
        
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError
        
        super().__init__(left, right, left._is_sparse, a, c)
    
    def get_coefficient(self, n):
        c = ZZ.zero()
        c = self._left[n] + self._right[n]
        return c

    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = ZZ.zero()
            c = self._left[n] + self._right[n]
            yield c
            n += 1


class LLS_sub(LLS_binary):
    """
    Operator for subtraction.
    """
    def __init__(self, left, right):
        a = min(left._approximate_valuation, right._approximate_valuation)

        if left._constant is not None and right._constant is not None:
            c = (left._constant[0] - right._constant[0],
                 max(left._constant[1], right._constant[1]))
        else:
            c = None
        
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError
        
        super().__init__(left, right, left._is_sparse, a, c)
    
    def get_coefficient(self, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = (1 + z)*(1 - z)
            sage: f.coefficient(2)
            -1
        """
        c = ZZ.zero()
        c = self._left[n] - self._right[n]
        return c

    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = ZZ.zero()
            c = self._left[n] - self._right[n]
            yield c
            n += 1


class LLS_rmul(LLS_unary):
    """
    Operator for multiplying with a scalar.
    """
    def __init__(self, series, scalar):
        self._scalar = scalar

        a = series._approximate_valuation
        if series._constant is not None:
            c = (scalar * series._constant[0], series._constant[1])
        else:
            c = None
        
        super().__init__(series, series._is_sparse, a, c)
    
    def get_coefficient(self, n):
        c = self._series[n] * self._scalar
        return c

    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = self._series[n] * self._scalar
            yield c
            n += 1


class LLS_neg(LLS_unary):
    """
    Operator for negative of the series.
    """
    def __init__(self, series):
        a = series._approximate_valuation

        if series._constant is not None:
            c = (-series._constant[0], series._constant[1])
        else:
            c = None
        
        super().__init__(series, series._is_sparse, a, c)
    
    def get_coefficient(self, n):
        c = -1 * self._series[n]
        return c

    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = -1 * self._series[n]
            yield c
            n += 1


class LLS_inv(LLS_unary):
    """
    Operator for multiplicative inverse of the series.
    """
    def __init__(self, series):
        v = series.valuation()
        # self._constant can be refined
        super().__init__(series, series._is_sparse, -v, None)

        if v is infinity:
            raise ZeroDivisionError('cannot invert zero')

        # self._ainv = ~series[v]
        self._ainv = series[v].inverse_of_unit()
        self._zero = ZZ.zero()
        
    def get_coefficient(self, n):
        v = self._approximate_valuation
        if n == v:
            return self._ainv
        c = self._zero
        for k in range(v, n):
            c += self[k] * self._series[n - v - k]
        return -c * self._ainv

    def iterate_coefficients(self):
        n = self._offset
        while True:
            v = self._approximate_valuation
            if n == v:
                yield self._ainv
                n += 1
                continue
            c = self._zero
            for k in range(v, n):
                c += self[k] * self._series[n - v - k]
            yield -c * self._ainv
            n += 1


class LLS_apply_coeff(LLS_unary):
    """
        Return the series with ``function`` applied to each coefficient of this series.
    """
    def __init__(self, series, function, ring):
        self._function = function
        self._ring = ring
        a = series._approximate_valuation

        if series._constant:
            c = (function(series._constant[0]), series._constant[1])
        else:
            c = None

        super().__init__(series, series._is_sparse, a, c)

    def get_coefficient(self, n):
        try:
            c = self._ring(self._function(self._series[n]))
            return c
        except TypeError:
            raise ValueError("The coefficients are not in the base ring.")
    
    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = self._function(self._series[n])
            try:
                c = self._ring(self._function(self._series[n]))
                yield c
                n += 1
            except TypeError:
                raise ValueError("The coefficients are not in the base ring.")


class LLS_trunc(LLS_unary):
    """
        Return this series with its terms of degree >= ``d`` truncated.
    """
    def __init__(self, series, d):
        self._d = d
        self._zero = ZZ.zero()
        a = series._approximate_valuation
        c = (ZZ.zero(), d)
        super().__init__(series, series._is_sparse, a, c)
    
    def get_coefficient(self, n):
        if n <= self._d:
            c = self._series[n]
            return c
        else:
            c = self._zero
            return c
    
    def iterate_coefficients(self):
        n = self._offset
        while True:
            if n <= self._d:
                c = self._series[n]
            else:
                c = self._zero
            yield c
            n += 1


class LLS_div(LLS_binary):
    """
        Return ``self`` divided by ``other``.
    """
    def __init__(self, left, right):
        
        super().__init__(left, right, left._is_sparse, left._approximate_valuation, None)

        self._approximate_valuation = left.valuation() - right.valuation()

        lv = left.valuation()
        rv = right.valuation()
        self._lv = lv
        self._rv = rv
        # self._ainv = ~right[rv]
        self._ainv = right[rv].inverse_of_unit()
    
    def get_coefficient(self, n):
        lv = self._lv
        rv = self._rv
        if n == lv - rv:
            return self._left[lv]/self._right[rv]
        c = self._left[n + rv]
        for k in range(lv - rv, n):
            c -= self[k] * self._right[n + rv - k]
        return c * self._ainv
    
    def iterate_coefficients(self):
        n = self._offset
        lv = self._lv
        rv = self._rv
        while True:
            if n == lv - rv:
                yield self._left[lv]/self._right[rv]
                n += 1
                continue
            c = self._left[n + rv]
            for k in range(lv - rv, n):
                c -= self[k] * self._right[n + rv - k]
            yield c * self._ainv
            n += 1
