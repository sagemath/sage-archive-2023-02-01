"""
Polynomial Compilers

AUTHOR:
    -- Tom Boothby

"""
################################################################################
#       Copyright (C) 2007 Tom Boothby <boothby@u.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################

from sage.misc.sagex_ds cimport BinaryTree
#from math import ceil, sqrt, log, floor
#from sage.rings.integer import Integer

cdef class CompiledPolynomialFunction:
    """
      Builds a reasonably optimized directed acyclic graph representation
    for a given polynomial.  A CompiledPolynomialFunction is callable from
    python, though it is a little faster to call the eval function from
    pyrex.

      This class is not intended to be called by a user, rather, it is
    intended to improve the performance of immutable polynomial objects.

      TODO:
        [ ] Recursive calling
        [ ] Faster casting of coefficients / argument
        [ ] Multivariate polynomials
        [ ] SageX implementation of Pippenger's Algorithm that doesn't
            depend heavily upon dicts.
        [ ] Computation of parameter sequence suggested by Pippenger
        [ ] Univariate exponentiation can use Brauer's method to improve
            extremely sparse polynomials of very high degree
    """


    def __init__(self, coeffs, method='quick'):
        cdef generic_pd max_gap
        self._coeffs = coeffs
        gaps, dag = self._parse_structure()
        max_gap = <generic_pd>(gaps.get_max())
        if max_gap.label > 1:
            if method == 'quick':
                self._fill_gaps_quick(gaps)
            elif method == 'pippenger':
                raise NotImplementedError, "Implementation of Pippenger's Algorithm is not ready for prime time."
                self._fill_gaps_pippenger()
            else:
                raise RuntimeError, "Method '%s' not supported."

        self._dag = dag.nodummies()

    cdef object eval(CompiledPolynomialFunction self, object x):
        cdef object temp
        try:
            pd_eval(self._dag, x, self._coeffs)
            temp = self._dag.value
            pd_clean(self._dag)
            return temp
        except:
            raise RuntimeError, "Polynomial evaluation error in val()!"

    cdef object _parse_structure(CompiledPolynomialFunction self):
        cdef BinaryTree gaps
        cdef int d, i
        cdef generic_pd s

        s = univar_pd()
        gaps = BinaryTree()
        gaps.insert(1, s)

        d = len(self._coeffs)-1

        s = coeff_pd(d)
        gap_width = 0

        d-= 1
        while d > 0:
            gap_width += 1
            if self._coeffs[d] != 0:
                s = abc_pd(s, self._get_gap(gaps, gap_width), d)
                gap_width = 0
            d-=1

        gap_width += 1
        if self._coeffs[0] != 0:
            s = abc_pd(s, self._get_gap(gaps, gap_width), 0)
        else:
            s = mul_pd(s, self._get_gap(gaps, gap_width))

        return gaps, s

    cdef generic_pd _get_gap(CompiledPolynomialFunction self, BinaryTree gaps, int gap):
        cdef generic_pd g
        g = gaps.get(gap)
        if g is not None:
            return g
        else:
            g = dummy_pd(gap)
            gaps.insert(gap,g)
            return g

    cdef void _fill_gaps_quick(CompiledPolynomialFunction self, BinaryTree gaps):
        cdef int m,n,k,r,half
        cdef generic_pd T,M,H
        cdef dummy_pd N
        T = gaps.pop_max()
        while T is not None:
            M = gaps.get_max()
            if M is None:
                return
            else:
                N = <dummy_pd>T
            m = M.label
            n = N.label
            k = n/m
            r = n%m
            half = n/2

            found = False
            if n % 2 == 0 and m >= half:
                H = gaps.get(half)
                if H is not None:
                    N.set_sqr(H)
                    found = True

            if not found:
                if r > 0:
                    N.set_mul(self._get_gap(gaps, m*k), self._get_gap(gaps, r))
                elif k % 2 == 0:
                    N.set_sqr(self._get_gap(gaps, m*k/2))
                else:
                    N.set_mul(self._get_gap(gaps, m*(k-1)), M)

            T = gaps.pop_max()





cdef inline void pd_eval(generic_pd pd, object vars, object coeffs):
    if pd.value is None:
        pd.eval(vars, coeffs)
    pd.hits += 1

cdef inline void pd_clean(generic_pd pd):
    if pd.hits >= pd.refs:
        pd.value = None
        pd.hits = 0

cdef class generic_pd:
    def __init__(generic_pd self):
        self.value = None
        self.hits = 0
        self.refs = 0
        self.label = -1

    cdef void eval(generic_pd self, object vars, object coeffs):
        raise NotImplementedError

    cdef generic_pd nodummies(generic_pd self):
        return self

cdef class dummy_pd(generic_pd):
    def __init__(dummy_pd self, int label):
        self.label = label

    cdef void set_sqr(dummy_pd self, generic_pd other):
        self.link = sqr_pd(other)

    cdef void set_mul(dummy_pd self, generic_pd left, generic_pd right):
        self.link = mul_pd(left, right)

    cdef generic_pd nodummies(dummy_pd self):
        #sorry guys, this is my stop
        self.link.refs = self.refs
        return self.link

cdef class var_pd(generic_pd):
    def __init__(var_pd self, int index):
        generic_pd.__init__(self)
        self.index = index
    cdef void eval(var_pd self, object vars, object coeffs):
        self.value = vars[self.index]

cdef class univar_pd(generic_pd):
    def __init__(univar_pd self):
        generic_pd.__init__(self)
        self.label = 1
    cdef void eval(univar_pd self, object var, object coeffs):
        self.value = var

cdef class coeff_pd(generic_pd):
    def __init__(coeff_pd self, int index):
        generic_pd.__init__(self)
        self.index = index
    cdef void eval(coeff_pd self, object vars, object coeffs):
        self.value = coeffs[self.index]

cdef class unary_pd(generic_pd):
    def __init__(unary_pd self, generic_pd operand):
        generic_pd.__init__(self)
        self.operand = operand
        self.operand.refs += 1

    def nodummies(unary_pd self):
        self.operand = self.operand.nodummies()


cdef class sqr_pd(unary_pd):
    cdef void eval(sqr_pd self, object vars, object coeffs):
        pd_eval(self.operand, vars, coeffs)
        self.value = self.operand.value * self.operand.value
        pd_clean(self.operand)

cdef class pow_pd(unary_pd):
    def __init__(unary_pd self, generic_pd base, object exponent):
        unary_pd.__init__(self, base)
        self.exponent = exponent

    cdef void eval(pow_pd self, object vars, object coeffs):
        pd_eval(self.operand, vars, coeffs)
        self.value = self.operand.value ** self.exponent
        pd_clean(self.operand)



cdef class binary_pd(generic_pd):
    def __init__(binary_pd self, generic_pd left, generic_pd right):
        generic_pd.__init__(self)
        self.left = left
        self.right = right
        self.left.refs+= 1
        self.right.refs+= 1

    def nodummies(binary_pd self):
        self.left = self.left.nodummies()
        self.right = self.right.nodummies()

cdef class add_pd(binary_pd):
    cdef void eval(add_pd self, object vars, object coeffs):
        pd_eval(self.left, vars, coeffs)
        pd_eval(self.right, vars, coeffs)
        self.value = self.left.value + self.right.value
        pd_clean(self.left)
        pd_clean(self.right)

cdef class mul_pd(binary_pd):
    cdef void eval(mul_pd self, object vars, object coeffs):
        pd_eval(self.left, vars, coeffs)
        pd_eval(self.right, vars, coeffs)
        self.value = self.left.value * self.right.value
        pd_clean(self.left)
        pd_clean(self.right)

cdef class abc_pd(binary_pd):
    def __init__(abc_pd self, generic_pd left, generic_pd right, int index):
        binary_pd.__init__(self, left, right)
        self.index = index

    cdef void eval(abc_pd self, object vars, object coeffs):
        pd_eval(self.left, vars, coeffs)
        pd_eval(self.right, vars, coeffs)
        self.value = self.left.value * self.right.value + coeffs[self.index]
        pd_clean(self.left)
        pd_clean(self.right)




