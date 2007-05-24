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
        self._coeffs = coeffs
        gaps, dag = self._parse_structure()
        if gaps.get_max().get_label() > 1:
            if method == 'quick':
                self._fill_gaps_quick(gaps)
            elif method == 'pippenger':
                raise NotImplementedError, "Implementation of Pippenger's Algorithm is not ready for prime time."
                self._fill_gaps_pippenger()
            else:
                raise RuntimeError, "Method '%s' not supported."

        self._dag = dag.accelerate()

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
        cdef poly_dag s

        s = univar_dag()
        s.set_label(1)
        gaps = BinaryTree()
        gaps.insert(1, s)

        d = len(self._coeffs)-1

        s = coeff_dag(d)
        gap_width = 0

        d-= 1
        while d > 0:
            gap_width += 1
            if self._coeffs[d] != 0:
                s = abc_dag(s, self._get_gap(gaps, gap_width), d)
                gap_width = 0
            d-=1

        gap_width += 1
        if self._coeffs[0] != 0:
            s = abc_dag(s, self._get_gap(gaps, gap_width), 0)
        else:
            s*= self._get_gap(gaps, gap_width)

        return gaps, s

    cdef poly_dag _get_gap(CompiledPolynomialFunction self, BinaryTree gaps, int gap):
        cdef poly_dag g
        g = gaps.get(gap)
        if g is not None:
            return g
        else:
            g = poly_dag()
            g.set_label(gap)
            gaps.insert(gap,g)
            return g

    cdef void _fill_gaps_quick(CompiledPolynomialFunction self, BinaryTree gaps):
        cdef int m,n,k,r,half
        cdef poly_dag N,M,H
        N = gaps.pop_max()
        while N is not None:
            M = gaps.get_max()
            if M is None:
                return
            m = M.get_label()
            n = N.get_label()
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

            N = gaps.pop_max()






cdef class poly_dag:
    def __init__(self, op=None, value=None, left=None, right=None, label = 0, index=0):
        self.dag = None
        self.label = label
        if op == pdCOEFF:
            self.dag = coeff_pd(index)
        elif op == pdSQR:
            self.dag = sqr_pd(left)
        elif op == pdMUL:
            self.dag = mul_pd(left,right)
        elif op == pdADD:
            self.dag = add_pd(left,right)
        elif op == pdPOW:
            self.dag = pow_pd(left, right)
        elif op == pdVAR:
            self.dag = var_pd(index)
        elif op == pdUNIVAR:
            self.dag = univar_pd()
        elif op == pdABC:
            self.dag = abc_pd(left, right, index)


    def set_label(self, label):
        self.label = label
    def get_label(self):
        return self.label

    def __add__(self,other):
        return add_dag(self,other);
    def __mul__(self,other):
        return mul_dag(self,other);
    def __sqr__(self):
        return sqr_dag(self)

    def value(self, vars, coeffs):
        return self.val(vars, coeffs)

    cdef object val(poly_dag self, object vars, object coeffs):
        cdef object temp
        try:
            pd_eval(self.dag, vars, coeffs)
            temp = self.dag.value
            pd_clean(self.dag)
            return temp
        except:
            raise RuntimeError, "Polynomial evaluation error in val()!"

    def set_mul(self, left, right):
        if self.dag is not None:
            raise RuntimeError, "Cannot re-cast a dag that's already been set"
        self.dag = mul_pd(left,right)

    def set_add(self,left, right):
        if self.dag is not None:
            raise RuntimeError, "Cannot re-cast a dag that's already been set"
        self.dag = add_pd(left, right)

    def set_sqr(self,left):
        if self.dag is not None:
            raise RuntimeError, "Cannot re-cast a dag that's already been set"
        self.dag = sqr_pd(left)

    def set_pow(self, base, pow):
        if self.dag is not None:
            raise RuntimeError, "Cannot re-cast a dag that's already been set"
        self.dag = pow_pd(base, pow)

    def clone(self, poly_dag other):
        if self.dag is not None:
            raise RuntimeError, "Cannot re-cast a dag that's already been set"
        self.dag = other.dag

    def accelerate(self):
        self.dag.accelerate()
        return self.dag


def add_dag(left, right):
    return poly_dag(op = pdADD, left = left, right = right)
def mul_dag(left, right):
    return poly_dag(op = pdMUL, left = left, right = right)
def abc_dag(left, right, index):
    return poly_dag(op = pdABC, left = left, right = right, index = index)
def pow_dag(left, right):
    return poly_dag(op = pdPOW, left = left, right = right)
def sqr_dag(left):
    return poly_dag(op = pdSQR, left = left)
def coeff_dag(index):
    return poly_dag(op = pdCOEFF, index = index)
def var_dag(index):
    return poly_dag(op = pdVAR, index = index)
def univar_dag():
    return poly_dag(op = pdUNIVAR)


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

    cdef void eval(generic_pd self, object vars, object coeffs):
        raise NotImplementedError

    cdef void accelerate(generic_pd self):
        pass

cdef class var_pd(generic_pd):
    def __init__(var_pd self, int index):
        generic_pd.__init__(self)
        self.index = index
    cdef void eval(var_pd self, object vars, object coeffs):
        self.value = vars[self.index]

cdef class univar_pd(generic_pd):
    cdef void eval(univar_pd self, object var, object coeffs):
        self.value = var

cdef class coeff_pd(generic_pd):
    def __init__(coeff_pd self, int index):
        generic_pd.__init__(self)
        self.index = index
    cdef void eval(coeff_pd self, object vars, object coeffs):
        self.value = coeffs[self.index]

cdef class unary_pd(generic_pd):
    def __init__(unary_pd self, poly_dag operand):
        generic_pd.__init__(self)
        self.operand_p = operand
    cdef void accelerate(unary_pd self):
        self.operand = self.operand_p.dag
        self.operand_p = None
        self.operand.refs += 1
        if self.operand.refs == 1:
            self.operand.accelerate()

cdef class sqr_pd(unary_pd):
    cdef void eval(sqr_pd self, object vars, object coeffs):
        pd_eval(self.operand, vars, coeffs)
        self.value = self.operand.value * self.operand.value
        pd_clean(self.operand)

cdef class pow_pd(unary_pd):
    def __init__(unary_pd self, poly_dag base, object exponent):
        unary_pd.__init__(self, base)
        self.exponent = exponent

    cdef void eval(pow_pd self, object vars, object coeffs):
        pd_eval(self.operand, vars, coeffs)
        self.value = self.operand.value ** self.exponent
        pd_clean(self.operand)



cdef class binary_pd(generic_pd):
    def __init__(binary_pd self, poly_dag left, poly_dag right):
        generic_pd.__init__(self)
        self.left_p = left
        self.right_p = right

    cdef void accelerate(binary_pd self):
        self.left = self.left_p.dag
        self.left_p = None
        self.left.refs += 1
        if self.left.refs == 1:
            self.left.accelerate()

        self.right = self.right_p.dag
        self.right_p = None
        self.right.refs += 1
        if self.right.refs == 1:
            self.right.accelerate()




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
    def __init__(abc_pd self, poly_dag left, poly_dag right, int index):
        binary_pd.__init__(self, left, right)
        self.index = index

    cdef void eval(abc_pd self, object vars, object coeffs):
        pd_eval(self.left, vars, coeffs)
        pd_eval(self.right, vars, coeffs)
        self.value = self.left.value * self.right.value + coeffs[self.index]
        pd_clean(self.left)
        pd_clean(self.right)




