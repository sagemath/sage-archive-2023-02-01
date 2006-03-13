"""
Elements of free modules
"""

import sage.misc.misc as misc
import sage.misc.latex as latex
from sage.rings.arith import lcm

from sage.rings.coerce import bin_op, cmp as coerce_cmp

import module_element

import sage.matrix.matrix

def is_FreeModuleElement(x):
    return isinstance(x, FreeModuleElement)

def Vector(R, elts):
    """
    Return a vector over R with given entries.

    INPUT:
        R -- ring
        elts -- entries of a vector
    OUTPUT:
        An element of the free module over R of rank len(elts).

    EXAMPLES:
        sage: v = Vector(Rationals(), [1,1]); v
        (1, 1)
    """
    return (R**len(elts))(elts)

class FreeModuleElement(module_element.ModuleElement):
    """
    An element of a generic free module.
    """
    def __init__(self, parent):
        module_element.ModuleElement.__init__(self, parent)

    def _vector_(self, R):
        return self.change_ring(R)

    def change_ring(self, R):
        P = self.parent()
        if P.base_ring() is R:
            return self
        return P.change_ring(R)(self)

    def __abs__(self):
        return self.norm()

    def __add__(self, right):
        """
        EXAMPLES:
            sage: V = Q^5
            sage: W = V.span([V.1, V.2])
            sage: W.0 + V.0
            (1, 1, 0, 0, 0)
            sage: V.0 + W.0
            (1, 1, 0, 0, 0)
        """
        if not isinstance(right, FreeModuleElement):
            right = self.parent()(right)
        if self.parent() is right.parent():
            V = self.parent()
        elif self.parent().ambient_vector_space() == right.parent().ambient_vector_space():
            V = self.parent().ambient_vector_space()
        else:
            raise TypeError, "self and right must have compatible parents"
        X = [self[i] + right[i] for i in xrange(self.degree())]
        return V(X, coerce_entries=False, copy=False, check_element=False)

    def __cmp__(self, right):
        if not isinstance(right, FreeModuleElement) or \
                  not (self.parent().ambient_vector_space() == right.parent().ambient_vector_space()):
            return coerce_cmp(self, right)
        for i in xrange(self.degree()):
            c = cmp(self[i], right[i])
            if c: return c
        return 0

    def __getitem__(self, i):
        raise NotImplementedError

    def __invert__(self):
        raise NotImplementedError

    def __len__(self):
        return self.parent().degree()

    def __mod__(self, p):
        raise NotImplementedError

    def __mul__(self, right):
        if isinstance(right, FreeModuleElement):
            # vector x vector -> dot product
            return self.dot_product(right)
        if isinstance(right, sage.matrix.matrix.Matrix):
            # vector x matrix multiply
            return self._matrix_multiply(right)

        global TypeError    # it's weird that this is necessary!
        try:
            right = self.base_ring()(right)
        except TypeError:
            try:
                right = self.parent().base_field()(right)
            except AttributeError, TypeError:
                raise TypeError, "unable to multiply %s by %s"%(self, right)
            else:
                return self.parent().ambient_vector_space()(self)._scalar_multiply(right)

        return self._scalar_multiply(right)

    def __div__(self, right):
        return self * ~(self.parent().base_ring()(right))

    def __neg__(self):
        """
        EXAMPLES:
            sage: V = QQ^3
            sage: v = V([1,2,3])
            sage: v.__neg__()
            (-1, -2, -3)
        """
        return self._scalar_multiply(-1)

    def __pos__(self):
        return self

    def __pow__(self, n):
        raise NotImplementedError

    def _repr_(self):
        d = self.degree()
        if d == 0: return "()"
        # compute column widths
        S = [str(x) for x in self.list()]
        #width = max([len(x) for x in S])
        s = "("
        for i in xrange(d):
            if i == d-1:
                sep = ""
            else:
                sep=", "
            entry = S[i]
            #if i > 0:
            #    entry = " "*(width-len(entry)) + entry
            s += entry + sep
        s = s + ")"
        return s

    def __setitem__(self, i, x):
        raise NotImplementedError

    def __sub__(self, right):
        if not isinstance(right, FreeModuleElement) or not self.parent() is right.parent():
            V = self.parent().ambient_vector_space()
            right = V(right)
        else:
            V = self.parent()
        X = [self[i] - right[i] for i in xrange(self.degree())]
        return V(X, coerce_entries=False, copy=False, check_element=False)

    def _matrix_multiply(self, A):
        """
        Return the product self*A.

        EXAMPLES:
            sage: MS = MatrixSpace(QQ,3)
            sage: A = MS([0,1,0,1,0,0,0,0,1])
            sage: V = QQ^3
            sage: v = V([1,2,3])
            sage: v._matrix_multiply(A)
            (2, 1, 3)

        The multiplication operator also just calls \code{_matrix_multiply}:
            sage: v*A
            (2, 1, 3)
        """
        return A.vector_matrix_multiply(self)

    def _scalar_multiply(self, s):
        """
        return the product s*self.
        """
        s = self.base_ring()(s)
        return self.parent()([x*s for x in self.list()], \
                    copy=False, coerce_entries=False, check_element=False)

    def copy(self):
        return self.parent()(self.entries(), \
                             coerce_entries=False, copy=True, check_element=False)
    def degree(self):
        return self.parent().degree()

    def denominator(self):
        R = self.base_ring()
        if self.degree() == 0: return 1
        x = self.list()
        d = x[0].denominator()
        for y in x:
            d = d.lcm(y.denominator())
        return d

    def dict(self):
        e = {}
        for i in xrange(self.degree()):
            c = self[i]
            if c != 0:
                e[i] = c
        return e

    def dot_product(self, right):
        """
        Return the dot product of self and right, which is the sum
        of the product of the corresponding entries.

        INPUT:
            right -- vector of the same degree as self.  it need not
                     be in the same vector space as self, as long as
                     the coefficients can be multiplied.

        EXAMPLES:
            sage: V = FreeModule(ZZ, 3)
            sage: v = V([1,2,3])
            sage: w = V([4,5,6])
            sage: v.dot_product(w)
            32

            sage: W = VectorSpace(GF(3),3)
            sage: w = W([0,1,2])
            sage: w.dot_product(v)
            2
            sage: w.dot_product(v).parent()
            Finite Field of size 3

        Implicit coercion is well defined (irregardless of order), so
        we get 2 even if we do the dot product in the other order.

            sage: v.dot_product(w)
            2
        """
        if not isinstance(right, FreeModuleElement):
            raise TypeError, "right must be a free module element"
        r = right.list()
        l = self.list()
        if len(r) != len(l):
            raise ArithmeticError, "degrees must be the same"%(len(l),len(r))
        zero = self.parent().base_ring()(0)
        return misc.add([l[i]*r[i] for i in xrange(len(l))], zero)

    def element(self):
        return self

    def get(self, i):
        """
        get is meant to be more efficient
        than getitem, because it does not do
        any error checking.
        """
        return self[i]

    def list(self):
        return [self[i] for i in range(self.degree())]

    def additive_order(self):
        return lcm([self[i].order() for i in range(self.degree())])

    def set(self, i, x):
        """
        set is meant to be more efficient
        than setitem, because it does not do
        any error checking or coercion.  Use with care.
        """
        self[i] = x


    def entries(self):
        # This *must* be defined in the derived class.
        # It is the underlying representation of elements.
        raise NotImplementedError

    def inner_product(self, right):
        """
        Returns the inner product of self and other, with respect to
        the inner product defined on the parent of self.

        EXAMPLES: todo
        """
        if self.parent().is_ambient() and self.parent()._inner_product_is_dot_product():
            return self.dot_product(right)
        if not isinstance(right, FreeModuleElement):
            raise TypeError, "right must be a free module element"
        M = self.parent()
        if M.is_ambient() or M.uses_ambient_inner_product():
            A = M.ambient_module().inner_product_matrix()
            return A.linear_combination_of_rows(self).dot_product(right)
        else:
            A = M.inner_product_matrix()
            v = M.coordinate_vector(self)
            w = M.coordinate_vector(right)
            return A.linear_combination_of_rows(v).dot_product(w)

    def is_dense(self):
        return self.parent().is_dense()

    def is_sparse(self):
        return self.parent().is_sparse()

    def is_vector(self):
        return True

    def nonzero_positions(self):
        """
        Return the sorted list of integers i such that self[i] != 0.
        """
        z = self.base_ring()(0)
        v = self.list()
        return [i for i in xrange(self.degree()) if v[i] != z]

    def _latex_(self):
        """
        Return a latex representation of self.  For example, if self is
        the free module element (1,2,3,4), then following latex is
        generated: "(1,2,3,4)"  (without the quotes).
        """
        s = '\\left('
        for a in self.list():
            s += latex.latex(a) + ','
        if len(self.list()) > 0:
            s = s[:-1]  # get rid of last comma
        return s + '\\right)'



#############################################
# Generic dense element
#############################################
class FreeModuleElement_generic_dense(FreeModuleElement):
    def __init__(self, parent, entries, coerce_entries=True, copy=True):
        FreeModuleElement.__init__(self, parent)
        R = self.parent().base_ring()
        if entries == 0:
            zero = R(0)
            entries = [zero for _ in xrange(self.degree())]
        else:
            if not isinstance(entries, list):
                raise TypeError, "entries (=%s) must be a list"%entries
            if len(entries) != self.degree():
                raise ArithmeticError, "entries must be a list of length %s"%\
                            self.degree()
            if coerce_entries:
                try:
                    entries = [R(x) for x in entries]
                except TypeError:
                    raise TypeError, "Unable to coerce entries (=%s) to %s"%(entries, R)
            elif copy:
                # Make a copy
                entries = list(entries)
        self.__entries = entries


    def __getitem__(self, i):
        """
        """
        i = int(i)
        #if not isinstance(i, int):
        #    raise TypeError, "index must an integer"
        if i < 0 or i >= self.degree():
            raise IndexError, "index (i=%s) must be between 0 and %s"%(i,
                            self.degree()-1)
        return self.__entries[i]

    def __setitem__(self, i, value):
        """
        Set entry i of self to value.
        """
        i = int(i)
        #if not isinstance(i, int):
        #    raise TypeError, "index must an integer"
        if i < 0 or i >= self.degree():
            raise IndexError, "index (i=%s) must be between 0 and %s"%(i,
                            self.degree()-1)
        self.__entries[i] = self.base_ring()(value)

    def entries(self):
        return self.__entries

    def list(self):
        return self.__entries


#############################################
# Generic sparse element
#############################################
def _sparse_dot_product(v, w):
    """
    v and w are dictionaries with integer keys.
    """
    x = set(v.keys()).intersection(set(w.keys()))
    return sum([v[k]*w[k] for k in x])

class FreeModuleElement_generic_sparse(FreeModuleElement):
    """
    A generic_sparse is a dictionary with keys ints i and entries in
    the base ring.
    """
    def __init__(self, parent,
                 entries=0,
                 coerce_entries=True,
                 copy=True):
        #WARNING: In creation, we do not check that the i pairs satisfy
        #     0 <= i < degree.
        FreeModuleElement.__init__(self, parent)
        R = self.base_ring()
        if entries == 0:
            entries = {}
        else:
            if isinstance(entries, list):
                if len(entries) != self.degree():
                    raise TypeError, "entries has the wrong length"
                x = entries
                entries = {}
                for i in xrange(self.degree()):
                    if x[i] != 0:
                        entries[i] = x[i]
                copy = False
            if not isinstance(entries, dict):
                raise TypeError, "entries must be a dict"
            if coerce_entries:
                try:
                    for k, x in entries.iteritems():
                        entries[k] = R(x)
                except TypeError:
                    raise TypeError, "Unable to coerce values of entries dict (=%s) to %s"%(entries, R)
            elif copy:
                # Make a copy
                entries = dict(entries)
        self.__entries = entries
        self.__zero = R(0)


    def __getitem__(self, i):
        #if not isinstance(i, int):
        i = int(i)
            #raise TypeError, "index must an integer"
        if i < 0 or i >= self.degree():
            raise IndexError, "index (i=%s) must be between 0 and %s"%(i,
                            self.degree()-1)
        if self.__entries.has_key(i):
            return self.__entries[i]
        return self.__zero

    def get(self, i):
        """
        Like __getitem__ but with no type or bounds checking.
        Returns 0 if access is out of bounds.
        """
        i = int(i)
        if self.__entries.has_key(i):
            return self.__entries[i]
        return self.__zero

    def set(self, i, x):
        """
        Like __setitem__ but with no type or bounds checking.
        """
        i = int(i)
        if x == 0:
            if self.__entries.has_key(i):
                del self.__entries[i]
            return
        self.__entries[i] = x

    def __setitem__(self, i, value):
        """
        """
        i = int(i)
        #if not isinstance(i, int):
        #    raise TypeError, "index must an integer"
        if i < 0 or i >= self.degree():
            raise IndexError, "index (i=%s) must be between 0 and %s"%(i,
                            self.degree()-1)
        self.set(i, value)

    def denominator(self):
        R = self.base_ring()
        x = self.entries()
        if len(x) == 0:
            return 1
        Z = x.iteritems()
        d = Z.next()[1].denominator()
        for _, y in Z:
            d = d.lcm(y.denominator())
        return d

    def entries(self):
        return self.__entries

    def dict(self):
        return self.__entries

    def nonzero_positions(self):
        """
        Returns the set of pairs (i,j) such that self[i,j] != 0.
        """
        K = self.__entries.keys()
        K.sort()
        return K


