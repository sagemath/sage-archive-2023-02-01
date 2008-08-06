#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
from sage.structure.element import ModuleElement
from sage.modules.free_module_element import vector
from sage.misc.misc import repr_lincomb
from sage.modules.module import Module
from sage.rings.all import Ring, Integer
import sage.structure.parent_base
from sage.combinat.family import Family
from sage.combinat.finite_class import FiniteCombinatorialClass
from sage.combinat.combinat import CombinatorialClass

# TODO:
# Rewrite all tests in a more self contained way

class CombinatorialFreeModuleElement(ModuleElement):
    def __init__(self, M, x):
        """
        Create a combinatorial module element x.  This should never
        be called directly, but only through the parent combinatorial
        module's __call__ method.

        """
        ModuleElement.__init__(self, M)
        self._monomial_coefficients = x

    def __iter__(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + s([3])
            sage: [i for i in sorted(a)]
            [([2, 1], 1), ([3], 1)]

        """
        return self._monomial_coefficients.iteritems()

    def __contains__(self, x):
        """
        Returns whether or not a combinatorial object x indexing a basis
        element is in the support of self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1]) + s([3])
            sage: Partition([2,1]) in a
            True
            sage: Partition([1,1,1]) in a
            False
        """
        return x in self._monomial_coefficients and self._monomial_coefficients[x] != 0

    def monomial_coefficients(self):
        """
        Return the internal dictionary which has the combinatorial
        objects indexing the basis as keys and their corresponding
        coefficients as values.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+2*s([3,2])
            sage: d = a.monomial_coefficients()
            sage: type(d)
            <type 'dict'>
            sage: d[ Partition([2,1]) ]
            1
            sage: d[ Partition([3,2]) ]
            2
        """
        return self._monomial_coefficients

    def __repr__(self):
        """
        EXAMPLES:
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='F')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b']
            F('a') + 2*F('b')

        """
        v = self._monomial_coefficients.items()
        v.sort()
        prefix = self.parent().prefix()
        mons = [ prefix + "(" + repr(m) + ")" for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def _latex_(self):
        """
        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: a = 2 + QS3([2,1,3])
            sage: latex(a) #indirect doctest
            2[1,2,3] + [2,1,3]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        prefix = self.parent().prefix()
        if prefix == "":
            mons = [ prefix + '[' + ",".join(map(str, m)) + ']' for (m, _) in v ]
        else:
            mons = [ prefix + '_{' + ",".join(map(str, m)) + '}' for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs, is_latex=True).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def __cmp__(left, right):
        """
        The ordering is the one on the underlying sorted list of
        (monomial,coefficients) pairs.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: b = s([1,1,1])
            sage: cmp(a,b) #indirect doctest
            1
        """
        nonzero = lambda mc: mc[1] != 0
        v = filter(nonzero, left._monomial_coefficients.items())
        v.sort()
        w = filter(nonzero, right._monomial_coefficients.items())
        w.sort()
        return cmp(v, w)

    def _add_(self, y):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]) + s([5,4]) # indirect doctest
            s[2, 1] + s[5, 4]
            sage: a = s([2,1]) + 0
            sage: len(a.monomial_coefficients())
            1
        """
        A = self.parent()
        BR = A.base_ring()
        z_elt = dict(self._monomial_coefficients)
        for m, c in y._monomial_coefficients.iteritems():
            if z_elt.has_key(m):
                cm = z_elt[m] + c
                if cm == 0:
                    del z_elt[m]
                else:
                    z_elt[m] = cm
            else:
                z_elt[m] = c


        #Remove all entries that are equal to 0
        del_list = []
        zero = BR(0)
        for m, c in z_elt.iteritems():
            if c == zero:
                del_list.append(m)
        for m in del_list:
            del z_elt[m]

        return A._from_dict(z_elt)


    def _neg_(self):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: -s([2,1]) # indirect doctest
            -s[2, 1]
        """
        return self.map_coefficients(lambda c: -c)


    def _sub_(self, y):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]) - s([5,4]) # indirect doctest
            s[2, 1] - s[5, 4]
        """
        A = self.parent()
        BR = A.base_ring()
        z_elt = dict(self._monomial_coefficients)
        for m, c in y._monomial_coefficients.iteritems():
            if z_elt.has_key(m):
                cm = z_elt[m] - c
                if cm == 0:
                    del z_elt[m]
                else:
                    z_elt[m] = cm
            else:
                z_elt[m] = -c

        #Remove all entries that are equal to 0
        zero = BR(0)
        del_list = []
        for m, c in z_elt.iteritems():
            if c == zero:
                del_list.append(m)
        for m in del_list:
            del z_elt[m]

        return A._from_dict(z_elt)


    def _coefficient_fast(self, m, default=None):
        """
        Returns the coefficient of m in self, where m is key
        in self._monomial_coefficients.

        EXAMPLES:
            sage: p = Partition([2,1])
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])
            sage: a._coefficient_fast([2,1])
            Traceback (most recent call last):
            ...
            TypeError: list objects are unhashable
            sage: a._coefficient_fast(p)
            1
            sage: a._coefficient_fast(p, 2)
            1
        """
        if default is None:
            default = self.base_ring()(0)
        return self._monomial_coefficients.get(m, default)

    def coefficient(self, m):
        """

        # NT: coefficient_fast should be the default, just with appropriate assertions
        # that can be turned on or of

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) - 2*s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficient([4])
            1
            sage: z.coefficient([2,1])
            -2
        """
        p = self.parent()
        if isinstance(m, p._combinatorial_class.object_class):
            return self._monomial_coefficients.get(m, p.base_ring().zero_element())
        if m in p._combinatorial_class:
            return self._monomial_coefficients.get(p._combinatorial_class.object_class(m), p.base_ring().zero_element())
        else:
            raise TypeError, "you must specify an element of %s"%p._combinatorial_class


    def is_zero(self):
        """
        Returns True if and only self == 0.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s([2,1]).is_zero()
            False
            sage: s(0).is_zero()
            True
            sage: (s([2,1]) - s([2,1])).is_zero()
            True
        """
        BR = self.parent().base_ring()
        for v in self._monomial_coefficients.values():
            if v != BR(0):
                return False
        return True

    def __len__(self):
        """
        Returns the number of basis elements of self with
        nonzero coefficients.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: len(z)
            4
        """
        return self.length()

    def length(self):
        """
        Returns the number of basis elements of self with
        nonzero coefficients.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.length()
            4
        """
        return len( filter(lambda x: self._monomial_coefficients[x] != 0, self._monomial_coefficients) )


    def support(self):
        """
        # NT: ???
        Returns a pair [mons, cffs] of lists of the monomials
        of self (mons) and their respective coefficients (cffs).

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.support()
            [[[1], [1, 1, 1], [2, 1], [4]], [1, 1, 1, 1]]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        mons = [ m for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        return [mons, cffs]

    def monomials(self):
        """
        Returns a list of the combinatorial objects indexing
        the basis elements of self which non-zero coefficients.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.monomials()
            [[1], [1, 1, 1], [2, 1], [4]]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        mons = [ m for (m, _) in v ]
        return mons

    def coefficients(self):
        """
        Returns a list of the coefficents appearing on the
        basiselements in self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficients()
            [1, 1, 1, 1]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        cffs = [ c for (_, c) in v ]
        return cffs

    def _vector_(self, new_BR=None):
        """
        Returns a vector version of self. If new_BR is specified,
        then in returns a vector over new_BR.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a._vector_()
            (2, 0, 0, 0, 0, 4)
            sage: vector(a)
            (2, 0, 0, 0, 0, 4)
            sage: a._vector_(RDF)
            (2.0, 0.0, 0.0, 0.0, 0.0, 4.0)
        """
        parent = self.parent()
        if parent.get_order() is None:
            cc = parent._combinatorial_class
        else:
            cc = parent.get_order()

        if new_BR is None:
            new_BR = parent.base_ring()

        return vector(new_BR, [new_BR(self._monomial_coefficients.get(m, 0)) for m in cc])

    def to_vector(self):
        """
        Returns a vector version of self.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a.to_vector()
            (2, 0, 0, 0, 0, 4)
        """
        return self._vector_()


    def map_coefficients(self, f):
        """
        Returns a new element of self.parent() obtained
        by applying the function f to all of the coefficients
        of self.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+2*s([3,2])
            sage: a.map_coefficients(lambda x: x*2)
            2*s[2, 1] + 4*s[3, 2]
        """
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            z_elt[m] = f(c)
        return self.parent()._from_dict(z_elt)


    def map_basis(self, f):
        """
        Returns a new element of self.parent() obtained
        by applying the function f to all of the combinatorial
        objects indexing the basis elements.

        FIXME: rename to map_support?

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+2*s([3,2])
            sage: a.map_basis(lambda x: x.conjugate())
            s[2, 1] + 2*s[2, 2, 1]
        """
        res = self.parent()(0)
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            z_elt[f(m)] = c
        res._monomial_coefficients = z_elt
        return res

    def map(self, f):
        """
        Returns a new element of self.parent() obtained
        by applying the function f to a monomial coefficient
        (m,c) pair.  f returns a (new_m, new_c) pair.

        FIXME: map_monomial?

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: f = lambda m,c: (m.conjugate(), 2*c)
            sage: a = s([2,1]) + s([1,1,1])
            sage: a.map_mc(f)
            2*s[2, 1] + 2*s[3]
        """
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            new_m, new_c = f(m,c)
            z_elt[new_m] = new_c
        return self.parent()._from_dict(z_elt)

    map_mc = map

    def _l_action(self, x):
        x = self.base_ring()(x)
        return self.map_coefficients(lambda c: x*c)

    def _r_action(self, x):
        x = self.base_ring()(x)
        return self.map_coefficients(lambda c: c*x)


# NT: There is nothing combinatorial here
# NT: It's really too bad that we can't have the ring as second optional argument
class CombinatorialFreeModuleInterface(): # Should not it inherit from ParentWithBase?
    def __init__(self, R, element_class):
        #Make sure R is a ring with unit element
        if not isinstance(R, Ring):
            raise TypeError, "Argument R must be a ring."
        try:
            # R._one_element?
            z = R(Integer(1))
        except:
            raise ValueError, "R must have a unit element"

        self._element_class = element_class
        self._order = None

        #Initialize the base structure
        sage.structure.parent_base.ParentWithBase.__init__(self, R)

    _prefix = ""
    _name   = "CombinatorialModule -- change me"


    # Should be an attribute?
    def basis(self):
        """
        Returns a list of the basis elements of self.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: list(QS3.basis())
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F.basis()
            Finite family {'a': B('a'), 'c': B('c'), 'b': B('b')}

        """
        return Family(self._combinatorial_class, self.term)

    def __call__(self, x):
        """
        Coerce x into self.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3(2)
            2*[1, 2, 3]
            sage: QS3([2,3,1])
            [2, 3, 1]
        """
        R = self.base_ring()
        eclass = self._element_class

        #Coerce ints to Integers
        if isinstance(x, int):
            x = Integer(x)


        if hasattr(self, '_coerce_start'):
            try:
                return self._coerce_start(x)
            except TypeError:
                pass

        # NT: coercion from the indexing set is louzy:
        # in FreeModule(ZZ, [1,2,3]), how do tell the difference between
        # an element of the ground ring and one of the indexing set?

        #x is an element of the same type of combinatorial algebra
        if hasattr(x, 'parent') and x.parent().__class__ is self.__class__:
            P = x.parent()
            #same base ring
            if P is self:
                return x
            #different base ring -- coerce the coefficients from into R
            else:
                return eclass(self, dict([ (e1,R(e2)) for e1,e2 in x._monomial_coefficients.items()]))
        #x is an element of the basis combinatorial class
        elif isinstance(self._combinatorial_class.object_class, type) and isinstance(x, self._combinatorial_class.object_class):
            return eclass(self, {x:R(1)})
        elif x in self._combinatorial_class:
            return eclass(self, {self._combinatorial_class.object_class(x):R(1)})
        else:
            if hasattr(self, '_coerce_end'):
                try:
                    return self._coerce_end(x)
                except TypeError:
                    pass
            raise TypeError, "do not know how to make x (= %s) an element of self (=%s)"%(x,self)

    def _an_element_impl(self):
        """
        Returns an element of self, namely the unit element.

        EXAMPLES:
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._an_element_impl()
            0
            sage: _.parent() is F
            True
        """
        return self._element_class(self, {})

    def __repr__(self):
        """
        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: print QS3.__repr__()
            Symmetric group algebra of order 3 over Rational Field
        """
        return self._name + " over %s"%self.base_ring()

    def combinatorial_class(self):
        """
        Returns the combinatorial class that indexes the basis
        elements.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s.combinatorial_class()
            Partitions
        """
        return self._combinatorial_class

    def _coerce_impl(self, x):
        """
        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s._coerce_impl(2)
            2*s[]
        """
        try:
            R = x.parent()
            if R.__class__ is self.__class__:
                #Only perform the coercion if we can go from the base
                #ring of x to the base ring of self
                if self.base_ring().has_coerce_map_from( R.base_ring() ):
                    return self(x)
        except AttributeError:
            pass

        # any ring that coerces to the base ring
        return self._coerce_try(x, [self.base_ring()])

    def dimension(self):
        """
        Returns the dimension of the combinatorial algebra (which is given
        by the number of elements in the associated combinatorial class).

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: s.dimension()
            +Infinity
        """
        return self._combinatorial_class.count()

    def set_order(self, order):
        """
        Sets the order of the elements of the combinatorial class.

        If .set_order() has not been called, then the ordering is
        the one used in the generation of the elements of self's
        associated combinatorial class.

        EXAMPLES:
            sage: QS2 = SymmetricGroupAlgebra(QQ,2)
            sage: b = list(QS2.basis().keys())
            sage: b.reverse()
            sage: QS2.set_order(b)
            sage: QS2.get_order()
            [[2, 1], [1, 2]]

        """
        self._order = order

    def get_order(self):
        """
        Returns the order of the elements in the basis.

        EXAMPLES:
            sage: QS2 = SymmetricGroupAlgebra(QQ,2)
            sage: QS2.get_order()
            [[1, 2], [2, 1]]
        """
        if self._order is None:
            self._order = self.combinatorial_class().list()
        return self._order

    def prefix(self):
        """
        Returns the prefix used when displaying elements of self.

        EXAMPLES:
            sage: X = SchubertPolynomialRing(QQ)
            sage: X.prefix()
            'X'
        """
        return self._prefix

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: XQ = SchubertPolynomialRing(QQ)
            sage: XZ = SchubertPolynomialRing(ZZ)
            sage: XQ == XZ #indirect doctest
            False
            sage: XQ == XQ
            True
        """
        if not isinstance(other, self.__class__):
            return -1
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        return 0

    def _apply_module_morphism(self, x, f):
        """
        Returns the image of x under the module morphism defined by
        extending f by linearity.

        INPUT:
            -- x : a element of self
            -- f : a function that takes in a combinatorial object
                   indexing a basis element and returns an element
                   of the target domain

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: a = s([3]) + s([2,1]) + s([1,1,1])
            sage: b = 2*a
            sage: f = lambda part: len(part)
            sage: s._apply_module_morphism(a, f) #1+2+3
            6
            sage: s._apply_module_morphism(b, f) #2*(1+2+3)
            12


        """
        res = 0
        for m, c in x._monomial_coefficients.iteritems():
            res += c*f(m)
        return res


    def _apply_module_endomorphism(self, a, f):
        """
        This takes in a function from the basis elements
        to the elements of self and applies it linearly
        to a. Note that _apply_module_endomorphism does not
        require multiplication on self to be defined.

        EXAMPLES:
            sage: s = SFASchur(QQ)
            sage: f = lambda part: 2*s(part.conjugate())
            sage: s._apply_module_endomorphism( s([2,1]) + s([1,1,1]), f)
            2*s[2, 1] + 2*s[3]

        """
        mcs = a.monomial_coefficients()
        base_ring = self.base_ring()
        zero = base_ring(0)

        z_elt = {}
        for basis_element in mcs:
            f_mcs = f(basis_element).monomial_coefficients()
            for f_basis_element in f_mcs:
                z_elt[ f_basis_element ] = z_elt.get(f_basis_element, zero) + mcs[basis_element]*f_mcs[f_basis_element]

        return self._from_dict(z_elt)

    def term(self, i):
        """
        EXAMPLES:
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.term('a')
            B('a')
        """
        return self._from_dict({i:self.base_ring().one_element()})

    def _from_dict(self, d, coerce=False):
        """
        Given a monomial coefficient dictionary d, return the element
        of self with the dictionary.

        EXAMPLES:
            sage: e = SFAElementary(QQ)
            sage: s = SFASchur(QQ)
            sage: a = e([2,1]) + e([1,1,1]); a
            e[1, 1, 1] + e[2, 1]
            sage: s._from_dict(a.monomial_coefficients())
            s[1, 1, 1] + s[2, 1]

            sage: part = Partition([2,1])
            sage: d = {part:1}
            sage: a = s._from_dict(d,coerce=True); a
            s[2, 1]
            sage: a.coefficient(part).parent()
            Rational Field
        """
        if coerce:
            R = self.base_ring()
            d = [ (m,R(c)) for m,c in d.iteritems() ]
            d = dict(d)

        return self._element_class(self, d)


class CombinatorialFreeModule(CombinatorialFreeModuleInterface, Module):
    r"""
    EXAMPLES:
        We construct a free module whose basis is indexed by the letters a,b,c:

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F
            Free module generated by ['a', 'b', 'c'] over Rational Field

        Its basis is a family, indexed by a,b,c:
        FIXME in family: we should preserve the order of the indices
	    sage: e = F.basis()

            sage: e.keys()
            ['a', 'b', 'c']
            sage: list(e)
            [B('a'), B('b'), B('c')]

        Let us construct some elements, and compute with them:
            sage: e['a']
            B('a')
            sage: 2*e['a']
            2*B('a')
            sage: e['a'] + 3*e['b']
            B('a') + 3*B('b')
    """


    def __init__(self, R, cc, prefix="B"):
        if isinstance(cc, list):
            cc = FiniteCombinatorialClass(cc)
        if not isinstance(cc, CombinatorialClass):
            raise TypeError, "cc = (%s) must be an instance of CombinatorialClass"%cc
        self._combinatorial_class = cc
        self._prefix = prefix
        self._name = "Free module generated by %s"%cc
        CombinatorialFreeModuleInterface.__init__(self, R, CombinatorialFreeModuleElement)
    pass
