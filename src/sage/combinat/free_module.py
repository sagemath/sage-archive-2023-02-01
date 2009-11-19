"""
Free modules
"""
#*****************************************************************************
#       Copyright (C) 2007      Mike Hansen <mhansen@gmail.com>,
#                     2007-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.modules.free_module_element import vector
from sage.misc.misc import repr_lincomb
from sage.modules.module import Module
from sage.rings.all import Integer
import sage.structure.parent_base
import sage.structure.element
from sage.combinat.family import Family
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.combinat.cartesian_product import CartesianProduct
from sage.categories.all import ModulesWithBasis
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.misc.all import lazy_attribute
from sage.categories.poor_man_map import PoorManMap

# TODO: move the content of this class to CombinatorialFreeModule.Element and ModulesWithBasis.Element
class CombinatorialFreeModuleElement(Element):
    def __init__(self, M, x):
        """
        Create a combinatorial module element. This should never be
        called directly, but only through the parent combinatorial
        free module's :meth:`__call__` method.

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']; f
            B['a'] + 3*B['c']
            sage: f == loads(dumps(f))
            True
        """
        Element.__init__(self, M)
        self._monomial_coefficients = x

    def __iter__(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: [i for i in sorted(f)]
            [('a', 1), ('c', 3)]

        ::

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

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: 'a' in f
            True
            sage: 'b' in f
            False

        ::

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
        Return the internal dictionary which has the combinatorial objects
        indexing the basis as keys and their corresponding coefficients as
        values.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: d = f.monomial_coefficients()
            sage: d['a']
            1
            sage: d['c']
            3

        ::

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

    def _repr_(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='F')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            F['a'] + 2*F['b']
        """
        v = self._monomial_coefficients.items()
        try:
            v = sorted(v)
        except StandardError: # Sorting the output is a plus, but if we can't, no big deal
            pass
        repr_term = self.parent()._repr_term
        mons = [ repr_term(m) for (m, _) in v ]
        cffs = [ x for (_, x) in v ]
        x = repr_lincomb(mons, cffs).replace("*1 "," ")
        if x[len(x)-2:] == "*1":
            return x[:len(x)-2]
        else:
            return x

    def _latex_(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: latex(f)
            B_{a} + 3B_{c}

        ::

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

        EXAMPLES::

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
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a'] + 3*B['c']
            B['a'] + 3*B['c']

        ::

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
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: -f
            -B['a'] - 3*B['c']

        ::

            sage: s = SFASchur(QQ)
            sage: -s([2,1]) # indirect doctest
            -s[2, 1]
        """
        return self.map_coefficients(lambda c: -c)


    def _sub_(self, y):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a'] - 3*B['c']
            B['a'] - 3*B['c']

        ::

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
        Returns the coefficient of m in self, where m is key in
        self._monomial_coefficients.

        EXAMPLES::

            sage: p = Partition([2,1])
            sage: q = Partition([1,1,1])
            sage: s = SFASchur(QQ)
            sage: a = s(p)
            sage: a._coefficient_fast([2,1])
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'

        ::

            sage: a._coefficient_fast(p)
            1
            sage: a._coefficient_fast(p, 2)
            1
            sage: a._coefficient_fast(q)
            0
            sage: a._coefficient_fast(q, 2)
            2
            sage: a[p]
            1
            sage: a[q]
            0
        """
        if default is None:
            default = self.base_ring()(0)
        return self._monomial_coefficients.get(m, default)

    __getitem__ = _coefficient_fast

    def coefficient(self, m):
        """
        EXAMPLES::

            sage: s = CombinatorialFreeModule(QQ, Partitions())
            sage: z = s([4]) - 2*s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficient([4])
            1
            sage: z.coefficient([2,1])
            -2
            sage: z.coefficient(Partition([2,1]))
            -2
            sage: z.coefficient([1,2])
            Traceback (most recent call last):
            ...
            AssertionError: [1, 2] should be an element of Partitions
            sage: z.coefficient(Composition([2,1]))
            Traceback (most recent call last):
            ...
            AssertionError: [2, 1] should be an element of Partitions

        Test that coefficient also works for those parents that do not yet have an element_class::

            sage: G = DihedralGroup(3)
            sage: F = CombinatorialFreeModule(QQ, G)
            sage: hasattr(G, "element_class")
            False
            sage: g = G.an_element()
            sage: (2*F.term(g)).coefficient(g)
            2
        """
        # NT: coefficient_fast should be the default, just with appropriate assertions
        # that can be turned on or off
        C = self.parent()._basis_keys
        assert m in C, "%s should be an element of %s"%(m, C)
        if hasattr(C, "element_class") and not isinstance(m, C.element_class):
            m = C(m)
        return self._coefficient_fast(m)


    def is_zero(self):
        """
        Returns True if and only self == 0.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.is_zero()
            False
            sage: F.zero().is_zero()
            True

        ::

            sage: s = SFASchur(QQ)
            sage: s([2,1]).is_zero()
            False
            sage: s(0).is_zero()
            True
            sage: (s([2,1]) - s([2,1])).is_zero()
            True
        """
        BR = self.parent().base_ring()
        return all(v == BR(0) for v in self._monomial_coefficients.values())

    def __len__(self):
        """
        Returns the number of basis elements of self with nonzero
        coefficients.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: len(f)
            2

        ::

            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: len(z)
            4
        """
        return self.length()

    def length(self):
        """
        Returns the number of basis elements of self with nonzero
        coefficients.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.length()
            2

        ::

            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.length()
            4
        """
        return len([mon for mon,coeff in self._monomial_coefficients.items() if coeff !=0 ])

    def support(self):
        """
        Returns a list of the combinatorial objects indexing the basis
        elements of self which non-zero coefficients.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.support()
            ['a', 'c']

        ::

            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.support()
            [[1], [1, 1, 1], [2, 1], [4]]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        mons = [ m for (m, _) in v ]
        # Usually we test the coeff for being non zero, why not here?
        # The same question arises in the following functions.
        return mons

    def monomials(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c']
            sage: f.monomials()
            [B['a'], B['c']]
        """
        P = self.parent()
        one = P.base_ring()(1)
        v = self._monomial_coefficients.items()
        v.sort()
        return [P._from_dict({key:one}) for key,value in v]

    def terms(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c']
            sage: f.terms()
            [B['a'], 2*B['c']]
        """
        P = self.parent()
        v = self._monomial_coefficients.items()
        v.sort()
        return [P._from_dict({key:value}) for key,value in v]

    def coefficients(self):
        """
        Returns a list of the coefficients appearing on the basiselements in
        self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.coefficients()
            [1, -3]

        ::

            sage: s = SFASchur(QQ)
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficients()
            [1, 1, 1, 1]
        """
        v = self._monomial_coefficients.items()
        v.sort()
        cffs = [ c for (_, c) in v ]
        return cffs

    def _vector_(self, new_base_ring=None):
        """
        Returns a vector version of self. If ''new_base_ring'' is specified,
        then in returns a vector over ''new_base_ring''.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: vector(f)
            (1, 0, -3)

        ::

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
        cc = parent.get_order()

        if new_base_ring is None:
            new_base_ring = parent.base_ring()
        # FIXME: should return a sparse vector
        return vector(new_base_ring,
                      [new_base_ring(self._monomial_coefficients.get(m, 0))
                       for m in list(cc)])

    def to_vector(self):
        """
        Returns a vector version of self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.to_vector()
            (1, 0, -3)

        ::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a.to_vector()
            (2, 0, 0, 0, 0, 4)
            sage: a == QS3.from_vector(a.to_vector())
            True
        """
        return self._vector_()


    def map_coefficients(self, f):
        """
        Returns a new element of self.parent() obtained by applying the
        function f to all of the coefficients of self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.map_coefficients(lambda x: x+5)
            6*B['a'] + 2*B['c']

        ::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+2*s([3,2])
            sage: a.map_coefficients(lambda x: x*2)
            2*s[2, 1] + 4*s[3, 2]
        """
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            z_elt[m] = f(c)
        return self.parent()._from_dict(z_elt)


    def map_support(self, f):
        """
        Returns a new element of self.parent() obtained by applying the
        function f to all of the combinatorial objects indexing the basis
        elements.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: a = s([2,1])+2*s([3,2])
            sage: a.map_support(lambda x: x.conjugate())
            s[2, 1] + 2*s[2, 2, 1]
        """
        res = self.parent()(0)
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            z_elt[f(m)] = c
        res._monomial_coefficients = z_elt
        return res

    def map_monomial(self, f):
        """
        Returns a new element of self.parent() obtained by applying the
        function f to a monomial coefficient (m,c) pair. f returns a
        (new_m, new_c) pair.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: f = lambda m,c: (m.conjugate(), 2*c)
            sage: a = s([2,1]) + s([1,1,1])
            sage: a.map_monomial(f)
            2*s[2, 1] + 2*s[3]
        """
        z_elt = {}
        for m,c in self.monomial_coefficients().iteritems():
            new_m, new_c = f(m,c)
            z_elt[new_m] = new_c
        return self.parent()._from_dict(z_elt)

    map_mc = map_monomial

    def _acted_upon_(self, scalar, self_on_left = False):
        """
        Returns the action of a scalar on self

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a']*(1/2)  # indirect doctest
            1/2*B['a']
            sage: B['a']/2
            1/2*B['a']
            sage: B['a']*2      # indirect doctest
            2*B['a']
            sage: B['a']*int(2) # indirect doctest
            2*B['a']

            sage: 1/2*B['a']
            1/2*B['a']
            sage: 2*B['a']      # indirect doctest
            2*B['a']
            sage: int(2)*B['a'] # indirect doctest
            2*B['a']

        TESTS::

            sage: F.get_action(QQ, operator.mul, True)
            Right action by Rational Field on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: F.get_action(QQ, operator.mul, False)
            Left action by Rational Field on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: F.get_action(ZZ, operator.mul, True)
            Right Integer Multiplication by Integer Ring on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: F.get_action(F, operator.mul, True)
            sage: F.get_action(F, operator.mul, False)

        TODO:
         - add non commutative tests
        """
        # With the current design, the coercion model does not have
        # enough information to detect apriori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        if scalar.parent() != self.base_ring():
            return None
        if self_on_left:
            return self.map_coefficients(lambda c: c*scalar)
        else:
            return self.map_coefficients(lambda c: scalar*c)

    # For backward compatibility
    _lmul_ = _acted_upon_
    _rmul_ = _acted_upon_

    def __div__(self, x):
        """
        Division by coefficients

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, [1,2,3])
            sage: x = F._from_dict({1:2, 2:3})
            sage: x/2
            B[1] + 3/2*B[2]

        ::

            sage: F = CombinatorialFreeModule(QQ, [1,2,3])
            sage: B = F.basis()
            sage: f = 2*B[2] + 4*B[3]
            sage: f/2
            B[2] + 2*B[3]
        """
        if self.base_ring().is_field():
            x = self.base_ring()(x)
            return self.map_coefficients(lambda c: c/x)
        else:
            return self.map_coefficients(lambda c: _divide_if_possible(c, x))

def _divide_if_possible(x, y):
    """
    EXAMPLES::

        sage: from sage.combinat.free_module import _divide_if_possible
        sage: _divide_if_possible(4, 2)
        2
        sage: _.parent()
        Integer Ring

    ::

        sage: _divide_if_possible(4, 3)
        Traceback (most recent call last):
        ...
        ValueError: 4 is not divisible by 3
    """
    q, r = x.quo_rem(y)
    if r != 0:
        raise ValueError, "%s is not divisible by %s"%(x, y)
    else:
        return q

class CombinatorialFreeModuleInterface(Parent):
    Element = CombinatorialFreeModuleElement

    def __init__(self, R, element_class = None, category = None):
        r"""
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'], category = FiniteDimensionalModulesWithBasis(QQ))
            sage: F.basis()
            Finite family {'a': B['a'], 'c': B['c'], 'b': B['b']}
            sage: F.category()
            Category of finite dimensional modules with basis over Rational Field
            sage: TestSuite(F).run()
        """
        #Make sure R is a ring with unit element
        from sage.categories.all import Rings
        if R not in Rings():
            raise TypeError, "Argument R must be a ring."
        try:
            # R._one_element?
            z = R(Integer(1))
        except:
            raise ValueError, "R must have a unit element"

        if category is None:
            category = ModulesWithBasis(R)
        if element_class is not None:
            self.Element = element_class
        Parent.__init__(self, base = R, category = category,
                        # Could we get rid of this?
                        element_constructor = self._element_constructor_)

        if not hasattr(self, "_name"):
            self._name = "CombinatorialModule -- change me"

        self._order = None

    #_prefix = ""

    # mostly for backward compatibility
    @lazy_attribute
    def _element_class(self):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F._element_class
            <class 'sage.combinat.free_module.CombinatorialFreeModule_with_category.element_class'>
        """
        return self.element_class

    def basis(self):
        """
        Returns the basis of self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F.basis()
            Finite family {'a': B['a'], 'c': B['c'], 'b': B['b']}

        ::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: list(QS3.basis())
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        return Family(self._basis_keys, self.term)

    def _an_element_(self):
        """
        EXAMPLES::

            sage: CombinatorialFreeModule(QQ, ("a", "b", "c")).an_element()
            2*B['a'] + 2*B['b'] + 3*B['c']
            sage: CombinatorialFreeModule(QQ, ("a", "b", "c"))._an_element_()
            2*B['a'] + 2*B['b'] + 3*B['c']
            sage: CombinatorialFreeModule(QQ, ()).an_element()
            0
            sage: CombinatorialFreeModule(QQ, ZZ).an_element()
            3*B[-1] + B[0] + 3*B[1]
            sage: CombinatorialFreeModule(QQ, RR).an_element()
            B[1.00000000000000]
        """
        # Try a couple heuristics to build a not completely trivial
        # element, while handling cases where R.an_element is not
        # implemented, or R has no iterator, or R has few elements.
        x = self.zero()
        I = self.basis().keys()
        R = self.base_ring()
        try:
            x = x + self.term(I.an_element())
        except:
            pass
        try:
            g = iter(self.basis().keys())
            for c in range(1,4):
                x = x + self.monomial(g.next(), R(c))
        except:
            pass
        return x

    # What semantic do we want for containment?
    # Accepting anything that can be coerced is not reasonnable, especially
    # if we allow coercion from the enumerated set.
    # Accepting anything that can be converted is an option, but that would
    # be expensive. So far, x in self if x.parent() == self

    def __contains__(self, x):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ,["a", "b"])
            sage: G = CombinatorialFreeModule(ZZ,["a", "b"])
            sage: F.term("a")in F
            True
            sage: G.term("a")in F
            False
            sage: "a" in F
            False
            sage: 5/3 in F
            False
        """
        return sage.structure.element.parent(x) == self # is self?

    def _element_constructor_(self, x):
        """
        Coerce x into self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ,["a", "b"])
            sage: F(F.term("a"))   # indirect doctest
            B['a']

        Do not rely on the following feature which may be removed in the future::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3([2,3,1])     # indirect doctest
            [2, 3, 1]

        instead, use::

            sage: P = QS3.basis().keys()
            sage: QS3.term(P([2,3,1]))   # indirect doctest
            [2, 3, 1]

        or:
            sage: B = QS3.basis()
            sage: B[P([2,3,1])]
            [2, 3, 1]

        TODO: The symmetric group algebra (and in general,
        combinatorial free modules on word-like object could instead
        provide an appropriate short-hand syntax QS3[2,3,1]).

        Rationale: this conversion is ambiguous in situations like::

            sage: F = CombinatorialFreeModule(QQ,[0,1])

        Is ``0`` the zero of the base ring, or the index of a basis
        element?  I.e. should the result be ``0`` or ``B[0]``?

            sage: F = CombinatorialFreeModule(QQ,[0,1])
            sage: F(0) # this feature may eventually disappear
            0
            sage: F(1)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= 1) an element of Free module generated by ... over Rational Field

        It is preferable not to rely either on the above, and instead, use::

            sage: F.zero()
            0

        Note that, on the other hand, conversions from the ground ring
        are systematically defined (and mathematically meaningful) for
        algebras.

        Conversions between distinct free modules are not allowed any
        more::

            sage: F = CombinatorialFreeModule(ZZ, ["a", "b"]);      F.rename("F")
            sage: G = CombinatorialFreeModule(QQ, ["a", "b"]);      G.rename("G")
            sage: H = CombinatorialFreeModule(ZZ, ["a", "b", "c"]); H.rename("H")
            sage: G(F.term("a"))
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= B['a']) an element of self (=G)
            sage: H(F.term("a"))
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= B['a']) an element of self (=H)

        Here is a real life example illustrating that this yielded
        mathematically wrong results::

            sage: S = SymmetricFunctions(QQ)
            sage: s = S.s(); p = S.p()
            sage: ss = tensor([s,s]); pp = tensor([p,p])
            sage: a = tensor((s[5],s[5]))
            sage: pp(a) # used to yield p[[5]] # p[[5]]
            Traceback (most recent call last):
               ...
            NotImplementedError

        Extensions of the ground ring should probably be reintroduced
        at some point, but via coercions, and with stronger sanity
        checks (ensuring that the codomain is really obtained by
        extending the scalar of the domain; checking that they share
        the same class is not sufficient).

        TESTS:

        Conversion from the ground ring is implemented for algebras::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3(2)
            2*[1, 2, 3]
        """
        R = self.base_ring()
        eclass = self.element_class

        #Coerce ints to Integers
        if isinstance(x, int):
            x = Integer(x)

#         if hasattr(self, '_coerce_start'):
#             try:
#                 return self._coerce_start(x)
#             except TypeError:
#                 pass

        # x is an element of the same type of combinatorial free module
        # (disabled: this yields mathematically wrong results)
        #if hasattr(x, 'parent') and x.parent().__class__ is self.__class__:
        #    P = x.parent()
        #    #same base ring
        #    if P is self:
        #        return x
        #    #different base ring -- coerce the coefficients from into R
        #    else:
        #        return eclass(self, dict([ (e1,R(e2)) for e1,e2 in x._monomial_coefficients.items()]))
        #x is an element of the ground ring (will be disabled at some point: not so useful)
        if x in R:
            if x == 0:
                return self.zero()
            else:
                raise TypeError, "do not know how to make x (= %s) an element of %s"%(x, self)
        #x is an element of the basis enumerated set;
        # This is a very ugly way of testing this
        elif ((hasattr(self._basis_keys, 'element_class') and
               isinstance(self._basis_keys.element_class, type) and
               isinstance(x, self._basis_keys.element_class))
              or (sage.structure.element.parent(x) == self._basis_keys)):
            return self.term(x)
        elif x in self._basis_keys:
            return self.term(self._basis_keys(x))
        else:
            if hasattr(self, '_coerce_end'):
                try:
                    return self._coerce_end(x)
                except TypeError:
                    pass
            raise TypeError, "do not know how to make x (= %s) an element of self (=%s)"%(x,self)

    def _an_element_impl(self):
        """
        Returns an element of self, namely the zero element.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._an_element_impl()
            0
            sage: _.parent() is F
            True
        """
        return self.element_class(self, {})

    def _repr_(self):
        """
        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3 # indirect doctest
            Symmetric group algebra of order 3 over Rational Field
            sage: QS3.__custom_name = "QS3"
            sage: QS3 # indirect doctest
            QS3

        TODO: allow for complete overriding by setting __custom_name or something similar
        """
        return self._name + " over %s"%self.base_ring()

    def combinatorial_class(self):
        """
        Returns the combinatorial class that indexes the basis elements.

        Deprecated: use self.basis().keys() instead.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.combinatorial_class()
            doctest:...: DeprecationWarning: "FM.combinatorial_class()" is deprecated. Use "F.basis().keys()" instead !
            {'a', 'b', 'c'}

        ::

            sage: s = SFASchur(QQ)
            sage: s.combinatorial_class()
            Partitions
        """
        from sage.misc.misc import deprecation
        deprecation('"FM.combinatorial_class()" is deprecated. Use "F.basis().keys()" instead !')
        return self._basis_keys

    def dimension(self):
        """
        Returns the dimension of the combinatorial algebra (which is given
        by the number of elements in the basis).

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.dimension()
            3
            sage: F.basis().cardinality()
            3
            sage: F.basis().keys().cardinality()
            3

        ::

            sage: s = SFASchur(QQ)
            sage: s.dimension()
            +Infinity
        """
        return self._basis_keys.cardinality()

    def set_order(self, order):
        """
        Sets the order of the elements of the basis.

        If :meth:`set_order` has not been called, then the ordering is
        the one used in the generation of the elements of self's
        associated enumerated set.

        EXAMPLES::

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

        EXAMPLES::

            sage: QS2 = SymmetricGroupAlgebra(QQ,2)
            sage: QS2.get_order()                    # note: order changed on 2009-03-13
            [[2, 1], [1, 2]]
        """
        if self._order is None:
            self._order = list(self.basis().keys())
        return self._order

    def from_vector(self, vector):
        """
        Build an element of self from a (sparse) vector

        See self.get_order and the method to_vector of the elements of self

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: b = QS3.from_vector(vector((2, 0, 0, 0, 0, 4))); b
            2*[1, 2, 3] + 4*[3, 2, 1]
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a == b
            True
        """
        cc = self.get_order()
        return self._from_dict(dict( (cc[index], coeff) for (index,coeff) in vector.iteritems()))


    def prefix(self):
        """
        Returns the prefix used when displaying elements of self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.prefix()
            'B'

        ::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X.prefix()
            'X'
        """
        return self._prefix

    _repr_option_bracket = True

    def _repr_term(self, m):
        """
        Returns a string representing the basis elements indexed by m.

        The output can be customized by mean of:
         - self.prefix()
         - self._repr_option_bracket # TODO: find a good name
           (suggestions: _repr_with_bracket or _repr_add_bracket)

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: e = F.basis()
            sage: e['a'] + 2*e['b']    # indirect doctest
            B['a'] + 2*B['b']

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix="C")
            sage: e = F.basis()
            sage: e['a'] + 2*e['b']    # indirect doctest
            C['a'] + 2*C['b']

            sage: QS3 = CombinatorialFreeModule(QQ, Permutations(3), prefix="")
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a                      # indirect doctest
            2*[[1, 2, 3]] + 4*[[3, 2, 1]]

            sage: QS3._repr_option_bracket = False
            sage: a              # indirect doctest
            2*[1, 2, 3] + 4*[3, 2, 1]

            sage: QS3._repr_option_bracket = True
            sage: a              # indirect doctest
            2*[[1, 2, 3]] + 4*[[3, 2, 1]]

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), ('c','d')])
            sage: e = F.basis()
            sage: e[('a','b')] + 2*e[('c','d')]    # indirect doctest
            B[('a', 'b')] + 2*B[('c', 'd')]
        """
        if self._repr_option_bracket:
            return self.prefix()+"["+repr(m)+"]" # mind the (m), to accept a tuple for m
        else:
            return self.prefix()+repr(m)


    def __cmp__(self, other):
        """
        EXAMPLES::

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


        - ``x`` : a element of self

        - ``f``f - a function that takes in a combinatorial
          object indexing a basis element and returns an element of the
          target domain


        EXAMPLES::

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
        This takes in a function from the basis elements to the elements of
        self and applies it linearly to a. Note that
        _apply_module_endomorphism does not require multiplication on
        self to be defined.

        EXAMPLES::

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

    def monomial(self, index, coeff):
        """
        INPUT:
         - index: the index of a basis element
         - coeff: an element of the coefficient ring

        Returns the corresponding monomial

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.monomial('a',3)
            3*B['a']

        Design: should this do coercion on the coefficient ring?
        """
        return self._from_dict({index: coeff})

    def _term(self, index):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._term('a')
            B['a']
        """
        return self.monomial(index, self.base_ring().one())

    # This is generic, and should be lifted into modules_with_basis
    @lazy_attribute
    def term(self):
        """
        INPUT:

         - ''i''

        Returns the basis element indexed by i

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.term('a')
            B['a']

        F.term is in fact (almost) a map::

            sage: F.term
            Term map from {'a', 'b', 'c'} to Free module generated by {'a', 'b', 'c'} over Rational Field
        """
        # Should use a real Map, as soon as combinatorial_classes are enumerated sets, and therefore parents
        return PoorManMap(self._term, domain = self._basis_keys, codomain = self, name = "Term map")

    def _sum_of_term(self, indices):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._sum_of_term(['a', 'b'])
            B['a'] + B['b']
        """
        #TODO: optimize
        return self.sum(self.term(index) for index in indices)

    @lazy_attribute
    def sum_of_terms(self):
        """
        INPUT:

         - ''indices'' -- an list (or iterable) of indices of basis elements

        Returns the sum of the corresponding basis elements

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.sum_of_terms(['a', 'b'])
            B['a'] + B['b']

            sage: F.sum_of_terms(['a', 'b', 'a'])
            2*B['a'] + B['b']

        F.sum_of_term is in fact (almost) a map::

            sage: F.sum_of_terms
            A map to Free module generated by {'a', 'b', 'c'} over Rational Field
        """
        # domain = sets of self.combinatorial_class(),
        return PoorManMap(self._sum_of_term, codomain = self)

    def sum_of_monomials(self, monomials):
        """
        INPUT:
         - monomials: an list (or iterable) of pairs (index, coeff)

        Returns the sum of the monomials

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.sum_of_monomials([('a',2), ('c',3)])
            2*B['a'] + 3*B['c']

        Extreme case::

            sage: F.sum_of_monomials([])
            0

        TODO: optimize, especially for the special case where all
        indices are distinct (should there be an option distinct=True
        for this, or a separate function sum_of_distinct_monomials)?
        """
        return self.sum(self.monomial(index, coeff) for (index, coeff) in monomials)

    def term_or_zero_if_none(self, i):
	"""
	EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.term_or_zero_if_none('a')
            B['a']
	    sage: F.term_or_zero_if_none(None)
	    0
        """
	if i == None:
	    return self.zero()
	return self.term(i)

    @cached_method
    def zero(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.zero()
            0
        """
        return self._from_dict({})

    def _from_dict(self, d, coerce=False):
        """
        Given a monomial coefficient dictionary ``d``, returns the
        element of self with those monomials

        EXAMPLES::

            sage: e = SFAElementary(QQ)
            sage: s = SFASchur(QQ)
            sage: a = e([2,1]) + e([1,1,1]); a
            e[1, 1, 1] + e[2, 1]
            sage: s._from_dict(a.monomial_coefficients())
            s[1, 1, 1] + s[2, 1]

        If the optional argument ``coerce`` is ``True``, then the
        coefficients are coerced into the base ring of ``self``::

            sage: part = Partition([2,1])
            sage: d = {part:1}
            sage: a = s._from_dict(d,coerce=True); a
            s[2, 1]
            sage: a.coefficient(part).parent()
            Rational Field
        """
        assert isinstance(d, dict)
        if coerce:
            R = self.base_ring()
            # FIXME: this should remove zero entries
            d = dict((m,R(c)) for m,c in d.iteritems())

        return self.element_class(self, d)



class CombinatorialFreeModule(UniqueRepresentation, CombinatorialFreeModuleInterface, Module):
    r"""
    EXAMPLES:

    We construct a free module whose basis is indexed by the letters a, b, c::

        sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
        sage: F
        Free module generated by {'a', 'b', 'c'} over Rational Field

    Its basis is a family, indexed by a, b, c::

        sage: e = F.basis()
        sage: e
        Finite family {'a': B['a'], 'c': B['c'], 'b': B['b']}

    ::

        sage: [x for x in e]
        [B['a'], B['b'], B['c']]
        sage: [k for k in e.keys()]
        ['a', 'b', 'c']

    Let us construct some elements, and compute with them::

        sage: e['a']
        B['a']
        sage: 2*e['a']
        2*B['a']
        sage: e['a'] + 3*e['b']
        B['a'] + 3*B['b']

    Some uses of :meth:`.summation` and :meth:`.sum`::

        sage: F = CombinatorialFreeModule(QQ, [1,2,3,4])
        sage: F.summation(F.term(1), F.term(3))
        B[1] + B[3]

        sage: F = CombinatorialFreeModule(QQ, [1,2,3,4])
        sage: F.sum(F.term(i) for i in [1,2,3])
        B[1] + B[2] + B[3]

    Note that the constructed free module depend on the order of the basis::

        sage: F1 = CombinatorialFreeModule(QQ, (1,2,3,4))
        sage: F1 is F
        True
        sage: F1 = CombinatorialFreeModule(QQ, [4,3,2,1])
        sage: F1 == F
        False
    """

    @staticmethod
    def __classcall__(cls, *args, **keywords):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: G = CombinatorialFreeModule(QQ, ('a','b','c'))
            sage: F is G
            True
        """
        # convert the argument cc into a tuple if it is a list.
        # What's the proper way to handle this? we want to support
        # subclasses whose constructor takes some completely different
        # set of arguments.
        if cls is CombinatorialFreeModule:
            # note: if args is too short, we still propagate it down
            # to __init__ to let it handle proper exception raising.
            if len(args) >= 2 and (type(args[1]) is list or type(args[1]) is tuple):
                args = (args[0], FiniteEnumeratedSet(args[1])) + args[2:]
        return super(CombinatorialFreeModule, cls).__classcall__(cls, *args, **keywords)

    def __init__(self, R, basis_keys, element_class = None, prefix="B", category = None):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])

            sage: F.category()
            Category of modules with basis over Rational Field

          One may specify the category this module belongs to:
            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'], category=AlgebrasWithBasis(QQ))
            sage: F.category()
            Category of algebras with basis over Rational Field

          This is mostly intended for subclasses.

        TESTS::

            sage: F == loads(dumps(F))
            True
        """
        if not hasattr(self, "_name"):
            self._name = "Free module generated by %s"%(basis_keys,)  # note: cc may be a tuple!
        if isinstance(basis_keys, list):
            basis_keys = FiniteEnumeratedSet(basis_keys)
        self._basis_keys = basis_keys
        if not hasattr(self, "_prefix"):
            self._prefix = prefix
        CombinatorialFreeModuleInterface.__init__(self, R, element_class, category = category)


class CombinatorialFreeModule_Tensor(CombinatorialFreeModule):
        """
        Tensor Product of Free Modules

        EXAMPLES:

        We construct two free modules, assign them short names, and construct their tensor product::

            sage: F = CombinatorialFreeModule(ZZ, [1,2]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [3,4]); G.__custom_name = "G"
            sage: T = tensor([F, G]); T
            F # G

            sage: T.category()
            Category of tensor products of modules with basis over Integer Ring

            sage: T.construction() # todo: not implemented
            [tensor, ]

        T is a free module, with same base ring as F and G::

            sage: T.base_ring()
            Integer Ring

        The basis of T is indexed by tuples of basis indices of F and G::

            sage: T.basis().keys()
            Image of Cartesian product of {1, 2}, {3, 4} by <type 'tuple'>
            sage: T.basis().keys().list()
            [(1, 3), (1, 4), (2, 3), (2, 4)]

        FIXME: Should elements of a CartesianProduct be tuples (making them hashable)?

        Here are the basis elements themselves::

            sage: T.basis().cardinality()
            4
            sage: list(T.basis())
            [B[1] # B[3], B[1] # B[4], B[2] # B[3], B[2] # B[4]]

        The tensor product is associative and flattens sub tensor products::

            sage: H = CombinatorialFreeModule(ZZ, [5,6]); H.rename("H")
            sage: tensor([F, tensor([G, H])])
            F # G # H
            sage: tensor([tensor([F, G]), H])
            F # G # H
            sage: tensor([F, G, H])
            F # G # H

        We now compute the tensor product of elements of free modules::

            sage: f =   F.term(1) + 2 * F.term(2)
            sage: g = 2*G.term(3) +     G.term(4)
            sage: h =   H.term(5) +     H.term(6)
            sage: tensor([f, g])
            2*B[1] # B[3] + B[1] # B[4] + 4*B[2] # B[3] + 2*B[2] # B[4]

        Again, the tensor product is associative on elements::

            sage: tensor([f, tensor([g, h])]) == tensor([f, g, h])
            True
            sage: tensor([tensor([f, g]), h]) == tensor([f, g, h])
            True

        Note further that the tensor product spaces need not preexist::

            sage: t = tensor([f, g, h])
            sage: t.parent()
            F # G # H


        TESTS::

            sage: tensor([tensor([F, G]), H]) == tensor([F, G, H])
            True
            sage: tensor([F, tensor([G, H])]) == tensor([F, G, H])
            True
        """
        @staticmethod
        def __classcall__(cls, modules, **options):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2])
                sage: G = CombinatorialFreeModule(ZZ, [3,4])
                sage: H = CombinatorialFreeModule(ZZ, [4])
                sage: tensor([tensor([F, G]), H]) == tensor([F, G, H])
                True
                sage: tensor([F, tensor([G, H])]) == tensor([F, G, H])
                True
            """
            assert(len(modules) > 0)
            R = modules[0].base_ring()
            assert(all(module in ModulesWithBasis(R)) for module in modules)
            # should check the base ring
            # flatten the list of modules so that tensor(A, tensor(B,C)) gets rewritten into tensor(A, B, C)
            modules = sum([module.modules if isinstance(module, CombinatorialFreeModule_Tensor) else (module,) for module in modules], ())
            return super(CombinatorialFreeModule.Tensor, cls).__classcall__(cls, modules, **options)


        def __init__(self, modules, **options):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2]); F
                F
            """
            self.modules = modules
            CombinatorialFreeModule.__init__(self, modules[0].base_ring(), CartesianProduct(*[module.basis().keys() for module in modules]).map(tuple), **options)

        def _repr_(self):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3])
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,8])
                sage: tensor([F, G]) # indirect doctest
                Free module generated by {1, 2, 3} over Integer Ring # Free module generated by {1, 2, 3, 8} over Integer Ring
            """
            from sage.categories.tensor import tensor
            return tensor.symbol.join(["%s"%module for module in self.modules])
            # TODO: make this overridable by setting _name

        def _repr_term(self, term):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3])
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4])
                sage: f =   F.term(1) + 2 * F.term(2)
                sage: g = 2*G.term(3) +     G.term(4)
                sage: tensor([f, g]) # indirect doctest
                2*B[1] # B[3] + B[1] # B[4] + 4*B[2] # B[3] + 2*B[2] # B[4]
            """
            from sage.categories.tensor import tensor
            return tensor.symbol.join(module._repr_term(t) for (module, t) in zip(self.modules, term))

        # TODO: latex

        @cached_method
        def tensor_constructor(self, modules):
            """
            INPUT:

             - ``modules`` -- a tuple `(F_1,\dots,F_n)` of
               free modules whose tensor product is self

            Returns the canonical multilinear morphism from
            `F_1 \times \dots \times F_n` to `F_1 \otimes \dots \otimes F_n`

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ, [1,2]); F.__custom_name = "F"
                sage: G = CombinatorialFreeModule(ZZ, [3,4]); G.__custom_name = "G"
                sage: H = CombinatorialFreeModule(ZZ, [5,6]); H.rename("H")

                sage: f =   F.term(1) + 2 * F.term(2)
                sage: g = 2*G.term(3) +     G.term(4)
                sage: h =   H.term(5) +     H.term(6)

                sage: FG  = tensor([F, G   ])
                sage: phi_fg = FG.tensor_constructor((F, G))
                sage: phi_fg(f,g)
                2*B[1] # B[3] + B[1] # B[4] + 4*B[2] # B[3] + 2*B[2] # B[4]

                sage: FGH = tensor([F, G, H])
                sage: phi_fgh = FGH.tensor_constructor((F, G, H))
                sage: phi_fgh(f, g, h)
                2*B[1] # B[3] # B[5] + 2*B[1] # B[3] # B[6] + B[1] # B[4] # B[5] + B[1] # B[4] # B[6] + 4*B[2] # B[3] # B[5] + 4*B[2] # B[3] # B[6] + 2*B[2] # B[4] # B[5] + 2*B[2] # B[4] # B[6]

                sage: phi_fg_h = FGH.tensor_constructor((FG, H))
                sage: phi_fg_h(phi_fg(f, g), h)
                2*B[1] # B[3] # B[5] + 2*B[1] # B[3] # B[6] + B[1] # B[4] # B[5] + B[1] # B[4] # B[6] + 4*B[2] # B[3] # B[5] + 4*B[2] # B[3] # B[6] + 2*B[2] # B[4] # B[5] + 2*B[2] # B[4] # B[6]
            """
            assert(module in ModulesWithBasis(self.base_ring()) for module in modules)
            assert(sage.categories.tensor.tensor(modules) == self)
            # a list l such that l[i] is True if modules[i] is readily a tensor product
            is_tensor = [isinstance(module, CombinatorialFreeModule_Tensor) for module in modules]
            # the tensor_constructor, on basis elements
            result = self.term * CartesianProductWithFlattening(is_tensor)
            # TODO: make this into an element of Hom( A x B, C ) when those will exist
            for i in range(0, len(modules)):
                result = modules[i]._module_morphism(result, position = i, codomain = self)
            return result

        def _tensor_of_elements(self, elements):
            """
            Returns the tensor product of the specified elements.
            The result should be in self.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ, [1,2]); F.__custom_name = "F"
                sage: G = CombinatorialFreeModule(ZZ, [3,4]); G.__custom_name = "G"
                sage: H = CombinatorialFreeModule(ZZ, [5,6]); H.rename("H")

                sage: f =   F.term(1) + 2 * F.term(2)
                sage: g = 2*G.term(3) +     G.term(4)
                sage: h =   H.term(5) +     H.term(6)

                sage: GH  = tensor([G, H])
                sage: gh = GH._tensor_of_elements([g, h]); gh
                2*B[3] # B[5] + 2*B[3] # B[6] + B[4] # B[5] + B[4] # B[6]

                sage: FGH = tensor([F, G, H])
                sage: FGH._tensor_of_elements([f, g, h])
                2*B[1] # B[3] # B[5] + 2*B[1] # B[3] # B[6] + B[1] # B[4] # B[5] + B[1] # B[4] # B[6] + 4*B[2] # B[3] # B[5] + 4*B[2] # B[3] # B[6] + 2*B[2] # B[4] # B[5] + 2*B[2] # B[4] # B[6]

                sage: FGH._tensor_of_elements([f, gh])
                2*B[1] # B[3] # B[5] + 2*B[1] # B[3] # B[6] + B[1] # B[4] # B[5] + B[1] # B[4] # B[6] + 4*B[2] # B[3] # B[5] + 4*B[2] # B[3] # B[6] + 2*B[2] # B[4] # B[5] + 2*B[2] # B[4] # B[6]
            """
            return self.tensor_constructor(tuple(element.parent() for element in elements))(*elements)

class CartesianProductWithFlattening(object):
    """
    A class for cartesian product constructor, with partial flattening
    """
    def __init__(self, flatten):
        """
        INPUT:

         - ``flatten`` -- a tuple of booleans

        This constructs a callable which accepts ``len(flatten)``
        arguments, and builds a tuple out them. When ``flatten[i]``,
        the i-th argument itself should be a tuple which is flattened
        in the result.

            sage: from sage.combinat.free_module import CartesianProductWithFlattening
            sage: CartesianProductWithFlattening([True, False, True, True])
            <sage.combinat.free_module.CartesianProductWithFlattening object at ...>

        """
        self._flatten = flatten

    def __call__(self, *indices):
        """
        EXAMPLES::

            sage: from sage.combinat.free_module import CartesianProductWithFlattening
            sage: cp = CartesianProductWithFlattening([True, False, True, True])
            sage: cp((1,2), (3,4), (5,6), (7,8))
            (1, 2, (3, 4), 5, 6, 7, 8)
            sage: cp((1,2,3), 4, (5,6), (7,8))
            (1, 2, 3, 4, 5, 6, 7, 8)

        """
        return sum( (i if flatten else (i,) for (i,flatten) in zip(indices, self._flatten) ), ())


# TODO: find a way to avoid this hack to allow for cross references
CombinatorialFreeModule.Tensor = CombinatorialFreeModule_Tensor


class CombinatorialFreeModule_CartesianProduct(CombinatorialFreeModule):
    """
    An implementation of cartesian products of modules with basis

    EXAMPLES:

    We construct two free modules, assign them short names, and construct their cartesian product::

        sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
        sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
        sage: H = CombinatorialFreeModule(ZZ, [4,7]); H.__custom_name = "H"
        sage: S = cartesian_product([F, G])
        sage: S
        F (+) G
        sage: S.basis()
        Lazy family (Term map from Disjoint union of Family ({4, 5}, {4, 6}) to F (+) G(i))_{i in Disjoint union of Family ({4, 5}, {4, 6})}

    Note that the indices of the basis elements of F and G intersect non
    trivially. This is handled by forcing the union to be disjoint::

        sage: list(S.basis())
        [B[(0, 4)], B[(0, 5)], B[(1, 4)], B[(1, 6)]]

    We now compute the cartesian product of elements of free modules::

        sage: f =   F.term(4) + 2 * F.term(5)
        sage: g = 2*G.term(4) +     G.term(6)
        sage: h =   H.term(4) +     H.term(7)
        sage: cartesian_product([f,g])
        B[(0, 4)] + 2*B[(0, 5)] + 2*B[(1, 4)] + B[(1, 6)]
        sage: cartesian_product([f,g,h])
        B[(0, 4)] + 2*B[(0, 5)] + 2*B[(1, 4)] + B[(1, 6)] + B[(2, 4)] + B[(2, 7)]
        sage: cartesian_product([f,g,h]).parent()
        F (+) G (+) H

    TODO: choose an appropriate semantic for cartesian products of cartesian products (associativity?)::

        sage: S = cartesian_product([cartesian_product([F, G]), H]) # todo: not implemented
        F (+) G (+) H
    """

    def __init__(self, modules, **options):
        r"""
        TESTS::

            sage: F = CombinatorialFreeModule(ZZ, [2,4,5])
            sage: G = CombinatorialFreeModule(ZZ, [2,4,7])
            sage: cartesian_product([F, G])
            Free module generated by {2, 4, 5} over Integer Ring (+) Free module generated by {2, 4, 7} over Integer Ring
        """
        assert(len(modules) > 0) # TODO: generalize to a family or tuple
        R = modules[0].base_ring()
        assert(all(module in ModulesWithBasis(R)) for module in modules)
        # should check the base ring
        self.modules = modules
        CombinatorialFreeModule.__init__(self, R,
            DisjointUnionEnumeratedSets(
                [module.basis().keys() for module in modules], keepkey=True),
            **options)

    def modules_keys(self):
        """
        In waiting for self.modules.keys()

        TESTS::

            sage: F = CombinatorialFreeModule(ZZ, [2,4,5])
            sage: G = CombinatorialFreeModule(ZZ, [2,4,7])
            sage: CP = cartesian_product([F, G])
            sage: CP.modules_keys()
            [0, 1]
        """
        return range(len(self.modules))

    def _repr_(self):
        """
        TESTS::

        sage: F = CombinatorialFreeModule(ZZ, [2,4,5])
        sage: CP = cartesian_product([F, F]); CP  # indirect doctest
        Free module generated by {2, 4, 5} over Integer Ring (+) Free module generated by {2, 4, 5} over Integer Ring
        sage: F.__custom_name = "F"; CP
        F (+) F
        """
        from sage.categories.cartesian_product import cartesian_product
        return cartesian_product.symbol.join(["%s"%module for module in self.modules])
        # TODO: make this overridable by setting _name

    @cached_method
    def summand_embedding(self, i):
        """
        Returns the natural embedding morphism of the i-th summand of self into self

        INPUTS:

         - ``i`` -- an integer

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: phi = S.summand_embedding(0)
            sage: phi(F.term(4) + 2 * F.term(5))
            B[(0, 4)] + 2*B[(0, 5)]
            sage: phi(F.term(4) + 2 * F.term(6)).parent() == S
            True
            sage: phi(G.term(4)) # not implemented Should raise an error! problem: G(F.term(4)) does not complain!!!!
        """
        assert i in self.modules_keys()
        return self.modules[i]._module_morphism(lambda t: self.term((i,t)), codomain = self)

    @cached_method
    def summand_projection(self, i):
        """
        Returns the natural projection onto the i-th summand of self

        INPUTS:

         - ``i`` -- an integer

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: x = S.term((0,4)) + 2 * S.term((0,5)) + 3 * S.term((1,6))
            sage: S.summand_projection(0)(x)
            B[4] + 2*B[5]
            sage: S.summand_projection(1)(x)
            3*B[6]
            sage: S.summand_projection(0)(x).parent() == F
            True
            sage: S.summand_projection(1)(x).parent() == G
            True
        """
        assert i in self.modules_keys()
        module = self.modules[i]
        return self._module_morphism(lambda (j,t): module.term(t) if i == j else module.zero(), codomain = module)

    def _cartesian_product_of_elements(self, elements):
        """
        Returns the cartesian product of the elements

        INPUT:

         - ``elements`` - a tuple with one element of each summand of self

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: f =   F.term(4) + 2 * F.term(5)
            sage: g = 2*G.term(4) +     G.term(6)
            sage: S._cartesian_product_of_elements([f, g])
            B[(0, 4)] + 2*B[(0, 5)] + 2*B[(1, 4)] + B[(1, 6)]
            sage: S._cartesian_product_of_elements([f, g]).parent() == S
            True

        """
        return self.sum(self.summand_embedding(i)(elements[i]) for i in self.modules_keys())

    class Element(CombinatorialFreeModule.Element): # TODO: get rid of this inheritance
        def summand_projection(self, i):
            """
            Returns the projection of self on the i-th summand of this cartesian product
            INPUTS:

             - ``i`` -- an integer

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
                sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
                sage: S = cartesian_product([F, G])
                sage: x = S.term((0,4)) + 2 * S.term((0,5)) + 3 * S.term((1,6))
                sage: x.summand_projection(0)
                B[4] + 2*B[5]
                sage: x.summand_projection(1)
                3*B[6]
            """
            return self.parent().summand_projection(i)(self)

        def summand_split(self):
            """
            Splits x into its summands

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
                sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
                sage: H = CombinatorialFreeModule(ZZ, [4,7]); H.__custom_name = "H"
                sage: S = cartesian_product([F, G, H])
                sage: x = S.term((0,4)) + 2 * S.term((0,5)) + 3 * S.term((1,6)) + 4 * S.term((2,4)) + 5 * S.term((2,7))
                sage: x.summand_split()
                (B[4] + 2*B[5], 3*B[6], 4*B[4] + 5*B[7])
                sage: [s.parent() for s in x.summand_split()]
                [F, G, H]
                sage: S.zero().summand_split()
                (0, 0, 0)
                sage: [s.parent() for s in S.zero().summand_split()]
                [F, G, H]
            """
            # TODO: optimize
            return tuple(self.summand_projection(i) for i in self.parent().modules_keys())
            #return Family(self.modules.keys(), self.projection)

CombinatorialFreeModule.CartesianProduct = CombinatorialFreeModule_CartesianProduct
