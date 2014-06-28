"""
Free modules
"""
#*****************************************************************************
#       Copyright (C) 2007      Mike Hansen <mhansen@gmail.com>,
#                     2007-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                     2010      Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.sage_object import have_same_parent
from sage.modules.free_module_element import vector
from sage.misc.misc import repr_lincomb
from sage.modules.module import Module
from sage.rings.all import Integer
import sage.structure.element
from sage.combinat.family import Family
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.combinat.cartesian_product import CartesianProduct
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.misc.all import lazy_attribute
from sage.categories.poor_man_map import PoorManMap
from sage.categories.all import ModulesWithBasis
from sage.combinat.dict_addition import dict_addition, dict_linear_combination
from sage.sets.family import Family
from sage.misc.ascii_art import AsciiArt, empty_ascii_art

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

            sage: s = SymmetricFunctions(QQ).schur()
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

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1]) + s([3])
            sage: Partition([2,1]) in a
            True
            sage: Partition([1,1,1]) in a
            False
        """
        return x in self._monomial_coefficients and self._monomial_coefficients[x] != 0

    @cached_method
    def __hash__(self):
        """
        Return the hash value for ``self``.

        The result is cached.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: hash(f)
            6429418278783588506           # 64-bit
            726440090                     # 32-bit

            sage: F = RootSystem(['A',2]).ambient_space()
            sage: f = F.simple_root(0)
            sage: hash(f)
            6920829894162680369           # 64-bit
            -528971215                    # 32-bit

        This uses the recipe that was proposed for frozendicts in `PEP
        0416 <http://legacy.python.org/dev/peps/pep-0416/>`_ (and adds
        the hash of the parent). This recipe relies on the hash
        function for frozensets which uses tricks to mix the hash
        values of the items in case they are similar.

        .. TODO::

            It would be desirable to make the hash value depend on the
            hash value of the parent. See :trac:`15959`.
        """
        return hash(frozenset(self._monomial_coefficients.items()))

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

        To run through the monomials of an element, it is better to
        use the idiom::

            sage: for (t,c) in f:
            ...       print t,c
            a 1
            c 3

        ::

            sage: s = SymmetricFunctions(QQ).schur()
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

    def _sorted_items_for_printing(self):
        """
        Returns the items (i.e terms) of ``self``, sorted for printing

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c'] + 3 * B['b']
            sage: f._sorted_items_for_printing()
            [('a', 1), ('b', 3), ('c', 2)]
            sage: F.print_options(monomial_cmp = lambda x,y: -cmp(x,y))
            sage: f._sorted_items_for_printing()
            [('c', 2), ('b', 3), ('a', 1)]
            sage: F.print_options(monomial_cmp=cmp) #reset to original state

        .. seealso:: :meth:`_repr_`, :meth:`_latex_`, :meth:`print_options`
        """
        print_options = self.parent().print_options()
        v = self._monomial_coefficients.items()
        try:
            v.sort(cmp = print_options['monomial_cmp'],
                   key = lambda monomial_coeff: monomial_coeff[0])
        except Exception: # Sorting the output is a plus, but if we can't, no big deal
            pass
        return v

    def _repr_(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='F')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            F['a'] + 2*F['b']
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            ['a'] + 2*['b']
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='', scalar_mult=' ', bracket=False)
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            'a' + 2 'b'

        Controling the order of terms by providing a comparison
        function on elements of the support::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'],
            ...                               monomial_cmp = lambda x,y: cmp(y,x))
            sage: e = F.basis()
            sage: e['a'] + 3*e['b'] + 2*e['c']
            2*B['c'] + 3*B['b'] + B['a']

            sage: F = CombinatorialFreeModule(QQ, ['ac', 'ba', 'cb'],
            ...                               monomial_cmp = lambda x,y: cmp(x[1],y[1]))
            sage: e = F.basis()
            sage: e['ac'] + 3*e['ba'] + 2*e['cb']
            3*B['ba'] + 2*B['cb'] + B['ac']
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self.parent()._print_options['scalar_mult'],
                            repr_monomial = self.parent()._repr_term,
                            strip_one = True)

    def _ascii_art_(self):
        """
        TESTS::

            sage: M = QuasiSymmetricFunctions(QQ).M()
            sage: ascii_art(M[1,3]**2)
            4*M      + 2*M       + 2*M      + 2*M       + 2*M       + M
                 ***      ******        ***         ***         ***     ******
               ***        *             *        ****         ***      **
               *          *           ***        *           **
               *                      *
        """
        from sage.misc.misc import coeff_repr
        terms = self._sorted_items_for_printing()
        scalar_mult = self.parent()._print_options['scalar_mult']
        repr_monomial = self.parent()._ascii_art_term
        strip_one = True

        if repr_monomial is None:
            repr_monomial = str

        s = empty_ascii_art # ""
        first = True
        i = 0

        if scalar_mult is None:
            scalar_mult = "*"

        all_atomic = True
        for (monomial,c) in terms:
            b = repr_monomial(monomial) # PCR
            if c != 0:
                break_points = []
                coeff = coeff_repr(c, False)
                if coeff != "0":
                    if coeff == "1":
                        coeff = ""
                    elif coeff == "-1":
                        coeff = "-"
                    elif b._l > 0:
                        if len(coeff) > 0 and monomial == 1 and strip_one:
                            b = empty_ascii_art # ""
                        else:
                            b = AsciiArt([scalar_mult]) + b
                    if not first:
                        if len(coeff) > 0 and coeff[0] == "-":
                            coeff = " - %s"%coeff[1:]
                        else:
                            coeff = " + %s"%coeff
                        break_points = [2]
                    else:
                        coeff = "%s"%coeff
                s += AsciiArt([coeff], break_points) + b
                first = False
        if first:
            return "0"
        elif s == empty_ascii_art:
            return AsciiArt(["1"])
        else:
            return s

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
            2[1, 2, 3] + [2, 1, 3]

       ::

            sage: F = CombinatorialFreeModule(QQ, ['a','b'], prefix='beta', latex_prefix='\\beta')
            sage: x = F.an_element()
            sage: x
            2*beta['a'] + 2*beta['b']
            sage: latex(x)
            2\beta_{a} + 2\beta_{b}

        Controling the order of terms by providing a comparison
        function on elements of the support::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'],
            ...                               monomial_cmp = lambda x,y: cmp(y,x))
            sage: e = F.basis()
            sage: latex(e['a'] + 3*e['b'] + 2*e['c'])
            2B_{c} + 3B_{b} + B_{a}

            sage: F = CombinatorialFreeModule(QQ, ['ac', 'ba', 'cb'],
            ...                               monomial_cmp = lambda x,y: cmp(x[1],y[1]))
            sage: e = F.basis()
            sage: latex(e['ac'] + 3*e['ba'] + 2*e['cb'])
            3B_{ba} + 2B_{cb} + B_{ac}
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult       = self.parent()._print_options['scalar_mult'],
                            latex_scalar_mult = self.parent()._print_options['latex_scalar_mult'],
                            repr_monomial = self.parent()._latex_term,
                            is_latex=True, strip_one = True)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: F1 = CombinatorialFreeModule(QQ, [1, 2, 3])
            sage: F2 = CombinatorialFreeModule(QQ, [1, 2, 3], prefix = "g")
            sage: F1.zero() == F1.zero()
            True
            sage: F1.zero() == F1.an_element()
            False
            sage: F1.an_element() == F1.an_element()
            True
            sage: F1.an_element() is None
            False

        .. TODO::

            Currently, if ``self`` and ``other`` do not have the same parent,
            seemingly equal elements do not evaluate equal, since conversions
            between different modules have not been established.

        ::

            sage: F1.zero() == 0
            True
            sage: F1(0)
            0

        ::

            sage: F1.zero() == F2.zero()
            False
            sage: F1(F2.zero())
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= 0) an element of self (=Free module generated by {1, 2, 3} over Rational Field)
            sage: F = AlgebrasWithBasis(QQ).example()
            sage: F.one() == 1
            True
            sage: 1 == F.one()
            True
            sage: 2 * F.one() == int(2)
            True
            sage: int(2) == 2 * F.one()
            True

            sage: S = SymmetricFunctions(QQ); s = S.s(); p = S.p()
            sage: p[2] == s[2] - s[1, 1]
            True
            sage: p[2] == s[2]
            False

        This feature is disputable, in particular since it can make
        equality testing costly. It may be removed at some point.

        Equality testing can be a bit tricky when the order of terms
        can vary because their indices are incomparable with
        ``cmp``. The following test did fail before :trac:`12489` ::

            sage: F = CombinatorialFreeModule(QQ, Subsets([1,2,3]))
            sage: x = F.an_element()
            sage: (x+F.zero()).terms()  # random
            [2*B[{1}], 3*B[{2}], B[{}]]
            sage: x.terms()             # random
            [2*B[{1}], B[{}], 3*B[{2}]]
            sage: x+F.zero() == x
            True

        TESTS::

            sage: TestSuite(F1).run()
            sage: TestSuite(F).run()
        """
        if have_same_parent(self, other):
            return self._monomial_coefficients == other._monomial_coefficients
        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.eq)
        except TypeError:
            return False

    def __ne__(left, right):
        """
        EXAMPLES::

            sage: F1 = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F1.an_element() != F1.an_element()
            False
            sage: F1.an_element() != F1.zero()
            True
        """
        return not left == right

    def __cmp__(left, right):
        """
        The ordering is the one on the underlying sorted list of
        (monomial,coefficients) pairs.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1])
            sage: b = s([1,1,1])
            sage: cmp(a,b) #indirect doctest
            1
        """
        if have_same_parent(left, right) and left._monomial_coefficients == right._monomial_coefficients:
            return 0
        v = sorted(mc for mc in left._monomial_coefficients.items() if mc[1] != 0)
        w = sorted(mc for mc in right._monomial_coefficients.items() if mc[1] != 0)
        return cmp(v, w)

    def _add_(self, other):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a'] + 3*B['c']
            B['a'] + 3*B['c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: s([2,1]) + s([5,4]) # indirect doctest
            s[2, 1] + s[5, 4]
            sage: a = s([2,1]) + 0
            sage: len(a.monomial_coefficients())
            1
        """

        assert hasattr( other, 'parent' ) and other.parent() == self.parent()

        F = self.parent()
        return F._from_dict( dict_addition( [ self._monomial_coefficients, other._monomial_coefficients ] ), remove_zeros=False )

    def _neg_(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: -f
            -B['a'] - 3*B['c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: -s([2,1]) # indirect doctest
            -s[2, 1]
        """
        F = self.parent()
        return F._from_dict( dict_linear_combination( [ ( self._monomial_coefficients, -1 ) ] ), remove_zeros=False )

    def _sub_(self, other):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a'] - 3*B['c']
            B['a'] - 3*B['c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: s([2,1]) - s([5,4]) # indirect doctest
            s[2, 1] - s[5, 4]
        """
        assert hasattr( other, 'parent' ) and other.parent() == self.parent()
        F = self.parent()
        return F._from_dict( dict_linear_combination( [ ( self._monomial_coefficients, 1 ), (other._monomial_coefficients, -1 ) ] ), remove_zeros=False )

    def _coefficient_fast(self, m, default=None):
        """
        Returns the coefficient of m in self, where m is key in
        self._monomial_coefficients.

        EXAMPLES::

            sage: p = Partition([2,1])
            sage: q = Partition([1,1,1])
            sage: s = SymmetricFunctions(QQ).schur()
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
            sage: (2*F.monomial(g)).coefficient(g)
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

            sage: s = SymmetricFunctions(QQ).schur()
            sage: s([2,1]).is_zero()
            False
            sage: s(0).is_zero()
            True
            sage: (s([2,1]) - s([2,1])).is_zero()
            True
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )
        return all( v == zero for v in self._monomial_coefficients.values() )

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

            sage: s = SymmetricFunctions(QQ).schur()
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

            sage: s = SymmetricFunctions(QQ).schur()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.length()
            4
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )
        return len( [ key for key, coeff in self._monomial_coefficients.iteritems() if coeff != zero ] )

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

            sage: s = SymmetricFunctions(QQ).schur()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.support()
            [[1], [1, 1, 1], [2, 1], [4]]
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )

        supp = sorted([ key for key, coeff in self._monomial_coefficients.iteritems() if coeff != zero ])

        return supp

    def monomials(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c']
            sage: f.monomials()
            [B['a'], B['c']]

            sage: (F.zero()).monomials()
            []
        """
        P = self.parent()
        BR = P.base_ring()
        zero = BR( 0 )
        one = BR( 1 )

        supp = sorted([ key for key, coeff in self._monomial_coefficients.iteritems() if coeff != zero ])

        return [ P._from_dict( { key : one }, remove_zeros=False ) for key in supp ]

    def terms(self):
        """
        Returns a list of the terms of ``self``

        .. seealso:: :meth:`monomials`

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c']
            sage: f.terms()
            [B['a'], 2*B['c']]
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )
        v = sorted([ ( key, value ) for key, value in self._monomial_coefficients.iteritems() if value != zero ])
        from_dict = self.parent()._from_dict
        return [ from_dict( { key : value } ) for key,value in v ]

    def coefficients(self):
        """
        Returns a list of the coefficients appearing on the basis elements in
        self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.coefficients()
            [1, -3]

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficients()
            [1, 1, 1, 1]
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )
        v = sorted([ ( key, value ) for key, value in self._monomial_coefficients.iteritems() if value != zero ])
        return [ value for key,value in v ]

    def _vector_(self, new_base_ring=None):
        """
        Returns ``self`` as a dense vector

        INPUT:

        - ``new_base_ring`` -- a ring (default: None)

        OUTPUT: a dense :func:`FreeModule` vector

        .. WARNING:: This will crash/run forever if ``self`` is infinite dimensional!

        .. SEEALSO::

            - :func:`vector`
            - :meth:`CombinatorialFreeModule.get_order`
            - :meth:`CombinatorialFreeModule.from_vector`
            - :meth:`CombinatorialFreeModule._dense_free_module`

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f._vector_()
            (1, 0, -3)

        One can use equivalently::

            sage: f.to_vector()
            (1, 0, -3)
            sage: vector(f)
            (1, 0, -3)

        More examples::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a._vector_()
            (2, 0, 0, 0, 0, 4)
            sage: a.to_vector()
            (2, 0, 0, 0, 0, 4)
            sage: vector(a)
            (2, 0, 0, 0, 0, 4)
            sage: a == QS3.from_vector(a.to_vector())
            True

        If ''new_base_ring'' is specified, then a vector over
        ''new_base_ring'' is returned::

            sage: a._vector_(RDF)
            (2.0, 0.0, 0.0, 0.0, 0.0, 4.0)

        .. NOTE::

            #13406: the current implementation has been optimized, at
            #the price of breaking the encapsulation for FreeModule
            elements creation, with the following use case as metric,
            on a 2008' Macbook Pro::

                sage: F = CombinatorialFreeModule(QQ, range(10))
                sage: f = F.an_element()
                sage: %timeit f._vector_()   # not tested
                625 loops, best of 3: 17.5 micros per loop

             Other use cases may call for different or further
             optimizations.
        """
        parent = self.parent()
        dense_free_module = parent._dense_free_module(new_base_ring)
        d = self._monomial_coefficients
        return dense_free_module._element_class(dense_free_module,
                                                [d.get(m, 0) for m in parent.get_order()],
                                                coerce=True, copy=False)

    to_vector = _vector_

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
            Right action by Integer Ring on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: F.get_action(F, operator.mul, True)
            sage: F.get_action(F, operator.mul, False)

        This also works when a coercion of the coefficient is needed, for
        example with polynomials or fraction fields #8832::

            sage: P.<q> = QQ['q']
            sage: V = CombinatorialFreeModule(P, Permutations())
            sage: el = V(Permutation([3,1,2]))
            sage: (3/2)*el
            3/2*B[[3, 1, 2]]

            sage: P.<q> = QQ['q']
            sage: F = FractionField(P)
            sage: V = CombinatorialFreeModule(F, Words())
            sage: w = Words()('abc')
            sage: (1+q)*V(w)
            (q+1)*B[word: abc]
            sage: ((1+q)/q)*V(w)
            ((q+1)/q)*B[word: abc]

        TODO:
         - add non commutative tests
        """
        # With the current design, the coercion model does not have
        # enough information to detect apriori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        if hasattr( scalar, 'parent' ) and scalar.parent() != self.base_ring():
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if self.base_ring().has_coerce_map_from(scalar.parent()):
                scalar = self.base_ring()( scalar )
            else:
                return None

        F = self.parent()
        D = self._monomial_coefficients
        if self_on_left:
            D = dict_linear_combination( [ ( D, scalar ) ], factor_on_left = False )
        else:
            D = dict_linear_combination( [ ( D, scalar ) ] )

        return F._from_dict( D, remove_zeros=False )

    # For backward compatibility
    _lmul_ = _acted_upon_
    _rmul_ = _acted_upon_

    def __div__(self, x, self_on_left=False ):
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
            F = self.parent()
            x = self.base_ring()( x )
            x_inv = x**-1
            D = self._monomial_coefficients
            if self_on_left:
                D = dict_linear_combination( [ ( D, x_inv ) ], factor_on_left=False )
            else:
                D = dict_linear_combination( [ ( D, x_inv ) ] )

            return F._from_dict( D, remove_zeros=False )
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
        raise ValueError("%s is not divisible by %s"%(x, y))
    else:
        return q

class CombinatorialFreeModule(UniqueRepresentation, Module):
    r"""
    Class for free modules with a named basis

    INPUT:

    - ``R`` - base ring

    - ``basis_keys`` - list, tuple, family, set, etc. defining the
      indexing set for the basis of this module

    - ``element_class`` - the class of which elements of this module
      should be instances (optional, default None, in which case the
      elements are instances of
      :class:`CombinatorialFreeModuleElement`)

    - ``category`` - the category in which this module lies (optional,
      default None, in which case use the "category of modules with
      basis" over the base ring ``R``)

    Options controlling the printing of elements:

    - ``prefix`` - string, prefix used for printing elements of this
      module (optional, default 'B').  With the default, a monomial
      indexed by 'a' would be printed as ``B['a']``.

    - ``latex_prefix`` - string or None, prefix used in the LaTeX
      representation of elements (optional, default None). If this is
      anything except the empty string, it prints the index as a
      subscript.  If this is None, it uses the setting for ``prefix``,
      so if ``prefix`` is set to "B", then a monomial indexed by 'a'
      would be printed as ``B_{a}``.  If this is the empty string, then
      don't print monomials as subscripts: the monomial indexed by 'a'
      would be printed as ``a``, or as ``[a]`` if ``latex_bracket`` is
      True.

    - ``bracket`` - None, bool, string, or list or tuple of
      strings (optional, default None): if None, use the value of the
      attribute ``self._repr_option_bracket``, which has default value
      True.  (``self._repr_option_bracket`` is available for backwards
      compatibility.  Users should set ``bracket`` instead.  If
      ``bracket`` is set to anything except None, it overrides
      the value of ``self._repr_option_bracket``.)  If False, do not
      include brackets when printing elements: a monomial indexed by
      'a' would be printed as ``B'a'``, and a monomial indexed by
      (1,2,3) would be printed as ``B(1,2,3)``.  If True, use "[" and
      "]" as brackets.  If it is one of "[", "(", or "{", use it and
      its partner as brackets.  If it is any other string, use it as
      both brackets.  If it is a list or tuple of strings, use the
      first entry as the left bracket and the second entry as the
      right bracket.

    - ``latex_bracket`` - bool, string, or list or tuple of strings
      (optional, default False): if False, do not include brackets in
      the LaTeX representation of elements.  This option is only
      relevant if ``latex_prefix`` is the empty string; otherwise,
      brackets are not used regardless.  If True, use "\\left[" and
      "\\right]" as brackets.  If this is one of "[", "(", "\\{", "|",
      or "||", use it and its partner, prepended with "\\left" and
      "\\right", as brackets.  If this is any other string, use it as
      both brackets.  If this is a list or tuple of strings, use the
      first entry as the left bracket and the second entry as the
      right bracket.

    - ``scalar_mult`` - string to use for scalar multiplication in
      the print representation (optional, default "*")

    - ``latex_scalar_mult`` - string or None (optional, default None),
      string to use for scalar multiplication in the latex
      representation.  If None, use the empty string if ``scalar_mult``
      is set to "*", otherwise use the value of ``scalar_mult``.

    - ``tensor_symbol`` - string or None (optional, default None),
      string to use for tensor product in the print representation. If
      None, use the ``sage.categories.tensor.symbol``.

    - ``monomial_cmp`` - a comparison function (optional, default cmp),
      to use for sorting elements in the output of elements

    .. note:: These print options may also be accessed and modified using the
       :meth:`print_options` method, after the module has been defined.

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

    Some uses of
    :meth:`sage.categories.commutative_additive_semigroups.CommutativeAdditiveSemigroups.ParentMethods.summation`
    and :meth:`.sum`::

        sage: F = CombinatorialFreeModule(QQ, [1,2,3,4])
        sage: F.summation(F.monomial(1), F.monomial(3))
        B[1] + B[3]

        sage: F = CombinatorialFreeModule(QQ, [1,2,3,4])
        sage: F.sum(F.monomial(i) for i in [1,2,3])
        B[1] + B[2] + B[3]

    Note that free modules with a given basis and parameters are unique::

        sage: F1 = CombinatorialFreeModule(QQ, (1,2,3,4))
        sage: F1 is F
        True

    The identity of the constructed free module depends on the order of the
    basis and on the other parameters, like the prefix. Note that :class:`CombinatorialFreeModule` is
    a :class:`~sage.structure.unique_representation.UniqueRepresentation`. Hence,
    two combinatorial free modules evaluate equal if and only if they are
    identical::

        sage: F1 = CombinatorialFreeModule(QQ, (1,2,3,4))
        sage: F1 is F
        True
        sage: F1 = CombinatorialFreeModule(QQ, [4,3,2,1])
        sage: F1 == F
        False
        sage: F2 = CombinatorialFreeModule(QQ, [1,2,3,4], prefix='F')
        sage: F2 == F
        False

    Because of this, if you create a free module with certain parameters and
    then modify its prefix or other print options, this affects all modules
    which were defined using the same parameters.
    ::

        sage: F2.print_options(prefix='x')
        sage: F2.prefix()
        'x'
        sage: F3 = CombinatorialFreeModule(QQ, [1,2,3,4], prefix='F')
        sage: F3 is F2   # F3 was defined just like F2
        True
        sage: F3.prefix()
        'x'
        sage: F4 = CombinatorialFreeModule(QQ, [1,2,3,4], prefix='F', bracket=True)
        sage: F4 == F2   # F4 was NOT defined just like F2
        False
        sage: F4.prefix()
        'F'

        sage: F2.print_options(prefix='F') #reset for following doctests

    The default category is the category of modules with basis over
    the base ring::

        sage: CombinatorialFreeModule(GF(3), ((1,2), (3,4))).category()
        Category of vector spaces with basis over Finite Field of size 3

    See :mod:`sage.categories.examples.algebras_with_basis` and
    :mod:`sage.categories.examples.hopf_algebras_with_basis` for
    illustrations of the use of the ``category`` keyword, and see
    :class:`sage.combinat.root_system.weight_space.WeightSpace` for an
    example of the use of ``element_class``.

    Customizing print and LaTeX representations of elements::

        sage: F = CombinatorialFreeModule(QQ, ['a','b'], prefix='x')
        sage: original_print_options = F.print_options()
        sage: sorted(original_print_options.items())
         [('bracket', None), ('latex_bracket', False), ('latex_prefix', None), ('latex_scalar_mult', None), ('monomial_cmp', <built-in function cmp>), ('prefix', 'x'), ('scalar_mult', '*'), ('tensor_symbol', None)]

        sage: e = F.basis()
        sage: e['a'] - 3 * e['b']
        x['a'] - 3*x['b']

        sage: F.print_options(prefix='x', scalar_mult=' ', bracket='{')
        sage: e['a'] - 3 * e['b']
        x{'a'} - 3 x{'b'}
        sage: latex(e['a'] - 3 * e['b'])
        x_{a} - 3 x_{b}

        sage: F.print_options(latex_prefix='y')
        sage: latex(e['a'] - 3 * e['b'])
        y_{a} - 3  y_{b}

        sage: F.print_options(monomial_cmp = lambda x,y: -cmp(x,y))
        sage: e['a'] - 3 * e['b']
        -3 x{'b'} + x{'a'}
        sage: F.print_options(**original_print_options) # reset print options

        sage: F = CombinatorialFreeModule(QQ, [(1,2), (3,4)])
        sage: e = F.basis()
        sage: e[(1,2)] - 3 * e[(3,4)]
        B[(1, 2)] - 3*B[(3, 4)]

        sage: F.print_options(bracket=['_{', '}'])
        sage: e[(1,2)] - 3 * e[(3,4)]
        B_{(1, 2)} - 3*B_{(3, 4)}

        sage: F.print_options(prefix='', bracket=False)
        sage: e[(1,2)] - 3 * e[(3,4)]
        (1, 2) - 3*(3, 4)

    TESTS:

    Before :trac:`14054`, combinatorial free modules violated the unique
    parent condition. That caused a problem. The tensor product construction
    involves maps, but maps check that their domain and the parent of a
    to-be-mapped element are identical (not just equal). However, the tensor
    product was cached by a :class:`~sage.misc.cachefunc.cached_method`, which
    involves comparison by equality (not identity). Hence, the last line of
    the following example used to fail with an assertion error::

        sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix="F")
        sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4], prefix="G")
        sage: f =   F.monomial(1) + 2 * F.monomial(2)
        sage: g = 2*G.monomial(3) +     G.monomial(4)
        sage: tensor([f, g])
        2*F[1] # G[3] + F[1] # G[4] + 4*F[2] # G[3] + 2*F[2] # G[4]
        sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix='x')
        sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4], prefix='y')
        sage: f =   F.monomial(1) + 2 * F.monomial(2)
        sage: g = 2*G.monomial(3) +     G.monomial(4)
        sage: tensor([f, g])
        2*x[1] # y[3] + x[1] # y[4] + 4*x[2] # y[3] + 2*x[2] # y[4]
    """

    @staticmethod
    def __classcall_private__(cls, base_ring, basis_keys, category = None, **keywords):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: G = CombinatorialFreeModule(QQ, ('a','b','c'))
            sage: F is G
            True

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'], latex_bracket=['LEFT', 'RIGHT'])
            sage: F.print_options()['latex_bracket']
            ('LEFT', 'RIGHT')

            sage: F is G
            False

        We check that the category is properly straightened::

            sage: F  = CombinatorialFreeModule(QQ, ['a','b'])
            sage: F1 = CombinatorialFreeModule(QQ, ['a','b'], category = ModulesWithBasis(QQ))
            sage: F2 = CombinatorialFreeModule(QQ, ['a','b'], category = [ModulesWithBasis(QQ)])
            sage: F3 = CombinatorialFreeModule(QQ, ['a','b'], category = (ModulesWithBasis(QQ),))
            sage: F4 = CombinatorialFreeModule(QQ, ['a','b'], category = (ModulesWithBasis(QQ),CommutativeAdditiveSemigroups()))
            sage: F5 = CombinatorialFreeModule(QQ, ['a','b'], category = (ModulesWithBasis(QQ),Category.join((LeftModules(QQ), RightModules(QQ)))))
            sage: F1 is F, F2 is F, F3 is F, F4 is F, F5 is F
            (True, True, True, True, True)

            sage: G  = CombinatorialFreeModule(QQ, ['a','b'], category = AlgebrasWithBasis(QQ))
            sage: F is G
            False
        """
        if isinstance(basis_keys, (list, tuple)):
            basis_keys = FiniteEnumeratedSet(basis_keys)
        category = ModulesWithBasis(base_ring).or_subcategory(category)
        # bracket or latex_bracket might be lists, so convert
        # them to tuples so that they're hashable.
        bracket = keywords.get('bracket', None)
        if isinstance(bracket, list):
            keywords['bracket'] = tuple(bracket)
        latex_bracket = keywords.get('latex_bracket', None)
        if isinstance(latex_bracket, list):
            keywords['latex_bracket'] = tuple(latex_bracket)
        return super(CombinatorialFreeModule, cls).__classcall__(cls, base_ring, basis_keys, category = category, **keywords)

    Element = CombinatorialFreeModuleElement

    def __init__(self, R, basis_keys, element_class = None, category = None, **kwds):
        r"""
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])

            sage: F.category()
            Category of vector spaces with basis over Rational Field

        One may specify the category this module belongs to::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'], category=AlgebrasWithBasis(QQ))
            sage: F.category()
            Category of algebras with basis over Rational Field

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'], category = FiniteDimensionalModulesWithBasis(QQ))
            sage: F.basis()
            Finite family {'a': B['a'], 'c': B['c'], 'b': B['b']}
            sage: F.category()
            Category of finite dimensional vector spaces with basis over Rational Field
            sage: TestSuite(F).run()

        TESTS:

        Regression test for #10127: ``self._basis_keys`` needs to be
        set early enough, in case the initialization of the categories
        use ``self.basis().keys()``. This occured on several occasions
        in non trivial constructions. In the following example,
        :class:`AlgebrasWithBasis` constructs ``Homset(self,self)`` to
        extend by bilinearity method ``product_on_basis``, which in
        turn triggers ``self._repr_()``::

            sage: class MyAlgebra(CombinatorialFreeModule):
            ...       def _repr_(self):
            ...           return "MyAlgebra on %s"%(self.basis().keys())
            ...       def product_on_basis(self,i,j):
            ...           pass
            sage: MyAlgebra(ZZ, ZZ, category = AlgebrasWithBasis(QQ))
            MyAlgebra on Integer Ring

        A simpler example would be welcome!

        We check that unknown options are caught::

            sage: CombinatorialFreeModule(ZZ, [1,2,3], keyy=2)
            Traceback (most recent call last):
            ...
            ValueError: keyy is not a valid print option.
        """
        #Make sure R is a ring with unit element
        from sage.categories.all import Rings
        if R not in Rings():
            raise TypeError("Argument R must be a ring.")

        if category is None:
            category = ModulesWithBasis(R)

        if element_class is not None:
            self.Element = element_class

        # The following is needed by e.g. root systems that don't call
        # the classcall and passes lists as basis_keys
        if isinstance(basis_keys, (list, tuple)):
            basis_keys = FiniteEnumeratedSet(basis_keys)
        self._basis_keys = basis_keys # Needs to be done early: #10127

        Parent.__init__(self, base = R, category = category,
                        # Could we get rid of this?
                        element_constructor = self._element_constructor_)


        self._order = None

        # printing options for elements (set when initializing self).
        # This includes self._repr_option_bracket (kept for backwards
        # compatibility, declared to be True by default, needs to be
        # overridden explicitly).
        self._print_options = {'prefix': "B",
                               'bracket': None,
                               'latex_bracket': False,
                               'latex_prefix': None,
                               'scalar_mult': "*",
                               'latex_scalar_mult': None,
                               'tensor_symbol': None,
                               'monomial_cmp': cmp}
        # 'bracket': its default value here is None, meaning that
        # the value of self._repr_option_bracket is used; the default
        # value of that attribute is True -- see immediately before
        # the method _repr_term.  If 'bracket' is any value
        # except None, then it overrides the value of
        # self._repr_option_bracket.  Future users might consider
        # using 'bracket' instead of _repr_option_bracket.

        # ignore the optional 'key' since it only affects CachedRepresentation
        kwds.pop('key', None)
        self.print_options(**kwds)

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
            x = x + self.monomial(I.an_element())
        except Exception:
            pass
        try:
            g = iter(self.basis().keys())
            for c in range(1,4):
                x = x + self.term(g.next(), R(c))
        except Exception:
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
            sage: F.monomial("a") in F
            True
            sage: G.monomial("a") in F
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
            sage: F(F.monomial("a"))   # indirect doctest
            B['a']

        Do not rely on the following feature which may be removed in the future::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3([2,3,1])     # indirect doctest
            [2, 3, 1]

        instead, use::

            sage: P = QS3.basis().keys()
            sage: QS3.monomial(P([2,3,1]))   # indirect doctest
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
            sage: G(F.monomial("a"))
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= B['a']) an element of self (=G)
            sage: H(F.monomial("a"))
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= B['a']) an element of self (=H)

        Here is a real life example illustrating that this yielded
        mathematically wrong results::

            sage: S = SymmetricFunctions(QQ)
            sage: s = S.s(); p = S.p()
            sage: ss = tensor([s,s]); pp = tensor([p,p])
            sage: a = tensor((s[2],s[2]))

        The following originally used to yield ``p[[2]] # p[[2]]``, and if
        there was no natural coercion between ``s`` and ``p``, this would
        raise a ``NotImplementedError``. Since :trac:`15305`, this takes the
        coercion between ``s`` and ``p`` and lifts it to the tensor product. ::

            sage: pp(a)
            1/4*p[1, 1] # p[1, 1] + 1/4*p[1, 1] # p[2] + 1/4*p[2] # p[1, 1] + 1/4*p[2] # p[2]

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
                raise TypeError("do not know how to make x (= %s) an element of %s"%(x, self))
        #x is an element of the basis enumerated set;
        # This is a very ugly way of testing this
        elif ((hasattr(self._basis_keys, 'element_class') and
               isinstance(self._basis_keys.element_class, type) and
               isinstance(x, self._basis_keys.element_class))
              or (sage.structure.element.parent(x) == self._basis_keys)):
            return self.monomial(x)
        elif x in self._basis_keys:
            return self.monomial(self._basis_keys(x))
        else:
            if hasattr(self, '_coerce_end'):
                try:
                    return self._coerce_end(x)
                except TypeError:
                    pass
            raise TypeError("do not know how to make x (= %s) an element of self (=%s)"%(x,self))

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

            sage: s = SymmetricFunctions(QQ).schur()
            sage: s.dimension()
            +Infinity
        """
        return self._basis_keys.cardinality()

    def gens(self):
        """
        Return a tuple consisting of the basis elements of ``self``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, ['a', 'b', 'c'])
            sage: F.gens()
            (B['a'], B['b'], B['c'])
        """
        return tuple(self.basis().values())

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

    @cached_method
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

    @cached_method
    def _dense_free_module(self, base_ring=None):
        """
        Returns a dense free module of the same dimension

        INPUT:

        - ``base_ring`` -- a ring or ``None``

        If ``base_ring`` is None, then the base ring of ``self`` is used.

        This method is mostly used by
        :meth:`CombinatorialFreeModule.Element._vector_`

        .. SEEALSO:: :meth:`from_vector`

        EXAMPLES::

            sage: C = CombinatorialFreeModule(QQ['x'], ['a','b','c']); C
            Free module generated by {'a', 'b', 'c'} over Univariate Polynomial Ring in x over Rational Field
            sage: C._dense_free_module()
            Ambient free module of rank 3 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: C._dense_free_module(QQ['x,y'])
            Ambient free module of rank 3 over the integral domain Multivariate Polynomial Ring in x, y over Rational Field
        """
        if base_ring is None:
            base_ring = self.base_ring()
        from sage.modules.free_module import FreeModule
        return FreeModule(base_ring, self.dimension())

    def from_vector(self, vector):
        """
        Build an element of ``self`` from a (sparse) vector

        .. SEEALSO:: :meth:`get_order`, :meth:`CombinatorialFreeModule.Element._vector_`

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
        return self._print_options['prefix']

    def print_options(self, **kwds):
        """
        Return the current print options, or set an option.

        INPUT: all of the input is optional; if present, it should be
        in the form of keyword pairs, such as
        ``latex_bracket='('``.  The allowable keywords are:

        - ``prefix``
        - ``latex_prefix``
        - ``bracket``
        - ``latex_bracket``
        - ``scalar_mult``
        - ``latex_scalar_mult``
        - ``tensor_symbol``
        - ``monomial_cmp``

        See the documentation for :class:`CombinatorialFreeModule` for
        descriptions of the effects of setting each of these options.

        OUTPUT: if the user provides any input, set the appropriate
        option(s) and return nothing.  Otherwise, return the
        dictionary of settings for print and LaTeX representations.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix='x')
            sage: F.print_options()
            {...'prefix': 'x'...}
            sage: F.print_options(bracket='(')
            sage: F.print_options()
            {...'bracket': '('...}

        TESTS::

            sage: sorted(F.print_options().items())
            [('bracket', '('), ('latex_bracket', False), ('latex_prefix', None), ('latex_scalar_mult', None), ('monomial_cmp', <built-in function cmp>), ('prefix', 'x'), ('scalar_mult', '*'), ('tensor_symbol', None)]
            sage: F.print_options(bracket='[') # reset
        """
        # don't just use kwds.get(...) because I want to distinguish
        # between an argument like "option=None" and the option not
        # being there altogether.
        if kwds:
            for option in kwds:
                if option in ['prefix', 'latex_prefix', 'bracket', 'latex_bracket',
                              'scalar_mult', 'latex_scalar_mult', 'tensor_symbol',
                              'monomial_cmp'
                             ]:
                    self._print_options[option] = kwds[option]
                else:
                    raise ValueError('%s is not a valid print option.' % option)
        else:
            return self._print_options

    _repr_option_bracket = True

    def _repr_term(self, m):
        """
        Returns a string representing the basis element indexed by m.

        The output can be customized by setting any of the following
        options when initializing the module:

        - prefix
        - bracket
        - scalar_mult

        Alternatively, one can use the :meth:`print_options` method
        to achieve the same effect.  To modify the bracket setting,
        one can also set ``self._repr_option_bracket`` as long as one
        has *not* set the ``bracket`` option: if the
        ``bracket`` option is anything but ``None``, it overrides
        the value of ``self._repr_option_bracket``.

        See the documentation for :class:`CombinatorialFreeModule` for
        details on the initialization options.

        .. todo:: rename to ``_repr_monomial``

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: e = F.basis()
            sage: e['a'] + 2*e['b']    # indirect doctest
            B['a'] + 2*B['b']

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix="F")
            sage: e = F.basis()
            sage: e['a'] + 2*e['b']    # indirect doctest
            F['a'] + 2*F['b']

            sage: QS3 = CombinatorialFreeModule(QQ, Permutations(3), prefix="")
            sage: original_print_options = QS3.print_options()
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a                      # indirect doctest
            2*[[1, 2, 3]] + 4*[[3, 2, 1]]

            sage: QS3.print_options(bracket = False)
            sage: a              # indirect doctest
            2*[1, 2, 3] + 4*[3, 2, 1]

            sage: QS3.print_options(prefix='')
            sage: a              # indirect doctest
            2*[1, 2, 3] + 4*[3, 2, 1]

            sage: QS3.print_options(bracket="|", scalar_mult=" *@* ")
            sage: a              # indirect doctest
            2 *@* |[1, 2, 3]| + 4 *@* |[3, 2, 1]|

            sage: QS3.print_options(**original_print_options) # reset

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), ('c','d')])
            sage: e = F.basis()
            sage: e[('a','b')] + 2*e[('c','d')]    # indirect doctest
            B[('a', 'b')] + 2*B[('c', 'd')]
        """
        bracket = self._print_options.get('bracket', None)
        bracket_d = {"{": "}", "[": "]", "(": ")"}
        if bracket is None:
            bracket = self._repr_option_bracket
        if bracket is True:
            left = "["
            right = "]"
        elif bracket is False:
            left = ""
            right = ""
        elif isinstance(bracket, (tuple, list)):
            left = bracket[0]
            right = bracket[1]
        elif bracket in bracket_d:
            left = bracket
            right = bracket_d[bracket]
        else:
            left = bracket
            right = bracket
        return self.prefix() + left + repr(m) + right # mind the (m), to accept a tuple for m

    def _ascii_art_term(self, el):
        r"""
        Return an ascii art representing of the term.

        TESTS::

            sage: R = NonCommutativeSymmetricFunctions(QQ).R()
            sage: ascii_art(R[1,2,2,4])
            R
               ****
              **
             **
             *
            sage: Partitions.global_options(diagram_str="#", convention="french")
            sage: ascii_art(R[1,2,2,4])
            R
             #
             ##
              ##
               ####
        """
        from sage.misc.ascii_art import ascii_art
        try:
            if el == self.one_basis():
                return AsciiArt(["1"])
        except Exception:
            pass
        pref = AsciiArt([self.prefix()])
        r = pref * (AsciiArt([" "**Integer(len(pref))]) + ascii_art(el))
        r._baseline = r._h - 1
        return r

    def _latex_term(self, m):
        r"""
        Returns a string for the LaTeX code for the basis element
        indexed by m.

        The output can be customized by setting any of the following
        options when initializing the module:

        - prefix
        - latex_prefix
        - latex_bracket

        (Alternatively, one can use the :meth:`print_options` method
        to achieve the same effect.)

        See the documentation for :class:`CombinatorialFreeModule` for
        details on the initialization options.

        .. todo:: rename to ``_latex_monomial``

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: e = F.basis()
            sage: latex(e['a'] + 2*e['b'])    # indirect doctest
            B_{a} + 2B_{b}

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix="C")
            sage: e = F.basis()
            sage: latex(e['a'] + 2*e['b'])    # indirect doctest
            C_{a} + 2C_{b}

            sage: QS3 = CombinatorialFreeModule(QQ, Permutations(3), prefix="", scalar_mult="*")
            sage: original_print_options = QS3.print_options()
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: latex(a)                     # indirect doctest
            2[1, 2, 3] + 4[3, 2, 1]
            sage: QS3.print_options(latex_bracket=True)
            sage: latex(a)                     # indirect doctest
            2\left[ [1, 2, 3] \right] + 4\left[ [3, 2, 1] \right]
            sage: QS3.print_options(latex_bracket="(")
            sage: latex(a)                     # indirect doctest
            2\left( [1, 2, 3] \right) + 4\left( [3, 2, 1] \right)
            sage: QS3.print_options(latex_bracket=('\\myleftbracket', '\\myrightbracket'))
            sage: latex(a)                     # indirect doctest
            2\myleftbracket [1, 2, 3] \myrightbracket + 4\myleftbracket [3, 2, 1] \myrightbracket
            sage: QS3.print_options(**original_print_options) # reset

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), (0,1,2)])
            sage: e = F.basis()
            sage: latex(e[('a','b')])    # indirect doctest
            B_{('a', 'b')}
            sage: latex(2*e[(0,1,2)])    # indirect doctest
            2B_{\left(0, 1, 2\right)}
            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), (0,1,2)], prefix="")
            sage: e = F.basis()
            sage: latex(2*e[(0,1,2)])    # indirect doctest
            2\left(0, 1, 2\right)
        """
        from sage.misc.latex import latex

        s = latex(m)
        if s.find('\\text{\\textt') != -1:
            # m contains "non-LaTeXed" strings, use string representation
            s = str(m)

        # dictionary with left-right pairs of "brackets".  put pairs
        # in here accept \\left and \\right as prefixes.
        bracket_d = {"{": "\\}", "[": "]", "(": ")", "\\{": "\\}",
                     "|": "|", "||": "||"}
        bracket = self._print_options.get('latex_bracket', False)
        if bracket is True:
            left = "\\left["
            right = "\\right]"
        elif bracket is False:
            left = ""
            right = ""
        elif isinstance(bracket, (tuple, list)):
            left = bracket[0]
            right = bracket[1]
        elif bracket in bracket_d:
            left = bracket
            right = bracket_d[bracket]
            if left == "{":
                left = "\\{"
            left = "\\left" + left
            right = "\\right" + right
        else:
            left = bracket
            right = bracket
        prefix = self._print_options.get('latex_prefix')
        if prefix is None:
            prefix = self._print_options.get('prefix')
        if prefix == "":
            return left + s + right
        return "%s_{%s}" % (prefix, s)

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

    def _apply_module_morphism( self, x, on_basis, codomain=False ):
        """
        Returns the image of ``x`` under the module morphism defined by
        extending :func:`on_basis` by linearity.

        INPUT:


        - ``x`` : a element of self

        - ``on_basis`` - a function that takes in a combinatorial
          object indexing a basis element and returns an element of the
          codomain

        - ``codomain`` (optional) - the codomain of the morphism, otherwise it is computed
          using :func:`on_basis`

        If ``codomain`` is not specified, then the function tries to compute the codomain
        of the module morphism by finding the image of one of the elements in the
        support, hence :func:`on_basis` should return an element whose parent is the
        codomain.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([3]) + s([2,1]) + s([1,1,1])
            sage: b = 2*a
            sage: f = lambda part: Integer( len(part) )
            sage: s._apply_module_morphism(a, f) #1+2+3
            6
            sage: s._apply_module_morphism(b, f) #2*(1+2+3)
            12
            sage: s._apply_module_morphism(s(0), f)
            0
            sage: s._apply_module_morphism(s(1), f)
            0
            sage: s._apply_module_morphism(s(1), lambda part: len(part), ZZ)
            0
            sage: s._apply_module_morphism(s(1), lambda part: len(part))
            Traceback (most recent call last):
            ...
            ValueError: Codomain could not be determined
        """

        if x == self.zero():
            if not codomain:
                B = Family(self.basis())
                try:
                    z = B.first()
                except StopIteration:
                    raise ValueError('Codomain could not be determined')
                codomain = on_basis(z).parent()
            return codomain.zero()
        else:
            if not codomain:
                keys = x.support()
                key = keys[0]
                if hasattr(on_basis( key ), 'parent'):
                    codomain = on_basis( key ).parent()
                else:
                    raise ValueError('Codomain could not be determined')

            if hasattr( codomain, 'linear_combination' ):
                return codomain.linear_combination( ( on_basis( key ), coeff ) for key, coeff in x._monomial_coefficients.iteritems() )
            else:
                return_sum = codomain.zero()
                for key, coeff in x._monomial_coefficients.iteritems():
                    return_sum += coeff * on_basis( key )
                return return_sum

    def _apply_module_endomorphism(self, x, on_basis):
        """
        This takes in a function from the basis elements to the elements of
        self and applies it linearly to a. Note that
        _apply_module_endomorphism does not require multiplication on
        self to be defined.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: f = lambda part: 2*s(part.conjugate())
            sage: s._apply_module_endomorphism( s([2,1]) + s([1,1,1]), f)
            2*s[2, 1] + 2*s[3]
        """
        return self.linear_combination( ( on_basis( key ), coeff ) for key, coeff in x._monomial_coefficients.iteritems() )

    def sum(self, iter_of_elements):
        """
        overrides method inherited from commutative additive monoid as it is much faster on dicts directly

        INPUT:

        - ``iter_of_elements``: iterator of elements of ``self``

        Returns the sum of all elements in ``iter_of_elements``

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ,[1,2])
            sage: f = F.an_element(); f
            2*B[1] + 2*B[2]
            sage: F.sum( f for _ in range(5) )
            10*B[1] + 10*B[2]
        """

        D = dict_addition( element._monomial_coefficients for element in iter_of_elements )
        return self._from_dict( D, remove_zeros=False )

    def linear_combination( self, iter_of_elements_coeff, factor_on_left=True ):
        """
        INPUT:

        - ``iter_of_elements_coeff`` -- iterator of pairs ``(element, coeff)``
          with ``element`` in ``self`` and ``coeff`` in ``self.base_ring()``

        - ``factor_on_left`` (optional) -- if ``True``, the coefficients are
          multiplied from the left if ``False``, the coefficients are
          multiplied from the right

        Returns the linear combination `\lambda_1 v_1 + ... + \lambda_k v_k`
        (resp.  the linear combination `v_1 \lambda_1 + ... + v_k \lambda_k`)
        where ``iter_of_elements_coeff`` iterates through the sequence
        `(\lambda_1, v_1) ... (\lambda_k, v_k)`.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ,[1,2])
            sage: f = F.an_element(); f
            2*B[1] + 2*B[2]
            sage: F.linear_combination( (f,i) for i in range(5) )
            20*B[1] + 20*B[2]
        """
        return self._from_dict( dict_linear_combination( ( ( element._monomial_coefficients, coeff ) for element, coeff in iter_of_elements_coeff ), factor_on_left=factor_on_left ), remove_zeros=False )

    def term(self, index, coeff=None):
        """
        Constructs a term in ``self``

        INPUT:

        - ``index`` -- the index of a basis element
        - ``coeff`` -- an element of the coefficient ring (default: one)

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.term('a',3)
            3*B['a']
            sage: F.term('a')
            B['a']

        Design: should this do coercion on the coefficient ring?
        """
        if coeff is None:
            coeff = self.base_ring().one()
        return self._from_dict( {index: coeff} )

    def _monomial(self, index):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._monomial('a')
            B['a']
        """
        return self._from_dict( {index: self.base_ring().one()}, remove_zeros = False )

    # This is generic, and should be lifted into modules_with_basis
    @lazy_attribute
    def monomial(self):
        """
        INPUT:

         - ''i''

        Returns the basis element indexed by i

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.monomial('a')
            B['a']

        F.monomial is in fact (almost) a map::

            sage: F.monomial
            Term map from {'a', 'b', 'c'} to Free module generated by {'a', 'b', 'c'} over Rational Field
        """
        # Should use a real Map, as soon as combinatorial_classes are enumerated sets, and therefore parents
        return PoorManMap(self._monomial, domain = self._basis_keys, codomain = self, name = "Term map")

    def _sum_of_monomials(self, indices):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._sum_of_monomials(['a', 'b'])
            B['a'] + B['b']
        """
        # TODO: optimize by calling directly _from_dict if we
        # know that all indices are distinct as sum_of_terms;
        # otherwise, maybe call dict_addition directly
        return self.sum(self.monomial(index) for index in indices)

    @lazy_attribute
    def sum_of_monomials(self):
        """
        INPUT:

         - ''indices'' -- an list (or iterable) of indices of basis elements

        Returns the sum of the corresponding basis elements

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.sum_of_monomials(['a', 'b'])
            B['a'] + B['b']

            sage: F.sum_of_monomials(['a', 'b', 'a'])
            2*B['a'] + B['b']

        F.sum_of_monomials is in fact (almost) a map::

            sage: F.sum_of_monomials
            A map to Free module generated by {'a', 'b', 'c'} over Rational Field
        """
        # domain = sets of self.combinatorial_class(),
        return PoorManMap(self._sum_of_monomials, codomain = self)

    def sum_of_terms(self, terms, distinct=False):
        """
        Constructs a sum of terms of ``self``

        INPUT:

        - ``terms`` -- a list (or iterable) of pairs (index, coeff)
        - ``distinct`` -- whether the indices are guaranteed to be distinct (default: ``False``)

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.sum_of_terms([('a',2), ('c',3)])
            2*B['a'] + 3*B['c']

        If ``distinct`` is True, then the construction is optimized::

            sage: F.sum_of_terms([('a',2), ('c',3)], distinct = True)
            2*B['a'] + 3*B['c']

        .. warning::

            Use ``distinct=True`` only if you are sure that the
            indices are indeed distinct::

                sage: F.sum_of_terms([('a',2), ('a',3)], distinct = True)
                3*B['a']

        Extreme case::

            sage: F.sum_of_terms([])
            0
        """
        if distinct:
            return self._from_dict(dict(terms))
        return self.sum(self.term(index, coeff) for (index, coeff) in terms)

    def monomial_or_zero_if_none(self, i):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.monomial_or_zero_if_none('a')
            B['a']
            sage: F.monomial_or_zero_if_none(None)
            0
        """
        if i is None:
            return self.zero()
        return self.monomial(i)

    @cached_method
    def zero(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.zero()
            0
        """
        return self._from_dict({}, remove_zeros=False)

    def _from_dict( self, d, coerce=False, remove_zeros=True ):
        """
        Construct an element of ``self`` from an `{index: coefficient}` dictionary

        INPUT:

        - ``d`` -- a dictionary ``{index: coeff}`` where each ``index`` is the
          index of a basis element and each ``coeff`` belongs to the
          coefficient ring ``self.base_ring()``

        - ``coerce`` -- a boolean (default: ``False``), whether to coerce the
          ``coeff``s to the coefficient ring

        - ``remove_zeros`` -- a boolean (default: ``True``), if some
          ``coeff``s may be zero and should therefore be removed

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).elementary()
            sage: s = SymmetricFunctions(QQ).schur()
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

        With ``remove_zeros=True``, zero coefficients are removed::

            sage: s._from_dict({part:0})
            0

        .. warning::

            With ``remove_zeros=True``, it is assumed that no
            coefficient of the dictionary is zero. Otherwise, this may
            lead to illegal results::

                sage: list(s._from_dict({part:0}, remove_zeros=False))
                [([2, 1], 0)]
        """
        assert isinstance(d, dict)
        if coerce:
            R = self.base_ring()
            d = dict( (key, R(coeff)) for key,coeff in d.iteritems())
        if remove_zeros:
            d = dict( (key, coeff) for key, coeff in d.iteritems() if coeff)
        return self.element_class( self, d )


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

            sage: f =   F.monomial(1) + 2 * F.monomial(2)
            sage: g = 2*G.monomial(3) +     G.monomial(4)
            sage: h =   H.monomial(5) +     H.monomial(6)
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
        def __classcall_private__(cls, modules, **options):
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
            modules = sum([module._sets if isinstance(module, CombinatorialFreeModule_Tensor) else (module,) for module in modules], ())
            return super(CombinatorialFreeModule.Tensor, cls).__classcall__(cls, modules, **options)


        def __init__(self, modules, **options):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2]); F
                F
            """
            from sage.categories.tensor import tensor
            self._sets = modules
            CombinatorialFreeModule.__init__(self, modules[0].base_ring(), CartesianProduct(*[module.basis().keys() for module in modules]).map(tuple), **options)
            # the following is not the best option, but it's better than nothing.
            self._print_options['tensor_symbol'] = options.get('tensor_symbol', tensor.symbol)

        def _repr_(self):
            """
            This is customizable by setting
            ``self.print_options('tensor_symbol'=...)``.

            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3])
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,8])
                sage: F.rename("F")
                sage: G.rename("G")
                sage: T = tensor([F, G])
                sage: T # indirect doctest
                F # G
                sage: T.print_options(tensor_symbol= ' @ ')  # note the spaces
                sage: T # indirect doctest
                F @ G

            To avoid a side\--effect on another doctest, we revert the change::

                sage: T.print_options(tensor_symbol= ' # ')
            """
            from sage.categories.tensor import tensor
            if hasattr(self, "_print_options"):
                symb = self._print_options['tensor_symbol']
                if symb is None:
                    symb = tensor.symbol
            else:
                symb = tensor.symbol
            return symb.join(["%s"%module for module in self._sets])
            # TODO: make this overridable by setting _name

        def _ascii_art_(self, term):
            """
            TESTS::

                sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                sage: Partitions.global_options(diagram_str="#", convention="french")
                sage: ascii_art(tensor((R[1,2], R[3,1,2])))
                R   # R
                 #     ###
                 ##      #
                         ##
            """
            from sage.categories.tensor import tensor
            if hasattr(self, "_print_options"):
                symb = self._print_options['tensor_symbol']
                if symb is None:
                    symb = tensor.symbol
            else:
                symb = tensor.symbol
            it = iter(zip(self._sets, term))
            module, t = it.next()
            rpr = module._ascii_art_term(t)
            for (module,t) in it:
                rpr += AsciiArt([symb], [len(symb)])
                rpr += module._ascii_art_term(t)
            return rpr

        _ascii_art_term = _ascii_art_

        def _latex_(self):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3])
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,8])
                sage: F.rename("F")
                sage: G.rename("G")
                sage: latex(tensor([F, F, G])) # indirect doctest
                \text{\texttt{F}} \otimes \text{\texttt{F}} \otimes \text{\texttt{G}}
                sage: F._latex_ = lambda : "F"
                sage: G._latex_ = lambda : "G"
                sage: latex(tensor([F, F, G])) # indirect doctest
                F \otimes F \otimes G
            """
            from sage.misc.latex import latex
            symb = " \\otimes "
            return symb.join(["%s"%latex(module) for module in self._sets])

        def _repr_term(self, term):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix="F")
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4], prefix="G")
                sage: f =   F.monomial(1) + 2 * F.monomial(2)
                sage: g = 2*G.monomial(3) +     G.monomial(4)
                sage: tensor([f, g]) # indirect doctest
                2*F[1] # G[3] + F[1] # G[4] + 4*F[2] # G[3] + 2*F[2] # G[4]
            """
            from sage.categories.tensor import tensor
            if hasattr(self, "_print_options"):
                symb = self._print_options['tensor_symbol']
                if symb is None:
                    symb = tensor.symbol
            else:
                symb = tensor.symbol
            return symb.join(module._repr_term(t) for (module, t) in zip(self._sets, term))

        def _latex_term(self, term):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix='x')
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4], prefix='y')
                sage: f =   F.monomial(1) + 2 * F.monomial(2)
                sage: g = 2*G.monomial(3) +     G.monomial(4)
                sage: latex(tensor([f, g])) # indirect doctest
                2x_{1} \otimes y_{3} + x_{1} \otimes y_{4} + 4x_{2} \otimes y_{3} + 2x_{2} \otimes y_{4}
            """
            symb = " \\otimes "
            return symb.join(module._latex_term(t) for (module, t) in zip(self._sets, term))

        @cached_method
        def tensor_constructor(self, modules):
            r"""
            INPUT:

             - ``modules`` -- a tuple `(F_1,\dots,F_n)` of
               free modules whose tensor product is self

            Returns the canonical multilinear morphism from
            `F_1 \times \dots \times F_n` to `F_1 \otimes \dots \otimes F_n`

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ, [1,2]); F.__custom_name = "F"
                sage: G = CombinatorialFreeModule(ZZ, [3,4]); G.__custom_name = "G"
                sage: H = CombinatorialFreeModule(ZZ, [5,6]); H.rename("H")

                sage: f =   F.monomial(1) + 2 * F.monomial(2)
                sage: g = 2*G.monomial(3) +     G.monomial(4)
                sage: h =   H.monomial(5) +     H.monomial(6)
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
            result = self.monomial * CartesianProductWithFlattening(is_tensor) #.
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

                sage: f =   F.monomial(1) + 2 * F.monomial(2)
                sage: g = 2*G.monomial(3) +     G.monomial(4)
                sage: h =   H.monomial(5) +     H.monomial(6)

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

        def _coerce_map_from_(self, R):
            """
            Return ``True`` if there is a coercion from ``R`` into ``self`` and
            ``False`` otherwise.  The things that coerce into ``self`` are:

            - Anything with a coercion into ``self.base_ring()``.

            - A tensor algebra whose factors have a coercion into the
              corresponding factors of ``self``.

            TESTS::

                sage: C = CombinatorialFreeModule(ZZ, ZZ)
                sage: C2 = CombinatorialFreeModule(ZZ, NN)
                sage: M = C.module_morphism(lambda x: C2.monomial(abs(x)), codomain=C2)
                sage: M.register_as_coercion()
                sage: C2(C.basis()[3])
                B[3]
                sage: C2(C.basis()[3] + C.basis()[-3])
                2*B[3]
                sage: S = C.tensor(C)
                sage: S2 = C2.tensor(C2)
                sage: S2.has_coerce_map_from(S)
                True
                sage: S.has_coerce_map_from(S2)
                False
                sage: S.an_element()
                3*B[0] # B[-1] + 2*B[0] # B[0] + 2*B[0] # B[1]
                sage: S2(S.an_element())
                2*B[0] # B[0] + 5*B[0] # B[1]

            ::

                sage: C = CombinatorialFreeModule(ZZ, Set([1,2]))
                sage: D = CombinatorialFreeModule(ZZ, Set([2,4]))
                sage: f = C.module_morphism(on_basis=lambda x: D.monomial(2*x), codomain=D)
                sage: f.register_as_coercion()
                sage: T = tensor((C,C))
                sage: p = D.an_element()
                sage: T(tensor((p,p)))
                Traceback (most recent call last):
                ...
                NotImplementedError
                sage: T = tensor((D,D))
                sage: p = C.an_element()
                sage: T(tensor((p,p)))
                4*B[2] # B[2] + 4*B[2] # B[4] + 4*B[4] # B[2] + 4*B[4] # B[4]
            """
            if R in ModulesWithBasis(self.base_ring()).TensorProducts() \
                    and isinstance(R, CombinatorialFreeModule_Tensor) \
                    and len(R._sets) == len(self._sets) \
                    and all(self._sets[i].has_coerce_map_from(M)
                            for i,M in enumerate(R._sets)):
                modules = R._sets
                vector_map = [self._sets[i]._internal_coerce_map_from(M)
                              for i,M in enumerate(modules)]
                return R.module_morphism(lambda x: self._tensor_of_elements(
                        [vector_map[i](M.monomial(x[i]))
                         for i,M in enumerate(modules)]),
                                         codomain=self)

            return super(CombinatorialFreeModule_Tensor, self)._coerce_map_from_(R)

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

        sage: f =   F.monomial(4) + 2 * F.monomial(5)
        sage: g = 2*G.monomial(4) +     G.monomial(6)
        sage: h =   H.monomial(4) +     H.monomial(7)
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
            sage: C = cartesian_product([F, G]) ; C
            Free module generated by {2, 4, 5} over Integer Ring (+) Free module generated by {2, 4, 7} over Integer Ring
            sage: TestSuite(C).run()
        """
        assert(len(modules) > 0) # TODO: generalize to a family or tuple
        R = modules[0].base_ring()
        assert(all(module in ModulesWithBasis(R)) for module in modules)
        # should check the base ring
        self._sets = modules
        CombinatorialFreeModule.__init__(self, R,
            DisjointUnionEnumeratedSets(
                [module.basis().keys() for module in modules], keepkey=True),
            **options)

    def _sets_keys(self):
        """
        In waiting for self._sets.keys()

        TESTS::

            sage: F = CombinatorialFreeModule(ZZ, [2,4,5])
            sage: G = CombinatorialFreeModule(ZZ, [2,4,7])
            sage: CP = cartesian_product([F, G])
            sage: CP._sets_keys()
            [0, 1]
        """
        return range(len(self._sets))

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
        return cartesian_product.symbol.join(["%s"%module for module in self._sets])
        # TODO: make this overridable by setting _name

    @cached_method
    def cartesian_embedding(self, i):
        """
        Return the natural embedding morphism of the ``i``-th
        cartesian factor (summand) of ``self`` into ``self``.

        INPUTS:

         - ``i`` -- an integer

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: phi = S.summand_embedding(0)
            sage: phi(F.monomial(4) + 2 * F.monomial(5))
            B[(0, 4)] + 2*B[(0, 5)]
            sage: phi(F.monomial(4) + 2 * F.monomial(6)).parent() == S
            True
            sage: phi(G.monomial(4)) # not implemented Should raise an error!  problem: G(F.monomial(4)) does not complain!!!!
        """
        assert i in self._sets_keys()
        return self._sets[i]._module_morphism(lambda t: self.monomial((i,t)), codomain = self)

    summand_embedding = cartesian_embedding

    @cached_method
    def cartesian_projection(self, i):
        """
        Return the natural projection onto the `i`-th cartesian factor
        (summand) of ``self``.

        INPUTS:

         - ``i`` -- an integer

        EXAMPLE::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: x = S.monomial((0,4)) + 2 * S.monomial((0,5)) + 3 * S.monomial((1,6))
            sage: S.cartesian_projection(0)(x)
            B[4] + 2*B[5]
            sage: S.cartesian_projection(1)(x)
            3*B[6]
            sage: S.cartesian_projection(0)(x).parent() == F
            True
            sage: S.cartesian_projection(1)(x).parent() == G
            True
        """
        assert i in self._sets_keys()
        module = self._sets[i]
        return self._module_morphism(lambda j_t: module.monomial(j_t[1]) if i == j_t[0] else module.zero(), codomain = module)

    summand_projection = cartesian_projection

    def _cartesian_product_of_elements(self, elements):
        """
        Return the cartesian product of the elements.

        INPUT:

        - ``elements`` -- a tuple with one element of each cartesian
          factor of ``self``

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: f =   F.monomial(4) + 2 * F.monomial(5)
            sage: g = 2*G.monomial(4) +     G.monomial(6)
            sage: S._cartesian_product_of_elements([f, g])
            B[(0, 4)] + 2*B[(0, 5)] + 2*B[(1, 4)] + B[(1, 6)]
            sage: S._cartesian_product_of_elements([f, g]).parent() == S
            True

        """
        return self.sum(self.summand_embedding(i)(elements[i]) for i in self._sets_keys())

    def cartesian_factors(self):
        """
        Return the factors of the cartesian product.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: S.cartesian_factors()
            (F, G)

        """
        return self._sets

    class Element(CombinatorialFreeModule.Element): # TODO: get rid of this inheritance
        pass

CombinatorialFreeModule.CartesianProduct = CombinatorialFreeModule_CartesianProduct
