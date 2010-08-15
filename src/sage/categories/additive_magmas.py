r"""
Additive Magmas
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.sets_cat import Sets
from sage.structure.sage_object import have_same_parent

class AdditiveMagmas(Category):
    """
    The category of additive magmas, i.e. sets with an binary
    operation ``+``.

    EXAMPLES::

        sage: AdditiveMagmas()
        Category of additive magmas
        sage: AdditiveMagmas().super_categories()
        [Category of sets]
        sage: AdditiveMagmas().all_super_categories()
        [Category of additive magmas, Category of sets, Category of sets with partial maps, Category of objects]

    TESTS::

        sage: C = AdditiveMagmas()
        sage: TestSuite(C).run()

    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: AdditiveMagmas().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class ParentMethods:

        def summation(self, x, y):
            """
            The binary addition operator of the semigroup

            INPUT:

             - ``x``, ``y`` -- elements of this additive semigroup

            Returns the sum of ``x`` and ``y``

            EXAMPLES::

                sage: S = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: S.summation(a, b)
                a + b

            A parent in ``AdditiveMagmas()`` must
            either implement :meth:`.summation` in the parent class or
            ``_add_`` in the element class. By default, the addition
            method on elements ``x._add_(y)`` calls
            ``S.summation(x,y)``, and reciprocally.


            As a bonus effect, ``S.summation`` by itself models the
            binary function from ``S`` to ``S``::

                sage: bin = S.summation
                sage: bin(a,b)
                a + b

            Here, ``S.summation`` is just a bound method. Whenever
            possible, it is recommended to enrich ``S.summation`` with
            extra mathematical structure. Lazy attributes can come
            handy for this.

            Todo: add an example.
            """
            return x._add_(y)

        summation_from_element_class_add = summation

        def __init_extra__(self):
            """
            TESTS::

                sage: S = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: a + b # indirect doctest
                a + b
                sage: a.__class__._add_ == a.__class__._add_parent
                True
            """
            # This should instead register the summation to the coercion model
            # But this is not yet implemented in the coercion model
            if (self.summation != self.summation_from_element_class_add) and hasattr(self, "element_class") and hasattr(self.element_class, "_add_parent"):
                self.element_class._add_ = self.element_class._add_parent


        def addition_table(self, names='letters', elements=None):
            r"""
            Returns a table describing the addition operation.

            .. note:: The order of the elements in the row and column
              headings is equal to the order given by the table's
              :meth:`~sage.matrix.operation_table.OperationTable.list`
              method.  The association can also be retrieved with the
              :meth:`~sage.matrix.operation_table.OperationTable.dict`
              method.

            INPUTS:

            - ``names`` - the type of names used

              * ``letters`` - lowercase ASCII letters are used
                for a base 26 representation of the elements'
                positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading 'a's.
              * ``digits`` - base 10 representation of the
                elements' positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading zeros.
              * ``elements`` - the string representations
                of the elements themselves.
              * a list - a list of strings, where the length
                of the list equals the number of elements.
            - ``elements`` - default = ``None``.  A list of
              elements of the set.  This may be used to impose an
              alternate ordering on the elements, perhaps
              when this is used in the context of a particular structure.
              The default is to use whatever ordering the
              ``S.list``
              method returns.  Or the ``elements`` can be a subset
              which is closed under the operation. In particular,
              this can be used when the base set is infinite.

            OUTPUT:
            The addition table as an object of the class
            :class:`~sage.matrix.operation_table.OperationTable`
            which defines several methods for manipulating and
            displaying the table.  See the documentation there
            for full details to supplement the documentation
            here.

            EXAMPLES:

            All that is required is that an algebraic structure
            has an addition defined.The default is to represent
            elements as lowercase ASCII letters.  ::

                sage: R=IntegerModRing(5)
                sage: R.addition_table()
                +  a b c d e
                 +----------
                a| a b c d e
                b| b c d e a
                c| c d e a b
                d| d e a b c
                e| e a b c d

            The ``names`` argument allows displaying the elements in different ways.  Requesting ``elements`` will use the representation of the elements of the set.  Requesting ``digits`` will include leading zeros as padding.  ::

                sage: R=IntegerModRing(11)
                sage: P=R.addition_table(names='elements')
                sage: P
                 +   0  1  2  3  4  5  6  7  8  9 10
                  +---------------------------------
                 0|  0  1  2  3  4  5  6  7  8  9 10
                 1|  1  2  3  4  5  6  7  8  9 10  0
                 2|  2  3  4  5  6  7  8  9 10  0  1
                 3|  3  4  5  6  7  8  9 10  0  1  2
                 4|  4  5  6  7  8  9 10  0  1  2  3
                 5|  5  6  7  8  9 10  0  1  2  3  4
                 6|  6  7  8  9 10  0  1  2  3  4  5
                 7|  7  8  9 10  0  1  2  3  4  5  6
                 8|  8  9 10  0  1  2  3  4  5  6  7
                 9|  9 10  0  1  2  3  4  5  6  7  8
                10| 10  0  1  2  3  4  5  6  7  8  9

                sage: T=R.addition_table(names='digits')
                sage: T
                +  00 01 02 03 04 05 06 07 08 09 10
                  +---------------------------------
                00| 00 01 02 03 04 05 06 07 08 09 10
                01| 01 02 03 04 05 06 07 08 09 10 00
                02| 02 03 04 05 06 07 08 09 10 00 01
                03| 03 04 05 06 07 08 09 10 00 01 02
                04| 04 05 06 07 08 09 10 00 01 02 03
                05| 05 06 07 08 09 10 00 01 02 03 04
                06| 06 07 08 09 10 00 01 02 03 04 05
                07| 07 08 09 10 00 01 02 03 04 05 06
                08| 08 09 10 00 01 02 03 04 05 06 07
                09| 09 10 00 01 02 03 04 05 06 07 08
                10| 10 00 01 02 03 04 05 06 07 08 09

            Specifying the elements in an alternative order can provide
            more insight into how the operation behaves.  ::

                sage: S=IntegerModRing(7)
                sage: elts = [0, 3, 6, 2, 5, 1, 4]
                sage: S.addition_table(elements=elts)
                +  a b c d e f g
                 +--------------
                a| a b c d e f g
                b| b c d e f g a
                c| c d e f g a b
                d| d e f g a b c
                e| e f g a b c d
                f| f g a b c d e
                g| g a b c d e f

            The ``elements`` argument can be used to provide
            a subset of the elements of the structure.  The subset
            must be closed under the operation.  Elements need only
            be in a form that can be coerced into the set.  The
            ``names`` argument can also be used to request that
            the elements be represented with their usual string
            representation.  ::

                sage: T=IntegerModRing(12)
                sage: elts=[0, 3, 6, 9]
                sage: T.addition_table(names='elements', elements=elts)
                +  0 3 6 9
                 +--------
                0| 0 3 6 9
                3| 3 6 9 0
                6| 6 9 0 3
                9| 9 0 3 6

            The table returned can be manipulated in various ways.  See
            the documentation for
            :class:`~sage.matrix.operation_table.OperationTable` for more
            comprehensive documentation. ::

                sage: R=IntegerModRing(3)
                sage: T=R.addition_table()
                sage: T.column_keys()
                (0, 1, 2)
                sage: sorted(T.translation().items())
                [('a', 0), ('b', 1), ('c', 2)]
                sage: T.change_names(['x', 'y', 'z'])
                sage: sorted(T.translation().items())
                [('x', 0), ('y', 1), ('z', 2)]
                sage: T
                +  x y z
                 +------
                x| x y z
                y| y z x
                z| z x y
            """
            from sage.matrix.operation_table import OperationTable
            import operator
            return OperationTable(self, operation=operator.add, names=names, elements=elements)

    class ElementMethods:

        # This could eventually be moved to SageObject
        def __add__(self, right):
            r"""
            Sum of two elements

            This calls the `_add_` method of ``self``, if it is
            available and the two elements have the same parent.

            Otherwise, the job is delegated to the coercion model.

            Do not override; instead implement an ``_add_`` method in the
            element class or a ``summation`` method in the parent class.

            EXAMPLES::

                sage: F = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = F.additive_semigroup_generators()
                sage: a + b
                a + b
            """
            if have_same_parent(self, right) and hasattr(self, "_add_"):
                return self._add_(right)
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(self, right, operator.add)

        def __radd__(self, left):
            r"""
            Handles the sum of two elements, when the left hand side
            needs to be coerced first.

            EXAMPLES::

                sage: F = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = F.additive_semigroup_generators()
                sage: a.__radd__(b)
                a + b
            """
            if have_same_parent(left, self) and hasattr(left, "_add_"):
                return left._add_(self)
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(left, self, operator.add)

        @abstract_method(optional = True)
        def _add_(self, right):
            """
            Sum of two elements

            INPUT:

             - ``self``, ``right`` -- two elements with the same parent

            OUTPUT:

             - an element of the same parent

            EXAMPLES::

                sage: F = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = F.additive_semigroup_generators()
                sage: a._add_(b)
                a + b
            """

        def _add_parent(self, other):
            r"""
            Returns the sum of the two elements, calculated using
            the ``summation`` method of the parent.

            This is the default implementation of _add_ if
            ``summation`` is implemented in the parent.

            INPUT:

             - ``other`` -- an element of the parent of ``self``

            OUTPUT:

            an element of the parent of ``self``

            EXAMPLES::

                sage: S = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: a._add_parent(b)
                a + b
            """
            return self.parent().summation(self, other)
