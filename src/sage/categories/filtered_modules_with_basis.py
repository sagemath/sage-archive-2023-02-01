r"""
Filtered modules with basis

A *filtered module with basis* over a commutative ring `R`
means (for the purpose of this code) a filtered `R`-module
`M` with filtration `(F_0, F_1, F_2, \ldots)` endowed with a
basis `(b_i)_{i \in I}` of `M` and a partition of the set
`I` into subsets `I_0, I_1, I_2, \ldots` (which can be
empty) such that for every `n \in \NN`, the subfamily
`(b_i)_{i \in I_0 \cup I_1 \cup \cdots \cup I_n}` is a basis
of the `R`-submodule `F_n`.

For every `n \in \NN`, the `R`-submodule of `M` spanned by
`(b_i)_{i \in I_n}` is called the `*n*-th graded component*
of the filtered-module-with-basis `M`; the elements of
this submodule are referred to as *homogeneous elements of
degree `n`*. The `R`-module `M` is the direct sum of its
`n`-th graded components over all `n \in \NN`, and thus
becomes a graded `R`-module with basis. Conversely, any
graded `R`-module with basis canonically becomes a filtered
`R`-module with basis (by defining `F_n` as the direct sum
of the `0`-th, `1`-st, ..., `n`-th graded components, and
`I_n` as the indexing set of the basis of the `n`-th graded
component). Hence, the notion of a filtered `R`-module with
basis is equivalent to the notion of a graded `R`-module
with basis. However, the *category* of filtered `R`-modules
with basis is not the category of graded `R`-modules with
basis. Indeed, the *morphisms* of filtered `R`-modules with
basis are defined to be morphisms of `R`-modules which send
each `F_n` of the domain to the corresponding `F_n` of the
target; in contrast, the morphisms of graded `R`-modules
with basis must preserve each homogeneous component. Also,
the notion of a filtered algebra with basis differs from
that of a graded algebra with basis.
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.filtered_modules import FilteredModulesCategory

class FilteredModulesWithBasis(FilteredModulesCategory):
    """
    The category of filtered modules with a distinguished basis.

    A *filtered module with basis* over a commutative ring `R`
    means (for the purpose of this code) a filtered `R`-module
    `M` with filtration `(F_0, F_1, F_2, \ldots)` endowed with a
    basis `(b_i)_{i \in I}` of `M` and a partition of the set
    `I` into subsets `I_0, I_1, I_2, \ldots` (which can be
    empty) such that for every `n \in \NN`, the subfamily
    `(b_i)_{i \in I_0 \cup I_1 \cup \cdots \cup I_n}` is a basis
    of the `R`-submodule `F_n`.

    For every `n \in \NN`, the `R`-submodule of `M` spanned by
    `(b_i)_{i \in I_n}` is called the `*n*-th graded component*
    of the filtered-module-with-basis `M`; the elements of
    this submodule are referred to as *homogeneous elements of
    degree `n`*.

    EXAMPLES::

        sage: C = ModulesWithBasis(ZZ).Filtered(); C
        Category of filtered modules with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of filtered modules over Integer Ring,
         Category of modules with basis over Integer Ring]
        sage: C is ModulesWithBasis(ZZ).Filtered()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:

        # TODO: which syntax do we prefer?
        # A.basis(degree = 3)
        # A.basis().subset(degree=3)

        # This is related to the following design question:
        # If F = (f_i)_{i\in I} is a family, should ``F.subset(degree = 3)``
        # be the elements of F of degree 3 or those whose index is of degree 3?

        def basis(self, d=None):
            r"""
            Return the basis for (the ``d``-th homogeneous component
            of) ``self``.

            INPUT:

            - ``d`` -- (optional, default ``None``) nonnegative integer
              or ``None``

            If ``d`` is ``None``, returns a basis of the module.
            Otherwise, returns the basis of the homogeneous component
            of degree ``d`` (i.e., the subfamily of the basis of the
            whole module which consists only of the basis vectors
            lying in `F_d \setminus F_{d-1}`).

            EXAMPLES::

                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: A.basis(4)
                Lazy family (Term map from Partitions to An example of a filtered module with basis: the free module on partitions over Integer Ring(i))_{i in Partitions of the integer 4}

            Without arguments, the full basis is returned::

                sage: A.basis()
                Lazy family (Term map from Partitions to An example of a filtered module with basis: the free module on partitions over Integer Ring(i))_{i in Partitions}
                sage: A.basis()
                Lazy family (Term map from Partitions to An example of a filtered module with basis: the free module on partitions over Integer Ring(i))_{i in Partitions}
            """
            from sage.sets.family import Family
            if d is None:
                return Family(self._indices, self.monomial)
            else:
                return Family(self._indices.subset(size=d), self.monomial)

    class ElementMethods:

        def is_homogeneous(self):
            r"""
            Return whether ``self`` is homogeneous.

            EXAMPLES::

                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x=A(Partition((3,2,1)))
                sage: y=A(Partition((4,4,1)))
                sage: z=A(Partition((2,2,2)))
                sage: (3*x).is_homogeneous()
                True
                sage: (x - y).is_homogeneous()
                False
                sage: (x+2*z).is_homogeneous()
                True
            """
            degree_on_basis = self.parent().degree_on_basis
            degree = None
            for m in self.support():
                if degree is None:
                    degree = degree_on_basis(m)
                else:
                    if degree != degree_on_basis(m):
                        return False
            return True

        def degree(self):
            r"""
            The degree of a nonzero homogeneous element ``self`` in the
            filtered module.

            .. NOTE::

                This raises an error if the element is not homogeneous.
                Another implementation option would be to return the
                maximum of the degrees of the homogeneous summands.

            EXAMPLES::

                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A(Partition((3,2,1)))
                sage: y = A(Partition((4,4,1)))
                sage: z = A(Partition((2,2,2)))
                sage: x.degree()
                6
                sage: (x + 2*z).degree()
                6
                sage: (y - x).degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous
            """
            if not self.support():
                raise ValueError("the zero element does not have a well-defined degree")
            if not self.is_homogeneous():
                raise ValueError("element is not homogeneous")
            return self.parent().degree_on_basis(self.leading_support())

        def homogeneous_component(self, n):
            """
            Return the homogeneous component of degree ``n`` of this
            element.

            Let `m` be an element of a filtered `R`-module `M` with
            basis. Then, `m` can be uniquely written in the form
            `m = m_0 + m_1 + m_2 + \ldots`, where each `m_i` is a
            homogeneous element of degree `i`. For `n \in \NN`, we
            define the homogeneous component of degree `n` of the
            element `m` to be `m_n`.

            EXAMPLES::

                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A.an_element(); x
                2*P[] + 2*P[1] + 3*P[2]
                sage: x.homogeneous_component(-1)
                0
                sage: x.homogeneous_component(0)
                2*P[]
                sage: x.homogeneous_component(1)
                2*P[1]
                sage: x.homogeneous_component(2)
                3*P[2]
                sage: x.homogeneous_component(3)
                0

            TESTS:

            Check that this really returns ``A.zero()`` and not a plain ``0``::

                sage: x.homogeneous_component(3).parent() is A
                True
            """
            degree_on_basis = self.parent().degree_on_basis
            return self.parent().sum_of_terms((i, c)
                                              for (i, c) in self
                                              if degree_on_basis(i) == n)

        def truncate(self, n):
            """
            Return the sum of the homogeneous components of degree
            strictly less than ``n`` of ``self``.

            See :meth:`homogeneous_component` for the notion of a
            homogeneous component.

            EXAMPLES::

                sage: A = ModulesWithBasis(ZZ).Filtered().example()
                sage: x = A.an_element(); x
                2*P[] + 2*P[1] + 3*P[2]
                sage: x.truncate(0)
                0
                sage: x.truncate(1)
                2*P[]
                sage: x.truncate(2)
                2*P[] + 2*P[1]
                sage: x.truncate(3)
                2*P[] + 2*P[1] + 3*P[2]

            TESTS:

            Check that this really return ``A.zero()`` and not a plain ``0``::

                sage: x.truncate(0).parent() is A
                True
            """
            degree_on_basis = self.parent().degree_on_basis
            return self.parent().sum_of_terms((i, c) for (i, c) in self
                                              if degree_on_basis(i) < n)

