r"""
Finite Weyl Groups
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.misc.cachefunc import cached_method, cached_in_parent_method

class FiniteWeylGroups(CategoryWithAxiom):
    """
    The category of finite Weyl groups.

    EXAMPLES::

        sage: C = FiniteWeylGroups()
        sage: C
        Category of finite weyl groups
        sage: C.super_categories()
        [Category of finite coxeter groups, Category of weyl groups]
        sage: C.example()
        The symmetric group on {0, ..., 3}

    TESTS::

        sage: W = FiniteWeylGroups().example()
        sage: TestSuite(W).run(verbose = "True")
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_has_descent() . . . pass
        running ._test_inverse() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_reduced_word() . . . pass
        running ._test_simple_projections() . . . pass
        running ._test_some_elements() . . . pass
    """

    class ParentMethods:

        @cached_method
        def m_cambrian_lattice(W,c,m):
            """
            INPUT:

            - ``c`` -- a Coxeter element of ``self`` (as a tuple, or as an element of ``self``)
            - ``m`` -- a positive integer (default: 1)

            Return the m-Cambrian lattice on ``m``-delta sequences, realized on roots rather than reflections (see arXiv:1503.00710 and arXiv:math/0611106).
            ``m``-delta sequences are certain ``m``-colored minimal factorizations of ``c`` into reflections.

            EXAMPLES::

                sage: WeylGroup(["A",2]).m_cambrian_lattice((1,2),1)
                Finite lattice containing 5 elements

                sage: WeylGroup(["A",2]).m_cambrian_lattice((1,2),2)
                Finite lattice containing 12 elements

            """
            from sage.combinat.posets.posets import Poset
            from sage.combinat.posets.lattices import LatticePoset
            if hasattr(c,"reduced_word"):
               c = c.reduced_word()
            elif not isinstance(c,list):
               c = list(c)
            inv_woc = [t.reflection_to_root() for t in W.inversion_sequence(W.long_element().coxeter_sorting_word(c))]
            S = [s.reflection_to_root() for s in W.simple_reflections()]
            PhiP = [t.reflection_to_root() for t in W.reflections().keys()]
            id = sorted([[t,0] for t in S])
            elements = []
            covers = []
            new = [id]
            while new != []:
                for new_element in new:
                    new.remove(new_element)
                    elements.append(new_element)
                    for t in new_element:
                        if t[1]<m:
                            cov_element = [s for s in new_element if s!=t]
                            cov_element.append([t[0],t[1]+1])
                            for t_conj in [[i,t[1]] for i in inv_woc[inv_woc.index(t[0]):]]+[[i,t[1]+1] for i in inv_woc[:inv_woc.index(t[0])]]:
                                if t_conj in cov_element:
                                    cov_element.remove(t_conj)
                                    tmp = t_conj[0].weyl_action(t[0].associated_reflection())
                                    if tmp in PhiP:
                                        cov_element.append([tmp,t_conj[1]])
                                    else:
                                        cov_element.append([-tmp,t_conj[1]-1])
                            cov_element = sorted(cov_element)
                            if cov_element not in elements and cov_element not in new:
                                new.append(cov_element)
                            covers.append([tuple(map(tuple,new_element)),tuple(map(tuple,cov_element))])
            return LatticePoset([[tuple(map(tuple,e)) for e in elements],covers],cover_relations=True)


    class ElementMethods:
        pass
