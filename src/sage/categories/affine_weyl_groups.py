r"""
Affine Weyl Groups
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.weyl_groups import WeylGroups

class AffineWeylGroups(Category):

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: AffineWeylGroups().super_categories()
            [Category of weyl groups]
        """
        return [WeylGroups()]

    class ParentMethods:

        @cached_method
        def special_node(self):
            """
            Returns the distinguished special node of the underlying
            Dynkin diagram

            EXAMPLES::

                sage: W=WeylGroup(['A',3,1])
                sage: W.special_node()
                0

            """
            return self.cartan_type().special_node()

        def affine_grassmannian_elements_of_given_length(self, k):
            """
            Returns the affine Grassmannian elements of length `k`, as a list.

            EXAMPLES::

                sage: W=WeylGroup(['A',3,1])
                sage: [x.reduced_word() for x in W.affine_grassmannian_elements_of_given_length(3)]
                [[2, 1, 0], [3, 1, 0], [2, 3, 0]]

            SEE ALSO: :meth:`AffineWeylGroups.ElementMethods.is_affine_grassmannian`.

            TODO: should return an enumerated set, with iterator, ...
            """
            if k == 0:
                return [self.unit()]
            w = []
            s = self.simple_reflections()
            for x in self.affine_grassmannian_elements_of_given_length(k-1):
                for i in x.descents(side="left", positive = True):
                    y = s[i]*x
                    if y not in w and y.is_affine_grassmannian():
                        w.append(y)
            return w

    class ElementMethods:
        def is_affine_grassmannian(self):
            """
            INPUT:
             - self: an element of an affine weyl group

            Tests whether self is Grassmannian, i.e. any of the following
            equivalent properties hold:
             - all reduced words for self end with 0.
             - self is the identity, or 0 is its single right descent.
             - self is a mimimal coset representative for W / cl W.

            EXAMPLES::

                sage: W=WeylGroup(['A',3,1])
                sage: w=W.from_reduced_word([2,1,0])
                sage: w.is_affine_grassmannian()
                True
                sage: w=W.from_reduced_word([2,0])
                sage: w.is_affine_grassmannian()
                False
                sage: W.unit().is_affine_grassmannian()
                True
            """
            return self.descents() in [[], [self.parent().special_node()]]
