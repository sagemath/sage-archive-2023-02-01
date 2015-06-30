r"""
Bialgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.all import Algebras, Coalgebras

class Bialgebras(Category_over_base_ring):
    """
    The category of bialgebras

    EXAMPLES::

        sage: Bialgebras(ZZ)
        Category of bialgebras over Integer Ring
        sage: Bialgebras(ZZ).super_categories()
        [Category of algebras over Integer Ring, Category of coalgebras over Integer Ring]

    TESTS::

        sage: TestSuite(Bialgebras(ZZ)).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: Bialgebras(QQ).super_categories()
            [Category of algebras over Rational Field, Category of coalgebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R), Coalgebras(R)]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of bialgebras defines no additional
        structure: a morphism of coalgebras and of algebras between
        two bialgebras is a bialgebra morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: This category should be a :class:`CategoryWithAxiom`.

        EXAMPLES::

            sage: Bialgebras(QQ).additional_structure()
        """
        return None

    class ParentMethods:
        pass

    class ElementMethods:

        def convolution_product(h,a,b):
            r"""
            input: h - an element of a Hopf algebra H
                   a,b - linear maps from H to H
            output: [a*b](h)
            """
            H = h.parent()
            out = 0
            for (bimonom,coef) in h.coproduct():
                out += coef*a(H(bimonom[0]))*b(H(bimonom[1]))
            return out

        def convolution_power(h, L, n):
            r"""
            input: h - an element of a Hopf algebra H
                   L - linear map from H to H
                   n - the convolution power to which to take 'L'
            output: [L^*n](h)
            """
            from sage.categories.tensor import tensor
            H = h.parent()
            def n_fold_coproduct(h, n):
                H = h.parent()
                if n == 0:
                    return H(h.counit())
                elif n == 1:
                    return h
                elif n == 2:
                    return h.coproduct()
                else:
                    # apply some kind of multilinear recursion
                    Hn = tensor([H]*n) # or: tensor([H for i in range(n)])
                    terms = []
                    hh = n_fold_coproduct(h, n-1)
                    for (monom,cof) in hh:
                        h0 = H(monom[0]).coproduct()
                        terms += [(tuple((h00, h01) + monom[1:]), cof0 * cof) for ((h00, h01), cof0) in h0]
                    return Hn.sum_of_terms(terms)
            hhh = n_fold_coproduct(h,n)
            out = H.zero()
            for term in hhh:
                out += H.prod(L(H(t)) for t in term[0]) * term[1]
            return out

        def hopf_power(h,n=2):
            r"""
            Input:
                h - an element of a Hopf algebra H
                n - the convolution power of the identity morphism to use
            Output:
                the nth convolution power of the identity morphism, applied to h., i.e., [id^*n](h)

            Remark: for historical reasons (see saga of Frobenius-Schur indicators), the second power deserves special attention.
            we use '2' as the default value for 'n'
            """
            H = h.parent()
            def Id(x):
                return x
            def S(x):
                return x.antipode()

            if n<0:
                L = S
            else:
                L = Id

            if n==0:
                return H(h.counit())
            elif abs(n)==1:
                return L(h)
            else:
                return h.convolution_power(L,abs(n))

