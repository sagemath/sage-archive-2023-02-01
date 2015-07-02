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
from sage.categories.tensor import tensor
from sage.functions.other import floor, ceil
from sage.rings.integer import Integer
from sage.categories.modules_with_basis import ModulesWithBasis

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

        def convolution_product(self, *maplist):
            """
            Given a maplist `(R, S, ..., T)` of length `n`, return the new map representing their convolution product.

            MATH::

                (R*S*\cdots *T)(h) := \mu^{(n-1)} \circ (R \otimes S \otimes\cdot\otimes T) \circ \Delta^{(n-1)}(h)

            .. SEEALSO::

                :meth:`sage.categories.bialgebras.ElementMethods.convolution_product`

            .. TODO::

                Remove dependency on modules_with_basis methods.

            TESTS::

                sage: Id = lambda x: x
                sage: Antipode = lambda x: x.antipode()
                sage: Proj2 = lambda x: x.parent().sum_of_terms([(m,c) for (m,c) in x if m.size()==2])

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: T = R.convolution_product([Id,Id])
                sage: [T(R(comp)) for comp in Compositions(3)]
                [4*R[1, 1, 1] + R[1, 2] + R[2, 1],
                 2*R[1, 1, 1] + 4*R[1, 2] + 2*R[2, 1] + 2*R[3],
                 2*R[1, 1, 1] + 2*R[1, 2] + 4*R[2, 1] + 2*R[3],
                 R[1, 2] + R[2, 1] + 4*R[3]]
                sage: T = R.convolution_product(Proj2,Id)
                sage: [T(R([i])) for i in range(1,5)]
                [0, R[2], R[2, 1] + R[3], R[2, 2] + R[4]]

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: T = m.convolution_product()
                sage: [T(m(lam)) for lam in SetPartitions(0).list()+SetPartitions(2).list()]
                [m{}, 0, 0]
                sage: T = m.convolution_product(Proj2,Id)
                sage: [T(m(lam)) for lam in SetPartitions(3)]
                [0,
                 m{{1, 2}, {3}} + m{{1, 2, 3}},
                 m{{1, 2}, {3}} + m{{1, 2, 3}},
                 m{{1, 2}, {3}} + m{{1, 2, 3}},
                 3*m{{1}, {2}, {3}} + 3*m{{1}, {2, 3}} + 3*m{{1, 3}, {2}}]

                sage: G = SymmetricGroup(3); QG = GroupAlgebra(G,QQ)
                sage: x = QG.sum_of_terms([(p,p.number_of_peaks() + p.number_of_inversions()) for p in Permutations(3)]); x
                2*[1, 3, 2] + [2, 1, 3] + 3*[2, 3, 1] + 2*[3, 1, 2] + 3*[3, 2, 1]
                sage: T = QG.convolution_product(Antipode,Antipode,Id)
                sage: T(x)
                2*[1, 3, 2] + [2, 1, 3] + 2*[2, 3, 1] + 3*[3, 1, 2] + 3*[3, 2, 1]
            """
            # Be flexible on how the maps are entered...
            if len(maplist)==1 and isinstance(maplist[0], (list,tuple)):
                T = tuple(maplist[0])
            else:
                T = maplist

            n = len(T)
            H = self

            # I do not know how to keep .convolution_product() from showing up in the
            # list of methods available to, e.g., SymmetricFunctions(QQ).
            # At present, the code below assumes Parent is something more like SymmetricFunctions(QQ).m()
            # The code should either be rewritten (I don't know how),
            # or we should do a check, e.g.,
            if not H in ModulesWithBasis:
                raise TypeError('`self` must belong to ModulesWithBasis. Try a basis of %s' % H._repr_())

            if n == 0:
                return H.module_morphism(on_basis=lambda x: H.one() * H.term(x).counit(), codomain = H)
            elif n == 1:
                return H.module_morphism(on_basis=lambda x: T[0](H.term(x)), codomain = H)
            else:
                HHH = tensor((H,)*n)
                Delta_n = H.module_morphism(on_basis=lambda x: H.term(x).coproduct_iterated(n-1), codomain = HHH)
                apply_T = HHH.module_morphism(on_basis=lambda args: tensor([T[i](H.term(args[i])) for i in range(n)]), codomain = HHH)
                Mu_n = HHH.module_morphism(on_basis=lambda args: H.prod([H.term(h) for h in args]), codomain = H)

                return Mu_n * apply_T * Delta_n

    class ElementMethods:

        def adams_operator(self, n):
            """
            Compute the n-th convolution power of the identity morphism `Id` on self.

            INPUT:

            - ``n`` -- a nonnegative integer.

            OUTPUT:

            - the element of self.parent() corresponding to `Id^{*n}(self)`.

            .. SEEALSO::

                :mod:`sage.categories.hopf_algebras.ElementMethods.convolution_product`
                :mod:`sage.categories.hopf_algebras.ElementMethods.convolution_product`

                (In the literature, this is also called a Hopf power or Sweedler power, cf. [AL2015]_.)

            REFERENCES:

            .. [AL2015] The characteristic polynomial of the Adams operators on graded connected Hopf algebras.
                Marcelo Aguiar and Aaron Lauve.
                Algebra Number Theory, v.9, 2015, n.3, 2015.

            .. TODO::

                Move to hopf_algebras.py (i.e., remove dependency on modules_with_basis methods).

            TESTS::

                sage: h = SymmetricFunctions(QQ).h()
                sage: h[5].adams_operator(2)
                2*h[3, 2] + 2*h[4, 1] + 2*h[5]
                sage: h[5].plethysm(2*h[1])
                2*h[3, 2] + 2*h[4, 1] + 2*h[5]
                sage: h([]).adams_operator(0), h([]).adams_operator(1)
                (h[], h[])
                sage: h[3,2].adams_operator(0), h[3,2].adams_operator(1)
                (0, h[3, 2])

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S[4].adams_operator(5)
                5*S[1, 1, 1, 1] + 10*S[1, 1, 2] + 10*S[1, 2, 1] + 10*S[1, 3] + 10*S[2, 1, 1] + 10*S[2, 2] + 10*S[3, 1] + 5*S[4]

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m[[1,3],[2]].adams_operator(-2)
                3*m{{1}, {2, 3}} + 3*m{{1, 2}, {3}} + 6*m{{1, 2, 3}} - 2*m{{1, 3}, {2}}

            """
            if n < 0:
                raise ValueError("cannot take less than 0 coproduct iterations: %s < 0" % str(n))
            return self.convolution_product([lambda x: x] * n)

        def convolution_product(self, *maplist):
            """
            Given a maplist `(R, S, ..., T)` of length `n`, compute the action of their convolution product on ``self.``

            MATH::

                (R*S*\cdots *T)(h) := \mu^{(n-1)} \circ (R \otimes S \otimes\cdot\otimes T) \circ \Delta^{(n-1)}(h)

            where `\Delta^{(k)} := \bigl(\Delta \otimes \mathrm{Id}^{\otimes(k-1)}\bigr) \circ \Delta^{(k-1)}`,
            with `\Delta^{(1)} = \Delta` (the ordinary coproduct) and `\Delta^{(0)} = \mathrm{Id}`;
            and with `\mu^{(k)} := \mu \circ \bigl(\mu^{(k-1)} \otimes \mathrm{Id})` and `\mu^{(1)} = \mu`
            (the ordinary product). See [Sw1969]_.

            (In the literature, one finds, e.g., `\Delta^{(2)}` for what we denote above as `\Delta^{(1)}`. See [KMN2012]_.)

            REFERENCES:

            .. [KMN2012] On the trace of the antipode and higher indicators.
                Yevgenia Kashina and Susan Montgomery and Richard Ng.
                Israel J. Math., v.188, 2012.

            .. [Sw1969] Hopf algebras.
                Moss Sweedler.
                W.A. Benjamin, Math Lec Note Ser., 1969.

            .. TODO::

                Remove dependency on modules_with_basis methods.

            EXAMPLES::

                sage: Id = lambda x: x
                sage: Antipode = lambda x: x.antipode()
                sage: s = SymmetricFunctions(QQ).schur()
                sage: s[3].convolution_product(Id,Id)
                2*s[2, 1] + 4*s[3]
                sage: s[3,2].convolution_product(Id,Id) == s[3,2].convolution_product([Id,Id])
                True
                sage: s[3,2].convolution_product(Id) == s[3,2]
                True
                sage: s[3,2].convolution_product() == s[3,2].convolution_product(Antipode,Id) == s[3,2].convolution_product(Id,Antipode)
                True
                sage: s[3,2].counit().parent() == s[3,2].convolution_product().parent()
                False

                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi[2,1].convolution_product(Id,Id,Id)
                3*Psi[1, 2] + 6*Psi[2, 1]
                sage: (Psi[5,1] - Psi[1,5]).convolution_product(Id,Id,Id)
                -3*Psi[1, 5] + 3*Psi[5, 1]

                sage: G = SymmetricGroup(3); QG = GroupAlgebra(G,QQ)
                sage: x = QG.sum_of_terms([(p,p.length()) for p in Permutations(3)]); x
                [1, 3, 2] + [2, 1, 3] + 2*[2, 3, 1] + 2*[3, 1, 2] + 3*[3, 2, 1]
                sage: x.convolution_product(Id,Id)
                5*[1, 2, 3] + 2*[2, 3, 1] + 2*[3, 1, 2]
                sage: x.convolution_product(Id,Id,Id)
                4*[1, 2, 3] + [1, 3, 2] + [2, 1, 3] + 3*[3, 2, 1]
                sage: x.convolution_product(Id,Id,Id,Id,Id,Id)
                9*[1, 2, 3]


            TESTS::

                sage: Id = lambda x: x
                sage: Antipode = lambda x: x.antipode()

                sage: h = SymmetricFunctions(QQ).h()
                sage: h[5].convolution_product([Id,Id])
                2*h[3, 2] + 2*h[4, 1] + 2*h[5]
                sage: h.one().convolution_product([Id,Antipode])
                h[]
                sage: h[3,2].convolution_product([Id,Antipode])
                0
                sage: h.one().convolution_product([Id,Antipode]) == h.one().convolution_product()
                True

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S[4].convolution_product([Id]*5)
                5*S[1, 1, 1, 1] + 10*S[1, 1, 2] + 10*S[1, 2, 1] + 10*S[1, 3] + 10*S[2, 1, 1] + 10*S[2, 2] + 10*S[3, 1] + 5*S[4]

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m[[1,3],[2]].convolution_product([Antipode,Antipode])
                3*m{{1}, {2, 3}} + 3*m{{1, 2}, {3}} + 6*m{{1, 2, 3}} - 2*m{{1, 3}, {2}}

                sage: m[[]].convolution_product([]), m[[1,3],[2]].convolution_product([])
                (m{}, 0)

                sage: QS = SymmetricGroupAlgebra(QQ,5)
                sage: x = QS.sum_of_terms(zip(Permutations(5)[3:6],[1,2,3])); x
                [1, 2, 4, 5, 3] + 2*[1, 2, 5, 3, 4] + 3*[1, 2, 5, 4, 3]
                sage: x.convolution_product([Antipode,Id])
                6*[1, 2, 3, 4, 5]
                sage: x.convolution_product(Id, Antipode, Antipode, Antipode)
                3*[1, 2, 3, 4, 5] + [1, 2, 4, 5, 3] + 2*[1, 2, 5, 3, 4]

                sage: G = SymmetricGroup(3); QG = GroupAlgebra(G,QQ)
                sage: x = QG.sum_of_terms([(p,p.length()) for p in Permutations(3)]); x
                [1, 3, 2] + [2, 1, 3] + 2*[2, 3, 1] + 2*[3, 1, 2] + 3*[3, 2, 1]
                sage: x.convolution_product(Antipode,Id)
                9*[1, 2, 3]
                sage: x.convolution_product([Id, Antipode, Antipode, Antipode])
                5*[1, 2, 3] + 2*[2, 3, 1] + 2*[3, 1, 2]
            """
            # Be flexible on how the maps are entered...
            if len(maplist)==1 and isinstance(maplist[0], (list,tuple)):
                T = tuple(maplist[0])
            else:
                T = maplist

            H = self.parent()
            n = len(T)
            if n == 0:
                return H.one() * self.counit()
            elif n == 1:
                return T[0](self)
            else:
                # We apply the maps T_i and products concurrently with coproducts, as this
                # seems to be faster than applying a composition of maps, e.g., (H.nfold_product) * tensor(T) * (H.nfold_coproduct).

                i = 0
                out = tensor((H.one(),self))
                HH = tensor((H,H))

                #ALGORITHM:
                #`split_convolve` moves terms of the form x # y to x*Ti(y1) # y2 in Sweedler notation.
                split_convolve = lambda (x,y): (((xy1,y2),c*d) for ((y1,y2),d) in H.term(y).coproduct() for (xy1,c) in H.term(x)*T[i](H.term(y1)))
                while i < n-1:
                    out = HH.module_morphism(on_basis=lambda t: HH.sum_of_terms(split_convolve(t)), codomain = HH)(out)
                    i += 1

                #Apply final map `T_n` to last term, `y`, and multiply.
                out = HH.module_morphism(on_basis=lambda (x,y): H.term(x)*T[n-1](H.term(y)), codomain=H)(out)

#                 #ALGORITHM:

                return out
                # ------------
                # IMPLEMENTATION NOTE:
                # In the `module_morphism()` version of this code (copied below), Sage sometimes throws a `TypeError`. E.g.,
                # -------
                # sage: Antipode = lambda x: x.antipode()
                #
                # sage: QS = SymmetricGroupAlgebra(QQ,3)
                # sage: x = QS.sum_of_terms([(p,p.length()) for p in Permutations(3)]); x
                # [1, 3, 2] + [2, 1, 3] + 2*[2, 3, 1] + 2*[3, 1, 2] + 3*[3, 2, 1]
                # sage: x.convolution_product([Antipode, Antipode])
                # 5*[1, 2, 3] + 2*[2, 3, 1] + 2*[3, 1, 2]
                #
                # sage: QG GroupAlgebra(SymmetricGroup(3),QQ)
                # sage: x = QG.sum_of_terms([(p,p.length()) for p in Permutations(3)]); x
                # [1, 3, 2] + [2, 1, 3] + 2*[2, 3, 1] + 2*[3, 1, 2] + 3*[3, 2, 1]
                # sage: x.convolution_product([Antipode,Antipode])
                # TypeError: Don't know how to create an element of Group algebra of group...
                # -------
                # #`split_convolve` moves terms of the form x # y to x*Ti(y1) # y2 in Sweedler notation.
                # split_convolve = lambda x,y: (x.tensor(y.coproduct())).apply_multilinear_morphism(convolve, codomain = HH)
                # convolve = lambda x,y1,y2: tensor([x * T[i](y1), y2])
                #
                # while i < n-1:
                #     out = out.apply_multilinear_morphism(split_convolve, codomain = HH)
                #     i += 1
                #
                # #Apply final map `T_n` to last term, `y`, and multiply.
                # out = out.apply_multilinear_morphism(lambda x,y: x * T[n-1](y), codomain = H)
                # -------
                # #`split_convolve` moves terms of the form x # y to x*Ti(y1) # y2 in Sweedler notation.
                # split_convolve = lambda (x,y): (((xy1,y2),c*d) for ((y1,y2),d) in H(y).coproduct() for (xy1,c) in H(x)*T[i](H(y1)))
                # while i < n-1:
                #     out = HH.module_morphism(on_basis=lambda t: HH.sum_of_terms(split_convolve(t)), codomain = HH)(out)
                #     i += 1

                # #Apply final map `T_n` to last term, `y`, and multiply.
                # out = HH.module_morphism(on_basis=lambda (x,y): H(x)*T[n-1](H(y)), codomain=H)(out)
                # ------------



        def coproduct_iterated(self, n=1):
            r"""
            Apply `k-1` coproducts to ``self``.

            EXAMPLES::

                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi[2,2].coproduct_iterated(0)
                Psi[2, 2]
                sage: Psi[2,2].coproduct_iterated(2)
                Psi[] # Psi[] # Psi[2, 2] + 2*Psi[] # Psi[2] # Psi[2] + Psi[] # Psi[2, 2] # Psi[] + 2*Psi[2] # Psi[] # Psi[2] + 2*Psi[2] # Psi[2] # Psi[] + Psi[2, 2] # Psi[] # Psi[]

            TESTS::

                sage: p = SymmetricFunctions(QQ).p()
                sage: p[5,2,2].coproduct_iterated()
                p[] # p[5, 2, 2] + 2*p[2] # p[5, 2] + p[2, 2] # p[5] + p[5] # p[2, 2] + 2*p[5, 2] # p[2] + p[5, 2, 2] # p[]
                sage: p([]).coproduct_iterated(3)
                p[] # p[] # p[] # p[]

                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi[2,2].coproduct_iterated(0)
                Psi[2, 2]
                sage: Psi[2,2].coproduct_iterated(3)
                Psi[] # Psi[] # Psi[] # Psi[2, 2] + 2*Psi[] # Psi[] # Psi[2] # Psi[2] + Psi[] # Psi[] # Psi[2, 2] # Psi[] + 2*Psi[] # Psi[2] # Psi[] # Psi[2] + 2*Psi[] # Psi[2] # Psi[2] # Psi[] + Psi[] # Psi[2, 2] # Psi[] # Psi[] + 2*Psi[2] # Psi[] # Psi[] # Psi[2] + 2*Psi[2] # Psi[] # Psi[2] # Psi[] + 2*Psi[2] # Psi[2] # Psi[] # Psi[] + Psi[2, 2] # Psi[] # Psi[] # Psi[]

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m[[1,3],[2]].coproduct_iterated(2)
                m{} # m{} # m{{1, 3}, {2}} + m{} # m{{1}} # m{{1, 2}} + m{} # m{{1, 2}} # m{{1}} + m{} # m{{1, 3}, {2}} # m{} + m{{1}} # m{} # m{{1, 2}} + m{{1}} # m{{1, 2}} # m{} + m{{1, 2}} # m{} # m{{1}} + m{{1, 2}} # m{{1}} # m{} + m{{1, 3}, {2}} # m{} # m{}

                sage: m[[]].coproduct_iterated(3), m[[1,3],[2]].coproduct_iterated(0)
                (m{} # m{} # m{} # m{}, m{{1, 3}, {2}})


            """
            if n < 0:
                raise ValueError("cannot take fewer than 0 coproduct iterations: %s < 0" % str(n))
            if n==0:
                return self
            elif n==1:
                return self.coproduct()
            else:
                # Use coassociativity of `\Delta` to perform many coproducts simultaneously.
                fn = floor(Integer(n-1)/2); cn = ceil(Integer(n-1)/2)
                def split(a,b): return tensor([a.coproduct_iterated(fn), b.coproduct_iterated(cn)])
                return (self.coproduct()).apply_multilinear_morphism(split)

