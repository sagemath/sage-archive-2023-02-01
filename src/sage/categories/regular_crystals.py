r"""
Regular Crystals
"""
#*****************************************************************************
#  Copyright (C) 2013    Anne Schilling <anne at math.ucdavis.edu>
#                        Travis Scrimshaw <tscrim at ucdavis.edu>
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
#****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals

class RegularCrystals(Category_singleton):
    r"""
    The category of regular crystals.

    A crystal is called *regular* if:

    .. MATH::

        \epsilon_i(b) = \max\{ k \mid e_i^k(b) \neq 0 \} \quad \text{and}
        \quad \phi_i(b) = \max\{ k \mid f_i^k(b) \neq 0 \}.

    .. NOTE::

        Regular crystals are sometimes referred to as *normal*. When only one
        of the conditions (on either `\phi_i` or `epsilon_i`) holds, these
        crystals are sometimes called *seminormal* or *semiregular*.

    EXAMPLES::

        sage: C = RegularCrystals()
        sage: C
        Category of regular crystals
        sage: C.super_categories()
        [Category of crystals]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = RegularCrystals().example()
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          running ._test_stembridge_local_axioms() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_stembridge_local_axioms() . . . pass
    """

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: RegularCrystals().super_categories()
            [Category of crystals]
        """
        return [Crystals()]

    def example(self, n = 3):
        """
        Returns an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = RegularCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.categories.crystals import Crystals
        return Crystals().example(n)

    class ParentMethods:

        # TODO: this could be a method in Crystals.Algebras.ElementMethods, so that
        # one could do:
        #   sage: C = CrystalOfTableaux(['A',2], shape=[2,1])
        #   sage: M = C.algebra(QQ)
        #   sage: m = M.an_element()
        #   sage: m.demazure_operator([1,4,2])
        def demazure_operator(self, element, reduced_word):
            r"""
            Returns the application of Demazure operators `D_i` for `i` from
            ``reduced_word`` on ``element``.

            INPUT:

            - ``element`` -- an element of a free module indexed by the
              underlying crystal
            - ``reduced_word`` -- a reduced word of the Weyl group of the
              same type as the underlying crystal

            OUTPUT:

            - an element of the free module indexed by the underlying crystal

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: C = CombinatorialFreeModule(QQ,T)
                sage: t = T.highest_weight_vector()
                sage: b = 2*C(t)
                sage: T.demazure_operator(b,[1,2,1])
                2*B[[[1, 1], [2]]] + 2*B[[[1, 2], [2]]] + 2*B[[[1, 3], [2]]] + 2*B[[[1, 1], [3]]]
                + 2*B[[[1, 2], [3]]] + 2*B[[[1, 3], [3]]] + 2*B[[[2, 2], [3]]] + 2*B[[[2, 3], [3]]]

            The Demazure operator is idempotent::

                sage: T = CrystalOfTableaux("A1",shape=[4])
                sage: C = CombinatorialFreeModule(QQ,T)
                sage: b = C(T.module_generators[0]); b
                B[[[1, 1, 1, 1]]]
                sage: e = T.demazure_operator(b,[1]); e
                B[[[1, 1, 1, 1]]] + B[[[1, 1, 1, 2]]] + B[[[1, 1, 2, 2]]] + B[[[1, 2, 2, 2]]] + B[[[2, 2, 2, 2]]]
                sage: e == T.demazure_operator(e,[1])
                True

                sage: all(T.demazure_operator(T.demazure_operator(C(t),[1]),[1]) == T.demazure_operator(C(t),[1]) for t in T)
                True
            """
            M = element.parent()
            for i in reversed(reduced_word):
                element = M.linear_combination((c.demazure_operator_simple(i), coeff)
                                               for c, coeff in element)
            return element

        def _test_stembridge_local_axioms(self, index_set=None, verbose=False, complete=False, **options):
            r"""
            This implements tests for the Stembridge local characterization
            on the finite crystal ``self``.

            The current implementation only uses the rules for simply-laced
            types.  Crystals of other types should still pass the test, but
            expansion of this test to non-simply laced type would be desirable.

            One can specify an index set smaller than the full index set of
            the crystal, using the option ``index_set``.

            Running with ``verbose=True`` will print each node for which a
            local axiom test applies.

            Running with ``complete=True`` will continue to run the test past
            the first failure of the local axioms.  This is probably only
            useful in conjunction with the verbose option, to see all places
            where the local axioms fail.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',3], shape=[2,1])
                sage: T._test_stembridge_local_axioms()
                True
                sage: T._test_stembridge_local_axioms(verbose=True)
                True
                sage: T._test_stembridge_local_axioms(index_set=[1,3])
                True

                sage: B=Crystals().example(choice='naive')
                sage: B._test_stembridge_local_axioms()
                Traceback (most recent call last):
                ...
                AssertionError: None
            """
            tester = self._tester(**options)
            goodness=True
            i = 0
            for x in self:
                goodness = x._test_stembridge_local_axioms(index_set, verbose)
                if goodness == False and not complete:
                    tester.fail()
                i += 1
                if i > tester._max_runs:
                    return
            tester.assertTrue(goodness)
            return goodness

    class ElementMethods:

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).epsilon(1)
                0
                sage: C(2).epsilon(1)
                1
            """
            assert i in self.index_set()
            x = self.e(i)
            eps = 0
            while x is not None:
                x = x.e(i)
                eps = eps + 1
            return eps

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).phi(1)
                1
                sage: C(2).phi(1)
                0
            """
            assert i in self.index_set()
            x = self.f(i)
            phi = 0
            while x is not None:
                x = x.f(i)
                phi = phi + 1
            return phi

        def weight(self):
            """
            Return the weight of this crystal element.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).weight()
                (1, 0, 0, 0, 0, 0)
            """
            return self.Phi() - self.Epsilon()

        def demazure_operator_simple(self, i, ring = None):
            r"""
            Return the Demazure operator `D_i` applied to ``self``.

            INPUT:

            - ``i`` -- an element of the index set of the underlying crystal
            - ``ring`` -- (default: ``QQ``) a ring

            OUTPUT:

            An element of the ``ring``-free module indexed by the underlying
            crystal.

            Let `r = \langle \mathrm{wt}(b), \alpha^{\vee}_i \rangle`, then
            `D_i(b)` is defined as follows:

            - If `r \geq 0`, this returns the sum of the elements obtained
              from ``self`` by application of `f_i^k` for `0 \leq k \leq r`.
            - If `r < 0`, this returns the opposite of the sum of the
              elements obtained by application of `e_i^k` for `0 < k < -r`.

            REFERENCES:

            .. [L1995] Peter Littelmann, Crystal graphs and Young tableaux,
               J. Algebra 175 (1995), no. 1, 65--87.

            .. [K1993] Masaki Kashiwara, The crystal base and Littelmann's
               refined Demazure character formula,
               Duke Math. J. 71 (1993), no. 3, 839--858.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: t = T(rows=[[1,2],[2]])
                sage: t.demazure_operator_simple(2)
                B[[[1, 2], [2]]] + B[[[1, 3], [2]]] + B[[[1, 3], [3]]]
                sage: t.demazure_operator_simple(2).parent()
                Free module generated by The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]] over Integer Ring

                sage: t.demazure_operator_simple(1)
                0

                sage: K = KirillovReshetikhinCrystal(['A',2,1],2,1)
                sage: t = K(rows=[[3],[2]])
                sage: t.demazure_operator_simple(0)
                B[[[2, 3]]] + B[[[1, 2]]]

            TESTS::

                sage: K = KirillovReshetikhinCrystal(['A',2,1],1,1)
                sage: x = K.an_element(); x
                [[1]]
                sage: x.demazure_operator_simple(0)
                0
                sage: x.demazure_operator_simple(0, ring = QQ).parent()
                Free module generated by Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1) over Rational Field
            """
            from sage.rings.integer_ring import ZZ
            if ring is None:
                ring = ZZ
            C = self.parent().algebra(ring)
            r = self.phi(i) - self.epsilon(i)
            if r >= 0:
                l = [self]
                element = self
                for k in range(r):
                    element = element.f(i)
                    l.append(element)
                return C.sum_of_monomials(l)
            else:
                l = []
                element = self
                for k in range(-r-1):
                    element = element.e(i)
                    l.append(element)
                return - C.sum_of_monomials(l)

        def stembridgeDelta_depth(self,i,j):
            r"""
            Return the difference in the `j`-depth of ``self`` and `e_i`
            of ``self``, where `i` and `j` are in the index set of the
            underlying crystal. This function is useful for checking the
            Stembridge local axioms for crystal bases.

            The `i`-depth of a crystal node `x` is `-\varepsilon_i(x)`.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: t=T(rows=[[1,2],[2]])
                sage: t.stembridgeDelta_depth(1,2)
                0
                sage: s=T(rows=[[2,3],[3]])
                sage: s.stembridgeDelta_depth(1,2)
                -1
            """
            if self.e(i) is None: return 0
            return -self.e(i).epsilon(j) + self.epsilon(j)

        def stembridgeDelta_rise(self,i,j):
            r"""
            Return the difference in the `j`-rise of ``self`` and `e_i` of
            ``self``, where `i` and `j` are in the index set of the
            underlying crystal. This function is useful for checking the
            Stembridge local axioms for crystal bases.

            The `i`-rise of a crystal node `x` is `\varphi_i(x)`.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: t=T(rows=[[1,2],[2]])
                sage: t.stembridgeDelta_rise(1,2)
                -1
                sage: s=T(rows=[[2,3],[3]])
                sage: s.stembridgeDelta_rise(1,2)
                0
            """
            if self.e(i) is None: return 0
            return self.e(i).phi(j) - self.phi(j)

        def stembridgeDel_depth(self,i,j):
            r"""
            Return the difference in the `j`-depth of ``self`` and `f_i` of
            ``self``, where `i` and `j` are in the index set of the
            underlying crystal. This function is useful for checking the
            Stembridge local axioms for crystal bases.

            The `i`-depth of a crystal node `x` is `\varepsilon_i(x)`.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: t=T(rows=[[1,1],[2]])
                sage: t.stembridgeDel_depth(1,2)
                0
                sage: s=T(rows=[[1,3],[3]])
                sage: s.stembridgeDel_depth(1,2)
                -1
            """
            if self.f(i) is None: return 0
            return -self.epsilon(j) + self.f(i).epsilon(j)

        def stembridgeDel_rise(self,i,j):
            r"""
            Return the difference in the `j`-rise of ``self`` and `f_i` of
            ``self``, where `i` and `j` are in the index set of the
            underlying crystal. This function is useful for checking the
            Stembridge local axioms for crystal bases.

            The `i`-rise of a crystal node `x` is `\varphi_i(x)`.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: t=T(rows=[[1,1],[2]])
                sage: t.stembridgeDel_rise(1,2)
                -1
                sage: s=T(rows=[[1,3],[3]])
                sage: s.stembridgeDel_rise(1,2)
                0
            """
            if self.f(i) is None: return 0
            return self.phi(j)-self.f(i).phi(j)

        def stembridgeTriple(self,i,j):
            r"""
            Let `A` be the Cartan matrix of the crystal, `x` a crystal element,
            and let `i` and `j` be in the index set of the crystal.
            Further, set
            ``b=stembridgeDelta_depth(x,i,j)``, and
            ``c=stembridgeDelta_rise(x,i,j))``.
            If ``x.e(i)`` is non-empty, this function returns the triple
            `( A_{ij}, b, c )`; otherwise it returns ``None``.
            By the Stembridge local characterization of crystal bases,
            one should have `A_{ij}=b+c`.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: t=T(rows=[[1,1],[2]])
                sage: t.stembridgeTriple(1,2)
                sage: s=T(rows=[[1,2],[2]])
                sage: s.stembridgeTriple(1,2)
                (-1, 0, -1)

                sage: T = CrystalOfTableaux(['B',2], shape=[2,1])
                sage: t=T(rows=[[1,2],[2]])
                sage: t.stembridgeTriple(1,2)
                (-2, 0, -2)
                sage: s=T(rows=[[-1,-1],[0]])
                sage: s.stembridgeTriple(1,2)
                (-2, -2, 0)
                sage: u=T(rows=[[0,2],[1]])
                sage: u.stembridgeTriple(1,2)
                (-2, -1, -1)
            """
            if self.e(i) is None: return None
            b=self.stembridgeDelta_depth(i,j)
            c=self.stembridgeDelta_rise(i,j)
            dd=self.cartan_type().dynkin_diagram()
            a=dd[j,i]
            return (a, b, c)

        def _test_stembridge_local_axioms(self, index_set=None, verbose=False, **options):
            r"""
            This implements tests for the Stembridge local characterization
            on the element of a crystal ``self``.

            The current implementation only uses the axioms for simply-laced
            types.  Crystals of other types should still pass the test, but
            in non-simply-laced types, passing is not a guarantee that the
            crystal arises from a representation.

            One can specify an index set smaller than the full index set of
            the crystal, using the option ``index_set``.

            Running with ``verbose=True`` will print warnings when a test fails.

            REFERENCES::

            .. [S2003] John R. Stembridge, A local characterization of
               simply-laced crystals,
               Transactions of the American Mathematical Society, Vol. 355,
               No. 12 (Dec., 2003), pp. 4807--4823

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: t=T(rows=[[1,1],[2]])
                sage: t._test_stembridge_local_axioms()
                True
                sage: t._test_stembridge_local_axioms(index_set=[1,3])
                True
                sage: t._test_stembridge_local_axioms(verbose=True)
                True
            """
            from sage.combinat.subset import Subsets
            tester = self._tester(**options)
            goodness=True
            if index_set is None: index_set=self.index_set()

            for (i,j) in Subsets(index_set, 2):
                if self.e(i) is not None and self.e(j) is not None:
                    triple=self.stembridgeTriple(i,j)
                    #Test axioms P3 and P4.
                    if not triple[0]==triple[1]+triple[2] or triple[1]>0 or triple[2]>0:
                        if verbose:
                            print 'Warning: Failed axiom P3 or P4 at vector ', self, 'i,j=', i, j, 'Stembridge triple:', self.stembridgeTriple(i,j)
                            goodness=False
                        else:
                            tester.fail()
                    if self.stembridgeDelta_depth(i,j)==0:
                        #check E_i E_j(x)= E_j E_i(x)
                        if self.e(i).e(j)!=self.e(j).e(i) or self.e(i).e(j).stembridgeDel_rise(j, i)!=0:
                            if verbose:
                                print 'Warning: Failed axiom P5 at: vector ', self, 'i,j=', i, j, 'Stembridge triple:', self.stembridgeTriple(i,j)
                                goodness=False
                            else:
                                tester.fail()
                    if self.stembridgeDelta_depth(i,j)==-1 and self.stembridgeDelta_depth(j,i)==-1:
                        #check E_i E_j^2 E_i (x)= E_j E_i^2 E_j (x)
                        y1=self.e(j).e(i).e(i).e(j)
                        y2=self.e(j).e(i).e(i).e(j)
                        a=y1.stembridgeDel_rise(j, i)
                        b=y2.stembridgeDel_rise(i, j)
                        if y1!=y2 or a!=-1 or b!=-1:
                            if verbose:
                                print 'Warning: Failed axiom P6 at: vector ', self, 'i,j=', i, j, 'Stembridge triple:', self.stembridgeTriple(i,j)
                                goodness=False
                            else:
                                tester.fail()
            tester.assertTrue(goodness)
            return goodness


