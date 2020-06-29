"""
Fusion Rings
"""
# ****************************************************************************
#  Copyright (C) 2019 Daniel Bump <bump at match.stanford.edu>
#                     Nicolas Thiery <nthiery at users.sf.net>
#                     Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.root_system.weyl_characters import WeylCharacterRing
from sage.combinat.q_analogues import q_int
from sage.matrix.special import diagonal_matrix
from sage.matrix.constructor import matrix
from sage.misc.misc import inject_variable
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import CyclotomicField
from sage.misc.cachefunc import cached_method

class FusionRing(WeylCharacterRing):
    r"""Return the Fusion Ring (Verlinde Algebra) of level ``k``.

    INPUT:

    - ``ct`` -- the Cartan type of a simple (finite-dimensional) Lie algebra
    - ``k`` -- a nonnegative integer

    This algebra has a basis (sometimes called *primary fields* but here
    called *simple objects*) indexed by the weights of level `\leq k`.
    These arise as the fusion algebras of WZW conformal field theories,
    or as Grothendieck groups of tilting modules for quantum groups at 
    roots of unity. The :class:`FusionRing` class is implemented as a
    variant of the :class:`WeylCharacterRing`.

    REFERENCES:

    - [BaKi2001]_ Chapter 3
    - [DFMS1996]_ Chapter 16
    - [EGNO2015]_ Chapter 8
    - [Feingold2004]_
    - [Fuchs1994]_
    - [Row2006]_
    - [Walton1990]_

    EXAMPLES::

        sage: A22 = FusionRing("A2",2)
        sage: [f1, f2] = A22.fundamental_weights()
        sage: M = [A22(x) for x in [0*f1, 2*f1, 2*f2, f1+f2, f2, f1]]
        sage: [M[3] * x for x in M]
        [A22(1,1),
         A22(0,1),
         A22(1,0),
         A22(0,0) + A22(1,1),
         A22(0,1) + A22(2,0),
         A22(1,0) + A22(0,2)]

    You may assign your own labels to the basis elements. In the next
    example, we create the `SO(5)` fusion ring of level `2`, check the
    weights of the basis elements, then assign new labels to them while
    injecting them into the global namespace::

        sage: B22 = FusionRing("B2", 2)
        sage: b = [B22(x) for x in B22.get_order()]; b
        [B22(0,0), B22(1,0), B22(0,1), B22(2,0), B22(1,1), B22(0,2)]
        sage: [x.weight() for x in b]
        [(0, 0), (1, 0), (1/2, 1/2), (2, 0), (3/2, 1/2), (1, 1)]
        sage: B22.fusion_labels(['I0','Y1','X','Z','Xp','Y2'], inject_variables=True)
        sage: b = [B22(x) for x in B22.get_order()]; b
        [I0, Y1, X, Z, Xp, Y2]
        sage: [(x, x.weight()) for x in b]
        [(I0, (0, 0)),
         (Y1, (1, 0)),
         (X, (1/2, 1/2)),
         (Z, (2, 0)),
         (Xp, (3/2, 1/2)),
         (Y2, (1, 1))]
        sage: X * Y1
        X + Xp
        sage: Z * Z
        I0

    A fixed order of the basis keys is available with :meth:`get_order`.
    This is the order used by methods such as :meth:`s_matrix`.
    You may use :meth:`CombinatorialFreeModule.set_order` to reorder the basis::

        sage: B22.set_order([x.weight() for x in [I0,Y1,Y2,X,Xp,Z]])
        sage: [B22(x) for x in B22.get_order()]
        [I0, Y1, Y2, X, Xp, Z]

    To reset the labels, you may run :meth:`fusion_labels` with no parameter::

        sage: B22.fusion_labels()
        sage: [B22(x) for x in B22.get_order()]
        [B22(0,0), B22(1,0), B22(0,2), B22(0,1), B22(1,1), B22(2,0)]

    To reset the order to the default, simply set it to the list of basis
    element keys::

        sage: B22.set_order(B22.basis().keys().list())
        sage: [B22(x) for x in B22.get_order()]
        [B22(0,0), B22(1,0), B22(0,1), B22(2,0), B22(1,1), B22(0,2)]

    The fusion ring has a number of methods that reflect its role
    as the Grothendieck ring of a modular tensor category. These
    include a twist method :meth:`Element.twist` for its elements related
    to the ribbon structure, and the S-matrix :meth:`s_ij`.

    There are two natural normalizations of the S-matrix. Both
    are explained in Chapter 3 of [BaKi2001]_. The one that is computed
    by the method :meth:`s_matrix`, or whose individual entries
    are computed by :meth:`s_ij` is denoted `\tilde{s}` in
    [BaKi2001]_. It is not unitary. To make it unitary, one would
    divide by the square root of `D = \sum_V d_i(V)^2` where the sum
    is over all simple objects `V` and `d_i(V)` is the quantum dimension
    of `V` computed by the method :meth:`Element.q_dimension`. The quantity `D`
    is computed by :meth:`global_q_dimension`.

    Let us check the Verlinde formula, which is [DFMS1996]_ (16.3). This
    famous identity states that

    .. MATH::

        N^k_{ij} = \sum_l \frac{s(i,\ell)\,s(j,\ell)\,\overline{s(k,\ell)}}{s(I,\ell)},

    where `N^k_{ij}` are the fusion coefficients, i.e. the structure
    constants of the fusion ring, and ``I`` is the unit object.
    The S-matrix has the property that if `i*` denotes the dual
    object of `i`, implemented in Sage as ``i.dual()``, then

    .. MATH::

        s(i*,j) = s(i,j*) = \overline{s(i,j)}.

    This is equation (16.5) in [DFMS1996]_. Thus with `N_{ijk}=N^{k*}_{ij}`
    the Verlinde formula is equivalent to

    .. MATH::

        N_{ijk} = \sum_l \frac{s(i,\ell)\,s(j,\ell)\,s(k,\ell)}{s(I,\ell)},

    In this formula `s` is the normalized unitary S-matrix
    denoted `s` in [BaKi2001]_. We may define a function that
    corresponds to the right-hand side, except using
    `\tilde{s}` instead of `s`::

        sage: def V(i,j,k):
        ....:     R = i.parent()
        ....:     return sum(R.s_ij(i,l) * R.s_ij(j,l) * R.s_ij(k,l) / R.s_ij(R.one(),l)
        ....:                for l in R.basis())

    This does not produce ``self.N_ijk(i,j,k)`` exactly, because of the
    missing normalization factor. The following code to check the
    Verlinde formula takes this into account::

       sage: def test_verlinde(R):
       ....:     b0 = R.one()
       ....:     c = R.global_q_dimension()
       ....:     return all(V(i,j,k) == c * R.N_ijk(i,j,k) for i in R.basis()
       ....:                for j in R.basis() for k in R.basis())

    Every fusion ring should pass this test::

        sage: test_verlinde(FusionRing("A2",1))
        True
        sage: test_verlinde(FusionRing("B4",2)) # long time (.56s)
        True

    As an exercise, the reader may verify the examples in
    Section 5.3 of [RoStWa2009]_. Here we check the example
    of the Ising modular tensor category, which is related
    to the BPZ minimal model `M(4,3)` or to an `E_8` coset
    model. See [DFMS1996]_ sections 7.4.2 and
    18.4.1. [RoStWa2009]_ Example 5.3.4 tells us how to
    construct it as the conjugate of the `E_8` level 2
    :class:`FusionRing`::

        sage: I = FusionRing("E8",2,conjugate=True)
        sage: I.fusion_labels(["i0","p","s"],inject_variables=True)
        sage: b = I.basis().list(); b
        [i0, p, s]
        sage: [[x*y for x in b] for y in b]
        [[i0, p, s], [p, i0, s], [s, s, i0 + p]]
        sage: [x.twist() for x in b]
        [0, 1, 1/8]
        sage: I.global_q_dimension()
        4
        sage: [x.q_dimension()^2 for x in b]
        [1, 1, 2]
        sage: I.s_matrix()
        [                       1                        1 -zeta128^48 + zeta128^16]
        [                       1                        1  zeta128^48 - zeta128^16]
        [-zeta128^48 + zeta128^16  zeta128^48 - zeta128^16                        0]
        sage: I.s_matrix().apply_map(lambda x:x^2)
        [1 1 2]
        [1 1 2]
        [2 2 0]

    """
    @staticmethod
    def __classcall__(cls, ct, k, base_ring=ZZ, prefix=None, style="coroots", conjugate=False):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: F1 = FusionRing('B3', 2)
            sage: F2 = FusionRing(CartanType('B3'), QQ(2), ZZ)
            sage: F3 = FusionRing(CartanType('B3'), int(2), style="coroots")
            sage: F1 is F2 and F2 is F3
            True

            sage: A23 = FusionRing('A2', 3)
            sage: TestSuite(A23).run()

            sage: B22 = FusionRing('B2', 2)
            sage: TestSuite(B22).run()

            sage: C31 = FusionRing('C3', 1)
            sage: TestSuite(C31).run()

            sage: D41 = FusionRing('D4', 1)
            sage: TestSuite(D41).run()

            sage: G22 = FusionRing('G2', 2)
            sage: TestSuite(G22).run()

            sage: F41 = FusionRing('F4', 1)
            sage: TestSuite(F41).run()

            sage: E61 = FusionRing('E6', 1)
            sage: TestSuite(E61).run()

            sage: E71 = FusionRing('E7', 1)
            sage: TestSuite(E71).run()

            sage: E81 = FusionRing('E8', 1)
            sage: TestSuite(E81).run()

        """
        return super(FusionRing, cls).__classcall__(cls, ct, base_ring=base_ring,
                                                    prefix=prefix, style=style, k=k, conjugate=conjugate)

    def _test_verlinde(self, **options):
        """
        Check the Verlinde formula for this :class:`FusionRing` instance.

        EXAMPLES::

            sage: G22=FusionRing("G2",2)
            sage: G22._test_verlinde()

        """
        tester = self._tester(**options)
        c = self.global_q_dimension()
        i0 = self.one()
        for x in self.basis():
            for y in self.basis():
                for z in self.basis():
                    v = sum(self.s_ij(x,w)*self.s_ij(y,w)*self.s_ij(z,w)/self.s_ij(i0,w) for w in self.basis())
                    tester.assertEqual(v,c*self.N_ijk(x,y,z))

    def field(self):
        r"""
        Return a cyclotomic field large enough to
        contain the `2\ell`-th roots of unity, as well as
        all the S-matrix entries.

        EXAMPLES::

            sage: FusionRing("A2",2).field()
            Cyclotomic Field of order 60 and degree 16
            sage: FusionRing("B2",2).field()
            Cyclotomic Field of order 40 and degree 16
        """
        return CyclotomicField(4 * self._fg * self._l)

    def get_order(self):
        r"""
        Return the weights of the basis vectors in a fixed order.

        You may change the order of the basis using :meth:`CombinatorialFreeModule.set_order`

        EXAMPLES::

            sage: A14=FusionRing("A1",4)
            sage: w = A14.get_order(); w
            [(0, 0), (1/2, -1/2), (1, -1), (3/2, -3/2), (2, -2)]
            sage: A14.set_order([w[k] for k in [0,4,1,3,2]])
            sage: [A14(x) for x in A14.get_order()]
            [A14(0), A14(4), A14(1), A14(3), A14(2)]

        .. WARNING::

            This duplicates :meth:`get_order` from
            :class:`CombinatorialFreeModule` except the result
            is *not* cached. Caching of
            :meth:`CombinatorialFreeModule.get_order` causes inconsistent
            results after calling :meth:`CombinatorialFreeModule.set_order`.

        """
        if self._order is None:
            self.set_order(self.basis().keys().list())
        return self._order

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: D41 = FusionRing('D4', 1)
            sage: D41.some_elements()
            [D41(1,0,0,0), D41(0,0,1,0), D41(0,0,0,1)]
        """
        return [self.monomial(x) for x in self.fundamental_weights()
                if self.level(x) <= self._k]

    def fusion_k(self):
        r"""
        Return the level of ``self``.

        EXAMPLES::

            sage: B22 = FusionRing('B2',2)
            sage: B22.fusion_k()
            2
        """
        return self._k

    def fusion_l(self):
        r"""
        Return the product `{\ell}=m_g(k + h^\vee)`, where `m_g` denotes the
        square of the ratio of the lengths of long to short roots of
        the underlying Lie algebra, `k` denotes the level of the FusionRing,
        and `h^\vee` denotes the dual Coxeter number of the underlying Lie
        algebra.

        This value is used to define the associated root `2\ell`-th
        of unity `q = e^{i\pi/\ell}`.

        EXAMPLES::

            sage: B22 = FusionRing('B2',2)
            sage: B22.fusion_l()
            10
            sage: D52 = FusionRing('D5',2)
            sage: D52.fusion_l()
            10
        """
        return self._l

    def twists_matrix(self):
        r"""
        Return a diagonal matrix describing the twist corresponding to
        each simple object in the ``FusionRing``.

        EXAMPLES::

            sage: B22 = FusionRing('B2',2)
            sage: B22.twists_matrix()
            [  0   0   0   0   0   0]
            [  0 4/5   0   0   0   0]
            [  0   0 1/2   0   0   0]
            [  0   0   0   2   0   0]
            [  0   0   0   0 3/2   0]
            [  0   0   0   0   0 6/5]
        """
        B = self.basis()
        return diagonal_matrix(B[x].twist() for x in self.get_order())

    def q_dims(self):
        r"""
        Return a list of quantum dimensions of the simple objects.

        EXAMPLES::

            sage: F41 = FusionRing('F4',1,conjugate=True)
            sage: F41.q_dims()
            [1, -zeta80^24 + zeta80^16 + 1]
            sage: B22 = FusionRing("B2",2)
            sage: B22.q_dims()
            [1, 2, -2*zeta40^12 + 2*zeta40^8 + 1, 1, -2*zeta40^12 + 2*zeta40^8 + 1, 2]
        """
        b = self.basis()
        return [b[x].q_dimension() for x in self.get_order()]

    @cached_method
    def N_ijk(self, elt_i, elt_j, elt_k):
        r"""
        Return the symmetric fusion coefficient `N_{ijk}`.

        INPUT:

        - ``elt_i``, ``elt_j``, ``elt_k`` -- elements of the fusion basis

        This is the same as `N_{ij}^{k\ast}`, where `N_{ij}^k` are
        the structure coefficients of the ring (see :meth:`Nk_ij`),
        and `k\ast`` denotes the dual element. The coefficient `N_{ijk}`
        is unchanged under permutations of the three basis vectors.

        EXAMPLES::

            sage: G23 = FusionRing("G2", 3)
            sage: G23.fusion_labels("g")
            sage: b = G23.basis().list(); b
            [g0, g1, g2, g3, g4, g5]
            sage: [(x,y,z) for x in b for y in b for z in b if G23.N_ijk(x,y,z) > 1]
            [(g3, g3, g3), (g3, g3, g4), (g3, g4, g3), (g4, g3, g3)]
            sage: all(G23.N_ijk(x,y,z)==G23.N_ijk(y,z,x) for x in b for y in b for z in b)
            True
            sage: all(G23.N_ijk(x,y,z)==G23.N_ijk(y,x,z) for x in b for y in b for z in b)
            True
        """
        return (elt_i * elt_j).monomial_coefficients().get(elt_k.dual().weight(), 0)

    @cached_method
    def Nk_ij(self, elt_i, elt_j, elt_k):
        r"""
        Return the fusion coefficient `N^k_{ij}`.

        These are the structure coefficients of the fusion ring, so

        .. MATH::

            i * j = \sum_{k} N_{ij}^k k.

        EXAMPLES::

            sage: A22 = FusionRing("A2", 2)
            sage: b = A22.basis().list()
            sage: all(x*y == sum(A22.Nk_ij(x,y,k)*k for k in b) for x in b for y in b)
            True
        """
        return (elt_i * elt_j).monomial_coefficients(copy=False).get(elt_k.weight(), 0)

    @cached_method
    def s_ij(self, elt_i, elt_j):
        r"""
        Return the element of the S-matrix of this FusionRing corresponding to
        the given elements.

        This is computed using the formula

        .. MATH::

            s_{i,j} = \frac{1}{\theta_i\theta_j} \sum_k N_{ik}^j d_k \theta_k,

        where `\theta_k` is the twist and `d_k` is the quantum
        dimension. See [Row2006]_ Equation (2.2) or [EGNO2015]_
        Proposition 8.13.8.

        INPUT:

        - ``elt_i``, ``elt_j`` -- elements of the fusion basis

        EXAMPLES::

            sage: G21 = FusionRing("G2", 1)
            sage: b = G21.basis()
            sage: [G21.s_ij(x, y) for x in b for y in b]
            [1, -zeta60^14 + zeta60^6 + zeta60^4, -zeta60^14 + zeta60^6 + zeta60^4, -1]
        """
        l = self.fusion_l() * self._fg
        K = self.field()
        q = K.gen()
        ijtwist = -2 * l * (elt_i.twist() + elt_j.twist())
        mc = (elt_i.dual() * elt_j).monomial_coefficients(copy=False)
        B = self.basis()
        return sum(B[k].q_dimension() * self.Nk_ij(elt_i, B[k], elt_j) * q**(2*l*B[k].twist() + ijtwist)
                   for k in mc)

    def s_matrix(self):
        r"""
        Return the S-matrix of this FusionRing.

        .. NOTE::

            This is the matrix denoted `\widetilde{s}` in [BaKi2001]_.
            It is not normalized to be unitary. To obtain a unitary
            matrix, divide by `\sqrt{D}` where `D` is computed
            by :meth:`global_q_dimension`.

        EXAMPLES::

            sage: D91 = FusionRing('D9', 1)
            sage: D91.s_matrix()
            [         1          1          1          1]
            [         1          1         -1         -1]
            [         1         -1 -zeta68^17  zeta68^17]
            [         1         -1  zeta68^17 -zeta68^17]

            sage: D41 = FusionRing('D4', 1)
            sage: D41.s_matrix()
            [ 1  1  1  1]
            [ 1  1 -1 -1]
            [ 1 -1  1 -1]
            [ 1 -1 -1  1]
        """
        b = self.basis()
        return matrix([[self.s_ij(b[x], b[y]) for x in self.get_order()] for y in self.get_order()])

    def fusion_labels(self, labels=None, inject_variables=False):
        r"""
        Set the labels of the basis.

        INPUT:

        - ``labels`` -- (default: ``None``) a list of strings or string
        - ``inject_variables`` -- (default: ``False``) if ``True``, then
          inject the variable names into the global namespace; note that
          this could override objects already defined

        If ``labels`` is a list, the length of the list must equal the
        number of basis elements. These become the names of
        the basis elements.

        If ``labels`` is a string, this is treated as a prefix and a
        list of names is generated.

        If ``labels`` is ``None``, then this resets the labels to the default.

        EXAMPLES::

            sage: A13 = FusionRing("A1", 3)
            sage: A13.fusion_labels("x")
            sage: fb = list(A13.basis()); fb
            [x0, x1, x2, x3]
            sage: Matrix([[x*y for y in A13.basis()] for x in A13.basis()])
            [     x0      x1      x2      x3]
            [     x1 x0 + x2 x1 + x3      x2]
            [     x2 x1 + x3 x0 + x2      x1]
            [     x3      x2      x1      x0]

        We give an example where the variables are injected into the
        global namepsace::

            sage: A13.fusion_labels("y", inject_variables=True)
            sage: y0
            y0
            sage: y0.parent() is A13
            True

        We reset the labels to the default::

            sage: A13.fusion_labels()
            sage: fb
            [A13(0), A13(1), A13(2), A13(3)]
            sage: y0
            A13(0)
        """
        if labels is None:
            # Remove the fusion labels
            self._fusion_labels = None
            return

        B = self.basis()
        if isinstance(labels, str):
            labels = [labels + str(k) for k in range(len(B))]
        elif len(labels) != len(B):
            raise ValueError('invalid data')

        d = {}
        ac = self.simple_coroots()
        for j, b in enumerate(self.get_order()):
            t = tuple([b.inner_product(x) for x in ac])
            d[t] = labels[j]
            if inject_variables:
                inject_variable(labels[j], B[b])
        self._fusion_labels = d

    def global_q_dimension(self):
        r"""
        Return `\sum d_i^2`, where the sum is over all simple objects
        and `d_i` is the quantum dimension.

        EXAMPLES::

            sage: FusionRing("E6",1).global_q_dimension()
            3
        """
        return sum(x.q_dimension()**2 for x in self.basis())

    class Element(WeylCharacterRing.Element):
        """
        A class for FusionRing elements.
        """
        def is_simple_object(self):
            r"""
            Determine whether ``self`` is a simple object of the fusion ring.

            EXAMPLES::

                sage: A22 = FusionRing("A2", 2)
                sage: x = A22(1,0); x
                A22(1,0)
                sage: x.is_simple_object()
                True
                sage: x^2
                A22(0,1) + A22(2,0)
                sage: (x^2).is_simple_object()
                False
            """
            return self.parent()._k is not None and len(self._monomial_coefficients) == 1

        def weight(self):
            r"""
            Return the parametrizing dominant weight in the level `k` alcove.

            This method is only available for basis elements.

            EXAMPLES::

                sage: A21 = FusionRing("A2",1)
                sage: sorted([x.weight() for x in A21.basis()])
                [(0, 0, 0), (1/3, 1/3, -2/3), (2/3, -1/3, -1/3)]
            """
            if len(self._monomial_coefficients) != 1:
                raise ValueError("fusion weight is valid for basis elements only")
            return next(iter(self._monomial_coefficients))

        def twist(self):
            r"""Compute the object's twist. This returns a rational number `h`
            such that  `\theta = e^{i \pi h}` is the twist of ``self``.

            This method is only available for simple objects. If
            `\lambda` is the weight of the object, then
            `h = \langle\lambda,\lambda+2\rho\rangle` where
            `\rho` is half the sum of the positive roots.
            As in [Row2006]_, this requires normalizing
            the invariant bilinear form so that `\langle \alpha, \alpha
            \rangle = 2` for short roots.

            EXAMPLES::

                sage: G21 = FusionRing('G2', 1)
                sage: G21.basis()
                Finite family {(0, 0, 0): G21(0,0), (1, 0, -1): G21(1,0)}
                sage: G21(1,0).twist()
                4/5
                sage: F41 = FusionRing('F4', 1, conjugate=True)
                sage: F41.basis()
                Finite family {(0, 0, 0, 0): F41(0,0,0,0), (1, 0, 0, 0): F41(0,0,0,1)}
                sage: F41(0,0,0,1).twist()
                4/5

            """
            if not self.is_simple_object():
                raise ValueError("quantum twist is only available for simple objects of a FusionRing")
            P = self.parent()
            rho = P.space().rho()
            lam = self.weight()
            inner = lam.inner_product(lam + 2*rho)
            twist = P._conj * P._nf * inner / P.fusion_l()
            # Reduce to canonical form
            while twist > 2:
                twist -= 2
            while twist < 0:
                twist += 2
            return twist

        @cached_method
        def q_dimension(self):
            r"""
            Return the quantum dimension as an element of the cyclotomic
            field of the `2\ell`-th roots of unity, where `l = m (k+h^\vee)`
            with `m=1,2,3` depending on whether type is simply, doubly or
            triply laced, `k` is the level and `h^\vee` is the dual
            Coxeter number.

            EXAMPLES::

                sage: B22 = FusionRing("B2",2)
                sage: [(b.q_dimension())^2 for b in B22.basis()]
                [1, 4, 5, 1, 5, 4]
            """
            if not self.is_simple_object():
                raise ValueError("quantum dimension is only available for simple objects of a FusionRing")
            P = self.parent()
            lam = self.weight()
            space = P.space()
            rho = space.rho()
            powers = {}
            for alpha in space.positive_roots():
                val = alpha.inner_product(lam + rho)
                if val in powers:
                    powers[val] += 1
                else:
                    powers[val] = 1
                val = alpha.inner_product(rho)
                if val in powers:
                    powers[val] -= 1
                else:
                    powers[val] = -1
            R = ZZ['q']
            q = R.gen()
            expr = R.fraction_field().one()
            for val in powers:
                exp = powers[val]
                if exp > 0:
                    expr *= q_int(P._nf * val, q)**exp
                elif exp < 0:
                    expr /= q_int(P._nf * val, q)**(-exp)
            expr = R(expr)
            expr = expr.substitute(q=q**4) / (q**(2*expr.degree()))
            zet = P.field().gen() ** P._fg
            return expr.substitute(q=zet)

