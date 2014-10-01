r"""
Hecke triangle group elements

AUTHORS:

- Jonas Jermann (2014): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import AA, ZZ, PolynomialRing, infinity, NumberField, AlgebraicField
from sage.misc.latex import latex
from sage.structure.parent_gens import localvars
from sage.rings.number_field.structure import AbsoluteFromRelative

from sage.groups.matrix_gps.group_element import MatrixGroupElement_generic
from sage.misc.cachefunc import cached_method


class HeckeTriangleGroupElement(MatrixGroupElement_generic):
    r"""
    Elements of HeckeTriangleGroup.
    """

    def __init__(self, parent, M, check=False, **kwargs):
        r"""
        An element of HeckeTriangle group given by a matrix ``M``.

        INPUT:

        - ``parent`` -- A ``HeckeTriangleGroup``.

        - ``M``      -- A matrix which coerces into the matrix space
                        of ``parent``. For example with entries in a
                        polynomial ring over ``ZZ`` with parameter ``lam``.

        - ``check``  -- ``True`` (default) or ``False``. If ``True``
                        then a (possibly long) check is performed
                        to see whether ``M`` really corresponds to a
                        group element of ``parent``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup, HeckeTriangleGroupElement
            sage: lam = PolynomialRing(ZZ, 'lam').gen()
            sage: M = matrix([[-1, 0], [-lam^4 + 5*lam^2 + lam - 5, -1]])
            sage: G = HeckeTriangleGroup(4)
            sage: G(M, check=True)
            Traceback (most recent call last):
            ...
            ValueError: The matrix is not an element of Hecke triangle group for n = 4, up to equivalence it identifies two nonequivalent points.

            sage: G = HeckeTriangleGroup(10)
            sage: HeckeTriangleGroupElement(G, M, check=True)
            [ -1   0]
            [lam  -1]
            sage: G(M).matrix().parent()
            Full MatrixSpace of 2 by 2 dense matrices over Maximal Order in Number Field in lam with defining polynomial x^4 - 5*x^2 + 5

            sage: M = matrix([[-1, lam], [0, 1]])
            sage: HeckeTriangleGroupElement(G, M, check=True)
            Traceback (most recent call last):
            ...
            ValueError: The matrix is not an element of Hecke triangle group for n = 10, it has determinant -1 != 1.

            sage: G.T().inverse()
            [   1 -lam]
            [   0    1]
            sage: G.U() == G.T()*G.S()
            True
            sage: G.U()^(-10) == -G.I()
            True
        """

        MatrixGroupElement_generic.__init__(self, parent, M, check=True, convert=True)

        if (check):
            if self._matrix.determinant() != 1:
                raise ValueError("The matrix is not an element of {}, it has determinant {} != 1.".format(parent, self._matrix.determinant()))
            self.decompose_basic()

    @cached_method
    def decompose_basic(self):
        r"""
        Decompose ``self`` into a product of the generators
        ``S`` and ``T`` of its parent,

        The function returns a tuple ``L`` consisting of either
        ``parent.S()`` or non-trivial integer powers of ``parent.T()``.
        Additionally ``L`` starts with +- the identity.
        It satisfies the property: ``prod(L) == self``.

        If this decomposition is not possible an ``ValueError``
        is raised. In particular this function can be used to
        check the membership in ``parent`` of an arbitrary matrix
        over the base ring.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=17)
            sage: L = (-G.V(2)).decompose_basic()
            sage: L
            (
            [  1 lam]  [ 0 -1]  [  1 lam]  [-1  0]
            [  0   1], [ 1  0], [  0   1], [ 0 -1]
            )
            sage: prod(L) == -G.V(2)
            True
            sage: L = G.U().decompose_basic()
            sage: L
            (
            [  1 lam]  [ 0 -1]  [1 0]
            [  0   1], [ 1  0], [0 1]
            )
            sage: prod(L) == G.U()
            True
        """

        def mshift(A):
            a = A[0][0]
            b = A[0][1]
            c = A[1][0]
            d = A[1][1]

            return (4*a*c + b*d) / (4*c*c + d*d)

        def mabs(A):
            a = A[0][0]
            b = A[0][1]
            c = A[1][0]
            d = A[1][1]

            return (4*a*a + b*b) / (4*c*c + d*d)

        res = []
        T   = self.parent().T()
        S   = self.parent().S()
        M   = self
        lam = self.parent().lam()

        while True:
            m = (AA(mshift(M) / lam) + ZZ(1)/ZZ(2)).floor()
            M = T**(-m) * M
            if (m != 0):
                res.append(T**m)

            abs_t = mabs(M)
            if AA(abs_t) < 1:
                M = (-S) * M
                res.append(S)
            elif M.is_identity():
                res.append(M)
                return tuple(res)
            else:
                raise ValueError("The matrix is not an element of {}, up to equivalence it identifies two nonequivalent points.".format(self.parent()))

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: V = HeckeTriangleGroup(17).V(13)
            sage: latex(V)
            \begin{pmatrix} \mathit{\lambda}^{3} - 2 \mathit{\lambda} & \mathit{\lambda}^{2} - 1 \\ \mathit{\lambda}^{4} - 3 \mathit{\lambda}^{2} + 1 & \mathit{\lambda}^{3} - 2 \mathit{\lambda} \end{pmatrix}
        """

        latex_out = r"\begin{pmatrix} %s & %s \\ %s & %s \end{pmatrix}"%(latex(self.a()), latex(self.b()), latex(self.c()), latex(self.d()))
        return latex_out.replace("lam", r"\lambda")

    def __neg__(self):
        r"""
        Return the group element corresponding to the negative of the underlying matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: -U
            [-lam    1]
            [  -1    0]
        """

        return self.parent().element_class(self.parent(), -self._matrix, check=False)

    def __getitem__(self, key):
        r"""
        Return the corresponding rows/entries of the underlying matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U[0]
            (lam, -1)
            sage: U[0].parent()
            Ambient free module of rank 2 over the principal ideal domain Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
            sage: U[1][0]
            1
            sage: U[1][0].parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
        """

        return self._matrix.__getitem__(key)

    def a(self):
        r"""
        Return the upper left entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.a()
            lam
            sage: U.a().parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
        """

        return self._matrix[0][0]

    def b(self):
        r"""
        Return the upper right entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.b()
            -1
            sage: U.b().parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
        """

        return self._matrix[0][1]

    def c(self):
        r"""
        Return the lower left entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.c()
            1
        """

        return self._matrix[1][0]

    def d(self):
        r"""
        Return the lower right of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.d()
            0
        """

        return self._matrix[1][1]

    def is_translation(self, exclude_one=False):
        r"""
        Return whether ``self`` is a translation. If ``exclude_one = True``,
        then the identity map is not considered as a translation.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: (-HeckeTriangleGroup(n=7).T(-4)).is_translation()
            True
            sage: (-HeckeTriangleGroup(n=7).I()).is_translation()
            True
            sage: (-HeckeTriangleGroup(n=7).I()).is_translation(exclude_one=True)
            False
        """

        exclude_one = bool(exclude_one)
        a = self.a()
        b = self.b()
        c = self.c()
        d = self.d()

        if (c != 0) or (a != d) or (a != 1 and a != -1):
            return False
        elif (exclude_one and b == 0):
            return False
        else:
            return True

    def is_reflection(self):
        r"""
        Return whether ``self`` is the usual reflection on the unit circle.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: (-HeckeTriangleGroup(n=7).S()).is_reflection()
            True
            sage: HeckeTriangleGroup(n=7).U().is_reflection()
            False
        """

        if (self == self.parent().S()) or (self == -self.parent().S()):
            return True
        else:
            return False

    def trace(self):
        r"""
        Return the trace of ``self``, which is the sum of the diagonal entries.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.U().trace()
            lam
            sage: G.S().trace()
            0
        """

        return self._matrix.trace()

    def discriminant(self):
        r"""
        Return the discriminant of ``self`` which corresponds to
        the discriminant of the corresponding quadratic form of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.V(3).discriminant()
            4*lam^2 + 4*lam - 4
            sage: AA(G.V(3).discriminant())
            16.19566935808922?
        """

        return self.trace()**2 - 4

    def is_hyperbolic(self):
        r"""
        Return whether ``self`` is a hyperbolic matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_hyperbolic() for k in range(1,8) ]
            [False, True, True, True, True, False, False]
            sage: G.U().is_hyperbolic()
            False
        """

        return AA(self.discriminant()) > 0

    def is_parabolic(self, exclude_one=False):
        r"""
        Return whether ``self`` is a parabolic matrix.

        If ``exclude_one`` is set, then +- the identity
        element is not considered parabolic.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_parabolic() for k in range(1,8) ]
            [True, False, False, False, False, True, False]
            sage: G.U().is_parabolic()
            False
            sage: G.V(6).is_parabolic(exclude_one=True)
            True
            sage: G.V(7).is_parabolic(exclude_one=True)
            False
        """

        if (exclude_one and self.is_identity()):
            return False

        return self.discriminant() == 0

    def is_identity(self):
        r"""
        Return whether ``self`` is the identity or minus the identity.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_identity() for k in range(1,8) ]
            [False, False, False, False, False, False, False]
            sage: G.U().is_identity()
            False
        """

        if self == self.parent().I() or self == -self.parent().I():
            return True
        else:
            return False

    def is_elliptic(self):
        r"""
        Return whether ``self`` is an elliptic matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_elliptic() for k in range(1,8) ]
            [False, False, False, False, False, False, True]
            sage: G.U().is_elliptic()
            True
        """

        return AA(self.discriminant()) < 0

    def root_extension_field(self):
        r"""
        Return the extension field of the base field of (the parent of)
        ``self`` in which the corresponding fixed points of ``self`` lie.

        The variable name of the corresponding generator is ``e``.
        It corresponds to ``e=sqrt(D)``, where ``D`` is the discriminant of ``self``.
        If the extension degree is ``1`` then the base field is returned.

        The correct embedding is the positive resp. positive imaginary one.
        Unfortunately no default embedding can be specified for relative
        number fields yet.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=infinity)
            sage: G.V(3).root_extension_field()
            Number Field in e with defining polynomial x^2 - 32
            sage: G.T().root_extension_field() == G.base_field()
            True
            sage: (G.S()).root_extension_field()
            Number Field in e with defining polynomial x^2 + 4

            sage: G = HeckeTriangleGroup(n=7)
            sage: G.V(3).root_extension_field()
            Number Field in e with defining polynomial x^2 - 4*lam^2 - 4*lam + 4 over its base field
            sage: G.V(1).root_extension_field() == G.base_field()
            True
            sage: G.U().root_extension_field()
            Number Field in e with defining polynomial x^2 - lam^2 + 4 over its base field
        """

        K = self.parent().base_field()
        x = PolynomialRing(K, 'x').gen()
        D = self.discriminant()

        if D.is_square():
            return K
        else:
            # unfortunately we can't set embeddings for relative extensions :-(
            # return K.extension(x**2 - D, 'e', embedding=AA(D).sqrt())

            L = K.extension(x**2 - D, 'e')

            #e = AA(D).sqrt()
            #emb = L.hom([e])
            #L._unset_embedding()
            #L.register_embedding(emb)

            #return NumberField(L.absolute_polynomial(), 'e', structure=AbsoluteFromRelative(L), embedding=(???))
            return L

    @cached_method
    def root_extension_embedding(self, K=AlgebraicField()):
        r"""
        Return the correct embedding from the root extension field to ``K`` (default: ``AlgebraicField()``).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=infinity)

            sage: fp = (-G.S()).fixed_points()[0]
            sage: alg_fp = (-G.S()).root_extension_embedding()(fp)
            sage: alg_fp
            1*I
            sage: alg_fp == (-G.S()).fixed_points(embedded=True)[0]
            True

            sage: fp = (-G.V(2)).fixed_points()[1]
            sage: alg_fp = (-G.V(2)).root_extension_embedding(AA)(fp)
            sage: alg_fp
            -1.732050807568...?
            sage: alg_fp == (-G.V(2)).fixed_points(embedded=True)[1]
            True

            sage: fp = (-G.V(2)).fixed_points()[0]
            sage: alg_fp = (-G.V(2)).root_extension_embedding(AA)(fp)
            sage: alg_fp
            1.732050807568...?
            sage: alg_fp == (-G.V(2)).fixed_points(embedded=True)[0]
            True

            sage: G = HeckeTriangleGroup(n=7)

            sage: fp = (-G.S()).fixed_points()[1]
            sage: alg_fp = (-G.S()).root_extension_embedding()(fp)
            sage: alg_fp
            0.?... - 1.000000000000...?*I
            sage: alg_fp == (-G.S()).fixed_points(embedded=True)[1]
            True

            sage: fp = (-G.U()^4).fixed_points()[0]
            sage: alg_fp = (-G.U()^4).root_extension_embedding()(fp)
            sage: alg_fp
            0.9009688679024...? + 0.4338837391175...?*I
            sage: alg_fp == (-G.U()^4).fixed_points(embedded=True)[0]
            True

            sage: (-G.U()^4).root_extension_embedding(CC)(fp)
            0.900968867902... + 0.433883739117...*I
            sage: (-G.U()^4).root_extension_embedding(CC)(fp).parent()
            Complex Field with 53 bits of precision

            sage: fp = (-G.V(5)).fixed_points()[1]
            sage: alg_fp = (-G.V(5)).root_extension_embedding(AA)(fp)
            sage: alg_fp
            -0.6671145837954...?
            sage: alg_fp == (-G.V(5)).fixed_points(embedded=True)[1]
            True
        """

        G = self.parent()
        D = self.discriminant()
        F = self.root_extension_field()
        n = G.n()

        if D.is_square():
            e   = D.sqrt()
            lam = G.lam()
        elif n in [3, infinity]:
            e   = F.gen(0)
            lam = G.lam()
        else:
            (e, lam) = F.gens()

        emb_lam  = K(G.lam())
        if self.is_elliptic():
            emb_e = K(AlgebraicField()(D).sqrt())
        else:
            emb_e = K(AA(D).sqrt())

        guess = ZZ(0)
        min_value = infinity
        index = ZZ(0)

        for (index, emb) in enumerate(self.root_extension_field().embeddings(K)):
            if K.is_exact():
                if emb(lam) == emb_lam and emb(e) == emb_e:
                    return emb
            else:
                value = (emb(lam) - emb_lam).n(K.prec()).abs() + (emb(e) - emb_e).n(K.prec()).abs()
                if (value < min_value):
                    guess = index
                    min_value = value

        if K.is_exact() or min_value == infinity:
            raise ValueError("No suitable embedding is available for K = {}!".format(K))
        else:
            return self.root_extension_field().embeddings(K)[guess]

    def fixed_points(self, embedded=False, order="default"):
        r"""
        Return a pair of (mutually conjugate) fixed points of ``self``
        in a possible quadratic extension of the base field.

        INPUT:

        - ``embedded`` -- If ``True`` the fixed points are embedded into
                          ``AlgebraicRealField`` resp. ``AlgebraicField``.
                          Default: ``False``.

        - ``order``    -- If ``order="none"`` the fixed points are choosen
                          and ordered according to a fixed formula.

                          If ``order="sign"`` the fixed points are always ordered
                          according to the sign in front of the square root.

                          If ``order="default"`` (default) then in case the fixed
                          points are hyperbolic they are ordered according to the
                          sign of the trace of ``self`` instead, such that the
                          attracting fixed point comes first.

                          If ``order="trace"`` the fixed points are always ordered
                          according to the sign of the trace of ``self``.
                          If the trace is zero they are ordered by the sign in
                          front of the square root. In particular the fixed_points
                          in this case remain the same for ``-self``.

        OUTPUT:

        If ``embedded=True`` an element of either ``AlgebraicRealField`` or
        ``AlgebraicField`` is returned. Otherwise an element of a relative field
        extension over the base field of (the parent of) ``self`` is returned.

        Warning: Relative field extensions don't support default embeddings.
        So the correct embedding (which is the positive resp. imaginary positive
        one) has to be choosen.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=infinity)
            sage: (-G.T(-4)).fixed_points()
            [+Infinity, +Infinity]
            sage: (-G.S()).fixed_points()
            [1/2*e, -1/2*e]
            sage: p = (-G.S()).fixed_points(embedded=True)[0]
            sage: p
            1*I
            sage: (-G.S()).acton(p) == p
            True
            sage: (-G.V(2)).fixed_points()
            [1/2*e, -1/2*e]
            sage: (-G.V(2)).fixed_points() == G.V(2).fixed_points()
            True
            sage: p = (-G.V(2)).fixed_points(embedded=True)[1]
            sage: p
            -1.732050807568878?
            sage: (-G.V(2)).acton(p) == p
            True

            sage: G = HeckeTriangleGroup(n=7)
            sage: (-G.S()).fixed_points()
            [1/2*e, -1/2*e]
            sage: p = (-G.S()).fixed_points(embedded=True)[1]
            sage: p
            -1*I
            sage: (-G.S()).acton(p) == p
            True
            sage: (G.U()^4).fixed_points()
            [(1/2*lam^2 - 1/2*lam - 1/2)*e + 1/2*lam, (-1/2*lam^2 + 1/2*lam + 1/2)*e + 1/2*lam]
            sage: pts = (G.U()^4).fixed_points(order="trace")
            sage: (G.U()^4).fixed_points() == [pts[1], pts[0]]
            True
            sage: (G.U()^4).fixed_points(order="trace") == (-G.U()^4).fixed_points(order="trace")
            True
            sage: (G.U()^4).fixed_points() == (G.U()^4).fixed_points(order="none")
            True
            sage: (-G.U()^4).fixed_points() == (G.U()^4).fixed_points()
            True
            sage: (-G.U()^4).fixed_points(order="none") == pts
            True
            sage: p = (G.U()^4).fixed_points(embedded=True)[1]
            sage: p
            0.9009688679024191? - 0.4338837391175581?*I
            sage: (G.U()^4).acton(p) == p
            True
            sage: (-G.V(5)).fixed_points()
            [(1/2*lam^2 - 1/2*lam - 1/2)*e, (-1/2*lam^2 + 1/2*lam + 1/2)*e]
            sage: (-G.V(5)).fixed_points() == G.V(5).fixed_points()
            True
            sage: p = (-G.V(5)).fixed_points(embedded=True)[0]
            sage: p
            0.6671145837954892?
            sage: (-G.V(5)).acton(p) == p
            True
        """

        if self.c() == 0:
            return [infinity, infinity]
        else:
            D = self.discriminant()
            if D.is_square():
                e = D.sqrt()
            else:
                e = self.root_extension_field().gen()

            a = self.a()
            d = self.d()
            c = self.c()

            if order == "none":
                sgn = ZZ(1)
            elif order == "sign":
                sgn = AA(c).sign()
            elif order == "default":
                if self.is_elliptic() or self.trace() == 0:
                    sgn = AA(c).sign()
                else:
                    sgn = AA(self.trace()).sign()
            elif order == "trace":
                if self.trace() == 0:
                    sgn = AA(c).sign()
                else:
                    sgn = AA(self.trace()).sign()
            else:
                raise NotImplementedError

            if embedded:
                e = AA(D).sqrt()
                a = AA(a)
                d = AA(d)
                c = AA(c)

            root1 = (a-d)/(2*c) + sgn*e/(2*c)
            root2 = (a-d)/(2*c) - sgn*e/(2*c)

            return [root1, root2]

    def acton(self, z):
        r"""
        Return the image of ``z`` under the action of ``self``
        by linear fractional transformations or by conjugation
        in case ``z`` is an element of the parent of ``self``.

        .. NOTE:

        There is a 1-1 correspondence between hyperbolic
        fixed points and the corresponding primitive element
        in the stabilizer. The action in the two cases above
        is compatible with this correspondence.

        INPUT:

        - ``z``     -- Either an element of ``self`` or any
                       element to which a linear fractional
                       transformation can be applied in
                       the usual way.

                       In particular ``infinity`` is a possible
                       argument and a possible return value.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(5)
            sage: G.S().acton(1 + i/2)
            2/5*I - 4/5
            sage: G.S().acton(1 + i/2).parent()
            Symbolic Ring
            sage: G.S().acton(i + exp(-2))
            -1/(e^(-2) + I)
            sage: G.S().acton(i + exp(-2)).parent()
            Symbolic Ring

            sage: G.T().acton(infinity) == infinity
            True
            sage: G.U().acton(infinity)
            lam
            sage: G.V(2).acton(-G.lam()) == infinity
            True

            sage: G.V(2).acton(G.U()) == G.V(2)*G.U()*G.V(2).inverse()
            True
            sage: G.V(2).inverse().acton(G.U())
            [  0  -1]
            [  1 lam]
        """

        if z.parent() == self.parent():
            return self*z*self.inverse()

        if z == infinity and self.c() == 0:
            return infinity
        elif z == infinity:
            return self.a()/self.c()
        elif self.c() != 0 and z == -self.d()/self.c():
            return infinity
        else:
            return (self.a()*z + self.b()) / (self.c()*z + self.d())

    # def _act_on_(self, other, self_on_left):
    #     TODO: implement default actions for "suitable" x
    #     if (self_on_left):
    #         return self.acton(other)
