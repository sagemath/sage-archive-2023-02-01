r"""
Hecke triangle groups

AUTHORS:

- Jonas Jermann (2013): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, QQ, AA, AlgebraicField, infinity, NumberField
from sage.functions.all import cos,exp,sec
from sage.functions.other import psi1
from sage.symbolic.all import pi,i
from sage.matrix.constructor import matrix
from sage.misc.latex import latex

from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_generic
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from hecke_triangle_group_element import HeckeTriangleGroupElement


class HeckeTriangleGroup(FinitelyGeneratedMatrixGroup_generic, UniqueRepresentation):
    r"""
    Hecke triangle group (2, n, infinity).
    """

    Element = HeckeTriangleGroupElement

    @staticmethod
    def __classcall__(cls, n=3):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(QQ(3)) == HeckeTriangleGroup(int(3))
            True
        """

        if (n == infinity):
            n = infinity
        else:
            n = ZZ(n)
            if (n < 3):
                raise AttributeError("n has to be infinity or an Integer >= 3.")

        return super(HeckeTriangleGroup, cls).__classcall__(cls, n)

    def __init__(self, n):
        r"""
        Hecke triangle group (2, n, infinity).
        Namely the von Dyck group corresponding to the triangle group
        with angles (pi/2, pi/n, 0).

        INPUT:

        - ``n``   - ``infinity`` or an integer greater or equal to ``3``.

        OUTPUT:

        The Hecke triangle group for the given parameter ``n``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(12)
            sage: G
            Hecke triangle group for n = 12
            sage: G.category()
            Category of groups
        """

        self._n = n

        if n in [3, infinity]:
            self._base_ring = ZZ
            self._lam = ZZ(1) if n==3 else ZZ(2)
        else:
            lam_symbolic = 2*cos(pi/n)
            K = NumberField(self.lam_minpoly(), 'lam', embedding = AA(lam_symbolic))
            #self._base_ring = K.order(K.gens())
            self._base_ring = K.maximal_order()
            self._lam = self._base_ring.gen(1)

        T = matrix(self._base_ring, [[1,self._lam],[0,1]])
        S = matrix(self._base_ring, [[0,-1],[1,0]])

        FinitelyGeneratedMatrixGroup_generic.__init__(self, ZZ(2), self._base_ring, [S, T])

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10)
            Hecke triangle group for n = 10
        """

        return "Hecke triangle group for n = {}".format(self._n)

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: a = HeckeTriangleGroup(5)
            sage: latex(a)
            \Gamma^{(5)}
        """

        return '\\Gamma^{(%s)}'%(latex(self._n))

    def one_element(self):
        r"""
        Return the identity element/matrix for ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10).one_element()
            [1 0]
            [0 1]

            sage: HeckeTriangleGroup(10).one_element().parent()
            Hecke triangle group for n = 10
        """

        return self.I()

    @cached_method
    def lam_minpoly(self):
        r"""
        Return the minimal polynomial of the corresponding lambda parameter of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10).lam_minpoly()
            x^4 - 5*x^2 + 5
            sage: HeckeTriangleGroup(17).lam_minpoly()
            x^8 - x^7 - 7*x^6 + 6*x^5 + 15*x^4 - 10*x^3 - 10*x^2 + 4*x + 1
            sage: HeckeTriangleGroup(infinity).lam_minpoly()
            x - 2
        """

        # TODO: Write an explicit (faster) implementation
        lam_symbolic = 2*cos(pi/self._n)
        return AA(lam_symbolic).minpoly()

    def base_ring(self):
        r"""
        Return the base ring of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(n=infinity).base_ring()
            Integer Ring
            sage: HeckeTriangleGroup(n=7).base_ring()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
        """

        return self._base_ring

    def base_field(self):
        r"""
        Return the base field of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(n=infinity).base_field()
            Rational Field
            sage: HeckeTriangleGroup(n=7).base_field()
            Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
        """

        if self._n in [3, infinity]:
            return QQ
        else:
            return self._base_ring.number_field()

    def n(self):
        r"""
        Return the parameter ``n`` of ``self``, where
        ``pi/n`` is the angle at ``rho`` of the corresponding
        basic hyperbolic triangle with vertices ``i``, ``rho``
        and ``infinity``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10).n()
            10
            sage: HeckeTriangleGroup(infinity).n()
            +Infinity
        """

        return self._n

    # TODO: rename this to a more descriptive lambda in a later update/patch
    def lam(self):
        r"""
        Return the parameter ``lambda`` of ``self``,
        where ``lambda`` is twice the real part of ``rho``,
        lying between ``1`` (when ``n=3``) and ``2`` (when ``n=infinity``).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).lam()
            1
            sage: HeckeTriangleGroup(4).lam()
            lam
            sage: HeckeTriangleGroup(4).lam()^2
            2
            sage: HeckeTriangleGroup(6).lam()^2
            3
            sage: AA(HeckeTriangleGroup(10).lam())
            1.9021130325903...?
            sage: HeckeTriangleGroup(infinity).lam()
            2
        """

        return self._lam

    def rho(self):
        r"""
        Return the vertex ``rho`` of the basic hyperbolic
        triangle which describes ``self``. ``rho`` has
        absolute value 1 and angle ``pi/n``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).rho() == 1/2 + sqrt(3)/2*i
            True
            sage: HeckeTriangleGroup(4).rho() == sqrt(2)/2*(1 + i)
            True
            sage: HeckeTriangleGroup(6).rho() == sqrt(3)/2 + 1/2*i
            True
            sage: HeckeTriangleGroup(10).rho()
            0.95105651629515...? + 0.30901699437494...?*I
            sage: HeckeTriangleGroup(infinity).rho()
            1
        """

        # TODO: maybe rho should be replaced by -rhobar
        if (self._n == infinity):
            return AA(1)
        else:
            return AlgebraicField()(exp(pi/self._n*i))

    def alpha(self):
        r"""
        Return the parameter ``alpha`` of ``self``.
        This is the first parameter of the hypergeometric series used
        in the calculation of the Hauptmodul of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).alpha()
            1/12
            sage: HeckeTriangleGroup(4).alpha()
            1/8
            sage: HeckeTriangleGroup(5).alpha()
            3/20
            sage: HeckeTriangleGroup(6).alpha()
            1/6
            sage: HeckeTriangleGroup(10).alpha()
            1/5
            sage: HeckeTriangleGroup(infinity).alpha()
            1/4
        """

        return ZZ(1)/ZZ(2) * (ZZ(1)/ZZ(2) - ZZ(1)/self._n)

    def beta(self):
        r"""
        Return the parameter ``beta`` of ``self``.
        This is the second parameter of the hypergeometric series used
        in the calculation of the Hauptmodul of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).beta()
            5/12
            sage: HeckeTriangleGroup(4).beta()
            3/8
            sage: HeckeTriangleGroup(5).beta()
            7/20
            sage: HeckeTriangleGroup(6).beta()
            1/3
            sage: HeckeTriangleGroup(10).beta()
            3/10
            sage: HeckeTriangleGroup(infinity).beta()
            1/4
        """

        return ZZ(1)/ZZ(2) * (ZZ(1)/ZZ(2) + ZZ(1)/self._n)

    @cached_method
    def I(self):
        r"""
        Return the identity element/matrix for ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(10).I()
            [1 0]
            [0 1]

            sage: HeckeTriangleGroup(10).I().parent()
            Hecke triangle group for n = 10
        """

        return self(matrix(self._base_ring, [[1,0],[0,1]]))

    @cached_method
    def T(self, m=1):
        r"""
        Return the element in ``self`` corresponding to the translation by ``m*self.lam()``.

        INPUT:

        - ``m`` -- An integer, default: ``1``, namely the second generator of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).T()
            [1 1]
            [0 1]
            sage: HeckeTriangleGroup(10).T(-4)
            [     1 -4*lam]
            [     0      1]
            sage: HeckeTriangleGroup(10).T().parent()
            Hecke triangle group for n = 10
        """

        return self(matrix(self._base_ring, [[1,self._lam*m],[0,1]]))

    @cached_method
    def S(self):
        r"""
        Return the generator of ``self`` corresponding to the
        conformal circle inversion.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).S()
            [ 0 -1]
            [ 1  0]
            sage: HeckeTriangleGroup(10).S()
            [ 0 -1]
            [ 1  0]
            sage: HeckeTriangleGroup(10).S()^2 == -HeckeTriangleGroup(10).I()
            True
            sage: HeckeTriangleGroup(10).S()^4 == HeckeTriangleGroup(10).I()
            True

            sage: HeckeTriangleGroup(10).S().parent()
            Hecke triangle group for n = 10
        """

        return self.gen(0)

    @cached_method
    def U(self):
        r"""
        Return an alternative generator of ``self`` instead of ``T``.
        ``U`` stabilizes ``rho`` and has order ``2*self.n()``.

        If ``n=infinity`` then ``U`` is parabolic and has infinite order,
        it then fixes the cusp ``[-1]``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).U()
            [ 1 -1]
            [ 1  0]
            sage: HeckeTriangleGroup(3).U()^3 == -HeckeTriangleGroup(3).I()
            True
            sage: HeckeTriangleGroup(3).U()^6 == HeckeTriangleGroup(3).I()
            True
            sage: HeckeTriangleGroup(10).U()
            [lam  -1]
            [  1   0]
            sage: HeckeTriangleGroup(10).U()^10 == -HeckeTriangleGroup(10).I()
            True
            sage: HeckeTriangleGroup(10).U()^20 == HeckeTriangleGroup(10).I()
            True

            sage: HeckeTriangleGroup(10).U().parent()
            Hecke triangle group for n = 10
        """

        return self.T() * self.S()

    def V(self, j):
        r"""
        Return the j'th generator for the usual representatives of
        conjugacy classes of ``self``. It is given by ``V=U^(j-1)*T``.

        INPUT:

        - ``j``  -- Any integer. To get the usual representatives
                    ``j`` should range from ``1`` to ``self.n()-1``.

        OUTPUT:

        The corresponding matrix/element.
        The matrix is parabolic if ``j`` is congruent to +-1 modulo ``self.n()``.
        It is elliptic if ``j`` is congruent to 0 modulo ``self.n()``.
        It is hyperbolic otherwise.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(3)
            sage: G.V(0) == -G.S()
            True
            sage: G.V(1) == G.T()
            True
            sage: G.V(2)
            [1 0]
            [1 1]
            sage: G.V(3) == G.S()
            True

            sage: G = HeckeTriangleGroup(5)
            sage: G.V(1)
            [  1 lam]
            [  0   1]
            sage: G.V(2)
            [lam lam]
            [  1 lam]
            sage: G.V(3)
            [lam   1]
            [lam lam]
            sage: G.V(4)
            [  1   0]
            [lam   1]
            sage: G.V(5) == G.S()
            True
        """

        return self.U()**(j-1) * self.T()

    @cached_method
    def dvalue(self):
        r"""
        Return a symbolic expression (or an exact value in case n=3, 4, 6)
        for the transfinite diameter (or capacity) of ``self``.
        I.e. the first nontrivial Fourier coefficient of the Hauptmodul
        for the Hecke triangle group in case it is normalized to ``J_inv(i)=1``.

        EXAMPLES:

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).dvalue()
            1/1728
            sage: HeckeTriangleGroup(4).dvalue()
            1/256
            sage: HeckeTriangleGroup(5).dvalue()
            e^(2*euler_gamma - 4*pi/(sqrt(5) + 1) + psi(17/20) + psi(13/20))
            sage: HeckeTriangleGroup(6).dvalue()
            1/108
            sage: HeckeTriangleGroup(10).dvalue()
            e^(2*euler_gamma - pi*sec(1/10*pi) + psi(4/5) + psi(7/10))
            sage: HeckeTriangleGroup(infinity).dvalue()
            1/64
        """

        n = self._n
        if (n==3):
            return ZZ(1)/ZZ(2**6*3**3)
        elif (n==4):
            return ZZ(1)/ZZ(2**8)
        elif (n==6):
            return ZZ(1)/ZZ(2**2*3**3)
        elif (n==infinity):
            return ZZ(1)/ZZ(2**6)
        else:
            return exp(-ZZ(2)*psi1(ZZ(1)) + psi1(ZZ(1)-self.alpha())+psi1(ZZ(1)-self.beta()) - pi*sec(pi/self._n))


    def is_arithmetic(self):
        r"""
        Return True if ``self`` is an arithmetic subgroup.

        EXAMPLES:

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(3).is_arithmetic()
            True
            sage: HeckeTriangleGroup(4).is_arithmetic()
            True
            sage: HeckeTriangleGroup(5).is_arithmetic()
            False
            sage: HeckeTriangleGroup(6).is_arithmetic()
            True
            sage: HeckeTriangleGroup(10).is_arithmetic()
            False
            sage: HeckeTriangleGroup(infinity).is_arithmetic()
            True
        """

        if (self._n in [ZZ(3), ZZ(4), ZZ(6), infinity]):
            return True
        else:
            return False

    def get_FD(self, z):
        r"""
        Return a tuple (A,w) which determines how to map ``z``
        to the usual (strict) fundamental domain of ``self``.

        INPUT:

        - ``z`` -- a complex number or an element of AlgebraicField().

        OUTPUT:

        A tuple ``(A, w)``.

        - ``A`` -- a matrix in ``self`` such that ``A.acton(w)==z``
          (if ``z`` is exact at least).

        - ``w`` -- a complex number or an element of AlgebraicField()
          (depending on the type ``z``) which lies inside the (strict)
          fundamental domain of ``self`` (``self.in_FD(w)==True``) and
          which is equivalent to ``z`` (by the above property).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(8)
            sage: z = AlgebraicField()(1+i/2)
            sage: (A, w) = G.get_FD(z)
            sage: A
            [-lam    1]
            [  -1    0]
            sage: A.acton(w) == z
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: z = (134.12 + 0.22*i).n()
            sage: (A, w) = G.get_FD(z)
            sage: A
            [-73*lam^3 + 74*lam       73*lam^2 - 1]
            [        -lam^2 + 1                lam]
            sage: w
            0.769070776942... + 0.779859114103...*I
            sage: z
            134.120000000... + 0.220000000000...*I
            sage: A.acton(w)
            134.1200000... + 0.2200000000...*I
        """

        ID = self.I()
        T  = self.T()
        S  = self.S()
        TI = self.T(-1)

        A = ID
        w = z
        while (abs(w) < ZZ(1) or abs(w.real()) > self.lam()/ZZ(2)):
            if (abs(w) < ZZ(1)):
                w = self.S().acton(w)
                A = S*A
            while (w.real() >= self.lam()/ZZ(2)):
                w = TI.acton(w)
                A = TI*A
            while (w.real() < -self.lam()/ZZ(2)):
                w = T.acton(w)
                A = T*A
        if (w.real() == self.lam()/ZZ(2)):
            w = TI.acton(w)
            A = TI*A
        if (abs(w) == ZZ(1) and w.real() > ZZ(0)):
            w = S.acton(w)
            A = S*A

        AI = A.inverse()

        return (AI, A.acton(z))

    def in_FD(self, z):
        r"""
        Returns ``True`` if ``z`` lies in the (strict) fundamental
        domain of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: HeckeTriangleGroup(5).in_FD(CC(1.5/2 + 0.9*i))
            True
            sage: HeckeTriangleGroup(4).in_FD(CC(1.5/2 + 0.9*i))
            False
        """

        return self.get_FD(z)[0] == self.I()
