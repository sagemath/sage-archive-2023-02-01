r"""
Graded rings of modular forms for Hecke triangle groups

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

from sage.rings.all import FractionField, PolynomialRing, ZZ, QQ, infinity
from sage.algebras.free_algebra import FreeAlgebra

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method

from hecke_triangle_groups import HeckeTriangleGroup
from constructor import FormsRing, FormsSpace, rational_type


# Maybe replace Parent by just SageObject?
class FormsRing_abstract(Parent):
    r"""
    Abstract (Hecke) forms ring.

    This should never be called directly. Instead one should
    instantiate one of the derived classes of this class.
    """

    from graded_ring_element import FormsRingElement
    Element = FormsRingElement

    from analytic_type import AnalyticType
    AT = AnalyticType()

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Abstract (Hecke) forms ring.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: `\Z).

        - ``red_hom``    -- If ``True`` then results of binary operations are considered
                            homogeneous whenever it makes sense (default: ``False``).
                            This is mainly used by the (Hecke) forms.

        OUTPUT:

        The corresponding abstract (Hecke) forms ring.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: MR = ModularFormsRing(n=5, base_ring=ZZ, red_hom=True)
            sage: MR
            ModularFormsRing(n=5) over Integer Ring
            sage: MR.group()
            Hecke triangle group for n = 5
            sage: MR.base_ring()
            Integer Ring
            sage: MR.has_reduce_hom()
            True
            sage: MR.is_homogeneous()
            False
        """

        #from graded_ring import canonical_parameters
        #(group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)

        if (group == infinity):
            raise NotImplementedError

        #if (not group.is_arithmetic() and base_ring.characteristic()>0):
        #    raise NotImplementedError
        #if (base_ring.characteristic().divides(2*group.n()*(group.n()-2))):
        #    raise NotImplementedError
        if (base_ring.characteristic() > 0):
            raise NotImplementedError
        self._group               = group
        self._red_hom             = red_hom
        self._base_ring           = base_ring
        self._coeff_ring          = FractionField(PolynomialRing(base_ring,'d'))
        self._pol_ring            = PolynomialRing(base_ring,'x,y,z,d')
        self._rat_field           = FractionField(self._pol_ring)

        # default values
        self._weight              = None
        self._ep                  = None
        self._analytic_type       = self.AT(["quasi", "mero"])

        self.default_prec(10)
        self.disp_prec(5)
        self.default_num_prec(53)

        #super(FormsRing_abstract, self).__init__(self.coeff_ring())

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: QuasiModularFormsRing(n=4)
            QuasiModularFormsRing(n=4) over Integer Ring
        """

        return "{}FormsRing(n={}) over {}".format(self._analytic_type.analytic_space_name(), self._group.n(), self._base_ring)

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiWeakModularFormsRing
            sage: latex(QuasiWeakModularFormsRing())
            \mathcal{ QM^! }_{n=3}(\Bold{Z})
        """

        from sage.misc.latex import latex
        return "\\mathcal{{ {} }}_{{n={}}}({})".format(self._analytic_type.latex_space_name(), self._group.n(), latex(self._base_ring))

    def _element_constructor_(self, x):
        r"""
        Return ``x`` coerced/converted into this forms ring.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: MR = ModularFormsRing()
            sage: (x,y,z,d) = MR.pol_ring().gens()

            sage: MR(x^3)
            f_rho^3

            sage: el = MR.Delta().full_reduce()
            sage: MR(el)
            f_rho^3*d - f_i^2*d
            sage: el.parent() == MR
            False
            sage: MR(el).parent() == MR
            True
        """

        from graded_ring_element import FormsRingElement
        if isinstance(x, FormsRingElement):
            x = self._rat_field(x._rat)
        else:
            x = self._rat_field(x)
        return self.element_class(self, x)

    def _coerce_map_from_(self, S):
        r"""
        Return whether or not there exists a coercion from ``S`` to ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiWeakModularFormsRing, ModularFormsRing, CuspFormsRing
            sage: MR1 = QuasiWeakModularFormsRing(base_ring=CC)
            sage: MR2 = ModularFormsRing()
            sage: MR3 = CuspFormsRing()
            sage: MR3.has_coerce_map_from(MR2)
            False
            sage: MR1.has_coerce_map_from(MR2)
            True
            sage: MR2.has_coerce_map_from(MR3)
            True
            sage: MR3.has_coerce_map_from(ZZ)
            False
            sage: MR1.has_coerce_map_from(ZZ)
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms, CuspForms
            sage: MF2 = ModularForms(k=6, ep=-1)
            sage: MF3 = CuspForms(k=12, ep=1)
            sage: MR1.has_coerce_map_from(MF2)
            True
            sage: MR2.has_coerce_map_from(MF3)
            True
        """

        from space import FormsSpace_abstract
        if (    isinstance(S, FormsRing_abstract)\
            and self._group         == S._group\
            and self._analytic_type >= S._analytic_type\
            and self.base_ring().has_coerce_map_from(S.base_ring()) ):
                return True
        elif isinstance(S, FormsSpace_abstract):
            raise RuntimeError( "This case should not occur." )
            # return self._coerce_map_from_(S.graded_ring())
        elif (self.AT("holo") <= self._analytic_type) and (self.coeff_ring().has_coerce_map_from(S)):
            return True
        else:
            return False

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import CuspFormsRing
            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: CuspFormsRing().an_element()
            f_rho^3*d - f_i^2*d
            sage: CuspFormsRing().an_element() == CuspFormsRing().Delta()
            True
            sage: WeakModularForms().an_element()
            O(q^5)
            sage: WeakModularForms().an_element() == WeakModularForms().zero()
            True
        """

        return self(self.Delta())

    def default_prec(self, prec = None):
        r"""
        Set the default precision ``prec`` for the Fourier expansion.
        If ``prec=None`` (default) then the current default precision is returned instead.

        INPUT:

        - ``prec`` -- An integer.

        NOTE:

        This is also used as the default precision for the Fourier
        expansion when evaluating forms.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MR = ModularFormsRing()
            sage: MR.default_prec(3)
            sage: MR.default_prec()
            3
            sage: MR.Delta().q_expansion_fixed_d()
            q - 24*q^2 + O(q^3)
            sage: MF = ModularForms(k=4)
            sage: MF.default_prec(2)
            sage: MF.E4()
            1 + 240*q + O(q^2)
            sage: MF.default_prec()
            2
        """

        if (prec is not None):
            self._prec = ZZ(prec)
        else:
            return self._prec

    def disp_prec(self, prec = None):
        r"""
        Set the maximal display precision to ``prec``.
        If ``prec="max"`` the precision is set to the default precision.
        If ``prec=None`` (default) then the current display precision is returned instead.

        NOTE:

        This is used for displaying/representing (elements of)
        ``self`` as Fourier expansions.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=4)
            sage: MF.default_prec(5)
            sage: MF.disp_prec(3)
            sage: MF.disp_prec()
            3
            sage: MF.E4()
            1 + 240*q + 2160*q^2 + O(q^3)
            sage: MF.disp_prec("max")
            sage: MF.E4()
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + O(q^5)
        """

        if (prec == "max"):
            self._disp_prec = self._prec;
        elif (prec is not None):
            self._disp_prec = ZZ(prec)
        else:
            return self._disp_prec

    def default_num_prec(self, prec = None):
        r"""
        Set the default numerical precision to ``prec`` (default: ``53``).
        If ``prec=None`` (default) the current default numerical
        precision is returned instead.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=6)
            sage: MF.default_prec(20)
            sage: MF.default_num_prec(10)
            sage: MF.default_num_prec()
            10
            sage: E6 = MF.E6()
            sage: E6(i)
            -0.002...
            sage: MF.default_num_prec(100)
            sage: E6(i)
            0.000000000000000000000000000...

            sage: MF = ModularForms(n=5, k=4/3)
            sage: f_rho = MF.f_rho()
            sage: f_rho.q_expansion(prec=2)[1]
            7/(100*d)
            sage: MF.default_num_prec(15)
            sage: f_rho.q_expansion_fixed_d(prec=2)[1]
            9.9...
            sage: MF.default_num_prec(100)
            sage: f_rho.q_expansion_fixed_d(prec=2)[1]
            9.92593243510795915276017782...
        """

        if (prec is not None):
            self._num_prec = ZZ(prec)
        else:
            return self._num_prec

    def change_ring(self, new_base_ring):
        r"""
        Return the same space as ``self`` but over a new base ring ``new_base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().change_ring(CC)
            ModularFormsRing(n=3) over Complex Field with 53 bits of precision
        """

        return self.__class__.__base__(self._group, new_base_ring, self._red_hom)

    def graded_ring(self):
        r"""
        Return the graded ring containing ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing, CuspFormsRing
            sage: from sage.modular.modform_hecketriangle.space import CuspForms

            sage: MR = ModularFormsRing(n=5)
            sage: MR.graded_ring() == MR
            True

            sage: CF=CuspForms(k=12)
            sage: CF.graded_ring() == CuspFormsRing()
            False
            sage: CF.graded_ring() == CuspFormsRing(red_hom=True)
            True

            sage: CF.subspace([CF.Delta()]).graded_ring() == CuspFormsRing(red_hom=True)
            True
        """

        return self.extend_type(ring=True)

    def extend_type(self, analytic_type=None, ring=False):
        r"""
        Return a new space which contains (elements of) ``self`` with the analytic type
        of ``self`` extended by ``analytic_type``, possibly extended to a graded ring
        in case ``ring`` is ``True``.

        INPUT:

        - ``analytic_type``  -- An ``AnalyticType`` or something which
                                coerces into it (default: ``None``).

        - ``ring``           -- Whether to extend to a graded ring (default: ``False``).

        OUTPUT:

        The new extended space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import CuspForms

            sage: MR = ModularFormsRing(n=5)
            sage: MR.extend_type(["quasi", "weak"])
            QuasiWeakModularFormsRing(n=5) over Integer Ring

            sage: CF=CuspForms(k=12)
            sage: CF.extend_type("holo")
            ModularForms(n=3, k=12, ep=1) over Integer Ring
            sage: CF.extend_type("quasi", ring=True)
            QuasiCuspFormsRing(n=3) over Integer Ring

            sage: CF.subspace([CF.Delta()]).extend_type()
            CuspForms(n=3, k=12, ep=1) over Integer Ring
        """

        if analytic_type == None:
            analytic_type = self._analytic_type
        else:
            analytic_type = self._analytic_type.extend_by(analytic_type)

        if (ring or not self.is_homogeneous()):
            return FormsRing(analytic_type, group=self.group(), base_ring=self.base_ring(), red_hom=self.has_reduce_hom())
        else:
            return FormsSpace(analytic_type, group=self.group(), base_ring=self.base_ring(), k=self.weight(), ep=self.ep())

    def reduce_type(self, analytic_type=None, degree=None):
        r"""
        Return a new space with analytic properties shared by both ``self`` and ``analytic_type``,
        possibly reduced to its homogeneous space of the given ``degree`` (if ``degree`` is set).
        Elements of the new space are contained in ``self``.

        INPUT:

        - ``analytic_type``   -- An ``AnalyticType`` or something which coerces into it (default: ``None``).

        - ``degree``          -- ``None`` (default) or the degree of the homogeneous component to which
                                 ``self`` should be reduced.

        OUTPUT:

        The new reduced space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms

            sage: MR = QuasiModularFormsRing()
            sage: MR.reduce_type(["quasi", "cusp"])
            QuasiCuspFormsRing(n=3) over Integer Ring

            sage: MR.reduce_type("cusp", degree=(12,1))
            CuspForms(n=3, k=12, ep=1) over Integer Ring

            sage: MF=QuasiModularForms(k=6)
            sage: MF.reduce_type("holo")
            ModularForms(n=3, k=6, ep=-1) over Integer Ring

            sage: MF.reduce_type([])
            ZeroForms(n=3, k=6, ep=-1) over Integer Ring
        """

        if analytic_type == None:
            analytic_type = self._analytic_type
        else:
            analytic_type = self._analytic_type.reduce_to(analytic_type)

        if (degree == None and not self.is_homogeneous()):
            return FormsRing(analytic_type, group=self.group(), base_ring=self.base_ring(), red_hom=self.has_reduce_hom())
        elif (degree == None):
            return FormsSpace(analytic_type, group=self.group(), base_ring=self.base_ring(), k=self.weight(), ep=self.ep())
        else:
            (weight, ep) = degree
            if (self.is_homogeneous() and (weight != self.weight() or ep!=self.ep())):
                analytic_type = self._analytic_type.reduce_to([])
            return FormsSpace(analytic_type, group=self.group(), base_ring=self.base_ring(), k=weight, ep=ep)

    @cached_method
    def contains_coeff_ring(self):
        r"""
        Return whether ``self`` contains its coefficient ring.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import CuspFormsRing, ModularFormsRing
            sage: CuspFormsRing(n=4).contains_coeff_ring()
            False
            sage: ModularFormsRing(n=5).contains_coeff_ring()
            True
        """

        return (self.AT("holo") <= self._analytic_type)

    def construction(self):
        r"""
        Return a functor that constructs ``self`` (used by the coercion machinery).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().construction()
            (ModularFormsRingFunctor(n=3), BaseFacade(Integer Ring))
        """

        from functors import FormsRingFunctor, BaseFacade
        return FormsRingFunctor(self._analytic_type, self._group, self._red_hom), BaseFacade(self._base_ring)

    @cached_method
    def group(self):
        r"""
        Return the (Hecke triangle) group of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: MR = ModularFormsRing(n=7)
            sage: MR.group()
            Hecke triangle group for n = 7

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: CF = CuspForms(n=7, k=4/5)
            sage: CF.group()
            Hecke triangle group for n = 7
        """

        return self._group

    @cached_method
    def hecke_n(self):
        r"""
        Return the parameter ``n`` of the
        (Hecke triangle) group of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: MR = ModularFormsRing(n=7)
            sage: MR.hecke_n()
            7

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: CF = CuspForms(n=7, k=4/5)
            sage: CF.hecke_n()
            7
        """

        return self._group.n()

    @cached_method
    def base_ring(self):
        r"""
        Return base ring of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().base_ring()
            Integer Ring

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: CuspForms(k=12, base_ring=AA).base_ring()
            Algebraic Real Field
        """

        return self._base_ring

    @cached_method
    def coeff_ring(self):
        r"""
        Return coefficient ring of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().coeff_ring()
            Fraction Field of Univariate Polynomial Ring in d over Integer Ring

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: CuspForms(k=12, base_ring=AA).coeff_ring()
            Fraction Field of Univariate Polynomial Ring in d over Algebraic Real Field
        """

        return self._coeff_ring

    @cached_method
    def pol_ring(self):
        r"""
        Return the underlying polynomial ring used
        by ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().pol_ring()
            Multivariate Polynomial Ring in x, y, z, d over Integer Ring

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: CuspForms(k=12, base_ring=AA).pol_ring()
            Multivariate Polynomial Ring in x, y, z, d over Algebraic Real Field
        """

        return self._pol_ring

    @cached_method
    def rat_field(self):
        r"""
        Return the underlying rational field used by
        ``self`` to construct/represent elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().rat_field()
            Fraction Field of Multivariate Polynomial Ring in x, y, z, d over Integer Ring

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: CuspForms(k=12, base_ring=AA).rat_field()
            Fraction Field of Multivariate Polynomial Ring in x, y, z, d over Algebraic Real Field
        """

        return self._rat_field

    @cached_method
    def diff_alg(self):
        r"""
        Return the algebra of differential operators
        (over QQ) which is used on rational functions
        representing elements of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().diff_alg()
            Noncommutative Multivariate Polynomial Ring in X, Y, Z, dX, dY, dZ over Rational Field, nc-relations: {dY*Y: Y*dY + 1, dZ*Z: Z*dZ + 1, dX*X: X*dX + 1}

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: CuspForms(k=12, base_ring=AA).diff_alg()
            Noncommutative Multivariate Polynomial Ring in X, Y, Z, dX, dY, dZ over Rational Field, nc-relations: {dY*Y: Y*dY + 1, dZ*Z: Z*dZ + 1, dX*X: X*dX + 1}
        """

        # We only use two operators for now which do not involve 'd', so for performance
        # reason and due to restrictions for possible rings that can be used with algebra
        # relations we choose FractionField(base_ring) instead of self.coeff_ring().
        # For our purposes it is currently enough to define the operators over ZZ resp. QQ.
        free_alg         = FreeAlgebra(FractionField(ZZ),6,'X,Y,Z,dX,dY,dZ')
        (X,Y,Z,dX,dY,dZ) = free_alg.gens()
        diff_alg         = free_alg.g_algebra({dX*X:1+X*dX,dY*Y:1+Y*dY,dZ*Z:1+Z*dZ})
        
        return diff_alg

    @cached_method
    def _derivative_op(self):
        r"""
        Return the differential operator in ``self.diff_alg()``
        corresponding to the derivative of forms.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing(n=7)._derivative_op()
            -1/2*X^6*dY - 5/28*X^5*dZ + 1/7*X*Z*dX + 1/2*Y*Z*dY + 5/28*Z^2*dZ - 1/7*Y*dX
        """

        (X,Y,Z,dX,dY,dZ) = self.diff_alg().gens()
        return   1/self._group.n() * (X*Z-Y)*dX\
               + ZZ(1)/ZZ(2) * (Y*Z-X**(self._group.n()-1))*dY\
               + (self._group.n()-2) / (4*self._group.n()) * (Z**2-X**(self._group.n()-2))*dZ

    @cached_method
    def _serre_derivative_op(self):
        r"""
        Return the differential operator in ``self.diff_alg()``
        corresponding to the Serre derivative of forms.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing(n=8)._serre_derivative_op()
            -1/2*X^7*dY - 3/16*X^6*dZ - 3/16*Z^2*dZ - 1/8*Y*dX
        """

        (X,Y,Z,dX,dY,dZ) = self.diff_alg().gens()
        return - 1/self._group.n() * Y*dX\
               - ZZ(1)/ZZ(2) * X**(self._group.n()-1)*dY\
               - (self._group.n()-2) / (4*self._group.n()) * (Z**2+X**(self._group.n()-2))*dZ

    @cached_method
    def has_reduce_hom(self):
        r"""
        Return whether the method ``reduce`` should reduce
        homogeneous elements to the corresponding homogeneous space.

        This is mainly used by binary operations on homogeneous
        spaces which temporarily produce an element of ``self``
        but want to consider it as a homogeneous element
        (also see ``reduce``).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().has_reduce_hom()
            False
            sage: ModularFormsRing(red_hom=True).has_reduce_hom()
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ModularForms(k=6).has_reduce_hom()
            True
            sage: ModularForms(k=6).graded_ring().has_reduce_hom()
            True
        """

        return self._red_hom

    def is_homogeneous(self):
        r"""
        Return whether ``self`` is homogeneous component.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().is_homogeneous()
            False

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ModularForms(k=6).is_homogeneous()
            True
        """

        return self._weight != None

    def is_modular(self):
        r"""
        Return whether ``self`` only contains modular elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiWeakModularFormsRing, CuspFormsRing
            sage: QuasiWeakModularFormsRing().is_modular()
            False
            sage: CuspFormsRing(n=7).is_modular()
            True

            sage: from sage.modular.modform_hecketriangle.space import QuasiWeakModularForms, CuspForms
            sage: QuasiWeakModularForms(k=10).is_modular()
            False
            sage: CuspForms(n=7, k=12, base_ring=AA).is_modular()
            True
        """

        return not (self.AT("quasi") <= self._analytic_type)

    def is_weakly_holomorphic(self):
        r"""
        Return whether ``self`` only contains weakly
        holomorphic modular elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, QuasiWeakModularFormsRing, CuspFormsRing
            sage: QuasiMeromorphicModularFormsRing().is_weakly_holomorphic()
            False
            sage: QuasiWeakModularFormsRing().is_weakly_holomorphic()
            True

            sage: from sage.modular.modform_hecketriangle.space import MeromorphicModularForms, CuspForms
            sage: MeromorphicModularForms(k=10).is_weakly_holomorphic()
            False
            sage: CuspForms(n=7, k=12, base_ring=AA).is_weakly_holomorphic()
            True
        """

        return (self.AT("weak", "quasi") >= self._analytic_type)

    def is_holomorphic(self):
        r"""
        Return whether ``self`` only contains holomorphic
        modular elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiWeakModularFormsRing, QuasiModularFormsRing
            sage: QuasiWeakModularFormsRing().is_holomorphic()
            False
            sage: QuasiModularFormsRing().is_holomorphic()
            True

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms, CuspForms
            sage: WeakModularForms(k=10).is_holomorphic()
            False
            sage: CuspForms(n=7, k=12, base_ring=AA).is_holomorphic()
            True
        """

        return (self.AT("holo", "quasi") >= self._analytic_type)

    def is_cuspidal(self):
        r"""
        Return whether ``self`` only contains cuspidal elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing, QuasiCuspFormsRing
            sage: QuasiModularFormsRing().is_cuspidal()
            False
            sage: QuasiCuspFormsRing().is_cuspidal()
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms, QuasiCuspForms
            sage: ModularForms(k=12).is_cuspidal()
            False
            sage: QuasiCuspForms(k=12).is_cuspidal()
            True
        """

        return (self.AT("cusp", "quasi") >= self._analytic_type)

    def is_zerospace(self):
        r"""
        Return whether ``self`` is the (`0`-dimensional) zero space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().is_zerospace()
            False

            sage: from sage.modular.modform_hecketriangle.space import ModularForms, CuspForms
            sage: ModularForms(k=12).is_zerospace()
            False
            sage: CuspForms(k=12).reduce_type([]).is_zerospace()
            True
        """

        return (self.AT(["quasi"]) >= self._analytic_type)

    def analytic_type(self):
        r"""
        Return the analytic type of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, QuasiWeakModularFormsRing
            sage: QuasiMeromorphicModularFormsRing().analytic_type()
            quasi meromorphic modular
            sage: QuasiWeakModularFormsRing().analytic_type()
            quasi weakly holomorphic modular

            sage: from sage.modular.modform_hecketriangle.space import MeromorphicModularForms, CuspForms
            sage: MeromorphicModularForms(k=10).analytic_type()
            meromorphic modular
            sage: CuspForms(n=7, k=12, base_ring=AA).analytic_type()
            cuspidal
        """

        return self._analytic_type

    def homogeneous_space(self, k, ep):
        r"""
        Return the homogeneous component of degree (``k``, ``e``) of ``self``.

        INPUT:

        - `k` -- An integer.

        - ``ep`` -- `+1` or `-1`.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, QuasiWeakModularFormsRing
            sage: QuasiMeromorphicModularFormsRing(n=7).homogeneous_space(k=2, ep=-1)
            QuasiMeromorphicModularForms(n=7, k=2, ep=-1) over Integer Ring
        """

        return self.reduce_type(degree = (k,ep))

    @cached_method
    def J_inv(self):
        r"""
        Return the J-invariant (Hauptmodul) of the group of ``self``.
        It is normalized such that ``J_inv(infinity) = infinity``,
        it has real Fourier coefficients starting with ``d > 0`` and ``J_inv(i) = 1``

        It lies in a (weak) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, WeakModularFormsRing, CuspFormsRing
            sage: MR = WeakModularFormsRing(n=7)
            sage: J_inv = MR.J_inv()
            sage: J_inv in MR
            True
            sage: CuspFormsRing(n=7).J_inv() == J_inv
            True
            sage: J_inv
            f_rho^7/(f_rho^7 - f_i^2)
            sage: QuasiMeromorphicModularFormsRing(n=7).J_inv() == QuasiMeromorphicModularFormsRing(n=7)(J_inv)
            True

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms, CuspForms
            sage: MF = WeakModularForms(n=5, k=0)
            sage: J_inv = MF.J_inv()
            sage: J_inv in MF
            True
            sage: WeakModularFormsRing(n=5, red_hom=True).J_inv() == J_inv
            True
            sage: CuspForms(n=5, k=12).J_inv() == J_inv
            True
            sage: MF.disp_prec(3)
            sage: J_inv
            d*q^-1 + 79/200 + 42877/(640000*d)*q + 12957/(2000000*d^2)*q^2 + O(q^3)

            sage: WeakModularForms().J_inv()
            1/1728*q^-1 + 31/72 + 1823/16*q + 335840/27*q^2 + 16005555/32*q^3 + 11716352*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type("weak", ring=True)(x**self._group.n()/(x**self._group.n()-y**2)).reduce()

    @cached_method
    def j_inv(self):
        r"""
        Return the j-invariant (Hauptmodul) of the group of ``self``.
        It is normalized such that ``j_inv(infinity) = infinity``,
        and such that it has real Fourier coefficients starting with ``1``.

        It lies in a (weak) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, WeakModularFormsRing, CuspFormsRing
            sage: MR = WeakModularFormsRing(n=7)
            sage: j_inv = MR.j_inv()
            sage: j_inv in MR
            True
            sage: CuspFormsRing(n=7).j_inv() == j_inv
            True
            sage: j_inv
            f_rho^7/(f_rho^7*d - f_i^2*d)
            sage: QuasiMeromorphicModularFormsRing(n=7).j_inv() == QuasiMeromorphicModularFormsRing(n=7)(j_inv)
            True

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms, CuspForms
            sage: MF = WeakModularForms(n=5, k=0)
            sage: j_inv = MF.j_inv()
            sage: j_inv in MF
            True
            sage: WeakModularFormsRing(n=5, red_hom=True).j_inv() == j_inv
            True
            sage: CuspForms(n=5, k=12).j_inv() == j_inv
            True
            sage: MF.disp_prec(3)
            sage: j_inv
            q^-1 + 79/(200*d) + 42877/(640000*d^2)*q + 12957/(2000000*d^3)*q^2 + O(q^3)

            sage: WeakModularForms().j_inv()
            q^-1 + 744 + 196884*q + 21493760*q^2 + 864299970*q^3 + 20245856256*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type("weak", ring=True)(1/d*x**self._group.n()/(x**self._group.n()-y**2)).reduce()

    @cached_method
    def f_rho(self):
        r"""
        Return the generator ``f_rho`` of the graded ring of ``self``.
        Up to the group action ``f_rho`` has exactly one simple zero at ``rho``. ``f_rho`` is
        normalized such that its first nontrivial Fourier coefficient is ``1``.

        The polynomial variable ``x`` exactly corresponds to ``f_rho``.

        It lies in a (cuspidal) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, ModularFormsRing, CuspFormsRing
            sage: MR = ModularFormsRing(n=7)
            sage: f_rho = MR.f_rho()
            sage: f_rho in MR
            True
            sage: CuspFormsRing(n=7).f_rho() == f_rho
            True
            sage: f_rho
            f_rho
            sage: QuasiMeromorphicModularFormsRing(n=7).f_rho() == QuasiMeromorphicModularFormsRing(n=7)(f_rho)
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms, CuspForms
            sage: MF = ModularForms(n=5, k=4/3)
            sage: f_rho = MF.f_rho()
            sage: f_rho in MF
            True
            sage: ModularFormsRing(n=5, red_hom=True).f_rho() == f_rho
            True
            sage: CuspForms(n=5, k=12).f_rho() == f_rho
            True
            sage: MF.disp_prec(3)
            sage: f_rho
            1 + 7/(100*d)*q + 21/(160000*d^2)*q^2 + O(q^3)

            sage: ModularForms(k=4).f_rho() == ModularForms(k=4).E4()
            True
            sage: ModularForms(k=4).f_rho()
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type("holo", ring=True)(x).reduce()

    @cached_method
    def f_i(self):
        r"""
        Return the generator ``f_i`` of the graded ring of ``self``.
        Up to the group action ``f_i`` has exactly one simple zero at ``i``. ``f_i`` is
        normalized such that its first nontrivial Fourier coefficient is ``1``.

        The polynomial variable ``y`` exactly corresponds to ``f_i``.

        It lies in a (holomorphic) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, ModularFormsRing, CuspFormsRing
            sage: MR = ModularFormsRing(n=7)
            sage: f_i = MR.f_i()
            sage: f_i in MR
            True
            sage: CuspFormsRing(n=7).f_i() == f_i
            True
            sage: f_i
            f_i
            sage: QuasiMeromorphicModularFormsRing(n=7).f_i() == QuasiMeromorphicModularFormsRing(n=7)(f_i)
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms, CuspForms
            sage: MF = ModularForms(n=5, k=10/3)
            sage: f_i = MF.f_i()
            sage: f_i in MF
            True
            sage: ModularFormsRing(n=5, red_hom=True).f_i() == f_i
            True
            sage: CuspForms(n=5, k=12).f_i() == f_i
            True
            sage: MF.disp_prec(3)
            sage: f_i
            1 - 13/(40*d)*q - 351/(64000*d^2)*q^2 + O(q^3)

            sage: ModularForms(k=6).f_i() == ModularForms(k=4).E6()
            True
            sage: ModularForms(k=6).f_i()
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type("holo", ring=True)(y).reduce()

    @cached_method
    def f_inf(self):
        r"""
        Return the first nontrivial cusp form ``f_inf`` of the graded ring of ``self``.
        Up to the group action ``f_inf`` has exactly one simple zero at ``infinity``.
        ``f_inf`` is normalized such that its first nontrivial Fourier coefficient is ``1``.

        It lies in a (holomorphic) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, CuspFormsRing
            sage: MR = CuspFormsRing(n=7)
            sage: f_inf = MR.f_inf()
            sage: f_inf in MR
            True
            sage: f_inf
            f_rho^7*d - f_i^2*d
            sage: QuasiMeromorphicModularFormsRing(n=7).f_inf() == QuasiMeromorphicModularFormsRing(n=7)(f_inf)
            True

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(n=5, k=20/3)
            sage: f_inf = MF.f_inf()
            sage: f_inf in MF
            True
            sage: CuspFormsRing(n=5, red_hom=True).f_inf() == f_inf
            True
            sage: CuspForms(n=5, k=0).f_inf() == f_inf
            True
            sage: MF.disp_prec(3)
            sage: f_inf
            q - 9/(200*d)*q^2 + O(q^3)

            sage: CuspForms(k=12).f_inf() == CuspForms(k=12).Delta()
            True
            sage: CuspForms(k=12).f_inf()
            q - 24*q^2 + 252*q^3 - 1472*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type("cusp", ring=True)(d*(x**self._group.n()-y**2)).reduce()

    @cached_method
    def G_inv(self):
        r"""
        If `2` divides `n`: Return the G-invariant of the group of ``self``.

        The G-invariant is analogous to the J-invariant but has multiplier `-1`.
        I.e. ``G_inv(-1/t) = -G_inv(t)``. It is a holomorphic square root
        of ``J_inv*(J_inv-1)`` with real Fourier coefficients.

        If `2` does not divide `n` the function does not exist and an exception is raised.

        The G-invariant lies in a (weak) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, WeakModularFormsRing, CuspFormsRing
            sage: MR = WeakModularFormsRing(n=8)
            sage: G_inv = MR.G_inv()
            sage: G_inv in MR
            True
            sage: CuspFormsRing(n=8).G_inv() == G_inv
            True
            sage: G_inv
            f_rho^4*f_i*d/(f_rho^8 - f_i^2)
            sage: QuasiMeromorphicModularFormsRing(n=8).G_inv() == QuasiMeromorphicModularFormsRing(n=8)(G_inv)
            True

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms, CuspForms
            sage: MF = WeakModularForms(n=8, k=0, ep=-1)
            sage: G_inv = MF.G_inv()
            sage: G_inv in MF
            True
            sage: WeakModularFormsRing(n=8, red_hom=True).G_inv() == G_inv
            True
            sage: CuspForms(n=8, k=12, ep=1).G_inv() == G_inv
            True
            sage: MF.disp_prec(3)
            sage: G_inv
            d^2*q^-1 - 15*d/128 - 15139/262144*q - 11575/(1572864*d)*q^2 + O(q^3)

            sage: WeakModularForms(n=4, k=0, ep=-1).G_inv()
            1/65536*q^-1 - 3/8192 - 955/16384*q - 49/32*q^2 - 608799/32768*q^3 - 659/4*q^4 + O(q^5)

        As explained above, the G-invariant exists only for even `n`::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(n=9)
            sage: MF.G_inv()
            Traceback (most recent call last):
            ...
            ArithmeticError: G_inv doesn't exists for n=9.
        """

        if (ZZ(2).divides(self._group.n())):
            (x,y,z,d) = self._pol_ring.gens()
            return self.extend_type("weak", ring=True)(d*y*x**(self._group.n()/ZZ(2))/(x**self._group.n()-y**2)).reduce()
        else:
           raise ArithmeticError("G_inv doesn't exists for n={}.".format(self._group.n()))

    @cached_method
    def g_inv(self):
        r"""
        If `2` divides `n`: Return the g-invariant of the group of ``self``.

        The g-invariant is analogous to the j-invariant but has multiplier ``-1``.
        I.e. ``g_inv(-1/t) = -g_inv(t)``. It is a (normalized) holomorphic square root
        of ``J_inv*(J_inv-1)``, normalized such that its first nontrivial Fourier coefficient is ``1``.

        If `2` does not divide ``n`` the function does not exist and an exception is raised.

        The g-invariant lies in a (weak) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, WeakModularFormsRing, CuspFormsRing
            sage: MR = WeakModularFormsRing(n=8)
            sage: g_inv = MR.g_inv()
            sage: g_inv in MR
            True
            sage: CuspFormsRing(n=8).g_inv() == g_inv
            True
            sage: g_inv
            f_rho^4*f_i/(f_rho^8*d - f_i^2*d)
            sage: QuasiMeromorphicModularFormsRing(n=8).g_inv() == QuasiMeromorphicModularFormsRing(n=8)(g_inv)
            True

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms, CuspForms
            sage: MF = WeakModularForms(n=8, k=0, ep=-1)
            sage: g_inv = MF.g_inv()
            sage: g_inv in MF
            True
            sage: WeakModularFormsRing(n=8, red_hom=True).g_inv() == g_inv
            True
            sage: CuspForms(n=8, k=12, ep=1).g_inv() == g_inv
            True
            sage: MF.disp_prec(3)
            sage: g_inv
            q^-1 - 15/(128*d) - 15139/(262144*d^2)*q - 11575/(1572864*d^3)*q^2 + O(q^3)

            sage: WeakModularForms(n=4, k=0, ep=-1).g_inv()
            q^-1 - 24 - 3820*q - 100352*q^2 - 1217598*q^3 - 10797056*q^4 + O(q^5)

        As explained above, the g-invariant exists only for even `n`::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(n=9)
            sage: MF.g_inv()
            Traceback (most recent call last):
            ...
            ArithmeticError: g_inv doesn't exists for n=9.
        """

        if (ZZ(2).divides(self._group.n())):
            (x,y,z,d) = self._pol_ring.gens()
            return self.extend_type("weak", ring=True)(1/d*y*x**(self._group.n()/ZZ(2))/(x**self._group.n()-y**2)).reduce()
        else:
           raise ArithmeticError("g_inv doesn't exists for n={}.".format(self._group.n()))

    @cached_method
    def E4(self):
        r"""
        Return the normalized Eisenstein series of weight 4` of the graded ring of ``self``.
        It is equal to ``f_rho^(n-2)``.

        It lies in a (holomorphic) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, ModularFormsRing, CuspFormsRing
            sage: MR = ModularFormsRing(n=7)
            sage: E4 = MR.E4()
            sage: E4 in MR
            True
            sage: CuspFormsRing(n=7).E4() == E4
            True
            sage: E4
            f_rho^5
            sage: QuasiMeromorphicModularFormsRing(n=7).E4() == QuasiMeromorphicModularFormsRing(n=7)(E4)
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms, CuspForms
            sage: MF = ModularForms(n=5, k=4)
            sage: E4 = MF.E4()
            sage: E4 in MF
            True
            sage: ModularFormsRing(n=5, red_hom=True).E4() == E4
            True
            sage: CuspForms(n=5, k=12).E4() == E4
            True
            sage: MF.disp_prec(3)
            sage: E4
            1 + 21/(100*d)*q + 483/(32000*d^2)*q^2 + O(q^3)

            sage: ModularForms(k=4).f_rho() == ModularForms(k=4).E4()
            True
            sage: ModularForms(k=4).E4()
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type("holo", ring=True)(x**(self._group.n()-2)).reduce()

    @cached_method
    def E6(self):
        r"""
        Return the normalized Eisenstein series of weight `6` of the graded ring of ``self``,
        It is equal to ``f_rho^(n-3) * f_i``.

        It lies in a (holomorphic) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, ModularFormsRing, CuspFormsRing
            sage: MR = ModularFormsRing(n=7)
            sage: E6 = MR.E6()
            sage: E6 in MR
            True
            sage: CuspFormsRing(n=7).E6() == E6
            True
            sage: E6
            f_rho^4*f_i
            sage: QuasiMeromorphicModularFormsRing(n=7).E6() == QuasiMeromorphicModularFormsRing(n=7)(E6)
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms, CuspForms
            sage: MF = ModularForms(n=5, k=6)
            sage: E6 = MF.E6()
            sage: E6 in MF
            True
            sage: ModularFormsRing(n=5, red_hom=True).E6() == E6
            True
            sage: CuspForms(n=5, k=12).E6() == E6
            True
            sage: MF.disp_prec(3)
            sage: E6
            1 - 37/(200*d)*q - 14663/(320000*d^2)*q^2 + O(q^3)

            sage: ModularForms(k=6).f_i() == ModularForms(k=6).E6()
            True
            sage: ModularForms(k=6).E6()
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type("holo", ring=True)(x**(self._group.n()-3)*y).reduce()

    @cached_method
    def Delta(self):
        r"""
        Return an analog of the Delta-function of the graded ring of ``self``.
        It is a cusp form of weight `12` and is equal to
        ``d*(E4^3 - E6^2)`` or (in terms of the generators) ``d*x^(2*n-6)*(x^n - y^2)``.

        It lies in a (cuspidal) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, CuspFormsRing
            sage: MR = CuspFormsRing(n=7)
            sage: Delta = MR.Delta()
            sage: Delta in MR
            True
            sage: Delta
            f_rho^15*d - f_rho^8*f_i^2*d
            sage: QuasiMeromorphicModularFormsRing(n=7).Delta() == QuasiMeromorphicModularFormsRing(n=7)(Delta)
            True

            sage: from sage.modular.modform_hecketriangle.space import CuspForms, ModularForms
            sage: MF = CuspForms(n=5, k=12)
            sage: Delta = MF.Delta()
            sage: Delta in MF
            True
            sage: CuspFormsRing(n=5, red_hom=True).Delta() == Delta
            True
            sage: CuspForms(n=5, k=0).Delta() == Delta
            True
            sage: MF.disp_prec(3)
            sage: Delta
            q + 47/(200*d)*q^2 + O(q^3)

            sage: d = ModularForms(n=5).coeff_ring().gen()
            sage: Delta == (d*(ModularForms(n=5).E4()^3-ModularForms(n=5).E6()^2))
            True

            sage: CuspForms(k=12).f_inf() == CuspForms(k=12).Delta()
            True
            sage: CuspForms(k=12).Delta()
            q - 24*q^2 + 252*q^3 - 1472*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type("cusp", ring=True)(d*x**(2*self._group.n()-6)*(x**self._group.n()-y**2)).reduce()

    @cached_method
    def E2(self):
        r"""
        Return the normalized quasi holomorphic Eisenstein series of weight `2` of the
        graded ring of ``self``. It is also a generator of the graded ring of
        ``self`` and  the polynomial variable ``z`` exactly corresponds to ``E2``.

        It lies in a (quasi holomorphic) extension of the graded ring of ``self``.
        In case ``has_reduce_hom`` is ``True`` it is given as an element of
        the corresponding homogeneous space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing, QuasiModularFormsRing, CuspFormsRing
            sage: MR = QuasiModularFormsRing(n=7)
            sage: E2 = MR.E2()
            sage: E2 in MR
            True
            sage: CuspFormsRing(n=7).E2() == E2
            True
            sage: E2
            E2
            sage: QuasiMeromorphicModularFormsRing(n=7).E2() == QuasiMeromorphicModularFormsRing(n=7)(E2)
            True

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms, CuspForms
            sage: MF = QuasiModularForms(n=5, k=2)
            sage: E2 = MF.E2()
            sage: E2 in MF
            True
            sage: QuasiModularFormsRing(n=5, red_hom=True).E2() == E2
            True
            sage: CuspForms(n=5, k=12, ep=1).E2() == E2
            True
            sage: MF.disp_prec(3)
            sage: E2
            1 - 9/(200*d)*q - 369/(320000*d^2)*q^2 + O(q^3)

            sage: f_inf = MF.f_inf()
            sage: E2 == f_inf.derivative() / f_inf
            True

            sage: QuasiModularForms(k=2).E2()
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 + O(q^5)
        """

        (x,y,z,d) = self._pol_ring.gens()
        return self.extend_type(["holo", "quasi"], ring=True)(z).reduce()
