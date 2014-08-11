r"""
Elements of graded rings of modular forms for Hecke triangle groups

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

from sage.rings.all import ZZ, infinity, LaurentSeries, O
from sage.functions.all import exp
from sage.symbolic.all import pi, i
from sage.structure.parent_gens import localvars

from sage.structure.element import CommutativeAlgebraElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from constructor import rational_type, FormsSpace, FormsRing
from series_constructor import MFSeriesConstructor


# Warning: We choose CommutativeAlgebraElement because we want the
# corresponding operations (e.g. __mul__) even though the category
# (and class) of the parent is in some cases not
# CommutativeAlgebras but Modules
class FormsRingElement(CommutativeAlgebraElement, UniqueRepresentation):
    r"""
    Element of a FormsRing.
    """

    from analytic_type import AnalyticType
    AT = AnalyticType()

    @staticmethod
    def __classcall__(cls, parent, rat):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring_element import FormsRingElement
            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: (x,d) = var("x","d")
            sage: el = FormsRingElement(ModularFormsRing(), x*d)
            sage: el.rat()
            x*d
            sage: el.rat().parent()
            Fraction Field of Multivariate Polynomial Ring in x, y, z, d over Integer Ring
        """

        rat = parent.rat_field()(rat)
        # rat.reduce() <- maybe add this for the nonexact case

        return super(FormsRingElement,cls).__classcall__(cls, parent, rat)

    def __init__(self, parent, rat):
        r"""
        Element of a FormsRing ``parent`` corresponding to the rational
        function ``rat`` evaluated at ``x=f_rho``, ``y=f_i``, ``z=E2``
        and ``d`` by the formal parameter from ``parent.coeff_ring()``.

        The functions ``f_rho, f_i, E2`` can be obtained from
        ``self.parent().graded_ring()``.

        INPUT:

        - ``parent`` -- An (non abstract) instance of ``FormsRing_abstract``.

        - ``rat``    -- A rational function in ``parent.rat_field()``, the
                        fraction field of the polynomial ring in ``x,y,z,d``
                        over the base ring of ``parent``.
        
        OUTPUT:

        An element of ``parent``. If ``rat`` does not correspond to such
        an element an exception is raised.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: MR = QuasiModularFormsRing(n=5)
            sage: el = MR(x^3*d + y*z)
            sage: el
            f_rho^3*d + f_i*E2
            sage: el.rat()
            x^3*d + y*z
            sage: el.parent()
            QuasiModularFormsRing(n=5) over Integer Ring
            sage: el.rat().parent()
            Fraction Field of Multivariate Polynomial Ring in x, y, z, d over Integer Ring
        """

        self._rat = rat
        (elem, homo, self._weight, self._ep, self._analytic_type) = rational_type(rat, parent.hecke_n(), parent.base_ring())

        if not (
            elem and\
            self._analytic_type <= parent.analytic_type() ):
                raise Exception("{} does not correspond to an element of the {}.".format(rat, parent))

        super(FormsRingElement, self).__init__(parent)

    # Unfortunately the polynomial ring does not give unique
    # representations of elements (with respect to ==)
    def __eq__(self, other):
        r"""
        Return whether ``self`` is equal to ``other``.
        They are considered equal if the corresponding rational
        functions are equal and the groups match up.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import MeromorphicModularFormsRing
            sage: (x,y,z,d) = MeromorphicModularFormsRing().pol_ring().gens()
            sage: MeromorphicModularFormsRing(n=3)(x) == MeromorphicModularFormsRing(n=4)(x)
            False
            sage: MeromorphicModularFormsRing()(-1/x) is MeromorphicModularFormsRing()(1/(-x))
            False
            sage: MeromorphicModularFormsRing()(-1/x) == MeromorphicModularFormsRing()(1/(-x))
            True
            sage: MeromorphicModularFormsRing(base_ring=CC)(-1/x) == MeromorphicModularFormsRing()(1/(-x))
            True
        """

        if (super(FormsRingElement, self).__eq__(other)):
            return True
        elif (isinstance(other, FormsRingElement)):
            if (self.group() == other.group()):
                if (self.group().is_arithmetic()):
                    return (self.rat().subs(d=self.group().dvalue()) == other.rat().subs(d=other.group().dvalue()))
                else:
                    return (self.rat() == other.rat())
            else:
                return False
        else:
            return False

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: QuasiModularFormsRing(n=5)(x^3*z-d*y)
            f_rho^3*E2 - f_i*d
        """

        return self._rat_repr()

    def _rat_repr(self):
        r"""
        Return a string representation of ``self`` as a rational function in the generators.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: QuasiModularForms(n=5, k=6, ep=-1)(x^3*z)._rat_repr()
            'f_rho^3*E2'
        """

        #return "{} in {}".format(str(self._rat), self.parent())
        with localvars(self.parent()._pol_ring, "f_rho, f_i, E2, d"):
            pol_str = str(self._rat)
        return "{}".format(pol_str)

    def _qexp_repr(self):
        r"""
        Return a string representation of ``self`` as a Fourier series.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: MR = QuasiModularFormsRing(n=5)
            sage: MR.disp_prec(3)
            sage: MR(x^3*z-d*y)._qexp_repr()
            '-d + 1 + ((65*d + 33)/(200*d))*q + ((1755*d + 1437)/(320000*d^2))*q^2 + O(q^3)'
        """

        n=self.hecke_n()

        # For now the series constructor doesn't behave well for non exact bases... :(
        if (self.group().is_arithmetic() or not self.base_ring().is_exact()):
            return str(self.q_expansion_fixed_d().add_bigoh(self.parent()._disp_prec))
        else:
            return str(self.q_expansion().add_bigoh(self.parent()._disp_prec))

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: (x,y,z,d)=var("x,y,z,d")
            sage: latex(QuasiModularFormsRing(n=5)(x^3*z-d*y))
            f_{\rho}^{3} E_{2} -  f_{i} d

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: latex(CuspForms(k=12)(x^3-y^2))
            f_{\rho}^{3} -  f_{i}^{2}
        """
        from sage.misc.latex import latex
        with localvars(self.parent()._pol_ring, "f_rho, f_i, E2, d"):
            latex_str = latex(self._rat)
        return latex_str

    def group(self):
        r"""
        Return the (Hecke triangle) group of ``self.parent()``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ModularForms(n=12, k=4).E4().group()
            Hecke triangle group for n = 12
        """

        return self.parent().group()

    def hecke_n(self):
        r"""
        Return the parameter ``n`` of the (Hecke triangle) group of ``self.parent()``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ModularForms(n=12, k=6).E6().hecke_n()
            12
        """

        return self.parent().hecke_n()

    def base_ring(self):
        r"""
        Return base ring of ``self.parent()``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ModularForms(n=12, k=4, base_ring=CC).E4().base_ring()
            Complex Field with 53 bits of precision
        """

        return self.parent().base_ring()

    def coeff_ring(self):
        r"""
        Return coefficient ring of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing().E6().coeff_ring()
            Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        return self.parent().coeff_ring()

    def rat(self):
        r"""
        Return the rational function representing ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing(n=12).Delta().rat()
            x^30*d - x^18*y^2*d
        """

        return self._rat

    def is_homogeneous(self):
        r"""
        Return whether ``self`` is homogeneous.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: QuasiModularFormsRing(n=12).Delta().is_homogeneous()
            True
            sage: QuasiModularFormsRing(n=12).Delta().parent().is_homogeneous()
            False
            sage: x,y,z,d=var("x,y,z,d")
            sage: QuasiModularFormsRing(n=12)(x^3+y^2+z+d).is_homogeneous()
            False
        """

        return self._weight != None

    def weight(self):
        r"""
        Return the weight of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiModularFormsRing()(x+y).weight() is None
            True
            sage: ModularForms(n=18).f_i().weight()
            9/4
        """

        return self._weight

    def ep(self):
        r"""
        Return the multiplier of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiModularFormsRing()(x+y).ep() is None
            True
            sage: ModularForms(n=18).f_i().ep()
            -1
        """

        return self._ep

    def degree(self):
        r"""
        Return the degree of ``self`` in the graded ring.
        If ``self`` is not homogeneous, then ``(None, None)``
        is returned.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiModularFormsRing()(x+y).degree() == (None, None)
            True
            sage: ModularForms(n=18).f_i().degree()
            (9/4, -1)
        """

        return (self._weight,self._ep)

    def is_modular(self):
        r"""
        Return whether ``self`` (resp. its homogeneous components)
        transform like modular forms.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiModularFormsRing(n=5)(x^2+y-d).is_modular()
            True
            sage: QuasiModularFormsRing(n=5)(x^2+y-d+z).is_modular()
            False
            sage: QuasiModularForms(n=18).f_i().is_modular()
            True
            sage: QuasiModularForms(n=18).E2().is_modular()
            False
        """

        return not (self.AT("quasi") <= self._analytic_type)

    def is_weakly_holomorphic(self):
        r"""
        Return whether ``self`` is weakly holomorphic
        in the sense that: ``self`` has at most a power of ``f_inf``
        in its denominator.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiMeromorphicModularFormsRing(n=5)(x/(x^5-y^2)+z).is_weakly_holomorphic()
            True
            sage: QuasiMeromorphicModularFormsRing(n=5)(x^2+y/x-d).is_weakly_holomorphic()
            False
            sage: QuasiMeromorphicModularForms(n=18).J_inv().is_weakly_holomorphic()
            True
        """

        return self.AT("weak", "quasi") >= self._analytic_type

    def is_holomorphic(self):
        r"""
        Return whether ``self`` is holomorphic
        in the sense that the denominator of ``self``
        is constant.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiMeromorphicModularFormsRing(n=5)((y^3-z^5)/(x^5-y^2)+y-d).is_holomorphic()
            False
            sage: QuasiMeromorphicModularFormsRing(n=5)(x^2+y-d+z).is_holomorphic()
            True
            sage: QuasiMeromorphicModularForms(n=18).J_inv().is_holomorphic()
            False
            sage: QuasiMeromorphicModularForms(n=18).f_i().is_holomorphic()
            True
        """

        return self.AT("holo", "quasi") >= self._analytic_type

    def is_cuspidal(self):
        r"""
        Return whether ``self`` is cuspidal
        in the sense that ``self`` is holomorphic and ``f_inf``
        divides the numerator.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiModularFormsRing(n=5)(y^3-z^5).is_cuspidal()
            False
            sage: QuasiModularFormsRing(n=5)(z*x^5-z*y^2).is_cuspidal()
            True
            sage: QuasiModularForms(n=18).Delta().is_cuspidal()
            True
            sage: QuasiModularForms(n=18).f_rho().is_cuspidal()
            False
        """

        return self.AT("cusp", "quasi") >= self._analytic_type

    def is_zero(self):
        r"""
        Return whether ``self`` is the zero function.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiModularFormsRing(n=5)(1).is_zero()
            False
            sage: QuasiModularFormsRing(n=5)(0).is_zero()
            True
            sage: QuasiModularForms(n=18).zero().is_zero()
            True
            sage: QuasiModularForms(n=18).Delta().is_zero()
            False
        """

        return self.AT(["quasi"]) >= self._analytic_type

    def analytic_type(self):
        r"""
        Return the analytic type of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiMeromorphicModularFormsRing(n=5)(x/z+d).analytic_type()
            quasi meromorphic modular
            sage: QuasiMeromorphicModularFormsRing(n=5)((y^3-z^5)/(x^5-y^2)+y-d).analytic_type()
            quasi weakly holomorphic modular
            sage: QuasiMeromorphicModularFormsRing(n=5)(x^2+y-d).analytic_type()
            modular
            sage: QuasiMeromorphicModularForms(n=18).J_inv().analytic_type()
            weakly holomorphic modular
            sage: QuasiMeromorphicModularForms(n=18).Delta().analytic_type()
            cuspidal
        """

        return self._analytic_type

    def numerator(self):
        r"""
        Return the numerator of ``self``.
        I.e. the (properly reduced) new form corresponding to
        the numerator of ``self.rat()``.
        
        Note that the parent of ``self`` might (probably will) change.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiMeromorphicModularFormsRing(n=5)((y^3-z^5)/(x^5-y^2)+y-d).numerator()
            f_rho^5*f_i - f_rho^5*d - E2^5 + f_i^2*d
            sage: QuasiMeromorphicModularFormsRing(n=5)((y^3-z^5)/(x^5-y^2)+y-d).numerator().parent()
            QuasiModularFormsRing(n=5) over Integer Ring
            sage: QuasiMeromorphicModularForms(n=5, k=-2, ep=-1)(x/y).numerator()
            1 + 7/(100*d)*q + 21/(160000*d^2)*q^2 + 1043/(192000000*d^3)*q^3 + 45479/(1228800000000*d^4)*q^4 + O(q^5)
            sage: QuasiMeromorphicModularForms(n=5, k=-2, ep=-1)(x/y).numerator().parent()
            QuasiModularForms(n=5, k=4/3, ep=1) over Integer Ring
        """

        res = self.parent().rat_field()(self._rat.numerator())
        # In general the numerator has a different weight than the original function...
        new_parent = self.parent().extend_type(ring=True).reduce_type(["holo", "quasi"])
        return new_parent(res).reduce()

    def denominator(self):
        r"""
        Return the denominator of ``self``.
        I.e. the (properly reduced) new form corresponding to
        the numerator of ``self.rat()``.
        
        Note that the parent of ``self`` might (probably will) change.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: x,y,z,d = var("x,y,z,d")
            sage: QuasiMeromorphicModularFormsRing(n=5).Delta().full_reduce().denominator()
            1 + O(q^5)
            sage: QuasiMeromorphicModularFormsRing(n=5)((y^3-z^5)/(x^5-y^2)+y-d).denominator()
            f_rho^5 - f_i^2
            sage: QuasiMeromorphicModularFormsRing(n=5)((y^3-z^5)/(x^5-y^2)+y-d).denominator().parent()
            QuasiModularFormsRing(n=5) over Integer Ring
            sage: QuasiMeromorphicModularForms(n=5, k=-2, ep=-1)(x/y).denominator()
            1 - 13/(40*d)*q - 351/(64000*d^2)*q^2 - 13819/(76800000*d^3)*q^3 - 1163669/(491520000000*d^4)*q^4 + O(q^5)
            sage: QuasiMeromorphicModularForms(n=5, k=-2, ep=-1)(x/y).denominator().parent()
            QuasiModularForms(n=5, k=10/3, ep=-1) over Integer Ring
        """

        res = self.parent().rat_field()(self._rat.denominator())
        # In general the denominator has a different weight than the original function...
        new_parent = self.parent().extend_type("holo", ring=True).reduce_type(["holo", "quasi"])
        return new_parent(res).reduce()

    def _add_(self,other):
        r"""
        Return the sum of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: MR = QuasiMeromorphicModularFormsRing(n=7)
            sage: E2 = MR.E2().full_reduce()
            sage: E4 = MR.E4().full_reduce()
            sage: E6 = MR.E6().full_reduce()
            sage: Delta = MR.Delta().full_reduce()
            sage: J_inv = MR.J_inv().full_reduce()
            sage: ring_el = MR(1/x+1).full_reduce()

            sage: (Delta^2*E2 + E6*E4^2).parent()
            QuasiModularFormsRing(n=7) over Integer Ring
            sage: E4 + Delta
            f_rho^15*d - f_rho^8*f_i^2*d + f_rho^5
            sage: (E4 + QQ(1) + ring_el).parent()
            MeromorphicModularFormsRing(n=7) over Integer Ring
            sage: (E4^3 + 1.1*Delta).parent()
            ModularForms(n=7, k=12, ep=1) over Real Field with 53 bits of precision
            sage: (E4 + FractionField(PolynomialRing(CC,'d')).gen()).parent()
            ModularFormsRing(n=7) over Complex Field with 53 bits of precision            

            sage: subspace = MR.reduce_type(["holo"], degree=(12,1)).subspace([Delta, E6^2])
            sage: gen0 = subspace.gen(0)
            sage: gen1 = subspace.gen(1)
            sage: subspace2 = MR.reduce_type(["holo"], degree=(12,1)).subspace([Delta, Delta + E6^2])
            sage: gen2 = subspace2.gen(0)
            sage: gen3 = subspace2.gen(1)

            sage: (gen0 + gen1).parent()
            Subspace of dimension 2 of ModularForms(n=7, k=12, ep=1) over Integer Ring
            sage: (gen0 + Delta*J_inv).parent()
            WeakModularForms(n=7, k=12, ep=1) over Integer Ring
            sage: gen0 + E2
            f_rho^15*d - f_rho^8*f_i^2*d + E2
            sage: (gen0 + E2).parent()
            QuasiModularFormsRing(n=7) over Integer Ring
            sage: gen2 + ring_el
            (f_rho^16*d - f_rho^9*f_i^2*d + f_rho + 1)/f_rho
            sage: (gen0 + int(1)).parent()
            ModularFormsRing(n=7) over Integer Ring
        """
        
        return self.parent()(self._rat+other._rat)

    def _sub_(self,other):
        r"""
        Return the difference of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: MR = QuasiMeromorphicModularFormsRing(n=7)
            sage: E2 = MR.E2().full_reduce()
            sage: E4 = MR.E4().full_reduce()
            sage: E6 = MR.E6().full_reduce()
            sage: Delta = MR.Delta().full_reduce()
            sage: J_inv = MR.J_inv().full_reduce()
            sage: ring_el = MR(1/x+1).full_reduce()

            sage: E6^2-E4^3
            -1/d*q - 17/(56*d^2)*q^2 - 88887/(2458624*d^3)*q^3 - 941331/(481890304*d^4)*q^4 + O(q^5)
            sage: (E6^2-E4^3).parent()
            ModularForms(n=7, k=12, ep=1) over Integer Ring
            sage: (E4 - int(1)).parent()
            ModularFormsRing(n=7) over Integer Ring
            sage: E4 - FractionField(PolynomialRing(CC,'d')).gen()
            f_rho^5 - d
            sage: ((E4+E6)-E6).parent()
            ModularFormsRing(n=7) over Integer Ring

            sage: subspace = MR.reduce_type(["holo"], degree=(12,1)).subspace([Delta, E6^2])
            sage: gen0 = subspace.gen(0)
            sage: gen1 = subspace.gen(1)
            sage: subspace2 = MR.reduce_type(["holo"], degree=(12,1)).subspace([Delta, Delta + E6^2])
            sage: gen2 = subspace2.gen(0)
            sage: gen3 = subspace2.gen(1)

            sage: (gen0 - gen2).parent()
            Subspace of dimension 2 of ModularForms(n=7, k=12, ep=1) over Integer Ring
            sage: (gen2 - ring_el).parent()
            MeromorphicModularFormsRing(n=7) over Integer Ring
            sage: (gen0 - 1.1).parent()
            ModularFormsRing(n=7) over Real Field with 53 bits of precision
        """

        #reduce at the end? See example "sage: ((E4+E6)-E6).parent()"
        return self.parent()(self._rat-other._rat)

    def _mul_(self,other):
        r"""
        Return the product of ``self`` and ``other``.

        Note that the parent might (probably will) change.

        If ``parent.has_reduce_hom() == True`` then
        the result is reduced to be an element of
        the corresponding forms space if possible.
        
        In particular this is the case if both ``self``
        and ``other`` are (homogeneous) elements of a
        forms space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: MR = QuasiMeromorphicModularFormsRing(n=8)
            sage: E2 = MR.E2().full_reduce()
            sage: E4 = MR.E4().full_reduce()
            sage: E6 = MR.E6().full_reduce()
            sage: Delta = MR.Delta().full_reduce()
            sage: J_inv = MR.J_inv().full_reduce()
            sage: ring_el = MR(1/x+1).full_reduce()

            sage: (1*Delta).parent()
            CuspForms(n=8, k=12, ep=1) over Integer Ring
            sage: (1.1*Delta).parent()
            ModularForms(n=8, k=12, ep=1) over Real Field with 53 bits of precision
            sage: (E2*E4).parent()
            QuasiModularForms(n=8, k=6, ep=-1) over Integer Ring
            sage: E4^2
            1 + 15/(32*d)*q + 3255/(32768*d^2)*q^2 + 105445/(8388608*d^3)*q^3 + 36379615/(34359738368*d^4)*q^4 + O(q^5)
            sage: (E4^2).parent()
            ModularForms(n=8, k=8, ep=1) over Integer Ring
            sage: (1.1*E4).parent()
            ModularForms(n=8, k=4, ep=1) over Real Field with 53 bits of precision
            sage: (J_inv*ring_el).parent()
            MeromorphicModularFormsRing(n=8) over Integer Ring
            sage: (E4*FractionField(PolynomialRing(CC,'d')).gen()).parent()
            ModularForms(n=8, k=4, ep=1) over Complex Field with 53 bits of precision

            sage: subspace = MR.reduce_type(["holo"], degree=(12,1)).subspace([Delta, E6^2])
            sage: gen0 = subspace.gen(0)
            sage: gen1 = subspace.gen(1)
            sage: subspace2 = MR.reduce_type(["holo"], degree=(12,1)).subspace([Delta, Delta + E6^2])
            sage: gen2 = subspace2.gen(0)
            sage: gen3 = subspace2.gen(1)

            sage: (gen0 * gen1).parent()
            ModularForms(n=8, k=24, ep=1) over Integer Ring
            sage: (gen0 * Delta*J_inv).parent()
            WeakModularForms(n=8, k=24, ep=1) over Integer Ring
            sage: (gen0 * Delta*J_inv).reduced_parent()
            CuspForms(n=8, k=24, ep=1) over Integer Ring
            sage: (gen0 * E2).parent()
            QuasiModularForms(n=8, k=14, ep=-1) over Integer Ring
            sage: gen2 * ring_el
            f_rho^18*d + f_rho^17*d - f_rho^10*f_i^2*d - f_rho^9*f_i^2*d
            sage: (gen2 * ring_el).parent()
            MeromorphicModularFormsRing(n=8) over Integer Ring
            sage: (1.1*gen0).parent()
            ModularForms(n=8, k=12, ep=1) over Real Field with 53 bits of precision
        """

        res = self.parent().rat_field()(self._rat*other._rat)
        new_parent = self.parent().extend_type(ring=True)
        # The product of two homogeneous elements is homogeneous
        return new_parent(res).reduce()

    def _div_(self,other):
        r"""
        Return the division of ``self`` by ``other``.

        Note that the parent might (probably will) change.

        If ``parent.has_reduce_hom() == True`` then
        the result is reduced to be an element of
        the corresponding forms space if possible.
        
        In particular this is the case if both ``self``
        and ``other`` are (homogeneous) elements of a
        forms space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: MR = QuasiMeromorphicModularFormsRing(n=8)
            sage: E2 = MR.E2().full_reduce()
            sage: E4 = MR.E4().full_reduce()
            sage: E6 = MR.E6().full_reduce()
            sage: Delta = MR.Delta().full_reduce()
            sage: J_inv = MR.J_inv().full_reduce()
            sage: ring_el = MR(1/x+1).full_reduce()

            sage: (1/Delta).parent()
            MeromorphicModularForms(n=8, k=-12, ep=1) over Integer Ring
            sage: (1.1/Delta).parent()
            MeromorphicModularForms(n=8, k=-12, ep=1) over Real Field with 53 bits of precision
            sage: (Delta/Delta).parent()
            MeromorphicModularForms(n=8, k=0, ep=1) over Integer Ring
            sage: 1/E4
            1 - 15/(64*d)*q + 2145/(65536*d^2)*q^2 - 59545/(16777216*d^3)*q^3 + 22622585/(68719476736*d^4)*q^4 + O(q^5)
            sage: (ring_el/J_inv).parent()
            MeromorphicModularFormsRing(n=8) over Integer Ring
            sage: (FractionField(PolynomialRing(CC,'d')).gen()/E4).parent()
            MeromorphicModularForms(n=8, k=-4, ep=1) over Complex Field with 53 bits of precision

            sage: (E4^(-2)).parent()
            MeromorphicModularForms(n=8, k=-8, ep=1) over Integer Ring
            sage: ((E4.as_ring_element())^(-2)).parent() == (E4^(-2)).parent()
            True
            sage: (MR(x)^(-3)).parent()
            QuasiMeromorphicModularFormsRing(n=8) over Integer Ring

            sage: subspace = MR.reduce_type(["holo"], degree=(12,1)).subspace([Delta, E6^2])
            sage: gen0 = subspace.gen(0)
            sage: gen1 = subspace.gen(1)
            sage: subspace2 = MR.reduce_type(["holo"], degree=(12,1)).subspace([Delta, Delta + E6^2])
            sage: gen2 = subspace2.gen(0)
            sage: gen3 = subspace2.gen(1)

            sage: (gen0 / gen1).parent()
            MeromorphicModularForms(n=8, k=0, ep=1) over Integer Ring
            sage: (gen0 / Delta).parent()
            MeromorphicModularForms(n=8, k=0, ep=1) over Integer Ring
            sage: ring_el / gen2
            (f_rho + 1)/(f_rho^19*d - f_rho^11*f_i^2*d)
            sage: (ring_el / gen2).parent()
            MeromorphicModularFormsRing(n=8) over Integer Ring
            sage: (1.1/gen0).parent()
            MeromorphicModularForms(n=8, k=-12, ep=1) over Real Field with 53 bits of precision
        """

        res = self.parent().rat_field()(self._rat/other._rat)
        new_parent = self.parent().extend_type("mero", ring=True)
        # The division of two homogeneous elements is homogeneous
        return new_parent(res).reduce()

    def sqrt(self):
        r"""
        Try to return the square root of ``self``.
        I.e. the element corresponding to ``sqrt(self.rat())``.

        Whether this works or not depends on whether
        ``sqrt(self.rat())`` works and coerces into
        ``self.parent().rat_field()``.

        Note that the parent might (probably will) change.

        If ``parent.has_reduce_hom() == True`` then
        the result is reduced to be an element of
        the corresponding forms space if possible.
        
        In particular this is the case if ``self``
        is a (homogeneous) element of a forms space.

        .. TODO::

            Make square root in the underlying rational field work.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: E2=QuasiModularForms(k=2, ep=-1).E2()
            sage: sqrt(E2^2)
            Traceback (most recent call last):
            ...
            NotImplementedError: is_square() not implemented for elements of Multivariate Polynomial Ring in x, y, z, d over Integer Ring
        """

        res = self.parent().rat_field()(self._rat.sqrt())
        new_parent = self.parent().extend_type(ring=True)
        # The sqrt of a homogeneous element is homogeneous if it exists
        return self.parent()(res).reduce()

    def diff_op(self, op, new_parent=None):
        r"""
        Return the differential operator ``op`` applied to ``self``.
        If ``parent.has_reduce_hom() == True`` then the result
        is reduced to be an element of the corresponding forms
        space if possible.

        INPUT:

        - ``op``          -- An element of ``self.parent().diff_alg()``.
                             I.e. an element of the algebra over ``QQ``
                             of differential operators generated
                             by ``X, Y, Z, dX, dY, DZ``, where e.g. ``X``
                             corresponds to the multiplication by ``x``
                             (resp. ``f_rho``) and ``dX`` corresponds to ``d/dx``.

                             To expect a homogeneous result after applying
                             the operator to a homogeneous element it should
                             should be homogeneous operator (with respect
                             to the the usual, special grading).

        - ``new_parent``  -- Try to convert the result to the specified
                             ``new_parent``. If ``new_parent == None`` (default)
                             then the parent is extended to a
                             "quasi meromorphic" ring.
                
        OUTPUT:

        The new element.

        EXAMPLES::
            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: MR = QuasiMeromorphicModularFormsRing(n=8, red_hom=True)
            sage: (X,Y,Z,dX,dY,dZ) = MR.diff_alg().gens()
            sage: n=MR.hecke_n()
            sage: mul_op = 4/(n-2)*X*dX + 2*n/(n-2)*Y*dY + 2*Z*dZ
            sage: der_op = MR._derivative_op()
            sage: ser_op = MR._serre_derivative_op()
            sage: der_op == ser_op + (n-2)/(4*n)*Z*mul_op
            True
            
            sage: Delta = MR.Delta().full_reduce()
            sage: E2 = MR.E2().full_reduce()
            sage: Delta.diff_op(mul_op) == 12*Delta
            True
            sage: Delta.diff_op(mul_op).parent()
            QuasiMeromorphicModularForms(n=8, k=12, ep=1) over Integer Ring
            sage: Delta.diff_op(mul_op, Delta.parent()).parent()
            CuspForms(n=8, k=12, ep=1) over Integer Ring
            sage: E2.diff_op(mul_op, E2.parent()) == 2*E2
            True
            sage: Delta.diff_op(Z*mul_op, Delta.parent().extend_type("quasi", ring=True)) == 12*E2*Delta
            True

            sage: ran_op = X + Y*X*dY*dX + dZ + dX^2
            sage: Delta.diff_op(ran_op)
            f_rho^19*d + 306*f_rho^16*d - f_rho^11*f_i^2*d - 20*f_rho^10*f_i^2*d - 90*f_rho^8*f_i^2*d
            sage: E2.diff_op(ran_op)
            f_rho*E2 + 1
        """

        (x,y,z,d) = self.parent().rat_field().gens()
        (X,Y,Z,dX,dY,dZ) = self.parent().diff_alg().gens()
        L = op.monomials()
        new_rat = 0
        for mon in L:
            mon_summand  = self._rat
            mon_summand  = mon_summand.derivative(x,mon.degree(dX))
            mon_summand  = mon_summand.derivative(y,mon.degree(dY))
            mon_summand  = mon_summand.derivative(z,mon.degree(dZ))
            mon_summand *= x**(mon.degree(X))
            mon_summand *= y**(mon.degree(Y))
            mon_summand *= z**(mon.degree(Z))
            new_rat     += op.monomial_coefficient(mon)*mon_summand
        res = self.parent().rat_field()(new_rat)
        if (new_parent == None):
            new_parent = self.parent().extend_type(["quasi", "mero"], ring=True)
        return new_parent(res).reduce()

    # note that this is qd/dq, resp 1/(2*pi*i)*d/dtau
    def derivative(self):
        r"""
        Return the derivative ``d/dq = 1/(2*pi*i) d/dtau`` of ``self``.

        Note that the parent might (probably will) change.
        In particular its analytic type will be extended
        to contain "quasi".

        If ``parent.has_reduce_hom() == True`` then
        the result is reduced to be an element of
        the corresponding forms space if possible.
        
        In particular this is the case if ``self``
        is a (homogeneous) element of a forms space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: MR = QuasiMeromorphicModularFormsRing(n=7, red_hom=True)
            sage: n = MR.hecke_n()
            sage: E2 = MR.E2().full_reduce()
            sage: E6 = MR.E6().full_reduce()
            sage: f_rho = MR.f_rho().full_reduce()
            sage: f_i = MR.f_i().full_reduce()
            sage: f_inf = MR.f_inf().full_reduce()

            sage: derivative(f_rho) == 1/n * (f_rho*E2 - f_i)
            True
            sage: derivative(f_i)   == 1/2 * (f_i*E2 - f_rho**(n-1))
            True
            sage: derivative(f_inf) == f_inf * E2
            True
            sage: derivative(f_inf).parent()
            QuasiCuspForms(n=7, k=38/5, ep=-1) over Integer Ring
            sage: derivative(E2)    == (n-2)/(4*n) * (E2**2 - f_rho**(n-2))
            True
            sage: derivative(E2).parent()
            QuasiModularForms(n=7, k=4, ep=1) over Integer Ring
        """

        return self.diff_op(self.parent()._derivative_op(), self.parent().extend_type("quasi", ring=True))

    def serre_derivative(self):
        r"""
        Return the Serre derivative of ``self``.

        Note that the parent might (probably will) change.
        However a modular element is returned if ``self``
        was already modular.

        If ``parent.has_reduce_hom() == True`` then
        the result is reduced to be an element of
        the corresponding forms space if possible.
        
        In particular this is the case if ``self``
        is a (homogeneous) element of a forms space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: MR = QuasiMeromorphicModularFormsRing(n=7, red_hom=True)
            sage: n = MR.hecke_n()
            sage: Delta = MR.Delta().full_reduce()
            sage: E2 = MR.E2().full_reduce()
            sage: E4 = MR.E4().full_reduce()
            sage: E6 = MR.E6().full_reduce()
            sage: f_rho = MR.f_rho().full_reduce()
            sage: f_i = MR.f_i().full_reduce()
            sage: f_inf = MR.f_inf().full_reduce()

            sage: f_rho.serre_derivative() == -1/n * f_i
            True
            sage: f_i.serre_derivative()   == -1/2 * E4 * f_rho
            True
            sage: f_inf.serre_derivative() == 0
            True
            sage: E2.serre_derivative()    == -(n-2)/(4*n) * (E2^2 + E4)
            True
            sage: E4.serre_derivative()    == -(n-2)/n * E6
            True
            sage: E6.serre_derivative()    == -1/2 * E4^2 - (n-3)/n * E6^2 / E4
            True
            sage: E6.serre_derivative().parent()
            ModularForms(n=7, k=8, ep=1) over Integer Ring
        """

        return self.diff_op(self.parent()._serre_derivative_op(), self.parent().extend_type(ring=True))

    @cached_method
    def order_inf(self):
        """
        Return the order at infinity of ``self``.

        If ``self`` is homogeneous and modular then
        only the rational function ``self.rat()`` is used.
        Otherwise the Fourier expansion is used
        (with increasing precision until the order
        can be determined).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: MR = QuasiMeromorphicModularFormsRing(red_hom=True)
            sage: (MR.Delta()^3).order_inf()
            3
            sage: MR.E2().order_inf()
            0
            sage: (MR.J_inv()^2).order_inf()
            -2
            sage: x,y,z,d = MR.pol_ring().gens()
            sage: el = MR((z^3-y)^2/(x^3-y^2)).full_reduce()
            sage: el
            108*q + 11664*q^2 + 502848*q^3 + 12010464*q^4 + O(q^5)
            sage: el.order_inf()
            1
            sage: el.parent()
            QuasiWeakModularForms(n=3, k=0, ep=1) over Integer Ring
            sage: el.is_holomorphic()
            False
            sage: MR((z-y)^2+(x-y)^3).order_inf()
            2
            sage: MR((x-y)^10).order_inf()
            10
            sage: MR.zero().order_inf()
            +Infinity
        """

        if self.is_zero():
            return infinity

        if (self.is_homogeneous() and self.is_modular()):
            rat       = self.parent().rat_field()(self._rat)
            R         = self.parent().pol_ring()
            numerator = R(rat.numerator())
            denom     = R(rat.denominator())
            (x,y,z,d) = R.gens()
            # We intentially leave out the d-factor!
            finf_pol  = (x**self.hecke_n()-y**2)
    
            order_inf = 0
            # There seems to be a bug in Singular, for now this "try, except" is a workaround
            # Also numerator /= finf_pol doesn't seem to return an element of R for non-exact rings...
            try:
                while (finf_pol.divides(numerator)):
                    numerator  = numerator.quo_rem(finf_pol)[0]
                    #numerator /= finf_pol
                    numerator  = R(numerator)
                    order_inf += 1
            except TypeError:
                pass
            try:
                while (finf_pol.divides(denom)):
                    denom      = denom.quo_rem(finf_pol)[0]
                    #denom     /= finf_pol
                    denom      = R(denom)
                    order_inf -= 1
            except TypeError:
                pass

            return order_inf

        else:
            num_val   = prec_num_bound   = 1 #(self.parent()._prec/ZZ(2)).ceil()
            denom_val = prec_denom_bound = 1 #(self.parent()._prec/ZZ(2)).ceil()

            while (num_val >= prec_num_bound):
                prec_num_bound   *= 2
                num_val   = self.numerator().q_expansion(prec=prec_num_bound, fix_prec=True).valuation()
            while (denom_val >= prec_denom_bound):
                prec_denom_bound *= 2
                denom_val = self.denominator().q_expansion(prec=prec_denom_bound, fix_prec=True).valuation()

            return (num_val-denom_val)

    def as_ring_element(self):
        r"""
        Coerce ``self`` into the graded ring of its parent.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: Delta = CuspForms(k=12).Delta()
            sage: Delta.parent()
            CuspForms(n=3, k=12, ep=1) over Integer Ring
            sage: Delta.as_ring_element()
            f_rho^3*d - f_i^2*d
            sage: Delta.as_ring_element().parent()
            CuspFormsRing(n=3) over Integer Ring
        """

        return self.parent().graded_ring()(self._rat)

    def reduce(self, force=False):
        r"""
        In case ``self.parent().has_reduce_hom() == True``
        (or ``force==True``) and ``self`` is homogeneous
        the converted element lying in the corresponding
        homogeneous_space is returned.
        
        Otherwise ``self`` is returned.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: E2 = ModularFormsRing(n=7).E2().reduce()
            sage: E2.parent()
            QuasiModularFormsRing(n=7) over Integer Ring
            sage: E2 = ModularFormsRing(n=7, red_hom=True).E2().reduce()
            sage: E2.parent()
            QuasiModularForms(n=7, k=2, ep=-1) over Integer Ring
            sage: ModularFormsRing(n=7)(x+1).reduce().parent()
            ModularFormsRing(n=7) over Integer Ring
            sage: E2 = ModularFormsRing(n=7).E2().reduce(force=True)
            sage: E2.parent()
            QuasiModularForms(n=7, k=2, ep=-1) over Integer Ring
            sage: ModularFormsRing(n=7)(x+1).reduce(force=True).parent()
            ModularFormsRing(n=7) over Integer Ring
        """

        if (force or self.parent().has_reduce_hom()) and self.is_homogeneous():
            return self.parent().homogeneous_space(self._weight, self._ep)(self._rat)
        else:
            return self

    def reduced_parent(self):
        r"""
        Return the space with the analytic type of ``self``.
        If ``self`` is homogeneous the corresponding ``FormsSpace`` is returned.
        
        I.e. return the smallest known ambient space of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: Delta = QuasiMeromorphicModularFormsRing(n=7).Delta()
            sage: Delta.parent()
            QuasiMeromorphicModularFormsRing(n=7) over Integer Ring
            sage: Delta.reduced_parent()
            CuspForms(n=7, k=12, ep=1) over Integer Ring
            sage: el = QuasiMeromorphicModularFormsRing()(x+1)
            sage: el.parent()
            QuasiMeromorphicModularFormsRing(n=3) over Integer Ring
            sage: el.reduced_parent()
            ModularFormsRing(n=3) over Integer Ring
        """

        if self.is_homogeneous():
            return FormsSpace(self.analytic_type(), self.group(), self.base_ring(), self.weight(), self.ep())
        else:
            return FormsRing(self.analytic_type(), self.group(), self.base_ring(), self.parent().has_reduce_hom())

    def full_reduce(self):
        r"""
        Convert ``self`` into its reduced parent.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: Delta = QuasiMeromorphicModularFormsRing().Delta()
            sage: Delta
            f_rho^3*d - f_i^2*d
            sage: Delta.full_reduce()
            q - 24*q^2 + 252*q^3 - 1472*q^4 + O(q^5)
            sage: Delta.full_reduce().parent() == Delta.reduced_parent()
            True
        """

        return self.reduced_parent()(self._rat)

    #precision is actually acuracy, maybe add "real precision" meaning number of rel. coef
    @cached_method
    def _q_expansion_cached(self, prec, fix_d, d, d_num_prec, fix_prec = False):
        """
        Returns the Fourier expansion of self (cached).
        Don't call this function, instead use ``q_expansion``.

        INPUT:

        - ``prec``       -- An integer, the desired output precision O(q^prec).

        - ``fix_d``      -- If ``True`` then the numerical value of d
                            corresponding to n will be used.
                            If n = 3, 4, 6 the used value is exact.
                            The base_ring will be changed accordingly (if possible).
                            Also see ``MFSeriesConstructor``.

        - ``d``         -- ``None`` or a value to substitute for d.
                            The base_ring will be changed accordingly (if possible).
                            Also see ``MFSeriesConstructor``.

        - ``d_num_prec`` -- An integer, namely the precision to be used if a
                            numerical value for d is substituted.

        - ``fix_prec``   -- If ``fix_prec`` is not ``False`` (default)
                            then the precision of the ``MFSeriesConstructor`` is
                            set such that the output has exactly the specified
                            precision O(q^prec).

        OUTPUT:

        The Fourier expansion of self as a ``FormalPowerSeries`` or ``FormalLaurentSeries``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import WeakModularFormsRing
            sage: J_inv = WeakModularFormsRing(red_hom=True).J_inv()
            sage: J_inv._q_expansion_cached(5, False, None, 53, False) == J_inv.q_expansion()
            True
            sage: J_inv._q_expansion_cached(5, True, None, 53, False) == J_inv.q_expansion_fixed_d()
            True
        """

        if (fix_prec == False):
            #if (prec <1):
            #    print "Warning: non-positiv precision!"
            if not (fix_d or (d is not None) or self.base_ring().is_exact()):
                from warnings import warn
                warn("For non-exact base rings it is strongly recommended to fix/set d!")
            if ((not self.is_zero()) and prec <= self.order_inf()):
                from warnings import warn
                warn("precision too low to determine any coefficient!")

            # This should _exactly_ ensure the given precision O(q^prec):
            prec += self.denominator().order_inf() + max(-self.order_inf(), 0)

            # The result will have "max(prec-self.order_inf(),0)" significant coefficients
            # So adding the following line to the above one
            # would instead give "prec" many significant coefficients:
            # prec += self.order_inf()

        SC = MFSeriesConstructor(self.group(), prec)
        (base_ring, coeff_ring, qseries_ring, D) = SC.series_data(self.base_ring(), fix_d, d, d_num_prec)
        X  = SC.f_rho(self.base_ring(), fix_d, d, d_num_prec)
        Y  = SC.f_i(self.base_ring(), fix_d, d, d_num_prec)

        if (self.parent().is_modular()):
            qexp = self._rat.subs(x=X,y=Y,d=D)
        else:
            Z = SC.E2(self.base_ring(), fix_d, d, d_num_prec)
            qexp = self._rat.subs(x=X,y=Y,z=Z,d=D)
 
        return (qexp + O(qseries_ring.gen()**prec)).parent()(qexp)

    def q_expansion(self, prec = None, fix_d = False, d=None, d_num_prec = None, fix_prec = False):
        """
        Returns the Fourier expansion of self.

        INPUT:

        - ``prec``       -- An integer, the desired output precision O(q^prec).
                            Default: ``None`` in which case the default precision
                            of ``self.parent()`` is used.

        - ``fix_d``      -- If ``True`` then the numerical value of d corresponding to n
                            will be used (default: ``False``).
                            If n = 3, 4, 6 the used value is exact.
                            The base_ring will be changed accordingly (if possible).
                            Alternatively also a value as in ``d`` can be specified.
                            Also see ``MFSeriesConstructor``.

        - ``d``      -- ``None`` (default) or a value to substitute for d.
                            The base_ring will be changed accordingly (if possible).
                            Also see ``MFSeriesConstructor``.

        - ``d_num_prec`` -- The precision to be used if a numerical value for d is substituted.
                            Default: ``None`` in which case the default
                            numerical precision of ``self.parent()`` is used.

        - ``fix_prec``   -- If ``fix_prec`` is not ``False`` (default)
                            then the precision of the ``MFSeriesConstructor`` is
                            increased such that the output has exactly the specified
                            precision O(q^prec).

        OUTPUT:

        The Fourier expansion of self as a ``FormalPowerSeries`` or ``FormalLaurentSeries``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import WeakModularFormsRing, QuasiModularFormsRing
            sage: j_inv = WeakModularFormsRing(red_hom=True).j_inv()
            sage: j_inv.q_expansion(prec=3)
            q^-1 + 31/(72*d) + 1823/(27648*d^2)*q + 10495/(2519424*d^3)*q^2 + O(q^3)

            sage: E2 = QuasiModularFormsRing(n=5, red_hom=True).E2()
            sage: E2.q_expansion(prec=3)
            1 - 9/(200*d)*q - 369/(320000*d^2)*q^2 + O(q^3)
            sage: E2.q_expansion(prec=3, fix_d=1)
            1 - 9/200*q - 369/320000*q^2 + O(q^3)

            sage: E6 = WeakModularFormsRing(n=5, red_hom=True).E6().full_reduce()
            sage: Delta = WeakModularFormsRing(n=5, red_hom=True).Delta().full_reduce()
            sage: E6.q_expansion(prec=3).prec() == 3
            True
            sage: (Delta/(E2^3-E6)).q_expansion(prec=3).prec() == 3
            True
            sage: (Delta/(E2^3-E6)^3).q_expansion(prec=3).prec() == 3
            True
            sage: ((E2^3-E6)/Delta^2).q_expansion(prec=3).prec() == 3
            True
            sage: ((E2^3-E6)^3/Delta).q_expansion(prec=3).prec() == 3
            True

            sage: x,y = var("x,y")
            sage: el = WeakModularFormsRing()((x+1)/(x^3-y^2))
            sage: el.q_expansion(prec=2, fix_prec = True)
            2*d*q^-1 + O(1)
            sage: el.q_expansion(prec=2)
            2*d*q^-1 + 1/6 + 119/(41472*d)*q + O(q^2)
        """

        if prec == None:
            prec = self.parent()._prec
        if d_num_prec == None:
            d_num_prec = self.parent()._num_prec
        
        if isinstance(fix_d,bool) and (fix_d == True):
            d = None
        elif d is None:
            d = fix_d
            fix_d = False
        else:
            fix_d = False

        return self._q_expansion_cached(prec, fix_d, d, d_num_prec, fix_prec)

    def q_expansion_fixed_d(self, prec = None, d_num_prec = None, fix_prec = False):
        """
        Returns the Fourier expansion of self.
        The numerical (or exact) value for d is substituted.


        INPUT:

        - ``prec``       -- An integer, the desired output precision O(q^prec).
                            Default: ``None`` in which case the default precision
                            of ``self.parent()`` is used.

        - ``d_num_prec`` -- The precision to be used if a numerical value for d is substituted.
                            Default: ``None`` in which case the default
                            numerical precision of ``self.parent()`` is used.

        - ``fix_prec``   -- If ``fix_prec`` is not ``False`` (default)
                            then the precision of the ``MFSeriesConstructor`` is
                            increased such that the output has exactly the specified
                            precision O(q^prec).

        OUTPUT:

        The Fourier expansion of self as a ``FormalPowerSeries`` or ``FormalLaurentSeries``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import WeakModularFormsRing, QuasiModularFormsRing
            sage: j_inv = WeakModularFormsRing(red_hom=True).j_inv()
            sage: j_inv.q_expansion_fixed_d(prec=3)
            q^-1 + 744 + 196884*q + 21493760*q^2 + O(q^3)

            sage: E2 = QuasiModularFormsRing(n=5, red_hom=True).E2()
            sage: E2.q_expansion_fixed_d(prec=3)
            1.000000000000... - 6.380956565426...*q - 23.18584547617...*q^2 + O(q^3)

            sage: x,y = var("x,y")
            sage: WeakModularFormsRing()((x+1)/(x^3-y^2)).q_expansion_fixed_d(prec=2, fix_prec = True)
            1/864*q^-1 + O(1)
            sage: WeakModularFormsRing()((x+1)/(x^3-y^2)).q_expansion_fixed_d(prec=2)
            1/864*q^-1 + 1/6 + 119/24*q + O(q^2)
        """
        return self.q_expansion(prec, True, None, d_num_prec, fix_prec)

    def __call__(self, tau, prec = None, num_prec = None):
        r"""
        Try too return ``self`` evaluated at a point ``tau``
        in the upper half plane, where ``self`` is interpreted
        as a function in ``tau``, where ``q=exp(2*pi*i*tau)``.
        
        Note that this interpretation might not make sense
        (and fail) for certain (many) choices of
        (``base_ring``, ``tau.parent()``).

        To obtain a precise and fast result the parameters
        ``prec`` and ``num_prec`` both have to be considered/balanced.
        A high ``prec`` value is usually quite costly.

        INPUT:
        
        - ``tau``       -- ``infinity`` or an element of the upper
                           half plane. E.g. with parent ``AA`` or ``CC``.

        - ``prec``      -- An integer, namely the precision used for the
                           Fourier expansion. If ``prec == None`` (default)
                           then the default precision of ``self.parent()``
                           is used.
                      
        - ``num_prec``  -- An integer, namely the minimal numerical precision
                           used for ``tau`` and ``d``. If `num_prec == None`
                           (default) then the default numerical precision of
                           ``self.parent()`` is used.

        OUTPUT:

        The (numerical) evaluated function value.


        ALGORITHM:

        #. If ``(tau == infinity)``:
        
           Return ``infinity`` unless the order at infinity is
           non-negative in which case the constant Fourier coefficient
           is returned.

        #. Else if ``self`` is homogeneous and modular:
        
           #. Determine the matrix ``A`` which sends ``tau``
              to ``w`` in the fundamental domain, together
              with ``aut_factor(A,w)``.
           
              Note: These values are determined by the method
              ``get_FD(tau, aut_factor)`` from ``self.group()``,
              where ``aut_factor = self.parent().aut_factor``
              (which is only defined on the basic generators).

           #. Because of the (modular) transformation property
              of ``self`` the evaluation at ``tau`` is given by
              the evaluation at ``w`` multiplied by ``aut_factor(A,w)``.

           #. The evaluation at ``w`` is calculated by
              evaluating the truncated Fourier expansion of
              self at ``q(w)``.
           
           Note that this is much faster and more precise
           than a direct evaluation at ``tau``.

        #. Else if ``self`` is exactly ``E2``:

           #. The same procedure as before is applied (with
              the aut_factor from the corresponding modular
              space).

           #. Except that at the end a correction term for
              the quasimodular form ``E2`` of the form 
              ``4*lambda/(2*pi*i)*n/(n-2) * c*(c*w + d)``
              has to be added, where ``lambda = 2*cos(pi/n)``
              and ``c,d`` are the lower entries of the matrix
              ``A``.

        #. Else:

           #. Evaluate ``f_rho, f_i, E2`` at ``tau``
              using the above procedures.

           #. Substitute ``x=f_rho(tau), y=f_i(tau), z=E2(tau)``
              and the numerical value of ``d`` for ``d``
              in ``self.rat()``.


        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: MR = QuasiMeromorphicModularFormsRing(n=5, red_hom=True)
            sage: f_rho = MR.f_rho().full_reduce()
            sage: f_i   = MR.f_i().full_reduce()
            sage: f_inf = MR.f_inf().full_reduce()
            sage: E2    = MR.E2().full_reduce()
            sage: E4    = MR.E4().full_reduce()
            sage: rho   = MR.group().rho()

            sage: f_rho(rho)
            7.924774170220...e-10 + 3.216846821336...e-16*I
            sage: f_i(i)
            6.072919944699...e-14
            sage: f_inf(infinity)
            0

            sage: z = -1/(-1/(2*i+30)-1)
            sage: z
            2/965*I + 934/965
            sage: E4(z)
            32288.05588811... - 118329.8566016...*I
            sage: E4(z, prec=30, num_prec=100)    # long time
            32288.0558872351130041311053... - 118329.856600349999751420381...*I
            sage: E2(z)
            409.3144737105... + 100.6926857489...*I
            sage: E2(z, prec=30, num_prec=100)    # long time
            409.314473710489761254584951... + 100.692685748952440684513866...*I
            sage: (E2^2-E4)(z)
            125111.2655383... + 200759.8039479...*I
            sage: (E2^2-E4)(z, prec=30, num_prec=100)    # long time
            125111.265538336196262200469... + 200759.803948009905410385699...*I

            sage: (E2^2-E4)(infinity)
            0
            sage: (1/(E2^2-E4))(infinity)
            +Infinity
            sage: ((E2^2-E4)/f_inf)(infinity)
            -3/(10*d)

            sage: G = HeckeTriangleGroup(n=8)
            sage: MR = QuasiMeromorphicModularFormsRing(group=G, red_hom=True)
            sage: f_rho = MR.f_rho().full_reduce()
            sage: f_i   = MR.f_i().full_reduce()
            sage: E2    = MR.E2().full_reduce()

            sage: z = AlgebraicField()(1/10+13/10*I)
            sage: A = G.V(4)
            sage: S = G.S()
            sage: T = G.T()
            sage: A == (T*S)**3*T
            True
            sage: az = G.act(A, z)
            sage: az == (A[0,0]*z + A[0,1]) / (A[1,0]*z + A[1,1])
            True

            sage: f_rho(z)
            1.03740476727... + 0.0131941034523...*I
            sage: f_rho(az)
            -2.29216470688... - 1.46235057536...*I
            sage: k = f_rho.weight()
            sage: aut_fact = f_rho.ep()^3 * (G.act((T*S)**2*T, z)/i)**k * (G.act((T*S)*T, z)/i)**k * (G.act(T, z)/i)**k
            sage: aut_fact * f_rho(z)
            -2.29216470688... - 1.46235057536...*I

            sage: f_rho.parent().default_num_prec(1000)
            sage: f_rho.parent().default_prec(300)
            sage: (f_rho.q_expansion_fixed_d().polynomial())(exp(2*pi*i*(z/G.lam()).n(1000))).n(100)    # long time
            1.0374047672719462149821251... + 0.013194103452368974597290332...*I
            sage: (f_rho.q_expansion_fixed_d().polynomial())(exp(2*pi*i*(az/G.lam()).n(1000))).n(100)    # long time
            -2.2921647068881834598616367... - 1.4623505753697635207183406...*I

            sage: f_i(z)
            0.667489320423... - 0.118902824870...*I
            sage: f_i(az)
            14.5845388476... - 28.4604652892...*I
            sage: k = f_i.weight()
            sage: aut_fact = f_i.ep()^3 * (G.act((T*S)**2*T, z)/i)**k * (G.act((T*S)*T, z)/i)**k * (G.act(T, z)/i)**k
            sage: aut_fact * f_i(z)
            14.5845388476... - 28.4604652892...*I

            sage: f_i.parent().default_num_prec(1000)
            sage: f_i.parent().default_prec(300)
            sage: (f_i.q_expansion_fixed_d().polynomial())(exp(2*pi*i*(z/G.lam()).n(1000))).n(100)    # long time
            0.66748932042300250077433252... - 0.11890282487028677063054267...*I
            sage: (f_i.q_expansion_fixed_d().polynomial())(exp(2*pi*i*(az/G.lam()).n(1000))).n(100)    # long time
            14.584538847698600875918891... - 28.460465289220303834894855...*I

            sage: f = f_rho*E2
            sage: f(z)
            0.966024386418... - 0.0138894699429...*I
            sage: f(az)
            -15.9978074989... - 29.2775758341...*I
            sage: k = f.weight()
            sage: aut_fact = f.ep()^3 * (G.act((T*S)**2*T, z)/i)**k * (G.act((T*S)*T, z)/i)**k * (G.act(T, z)/i)**k
            sage: k2 = f_rho.weight()
            sage: aut_fact2 = f_rho.ep() * (G.act((T*S)**2*T, z)/i)**k2 * (G.act((T*S)*T, z)/i)**k2 * (G.act(T, z)/i)**k2
            sage: cor_term = (4 * G.n() / (G.n()-2) * A[1][0] * (A[1][0]*z+A[1][1])) / (2*pi*i).n(1000) * G.lam()
            sage: aut_fact*f(z) + cor_term*aut_fact2*f_rho(z)
            -15.9978074989... - 29.2775758341...*I

            sage: f.parent().default_num_prec(1000)
            sage: f.parent().default_prec(300)
            sage: (f.q_expansion_fixed_d().polynomial())(exp(2*pi*i*(z/G.lam()).n(1000))).n(100)    # long time
            0.96602438641867296777809436... - 0.013889469942995530807311503...*I
            sage: (f.q_expansion_fixed_d().polynomial())(exp(2*pi*i*(az/G.lam()).n(1000))).n(100)    # long time
            -15.997807498958825352887040... - 29.277575834123246063432206...*I
        """

        if (prec == None):
            prec = self.parent().default_prec()
        if (num_prec == None):
            num_prec = self.parent().default_num_prec()
    
        if (tau==infinity):
            if (self.AT(["cusp", "quasi"]) >= self.analytic_type()):
                return ZZ(0)
            #if (self.AT(["holo", "quasi"]) >= self.analytic_type()):
            elif self.order_inf() > 0:
                return ZZ(0)
            elif self.order_inf() < 0:
                return infinity
            else:
                return self.q_expansion(prec=1)[0]

        num_prec = max(\
            ZZ(getattr(tau,'prec',lambda: num_prec)()),\
            num_prec\
        )
        tau = tau.n(num_prec)
        (x,y,z,d) = self.parent().rat_field().gens()

        if (self.is_homogeneous() and self.is_modular()):
            q_exp = self.q_expansion_fixed_d(prec=prec, d_num_prec=num_prec)
            A, w, aut_factor = self.group().get_FD(tau,self.reduce(force=True).parent().aut_factor)
            if (type(q_exp) == LaurentSeries):
                return q_exp.laurent_polynomial()(exp((2 * pi * i).n(num_prec) / self.group().lam() * w)) * aut_factor
            else:
                return q_exp.polynomial()(exp((2 * pi * i).n(num_prec) / self.group().lam() * w)) * aut_factor
        elif (self._rat == z):
            E2 = self.parent().graded_ring().E2().reduce(force=True)
            A, w, aut_factor = self.group().get_FD(tau,E2.parent().aut_factor)
            E2_wvalue = E2.q_expansion_fixed_d(prec=prec, d_num_prec=num_prec).polynomial()(exp((2 * pi * i).n(num_prec) / self.group().lam() * w))
            E2_cor_term = 4 * self.group().lam() / (2*pi*i).n(num_prec) * self.hecke_n() / (self.hecke_n()-2) * A[1][0] * (A[1][0]*w + A[1][1])
            return E2_wvalue*aut_factor + E2_cor_term
        else:
            f_rho = self.parent().graded_ring().f_rho()
            f_i   = self.parent().graded_ring().f_i()
            E2    = self.parent().graded_ring().E2()
            dval  = self.parent().group().dvalue().n(num_prec)
            return self._rat.subs(x=f_rho(tau), y=f_i(tau), z=E2(tau), d=dval)
