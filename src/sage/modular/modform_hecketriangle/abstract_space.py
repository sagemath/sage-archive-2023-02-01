r"""
Modular forms for Hecke triangle groups

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

from sage.symbolic.all import i
from sage.rings.all import ZZ, QQ, infinity, AlgebraicField
from sage.rings.power_series_ring import is_PowerSeriesRing
from sage.rings.laurent_series_ring import is_LaurentSeriesRing
from sage.modules.free_module_element import is_FreeModuleElement

from sage.misc.cachefunc import cached_method

from abstract_ring import FormsRing_abstract
from series_constructor import MFSeriesConstructor


class FormsSpace_abstract(FormsRing_abstract):
    r"""
    Abstract (Hecke) forms space.

    This should never be called directly. Instead one should
    instantiate one of the derived classes of this class. 
    """

    from element import FormsElement
    Element = FormsElement

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Abstract (Hecke) forms space.

        INPUT:

        - ``group``       -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``k``           -- The weight (default: `0`)

        - ``ep``          -- The epsilon (default: ``None``).
                             If ``None``, then k*(n-2) has to be divisible by `2` and
                             ``ep=(-1)^(k*(n-2)/2)`` is used.

        - ``base_ring``   -- The base_ring (default: `\Z`).

        OUTPUT:

        The corresponding abstract (Hecke) forms space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=5, base_ring=ZZ, k=6, ep=-1)
            sage: MF
            ModularForms(n=5, k=6, ep=-1) over Integer Ring
            sage: MF.group()
            Hecke triangle group for n = 5
            sage: MF.base_ring()
            Integer Ring
            sage: MF.weight()
            6
            sage: MF.ep()
            -1
            sage: MF.has_reduce_hom()
            True
            sage: MF.is_homogeneous()
            True
        """

        #from space import canonical_parameters
        #(group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)

        super(FormsSpace_abstract, self).__init__(group=group, base_ring=base_ring, red_hom=True, n=n)
        #self.register_embedding(self.hom(lambda f: f.parent().graded_ring()(f), codomain=self.graded_ring()))

        self._weight = k
        self._ep = ep
        (self._l1,self._l2) = self.weight_parameters()
        self._module = None
        self._ambient_space = self

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: QuasiModularForms(n=4, k=2, ep=-1)
            QuasiModularForms(n=4, k=2, ep=-1) over Integer Ring
        """

        return "{}Forms(n={}, k={}, ep={}) over {}".format(self._analytic_type.analytic_space_name(), self._group.n(), self._weight, self._ep, self._base_ring)

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiWeakModularForms
            sage: latex(QuasiWeakModularForms())
            QM^!_{ n=3 }(0,\ 1)(\Bold{Z})
        """

        from sage.misc.latex import latex
        return "{}_{{ n={} }}({},\ {})({})".format(self._analytic_type.latex_space_name(), self._group.n(), self._weight, self._ep, latex(self._base_ring))

    def _element_constructor_(self, x):
        r"""
        Return ``x`` coerced into this forms space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import MeromorphicModularFormsRing
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=12, ep=1)
            sage: (x,y,z,d) = MF.pol_ring().gens()

            sage: Delta = MeromorphicModularFormsRing().Delta()
            sage: Delta.parent()
            MeromorphicModularFormsRing(n=3) over Integer Ring
            sage: MF(Delta)
            q - 24*q^2 + 252*q^3 - 1472*q^4 + O(q^5)
            sage: MF(Delta).parent() == MF
            True

            sage: MF(x^3)
            1 + 720*q + 179280*q^2 + 16954560*q^3 + 396974160*q^4 + O(q^5)
            sage: MF(x^3).parent() == MF
            True

            sage: qexp = Delta.q_expansion(prec=2)
            sage: qexp
            q + O(q^2)
            sage: qexp.parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF(qexp)
            q - 24*q^2 + 252*q^3 - 1472*q^4 + O(q^5)
            sage: MF(qexp) == MF(Delta)
            True

            sage: MF([0,1]) == MF(Delta)
            True
            sage: MF([1,0]) == MF(x^3) - 720*MF(Delta)
            True

            sage: vec = MF(Delta).coordinate_vector()
            sage: vec
            (0, 1)
            sage: vec.parent()
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: vec in MF.module()
            True
            sage: MF(vec) == MF(Delta)
            True

            sage: subspace = MF.subspace([MF(Delta)])
            sage: subspace
            Subspace of dimension 1 of ModularForms(n=3, k=12, ep=1) over Integer Ring
            sage: subspace(MF(Delta)) == subspace(d*(x^3-y^2)) == subspace(qexp) == subspace([0,1]) == subspace(vec) == subspace.gen()
            True
            sage: subspace(MF(Delta)).parent() == subspace(d*(x^3-y^2)).parent() == subspace(qexp).parent() == subspace([0,1]).parent() == subspace(vec).parent()
            True
            sage: subspace([1]) == subspace.gen()
            True
            sage: ssvec = subspace(vec).coordinate_vector()
            sage: ssvec
            (1)
            sage: ssvec.parent()
            Vector space of dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: ambvec = subspace(vec).ambient_coordinate_vector()
            sage: ambvec
            (0, 1)
            sage: ambvec.parent()
            Vector space of degree 2 and dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            Basis matrix:
            [0 1]
            sage: subspace(ambvec) == subspace(vec) and subspace(ambvec).parent() == subspace(vec).parent()
            True
        """

        from graded_ring_element import FormsRingElement
        if isinstance(x, FormsRingElement):
            return self.element_class(self, x._rat)
        if hasattr(x, 'parent') and (is_LaurentSeriesRing(x.parent()) or is_PowerSeriesRing(x.parent())):
            # This assumes that the series corresponds to a weakly holomorphic modular form!
            # But the construction method (with the assumption) may also be used for more general form spaces...
            return self.construct_form(x)
        if is_FreeModuleElement(x) and (self.module() == x.parent() or self.ambient_module() == x.parent()):
            return self.element_from_ambient_coordinates(x)
        if (not self.is_ambient()) and (isinstance(x, list) or isinstance(x, tuple) or is_FreeModuleElement(x)) and len(x) == self.rank():
            try:
                return self.element_from_coordinates(x)
            except ArithmeticError, TypeError:
                pass
        if hasattr(x, 'parent') and self.ambient_module() and self.ambient_module().has_coerce_map_from(x.parent()):
            return self.element_from_ambient_coordinates(self.ambient_module()(x))
        if (isinstance(x,list) or isinstance(x, tuple)) and len(x) == self.degree():
            try:
                return self.element_from_ambient_coordinates(x)
            except ArithmeticError, TypeError:
                pass

        return self.element_class(self, x)

    def _coerce_map_from_(self, S):
        r"""
        Return whether or not there exists a coercion from ``S`` to ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiWeakModularForms, ModularForms, CuspForms, ZeroForm
            sage: MF1 = QuasiWeakModularForms(n=4, base_ring=CC, k=0, ep=1)
            sage: MF2 = ModularForms(n=4, k=24, ep=1)
            sage: MF3 = ModularForms(n=4, k=24, ep=-1)
            sage: MF4 = CuspForms(n=4, k=0, ep=1)
            sage: MF5 = ZeroForm(n=4, k=10, ep=-1)
            sage: subspace1 = MF3.subspace([MF3.gen(0), MF3.gen(1)])
            sage: subspace2 = MF3.subspace([MF3.gen(2)])
            sage: subspace3 = MF3.subspace([MF3.gen(0), MF3.gen(0)+MF3.gen(2)])

            sage: MF2.has_coerce_map_from(MF3)
            False
            sage: MF1.has_coerce_map_from(MF4)
            True
            sage: MF4.has_coerce_map_from(MF5)
            True
            sage: MF4.has_coerce_map_from(ZZ)
            False
            sage: MF1.has_coerce_map_from(ZZ)
            True
            sage: MF3.has_coerce_map_from(subspace1)
            True
            sage: subspace1.has_coerce_map_from(MF3)
            False
            sage: subspace3.has_coerce_map_from(subspace1)
            False
            sage: subspace3.has_coerce_map_from(subspace2)
            True
        """

        from space import ZeroForm
        from subspace import SubSpaceForms
        if   (  isinstance(S, ZeroForm)):
            return True
        elif (  isinstance(S, SubSpaceForms)\
            and isinstance(self, SubSpaceForms) ):
                if (self.ambient_space().has_coerce_map_from(S.ambient_space())\
                    and self.module().has_coerce_map_from(S.module()) ):
                        return True
                else:
                        return False
        elif (  isinstance(S, FormsSpace_abstract)\
            and self.graded_ring().has_coerce_map_from(S.graded_ring())\
            and S.weight()    == self._weight\
            and S.ep()        == self._ep\
            and not isinstance(self, SubSpaceForms)):
                return True
        elif (self.contains_coeff_ring() and self.coeff_ring().has_coerce_map_from(S) ):
            return True
        else:
            return False

    # Since forms spaces are modules instead of rings
    # we have to manually define the one element.
    # This makes it possible to take negative powers of elements.
    @cached_method
    def one_element(self):
        r"""
        Return the one element from the corresponding space
        of constant forms.
        
        Note: The one element does not lie in ``self`` in general.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(k=12)
            sage: MF.one_element()
            1 + O(q^5)
            sage: MF.one_element().parent()
            ModularForms(n=3, k=0, ep=1) over Integer Ring
        """

        return self.extend_type("holo", ring=True)(1).reduce()

    def is_ambient(self):
        r"""
        Return whether ``self`` is an ambient space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=12)
            sage: MF.is_ambient()
            True
            sage: MF.subspace([MF.gen(0)]).is_ambient()
            False
        """

        return self._ambient_space == self

    def ambient_space(self):
        r"""
        Return the ambient space of self.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=12)
            sage: MF.ambient_space()
            ModularForms(n=3, k=12, ep=1) over Integer Ring
            sage: MF.ambient_space() == MF
            True
            sage: subspace = MF.subspace([MF.gen(0)])
            sage: subspace
            Subspace of dimension 1 of ModularForms(n=3, k=12, ep=1) over Integer Ring
            sage: subspace.ambient_space() == MF
            True
        """

        return self._ambient_space

    def module(self):
        r"""
        Return the module associated to self.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=12)
            sage: MF.module()
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: subspace = MF.subspace([MF.gen(0)])
            sage: subspace.module()
            Vector space of degree 2 and dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            Basis matrix:
            [1 0]
        """

        return self._module

    def ambient_module(self):
        r"""
        Return the module associated to the ambient space of self.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=12)
            sage: MF.ambient_module()
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_module() == MF.module()
            True
            sage: subspace = MF.subspace([MF.gen(0)])
            sage: subspace.ambient_module() == MF.module()
            True
        """

        return self._ambient_space._module

    def subspace(self, basis):
        r"""
        Return the subspace of ``self`` generated by ``basis``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=24)
            sage: MF.dimension()
            3
            sage: subspace = MF.subspace([MF.gen(0), MF.gen(1)])
            sage: subspace
            Subspace of dimension 2 of ModularForms(n=3, k=24, ep=1) over Integer Ring
        """

        from subspace import SubSpaceForms
        return SubSpaceForms(self, basis)

    def change_ring(self, new_base_ring):
        r"""
        Return the same space as ``self`` but over a new base ring ``new_base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: CuspForms(n=5, k=24).change_ring(CC)
            CuspForms(n=5, k=24, ep=1) over Complex Field with 53 bits of precision
        """

        return self.__class__.__base__(self.group(), new_base_ring, self.weight(), self.ep())

    def construction(self):
        r"""
        Return a functor that constructs ``self`` (used by the coercion machinery).

        EXAMPLE:: 

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: QuasiModularForms(n=4, k=2, ep=1, base_ring=CC).construction()
            (QuasiModularFormsFunctor(n=4, k=2, ep=1),
             BaseFacade(Complex Field with 53 bits of precision))

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF=ModularForms(k=12)
            sage: MF.subspace([MF.gen(1)]).construction()
            (FormsSubSpaceFunctor with 1 generator for the ModularFormsFunctor(n=3, k=12, ep=1), BaseFacade(Integer Ring))
        """

        from functors import FormsSubSpaceFunctor, FormsSpaceFunctor, BaseFacade
        ambient_space_functor = FormsSpaceFunctor(self._analytic_type, self._group, self._weight, self._ep)

        if (self.is_ambient()):
            return (ambient_space_functor, BaseFacade(self._base_ring))
        else:
            return (FormsSubSpaceFunctor(ambient_space_functor, self._basis), BaseFacade(self._base_ring))

    @cached_method
    def weight(self):
        r"""
        Return the weight of (elements of) ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: QuasiModularForms(n=16, k=16/7, ep=-1).weight()
            16/7
        """

        return self._weight

    @cached_method
    def ep(self):
        r"""
        Return the multiplier of (elements of) ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: QuasiModularForms(n=16, k=16/7, ep=-1).ep()
            -1
        """

        return self._ep

    @cached_method
    def contains_coeff_ring(self):
        r"""
        Return whether ``self`` contains its coefficient ring.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: QuasiModularForms(k=0, ep=1, n=8).contains_coeff_ring()
            True
            sage: QuasiModularForms(k=0, ep=-1, n=8).contains_coeff_ring()
            False
        """

        return ((self.AT("holo") <= self._analytic_type) and (self.weight()==QQ(0)) and (self.ep()==ZZ(1)))

    def element_from_coordinates(self, vec):
        r"""
        If ``self`` has an associated free module, then return the element of ``self``
        corresponding to the given coordinate vector ``vec``. Otherwise raise an exception.

        INPUT:

        - ``vec``     - A coordinate vector with respect to ``self.gens()``.

        OUTPUT:

        An element of ``self`` corresponding to the coordinate vector ``vec``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=24)
            sage: MF.dimension()
            3
            sage: el = MF.element_from_coordinates([1,1,1])
            sage: el
            1 + q + q^2 + 52611612*q^3 + 39019413208*q^4 + O(q^5)
            sage: el == MF.gen(0) + MF.gen(1) + MF.gen(2)
            True
            sage: el.parent() == MF
            True
            
            sage: subspace = MF.subspace([MF.gen(0), MF.gen(1)])
            sage: el = subspace.element_from_coordinates([1,1])
            sage: el
            1 + q + 52611660*q^3 + 39019412128*q^4 + O(q^5)
            sage: el == subspace.gen(0) + subspace.gen(1)
            True
            sage: el.parent() == subspace
            True
        """

        if not self.module():
            raise ValueError("No free module defined for {}".format(self))
        basis = self.gens()
        assert(len(basis) == len(vec))
        # vec = self.module()(self.module().linear_combination_of_basis(vec))
        # this also handles the trivial case (dimension 0)
        return self(sum([vec[k]*basis[k] for k in range(0, len(vec))]))

    def element_from_ambient_coordinates(self, vec):
        r"""
        If ``self`` has an associated free module, then return the element of ``self``
        corresponding to the given ``vec``. Otherwise raise an exception.

        INPUT:

        - ``vec``     - An element of ``self.module()`` or ``self.ambient_module()``.

        OUTPUT:

        An element of ``self`` corresponding to ``vec``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(k=24)
            sage: MF.dimension()
            3
            sage: el = MF.element_from_ambient_coordinates([1,1,1])
            sage: el == MF.element_from_coordinates([1,1,1])
            True
            sage: el.parent() == MF
            True

            sage: subspace = MF.subspace([MF.gen(0), MF.gen(1)])
            sage: el = subspace.element_from_ambient_coordinates([1,1,0])
            sage: el
            1 + q + 52611660*q^3 + 39019412128*q^4 + O(q^5)
            sage: el.parent() == subspace
            True
        """

        return self(self.ambient_space().element_from_coordinates(vec))

    def homogeneous_space(self, k, ep):
        r"""
        Since ``self`` already is a homogeneous component return ``self``
        unless the degree differs in which case an Exception is raised.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: MF = QuasiMeromorphicModularForms(n=6, k=4)
            sage: MF == MF.homogeneous_space(4,1)
            True
            sage: MF.homogeneous_space(5,1)
            Traceback (most recent call last):
            ...
            ValueError: QuasiMeromorphicModularForms(n=6, k=4, ep=1) over Integer Ring already is homogeneous with degree (4, 1) != (5, 1)!
        """

        if (k==self._weight and ep==self._ep):
            return self
        else:
            raise ValueError("{} already is homogeneous with degree ({}, {}) != ({}, {})!".format(self, self._weight, self._ep, k, ep))

    def weight_parameters(self):
        r"""
        Check whether ``self`` has a valid weight and multiplier.
        If not then an exception is raised. Otherwise the two weight
        paramters corresponding to the weight and multiplier of ``self``
        are returned.

        The weight parameters are e.g. used to calculate dimensions
        or precisions of Fourier expansion.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import MeromorphicModularForms
            sage: MF = MeromorphicModularForms(n=18, k=-7, ep=-1)
            sage: MF.weight_parameters()
            (-3, 17)
            sage: (MF._l1, MF._l2) == MF.weight_parameters()
            True
            sage: (k, ep) = (MF.weight(), MF.ep())
            sage: n = MF.hecke_n()
            sage: k == 4*(n*MF._l1 + MF._l2)/(n-2) + (1-ep)*n/(n-2)
            True

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=5, k=12, ep=1)
            sage: MF.weight_parameters()
            (1, 4)
            sage: (MF._l1, MF._l2) == MF.weight_parameters()
            True
            sage: (k, ep) = (MF.weight(), MF.ep())
            sage: n = MF.hecke_n()
            sage: k == 4*(n*MF._l1 + MF._l2)/(n-2) + (1-ep)*n/(n-2)
            True
            sage: MF.dimension() == MF._l1 + 1
            True
        """

        n = self._group.n()
        k = self._weight
        ep = self._ep
        num = (k-(1-ep)*ZZ(n)/ZZ(n-2))*ZZ(n-2)/ZZ(4)

        if (num.is_integral()):
            num = ZZ(num)
            l2 = num%n
            l1 = ((num-l2)/n).numerator()
        else:
            raise ValueError('Invalid or non-occuring weight!')
        return (l1,l2)

    # TODO: this only makes sense for modular forms,
    # resp. needs a big adjustment for quasi modular forms
    def aut_factor(self,gamma,t):
        r"""
        The automorphy factor of ``self``.

        For now it is only defined on the two basic generators of the
        Hecke group of ``self`` and their inverses.

        However, when determening the map which sends an element ``t``
        of the upper half plane to the fundamental domain, the
        function ``self.group().get_FD(t, self.aut_factor)`` can be used.
        It returns the full automorphy factor of the transformation matrix
        applied to ``t``.        

        INPUT:

        - ``gamma``   -- An element of the group of ``self``. For now only
                         the basic generators (and their inverses) are supported.

        - ``t``       -- An element of the upper half plane.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=8, k=4, ep=1)
            sage: full_factor = lambda mat, t: (mat[1][0]*t+mat[1][1])**4
            sage: T = MF.group().T()
            sage: S = MF.group().S()
            sage: z = AlgebraicField()(1+i/2)

            sage: MF.aut_factor(S, z)
            3/2*I - 7/16
            sage: MF.aut_factor(T, z)
            1
            sage: MF.aut_factor(S, z) == full_factor(S, z)
            True
            sage: MF.aut_factor(T, z) == full_factor(T, z)
            True

            sage: (A, w, fact) = MF.group().get_FD(z, MF.aut_factor)
            sage: fact == full_factor(A,w)
            True

            sage: MF = ModularForms(n=7, k=14/5, ep=-1)
            sage: T = MF.group().T()
            sage: S = MF.group().S()
            sage: z = AlgebraicField()(1+i/2)

            sage: MF.aut_factor(S, z)
            1.3655215324256...? + 0.056805991182877...?*I
            sage: MF.aut_factor(T, z)
            1
            sage: MF.aut_factor(S, z) == MF.ep() * (z/i)^MF.weight()
            True
        """

        if (gamma == self._group.T() or gamma == self._group.T().inverse()):
            return 1
        elif (gamma == self._group.S() or gamma == -self._group.S()):
            return self._ep*(t/AlgebraicField()(i))**self._weight
        else:
            raise NotImplementedError("Factor of autormorphy is implemented only for some group elements.")

    @cached_method
    def F_simple(self):
        r"""
        Return a (the most) simple normalized element of ``self``
        corresponding to the weight parameters ``l1=self._l1`` and
        ``l2=self._l2``. If the element does not lie in ``self`` the
        type of its parent is extended accordingly.

        The main part of the element is given by the ``l1``-th power of ``f_inf``,
        up to a small holomorphic correction factor.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(n=18, k=-7, ep=-1)
            sage: MF.disp_prec(1)
            sage: MF.F_simple()
            q^-3 + 16/(81*d)*q^-2 - 4775/(104976*d^2)*q^-1 - 14300/(531441*d^3) + O(q)
            sage: MF.F_simple() == MF.f_inf()^MF._l1 * MF.f_rho()^MF._l2 * MF.f_i()
            True

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(n=5, k=2, ep=-1)
            sage: MF._l1
            -1
            sage: MF.F_simple().parent()
            WeakModularForms(n=5, k=2, ep=-1) over Integer Ring
        """

        (x,y,z,d) = self.rat_field().gens()

        finf_pol = d*(x**self._group.n() - y**2)
        rat = finf_pol**self._l1 * x**self._l2 * y**(ZZ(1-self._ep)/ZZ(2))

        if (self._l1 > 0):
            new_space = self.extend_type("cusp")
        elif (self._l1 == 0):
            new_space = self.extend_type("holo")
        else:
            new_space = self.extend_type("weak")

        return new_space(rat)

    def Faber_pol(self, m, fix_d=False, d=None, d_num_prec=None):
        r"""
        Return the ``m``'th Faber polynomial of ``self``.

        Namely a polynomial `P(q)` such that ``P(J_inv)*F_simple()``
        has a Fourier expansion of the form ``q^(-m) + O(q^(self._l1 + 1))``.
        ``d^(self._l1 + m)*P(q)`` is a monic polynomial of degree ``self._l1 + m``.

        The Faber polynomials are used to construct a basis of weakly holomorphic forms
        and to recover such forms from their initial Fourier coefficients.

        INPUT:

        - ``m``           -- An integer ``m >= -self._l1``.

        - ``fix_d``       -- ``True`` if the value of ``d`` should be
                             (numerically) substituted for the coefficients
                             of the polynomial (default: ``False``).

        - ``d``           -- The value which should be substituted for ``d`` (default: ``None``).
                             The value is ignored if ``fix_d=True``.                             

        - ``d_num_prec``  -- The numerical precision to be used for ``d``
                             in case ``fix_d=True`` or ``d`` is set,
                             Default: ``None``, in which case the default
                             numerical precision ``self.num_prec()`` is used.

        OUTPUT:

        The corresponding Faber polynomial `P(q)` with coefficients in ``self.coeff_ring()``
        resp. a numerical ring in case ``fix_d=True`` or ``d`` is set
        (and the group of ``self`` is not arithmetic).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(n=5, k=62/3, ep=-1)
            sage: MF.weight_parameters()
            (2, 3)

            sage: MF.Faber_pol(-2)
            1
            sage: MF.Faber_pol(-1)
            1/d*q - 19/(100*d)
            sage: MF.Faber_pol(0)
            1/d^2*q^2 - 117/(200*d^2)*q + 9113/(320000*d^2)
            sage: MF.Faber_pol(2)
            1/d^4*q^4 - 11/(8*d^4)*q^3 + 41013/(80000*d^4)*q^2 - 2251291/(48000000*d^4)*q + 1974089431/(4915200000000*d^4)            
            sage: (MF.Faber_pol(-2)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2)
            q^2 - 41/(200*d)*q^3 + O(q^4)
            sage: (MF.Faber_pol(-1)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+1)
            q + O(q^3)
            sage: (MF.Faber_pol(0)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+1)
            1 + O(q^3)
            sage: (MF.Faber_pol(2)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+1)
            q^-2 + O(q^3)

            sage: MF.Faber_pol(-2, d=1)
            1
            sage: MF.Faber_pol(-1, d=1)
            q - 19/100
            sage: MF.Faber_pol(2, d=1)
            q^4 - 11/8*q^3 + 41013/80000*q^2 - 2251291/48000000*q + 1974089431/4915200000000
            sage: (MF.Faber_pol(-2, d=1)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2, d=1)
            q^2 - 41/200*q^3 + O(q^4)
            sage: (MF.Faber_pol(2)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+1, d=1)
            q^-2 + O(q^3)

            sage: MF = WeakModularForms(n=4, k=-2)
            sage: MF.weight_parameters()
            (-1, 3)

            sage: MF.Faber_pol(1)
            1
            sage: MF.Faber_pol(2, fix_d=True)
            256*q - 184
            sage: MF.Faber_pol(3, fix_d=True)
            65536*q^2 - 73728*q + 14364
            sage: (MF.Faber_pol(1, fix_d=True)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2, fix_d=True)
            q^-1 + 80 + O(q)
            sage: (MF.Faber_pol(2, fix_d=True)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2, fix_d=True)
            q^-2 + 400 + O(q)
            sage: (MF.Faber_pol(3)(MF.J_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2, fix_d=True)
            q^-3 + 2240 + O(q)
        """

        m = ZZ(m)
        if (m < -self._l1):
            raise ValueError("Invalid basis index: {}<{}!".format(m,-self._l1))
        if (d_num_prec == None):
            d_num_prec = self._num_prec

        prec          = 2*self._l1+m+1
        (base_ring, coeff_ring, qseries_ring, d) = MFSeriesConstructor(self._group, prec).series_data(self.base_ring(), fix_d, d, d_num_prec)
        q             = qseries_ring.gen()

        simple_qexp   = self.F_simple().q_expansion(prec=prec, fix_d=fix_d, d=d, d_num_prec=d_num_prec)
        J_qexp        = self.J_inv().q_expansion(prec=m + self._l1, fix_d=fix_d, d=d, d_num_prec=d_num_prec)

        # The precision could be infinity, otherwise we could do this:
        #assert(temp_reminder.prec() == 1)
        temp_reminder = (1 / simple_qexp / q**m).add_bigoh(1)

        fab_pol       = qseries_ring([])
        while (len(temp_reminder.coefficients()) > 0):
            temp_coeff     = temp_reminder.coefficients()[0]
            temp_exp       = -temp_reminder.exponents()[0]
            fab_pol       += temp_coeff*(q/d)**temp_exp

            temp_reminder -= temp_coeff*(J_qexp/d)**temp_exp
            # The first term is zero only up to numerical errors,
            # so we manually have to remove it
            if (not base_ring.is_exact()):
                temp_reminder=temp_reminder.truncate_neg(-temp_exp+1)

        return fab_pol.polynomial()

    # very similar to Faber_pol: faber_pol(q)=Faber_pol(d*q)
    def faber_pol(self, m, fix_d=False, d=None, d_num_prec=None):
        r"""
        Return the `m`'th Faber polynomial of ``self``
        with a different normalization based on ``j_inv``
        instead of ``J_inv``.

        Namely a polynomial p(q) such that ``p(j_inv)*F_simple()``
        has a Fourier expansion of the form ``q^(-m) + O(q^(self._l1 + 1))``.
        ``p(q)`` is a monic polynomial of degree ``self._l1 + m``.

        The relation to ``Faber_pol`` is: ``faber_pol(q) = Faber_pol(d*q)``.

        INPUT:

        - ``m``           -- An integer ``m >= -self._l1``.

        - ``fix_d``       -- ``True`` if the value of ``d`` should be
                             (numerically) substituted for the coefficients
                             of the polynomial (default: ``False``).

        - ``d``           -- The value which should be substituted for ``d`` (default: ``None``).
                             The value is ignored if ``fix_d=True``.                             

        - ``d_num_prec``  -- The numerical precision to be used for ``d``
                             in case ``fix_d=True`` or ``d`` is set,
                             Default: ``None``, in which case the default
                             numerical precision ``self.num_prec()`` is used.

        OUTPUT:

        The corresponding Faber polynomial p(q) with coefficients in ``self.coeff_ring()``
        resp. a numerical ring in case ``fix_d=True`` or ``d`` is set
        (and the group of ``self`` is not arithmetic).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(n=5, k=62/3, ep=-1)
            sage: MF.weight_parameters()
            (2, 3)

            sage: MF.faber_pol(-2)
            1
            sage: MF.faber_pol(-1)
            q - 19/(100*d)
            sage: MF.faber_pol(0)
            q^2 - 117/(200*d)*q + 9113/(320000*d^2)
            sage: MF.faber_pol(2)
            q^4 - 11/(8*d)*q^3 + 41013/(80000*d^2)*q^2 - 2251291/(48000000*d^3)*q + 1974089431/(4915200000000*d^4)
            sage: (MF.faber_pol(-2)(MF.j_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2)
            q^2 - 41/(200*d)*q^3 + O(q^4)
            sage: (MF.faber_pol(-1)(MF.j_inv())*MF.F_simple()).q_expansion(prec=MF._l1+1)
            q + O(q^3)
            sage: (MF.faber_pol(0)(MF.j_inv())*MF.F_simple()).q_expansion(prec=MF._l1+1)
            1 + O(q^3)
            sage: (MF.faber_pol(2)(MF.j_inv())*MF.F_simple()).q_expansion(prec=MF._l1+1)
            q^-2 + O(q^3)

            sage: MF = WeakModularForms(n=4, k=-2)
            sage: MF.weight_parameters()
            (-1, 3)

            sage: MF.faber_pol(1)
            1
            sage: MF.faber_pol(2, fix_d=True)
            q - 184
            sage: MF.faber_pol(3, fix_d=True)
            q^2 - 288*q + 14364
            sage: (MF.faber_pol(1, fix_d=True)(MF.j_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2, fix_d=True)
            q^-1 + 80 + O(q)
            sage: (MF.faber_pol(2, fix_d=True)(MF.j_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2, fix_d=True)
            q^-2 + 400 + O(q)
            sage: (MF.faber_pol(3)(MF.j_inv())*MF.F_simple()).q_expansion(prec=MF._l1+2, fix_d=True)
            q^-3 + 2240 + O(q)
        """

        m = ZZ(m)
        if (m < -self._l1):
            raise ValueError("Invalid basis index: {}<{}!".format(m,-self._l1))
        if (d_num_prec == None):
            d_num_prec = self._num_prec

        prec          = 2*self._l1+m+1
        (base_ring, coeff_ring, qseries_ring, d) = MFSeriesConstructor(self._group, prec).series_data(self.base_ring(), fix_d, d, d_num_prec)
        q             = qseries_ring.gen()

        simple_qexp   = self.F_simple().q_expansion(prec=prec, fix_d=fix_d, d=d, d_num_prec=d_num_prec)
        j_qexp        = self.j_inv().q_expansion(prec=m + self._l1, fix_d=fix_d, d=d, d_num_prec=d_num_prec)

        # The precision could be infinity, otherwise we could do this:
        #assert(temp_reminder.prec() == 1)
        temp_reminder = (1 / simple_qexp / q**m).add_bigoh(1)

        fab_pol       = qseries_ring([])
        while (len(temp_reminder.coefficients()) > 0):
            temp_coeff     = temp_reminder.coefficients()[0]
            temp_exp       = -temp_reminder.exponents()[0]
            fab_pol       += temp_coeff*q**temp_exp

            temp_reminder -= temp_coeff*j_qexp**temp_exp
            # The first term is zero only up to numerical errors,
            # so we manually have to remove it
            if (not base_ring.is_exact()):
                temp_reminder=temp_reminder.truncate_neg(-temp_exp+1)

        return fab_pol.polynomial()

    def F_basis_pol(self, m):
        r"""
        Returns a polynomial corresponding to the basis
        element of the correponding space of weakly holomorphic
        forms of the same degree as ``self``. The basis element
        is determined by the property that the Fourier expansion
        is of the form ``q^(-m) + O(q^(self._l1 + 1))``.

        INPUT:
        
        - ``m``           -- An integer ``m >= -self._l1``.

        OUTPUT:

        A polynomial in ``x,y,z,d``, corresponding to ``f_rho, f_i, E2``
        and the (possibly) transcendental parameter ``d``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(n=5, k=62/3, ep=-1)
            sage: MF.weight_parameters()
            (2, 3)

            sage: MF.F_basis_pol(-2)
            x^13*y*d^2 - 2*x^8*y^3*d^2 + x^3*y^5*d^2
            sage: MF.F_basis_pol(-1)
            (-81*x^13*y*d + 62*x^8*y^3*d + 19*x^3*y^5*d)/(-100)
            sage: MF.F_basis_pol(0)
            (141913*x^13*y + 168974*x^8*y^3 + 9113*x^3*y^5)/320000
            
            sage: MF(MF.F_basis_pol(-2)).q_expansion(prec=MF._l1+2)
            q^2 - 41/(200*d)*q^3 + O(q^4)
            sage: MF(MF.F_basis_pol(-1)).q_expansion(prec=MF._l1+1)
            q + O(q^3)
            sage: MF(MF.F_basis_pol(0)).q_expansion(prec=MF._l1+1)
            1 + O(q^3)
            sage: MF(MF.F_basis_pol(2)).q_expansion(prec=MF._l1+1)
            q^-2 + O(q^3)
            sage: MF(MF.F_basis_pol(2)).parent()
            WeakModularForms(n=5, k=62/3, ep=-1) over Integer Ring

            sage: MF = WeakModularForms(n=4, k=-2)
            sage: MF.weight_parameters()
            (-1, 3)

            sage: MF.F_basis_pol(1)
            x^3/(x^4*d - y^2*d)
            sage: MF.F_basis_pol(2)
            (9*x^7 + 23*x^3*y^2)/(32*x^8*d^2 - 64*x^4*y^2*d^2 + 32*y^4*d^2)

            sage: MF(MF.F_basis_pol(1)).q_expansion(prec=MF._l1+2)
            q^-1 + 5/(16*d) + O(q)
            sage: MF(MF.F_basis_pol(2)).q_expansion(prec=MF._l1+2)
            q^-2 + 25/(4096*d^2) + O(q)
        """

        (x,y,z,d) = self.rat_field().gens()

        n        = self._group.n()
        finf_pol = d*(x**n-y**2)
        jinv_pol = x**n/(x**n-y**2)
        rat      = finf_pol**self._l1 * x**self._l2 * y**(ZZ(1-self._ep)/ZZ(2)) * self.Faber_pol(m)(jinv_pol)

        #return self(self.F_simple()*self.Faber_pol(m)(self.graded_ring().J_inv()))
        return rat

    def F_basis(self, m):
        r"""
        Returns a weakly holomorphic element of ``self``
        (extended if necessarily) determined by the property that
        the Fourier expansion is of the form is of the form
        ``q^(-m) + O(q^(self._l1 + 1))``.

        In particular for all ``m >= -self._l1`` these elements form
        a basis of the space of weakly holomorphic modular forms
        of the corresponding degree.

        INPUT:
        
        - ``m``            - An integer ``m >= -self._l1``.

        OUTPUT:

        The corresponding element in (possibly an extension of) ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms, CuspForms
            sage: MF = WeakModularForms(n=5, k=62/3, ep=-1)
            sage: MF.disp_prec(MF._l1+2)
            sage: MF.weight_parameters()
            (2, 3)

            sage: MF.F_basis(-2)
            q^2 - 41/(200*d)*q^3 + O(q^4)
            sage: MF.F_basis(-1)
            q - 13071/(640000*d^2)*q^3 + O(q^4)
            sage: MF.F_basis(0)
            1 - 277043/(192000000*d^3)*q^3 + O(q^4)
            sage: MF.F_basis(2)
            q^-2 - 162727620113/(40960000000000000*d^5)*q^3 + O(q^4)
            sage: MF.F_basis(2).parent() == MF
            True

            sage: MF = CuspForms(n=4, k=-2)
            sage: MF.weight_parameters()
            (-1, 3)

            sage: MF.F_basis(1).parent()
            WeakModularForms(n=4, k=-2, ep=1) over Integer Ring
            sage: MF.F_basis(1).parent().disp_prec(MF._l1+2)
            sage: MF.F_basis(1)
            q^-1 + 80 + O(q)
            sage: MF.F_basis(2)
            q^-2 + 400 + O(q)
        """

        if (m < 0):
            new_space = self.extend_type("cusp")
        elif (m == 0):
            new_space = self.extend_type("holo")
        else:
            new_space = self.extend_type("weak")

        return new_space(self.F_basis_pol(m))

    # TODO: This only works for weakly holomorphic modular forms!
    def construct_form(self, laurent_series):
        r"""
        Tries to construct an element of self with the given Fourier
        expansion. The assumption is made that the specified Fourier
        expansion corresponds to a weakly holomorphic modular form.

        If the precision is too low to determine the
        element an exception is raised.

        INPUT:

        - ``laurent_series``  -- A Laurent or Power series.        

        OUTPUT:

        If possible: An element of self with the same initial
        Fourier expansion as ``laurent_series``.

        Note: Instead of calling ``construct_form`` directly it is also
        possible to use ``self(laurent_series)`` which then calls this
        function.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: Delta = CuspForms(k=12).Delta()
            sage: qexp = Delta.q_expansion(prec=2)
            sage: qexp.parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: qexp
            q + O(q^2)
            sage: CuspForms(k=12).construct_form(qexp) == Delta
            True

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: J_inv = WeakModularForms(n=7).J_inv()
            sage: qexp = J_inv.q_expansion(prec=1)
            sage: qexp.parent()
            Laurent Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: qexp
            d*q^-1 + 151/392 + O(q)
            sage: WeakModularForms(n=7).construct_form(qexp) == J_inv
            True

            sage: MF = WeakModularForms(n=5, k=62/3, ep=-1)
            sage: MF.default_prec(MF._l1+1)
            sage: d = MF.coeff_ring().gen()
            sage: MF.weight_parameters()
            (2, 3)
            sage: fun = d*MF.F_basis(-2) + 2*MF.F_basis(-1) + MF.F_basis(2)
            sage: qexp = fun.q_expansion()
            sage: qexp.parent()
            Laurent Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: qexp
            q^-2 + 2*q + d*q^2 + O(q^3)
            sage: WeakModularForms(n=5, k=62/3, ep=-1).construct_form(qexp) == fun
            True
        """

        if (laurent_series.prec() < self._l1+1):
            raise ValueError('Insufficient precision!')

        laurent_series = laurent_series.add_bigoh(self._l1+1)
        coefficients   = laurent_series.coefficients()
        exponents      = laurent_series.exponents()

        if (len(coefficients) == 0):
            return self(0)

        rat = sum([\
                  coefficients[j] * self.F_basis_pol(-exponents[j])\
                  for j in range(ZZ(0), ZZ(len(coefficients)))
              ])

        return self(rat)


    # DEFAULT METHODS (should be overwritten in concrete classes)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: el = QuasiMeromorphicModularForms(k=2, ep=-1).an_element()
            sage: el.parent()
            QuasiMeromorphicModularForms(n=3, k=2, ep=-1) over Integer Ring
            sage: el.is_zero()
            True
            sage: el
            O(q^5)
        """

        # this seems ok, so might as well leave it as is for everything
        return self(ZZ(0))
        #return self.F_simple()

    @cached_method
    def dimension(self):
        r"""
        This method should be overloaded by subclasses.

        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: QuasiMeromorphicModularForms(k=2, ep=-1).dimension()
            +Infinity
        """

        return infinity
 
    def rank(self):
        r"""
        Return the rank of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=4, k=24, ep=-1)
            sage: MF.rank()
            3
            sage: MF.subspace([MF.gen(0), MF.gen(2)]).rank()
            2
        """

        return self.dimension()

    def degree(self):
        r"""
        Return the degree of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=4, k=24, ep=-1)
            sage: MF.degree()
            3
            sage: MF.subspace([MF.gen(0), MF.gen(2)]).degree() # defined in subspace.py
            3
        """

        return self.dimension()
 
    def coordinate_vector(self, v):
        r"""
        This method should be overloaded by subclasses.

        Return the coordinate vector of the element ``v``
        with respect to ``self.gens()``.

        NOTE:

        Elements use this method (from their parent)
        to calculate their coordinates.

        INPUT:

        - `v` -- An element of ``self``.
        
        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=4, k=24, ep=-1)
            sage: MF.coordinate_vector(MF.gen(0)).parent() # defined in space.py
            Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.coordinate_vector(MF.gen(0))          # defined in space.py
            (1, 0, 0)
            sage: subspace = MF.subspace([MF.gen(0), MF.gen(2)])
            sage: subspace.coordinate_vector(subspace.gen(0)).parent()  # defined in subspace.py
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: subspace.coordinate_vector(subspace.gen(0))           # defined in subspace.py
            (1, 0)
        """

        raise NotImplementedError()

    @cached_method
    def ambient_coordinate_vector(self, v):
        r"""
        Return the coordinate vector of the element ``v``
        in ``self.module()`` with respect to the basis
        from ``self.ambient_space``.
        
        NOTE:

        Elements use this method (from their parent)
        to calculate their coordinates.

        INPUT:

        - `v` -- An element of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=4, k=24, ep=-1)
            sage: MF.ambient_coordinate_vector(MF.gen(0)).parent()
            Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_coordinate_vector(MF.gen(0))
            (1, 0, 0)
            sage: subspace = MF.subspace([MF.gen(0), MF.gen(2)])
            sage: subspace.ambient_coordinate_vector(subspace.gen(0)).parent()
            Vector space of degree 3 and dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            Basis matrix:
            [1 0 0]
            [0 0 1]
            sage: subspace.ambient_coordinate_vector(subspace.gen(0))
            (1, 0, 0)
        """

        return self.module()(self.ambient_space().coordinate_vector(v))
 
    def gens(self):
        r"""
        This method should be overloaded by subclasses.

        Return a basis of ``self``.
        
        Note that the coordinate vector of elements of ``self``
        are with respect to this basis.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ModularForms(k=12).gens() # defined in space.py
            [1 + 196560*q^2 + 16773120*q^3 + 398034000*q^4 + O(q^5),
             q - 24*q^2 + 252*q^3 - 1472*q^4 + O(q^5)]
        """

        raise NotImplementedError()
 
    def gen(self, k=0):
        r"""
        Return the ``k``'th basis element of ``self``
        if possible (default: ``k=0``).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ModularForms(k=12).gen(1).parent()
            ModularForms(n=3, k=12, ep=1) over Integer Ring
            sage: ModularForms(k=12).gen(1)
            q - 24*q^2 + 252*q^3 - 1472*q^4 + O(q^5)
        """

        k = ZZ(k)
        if k>=0 and k < self.dimension():
            return self.gens()[k]
        else:
            raise ValueError("Invalid index: k={} does not satisfy 0 <= k <= {}!".format(k, self.dimension()))
