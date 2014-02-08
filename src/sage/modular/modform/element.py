# -*- coding: utf-8 -*-
"""
Elements of modular forms spaces
"""

#########################################################################
#       Copyright (C) 2004--2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import space
import sage.modular.hecke.element as element
import sage.rings.all as rings
from sage.modular.modsym.space import is_ModularSymbolsSpace
from sage.modules.module_element import ModuleElement
from sage.modules.free_module_element import vector
from sage.misc.misc import verbose
from sage.modular.dirichlet import DirichletGroup

def is_ModularFormElement(x):
    """
    Return True if x is a modular form.

    EXAMPLES::

        sage: from sage.modular.modform.element import is_ModularFormElement
        sage: is_ModularFormElement(5)
        False
        sage: is_ModularFormElement(ModularForms(11).0)
        True
    """
    return isinstance(x, ModularFormElement)

def delta_lseries(prec=53,
                 max_imaginary_part=0,
                 max_asymp_coeffs=40):
    r"""
    Return the L-series of the modular form Delta.

    This actually returns an interface to Tim Dokchitser's program
    for computing with the L-series of the modular form `\Delta`.

    INPUT:

    - ``prec`` - integer (bits precision)

    - ``max_imaginary_part`` - real number

    - ``max_asymp_coeffs`` - integer

    OUTPUT:

    The L-series of `\Delta`.

    EXAMPLES::

        sage: L = delta_lseries()
        sage: L(1)
        0.0374412812685155
    """
    from sage.lfunctions.all import Dokchitser
    key = (prec, max_imaginary_part, max_asymp_coeffs)
    L = Dokchitser(conductor = 1,
                   gammaV = [0,1],
                   weight = 12,
                   eps = 1,
                   prec = prec)
    s = 'tau(n) = (5*sigma(n,3)+7*sigma(n,5))*n/12-35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5));'
    L.init_coeffs('tau(k)',pari_precode = s,
                  max_imaginary_part=max_imaginary_part,
                  max_asymp_coeffs=max_asymp_coeffs)
    L.set_coeff_growth('2*n^(11/2)')
    L.rename('L-series associated to the modular form Delta')
    return L

class ModularForm_abstract(ModuleElement):
    """
    Constructor for generic class of a modular form. This
    should never be called directly; instead one should
    instantiate one of the derived classes of this
    class.
    """
    def group(self):
        """
        Return the group for which self is a modular form.

        EXAMPLES::

            sage: ModularForms(Gamma1(11), 2).gen(0).group()
            Congruence Subgroup Gamma1(11)
        """
        return self.parent().group()

    def weight(self):
        """
        Return the weight of self.

        EXAMPLES::

            sage: (ModularForms(Gamma1(9),2).6).weight()
            2
        """
        return self.parent().weight()

    def level(self):
        """
        Return the level of self.

        EXAMPLES::

            sage: ModularForms(25,4).0.level()
            25
        """
        return self.parent().level()

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: ModularForms(25,4).0._repr_()
            'q + O(q^6)'

            sage: ModularForms(25,4).3._repr_()
            'q^4 + O(q^6)'
        """
        return str(self.q_expansion())

    def _ensure_is_compatible(self, other):
        """
        Make sure self and other are compatible for arithmetic or
        comparison operations. Raise an error if incompatible,
        do nothing otherwise.

        EXAMPLES::

            sage: f = ModularForms(DirichletGroup(17).0^2,2).2
            sage: g = ModularForms(DirichletGroup(17).0^2,2).1
            sage: h = ModularForms(17,4).0

            sage: f._ensure_is_compatible(g)

            sage: f._ensure_is_compatible(h)
            Traceback (most recent call last):
            ...
            ArithmeticError: Modular forms must be in the same ambient space.
        """
        if not isinstance(other, ModularForm_abstract):
            raise TypeError, "Second argument must be a modular form."
        if self.parent().ambient() != other.parent().ambient():
            raise ArithmeticError, "Modular forms must be in the same ambient space."

    def __call__(self, x, prec=None):
        """
        Evaluate the q-expansion of this modular form at x.

        EXAMPLES::

            sage: f = ModularForms(DirichletGroup(17).0^2,2).2

            sage: q = f.q_expansion().parent().gen()
            sage: f(q^2 + O(q^7))
            q^2 + (-zeta8^2 + 2)*q^4 + (zeta8 + 3)*q^6 + O(q^7)

            sage: f(0)
            0
        """
        return self.q_expansion(prec)(x)

    def valuation(self):
        """
        Return the valuation of self (i.e. as an element of the power
        series ring in q).

        EXAMPLES::

            sage: ModularForms(11,2).0.valuation()
            1
            sage: ModularForms(11,2).1.valuation()
            0
            sage: ModularForms(25,6).1.valuation()
            2
            sage: ModularForms(25,6).6.valuation()
            7
        """
        try:
            return self.__valuation
        except AttributeError:
            v = self.qexp().valuation()
            if v != self.qexp().prec():
                self.__valuation = v
                return v
            v = self.qexp(self.parent().sturm_bound()).valuation()
            self.__valuation = v
            return v

    def qexp(self, prec=None):
        """
        Same as ``self.q_expansion(prec)``.

        .. seealso:: :meth:`q_expansion`

        EXAMPLES::

            sage: CuspForms(1,12).0.qexp()
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
        """
        return self.q_expansion(prec)


    def __eq__(self, other):
        """
        Compare self to other.

        EXAMPLES::

            sage: f = ModularForms(6,4).0
            sage: g = ModularForms(23,2).0
            sage: f == g ## indirect doctest
            False
            sage: f == f
            True
            sage: f == loads(dumps(f))
            True
        """
        if not isinstance(other, ModularFormElement) or \
           self.ambient_module() != other.ambient_module():
            return False
        else:
            return self.element() == other.element()

    def __cmp__(self, other):
        """
        Compare self to other. If they are not the same object, but
        are of the same type, compare them as vectors.

        EXAMPLES::

            sage: f = ModularForms(DirichletGroup(17).0^2,2).2
            sage: g = ModularForms(DirichletGroup(17).0^2,2).1
            sage: f == g ## indirect doctest
            False
            sage: f == f
            True
        """
        try:
            self._ensure_is_compatible(other)
        except Exception:
            return self.parent().__cmp__(other.parent())
        if self.element() == other.element():
            return 0
        else:
            return -1

    def _compute(self, X):
        """
        Compute the coefficients of `q^n` of the power series of self,
        for `n` in the list `X`.  The results are not cached.  (Use
        coefficients for cached results).

        EXAMPLES::

            sage: f = ModularForms(18,2).1
            sage: f.q_expansion(20)
            q + 8*q^7 + 4*q^10 + 14*q^13 - 4*q^16 + 20*q^19 + O(q^20)
            sage: f._compute([10,17])
            [4, 0]
            sage: f._compute([])
            []
        """
        if not isinstance(X, list) or len(X) == 0:
            return []
        bound = max(X)
        q_exp = self.q_expansion(bound+1)
        return [q_exp[i] for i in X]

    def coefficients(self, X):
        """
        The coefficients a_n of self, for integers n>=0 in the list
        X. If X is an Integer, return coefficients for indices from 1
        to X.

        This function caches the results of the compute function.

        TESTS::

            sage: e = DirichletGroup(11).gen()
            sage: f = EisensteinForms(e, 3).eisenstein_series()[0]
            sage: f.coefficients([0,1])
            [15/11*zeta10^3 - 9/11*zeta10^2 - 26/11*zeta10 - 10/11,
            1]
            sage: f.coefficients([0,1,2,3])
            [15/11*zeta10^3 - 9/11*zeta10^2 - 26/11*zeta10 - 10/11,
            1,
            4*zeta10 + 1,
            -9*zeta10^3 + 1]
            sage: f.coefficients([2,3])
            [4*zeta10 + 1,
            -9*zeta10^3 + 1]

        Running this twice once revealed a bug, so we test it::

            sage: f.coefficients([0,1,2,3])
            [15/11*zeta10^3 - 9/11*zeta10^2 - 26/11*zeta10 - 10/11,
            1,
            4*zeta10 + 1,
            -9*zeta10^3 + 1]
        """
        try:
            self.__coefficients
        except AttributeError:
            self.__coefficients = {}
        if isinstance(X, rings.Integer):
            X = range(1,X+1)
        Y = [n for n in X   if  not (n in self.__coefficients.keys())]
        v = self._compute(Y)
        for i in range(len(v)):
            self.__coefficients[Y[i]] = v[i]
        return [ self.__coefficients[x] for x in X ]

    def __getitem__(self, n):
        """
        Returns the `q^n` coefficient of the `q`-expansion of self or
        returns a list containing the `q^i` coefficients of self
        where `i` is in slice `n`.

        EXAMPLES::

            sage: f = ModularForms(DirichletGroup(17).0^2,2).2
            sage: f.__getitem__(10)
            zeta8^3 - 5*zeta8^2 - 2*zeta8 + 10
            sage: f[30]
            -2*zeta8^3 - 17*zeta8^2 + 4*zeta8 + 29
            sage: f[10:15]
            [zeta8^3 - 5*zeta8^2 - 2*zeta8 + 10,
            -zeta8^3 + 11,
            -2*zeta8^3 - 6*zeta8^2 + 3*zeta8 + 9,
            12,
            2*zeta8^3 - 7*zeta8^2 + zeta8 + 14]
        """
        if isinstance(n, slice):
            if n.stop is None:
                return self.q_expansion().list()[n]
            else:
                return self.q_expansion(n.stop+1).list()[n]
        else:
            return self.q_expansion(n+1)[int(n)]

    def padded_list(self, n):
        """
        Return a list of length n whose entries are the first n
        coefficients of the q-expansion of self.

        EXAMPLES::

            sage: CuspForms(1,12).0.padded_list(20)
            [0, 1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920, 534612, -370944, -577738, 401856, 1217160, 987136, -6905934, 2727432, 10661420]
        """
        return self.q_expansion(n).padded_list(n)


    def _latex_(self):
        """
        Return the LaTeX expression of self.

        EXAMPLES:
            sage: ModularForms(25,4).0._latex_()
            'q + O(q^{6})'

            sage: ModularForms(25,4).4._latex_()
            'q^{5} + O(q^{6})'
        """
        return self.q_expansion()._latex_()

    def base_ring(self):
        """
        Return the base_ring of self.

        EXAMPLES::

            sage: (ModularForms(117, 2).13).base_ring()
            Rational Field
            sage: (ModularForms(119, 2, base_ring=GF(7)).12).base_ring()
            Finite Field of size 7
        """
        return self.parent().base_ring()

    def character(self, compute=True):
        """
        Return the character of self. If ``compute=False``, then this will
        return None unless the form was explicitly created as an element of a
        space of forms with character, skipping the (potentially expensive)
        computation of the matrices of the diamond operators.

        EXAMPLES::

            sage: ModularForms(DirichletGroup(17).0^2,2).2.character()
            Dirichlet character modulo 17 of conductor 17 mapping 3 |--> zeta8

            sage: CuspForms(Gamma1(7), 3).gen(0).character()
            Dirichlet character modulo 7 of conductor 7 mapping 3 |--> -1
            sage: CuspForms(Gamma1(7), 3).gen(0).character(compute = False) is None
            True
            sage: M = CuspForms(Gamma1(7), 5).gen(0).character()
            Traceback (most recent call last):
            ...
            ValueError: Form is not an eigenvector for <3>
        """
        chi = self.parent().character()
        if (chi is not None) or (not compute):
            return chi
        else: # do the expensive computation
            G = DirichletGroup(self.parent().level(), base_ring = self.parent().base_ring())
            gens = G.unit_gens()
            i = self.valuation()
            vals = []
            for g in gens:
                df = self.parent().diamond_bracket_operator(g)(self)
                if df != (df[i] / self[i]) * self:
                    raise ValueError, "Form is not an eigenvector for <%s>" % g
                vals.append(df[i] / self[i])
            return G(vals)

    def __nonzero__(self):
        """
        Return True if self is nonzero, and False if not.

        EXAMPLES::

            sage: ModularForms(25,6).6.__nonzero__()
            True
        """
        return not self.element().is_zero()

    def prec(self):
        """
        Return the precision to which self.q_expansion() is
        currently known. Note that this may be 0.

        EXAMPLES::

            sage: M = ModularForms(2,14)
            sage: f = M.0
            sage: f.prec()
            0

            sage: M.prec(20)
            20
            sage: f.prec()
            0
            sage: x = f.q_expansion() ; f.prec()
            20
        """
        try:
            return self.__q_expansion[0]
        except AttributeError:
            return 0

    def q_expansion(self, prec=None):
        r"""
        The `q`-expansion of the modular form to precision `O(q^\text{prec})`.
        This function takes one argument, which is the integer prec.

        EXAMPLES:

        We compute the cusp form `\Delta`::

            sage: delta = CuspForms(1,12).0
            sage: delta.q_expansion()
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)

        We compute the `q`-expansion of one of the cusp forms of level 23::

            sage: f = CuspForms(23,2).0
            sage: f.q_expansion()
            q - q^3 - q^4 + O(q^6)
            sage: f.q_expansion(10)
            q - q^3 - q^4 - 2*q^6 + 2*q^7 - q^8 + 2*q^9 + O(q^10)
            sage: f.q_expansion(2)
            q + O(q^2)
            sage: f.q_expansion(1)
            O(q^1)
            sage: f.q_expansion(0)
            O(q^0)
        """
        if prec is None:
            prec = self.parent().prec()
        prec = rings.Integer(prec)
        if prec < 0:
            raise ValueError, "prec (=%s) must be at least 0"%prec
        try:
            current_prec, f = self.__q_expansion
        except AttributeError:
            current_prec = 0
            f = self.parent()._q_expansion_ring()(0, -1)

        if current_prec == prec:
            return f
        elif current_prec > prec:
            return f.add_bigoh(prec)
        else:
            f = self._compute_q_expansion(prec)
            self.__q_expansion = (prec, f)
            return f

    def atkin_lehner_eigenvalue(self, d=None):
        r"""
        Return the eigenvalue of the Atkin-Lehner operator W_d acting on self
        (which is either 1 or -1), or None if this form is not an eigenvector
        for this operator. If d is not given or is None, use d = the level.

        EXAMPLES::

            sage: sage.modular.modform.element.ModularForm_abstract.atkin_lehner_eigenvalue(CuspForms(2, 8).0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    # this function lives here so it is inherited by Newform (which does *not* derive from ModularFormElement)
    def cuspform_lseries(self, prec=53,
                         max_imaginary_part=0,
                         max_asymp_coeffs=40):
        r"""
        Return the L-series of the weight k cusp form
        f on `\Gamma_0(N)`.

        This actually returns an interface to Tim Dokchitser's program
        for computing with the L-series of the cusp form.

        INPUT:

        - ``prec`` - integer (bits precision)

        - ``max_imaginary_part`` - real number

        - ``max_asymp_coeffs`` - integer

        OUTPUT:

        The L-series of the cusp form.

        EXAMPLES::

           sage: f = CuspForms(2,8).newforms()[0]
           sage: L = f.cuspform_lseries()
           sage: L(1)
           0.0884317737041015
           sage: L(0.5)
           0.0296568512531983

        Consistency check with delta_lseries (which computes coefficients in pari)::

           sage: delta = CuspForms(1,12).0
           sage: L = delta.cuspform_lseries()
           sage: L(1)
           0.0374412812685155
           sage: L = delta_lseries()
           sage: L(1)
           0.0374412812685155

        We check that #5262 is fixed::

            sage: E=EllipticCurve('37b2')
            sage: h=Newforms(37)[1]
            sage: Lh = h.cuspform_lseries()
            sage: LE=E.lseries()
            sage: Lh(1), LE(1)
            (0.725681061936153, 0.725681061936153)
            sage: CuspForms(1, 30).0.cuspform_lseries().eps
            -1
        """
        if self.q_expansion().list()[0] !=0:
            raise TypeError,"f = %s is not a cusp form"%self
        from sage.lfunctions.all import Dokchitser
        key = (prec, max_imaginary_part, max_asymp_coeffs)
        l = self.weight()
        N = self.level()
        w = self.atkin_lehner_eigenvalue()
        if w is None:
            raise ValueError, "Form is not an eigenform for Atkin-Lehner"
        e = (-1)**(l/2)*w
        L = Dokchitser(conductor = N,
                       gammaV = [0,1],
                       weight = l,
                       eps = e,
                       prec = prec)
        s = 'coeff = %s;'%self.q_expansion(prec).list()
        L.init_coeffs('coeff[k+1]',pari_precode = s,
                      max_imaginary_part=max_imaginary_part,
                      max_asymp_coeffs=max_asymp_coeffs)
        L.check_functional_equation()
        L.rename('L-series associated to the cusp form %s'%self)
        return L

class Newform(ModularForm_abstract):
    def __init__(self, parent, component, names, check=True):
        r"""
        Initialize a Newform object.

        INPUT:

        - ``parent`` - An ambient cuspidal space of modular forms for
          which self is a newform.

        - ``component`` - A simple component of a cuspidal modular
          symbols space of any sign corresponding to this newform.

        - ``check`` - If check is ``True``, check that parent and
          component have the same weight, level, and character, that
          component has sign 1 and is simple, and that the types are
          correct on all inputs.

        EXAMPLES::

            sage: sage.modular.modform.element.Newform(CuspForms(11,2), ModularSymbols(11,2,sign=1).cuspidal_subspace(), 'a')
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)

            sage: f = Newforms(DirichletGroup(5).0, 7,names='a')[0]; f[2].trace(f.base_ring().base_field())
            -5*zeta4 - 5
        """
        if check:
            if not space.is_ModularFormsSpace(parent):
                raise TypeError, "parent must be a space of modular forms"
            if not is_ModularSymbolsSpace(component):
                raise TypeError, "component must be a space of modular symbols"
            if parent.group() != component.group():
                raise ValueError, "parent and component must be defined by the same congruence subgroup"
            if parent.weight() != component.weight():
                raise ValueError, "parent and component must have the same weight"
            if not component.is_cuspidal():
                raise ValueError, "component must be cuspidal"
            if not component.is_simple():
                raise ValueError, "component must be simple"
        extension_field = component.eigenvalue(1,name=names).parent()
        if extension_field != parent.base_ring(): # .degree() != 1 and rings.is_NumberField(extension_field):
            assert extension_field.base_field() == parent.base_ring()
            extension_field = parent.base_ring().extension(extension_field.relative_polynomial(), names=names)
        self.__name = names
        ModuleElement.__init__(self, parent.base_extend(extension_field))
        self.__modsym_space = component
        self.__hecke_eigenvalue_field = extension_field

    def _name(self):
        """
        Return the name of the generator of the Hecke eigenvalue field
        of self. Note that a name exists even when this field is QQ.

        EXAMPLES::

            sage: [ f._name() for f in Newforms(38,4,names='a') ]
            ['a0', 'a1', 'a2']
        """
        return self.__name

    def _compute_q_expansion(self, prec):
        """
        Return the q-expansion of self to precision prec.

        EXAMPLES::

            sage: forms = Newforms(31, 6, names='a')
            sage: forms[0]._compute_q_expansion(10)
            q + a0*q^2 + (5/704*a0^4 + 43/704*a0^3 - 61/88*a0^2 - 197/44*a0 + 717/88)*q^3 + (a0^2 - 32)*q^4 + (-31/352*a0^4 - 249/352*a0^3 + 111/22*a0^2 + 218/11*a0 - 2879/44)*q^5 + (-1/352*a0^4 - 79/352*a0^3 - 67/44*a0^2 + 13/22*a0 - 425/44)*q^6 + (17/88*a0^4 + 133/88*a0^3 - 405/44*a0^2 - 1005/22*a0 - 35/11)*q^7 + (a0^3 - 64*a0)*q^8 + (39/352*a0^4 + 441/352*a0^3 - 93/44*a0^2 - 441/22*a0 - 5293/44)*q^9 + O(q^10)
            sage: forms[0]._compute_q_expansion(15)
            q + a0*q^2 + (5/704*a0^4 + 43/704*a0^3 - 61/88*a0^2 - 197/44*a0 + 717/88)*q^3 + (a0^2 - 32)*q^4 + (-31/352*a0^4 - 249/352*a0^3 + 111/22*a0^2 + 218/11*a0 - 2879/44)*q^5 + (-1/352*a0^4 - 79/352*a0^3 - 67/44*a0^2 + 13/22*a0 - 425/44)*q^6 + (17/88*a0^4 + 133/88*a0^3 - 405/44*a0^2 - 1005/22*a0 - 35/11)*q^7 + (a0^3 - 64*a0)*q^8 + (39/352*a0^4 + 441/352*a0^3 - 93/44*a0^2 - 441/22*a0 - 5293/44)*q^9 + (15/176*a0^4 - 135/176*a0^3 - 185/11*a0^2 + 311/11*a0 + 2635/22)*q^10 + (-291/704*a0^4 - 3629/704*a0^3 + 1139/88*a0^2 + 10295/44*a0 - 21067/88)*q^11 + (-75/176*a0^4 - 645/176*a0^3 + 475/22*a0^2 + 1503/11*a0 - 5651/22)*q^12 + (207/704*a0^4 + 2977/704*a0^3 + 581/88*a0^2 - 3307/44*a0 - 35753/88)*q^13 + (-5/22*a0^4 + 39/11*a0^3 + 763/22*a0^2 - 2296/11*a0 - 2890/11)*q^14 + O(q^15)
        """
        return self.modular_symbols(1).q_eigenform(prec, names=self._name())

    def __eq__(self, other):
        """
        Return True if self equals other, and False otherwise.

        EXAMPLES::

            sage: f1, f2 = Newforms(17,4,names='a')
            sage: f1.__eq__(f1)
            True
            sage: f1.__eq__(f2)
            False
        """
        try:
            self._ensure_is_compatible(other)
        except Exception:
            return False
        if isinstance(other, Newform):
            if self.q_expansion(self.parent().sturm_bound()) == other.q_expansion(other.parent().sturm_bound()):
                return True
            else:
                return False
        if is_ModularFormElement(other):
            if self.element() == other.element():
                return True
            else:
                return False

    def __cmp__(self, other):
        """
        Compare self with other.

        EXAMPLES::

            sage: f1, f2 = Newforms(19,4,names='a')
            sage: f1.__cmp__(f1)
            0
            sage: f1.__cmp__(f2)
            -1
            sage: f2.__cmp__(f1)
            -1
        """
        try:
            self._ensure_is_compatible(other)
        except Exception:
            return self.parent().__cmp__(other.parent())
        if isinstance(other, Newform):
            if self.q_expansion(self.parent().sturm_bound()) == other.q_expansion(other.parent().sturm_bound()):
                return 0
            else:
                return -1
        if is_ModularFormElement(other):
            if self.element() == other.element():
                return 0
            else:
                return -1

    def abelian_variety(self):
        """
        Return the abelian variety associated to self.

        EXAMPLES::

            sage: Newforms(14,2)[0]
            q - q^2 - 2*q^3 + q^4 + O(q^6)
            sage: Newforms(14,2)[0].abelian_variety()
            Newform abelian subvariety 14a of dimension 1 of J0(14)
        """
        try:
            return self.__abelian_variety
        except AttributeError:
            from sage.modular.abvar.abvar_newform import ModularAbelianVariety_newform
            self.__abelian_variety = ModularAbelianVariety_newform(self)
            return self.__abelian_variety

    def hecke_eigenvalue_field(self):
        r"""
        Return the field generated over the rationals by the
        coefficients of this newform.

        EXAMPLES::

            sage: ls = Newforms(35, 2, names='a') ; ls
            [q + q^3 - 2*q^4 - q^5 + O(q^6),
            q + a1*q^2 + (-a1 - 1)*q^3 + (-a1 + 2)*q^4 + q^5 + O(q^6)]
            sage: ls[0].hecke_eigenvalue_field()
            Rational Field
            sage: ls[1].hecke_eigenvalue_field()
            Number Field in a1 with defining polynomial x^2 + x - 4
        """
        return self.__hecke_eigenvalue_field

    def _compute(self, X):
        """
        Compute the coefficients of `q^n` of the power series of self,
        for `n` in the list `X`.  The results are not cached.  (Use
        coefficients for cached results).

        EXAMPLES::

            sage: f = Newforms(39,4,names='a')[1] ; f
            q + a1*q^2 - 3*q^3 + (2*a1 + 5)*q^4 + (-2*a1 + 14)*q^5 + O(q^6)
            sage: f._compute([2,3,7])
            [alpha, -3, -2*alpha + 2]
            sage: f._compute([])
            []
        """
        M = self.modular_symbols(1)
        return [M.eigenvalue(x) for x in X]

    def element(self):
        """
        Find an element of the ambient space of modular forms which
        represents this newform.

        .. note::

           This can be quite expensive. Also, the polynomial defining
           the field of Hecke eigenvalues should be considered random,
           since it is generated by a random sum of Hecke
           operators. (The field itself is not random, of course.)

        EXAMPLES::

            sage: ls = Newforms(38,4,names='a')
            sage: ls[0]
            q - 2*q^2 - 2*q^3 + 4*q^4 - 9*q^5 + O(q^6)
            sage: ls # random
            [q - 2*q^2 - 2*q^3 + 4*q^4 - 9*q^5 + O(q^6),
            q - 2*q^2 + (-a1 - 2)*q^3 + 4*q^4 + (2*a1 + 10)*q^5 + O(q^6),
            q + 2*q^2 + (1/2*a2 - 1)*q^3 + 4*q^4 + (-3/2*a2 + 12)*q^5 + O(q^6)]
            sage: type(ls[0])
            <class 'sage.modular.modform.element.Newform'>
            sage: ls[2][3].minpoly()
            x^2 - 9*x + 2
            sage: ls2 = [ x.element() for x in ls ]
            sage: ls2 # random
            [q - 2*q^2 - 2*q^3 + 4*q^4 - 9*q^5 + O(q^6),
            q - 2*q^2 + (-a1 - 2)*q^3 + 4*q^4 + (2*a1 + 10)*q^5 + O(q^6),
            q + 2*q^2 + (1/2*a2 - 1)*q^3 + 4*q^4 + (-3/2*a2 + 12)*q^5 + O(q^6)]
            sage: type(ls2[0])
            <class 'sage.modular.modform.element.ModularFormElement'>
            sage: ls2[2][3].minpoly()
            x^2 - 9*x + 2
        """
        S = self.parent()
        return S(self.q_expansion(S.sturm_bound()))

    def modular_symbols(self, sign=0):
        """
        Return the subspace with the specified sign of the space of
        modular symbols corresponding to this newform.

        EXAMPLES::

            sage: f = Newforms(18,4)[0]
            sage: f.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_0(18) of weight 4 with sign 0 over Rational Field
            sage: f.modular_symbols(1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 11 for Gamma_0(18) of weight 4 with sign 1 over Rational Field
        """
        return self.__modsym_space.modular_symbols_of_sign(sign)

    def _defining_modular_symbols(self):
        """
        Return the modular symbols space corresponding to self.

        EXAMPLES::

            sage: Newforms(43,2,names='a')
            [q - 2*q^2 - 2*q^3 + 2*q^4 - 4*q^5 + O(q^6),
            q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6)]
            sage: [ x._defining_modular_symbols() for x in Newforms(43,2,names='a') ]
            [Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field]
            sage: ModularSymbols(43,2,sign=1).cuspidal_subspace().new_subspace().decomposition()
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
            ]
        """
        return self.__modsym_space

    def number(self):
        """
        Return the index of this space in the list of simple, new,
        cuspidal subspaces of the full space of modular symbols for
        this weight and level.

        EXAMPLES::

            sage: Newforms(43, 2, names='a')[1].number()
            1
        """
        return self._defining_modular_symbols().ambient().cuspidal_subspace().new_subspace().decomposition().index(self._defining_modular_symbols())

    def __nonzero__(self):
        """
        Return True, as newforms are never zero.

        EXAMPLES::

            sage: Newforms(14,2)[0].__nonzero__()
            True
        """
        return True

    def character(self):
        r"""
        The nebentypus character of this newform (as a Dirichlet character with
        values in the field of Hecke eigenvalues of the form).

        EXAMPLES::

            sage: Newforms(Gamma1(7), 4,names='a')[1].character()
            Dirichlet character modulo 7 of conductor 7 mapping 3 |--> 1/2*a1
            sage: chi = DirichletGroup(3).0; Newforms(chi, 7)[0].character() == chi
            True
        """
        return self._defining_modular_symbols().q_eigenform_character(self._name())


    def atkin_lehner_eigenvalue(self, d=None):
        r"""
        Return the eigenvalue of the Atkin-Lehner operator W_d acting on this newform
        (which is either 1 or -1). A ValueError will be raised if the character
        of this form is not either trivial or quadratic. If d is not given or
        is None, then d defaults to the level of self.

        EXAMPLE::

            sage: [x.atkin_lehner_eigenvalue() for x in ModularForms(53).newforms('a')]
            [1, -1]
            sage: CuspForms(DirichletGroup(5).0, 5).newforms()[0].atkin_lehner_eigenvalue()
            Traceback (most recent call last):
            ...
            ValueError: Atkin-Lehner only leaves space invariant when character is trivial or quadratic.  In general it sends M_k(chi) to M_k(1/chi)
        """
        return self.modular_symbols(sign=1).atkin_lehner_operator(d).matrix()[0,0]

class ModularFormElement(ModularForm_abstract, element.HeckeModuleElement):
    def __init__(self, parent, x, check=True):
        r"""
        An element of a space of modular forms.

        INPUT:

        - ``parent`` - ModularForms (an ambient space of modular forms)

        - ``x`` - a vector on the basis for parent

        - ``check`` - if check is ``True``, check the types of the
          inputs.

        OUTPUT:

        - ``ModularFormElement`` - a modular form

        EXAMPLES::

            sage: M = ModularForms(Gamma0(11),2)
            sage: f = M.0
            sage: f.parent()
            Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
        """
        if not isinstance(parent, space.ModularFormsSpace):
            raise TypeError, "First argument must be an ambient space of modular forms."
        element.HeckeModuleElement.__init__(self, parent, x)

    def _compute_q_expansion(self, prec):
        """
        Computes the q-expansion of self to precision prec.

        EXAMPLES::

            sage: f = EllipticCurve('37a').modular_form()
            sage: f.q_expansion() ## indirect doctest
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)

            sage: f._compute_q_expansion(10)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + O(q^10)
        """
        return self.parent()._q_expansion(element = self.element(), prec=prec)

    def _add_(self, other):
        """
        Add self to other.

        EXAMPLES::

            sage: f = ModularForms(DirichletGroup(17).0^2,2).2
            sage: g = ModularForms(DirichletGroup(17).0^2,2).1
            sage: f
            q + (-zeta8^2 + 2)*q^2 + (zeta8 + 3)*q^3 + (-2*zeta8^2 + 3)*q^4 + (-zeta8 + 5)*q^5 + O(q^6)

            sage: g
            1 + (-14/73*zeta8^3 + 57/73*zeta8^2 + 13/73*zeta8 - 6/73)*q^2 + (-90/73*zeta8^3 + 64/73*zeta8^2 - 52/73*zeta8 + 24/73)*q^3 + (-81/73*zeta8^3 + 189/73*zeta8^2 - 3/73*zeta8 + 153/73)*q^4 + (72/73*zeta8^3 + 124/73*zeta8^2 + 100/73*zeta8 + 156/73)*q^5 + O(q^6)

            sage: f+g ## indirect doctest
            1 + q + (-14/73*zeta8^3 - 16/73*zeta8^2 + 13/73*zeta8 + 140/73)*q^2 + (-90/73*zeta8^3 + 64/73*zeta8^2 + 21/73*zeta8 + 243/73)*q^3 + (-81/73*zeta8^3 + 43/73*zeta8^2 - 3/73*zeta8 + 372/73)*q^4 + (72/73*zeta8^3 + 124/73*zeta8^2 + 27/73*zeta8 + 521/73)*q^5 + O(q^6)
        """
        return ModularFormElement(self.parent(), self.element() + other.element())

    def __mul__(self, other):
        r"""
        Calculate the product self * other.

        This tries to determine the
        characters of self and other, in order to avoid having to compute a
        (potentially very large) Gamma1 space. Note that this might lead to
        a modular form that is defined with respect to a larger subgroup than
        the factors are.

        An example with character::

            sage: f = ModularForms(DirichletGroup(3).0, 3).0
            sage: f * f
            1 + 108*q^2 + 144*q^3 + 2916*q^4 + 8640*q^5 + O(q^6)
            sage: (f*f).parent()
            Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(3) of weight 6 over Rational Field
            sage: (f*f*f).parent()
            Modular Forms space of dimension 4, character [-1] and weight 9 over Rational Field

        An example where the character is computed on-the-fly::

            sage: f = ModularForms(Gamma1(3), 5).0
            sage: f*f
            1 - 180*q^2 - 480*q^3 + 8100*q^4 + 35712*q^5 + O(q^6)
            sage: (f*f).parent()
            Modular Forms space of dimension 4 for Congruence Subgroup Gamma0(3) of weight 10 over Rational Field

            sage: f = ModularForms(Gamma1(3), 7).0
            sage: f*f
            q^2 - 54*q^4 + 128*q^5 + O(q^6)
            sage: (f*f).parent()
            Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(3) of weight 14 over Rational Field

        An example with no character::

            sage: f = ModularForms(Gamma1(5), 2).0
            sage: f*f
            1 + 120*q^3 - 240*q^4 + 480*q^5 + O(q^6)
            sage: (f*f).parent()
            Modular Forms space of dimension 5 for Congruence Subgroup Gamma1(5) of weight 4 over Rational Field

        TESTS:

        This shows that the issue at trac ticket #7548 is fixed::

            sage: M = CuspForms(Gamma0(5*3^2), 2)
            sage: f = M.basis()[0]
            sage: 2*f
            2*q - 2*q^4 + O(q^6)
            sage: f*2
            2*q - 2*q^4 + O(q^6)
        """

        # boring case: scalar multiplication
        if not isinstance(other, ModularFormElement):
            return element.HeckeModuleElement.__mul__(self, other)

        # first ensure the levels are equal
        if self.level() != other.level():
            raise NotImplementedError, "Cannot multiply forms of different levels"

        # find out about characters
        try:
            eps1 = self.character()
            verbose("character of left is %s" % eps1)
            eps2 = other.character()
            verbose("character of right is %s" % eps2)
            newchar = eps1 * eps2
            verbose("character of product is %s" % newchar)
        except (NotImplementedError, ValueError):
            newchar = None
            verbose("character of product not determined")

        # now do the math
        from constructor import ModularForms
        if newchar is not None:
            verbose("creating a parent with char")
            newparent = ModularForms(newchar, self.weight() + other.weight(), base_ring = newchar.base_ring())
            verbose("parent is %s" % newparent)
        else:
            newparent = ModularForms(self.group(), self.weight() + other.weight(), base_ring = rings.ZZ)
        m = newparent.sturm_bound()
        newqexp = self.qexp(m) * other.qexp(m)

        return newparent.base_extend(newqexp.base_ring())(newqexp)

    def modform_lseries(self, prec=53,
                        max_imaginary_part=0,
                        max_asymp_coeffs=40):
        r"""
        Return the L-series of the weight `k` modular form
        `f` on `\mathrm{SL}_2(\ZZ)`.

        This actually returns an interface to Tim Dokchitser's program
        for computing with the L-series of the modular form.

        INPUT:

        - ``prec`` - integer (bits precision)

        - ``max_imaginary_part`` - real number

        - ``max_asymp_coeffs`` - integer

        OUTPUT:

        The L-series of the modular form.

        EXAMPLES:

        We compute with the L-series of the Eisenstein series `E_4`::

           sage: f = ModularForms(1,4).0
           sage: L = f.modform_lseries()
           sage: L(1)
           -0.0304484570583933
        """
        a = self.q_expansion(prec).list()
        if a[0] == 0:
            raise TypeError,"f = %s is a cusp form; please use f.cuspform_lseries() instead!"%self
        if self.level() != 1:
            raise TypeError, "f = %s is not a modular form for SL_2(Z)"%self
        from sage.lfunctions.all import Dokchitser
        key = (prec, max_imaginary_part, max_asymp_coeffs)
        l = self.weight()
        L = Dokchitser(conductor = 1,
                       gammaV = [0,1],
                       weight = l,
                       eps = (-1)**l,
                       poles = [l],
                       prec = prec)
        b = a[1]
        for i in range(len(a)):    ##to renormalize so that coefficient of q is 1
            a[i] =(1/b)*a[i]
        s = 'coeff = %s;'%a
        L.init_coeffs('coeff[k+1]',pari_precode = s,
                      max_imaginary_part=max_imaginary_part,
                      max_asymp_coeffs=max_asymp_coeffs)
        L.check_functional_equation()
        L.rename('L-series associated to the weight %s modular form on SL_2(Z)'%l)
        return L

    def atkin_lehner_eigenvalue(self, d=None):
        r"""
        Return the eigenvalue of the Atkin-Lehner operator W_d acting on this
        modular form (which is either 1 or -1), or None if this form is not an
        eigenvector for this operator.

        EXAMPLE::

             sage: CuspForms(1, 30).0.atkin_lehner_eigenvalue()
             1
             sage: CuspForms(2, 8).0.atkin_lehner_eigenvalue()
             Traceback (most recent call last):
             ...
             NotImplementedError: Don't know how to compute Atkin-Lehner matrix acting on this space (try using a newform constructor instead)
        """
        try:
            f = self.parent().atkin_lehner_operator(d)(self)
        except NotImplementedError:
            raise NotImplementedError, "Don't know how to compute Atkin-Lehner matrix acting on this space" \
                + " (try using a newform constructor instead)"
        if f == self:
            return 1
        elif f == -self:
            return -1
        else:
            return None

class ModularFormElement_elliptic_curve(ModularFormElement):
    """
    A modular form attached to an elliptic curve.
    """
    def __init__(self, parent, E):
        """
        Modular form attached to an elliptic curve as an element
        of a space of modular forms.

        EXAMPLES::

            sage: E = EllipticCurve('5077a')
            sage: f = E.modular_form()
            sage: f
            q - 2*q^2 - 3*q^3 + 2*q^4 - 4*q^5 + O(q^6)
            sage: f.q_expansion(10)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 4*q^5 + 6*q^6 - 4*q^7 + 6*q^9 + O(q^10)
            sage: f.parent()
            Modular Forms space of dimension 423 for Congruence Subgroup Gamma0(5077) of weight 2 over Rational Field

            sage: E = EllipticCurve('37a')
            sage: f = E.modular_form() ; f
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)
            sage: f == loads(dumps(f))
            True
        """
        ModularFormElement.__init__(self, parent, None)
##                                    parent.find_in_space( E.q_expansion(parent.hecke_bound()) ))
        self.__E = E


    def elliptic_curve(self):
        """
        Return elliptic curve associated to self.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: f = E.modular_form()
            sage: f.elliptic_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: f.elliptic_curve() is E
            True
        """
        return self.__E

    def _compute_element(self):
        """
        Compute self as a linear combination of the basis elements
        of parent.

        EXAMPLES::

            sage: EllipticCurve('11a1').modular_form()._compute_element()
            (1, 0)
            sage: EllipticCurve('389a1').modular_form()._compute_element()
            (1, -2, -2, 2, -3, 4, -5, 0, 1, 6, -4, -4, -3, 10, 6, -4, -6, -2, 5, -6, 10, 8, -4, 0, 4, 6, 4, -10, -6, -12, 4, 8, 0)
        """
        M = self.parent()
        S = M.cuspidal_subspace()
##        return S.find_in_space( self.__E.q_expansion( S.q_expansion_basis()[0].prec() ) ) + [0] * ( M.dimension() - S.dimension() )
        return vector(S.find_in_space( self.__E.q_expansion( S.sturm_bound() ) ) + [0] * ( M.dimension() - S.dimension() ))

    def _compute_q_expansion(self, prec):
        r"""
        The `q`-expansion of the modular form to precision `O(q^\text{prec})`.
        This function takes one argument, which is the integer prec.

        EXAMPLES::

            sage: E = EllipticCurve('11a') ; f = E.modular_form()
            sage: f._compute_q_expansion(10)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)

            sage: f._compute_q_expansion(30)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 - 2*q^10 + q^11 - 2*q^12 + 4*q^13 + 4*q^14 - q^15 - 4*q^16 - 2*q^17 + 4*q^18 + 2*q^20 + 2*q^21 - 2*q^22 - q^23 - 4*q^25 - 8*q^26 + 5*q^27 - 4*q^28 + O(q^30)

            sage: f._compute_q_expansion(10)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)
        """
        return self.__E.q_expansion(prec)

    def atkin_lehner_eigenvalue(self, d=None):
        r"""
        Calculate the eigenvalue of the Atkin-Lehner operator W_d acting on
        this form. If d is None, default to the level of the form. As this form
        is attached to an elliptic curve, we can read this off from the root
        number of the curve if d is the level.

        EXAMPLE::

            sage: EllipticCurve('57a1').newform().atkin_lehner_eigenvalue()
            1
            sage: EllipticCurve('57b1').newform().atkin_lehner_eigenvalue()
            -1
            sage: EllipticCurve('57b1').newform().atkin_lehner_eigenvalue(19)
            1
        """
        if d is None:
            return -self.__E.root_number()
        else:
            return self.__E.modular_symbol_space().atkin_lehner_operator(d).matrix()[0,0]

######################################################################

class EisensteinSeries(ModularFormElement):
    """
    An Eisenstein series.

    EXAMPLES::

        sage: E = EisensteinForms(1,12)
        sage: E.eisenstein_series()
        [
        691/65520 + q + 2049*q^2 + 177148*q^3 + 4196353*q^4 + 48828126*q^5 + O(q^6)
        ]
        sage: E = EisensteinForms(11,2)
        sage: E.eisenstein_series()
        [
        5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
        ]
        sage: E = EisensteinForms(Gamma1(7),2)
        sage: E.set_precision(4)
        sage: E.eisenstein_series()
        [
        1/4 + q + 3*q^2 + 4*q^3 + O(q^4),
        1/7*zeta6 - 3/7 + q + (-2*zeta6 + 1)*q^2 + (3*zeta6 - 2)*q^3 + O(q^4),
        q + (-zeta6 + 2)*q^2 + (zeta6 + 2)*q^3 + O(q^4),
        -1/7*zeta6 - 2/7 + q + (2*zeta6 - 1)*q^2 + (-3*zeta6 + 1)*q^3 + O(q^4),
        q + (zeta6 + 1)*q^2 + (-zeta6 + 3)*q^3 + O(q^4)
        ]
    """
    def __init__(self, parent, vector, t, chi, psi):
        """
        An Eisenstein series.

        EXAMPLES::

            sage: E = EisensteinForms(1,12) ## indirect doctest
            sage: E.eisenstein_series()
            [
            691/65520 + q + 2049*q^2 + 177148*q^3 + 4196353*q^4 + 48828126*q^5 + O(q^6)
            ]
            sage: E = EisensteinForms(11,2)
            sage: E.eisenstein_series()
            [
            5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
            ]
            sage: E = EisensteinForms(Gamma1(7),2)
            sage: E.set_precision(4)
            sage: E.eisenstein_series()
            [
            1/4 + q + 3*q^2 + 4*q^3 + O(q^4),
            1/7*zeta6 - 3/7 + q + (-2*zeta6 + 1)*q^2 + (3*zeta6 - 2)*q^3 + O(q^4),
            q + (-zeta6 + 2)*q^2 + (zeta6 + 2)*q^3 + O(q^4),
            -1/7*zeta6 - 2/7 + q + (2*zeta6 - 1)*q^2 + (-3*zeta6 + 1)*q^3 + O(q^4),
            q + (zeta6 + 1)*q^2 + (-zeta6 + 3)*q^3 + O(q^4)
            ]
        """
        N = parent.level()
        K = parent.base_ring()
        if chi.parent().modulus() != N or psi.parent().modulus() != N:
            raise ArithmeticError, "Incompatible moduli"
        if chi.parent().base_ring() != K or psi.parent().base_ring() != K:
            raise ArithmeticError, "Incompatible base rings"
        t = int(t)
        #if not isinstance(t, int): raise TypeError, "weight must be an int"
        if parent.weight() == 2 and chi.is_trivial() \
               and psi.is_trivial() and t==1:
            raise ArithmeticError, "If chi and psi are trivial and k=2, then t must be >1."
        ModularFormElement.__init__(self, parent, vector)
        self.__chi = chi
        self.__psi = psi
        self.__t   = t

    def _compute_q_expansion(self, prec=None):
        """
        Compute the q-expansion of self to precision prec.

        EXAMPLES::

            sage: EisensteinForms(11,2).eisenstein_series()[0]._compute_q_expansion(10)
            5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + 12*q^6 + 8*q^7 + 15*q^8 + 13*q^9 + O(q^10)
        """
        if prec is None:
            prec = self.parent().prec()
        F = self._compute(range(prec))
        R = self.parent()._q_expansion_ring()
        return R(F, prec)

    def _compute(self, X):
        r"""
        Compute the coefficients of `q^n` of the power series of self,
        for `n` in the list `X`.  The results are not cached.  (Use
        coefficients for cached results).

        EXAMPLES::

            sage: e = DirichletGroup(11).gen()
            sage: f = EisensteinForms(e, 3).eisenstein_series()[0]
            sage: f._compute([3,4,5])
            [-9*zeta10^3 + 1,
             16*zeta10^2 + 4*zeta10 + 1,
             25*zeta10^3 - 25*zeta10^2 + 25*zeta10 - 24]

        """
        if self.weight() == 2 and (self.__chi.is_trivial() and self.__psi.is_trivial()):
            return self.__compute_weight2_trivial_character(X)
        else: # general case
            return self.__compute_general_case(X)

    def __compute_weight2_trivial_character(self, X):
        r"""
        Compute coefficients for self an Eisenstein series of the form
        `E_2 - t*E_2(q^t)`. Computes `a_n` for each `n \in X`.

        EXAMPLES::

            sage: EisensteinForms(14,2).eisenstein_series()[0]._EisensteinSeries__compute_weight2_trivial_character([0])
            [1/24]
            sage: EisensteinForms(14,2).eisenstein_series()[0]._EisensteinSeries__compute_weight2_trivial_character([0,4,11,38])
            [1/24, 1, 12, 20]
        """
        F = self.base_ring()
        v = []
        t = self.__t
        for n in X:
            if n < 0:
                pass
            elif n == 0:
                v.append(F(t-1)/F(24))
            else:
                an = rings.sigma(n,1)
                if n%t==0:
                    an -= t*rings.sigma(n/t,1)
                v.append(an)
        return v

    def __compute_general_case(self, X):
        """
        Returns the list coefficients of `q^n` of the power series of self,
        for `n` in the list `X`.  The results are not cached.  (Use
        coefficients for cached results).

        General case (except weight 2, trivial character, where this
        is wrong!)  `\chi` is a primitive character of conductor `L`
        `\psi` is a primitive character of conductor `M` We have
        `MLt \mid N`, and

        .. math::

          E_k(chi,psi,t) =
           c_0 + sum_{m \geq 1}[sum_{n|m} psi(n) * chi(m/n) * n^(k-1)] q^{mt},

        with `c_0=0` if `L>1`, and `c_0=L(1-k,psi)/2` if `L=1` (that
        second `L` is an `L`-function `L`).

        EXAMPLES::

            sage: e = DirichletGroup(11).gen()
            sage: f = EisensteinForms(e, 3).eisenstein_series()[0]
            sage: f._EisensteinSeries__compute_general_case([1])
            [1]
            sage: f._EisensteinSeries__compute_general_case([2])
            [4*zeta10 + 1]
            sage: f._EisensteinSeries__compute_general_case([0,1,2])
            [15/11*zeta10^3 - 9/11*zeta10^2 - 26/11*zeta10 - 10/11, 1, 4*zeta10 + 1]
        """
        c0, chi, psi, K, n, t, L, M = self.__defining_parameters()
        zero = K(0)
        k = self.weight()
        v = []
        for i in X:
            if i==0:
                v.append(c0)
                continue
            if i%t != 0:
                v.append(zero)
            else:
                m = i//t
                v.append(sum([psi(n)*chi(m/n)*n**(k-1) for \
                               n in rings.divisors(m)]))
        return v

    def __defining_parameters(self):
        """
        Return defining parameters for self.

        EXAMPLES::

            sage: EisensteinForms(11,2).eisenstein_series()[0]._EisensteinSeries__defining_parameters()
            (-1/24, Dirichlet character modulo 1 of conductor 1 mapping 0 |--> 1, Dirichlet character modulo 1 of conductor 1 mapping 0 |--> 1, Rational Field, 2, 11, 1, 1)
        """
        try:
            return self.__defining_params
        except AttributeError:
            chi = self.__chi.primitive_character()
            psi = self.__psi.primitive_character()
            k = self.weight()
            t = self.__t
            L = chi.conductor()
            M = psi.conductor()
            K = chi.base_ring()
            n = K.zeta_order()
            if L == 1:
                c0 = K(-psi.bernoulli(k))/K(2*k)
            else:
                c0 = K(0)
            self.__defining_params = (c0, chi, psi, K, n, t, L, M)
        return self.__defining_params

    def chi(self):
        """
        Return the parameter chi associated to self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].chi()
            Dirichlet character modulo 17 of conductor 17 mapping 3 |--> zeta16
        """
        return self.__chi

    def psi(self):
        """
        Return the parameter psi associated to self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].psi()
             Dirichlet character modulo 17 of conductor 1 mapping 3 |--> 1
        """
        return self.__psi

    def t(self):
        """
        Return the parameter t associated to self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].t()
            1
        """
        return self.__t

    def parameters(self):
        """
        Return chi, psi, and t, which are the defining parameters of self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].parameters()
            (Dirichlet character modulo 17 of conductor 17 mapping 3 |--> zeta16, Dirichlet character modulo 17 of conductor 1 mapping 3 |--> 1, 1)
        """
        return self.__chi, self.__psi, self.__t

    def L(self):
        """
        Return the conductor of self.chi().

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].L()
            17
        """
        return self.__chi.conductor()

    def M(self):
        """
        Return the conductor of self.psi().

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].M()
            1
        """
        return self.__psi.conductor()

    def character(self):
        """
        Return the character associated to self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].character()
            Dirichlet character modulo 17 of conductor 17 mapping 3 |--> zeta16

            sage: chi = DirichletGroup(7)[4]
            sage: E = EisensteinForms(chi).eisenstein_series() ; E
            [
            -1/7*zeta6 - 2/7 + q + (2*zeta6 - 1)*q^2 + (-3*zeta6 + 1)*q^3 + (-2*zeta6 - 1)*q^4 + (5*zeta6 - 4)*q^5 + O(q^6),
            q + (zeta6 + 1)*q^2 + (-zeta6 + 3)*q^3 + (zeta6 + 2)*q^4 + (zeta6 + 4)*q^5 + O(q^6)
            ]
            sage: E[0].character() == chi
            True
            sage: E[1].character() == chi
            True

        TESTS::

            sage: [ [ f.character() == chi for f in EisensteinForms(chi).eisenstein_series() ] for chi in DirichletGroup(17) ]
            [[True], [], [True, True], [], [True, True], [], [True, True], [], [True, True], [], [True, True], [], [True, True], [], [True, True], []]

            sage: [ [ f.character() == chi for f in EisensteinForms(chi).eisenstein_series() ] for chi in DirichletGroup(16) ]
            [[True, True, True, True, True], [], [True, True], [], [True, True, True, True], [], [True, True], []]
        """
        try:
            return self.__character
        except AttributeError:
            self.__character = self.__chi * self.__psi
        return self.__character

    def new_level(self):
        """
        Return level at which self is new.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].level()
            17
            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].new_level()
            17
            sage: [ [x.level(), x.new_level()] for x in EisensteinForms(DirichletGroup(60).0^2,2).eisenstein_series() ]
            [[60, 2], [60, 3], [60, 2], [60, 5], [60, 2], [60, 2], [60, 2], [60, 3], [60, 2], [60, 2], [60, 2]]
        """
        if self.__chi.is_trivial() and self.__psi.is_trivial() and self.weight() == 2:
            return rings.factor(self.__t)[0][0]
        return self.L()*self.M()



