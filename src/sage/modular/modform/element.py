"""
Elements of modular forms spaces.
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import space
import sage.modular.hecke.element as element
import sage.rings.all as rings

def is_ModularFormElement(x):
    """
    Return True if x is a modular form.

    EXAMPLES:
        sage: is_ModularFormElement(5)
        False
        sage: is_ModularFormElement(ModularForms(11).0)
        True
    """
    return isinstance(x, ModularFormElement)

class ModularFormElement(element.HeckeModuleElement):
    """
    An element of a space of modular forms.
    """
    def __init__(self, parent, x):
        """
        INPUT:
            parent -- ModularForms (an ambient space of modular forms)
            x -- a vector on the basis for parent

        OUTPUT:
            ModularFormElement -- a modular form

        EXAMPLES:
            sage: M = ModularForms(Gamma0(11),2)
            sage: f = M.0
            sage: f.parent()
            Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
        """
        if not isinstance(parent, space.ModularFormsSpace):
            raise TypeError, "First argument must be an ambient space of modular forms."
        element.HeckeModuleElement.__init__(self, parent, x)

    def __ensure_is_compatible(self, other):
        if not isinstance(other, ModularFormElement):
            raise TypeError, "Second argument must be a modular form."
        if self.ambient_module() != other.ambient_module():
            raise ArithmeticError, "Modular forms must be in the same ambient space."

    def __call__(self, x, prec=None):
        """
        Evaluate the q-expansion of this modular form at x.
        """
        return self.q_expansion(prec)(x)

    def _add_(self, other):
        return ModularFormElement(self.ambient_module(), self.element() + other.element())

    def __cmp__(self, other):
        self.__ensure_is_compatible(other)
        return self.element() == other.element()

    def coefficients(self, X):
        """
        The coefficients a_n of self, for integers n>=0 in the list X.

        This function caches the results of the compute function.
        """
        try:
            self.__coefficients
        except AttributeError:
            self.__coefficients = {}
        Y = [n for n in X   if    not (n in self.__coefficients.keys())]
        v = self._compute(Y)
        for i in range(len(v)):
            self.__coefficients[X[i]] = v[i]
        return v

##     def _compute(self, X):
##         r"""
##         Compute the coefficients $a_n$ of self, for integers $n \geq
##         0$ in the list $X$.

##         NOTES: The results need not be cached; use the coefficients
##         method instead for cached results.
##         """
##         return self.parent()._compute_coefficients(self.element(), X)

    def __getitem__(self, n):
        return self.q_expansion(n+1)[int(n)]

    def __getslice__(self, i, j):
        return self.q_expansion(j+1)[int(i):int(j)]

    def padded_list(self, n):
        return self.q_expansion(n).padded_list(n)


    def _repr_(self):
        return str(self.q_expansion())

    def _latex_(self):
        return self.q_expansion()._latex_()

    def base_ring(self):
        return self.parent().base_ring()

    def character(self):
        chi = self.parent().character()
        if chi is None:
            raise NotImplementedError, "Determination of character in this " + \
                  "case not implemented yet."
        return chi

    def __nonzero__(self):
        return not self.element().is_zero()

    def level(self):
        return self.parent().level()

    def prec(self):
        try:
            self.__qexp
        except AttributeError:
            return self.parent().prec()
        return self.__qexp.prec()

    def q_expansion(self, prec=None):
        r"""
        The $q$-expansion of the modular form to precision $O(q^\text{prec})$.
        This function takes one argument, which is the integer prec.

        EXAMPLES:
        We compute the cusp form $\Delta$.
            sage: delta = CuspForms(1,12).0
            sage: delta.q_expansion()
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)

        We compute the $q$-expansion of one of the cusp forms of level 23:
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
            prec = self.prec()
        prec = rings.Integer(prec)
        if prec < 0:
            raise ValueError, "prec (=%s) must be at least 0"%prec
        try:
            current_prec, f = self.__q_expansion
        except AttributeError:
            zero = self.parent()._q_expansion_ring()(0, -1)
            current_prec, f = 0, zero
        if current_prec == prec:
            return f
        elif current_prec > prec:
            return f.add_bigoh(prec)
        f = self._compute_q_expansion(prec)
        self.__q_expansion = (prec, f)
        return f

    def _compute_q_expansion(self, prec=None):
        return self.parent()._q_expansion(element = self.element(), prec=prec)

    def qexp(self, prec=None):
        """
        Same as self.q_expansion(prec).
        """
        return self.q_expansion(prec)

    def weight(self):
        return self.parent().weight()

    def valuation(self):
        if self.__valuation != None:
            return self.__valuation
        v = self.qexp().valuation()
        if not (v is rings.infinity):
            self.__valuation = v
            return v
        v = self.qexp(self.parent().sturm_bound()).valuation()
        self.__valuation = v
        return v



class ModularFormElement_elliptic_curve(ModularFormElement):
    """
    A modular form attached to an elliptic curve.

    EXAMPLES:
        sage: E = EllipticCurve('5077a')
        sage: f = E.modular_form()
        sage: f
        q - 2*q^2 - 3*q^3 + 2*q^4 - 4*q^5 + O(q^6)
        sage: f.q_expansion(10)
        q - 2*q^2 - 3*q^3 + 2*q^4 - 4*q^5 + 6*q^6 - 4*q^7 + 6*q^9 + O(q^10)
        sage: f.parent()
        Modular Forms space of dimension 423 for Congruence Subgroup Gamma0(5077) of weight 2 over Rational Field
    """
    def __init__(self, parent, E, x=None):
        """
        Modular form attached to an elliptic curve as an element
        of a space of modular forms.
        """
        ModularFormElement.__init__(self, parent, x)
        self.__E = E

    def elliptic_curve(self):
        return self.__E

    def _compute_element(self):
        # TODO
        raise NotImplementedError, "todo -- compute q-exp, find element of space, etc."

    def _compute_q_expansion(self, prec=None):
        r"""
        The $q$-expansion of the modular form to precision $O(q^\text{prec})$.
        This function takes one argument, which is the integer prec.

        EXAMPLES:
        """
        if prec is None:
            prec = self.parent().prec()
        return self.__E.q_expansion(prec)

######################################################################

class EisensteinSeries(ModularFormElement):
    """
    An Eisenstein series.

    EXAMPLES:
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

    def _compute_q_expansion(self, prec):
        F = self._compute(range(prec))
        R = self.parent()._q_expansion_ring()
        return R(F, prec)

    def _compute(self, X):
        """
        Compute the coefficients of $q^n$ of the power series of self,
        for $n$ in the list $X$.  The results are not cached.  (Use
        coefficients for cached results).
        """
        if self.weight() == 2 and (self.__chi.is_trivial() and self.__psi.is_trivial()):
            return self.__compute_weight2_trivial_character(X)
        else: # general case
            return self.__compute_general_case(X)

    def __compute_weight2_trivial_character(self, X):
        """
        Compute $E_2 - t*E_2(q^t)$.
        """
        F = self.base_ring()
        v = []
        t = self.__t
        for n in X:
            if n <= 0:
                v.append(F(t-1)/F(24))
            else:
                an = rings.sigma(n,1)
                if n%t==0:
                    an -= t*rings.sigma(n/t,1)
                v.append(an)
        return v

    def __compute_general_case(self, X):
        """
        General case (except weight 2, trivial character, where this is wrong!)
        $\chi$ is a primitive character of conductor $L$
        $\psi$ is a primitive character of conductor $M$
        We have $MLt \mid N$, and
        $$
          E_k(chi,psi,t) =
           c_0 + sum_{m \geq 1}[sum_{n|m} psi(n) * chi(m/n) * n^(k-1)] q^{mt},
        $$
        with $c_0=0$ if $L>1$,
         and
        $c_0=L(1-k,psi)/2$ if $L=1$ (that second $L$ is an $L$-function $L$).
        """
        c0, chi, psi, K, n, t, L, M = self.__defining_parameters()
        zero = K(0)
        k = self.weight()
        v = [c0]
        for i in X:
            if i==0: continue
            if i%t != 0:
                v.append(zero)
            else:
                m = i//t
                v.append(sum([psi(n)*chi(m/n)*n**(k-1) for \
                               n in rings.divisors(m)]))
        return v

    def __defining_parameters(self):
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
        return self.__chi

    def psi(self):
        return self.__psi

    def t(self):
        return self.__t

    def parameters(self):
        """
        Return chi, psi, and t, which are the defining parameters of self.
        """
        return self.__chi, self.__psi, self.__t

    def L(self):
        return self.__chi.conductor()

    def M(self):
        return self.__psi.conductor()

    def character(self):
        try:
            return self.__character
        except AttributeError:
            self.__character = self.__chi * (~self.__psi)
        return self.__character

    def new_level(self):
        if self.__chi.is_trivial() and self.__psi.is_trivial() and self.weight() == 2:
            return rings.factor(self.__t)[0][0]
        return self.L()*self.M()



