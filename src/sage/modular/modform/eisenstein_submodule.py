"""
The Eisenstein Subspace

EXAMPLES:

"""

from sage.structure.all import Sequence
from sage.misc.all import verbose
import sage.rings.all as rings
from sage.categories.all import Objects

import eis_series
import element
import submodule

class EisensteinSubmodule(submodule.ModularFormsSubmodule):
    """
    The Eisenstein submodule of an ambient space of modular forms.
    """
    def __init__(self, ambient_space):
        verbose('creating eisenstein submodule of %s'%ambient_space)
        d = ambient_space._dim_eisenstein()
        V = ambient_space.module()
        n = V.dimension()
        self._start_position = int(n - d)
        S = V.submodule([V.gen(i) for i in range(n-d,n)], check=False,
                        already_echelonized=True)
        submodule.ModularFormsSubmodule.__init__(self, ambient_space, S)

    def _repr_(self):
        return "Eisenstein subspace of dimension %s of %s"%(self.dimension(), self.ambient_module())

    def eisenstein_submodule(self):
        return self

    def modular_symbols(self, sign=0):
        r"""
        Return the corresponding space of modular symbols with given sign.

        WARNING: If sign != 0, then the space of modular symbols will,
        in general, only correspond to a \emph{subspace} of this space
        of modular forms.  This can be the case for both sign +1 or -1.

        EXAMPLES:
            sage: E = ModularForms(11,2).eisenstein_submodule()
            sage: M = E.modular_symbols(); M
            Modular Symbols subspace of dimension 1 of Modular Symbols space
            of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: M.sign()
            0

            sage: M = E.modular_symbols(sign=-1); M
            Modular Symbols subspace of dimension 0 of Modular Symbols space of
            dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

            sage: E = ModularForms(1,12).eisenstein_submodule()
            sage: E.modular_symbols()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of
            dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field

            sage: eps = DirichletGroup(13).0
            sage: E = EisensteinForms(eps^2, 2)
            sage: E.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
        """
        try:
            return self.__modular_symbols[sign]
        except AttributeError:
            self.__modular_symbols = {}
        except KeyError:
            pass
        A = self.ambient_module()
        S = A.modular_symbols(sign).eisenstein_submodule()
        self.__modular_symbols[sign] = S
        return S

class EisensteinSubmodule_params(EisensteinSubmodule):
    def parameters(self):
        try:
            return self.__parameters
        except AttributeError:
            P = eis_series.compute_eisenstein_params(self._parameters_character(), self.weight())
            self.__parameters = P
            return P

    def _parameters_character(self):
        return self.character()

    def change_ring(self, base_ring):
        if base_ring == self.base_ring():
            return self
        A = self.ambient_module()
        B = A.change_ring(base_ring)
        return B.eisenstein_submodule()

    def eisenstein_series(self):
        """
        Return the Eisenstein series that span this space (over the
        algebraic closure).

        EXAMPLES:
            sage: EisensteinForms(11,2).eisenstein_series()
            [
            5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
            ]
            sage: EisensteinForms(1,4).eisenstein_series()
            [
            1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + 126*q^5 + O(q^6)
            ]
            sage: EisensteinForms(1,24).eisenstein_series()
            [
            236364091/131040 + q + 8388609*q^2 + 94143178828*q^3 + 70368752566273*q^4 + 11920928955078126*q^5 + O(q^6)
            ]
            sage: EisensteinForms(5,4).eisenstein_series()
            [
            1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + 126*q^5 + O(q^6),
            1/240 + q^5 + O(q^6)
            ]
            sage: EisensteinForms(13,2).eisenstein_series()
            [
            1/2 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
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

            sage: eps = DirichletGroup(13).0^2
            sage: ModularForms(eps,2).eisenstein_series()
            [
            -7/13*zeta6 - 11/13 + q + (2*zeta6 + 1)*q^2 + (-3*zeta6 + 1)*q^3 + (6*zeta6 - 3)*q^4 - 4*q^5 + O(q^6),
            q + (zeta6 + 2)*q^2 + (-zeta6 + 3)*q^3 + (3*zeta6 + 3)*q^4 + 4*q^5 + O(q^6)
            ]

	    sage: M = ModularForms(19,3).eisenstein_subspace()
	    sage: M.eisenstein_series()
            [
            ]
        """
        try:
            return self.__eisenstein_series
        except AttributeError:
            P = self.parameters()
            E = Sequence([element.EisensteinSeries(self.change_ring(chi.base_ring()),
                                                      None, t, chi, psi) for \
                                              chi, psi, t in P], immutable=True,
                         cr = True, universe=Objects())
            assert len(E) == self.dimension(), "bug in enumeration of Eisenstein series."
            self.__eisenstein_series = E
            return E

    def xxx_compute_q_expansion_basis(self, prec):
        E = self.eisenstein_series()
        QQ = rings.QQ
        V = QQ**prec
        G = []
        z = QQ(0)
        # TODO: Possibly optimization: include one form from each Galois orbit.
        for e in E:
            f = e.q_expansion(prec)
            w = f.padded_list(prec)
            if f.base_ring() == QQ:
                G.append(V(w))
            else:
                assert False, 'the following is just wrong in general'
                for i in range(f.base_ring().degree()):
                    G.append(V([x[i] for x in w]))
        #print Sequence(G,cr=True)
        W = V.submodule(G, check=False)
        R = self._q_expansion_ring()
        X = [R(f.list(), prec) for f in W.basis()]
        return X + [R(0,prec)]*(self.dimension() - len(X))

    def _compute_q_expansion_basis(self, prec):
        E = self.eisenstein_series()
        K = self.base_ring()
        V = K**prec
        G = []
        for e in E:
            f = e.q_expansion(prec)
            w = f.padded_list(prec)
            L = f.base_ring()
            if K.has_coerce_map_from(L):
                G.append(V(w))
            else:
                # restrict scalars from L to K
                r,d = cyclotomic_restriction(L,K)
                s = [r(x) for x in w]
                for i in range(d):
                    G.append(V([x[i] for x in s]))
            # end if
        # end for
        W = V.submodule(G, check=False)
        R = self._q_expansion_ring()
        X = [R(f.list(), prec) for f in W.basis()]
        return X + [R(0,prec)]*(self.dimension() - len(X))

    def _q_expansion(self, element, prec):
        B = self.q_expansion_basis(prec)
        f = self._q_expansion_zero()
        for i in range(self._start_position, len(element)):
            if element[i] != 0:
                f += element[i] * B[i - self._start_position]
        return f


class EisensteinSubmodule_g0_Q(EisensteinSubmodule_params):
    """
    Space of Eisenstein forms for Gamma0(N).
    """

class EisensteinSubmodule_g1_Q(EisensteinSubmodule_params):
    """
    Space of Eisenstein forms for Gamma1(N).
    """
    def _parameters_character(self):
        return self.level()


class EisensteinSubmodule_eps(EisensteinSubmodule_params):
    """
    Space of Eisenstein forms with given Dirichlet character.

    EXAMPLES:
	sage: e = DirichletGroup(27,CyclotomicField(3)).0
	sage: M = ModularForms(e,2,prec=10).eisenstein_subspace()
	sage: M.dimension()
	6

        sage: M.eisenstein_series()
        [
        -1/3*zeta6 - 1/3 + q + (2*zeta6 - 1)*q^2 + q^3 + (-2*zeta6 - 1)*q^4 + (-5*zeta6 + 1)*q^5 + O(q^6),
        -1/3*zeta6 - 1/3 + q^3 + O(q^6),
        q + (-2*zeta6 + 1)*q^2 + (-2*zeta6 - 1)*q^4 + (5*zeta6 - 1)*q^5 + O(q^6),
        q + (zeta6 + 1)*q^2 + 3*q^3 + (zeta6 + 2)*q^4 + (-zeta6 + 5)*q^5 + O(q^6),
        q^3 + O(q^6),
        q + (-zeta6 - 1)*q^2 + (zeta6 + 2)*q^4 + (zeta6 - 5)*q^5 + O(q^6)
        ]
        sage: M.eisenstein_subspace().T(2).matrix().fcp()
        (x + zeta3 + 2) * (x + 2*zeta3 + 1) * (x - 2*zeta3 - 1)^2 * (x - zeta3 - 2)^2
        sage: ModularSymbols(e,2).eisenstein_subspace().T(2).matrix().fcp()
        (x + zeta3 + 2) * (x + 2*zeta3 + 1) * (x - 2*zeta3 - 1)^2 * (x - zeta3 - 2)^2


	sage: M.basis()
        [
        1 - 3*zeta3*q^6 + (-2*zeta3 + 2)*q^9 + O(q^10),
        q + (5*zeta3 + 5)*q^7 + O(q^10),
        q^2 - 2*zeta3*q^8 + O(q^10),
        q^3 + (zeta3 + 2)*q^6 + 3*q^9 + O(q^10),
        q^4 - 2*zeta3*q^7 + O(q^10),
        q^5 + (zeta3 + 1)*q^8 + O(q^10)
        ]


    """
    # TODO
    #def _compute_q_expansion_basis(self, prec):
	#B = EisensteinSubmodule_params._compute_q_expansion_basis(self, prec)
        #raise NotImplementedError, "must restrict scalars down correctly."


from sage.rings.all import CyclotomicField, lcm, euler_phi

def cyclotomic_restriction(L,K):
    if not L.has_coerce_map_from(K):
        M = CyclotomicField(lcm(L.zeta_order(), K.zeta_order()))
        f = cyclotomic_restriction_tower(M,K)
        def g(x):
            return f(M(x)), euler_phi(M.zeta_order())//euler_phi(K.zeta_order())
    else:
        return cyclotomic_restriction_tower(L,K), \
               euler_phi(L.zeta_order())//euler_phi(K.zeta_order())


def cyclotomic_restriction_tower(L,K):
    """
    Suppose L/K is an extension of cyclotomic fields and L=Q(zeta_m).
    This function computes a map with the following property:

       INPUT: element alpha in L
       OUTPUT: a list that gives the coefficients of a polynomial f(x) in K[x]
               such that f(zeta_m) = alpha.
    """
    if not L.has_coerce_map_from(K):
        raise ValueError, "K must be contained in L"
    f = L.defining_polynomial()
    R = K['x']
    x = R.gen()
    g = R(f)
    h = g.factor()[0][0]    # all factors are the same, since everybody is Galois
    def z(a):
        return R(a.polynomial()) % h
    return z

