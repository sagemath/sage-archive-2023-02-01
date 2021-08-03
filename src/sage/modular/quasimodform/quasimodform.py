# r"""
# Graded quasi-modular forms ring

# .. NOTE:

#     Only the space of quasimodular forms for the full modular group have been implemented.

# AUTHORS:

# - DAVID AYOTTE (2021-03-18): initial version

# """

# # ****************************************************************************
# #       Copyright (C) 2021 DAVID AYOTTE
# #
# # This program is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 2 of the License, or
# # (at your option) any later version.
# #                  https://www.gnu.org/licenses/
# # ****************************************************************************


# from sage.modular.modform.eis_series import eisenstein_series_qexp
# from sage.modular.modform.constructor import ModularForms
# from sage.modular.arithgroup.all import Gamma0, is_CongruenceSubgroup
# from sage.modular.modform.element import GradedModularFormElement, ModularFormElement
# from sage.modular.modform.space import ModularFormsSpace

# from sage.rings.all import Integer, QQ, ZZ

# from sage.structure.parent import Parent
# from sage.structure.element import Element
# from sage.structure.unique_representation import UniqueRepresentation
# from sage.structure.richcmp import richcmp, op_NE, op_EQ

# from sage.categories.graded_modules import GradedModules
# from sage.categories.hecke_modules import HeckeModules
# from sage.categories.rings import Rings
# from sage.rings.find_generators import ModularFormsRing
# from sage.rings.polynomial.polynomial_element import Polynomial

# class QuasiModularFormsElement(Element):
#     r"""
#     A quasimodular forms ring element

#     .. TODO::

#     Move this class somewhere else?
#     """
#     def __init__(self, parent, polynomial):
#         r"""
#         An element of a graded ring of quasimodular form.

#         INPUTS:

#         - ``parent`` - QuasiModularForms
#         - ``polynomial`` - a polynomial `f_0 + f_1 E_2 + ... + f_n E_2^n` where each `f_i`
#           are modular forms or base ring elements and `E_2` correspond to the weight 2 Eisenstein series.

#         OUTPUT:

#         - ``QuasiModularFormsElement``

#         EXAMPLES::

#             sage: QM = QuasiModularForms()
#             sage: QM.gens()
#             [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
#             1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
#             1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
#             sage: QM.0 + QM.1
#             2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
#         """
#         if not isinstance(polynomial, Polynomial):
#             raise TypeError("'polynomial' argument should be of type 'Polynomial'")
#         for f in polynomial.coefficients():
#             if not isinstance(f, GradedModularFormElement):
#                 raise ValueError("at least one coefficient is not a 'GradedModularFormElement'")
#         self._polynomial = polynomial
#         self._coefficients = polynomial.coefficients(sparse=False)
#         self.__base_ring = parent.base_ring()
#         Element.__init__(self, parent)

#     def q_expansion(self, prec=6):
#         r"""
#         Computes the `q`-expansion of self to precision `prec`.

#         EXAMPLES:::

#             sage: QM = QuasiModularForms()
#             sage: E2 = QM.0
#             sage: E2.q_expansion()
#             1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
#             sage: E2.q_expansion(prec=10)
#             1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 - 288*q^6 - 192*q^7 - 360*q^8 - 312*q^9 + O(q^10)
#         """
#         E2 = eisenstein_series_qexp(2, prec=prec, K=self.__base_ring, normalization='constant') #normalization -> to force integer coefficients
#         return sum(f.q_expansion(prec=prec)*E2**idx for idx, f in enumerate(self._coefficients))

#     def _repr_(self):
#         r"""
#         String representation of self.

#         TESTS::

#             sage: QM = QuasiModularForms()
#             sage: QM.0
#             1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
#             sage: QM.1
#             1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
#             sage: QM.2
#             1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
#         """
#         return "%s" % (self.q_expansion())

#     def _richcmp_(self, other, op):
#         r"""
#         Compare self with other.

#         TESTS::

#             sage: QM = QuasiModularForms(1)
#             sage: QM.0 == QM.1
#             False
#             sage: QM.0 == QM.0
#             True
#             sage: QM.0 != QM.1
#             True
#             sage: QM.0 != QM.0
#             False
#             sage: QM.0 < QM.1
#             Traceback (most recent call last):
#             ...
#             TypeError: invalid comparison between quasimodular forms ring elements
#         """
#         if op != op_EQ and op != op_NE:
#             raise TypeError('invalid comparison between quasimodular forms ring elements')
#         return richcmp(self._polynomial, other._polynomial, op)

#     def _add_(self, other):
#         r"""
#         Addition of two ``QuasiModularFormElement``.

#         INPUT:

#         - ``other`` - ``QuasiModularFormElement``

#         OUTPUT: a ``QuasiModularFormElement``

#         TESTS::

#             sage: QM = QuasiModularForms(1)
#             sage: QM.0 + QM.1
#             2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
#             sage: QM.0 + (QM.1 + QM.2) == (QM.0 + QM.1) + QM.2
#             True
#         """
#         return self.__class__(self.parent(), self._polynomial + other._polynomial)

#     def __neg__(self):
#         r"""
#         The negation of ``self```

#         TESTS::

#             sage: -QuasiModularForms(1).0
#             -1 + 24*q + 72*q^2 + 96*q^3 + 168*q^4 + 144*q^5 + O(q^6)
#             sage: QuasiModularForms(1).0 - QuasiModularForms(1).0
#             0
#         """
#         return self.__class__(self.parent(), -self._polynomial)

#     def _mul_(self, other):
#         r"""
#         The multiplication of two ``QuasiModularFormElement``

#         INPUT:

#         - ``other`` - ``QuasiModularFormElement``

#         OUTPUT: a ``QuasiModularFormElement``

#         TESTS::

#             sage: QM = QuasiModularForms(1)
#             sage: QM.0 * QM.1
#             1 + 216*q - 3672*q^2 - 62496*q^3 - 322488*q^4 - 1121904*q^5 + O(q^6)
#             sage: (QM.0 * QM.1) * QM.2 == QM.0 * (QM.1 * QM.2)
#             True
#         """
#         return self.__class__(self.parent(), self._polynomial * other._polynomial)

#     def _lmul_(self, c):
#         r"""
#         The left action of the base ring on self.

#         INPUT:

#         - ``other`` - ``QuasiModularFormElement``

#         OUTPUT: a ``QuasiModularFormElement``

#         TESTS::

#             sage: QM = QuasiModularForms(1)
#             sage: (1/2) * QM.0
#             1/2 - 12*q - 36*q^2 - 48*q^3 - 84*q^4 - 72*q^5 + O(q^6)
#             sage: QM.0 * (3/2)
#             3/2 - 36*q - 108*q^2 - 144*q^3 - 252*q^4 - 216*q^5 + O(q^6)
#         """
#         return self.__class__(self.parent(), c * self._polynomial)


# class QuasiModularForms(Parent, UniqueRepresentation):
#     r"""
#     Base class of quasimodular forms ring.

#     ..TODO::

#     Move this class in sage.rings
#     """
#     Element = QuasiModularFormsElement
#     def __init__(self, group=1, base_ring=QQ, name='E2'):
#         r"""
#         The graded ring of quasimodular forms for the full modular group `{\rm SL}_2(\ZZ)`, with
#         coefficients in a ring.

#         INPUT:

#         - ``group`` (default: `{\rm SL}_2(\ZZ)`) -- a congruence subgroup of `{\rm SL}_2(\ZZ)`, or a
#           positive integer `N` (interpreted as `\Gamma_0(N)`)

#         - ``base_ring`` (ring, default: `\QQ`) -- a base ring, which should be
#           `\QQ`, `\ZZ`, or the integers mod `p` for some prime `p`.

#         EXAMPLES::

#             sage: QM = QuasiModularForms(); QM
#             Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field
#             sage: QM.gens()
#             [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
#             1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
#             1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]

#         .. TESTS:

#             sage: M = QuasiModularForms(1)
#             sage: M.group()
#             Modular Group SL(2,Z)
#             sage: M.base_ring()
#             Rational Field
#             sage: M = QuasiModularForms(1, Integers(5))
#             sage: M.base_ring()
#             Ring of integers modulo 5
#             sage: QuasiModularForms(2)
#             Traceback (most recent call last):
#             ...
#             NotImplementedError: space of quasimodular forms are only implemented for the full modular group
#             sage: QuasiModularForms(Integers(5))
#             Traceback (most recent call last):
#             ...
#             ValueError: Group (=Ring of integers modulo 5) should be a congruence subgroup
#             sage: M2 = QuasiModularForms(1, GF(7))
#             sage: M == M2
#             False
#         """
#         #check if the group is SL2(Z)
#         if isinstance(group, (int, Integer)):
#             if group>1:
#                 raise NotImplementedError("space of quasimodular forms are only implemented for the full modular group")
#             group = Gamma0(1)
#         elif not is_CongruenceSubgroup(group):
#             raise ValueError("Group (=%s) should be a congruence subgroup" % group)
#         elif group is not Gamma0(1):
#             raise NotImplementedError("space of quasimodular forms are implemented for the full modular group")

#         #Check if the base ring is a field
#         #For some reasons, there is a problem when computing a basis of ModularForms
#         if not base_ring.is_field():
#             raise ValueError("The base ring must be a field")

#         self.__group = group
#         self.__base_ring = base_ring
#         self.__modular_forms_subring = ModularFormsRing(group, base_ring)
#         self.__polynomial_subring = self.__modular_forms_subring[name]
#         Parent.__init__(self, base=base_ring, category=GradedModules(base_ring))

#     def group(self):
#         r"""
#         Return the congruence subgroup for which this is the ring of quasimodular forms.

#         EXAMPLES::

#             sage: QM = QuasiModularForms(1)
#             sage: QM.group() is SL2Z
#             True
#             sage: QM = QuasiModularForms(Gamma0(1)); QM
#             Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field

#         Higher level congruence subgroups are not yet implemented::

#             sage: QuasiModularForms(2)
#             Traceback (most recent call last):
#             ...
#             NotImplementedError: space of quasimodular forms are only implemented for the full modular group
#         """
#         return self.__group

#     def base_ring(self):
#         r"""
#         Return the coefficient ring of this quasimodular forms ring.

#         EXAMPLES::

#             sage: QuasiModularForms(1).base_ring()
#             Rational Field
#             sage: QuasiModularForms(1, base_ring=Integers(5)).base_ring()
#             Ring of integers modulo 5
#         """
#         return self.__base_ring

#     def modular_forms_subring(self):
#         r"""
#         Return the subring of modular forms of this ring of quasimodular forms.

#         EXAMPLES::

#             sage: QuasiModularForms(1).modular_forms_subring()
#             Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
#         """
#         return self.__modular_forms_subring

#     def _repr_(self):
#         r"""
#         String representation of self.

#         EXAMPLES::

#             sage: QuasiModularForms(1)._repr_()
#             'Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field'
#             sage: QuasiModularForms(1, base_ring=Integers(13))._repr_()
#             'Ring of Quasimodular Forms for Modular Group SL(2,Z) over Ring of integers modulo 13'
#         """
#         return "Ring of Quasimodular Forms for %s over %s" % (self.group(), self.base_ring())

#     def _coerce_map_from_(self, M):
#         r"""
#         Code to make QuasiModularForms work well with coercion framework.

#         TESTS::

#             sage: E2 = QuasiModularForms(1).0
#             sage: M = ModularFormsRing(1)
#             sage: E2 + M.0
#             2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
#             sage: M.0 + E2
#             2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
#             sage: 1 + E2
#             2 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
#             sage: E2 + 1
#             2 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
#             sage: f = ModularForms(1, 12).0
#             sage: E2 + f
#             1 - 23*q - 96*q^2 + 156*q^3 - 1640*q^4 + 4686*q^5 + O(q^6)
#         """
#         if isinstance(M, (ModularFormsRing, ModularFormsSpace)):
#             if M.group() == self.group() and self.has_coerce_map_from(M.base_ring()):
#                 return True
#         if self.base_ring().has_coerce_map_from(M):
#             return True
#         return False

#     def _element_constructor_(self, datum):
#         r"""
#         The call method of self.

#         INPUT:

#         - datum - list, GradedModularFormElement, ModularFormElement, Polynomial, base ring element

#         OUTPUT: QuasiModularFormElement
#         """
#         if isinstance(datum, list):
#             if len(datum) == 0:
#                 raise ValueError("the given list should be non-empty")
#             for idx, f in enumerate(datum):
#                 if not isinstance(f, (GradedModularFormElement, ModularFormElement)):
#                     raise ValueError("one list element is not a modular form")
#                 datum[idx] = self.__modular_forms_subring(f) #to ensure that every forms is a GradedModularFormElement
#             datum = self.__polynomial_subring(datum)
#         elif isinstance(datum, (GradedModularFormElement, ModularFormElement)):
#             datum = self.__modular_forms_subring(datum)
#             datum = self.__polynomial_subring(datum)
#         elif isinstance(datum, Polynomial):
#             datum = self.__polynomial_subring(datum.coefficients(sparse=False))
#         elif self.base_ring().has_coerce_map_from(datum.base_ring()):
#             datum = self.__polynomial_subring(datum)
#         return self.element_class(self, datum)

#     def weigt_2_eisenstein_series(self):
#         r"""
#         Return the weight 2 Eisenstein series.

#         EXAMPLES::

#             sage: QM = QuasiModularForms(1)
#             sage: E2 = QM.weigt_2_eisenstein_series(); E2
#             1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
#             sage: E2.parent()
#             Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field
#         """
#         return self(self.__polynomial_subring.gen())

#     def gens(self):
#         r"""
#         Return a list of generators of the quasimodular forms ring. Note that the generators
#         of the modular forms subring is given are the one given by the method
#         :meth: `sage.rings.find_generators.ModularFormsRing.gen_forms`

#         EXAMPLES::

#             sage: QM = QuasiModularForms(1)
#             sage: QM.gens()
#             [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
#             1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
#             1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
#             sage: QM.modular_forms_subring().gen_forms()
#             [1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
#             1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
#         """
#         gen_list = [self.weigt_2_eisenstein_series()]
#         for f in self.__modular_forms_subring.gen_forms():
#             gen_list.append(self(f))
#         return gen_list

#     generators = gens # alias

#     def gen(self, n):
#         r"""
#         Return the `n`-th generator of the quasimodular forms ring.

#         EXAMPLES::

#             sage: QM = QuasiModularForms(1)
#             sage: QM.0
#             1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
#             sage: QM.1
#             1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
#             sage: QM.2
#             1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
#         """
#         return self.gens()[n]

#     def polygen(self):
#         r"""
#         Return the generator of this quasimodular form space as a polynomial ring over the modular
#         form subring. Note that this generator correspond to the weight-2 Eisenstein series. The default
#         name of this generator is 'E2'.

#         EXAMPLES::

#             sage: QM = QuasiModularForms(1)
#             sage: QM.polygen()
#             E2
#             sage: QuasiModularForms(1, name='X').polygen()
#             X
#             sage: QM.polygen().parent()
#             Univariate Polynomial Ring in E2 over Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
#         """
#         return self.__polynomial_subring.gen()

#     def differentiation_operator(self, f):
#         r"""
#         Compute the formal derivative `q\frac{d}{dq}` of the q-expansion of a quasimodular form `f`

#         INPUT:

#         - ``f`` -- a power serie in corresponding to the q-expansion of a quasimodular form.

#         OUTPUT:

#         The power series `q\frac{d}{dq}(f)`

#         EXAMPLES::

#             sage: M = QuasiModularForms()
#             sage: D = M.differentiation_operator
#             sage: B = M.gens()
#             sage: P = B[0]; Q = B[1]; R = B[2]
#             sage: D(P)
#             -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 + O(q^6)
#             sage: (P.q_expansion()^2 - Q.q_expansion())/12
#             -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 + O(q^6)
#             sage: D(Q)
#             240*q + 4320*q^2 + 20160*q^3 + 70080*q^4 + 151200*q^5 + O(q^6)
#             sage: (P.q_expansion()*Q.q_expansion() - R.q_expansion())/3
#             240*q + 4320*q^2 + 20160*q^3 + 70080*q^4 + 151200*q^5 + O(q^6)
#             sage: D(R)
#             -504*q - 33264*q^2 - 368928*q^3 - 2130912*q^4 - 7877520*q^5 + O(q^6)
#             sage: (P.q_expansion()*R.q_expansion() - Q.q_expansion()^2)/2
#             -504*q - 33264*q^2 - 368928*q^3 - 2130912*q^4 - 7877520*q^5 + O(q^6)

#         TODO:: This method need some work. It should return a QuasiModularFormsElement (not a power series in q)
#         """
#         q = f.q_expansion().parent().gen()
#         return q*f.q_expansion().derivative()
