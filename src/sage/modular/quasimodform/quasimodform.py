r"""
Graded quasi-modular forms ring

TODO: add more info

.. NOTE:

    Currently, all the methods are implemented only for the full modular groups. Congruence subgroups
    of higher level are not yet supported

AUTHORS:

- DAVID AYOTTE (2021-03-18): initial version

"""

# ****************************************************************************
#       Copyright (C) 2021 DAVID AYOTTE
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.modular.modform.eis_series import eisenstein_series_qexp
from sage.modular.arithgroup.all import Gamma0, is_CongruenceSubgroup
from sage.rings.all import Integer, QQ, ZZ
from sage.structure.sage_object  import SageObject

class QuasiModularFormsRing(SageObject):
    def __init__(self, group=1, base_ring=QQ):
        r"""
        The graded ring of quasimodular forms for the full modular group `{\rm SL}_2(\ZZ)`, with
        coefficients in a ring.

        INPUT:

        - ``group`` (default: `{\rm SL}_2(\ZZ)`) -- a congruence subgroup of `{\rm SL}_2(\ZZ)`, or a
          positive integer `N` (interpreted as `\Gamma_0(N)`)

        - ``base_ring`` (ring, default: `\QQ`) -- a base ring, which should be
          `\QQ`, `\ZZ`, or the integers mod `p` for some prime `p`.

        EXAMPLES::

            sage: M = QuasiModularFormsRing(); M
            Ring of quasimodular forms for Modular Group SL(2,Z) with coefficients in Rational Field
            sage: B = M.generators(); B
            [(2,
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 - 288*q^6 - 192*q^7 - 360*q^8 - 312*q^9 + O(q^10)),
            (4,
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + 82560*q^7 + 140400*q^8 + 181680*q^9 + O(q^10)),
            (6,
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 - 8471232*q^7 - 17047800*q^8 - 29883672*q^9 + O(q^10))]
            sage: P = B[0][1]; Q = B[1][1]
            sage: D = M.differentiation_operator
            sage: D(P)
            -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 - 1728*q^6 - 1344*q^7 - 2880*q^8 - 2808*q^9 + O(q^10)
            sage: (P^2 - Q)/12
            -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 - 1728*q^6 - 1344*q^7 - 2880*q^8 - 2808*q^9 + O(q^10)
            sage: M = QuasiModularFormsRing(1, Integers(5)); M
            Ring of quasimodular forms for Modular Group SL(2,Z) with coefficients in Ring of integers modulo 5
            sage: B = M.generators(); B
            [(2, 1 + q + 3*q^2 + 4*q^3 + 2*q^4 + q^5 + 2*q^6 + 3*q^7 + 3*q^9 + O(q^10)),
            (4, 1 + O(q^10)),
            (6, 1 + q + 3*q^2 + 4*q^3 + 2*q^4 + q^5 + 2*q^6 + 3*q^7 + 3*q^9 + O(q^10))]

        .. TESTS:

            sage: M = QuasiModularFormsRing(1)
            sage: M.group()
            Modular Group SL(2,Z)
            sage: M.base_ring()
            Rational Field
            sage: M = QuasiModularFormsRing(1, ZZ)
            sage: M.base_ring()
            Integer Ring
            sage: M = QuasiModularFormsRing(1, Integers(5))
            sage: M.base_ring()
            Ring of integers modulo 5
            sage: QuasiModularFormsRing(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: The space of quasimodular forms for higher levels are not yet implemented
            sage: QuasiModularFormsRing(Integers(5))
            Traceback (most recent call last):
            ...
            ValueError: Group (=Ring of integers modulo 5) should be a congruence subgroup           


        """
        if isinstance(group, (int, Integer)):
            if group>1:
                raise NotImplementedError("The space of quasimodular forms for higher levels are not yet implemented")
            group = Gamma0(1)
        elif not is_CongruenceSubgroup(group):
            raise ValueError("Group (=%s) should be a congruence subgroup" % group)
        elif group is not Gamma0(1):
            raise NotImplementedError("The space of quasimodular forms for higher levels are not yet implemented")
        
        self.__group = group
        self.__base_ring = base_ring

    def group(self):
        r"""
        Return the congruence subgroup for which this is the ring of quasimodular forms.

        EXAMPLES::

            sage: M = QuasiModularFormsRing(1)
            sage: M.group() is SL2Z
            True
            sage: M = QuasiModularFormsRing(Gamma0(1)); M
            Ring of quasimodular forms for Modular Group SL(2,Z) with coefficients in Rational Field
            
        Higher level congruence subgroups are not yet implemented::

            sage: QuasiModularFormsRing(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: The space of quasimodular forms for higher levels are not yet implemented

        """
        return self.__group

    def base_ring(self):
        r"""
        Return the coefficient ring of this quasimodular forms ring.

        EXAMPLES::

            sage: QuasiModularFormsRing(1).base_ring()
            Rational Field
            sage: QuasiModularFormsRing(1, base_ring = ZZ).base_ring()
            Integer Ring
            sage: QuasiModularFormsRing(1, base_ring = Integers(5)).base_ring()
            Ring of integers modulo 5
        """
        return self.__base_ring

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: QuasiModularFormsRing(1)._repr_()
            'Ring of quasimodular forms for Modular Group SL(2,Z) with coefficients in Rational Field'
            sage: QuasiModularFormsRing(1, base_ring=Integers(13))._repr_()
            'Ring of quasimodular forms for Modular Group SL(2,Z) with coefficients in Ring of integers modulo 13'
        """
        return "Ring of quasimodular forms for %s with coefficients in %s" % (self.group(), self.base_ring())

    def generators(self, prec=10):
        r"""
        If `R` is the base ring of self, then this method returns a set of
        quasimodular forms which generate the `R`-algebra of all quasimodular forms.

        INPUT:

        - ``prec`` (integer, default: 10) -- return `q`-expansions to this
          precision

        OUPUT:

        a list of pairs (k, f), where f is the q-expansion to precision
        ``prec`` of a quasimodular form of weight k. For the full modular group, these
        forms are precisely the normalized eisenstein series of weight 2, 4 and 6 respectively.

        EXAMPLES::

            sage: M = QuasiModularFormsRing(1)
            sage: M.generators()
            [(2,
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 - 288*q^6 - 192*q^7 - 360*q^8 - 312*q^9 + O(q^10)),
            (4,
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + 82560*q^7 + 140400*q^8 + 181680*q^9 + O(q^10)),
            (6,
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 - 8471232*q^7 - 17047800*q^8 - 29883672*q^9 + O(q^10))]
            sage: QuasiModularFormsRing(1, Integers(17)).generators(prec=6)
            [(2, 1 + 10*q + 13*q^2 + 6*q^3 + 2*q^4 + 9*q^5 + O(q^6)),
            (4, 1 + 2*q + q^2 + 5*q^3 + 10*q^4 + 14*q^5 + O(q^6)),
            (6, 1 + 6*q + 11*q^2 + 2*q^3 + q^4 + 5*q^5 + O(q^6))]
        """
        E2 = eisenstein_series_qexp(2, prec=prec, K=self.base_ring(), normalization='constant')
        E4 = eisenstein_series_qexp(4, prec=prec, K=self.base_ring(), normalization='constant')
        E6 = eisenstein_series_qexp(6, prec=prec, K=self.base_ring(), normalization='constant')
        return [(2, E2), (4, E4), (6, E6)]

    def differentiation_operator(self, f):
        r"""
        Compute the formal derivative `q\frac{d}{dq}` of the q-expansion of a quasimodular form `f`

        INPUT:

        - ``f`` -- a power serie in corresponding to the q-expansion of a quasimodular form.

        OUTPUT:

        The power serie `q\frac{d}{dq}(f)`

        EXAMPLES::

            sage: M = QuasiModularFormsRing()
            sage: D = M.differentiation_operator
            sage: B = M.generators()
            sage: P = B[0][1]; Q = B[1][1]; R = B[2][1]
            sage: D(P)
            -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 - 1728*q^6 - 1344*q^7 - 2880*q^8 - 2808*q^9 + O(q^10)
            sage: (P^2 - Q)/12
            -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 - 1728*q^6 - 1344*q^7 - 2880*q^8 - 2808*q^9 + O(q^10)
            sage: D(Q)
            240*q + 4320*q^2 + 20160*q^3 + 70080*q^4 + 151200*q^5 + 362880*q^6 + 577920*q^7 + 1123200*q^8 + 1635120*q^9 + O(q^10)
            sage: (P*Q - R)/3
            240*q + 4320*q^2 + 20160*q^3 + 70080*q^4 + 151200*q^5 + 362880*q^6 + 577920*q^7 + 1123200*q^8 + 1635120*q^9 + O(q^10)
            sage: D(R)
            -504*q - 33264*q^2 - 368928*q^3 - 2130912*q^4 - 7877520*q^5 - 24349248*q^6 - 59298624*q^7 - 136382400*q^8 - 268953048*q^9 + O(q^10)
            sage: (P*R - Q^2)/2
            -504*q - 33264*q^2 - 368928*q^3 - 2130912*q^4 - 7877520*q^5 - 24349248*q^6 - 59298624*q^7 - 136382400*q^8 - 268953048*q^9 + O(q^10)
        """
        q = f.parent().gen()
        return q*f.derivative()
    

    
