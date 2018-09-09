from sage.rings.ideal import Ideal_generic
from sage.misc.cachefunc import cached_method

from sage.structure.richcmp import op_EQ, op_NE, op_LT, op_GT, op_LE, op_GE

class TateAlgebraIdeal(Ideal_generic):

    @cached_method
    def groebner_basis(self, algorithm=None):
        if algorithm is None:
            algorithm = "buchberger"
        if algorithm == "buchberger":
            gb = self._groebner_basis_buchberger()
        elif algorithm == "F4":
            gb = self._groebner_basis_F4()
        else:
            raise ValueError("Algorithm must be 'buchberger', 'F4' or None")
        # We minimize the Groebner basis
        i = 0
        while i < len(gb):
            ti = gb[i].leading_term()
            for j in range(len(gb)):
                if j != i and ti.is_divisible_by(gb[j].leading_term()):
                    del gb[i]
                    break
            else:
                i += 1
        return gb

    def _groebner_basis_buchberger(self):
        from heapq import heappush, heappop
        gb = list(self.gens())
        S = [ ]; ind = 0
        for i in range(len(gb)):
            for j in range(i+1, len(gb)):
                if not gb[i].leading_term().is_coprime_with(gb[j].leading_term()):
                    s = gb[i].Spoly(gb[j])
                    heappush(S,(s.leading_term(), ind, s))
                    ind += 1
        while S:
            _, _, f = heappop(S)
            _, r = f.quo_rem(gb, integral=True)
            #print len(S), ":", r
            if not r.is_zero():
                for g in gb:
                    if not g.leading_term().is_coprime_with(r.leading_term()):
                        s = g.Spoly(r)
                        heappush(S,(s.leading_term(), ind, s))
                        ind += 1
                gb.append(r)
        return gb

    def _groebner_basis_F4(self):
        raise NotImplementedError

    def _contains_(self, x):
        gb = self.groebner_basis()
        return (x % gb).is_zero()

    def _contains_ideal(self, I):
        gb = self.groebner_basis()
        for f in I.gens():
            if not (f % gb).is_zero():
                return False
        return True

    def _richcmp_(self, other, op):
        if not isinstance(other, TateAlgebraIdeal):
            raise NotImplementedError
        if op == op_GT:
            return self._contains_ideal(other) and not other._contains_ideal(self)
        elif op == op_GE:
            return self._contains_ideal(other)
        elif op == op_EQ:
            return self._contains_ideal(other) and other._contains_ideal(self)
        elif op == op_NE:
            return not(self._contains_ideal(other) and other._contains_ideal(self))
        elif op == op_LE:
            return other._contains_ideal(self)
        elif op == op_LT:
            return other._contains_ideal(self) and not self._contains_ideal(other)

    def is_saturated(self):
        if self.ring().base_ring().is_field():
            return True
        gb = self.groebner_basis()
        for g in gb:
            if g.valuation() > 0:
                return False
        return True

    def saturate(self):
        gb = self.groebner_basis()
        gens = [ g.monic() for g in gb ]
        return self.ring().ideal(gens)
