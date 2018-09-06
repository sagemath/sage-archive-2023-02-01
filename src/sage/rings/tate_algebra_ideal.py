from sage.rings.ideal import Ideal_generic
from sage.misc.cachefunc import cached_method

class TateAlgebraIdeal(Ideal_generic):

    @cached_method
    def groebner_basis(self, algorithm=None):
        if algorithm is None:
            algorithm = "buchberger"
        if algorithm == "buchberger":
            return self._groebner_basis_buchberger()
        elif algorithm == "F4":
            return self._groebner_basis_F4()
        else:
            raise ValueError("Algorithm must be 'buchberger', 'F4' or None")

    def _groebner_basis_buchberger(self):
        gb = list(self.gens())
        S = [ ]
        for i in range(len(gb)):
            for j in range(i+1, len(gb)):
                S.append(gb[i].Spoly(gb[j]))
        while S:
            f = S.pop()
            r = f % gb
            if not r.is_zero():
                for g in gb:
                    S.append(g.Spoly(r))
                gb.append(r)
        return gb

    def _groebner_basis_F4(self):
        raise NotImplementedError
