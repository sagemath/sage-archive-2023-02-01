
import sage.modular.hecke.hecke_operator
from sage.rings.arith import is_prime
import heilbronn

class HeckeOperator(sage.modular.hecke.hecke_operator.HeckeOperator):
    def apply_sparse(self, x):
        if x not in self.domain():
            raise TypeError, "x (=%s) must be in %s"%(x, self.domain())

        # old version just to check for correctness
        #return self.hecke_module_morphism()(x)

        p = self.index()
        if is_prime(p):
            H = heilbronn.HeilbronnCremona(p)
        else:
            H = heilbronn.HeilbronnMerel(p)

        M = self.parent().module()
        mod2term = M._mod2term
        syms = M.manin_symbols()
        K = M.base_ring()
        R = M.manin_gens_to_basis()

        W = R.new_matrix(nrows=1, ncols = R.nrows())

        B = M.manin_basis()

        #from sage.all import cputime
        #t = cputime()
        v = x.element()
        for i in v.nonzero_positions():
            for h in H:
                entries = syms.apply(B[i], h)
                for k, w in entries:
                    f, s = mod2term[k]
                    if s:
                        W[0,f] += s*K(w)*v[i]

        #print 'sym', cputime(t)
        #t = cputime()
        ans = M( v.parent()((W * R).row(0)) )
        #print 'mul', cputime(t)
        #print 'density: ', len(W.nonzero_positions())/(W.nrows()*float(W.ncols()))

        return ans
