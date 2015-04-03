"""
Integrable Representations of Affine Lie Algebras
"""
#*****************************************************************************
#  Copyright (C) 2014 Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from cartan_type import CartanType
from sage.rings.all import ZZ, QQ
from sage.misc.all import cached_method
from root_space import RootSpace
from weight_space import WeightSpace

class IntegrableRepresentation():
    """
    This is a class for an irreducible representation of affine Lie algebras.

    INPUT:

        - ``Lam`` - a dominant weight in an extended weight lattice of affine type.

    OPTIONAL:

        - ``depth`` - a parameter indicating how far to push computations.

    REFERENCES:

    .. [Kac] Kac, *Infinite-dimensional Lie algebras*, Third Edition.
       Cambridge, 1990.

    .. [KMPS] Kass, Moody, Patera and Slansky, Affine Lie algebras,
       weight multiplicities, and branching rules. Vols. 1, 2. University of
       California Press, Berkeley, CA, 1990.

    .. [KacPeterson] Kac and Peterson. Infinite-dimensional Lie algebras, theta
       functions and modular forms. Adv. in Math. 53 (1984), no. 2, 125-264.

    If `\Lambda` is a dominant integral weight for an affine root system, there exists a unique
    integrable representation of highest weight `\Lambda`. If `\mu` is another weight such that
    `\Lambda-\mu` is in the root lattice, then multiplicity of `\mu` in this representation will
    be denoted `m(\mu)`.

    Let `\delta` be the nullroot. Then for fixed `\mu` the function `m(\mu-k\delta)` is
    a monotone increasing function of `\mu`. It is useful to take `\mu` to be such that this
    function is nonzero if and only if `k\ge 0`. Therefore we make the following definition.
    If `\mu` is such that `m(\mu)\\ne 0` but `m(\mu+\delta)=0` then `\mu` is called *maximal*.

    Since `\delta` is fixed under the action of the affine Weyl group,
    and since the weight multiplicities are Weyl group invariant, the *string
    function* `k \mapsto m(\mu-k\delta)` is unchanged if `\mu` is replaced by an equivalent
    weight. Therefore in tabulating the string functions, we may assume that `\mu`
    is dominant. There are only a finite number of dominant maximal weights.

    Since every nonzero weight multiplicity appears in the string `\mu-k\delta` for
    one of the finite number of dominant maximal weights `\mu`, it is important to
    be able to compute these. We may do this as follows.

    EXAMPLE::

         sage: Lambda = RootSystem(['A',3,1]).weight_lattice(extended=true).fundamental_weights()
         sage: IntegrableRepresentation(Lambda[1]+Lambda[2]+Lambda[3]).strings()
         3*Lambda[2] - delta: 3 21 107 450 1638 5367 16194 45687 121876 310056 757056 1783324
         2*Lambda[0] + Lambda[2]: 4 31 161 665 2380 7658 22721 63120 166085 417295 1007601 2349655
         Lambda[1] + Lambda[2] + Lambda[3]: 1 10 60 274 1056 3601 11199 32354 88009 227555 563390 1343178
         Lambda[0] + 2*Lambda[3]: 2 18 99 430 1593 5274 16005 45324 121200 308829 754884 1779570
         Lambda[0] + 2*Lambda[1]: 2 18 99 430 1593 5274 16005 45324 121200 308829 754884 1779570

    In this example, we construct the extended weight lattice of Cartan type ['A',3,1],
    then define ``Lambda`` to be the fundamental weights. We find there are 5 maximal
    dominant weights in irreducible representation of highest weight ``Lambda[1]+Lambda[2]+Lambda[3]``,
    and we determine their string functions. If you want more values, give IntegrableRepresentation
    the optional parameter ``depth``, which defaults to 12.

    It was shown by Kac and Peterson that each string function is the set of
    Fourier coefficients of a modular form.

    Every weight `\mu` such that the weight multiplicity `m(\mu)` is nonzero has the
    form

      .. MATH::
          \Lambda - n_0\\alpha_0 - n_1\\alpha_1 - \cdots,

    where the `n_i` are nonnegative integers. This is represented internally as a
    tuple ``(n0, n1, n2, ...)``. If you want an individual multiplicity you use
    the method ``m`` and supply it with this tuple. The ``IntegrableRepresentation`` 
    class has methods ``to_weight`` and ``from_weight`` to convert between this
    internal representation and the weight lattice::

        sage: v = IntegrableRepresentation(2*Lambda[0]); v
        The IntegrableRepresentation of ['A', 3, 1] with highest weight 2*Lambda[0]
        sage: v.m((5,4,3,2))
        11
        sage: delta = v._delta
        sage: v.from_weight(-2*Lambda[0] + 4*Lambda[3] - 5*delta)
        (5, 4, 3, 2)

    To get more values, use the depth parameter::

        sage: sage: L0 = RootSystem(["A",1,1]).weight_lattice(extended=true).fundamental_weight(0); L0
        Lambda[0]
        sage: IntegrableRepresentation(4*L0, depth=20).strings()
        4*Lambda[1] - 2*delta: 1 2 6 11 23 41 75 126 215 347 561 878 1368 2082 3153 4690 6936 10121 14677 21055
        2*Lambda[0] + 2*Lambda[1] - delta: 1 2 5 10 20 36 66 112 190 310 501 788 1230 1880 2850 4256 6303 9222 13396 19262
        4*Lambda[0]: 1 1 3 6 13 23 44 75 131 215 354 561 889 1368 2097 3153 4712 6936 10151 14677

    .. WARNING:

        Currently this code works only for untwisted Type A.

    """
    def __init__(self, Lam, depth=12, debug=False, new=True):
        self._depth = depth
        self._P = Lam.parent()
        self._RS = self._P.root_system
        self._Q = self._RS.root_lattice()
        self._Lam = self._P(Lam)
        self._cartan_matrix = self._RS.cartan_matrix()
        self._cartan_type = self._RS.cartan_type()
        self._classical_rank = self._cartan_type.classical().rank()
        self._Lambda = self._P.fundamental_weights()
        self._W = self._P.weyl_group(prefix="s")
        self._s = self._W.simple_reflections()
        self._index_set = self._P.index_set()
        self._alpha = {}
        for i in self._index_set:
            self._alpha[i]=self._P(self._RS.root_lattice().simple_root(i))
        self._rho = self._P.rho()
        self._delta = self._P.null_root()
        self._shift = {}
        self._smat = {}
        for i in self._index_set:
            for j in self._index_set:
                self._smat[(i,j)] = -self._cartan_matrix[i][j]
            self._smat[(i,i)] += 1
        for i in self._index_set:
            self._shift[i] = tuple([ZZ((self._Lam-self._Lam.weyl_action(self._s[i])).symmetric_form(self._Lambda[j])) for j in self._index_set])
        self._ddict = {}
        self._mdict = {tuple(0 for i in self._index_set):1}
        if debug:
            self._debug = 1
        else:
            self._debug = None
        self._new = new
        self._den0 = (self._Lam+self._rho).symmetric_form(self._Lam+self._rho)
        def from_classical_root(h):
            """
            Coerces a classical root into P.
            """
            index_set_classical = [i for i in self._index_set if i != 0]
            return sum(self._alpha[i]*h.monomial_coefficients().get(i,0) for i in index_set_classical)
        self._classical_roots = [from_classical_root(al) for al in self._Q.classical().roots()]
        self._classical_positive_roots = [from_classical_root(al) for al in self._Q.classical().positive_roots()]
        self.string(Lam)

    def __repr__(self):
        return "The IntegrableRepresentation of %s with highest weight %s"%(self._cartan_type, self._Lam)

    def to_weight(self, n):
        """
        returns `Lam - \sum n[i]*alpha[i]`
        """
        return self._Lam - sum(ZZ(n[i])*self._alpha[i] for i in self._index_set)

    def from_weight(self, mu):
        """
        returns `(n[0], n[1], ...)` such that `mu = Lam - \sum n[i]*alpha[i]`
        """
        return tuple([ZZ((self._Lam-mu).symmetric_form(self._Lambda[i])) for i in self._index_set])

    def s(self, n, i):
        """
        Implements the `i`-th simple reflection in the internal representation of
        weights by tuples.
        """
        ret = [n[j] for j in self._index_set]
        for j in self._index_set:
            ret[j] += self._shift[i][j]
        ret[i] -= sum(n[j]*self._cartan_matrix[i][j] for j in self._index_set)
        return tuple(ret)

    def to_dominant(self, n):
        """
        Returns the dominant weight equivalent to ``n`` in the internal
        representation of weights by tuples.
        """

        if self._ddict.has_key(n):
            return self._ddict[n]
        for i in self._index_set:
            m = self.s(n, i)
            if m[i] < n[i]:
                v = self.to_dominant(m)
                self._ddict[n] = v
                return v
        return n

    def _freudenthal_roots_imaginary(self, nu):
        """
        It is assumed that `\\nu` is in the root lattice `Q`. Returns the set of
        imaginary roots `\\alpha\in\Delta^+` such that `\\nu+\\alpha\in Q^+`.
        """
        kp = min(ZZ(nu.symmetric_form(self._Lambda[i])) for i in self._index_set)
        return [u*self._delta for u in range(1,kp+1)]

    def _freudenthal_roots_real(self, nu):
        """
        It is assumed that `nu` is in `Q`. Returns the set of realroots
        `\alpha\in\Delta^+` such that `nu+\alpha\in Q^+`.
        """
        ret = []
        for al in self._classical_positive_roots:
            if all((nu-al).symmetric_form(self._Lambda[i])>=0 for i in self._index_set):
                ret.append(al)
        for al in self._classical_roots:
            for ir in self._freudenthal_roots_imaginary(nu-al):
                ret.append(al+ir)
        return ret

    def _freudenthal_accum(self, nu, al):
        """
        Helper function for ``self.m_freudenthal``.
        """
        if self._debug:
            print "g: nu=%s, a=%s"%(nu, al)
        ret = 0
        k = 1
        while 1:
            mk = self.m(self.from_weight(nu+k*al))
            sf = al.symmetric_form(nu+k*al)
            if self._debug:
                print "h:   k=%s, m(nu+k*al)=%s, (a|nu+k*al)=%s"%(k, mk, sf)
            if mk == 0:
                break
            ret += 2*mk*sf
            k += 1
        return ret

    def _m_freudenthal(self, n):
        """
        Computes the weight multiplicity using the Freudenthal multiplicity formula.
        """
        if min(n)<0:
            return 0
        if self._mdict.has_key(n):
            return self._mdict[n]
        mu = self.to_weight(n)
        den = self._den0-(mu+self._rho).symmetric_form(mu+self._rho)
        num = 0
        if self._debug:
            print "f:", self._Lam-mu
        for al in self._freudenthal_roots_real(self._Lam-mu):
            num += self._freudenthal_accum(mu, al)
        for al in self._freudenthal_roots_imaginary(self._Lam-mu):
            num += self._classical_rank*self._freudenthal_accum(mu, al)
        if den == 0 or self._debug:
            print "i: num=%s, den=%s"%(num,den)
        if self._debug:
            print " e:%s,%s,%s"%(n.__repr__(), num, den)
        return num/den

    def m(self, n):
        """
        Returns the multiplicity
        """
        if self._mdict.has_key(n):
            return self._mdict[n]
        elif self._ddict.has_key(n):
            self._mdict[n] = self.m(self._ddict[n])
        m = self.to_dominant(n)
        ret = self._m_freudenthal(m)
        self._mdict[n] = ret
        return ret

    def dominant_maximal(self):
        """
        Returns the finite set of dominant maximal weights.
        """
        ret = set()
        for x in self._ddict.values():
            if self.m(x) > 0:
                if min(x) == 0:
                    ret.add(x)
                else:
                    y = tuple([u-1 for u in x])
                    if self.m(y) == 0:
                        ret.add(y)
        return [self.to_weight(x) for x in ret]

    def string(self, max):
        """
        INPUT:

            - ``max`` - a dominant maximal weight.

        Returns the list of multiplicities ``m(max-k*delta)`
        for ``k = 0,1,2,`` up to ``self._depth``.
            """
        
        ret = []
        k = self._depth
        for j in range(k):
            ret.append(self.m(self.from_weight(max-j*self._delta)))
        return ret

    def strings(self):
        """
        Returns the set of dominant maximal weights of self, together with the string
        coefficients for each.
        """
        for max in self.dominant_maximal():
            s = self.string(max)
            print "%s:"%max,
            for j in s:
                print j,
            print

