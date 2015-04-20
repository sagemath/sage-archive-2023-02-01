"""
Integrable Representations of Affine Lie Algebras
"""

#*****************************************************************************
#  Copyright (C) 2014, 2105 Daniel Bump <bump at match.stanford.edu>
#                           Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.category_object import CategoryObject
from sage.categories.modules import Modules
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_space import RootSpace
from sage.combinat.root_system.weight_space import WeightSpace
from sage.rings.all import ZZ, QQ
from sage.misc.all import cached_method
from sage.matrix.constructor import Matrix
from sage.functions.other import floor
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet

# TODO: Make this a proper parent and implement actions
class IntegrableRepresentation(CategoryObject, UniqueRepresentation):
    r"""
    An irreducible highest weight representation of an affine Lie algebra.

    INPUT:

    - ``Lam`` -- a dominant weight in an extended weight lattice
      of affine type

    REFERENCES:

    .. [Kac] Kac, *Infinite-dimensional Lie algebras*, Third Edition.
       Cambridge, 1990.

    .. [KMPS] Kass, Moody, Patera and Slansky, Affine Lie algebras,
       weight multiplicities, and branching rules. Vols. 1, 2. University of
       California Press, Berkeley, CA, 1990.

    .. [KacPeterson] Kac and Peterson. Infinite-dimensional Lie algebras, theta
       functions and modular forms. Adv. in Math. 53 (1984), no. 2, 125-264.

    If `\Lambda` is a dominant integral weight for an affine root system,
    there exists a unique integrable representation `V=V_\Lambda` of highest weight
    `\Lambda`. If `\mu` is another weight, let `m(\mu)` denote the
    multiplicity of the weight `\mu` in this representation. The set
    `\text{supp}(V)` of `\mu` such that `m(\mu)>0` is contained in the
    paraboloid

    .. MATH::
             (\Lambda+\rho|\Lambda+\rho) - (\mu+\rho|\mu+\rho) \ge 0

    where `(\,|\,)` is the invariant inner product on the weight
    lattice and `\rho` is the Weyl vector. Moreover every `\mu` differs from
    `\Lambda` by an element of the root lattice. ([Kac], Propositions 11.3 and
    11.4) such that `\Lambda - \mu` is in the root lattice, then multiplicity
    of `\mu` in this representation will be denoted `m(\mu)`. The set of `\mu`
    such that the multiplicity
    
    Let `\delta` be the nullroot, which is the lowest positive imaginary
    root. Then by [Kac], Proposition 11.3 or Corollary 11.9, for fixed `\mu`
    the function `m(\mu - k\delta)` is a monotone increasing function of
    `k`. It is useful to take `\mu` to be such that this function is nonzero
    if and only if `k \geq 0`. Therefore we make the following definition.  If
    `\mu` is such that `m(\mu) \neq 0` but `m(\mu + \delta) = 0` then `\mu` is
    called *maximal*.

    Since `\delta` is fixed under the action of the affine Weyl group,
    and since the weight multiplicities are Weyl group invariant, the
    function `k \mapsto m(\mu-k\delta)` is unchanged if `\mu` is replaced by
    an equivalent weight. Therefore in tabulating these functions, we may
    assume that `\mu` is dominant. There are only a finite number of dominant
    maximal weights.

    Since every nonzero weight multiplicity appears in the string
    `\mu - k\delta` for one of the finite number of dominant maximal
    weights `\mu`, it is important to be able to compute these. We may
    do this as follows.

    EXAMPLE::

         sage: Lambda = RootSystem(['A',3,1]).weight_lattice(extended=true).fundamental_weights()
         sage: IntegrableRepresentation(Lambda[1]+Lambda[2]+Lambda[3]).strings()
         2*Lambda[0] + Lambda[2]: 4 31 161 665 2380 7658 22721 63120 166085 417295 1007601 2349655
         Lambda[0] + 2*Lambda[1]: 2 18 99 430 1593 5274 16005 45324 121200 308829 754884 1779570
         Lambda[0] + 2*Lambda[3]: 2 18 99 430 1593 5274 16005 45324 121200 308829 754884 1779570
         Lambda[1] + Lambda[2] + Lambda[3]: 1 10 60 274 1056 3601 11199 32354 88009 227555 563390 1343178
         3*Lambda[2] - delta: 3 21 107 450 1638 5367 16194 45687 121876 310056 757056 1783324
         sage: Lambda = RootSystem(['D',4,1]).weight_lattice(extended=true).fundamental_weights()
         sage: IntegrableRepresentation(Lambda[0]+Lambda[1]).strings()                        # long time
         Lambda[0] + Lambda[1]: 1 10 62 293 1165 4097 13120 38997 109036 289575 735870 1799620
         Lambda[3] + Lambda[4] - delta: 3 25 136 590 2205 7391 22780 65613 178660 463842 1155717 2777795

    In this example, we construct the extended weight lattice of Cartan
    type `A_3^{(1)}`, then define ``Lambda`` to be the fundamental
    weights. We find there are 5 maximal dominant weights in irreducible
    representation of highest weight ``Lambda[1]+Lambda[2]+Lambda[3]``,
    and we determine their string functions.

    It was shown in [KacPeterson] that each string is the set of Fourier
    coefficients of a modular form.

    Every weight `\mu` such that the weight multiplicity `m(\mu)` is
    nonzero has the form

      .. MATH::

          \Lambda - n_0 \alpha_0 - n_1 \alpha_1 - \cdots,

    where the `n_i` are nonnegative integers. This is represented internally
    as a tuple ``(n0, n1, n2, ...)``. If you want an individual multiplicity
    you use the method ``m`` and supply it with this tuple::

        sage: Lambda = RootSystem(['C',2,1]).weight_lattice(extended=true).fundamental_weights()
        sage: v = IntegrableRepresentation(2*Lambda[0]); v
        The Integrable representation of ['C', 2, 1] with highest weight 2*Lambda[0]
        sage: v.m((3,5,3))
        18

    The :class:`IntegrableRepresentation` class has methods :meth:`to_weight`
    and :meth:`from_weight` to convert between this internal representation
    and the weight lattice::

        sage: delta = v.null_root()
        sage: v.to_weight((4,3,2))
        -3*Lambda[0] + 6*Lambda[1] - Lambda[2] - 4*delta
        sage: v.from_weight(-3*Lambda[0] + 6*Lambda[1] - Lambda[2] - 4*delta)
        (4, 3, 2)

    To get more values, use the depth parameter::

        sage: L0 = RootSystem(["A",1,1]).weight_lattice(extended=true).fundamental_weight(0); L0
        Lambda[0]
        sage: IntegrableRepresentation(4*L0).strings(depth=20)
        4*Lambda[0]: 1 1 3 6 13 23 44 75 131 215 354 561 889 1368 2097 3153 4712 6936 10151 14677
        2*Lambda[0] + 2*Lambda[1] - delta: 1 2 5 10 20 36 66 112 190 310 501 788 1230 1880 2850 4256 6303 9222 13396 19262
        4*Lambda[1] - 2*delta: 1 2 6 11 23 41 75 126 215 347 561 878 1368 2082 3153 4690 6936 10121 14677 21055

    An example in type `C_2^{(1)}`::

        sage: Lambda = RootSystem(['C',2,1]).weight_lattice(extended=true).fundamental_weights()
        sage: v = IntegrableRepresentation(2*Lambda[0])
        sage: v.strings()    # long time
        2*Lambda[0]: 1 2 9 26 77 194 477 1084 2387 5010 10227 20198
        Lambda[0] + Lambda[2] - delta: 1 5 18 55 149 372 872 1941 4141 8523 17005 33019
        2*Lambda[1] - delta: 1 4 15 44 122 304 721 1612 3469 7176 14414 28124
        2*Lambda[2] - 2*delta: 2 7 26 72 194 467 1084 2367 5010 10191 20198 38907
    """
    def __init__(self, Lam):
        """
        Initialize ``self``.
        """
        CategoryObject.__init__(self, base=ZZ, category=Modules(ZZ))

        if not Lam.parent()._extended:
            raise ValueError("The parent of %s must be an extended affine root lattice."%Lam)
        self._Lam = Lam
        self._P = Lam.parent()
        self._Q = self._P.root_system.root_lattice()

        self._cartan_matrix = self._P.root_system.cartan_matrix()
        self._cartan_type = self._P.root_system.cartan_type()
        if not self._cartan_type.is_untwisted_affine():
            raise NotImplementedError("Integrable Representations are only implemented for Untwisted Affine Types.")
        self._classical_rank = self._cartan_type.classical().rank()
        self._index_set = self._P.index_set()
        self._index_set_classical = self._cartan_type.classical().index_set()
        self._cminv = self._cartan_type.classical().cartan_matrix().inverse()

        self._ddict = {}
        self._mdict = {tuple(0 for i in self._index_set): 1}
        self._rho = self._P.rho()
        self._delta = self._Q.null_root()
        # Coerce a classical root into Q
        from_cl_root = lambda h: self._Q._from_dict(h._monomial_coefficients)
        self._classical_roots = [from_cl_root(al)
                                 for al in self._Q.classical().roots()]
        self._classical_positive_roots = [from_cl_root(al)
                                          for al in self._Q.classical().positive_roots()]
        self._eps = {i: self._cartan_type.a()[i] / self._cartan_type.dual().a()[i]
                     for i in self._index_set}
        self._a = self._cartan_type.a()
        self._ac = self._cartan_type.dual().a()
        self._coxeter_number = sum(self._a[i] for i in self._index_set)
        self._dual_coxeter_number = sum(self._ac[i] for i in self._index_set)
        E = Matrix.diagonal([self._eps[i] for i in self._index_set_classical])
        self._ip = (self._cartan_type.classical().cartan_matrix()*E).inverse()
        self._den0 = self._inner_pp(self._Lam+self._rho, self._Lam+self._rho)

    def highest_weight(self):
        """
        Returns the highest weight of the IntegrableRepresentation ``self``.

        EXAMPLES::

            sage: Lambda = RootSystem(['D',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: IntegrableRepresentation(Lambda[0]+2*Lambda[2]).highest_weight()
            Lambda[0] + 2*Lambda[2]

        """
        return self._Lam

    def weight_lattice(self):
        """
        Returns the Weight Lattice of ``self``.

        EXAMPLES::

            sage: v=IntegrableRepresentation(RootSystem(['E',6,1]).weight_lattice(extended=true).fundamental_weight(0))
            sage: v.weight_lattice()
            Extended weight lattice of the Root system of type ['E', 6, 1]
        """

        return self._P

    def root_lattice(self):
        """
        Returns the Root Lattice of ``self``.

        EXAMPLES::

            sage: v=IntegrableRepresentation(RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weight(0))
            sage: v.root_lattice()
            Root lattice of the Root system of type ['F', 4, 1]
        """
        return self._Q

    def null_root(self):
        """
        Returns the nullroot of ``self``.

        EXAMPLES::
        
            sage: Lambda = RootSystem(['G',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: v=IntegrableRepresentation(Lambda[0])
            sage: delta = v.null_root(); delta
            alpha[0] + 3*alpha[1] + 2*alpha[2]
        """
        return self._delta

    def rho(self):
        """
        Returns the Weyl vector `\\rho`.

        EXAMPLES::

            sage: Lambda = RootSystem(['G',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: IntegrableRepresentation(Lambda[0]).rho()
            Lambda[0] + Lambda[1] + Lambda[2]
        """
        return self._rho

    def level(self):
        """
        The level of the representation equals `(\Lambda|\delta)`. See [Kac] section 12.4.

        EXAMPLES::

            sage: Lambda = RootSystem(['G',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: [IntegrableRepresentation(Lambda[i]).level() for i in [0,1,2]]
            [1, 1, 2]

        """
        return ZZ(self._inner_pq(self._Lam, self._delta))

    def coxeter_number(self):
        """
        The dual Coxeter number is defined in [Kac] Chapter 6, and commonly
        denoted `h`.

        EXAMPLES::

            sage: Lambda = RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: v = IntegrableRepresentation(Lambda[0])
            sage: v.coxeter_number()
            12

        """
        return self._coxeter_number

    def dual_coxeter_number(self):
        """
        The dual Coxeter number is defined in [Kac] Chapter 6, and commonly
        denoted `h^\\vee`.

        EXAMPLES::

            sage: Lambda = RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: v = IntegrableRepresentation(Lambda[0])
            sage: v.dual_coxeter_number()
            9

        """
        return self._dual_coxeter_number

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return "The Integrable representation of %s with highest weight %s"%(self._cartan_type, self._Lam)

    def _inner_qq(self, qelt1, qelt2):
        """
        Symmetric form between two elements of the root lattice.
        
        EXAMPLES::
            
            sage: P = RootSystem(['F',4,1]).weight_lattice(extended=true)
            sage: Lambda = P.fundamental_weights()
            sage: v = IntegrableRepresentation(Lambda[0])
            sage: alpha = v.root_lattice().simple_roots()
            sage: Matrix([[v._inner_qq(alpha[i], alpha[j]) for j in v._index_set] for i in v._index_set])
            [   2   -1    0    0    0]
            [  -1    2   -1    0    0]
            [   0   -1    2   -1    0]
            [   0    0   -1    1 -1/2]
            [   0    0    0 -1/2    1]

        .. WARNING:

            If ``qelt1`` or ``qelt1`` accidentally gets coerced into
            the extended weight lattice, this will return an answer,
            and it will be wrong. To make this code robust, parents
            should be checked. This is not done since in the application
            the parents are known, so checking would unnecessarily slow
=            us down.
        """
        mc1 = qelt1.monomial_coefficients()
        mc2 = qelt2.monomial_coefficients()
        zero = ZZ.zero()
        return sum(mc1.get(i, zero) * mc2.get(j, zero)
                   * self._cartan_matrix[i,j] / self._eps[i]
                   for i in self._index_set for j in self._index_set)

    def _inner_pq(self, pelt, qelt):
        """
        Symmetric form between an element of the weight and root lattices.

        .. WARNING:

            If ``qelt`` accidentally gets coerced into the extended weight
            lattice, this will return an answer, and it will be wrong. To make
            this code robust, parents should be checked. This is not done
            since in the application the parents are known, so checking would
            unnecessarily slow us down.

        EXAMPLES::

            sage: P = RootSystem(['F',4,1]).weight_lattice(extended=true)
            sage: Lambda = P.fundamental_weights()
            sage: v = IntegrableRepresentation(Lambda[0])
            sage: alpha = v.root_lattice().simple_roots()
            sage: Matrix([[v._inner_pq(P(alpha[i]), alpha[j]) for j in v._index_set] for i in v._index_set])
            [   2   -1    0    0    0]
            [  -1    2   -1    0    0]
            [   0   -1    2   -1    0]
            [   0    0   -1    1 -1/2]
            [   0    0    0 -1/2    1]
            sage: P = RootSystem(['G',2,1]).weight_lattice(extended=true)
            sage: P = RootSystem(['G',2,1]).weight_lattice(extended=true)
            sage: Lambda = P.fundamental_weights()
            sage: v = IntegrableRepresentation(Lambda[0])
            sage: alpha = v.root_lattice().simple_roots()
            sage: Matrix([[v._inner_pq(Lambda[i],alpha[j]) for j in v._index_set] for i in v._index_set])
            [  1   0   0]
            [  0 1/3   0]
            [  0   0   1]
        """
        mcp = pelt.monomial_coefficients()
        mcq = qelt.monomial_coefficients()
        zero = ZZ.zero()
        return sum(mcp.get(i, zero) * mcq[i] / self._eps[i] for i in mcq) 

    def _inner_pp(self, pelt1, pelt2):
        """
        Symmetric form between an two elements of the weight lattice.

        EXAMPLES::

            sage: P = RootSystem(['G',2,1]).weight_lattice(extended=true)
            sage: Lambda = P.fundamental_weights()
            sage: v = IntegrableRepresentation(Lambda[0])
            sage: alpha = v.root_lattice().simple_roots()
            sage: Matrix([[v._inner_pp(Lambda[i],P(alpha[j])) for j in v._index_set] for i in v._index_set])
            [  1   0   0]
            [  0 1/3   0]
            [  0   0   1]
            sage: Matrix([[v._inner_pp(Lambda[i],Lambda[j]) for j in v._index_set] for i in v._index_set])
            [  0   0   0]
            [  0 2/3   1]
            [  0   1   2]

        """
        mc1 = pelt1.monomial_coefficients()
        mc2 = pelt2.monomial_coefficients()
        return sum(mc1.get(i,0)*self._ac[i]*mc2.get('delta',0) for i in self._index_set) +\
           sum(mc2.get(i,0)*self._ac[i]*mc1.get('delta',0) for i in self._index_set) +\
           sum(mc1.get(i,0)*mc2.get(j,0)*self._ip[i-1][j-1] for i in self._index_set_classical for j in self._index_set_classical)

    def to_weight(self, n):
        """
        INPUT:

        - ``n`` -- a tuple representing a weight.

        Returns the weight associated to the tuple ``n``.  If ``n``
        is the tuple `(n_1, n_2, \ldots)`, then the associated
        weight is `\Lambda - \sum_i n_i \\alpha_i`, where `\Lambda`
        is the weight of the representation.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: v = IntegrableRepresentation(2*Lambda[2])
            sage: v.to_weight((1,0,0))
            -2*Lambda[0] + Lambda[1] + 3*Lambda[2] - delta

        """
        alpha = self._P.simple_roots()
        I = self._index_set
        return self._Lam - self._P.sum(val * alpha[I[i]]
                                       for i,val in enumerate(n))

    def _from_weight_helper(self, mu, check=False):
        r"""
        Return ``(n[0], n[1], ...)`` such that ``mu = sum n[i]*alpha[i]``.

        INPUT:

        - ``mu`` -- an element in the root lattice

        .. TODO::

            Implement this as a section map of the inverse of the
            coercion from `Q \to P`.
        """
        mu = self._P(mu)
        n0 = mu.monomial_coefficients().get('delta', 0)
        mu0 = mu - n0*self._P.simple_root(self._cartan_type.special_node())
        ret = [n0] # This should be in ZZ because it is in the weight lattice
        mc_mu0 = mu0.monomial_coefficients()
        zero = ZZ.zero()
        for i in self._index_set_classical:
            # -1 for indexing
            ret.append( sum(self._cminv[i-1,j-1] * mc_mu0.get(j, zero)
                               for j in self._index_set_classical) )
        if check:
            return all(x in ZZ for x in ret)
        else:
            return tuple(ZZ(x) for x in ret)
    
    def from_weight(self, mu):
        """
        Return ``(n[0], n[1], ...)`` such that ``mu = Lam - sum n[i]*alpha[i]``

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: v = IntegrableRepresentation(2*Lambda[2])
            sage: v.to_weight((1,0,0))
            -2*Lambda[0] + Lambda[1] + 3*Lambda[2] - delta
            sage: delta = v.null_root()
            sage: v.from_weight(-2*Lambda[0] + Lambda[1] + 3*Lambda[2] - delta)
            (1, 0, 0)

        """
        return self._from_weight_helper(self._Lam - mu)

    def s(self, n, i):
        """
        Implements the `i`-th simple reflection in the internal
        representation of weights by tuples.

        EXAMPLES::

            sage: v = IntegrableRepresentation(RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weight(0))
            sage: [v.s((0,0,0),i) for i in v._index_set]
            [(1, 0, 0), (0, 0, 0), (0, 0, 0)]

        """
        ret = list(n) # This makes a copy
        ret[i] += self._Lam._monomial_coefficients.get(i, ZZ.zero())
        ret[i] -= sum(val * self._cartan_matrix[i,j] for j,val in enumerate(n))
        return tuple(ret)

    def to_dominant(self, n):
        """
        Return the dominant weight equivalent to ``n`` under the affine Weyl group.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: v = IntegrableRepresentation(3*Lambda[0])
            sage: n = v.to_dominant((13,11,7)); n
            (4, 3, 3)
            sage: v.to_weight(n)
            Lambda[0] + Lambda[1] + Lambda[2] - 4*delta
        """
        if n in self._ddict:
            return self._ddict[n]

        path = [n]
        alpha = self._P.simple_roots()
        next = True
        cur_wt = self.to_weight(n)

        while next:
            if path[-1] in self._ddict:
                path.append( self._ddict[path[-1]] )
                break

            next = False
            mc = cur_wt.monomial_coefficients()
            # Most weights are dense over the index set
            for i in self._index_set:
                if mc.get(i, 0) < 0:
                    m = self.s(path[-1], i)
                    if m in self._ddict:
                        path.append(self._ddict[m])
                    else:
                        cur_wt -= (m[i] - path[-1][i]) * alpha[i]
                        path.append(m)
                        next = True
                    break

        # We don't want any dominant weight to refer to itself in self._ddict
        #   as this leads to an infinite loop with self.m() when the dominant
        #   weight does not have a known multiplicity.
        v = path.pop()
        for m in path:
            self._ddict[m] = v
        return v

    def _freudenthal_roots_imaginary(self, nu):
        r"""
        It is assumed that ``nu`` is in `Q`. Returns the set of imaginary
        roots `\alpha \in \Delta^+` such that `\\nu - \alpha \in Q^+`.
        """
        l = self._from_weight_helper(nu)
        kp = min(floor(l[i] / self._a[i]) for i in self._index_set)
        delta = self._Q.null_root()
        return [u*delta for u in range(1, kp+1)]

    def _freudenthal_roots_real(self, nu):
        r"""
        It is assumed that ``nu`` is in `Q`. Returns the set of real positive
        roots `\alpha \in \Delta^+` such that `\\nu - \alpha \in Q^+`.
        """
        ret = []
        for al in self._classical_positive_roots:
            if all(x >= 0 for x in self._from_weight_helper(nu-al)):
                ret.append(al)
        for al in self._classical_roots:
            for ir in self._freudenthal_roots_imaginary(nu-al):
                ret.append(al+ir)
        return ret

    def _freudenthal_accum(self, nu, al):
        """
        Helper method for computing the Freudenthal formula.
        """
        ret = 0
        n = list(self._from_weight_helper(self._Lam - nu))
        ip = self._inner_pq(nu, al)
        n_shift = self._from_weight_helper(al)
        ip_shift = self._inner_qq(al, al)

        while all(val >= 0 for val in n):
            # Change in data by adding ``al`` to our current weight
            ip += ip_shift
            for i,val in enumerate(n_shift):
                n[i] -= val
            # Compute the multiplicity
            mk = self.m(tuple(n))
            ret += 2*mk*ip
        return ret

    def _m_freudenthal(self, n):
        """
        Compute the weight multiplicity using the Freudenthal multiplicity formula.
        """
        if min(n) < 0:
            return 0
        mu = self.to_weight(n)
        I = self._index_set
        al = self._Q._from_dict({I[i]: val for i,val in enumerate(n) if val},
                                remove_zeros=False)
        den = 2*self._inner_pq(self._Lam+self._P.rho(), al) - self._inner_qq(al, al)
        num = 0
        for al in self._freudenthal_roots_real(self._Lam - mu):
            num += self._freudenthal_accum(mu, al)
        for al in self._freudenthal_roots_imaginary(self._Lam - mu):
            num += self._classical_rank * self._freudenthal_accum(mu, al)
        if den == 0 and num == 0:
            print "_m_freudenthal","m: n=%s, num=den=0"%n.__repr__()
        try:
            return ZZ(num / den)
        except:
            return None

    def m(self, n):
        """
        INPUT:
        
        - ``n`` -- a tuple representing a weight `\mu`.

        Return the multiplicity of the weight `\mu` in ``self``, where
        `\mu = \Lambda - \sum_i n_i \\alpha_i`.

        EXAMPLES::

            sage: Lambda = RootSystem(['E',6,1]).weight_lattice(extended=true).fundamental_weights()
            sage: v = IntegrableRepresentation(Lambda[0])
            sage: u = v.highest_weight()-v.null_root()
            sage: v.from_weight(u)
            (1, 1, 2, 2, 3, 2, 1)
            sage: v.m(v.from_weight(u))
            6
        """
        # TODO: Make this non-recursive by implementing our own stack
        # The recursion follows:
        #   - m
        #   - _m_freudenthal
        #   - _freudenthal_accum
        if self._mdict.has_key(n):
            return self._mdict[n]
        elif self._ddict.has_key(n):
            self._mdict[n] = self.m(self._ddict[n])
        m = self.to_dominant(n)
        if self._mdict.has_key(m):
            return self._mdict[m]
        ret = self._m_freudenthal(m)
        if ret is not None:
            self._mdict[n] = ret
        else:
            print "m: error - failed to compute m%s"%n.__repr__()
        return ret

    def dominant_maximal_weights(self):
        """
        A weight `\\mu` is *maximal* if it has nonzero multiplicity but
        `\\mu+\\delta`` has multiplicity zero. There are a finite number
        of dominant maximal weights. Indeed, [Kac] Proposition 12.6
        shows that the dominant maximal weights are in bijection with
        the classical weights in `k\cdot F` where `F` is the fundamental
        alcove and `k` is the level. The construction used in this
        method is based on that Proposition.

        EXAMPLES::

            sage: Lambda = RootSystem(['C',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: IntegrableRepresentation(2*Lambda[0]).dominant_maximal_weights()
            [2*Lambda[0],
            Lambda[0] + Lambda[2] - delta,
            2*Lambda[1] - delta,
            Lambda[1] + Lambda[3] - 2*delta,
            2*Lambda[2] - 2*delta,
            2*Lambda[3] - 3*delta]
     
        """
        k = self.level()
        Lambda = self._P.fundamental_weights()
        def next_level(wt):
            return [wt + Lambda[i] for i in self._index_set_classical if (wt + Lambda[i]).level() <= k]
        R = RecursivelyEnumeratedSet([self._P.zero()], next_level)
        candidates = [x + (k - x.level())*Lambda[0] for x in list(R)]
        ret = []
        for x in candidates:
            if self._from_weight_helper(self._Lam-x, check=True):
                t = 0
                while self.m(self.from_weight(x-t*self._delta)) == 0:
                    t += 1
                ret.append(x-t*self._delta)
        return ret

    def string(self, max_weight, depth=12):
        """
        Return the list of multiplicities `m(\Lambda - k \delta)` where
        `\Lambda` is ``max_weight`` and `k` runs from `0` to ``depth``.

        INPUT:

        - ``max_weight`` -- a dominant maximal weight
        - ``depth`` -- (default: 12) the maximum value of `k`
            """
        ret = []
        delta = self._Q.null_root()
        cur_weight = max_weight
        for k in range(depth):
            ret.append(self.m( self.from_weight(cur_weight) ))
            cur_weight -= delta
        return ret

    # FIXME: Return a dictionary
    def strings(self, depth=12):
        """
        Return the set of dominant maximal weights of ``self``, together
        with the string coefficients for each.

        OPTIONAL:

        - ``depth`` -- a parameter indicating how far to push computations

        The depth parameter defaults to 12.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',1,1]).weight_lattice(extended=true).fundamental_weights()
            sage: IntegrableRepresentation(2*Lambda[0]).strings(depth=25)
            2*Lambda[0]: 1 1 3 5 10 16 28 43 70 105 161 236 350 501 722 1016 1431 1981 2741 3740 5096 6868 9233 12306 16357
            2*Lambda[1] - delta: 1 2 4 7 13 21 35 55 86 130 196 287 420 602 858 1206 1687 2331 3206 4368 5922 7967 10670 14193 18803

        """
        for max_weight in self.dominant_maximal_weights():
            s = self.string(max_weight, depth)
            print "%s:"%max_weight,
            for j in s:
                print j,
            print

    def modular_characteristic(self, mu=None):
        """
        The modular characteristic is a rational number introduced
        by Kac and Peterson (1984), required to interpret the string
        functions as Fourier coefficients of modular forms. See
        [Kac] Section 12.7. Let `k` be the level, and let `h^\\vee`
        be the dual Coxeter number. Then

        .. MATH::

            m_\Lambda = \\frac{|\Lambda+\\rho|^2}{2(k+h^\\vee)} - \\frac{|\\rho|^2}{2h^\\vee}

        If $\mu$ is a weight, then

        .. MATH::

            m_{\Lambda,\mu} = m_\Lambda - \\frac{|\mu|^2}{2k}.
            
        OPTIONAL:

        - ``mu`` -- a weight; or alternatively:
        - ``n`` -- a tuple representing a weight `\mu`.

        If not optional parameter is specified, The function returns `m_\Lambda`.
        If ``mu`` is specified, it returns `m_{\Lambda,\mu}`. You may use the
        tuple ``n`` to specify `\mu`. If you do this, `\mu` is
        `\Lambda - \sum_i n_i \\alpha_i`.
        
        EXAMPLES::

            sage: Lambda = RootSystem(['A',1,1]).weight_lattice(extended=true).fundamental_weights()
            sage: v = IntegrableRepresentation(3*Lambda[0]+2*Lambda[1])
            sage: [v.modular_characteristic(x) for x in v.dominant_maximal_weights()]
            [11/56, -1/280, 111/280]
        """
        if type(mu) is tuple:
            n = mu
        else:
            n = self.from_weight(mu)
        k = self.level()
        hd = self._dual_coxeter_number
        m_Lambda = self._inner_pp(self._Lam+self._rho, self._Lam+self._rho)/(2*(k+hd)) - \
                   self._inner_pp(self._rho,self._rho)/(2*hd)
        if n is None:
            return m_Lambda
        else:
            mu = self.to_weight(n)
            return m_Lambda-self._inner_pp(mu,mu)/(2*k)
