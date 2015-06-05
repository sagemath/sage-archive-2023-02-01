"""
Rational Cherednik Algebras
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class RationalCherednikAlgebra(CombinatorialFreeModule):
    r"""
    Rational Cherednik algebra.

    INPUT:

    - ``ct`` -- a Cartan type
    - ``c`` -- the parameter `c`
    - ``t`` -- the parameter `t`
    - ``base_ring`` -- (optional) the base ring
    - ``prefix`` -- (default: ``'R'``) the prefix
    """
    @staticmethod
    def __classcall_private__(cls, ct, c, t, base_ring=None, prefix=('a', 's', 'ac')):
        """
        Normalize input to ensure a unique representation.
        """
        ct = CartanType(ct)
        if base_ring is None:
            base_ring = t.parent()
        t = base_ring(t)
        c = base_ring(c)
        return super(RationalCherednikAlgebra, cls).__classcall__(cls, ct, c, t, base_ring, tuple(prefix))

    def __init__(self, ct, c, t, base_ring, prefix):
        r"""
        Initialize ``self``.
        """
        self._c = c
        self._t = t
        self._cartan_type = ct
        self._weyl = WeylGroup(ct, prefix=prefix[1])
        self._hd = IndexedFreeAbelianMonoid(ct.index_set(), prefix=prefix[0],
                                            bracket=False)
        self._h = IndexedFreeAbelianMonoid(ct.index_set(), prefix=prefix[2],
                                           bracket=False)
        indices = DisjointUnionEnumeratedSets([self._hd, self._weyl, self._h])
        CombinatorialFreeModule.__init__(self, base_ring, indices,
                                         category=Algebras(base_ring).WithBasis().Graded())

    @lazy_attribute
    def _reflections(self):
        """
        A dictionary of reflections to a pair of the associated root
        and coroot.
        """
        d = {}
        for s in self._weyl.reflections().keys():
            r = s.reflection_to_root()
            d[s] = (r, r.associated_coroot())
        return d

    def _repr_(self):
        r"""
        EXAMPLES ::

            sage: R = RationalCherednikAlgebra(['A',4], 2, 1, QQ)
            Rational Cherednik Algebra of type ['A', 4] with c=2 and t=1
             over Rational Field
        """
        return "Rational Cherednik Algebra of type {} with c={} and t={} over {}".format(
                        self._cartan_type, self._c, self._t, self.base_ring())

    def _repr_term(self, t):
        """
        Return a string representation of the term ``t``.
        """
        r = []
        if t[0] != self._hd.one():
            r.append(t[0])
        if t[1] != self._weyl.one():
            r.append(t[1])
        if t[2] != self._h.one():
            r.append(t[2])
        if not r:
            return 'I'
        return '*'.join(repr(x) for x in r)

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.
        """
        p = self
        keys  = ['a'+str(i) for i in self._cartan_type.index_set()]
        keys += ['s'+str(i) for i in self._cartan_type.index_set()]
        keys += ['ac'+str(i) for i in self._cartan_type.index_set()]
        one = self.base_ring().one()
        def gen_map(k):
            if k[0] == 's':
                i = int(k[1:])
                return self.monomial( (self._hd.one(),
                                       self._weyl.group_generators()[i],
                                       self._h.one()) )
            if k[1] == 'c':
                i = int(k[2:])
                return self.monomial( (self._hd.one(),
                                       self._weyl.one(),
                                       self._h.monoid_generators()[i]) )

            i = int(k[1:])
            return self.monomial( (self._hd.monoid_generators()[i],
                                   self._weyl.one(),
                                   self._h.one()) )
        return Family(keys, gen_map)

    @cached_method
    def one_basis(self):
        """
        Return the index of the element `1`.
        """
        return (self._hd.one(), self._weyl.one(), self._h.one())

    def product_on_basis(self, left, right):
        r"""
        Return ``left`` multiplied by ``right`` in ``self``.

        Let `\alpha \in \mathfrak{h}` and `\alpha^{\vee} \in \mathfrak{h}^*`.
        Let `w \in W` and `\cdot denote the natural action of `w` on
        `\mathfrak{h}` and `\mathfrak{h}^*`. We have the following relations:

        .. MATH::

            \begin{aligned}
            w \alpha & = (w \cdot \alpha) w,
            \\ \alpha^{\vee} w & = w (w^{-1} \cdot \alpha^{\vee}),
            \\ \alpha \alpha^{\vee} & = \alpha^{\vee} \alpha
            + t \langle \alpha^{\vee}, \alpha \rangle
            + \sum_{s \in \mathcal{S}} c_s \frac{\langle \alpha^{\vee},
            \alpha_s \rangle \langle \alpha^{\vee}_s, \alpha \rangle}{
            \langle \alpha^{\vee}, \alpha \rangle} s,

        where `\mathcal{S}` are the set of all reflections of `W` and
        `\alpha_s` and `\alpha_s^{\vee}` are the associated root and
        coroot of `s`.
        """
        # Make copies of the internal dictionaries
        dl = dict(left[2]._monomial)
        dr = dict(right[0]._monomial)

        # If there is nothing to commute
        if not dl and not dr:
            return self.monomial((left[0], left[1] * right[1], right[2]))

        hd_one = self._hd.one() # root - a
        h_one = self._h.one() # coroot - ac
        R = self.base_ring()
        I = self._cartan_type.index_set()
        P = PolynomialRing(R, 'x', len(I))
        G = P.gens()
        gens_dict = {a:G[i] for i,a in enumerate(I)}
        Q = RootSystem(self._cartan_type).root_lattice()
        alpha = Q.simple_roots()
        alphacheck = Q.simple_coroots()

        def commute_w_hd(w, al): # al is given as a dictionary
            ret = P.one()
            for k in al:
                x = sum(c * gens_dict[i] for i,c in alpha[k].weyl_action(w))
                ret *= x**al[k]
            ret = ret.dict()
            for k in ret:
                yield (self._hd({I[i]: e for i,e in enumerate(k) if e != 0}), ret[k])

        # Do Lac Ra if they are both non-trivial
        if dl and dr:
            il = dl.keys()[0]
            ir = dr.keys()[0]

            # Compute the commutator
            terms = self._product_coroot_root(il, ir)

            # remove the generator from the elements
            dl[il] -= 1
            if dl[il] == 0:
                del dl[il]
            dr[ir] -= 1
            if dr[ir] == 0:
                del dr[ir]

            # We now commute right roots past the left reflections: s Ra = Ra' s
            cur = self._from_dict({ (hd, s*right[1], right[2]): c * cc
                                    for s,c in terms
                                    for hd, cc in commute_w_hd(s, dr) })
            cur = self.monomial( (left[0], left[1], self._h(dl)) ) * cur

            # Add back in the commuted h and hd elements
            rem = self.monomial( (left[0], left[1], self._h(dl)) )
            rem = rem * self.monomial( (self._hd({ir:1}), self._weyl.one(),
                                        self._h({il:1})) )
            rem = rem * self.monomial( (self._hd(dr), right[1], right[2]) )

            return cur + rem

        if dl:
            # We have La Ls Lac Rs Rac,
            #   so we must commute Lac Rs = Rs Lac'
            #   and obtain La (Ls Rs) (Lac' Rac)
            ret = P.one()
            for k in dl:
                x = sum(c * gens_dict[i]
                        for i,c in alphacheck[k].weyl_action(right[1].reduced_word(),
                                                             inverse=True))
                ret *= x**dl[k]
            ret = ret.dict()
            w = left[1]*right[1]
            return self._from_dict({ (left[0], w,
                                      self._h({I[i]: e for i,e in enumerate(k)
                                               if e != 0}) * right[2]
                                     ): ret[k]
                                     for k in ret })

        # Otherwise dr is non-trivial and we have La Ls Ra Rs Rac,
        #   so we must commute Ls Ra = Ra' Ls
        w = left[1]*right[1]
        return self._from_dict({ (left[0] * hd, w, right[2]): c
                                 for hd, c in commute_w_hd(left[1], dr) })

    @cached_method
    def _product_coroot_root(self, i, j):
        r"""
        Return the product `\alpha^{\vee}_i \alpha_j`
        """
        Q = RootSystem(self._cartan_type).root_lattice()
        ac = Q.simple_coroot(i)
        al = Q.simple_root(j)

        R = self.base_ring()
        terms = [( self._weyl.one(), self._t * R(ac.scalar(al)) )]
        for s in self._reflections:
            p = self._reflections[s] # p[0] is the root, p[1] is the coroot
            terms.append(( s, self._c * R(ac.scalar(p[0]) * p[1].scalar(al)
                                           / p[1].scalar(p[0])) ))
        return tuple(terms)

    def degree_on_basis(self, m):
        """
        Return the degree on the monomial indexed by ``m``.
        """
        return m[0].length() - m[2].length()

    @cached_method
    def trivial_idempotent(self):
        """
        Return the trivial idempotent of ``self``.

        Let `e = |W|^{-1} \sum_{w \in W} w` denote the trivial idempotent.
        We construct the spherical Cherednik algebra by `U(W) = e H(W) e`,
        where `H(W)` is the rational Cherednik algebra.
        """
        coeff = self.base_ring()(~self._weyl.cardinality())
        hd_one = self._hd.one() # root - a
        h_one = self._h.one() # coroot - ac
        return self._from_dict({(hd_one, w, h_one): coeff for w in self._weyl},
                               remove_zeros=False)

    @cached_method
    def deformed_euler(self):
        """
        Return the element `eu_k`.
        """
        I = self._cartan_type.index_set()
        G = self.algebra_generators()
        cm = ~CartanMatrix(self._cartan_type)
        n = len(I)
        ac = [G['ac'+str(i)] for i in I]
        la = [sum(cm[i,j]*G['a'+str(I[i])] for i in range(n)) for j in range(n)]
        return self.sum(ac[i]*la[i] for i in range(n))

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.
        """
        G = self.algebra_generators()
        i = str(self._cartan_type.index_set()[0])
        return G['a'+i] + 2*G['s'+i] + 3*G['ac'+i]

    def some_elements(self):
        """
        Return some elements of ``self``.
        """
        ret = [self.zero(), self.one(), self.an_element()]
        ret += list(self.algebra_generators())
        return ret

    #class Element(CombinatorialFreeModuleElement):

