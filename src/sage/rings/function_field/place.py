"""
Places

Sage can find and compute with places on global function fields.

EXAMPLES:

All rational places of the function field can be computed::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
    sage: L.places()
    [Place (1/x, 1/x^3*y^2 + 1/x),
     Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1),
     Place (x, y)]

The residue field associated with a place is given as an extension of the
constant field::

    sage: F.<x> = FunctionField(GF(2))
    sage: O = F.maximal_order()
    sage: p = O.ideal(x^2 + x + 1).place()
    sage: k, fr_k, to_k = p.residue_field()
    sage: k
    Finite Field in z2 of size 2^2

The isomorphisms are between the valuation ring and the residue field::

    sage: fr_k
    Ring morphism:
      From: Finite Field in z2 of size 2^2
      To:   Valuation ring at Place (x^2 + x + 1)
    sage: to_k
    Ring morphism:
      From: Valuation ring at Place (x^2 + x + 1)
      To:   Finite Field in z2 of size 2^2

AUTHORS:

- Kwankyu Lee (2017-04-30): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

from sage.arith.all import lcm

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.sage_object import richcmp

from sage.modules.free_module_element import vector

lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.rings.function_field.divisor', 'prime_divisor')

def is_Place(x):
    """
    Return True if x is a place of a function field.
    """
    return isinstance(x, FunctionFieldPlace)

class FunctionFieldPlace(Element):
    """
    Places of function fields.
    """
    def __init__(self, field, prime):
        """
        Initialize the place.

        INPUT:

        - ``field`` -- a function field

        - ``prime`` -- a prime ideal associated with the place

        """
        Element.__init__(self, field.place_set())

        self._field = field
        self._prime = prime

    def __hash__(self):
        """
        Return the hash of the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y>=K.extension(Y^3 + x + x^3*Y)
            sage: p = L.places_finite()[0]
            sage: {p: 1}
            {Place (x, y): 1}
        """
        return hash((self._field, self._prime))

    def _repr_(self):
        """
        Return the string representation of the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p = L.places_finite()[0]
            sage: p
            Place (x, y)
        """
        try:
            gens = self._prime.gens_two()
        except AttributeError:
            gens = self._prime.gens()
        gens_str = ', '.join(repr(g) for g in gens)
        return "Place ({})".format(gens_str)

    def _richcmp_(self, other, op):
        """
        Compare the place with ``other`` place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: p1, p2, p3 = L.places()[:3]
            sage: p1 < p2
            True
            sage: p2 < p1
            False
            sage: p1 == p3
            False
        """
        # First compare the rings. In effect, infinite
        # places are compared less than finite places so that
        # they are listed first.
        return richcmp((self._prime.ring(), self._prime),
                       (other._prime.ring(), other._prime), op)

    def _acted_upon_(self, other, self_on_left):
        """
        Define integer multiplication upon the prime divisor
        of the place on the left.

        The output is a divisor.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x+1,y)
            sage: P = I.place()
            sage: -3*P + 5*P
            2*Place (x + 1, y)
        """
        if self_on_left:
            raise TypeError("only left multiplication by integers is allowed")
        return other * self.divisor()

    def _add_(self, other):
        """
        Return the divisor that is the sum of the place and ``other``.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p1, p2, p3 = L.places()[:3]
            sage: p1 + p2 + p3
            Place (1/x, 1/x^3*y^2 + 1/x)
             + Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             + Place (x, y)
        """
        return prime_divisor(self._field, self) + other

    def __radd__(self, other):
        """
        Return the prime divisor of the place if ``other`` is zero.

        This is only to support the ``sum`` function, that adds
        the argument to initial (int) zero.

        EXAMPLES::

            sage: k.<a>=GF(2)
            sage: K.<x>=FunctionField(k)
            sage: sum(K.places_finite())
            Place (x) + Place (x + 1)

        .. NOTE:

        This does not work though::

            sage: 0 + K.place_infinite()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +: ...

        The reason is that the ``0`` is a Sage integer, for which
        the coercion system applies.
        """
        if other == 0:
            return prime_divisor(self._field, self)
        raise NotImplementedError

    def function_field(self):
        """
        Return the function field to which the place belongs.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p = L.places()[0]
            sage: p.function_field() == L
            True
        """
        return self._field

    def prime_ideal(self):
        """
        Return the prime ideal associated with the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: p = L.places()[0]
            sage: p.prime_ideal()
            Ideal (1/x,1/x^3*y^2 + 1/x) of Maximal infinite order of Function field
            in y defined by y^3 + x^3*y + x
        """
        return self._prime

    def divisor(self, multiplicity=1):
        """
        Return the prime divisor corresponding to the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x+1,y)
            sage: P = I.place()
            sage: P.divisor()
            Place (x + 1, y)
        """
        return prime_divisor(self._field, self, multiplicity)

class FunctionFieldPlace_rational(FunctionFieldPlace):
    """
    Places of rational function field.
    """
    def degree(self):
        """
        Return the degree of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: i = O.ideal(x^2+x+1)
            sage: p = i.place()
            sage: p.degree()
            2
        """
        if self.is_infinite_place():
            return 1
        else:
            return self._prime.gen().numerator().degree()

    def is_infinite_place(self):
        """
        Return True if the place is at infinite.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: F.places()
            [Place (1/x), Place (x), Place (x + 1)]
            sage: [p.is_infinite_place() for p in F.places()]
            [True, False, False]
        """
        F = self.function_field()
        return self.prime_ideal().ring() == F.maximal_order_infinite()

    def local_uniformizer(self):
        """
        Return a local uniformizer of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: F.places()
            [Place (1/x), Place (x), Place (x + 1)]
            sage: [p.local_uniformizer() for p in F.places()]
            [1/x, x, x + 1]
        """
        return self.prime_ideal().gen()

    def residue_field(self, name=None):
        """
        Return the residue field of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: p = O.ideal(x^2 + x + 1).place()
            sage: k, fr_k, to_k = p.residue_field()
            sage: k
            Finite Field in z2 of size 2^2
            sage: fr_k
            Ring morphism:
              From: Finite Field in z2 of size 2^2
              To:   Valuation ring at Place (x^2 + x + 1)
            sage: to_k
            Ring morphism:
              From: Valuation ring at Place (x^2 + x + 1)
              To:   Finite Field in z2 of size 2^2
        """
        return self.valuation_ring().residue_field(name=name)

    def _residue_field(self, name=None):
        """
        Return the residue field of the place along with the maps from
        and to it.

        INPUT:

        - ``name`` -- a string; name of the generator of the residue field

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: i = O.ideal(x^2+x+1)
            sage: p = i.place()
            sage: R, fr, to = p._residue_field()
            sage: R
            Finite Field in z2 of size 2^2
            sage: [fr(e) for e in R.list()]
            [0, x, x + 1, 1]
            sage: to(x*(x+1)) == to(x) * to(x+1)
            True
        """
        F = self._field
        prime = self.prime_ideal()

        if self.is_infinite_place():
            K = F.constant_base_field()

            def from_K(e):
                return F(e)

            def to_K(f):
                n = f.numerator()
                d = f.denominator()

                n_deg = n.degree()
                d_deg =d.degree()

                if n_deg < d_deg:
                    return K(0)
                elif n_deg == d_deg:
                    return n.lc() / d.lc()
                else:
                    raise TypeError("not in the valuation ring")
        else:
            O = F.maximal_order()
            K, from_K, _to_K = O._residue_field(prime, name=name)

            def to_K(f):
                if f in O: # f.denominator() is 1
                    return _to_K(f.numerator())
                else:
                    d = F(f.denominator())
                    n = d * f

                    nv = prime.valuation(O.ideal(n))
                    dv = prime.valuation(O.ideal(d))

                    if nv > dv:
                        return K(0)
                    elif dv > nv:
                        raise TypeError("not in the valuation ring")

                    s = ~prime.gen()
                    rd = d * s**dv # in O but not in prime
                    rn = n * s**nv # in O but not in prime
                    return to_K(rn) / to_K(rd)

        return K, from_K, to_K

    @cached_method
    def valuation_ring(self):
        """
        Return the valuation ring at the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: p.valuation_ring()
            Valuation ring at Place (x, x*y)
        """
        from .valuation_ring import FunctionFieldValuationRing_global
        return FunctionFieldValuationRing_global(self._field, self)

class FunctionFieldPlace_global(FunctionFieldPlace):
    """
    Places of function fields
    """
    def place_below(self):
        """
        Return the place lying below the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: OK = K.maximal_order()
            sage: OL = L.maximal_order()
            sage: p = OK.ideal(x^2 + x + 1)
            sage: dec = OL.decomposition(p)
            sage: q = dec[0][0].place()
            sage: q.place_below()
            Place (x^2 + x + 1)
        """
        return self.prime_ideal().prime_below().place()

    def relative_degree(self):
        """
        Return the relative degree of the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: OK = K.maximal_order()
            sage: OL = L.maximal_order()
            sage: p = OK.ideal(x^2 + x + 1)
            sage: dec = OL.decomposition(p)
            sage: q = dec[0][0].place()
            sage: q.relative_degree()
            1
        """
        return self._prime._relative_degree

    def degree(self):
        """
        Return the degree of the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: OK = K.maximal_order()
            sage: OL = L.maximal_order()
            sage: p = OK.ideal(x^2 + x + 1)
            sage: dec = OL.decomposition(p)
            sage: q = dec[0][0].place()
            sage: q.degree()
            2
        """
        return self.relative_degree() * self.place_below().degree()

    def is_infinite_place(self):
        """
        Return True if the place is above the unique infinite place
        of the underlying rational function field.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: pls = L.places()
            sage: [p.is_infinite_place() for p in pls]
            [True, True, False]
            sage: [p.place_below() for p in pls]
            [Place (1/x), Place (1/x), Place (x)]
        """
        return self.place_below().is_infinite_place()

    def local_uniformizer(self):
        """
        Return an element of the function field that has a simple zero
        at the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: pls = L.places()
            sage: [p.local_uniformizer().valuation(p) for p in pls]
            [1, 1, 1, 1, 1]
        """
        gens = self._prime.gens()
        for g in gens:
            if g.valuation(self) == 1:
                return g
        assert False, "Internal error"

    def gaps(self):
        """
        Return the gap sequence for the place.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x,y).place()
            sage: p.gaps() # a Weierstrass place
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: [p.gaps() for p in L.places()]
            [[1, 2, 4], [1, 2, 4], [1, 2, 4]]
        """
        if self.degree() == 1:
            return self._gaps_rational()
        else:
            return self._gaps_wronskian()

    def _gaps_rational(self):
        """
        Return the gap sequence for the rational place.

        This method computes the gap numbers using the definition of gap
        numbers. The dimension of the multiple of the prime divisor
        supported at the place is computed by Hess' algorithm.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x,y).place()
            sage: p.gaps() # a Weierstrass place
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: [p.gaps() for p in L.places()]
            [[1, 2, 4], [1, 2, 4], [1, 2, 4]]
        """
        F = self._field
        n = F.degree()
        O = F.maximal_order()
        Oinf = F.maximal_order_infinite()

        R = O._module_base_ring._ring
        one = R.one()

        # Hess' Riemann-Roch basis algorithm stripped down for gaps computation
        def dim_RR(M):
            den = lcm([e.denominator() for e in M.list()])
            mat = matrix(R, M.nrows(), [(den*e).numerator() for e in M.list()])

            # initialise pivot_row and conflicts list
            pivot_row = [[] for i in range(n)]
            conflicts = []
            for i in range(n):
                bestp = -1
                best = -1
                for c in range(n):
                    d = mat[i,c].degree()
                    if d >= best:
                        bestp = c
                        best = d

                if best >= 0:
                    pivot_row[bestp].append((i,best))
                    if len(pivot_row[bestp]) > 1:
                        conflicts.append(bestp)

            # while there is a conflict, do a simple transformation
            while conflicts:
                c = conflicts.pop()
                row = pivot_row[c]
                i,ideg = row.pop()
                j,jdeg = row.pop()

                if jdeg > ideg:
                    i,j = j,i
                    ideg,jdeg = jdeg,ideg

                coeff = - mat[i,c].lc() / mat[j,c].lc()
                s = coeff * one.shift(ideg - jdeg)

                mat.add_multiple_of_row(i, j, s)

                row.append((j,jdeg))

                bestp = -1
                best = -1
                for c in range(n):
                    d = mat[i,c].degree()
                    if d >= best:
                        bestp = c
                        best = d

                if best >= 0:
                    pivot_row[bestp].append((i,best))
                    if len(pivot_row[bestp]) > 1:
                        conflicts.append(bestp)

            dim = 0
            for j in range(n):
                i,ideg = pivot_row[j][0]
                k = den.degree() - ideg + 1
                if k > 0:
                    dim += k
            return dim

        V,fr,to = F.vector_space()

        prime_inv = ~ self.prime_ideal()
        I = O.ideal(1)
        J = Oinf.ideal(1)

        B = matrix([to(b) for b in J.gens_over_base()])
        C = matrix([to(v) for v in I.gens_over_base()])

        prev = dim_RR(C * B.inverse())
        gaps = []
        g = F.genus()
        i = 1
        if self.is_infinite_place():
            while len(gaps) < g:
                J = J * prime_inv
                B = matrix([to(b) for b in J.gens_over_base()])
                dim = dim_RR(C * B.inverse())
                if dim == prev:
                    gaps.append(i)
                else:
                    prev = dim
                i += 1
        else: # self is a finite place
            Binv = B.inverse()
            while len(gaps) < g:
                I = I * prime_inv
                C = matrix([to(v) for v in I.gens_over_base()])
                dim = dim_RR(C * Binv)
                if dim == prev:
                    gaps.append(i)
                else:
                    prev = dim
                i += 1

        return gaps

    def _gaps_wronskian(self):
        """
        Return the gap sequence for the place.

        This method implements the local version of Hess' Algorithm 30 of [Hes2002b]_
        based on the Wronskian determinant.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: O = L.maximal_order()
            sage: p = O.ideal(x,y).place()
            sage: p._gaps_wronskian() # a Weierstrass place
            [1, 2, 4]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3 * Y + x)
            sage: [p._gaps_wronskian() for p in L.places()]
            [[1, 2, 4], [1, 2, 4], [1, 2, 4]]
        """
        F = self.function_field()
        R,fr_R,to_R = self._residue_field()
        hasse = F.hasse_derivation()

        sep = self.local_uniformizer()

        # a differential divisor satisfying
        # v_p(W) = 0 for the place p
        W = sep.differential().divisor()

        # Step 3:
        basis = W._basis()
        d = len(basis)
        M = matrix([to_R(b) for b in basis])
        if M.rank() == 0:
            return []

        # Steps 4, 5, 6, 7:
        e = 1
        gaps = [1]
        while M.nrows() < d:
            row = vector([to_R(hasse._derive(basis[i], e, sep)) for i in range(d)])
            if not row in M.row_space():
                M = matrix(M.rows() + [row])
                M.echelonize()
                gaps.append(e+1)
            e += 1

        return gaps

    def residue_field(self, name=None):
        """
        Return the residue field of the place.

        INPUT:

        - ``name`` -- a string; name of the generator of the residue field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: k, fr_k, to_k = p.residue_field()
            sage: k
            Finite Field of size 2
            sage: fr_k
            Ring morphism:
              From: Finite Field of size 2
              To:   Valuation ring at Place (x, x*y)
            sage: to_k
            Ring morphism:
              From: Valuation ring at Place (x, x*y)
              To:   Finite Field of size 2
            sage: to_k(y)
            Traceback (most recent call last):
            ...
            TypeError: y fails to convert into the map's domain
            Valuation ring at Place (x, x*y)...
            sage: to_k(1/y)
            0
            sage: to_k(y/(1+y))
            1
        """
        return self.valuation_ring().residue_field(name=name)

    @cached_method
    def _residue_field(self, name=None):
        """
        Return the residue field of the place along with the functions
        mapping from and to it.

        INPUT:

        - ``name`` -- a string; name of the generator of the residue field

        If name is not given, it defaults to ``'a'``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: k,fr_k,to_k = p._residue_field()
            sage: k
            Finite Field of size 2
            sage: [fr_k(e) for e in k]
            [0, 1]

            sage: K.<x> = FunctionField(GF(9)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + Y - x^4)
            sage: p = L.places()[-1]
            sage: p.residue_field()
            (Finite Field in z2 of size 3^2, Ring morphism:
               From: Finite Field in z2 of size 3^2
               To:   Valuation ring at Place (x + 1, y + 2*z2), Ring morphism:
               From: Valuation ring at Place (x + 1, y + 2*z2)
               To:   Finite Field in z2 of size 3^2)
        """
        F = self._field
        prime = self.prime_ideal()

        if self.is_infinite_place():
            _F, from_F, to_F  = F._inversion_isomorphism()
            _prime = prime._ideal
            _place = _prime.place()

            K, _from_K, _to_K = _place._residue_field(name=name)

            from_K = lambda e: from_F(_from_K(e))
            to_K = lambda f: _to_K(to_F(f))
            return K, from_K, to_K

        O = F.maximal_order()
        Obasis = O.basis()

        M = prime.hnf()
        R = M.base_ring() # univariate polynomial ring
        n = M.nrows() # extension degree of the function field

        degs = [M[i,i].degree() for i in range(n)]
        deg = sum(degs) # degree of the place

        # Step 1: construct a vector space representing the residue field
        k = F.constant_base_field()
        V = k**deg

        def to_V(e):
            """
            An example to show the idea: Suppose that

                    [x 0 0]
                M = [0 1 0] and v = (x^10, x^7 + x^3, x^7 + x^4 + x^3 + 1)
                    [1 0 1]

            Then to_V(e) = [1]
            """
            v = O._coordinate_vector(e)
            vec = []
            for i in reversed(range(n)):
                q,r = v[i].quo_rem(M[i,i])
                v -= q * M[i]
                for j in range(degs[i]):
                    vec.append(r[j])
            return vector(vec)

        def fr_V(vec): # to_O
            vec = vec.list()
            pos = 0
            e = F(0)
            for i in reversed(range(n)):
                if degs[i] == 0: continue
                else:
                    end = pos + degs[i]
                    e += R(vec[pos:end]) * Obasis[i]
                    pos = end
            return e

        # Step 2: find a generator of the residue field
        def candidates():
            # Trial 1: this suffices for places obtained from Kummers' theorem

            # Note that a = O._kummer_gen is a simple generator of O/prime over
            # o/p. If b is a simple generator of o/p over the constant base field
            # k, then the set a + k * b contains a simple generator of O/prime
            # over k (as there are finite number of intermediate fields).
            a = O._kummer_gen
            if a is not None:
                K,fr_K,_ = self.place_below().residue_field()
                b = fr_K(K.gen())
                for c in reversed(k.list()):
                    yield a + c * b

            # Trial 2: basis elements of the maximal order
            for gen in reversed(Obasis):
                yield gen

            import itertools

            # Trial 3: exhaustive search in O using only polynomials
            # with coefficients 0 or 1
            for d in range(deg):
                G = itertools.product(itertools.product([0,1],repeat=d+1), repeat=n)
                for g in G:
                    gen = sum([R(c1)*c2 for c1,c2 in zip(g, Obasis)])
                    yield gen

            # Trial 4: exhaustive search in O using all polynomials
            for d in range(deg):
                G = itertools.product(R.polynomials(max_degree=d), repeat=n)
                for g in G:
                    # discard duplicate cases
                    if max(c.degree() for c in g) != d: continue
                    for j in range(n):
                        if g[j] != 0: break
                    if g[j].leading_coefficient() != 1: continue

                    gen = sum([c1*c2 for c1,c2 in zip(g, Obasis)])
                    yield gen

        for gen in candidates():
            g = F.one()
            m = []
            for i in range(deg):
                m.append(to_V(g))
                g *= gen
            mat = matrix(m)
            if mat.rank() == deg:
                break

        # Step 3: compute the minimal polynomial of g
        min_poly = R((-mat.solve_left(to_V(g))).list() + [1])

        if deg > 1:
            if name is None:
                name = 'a' # default

            # Step 4: construct the finite field
            K = k.extension(min_poly, name=name)

            # Step 5: compute the matrix of change of basis
            C = mat.inverse()

            # Step 6: construct an isomorphism
            def from_K(e):
                return fr_V(vector(e.list()) * mat)
        else: # deg == 1
            # Step 4: construct the finite field
            K = k

            # Step 5: compute the matrix of change of basis
            C = mat.inverse()

            # Step 6: construct an isomorphism
            def from_K(e):
                return fr_V(vector([e]) * mat)

        p = prime.prime_below().gen().numerator()
        beta = prime._beta
        alpha = ~p * sum(c1*c2 for c1,c2 in zip(beta, O.basis()))
        alpha_powered_by_ramification_index = alpha ** prime._ramification_index

        def to_K(f):
            if not f in O:
                den = O.coordinate_vector(f).denominator()
                num = den * f

                # s powered by the valuation of den at the prime
                alpha_power = alpha_powered_by_ramification_index ** den.valuation(p)
                rn = num * alpha_power # in O
                rd = den * alpha_power # in O but not in prime

                # Note that rn is not in O if and only if f is
                # not in the valuation ring. Hence f is in the
                # valuation ring if and only if this procedure
                # does not fall into an infinite loop.
                return to_K(rn) / to_K(rd)

            e = (to_V(f)*C).list()
            if deg > 1:
                return K(e)
            else: # len(e) == 1
                return K(e[0])

        return K, from_K, to_K

    @cached_method
    def valuation_ring(self):
        """
        Return the valuation ring at the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: p.valuation_ring()
            Valuation ring at Place (x, x*y)
        """
        from .valuation_ring import FunctionFieldValuationRing_global

        return FunctionFieldValuationRing_global(self._field, self)

class PlaceSet(Parent):
    """
    Sets of Places of function fields.
    """
    def __init__(self, field):
        """
        Initialize the set of places of the function ``field``.

        INPUT:

        - ``field`` -- a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: L.place_set()
            Set of places of Function field in y defined by y^3 + x^3*y + x
        """
        Parent.__init__(self)

        self._field = field

    def _repr_(self):
        """
        Return the string representation of the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: L.place_set()
            Set of places of Function field in y defined by y^3 + x^3*y + x
        """
        return "Set of places of {}".format(self._field)


