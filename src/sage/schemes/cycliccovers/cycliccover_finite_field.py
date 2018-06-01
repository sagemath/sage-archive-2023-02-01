"""
Cyclic covers curves over a finite field

EXAMPLES::
    sage: p = 5;
    sage: x = PolynomialRing(GF(p),"x").gen();
    sage: C = CyclicCover(6, x^6 + 1);
    sage: C.frobenius_polynomial()
    x^20 + 50*x^18 + 1125*x^16 + 15000*x^14 + 131250*x^12 + 787500*x^10 + 3281250*x^8 + 9375000*x^6 + 17578125*x^4 + 19531250*x^2 + 9765625
    sage: R.<t> = PowerSeriesRing(Integers());
    sage: C.projective_closure().zeta_series(4,t)
    1 + 6*t + 81*t^2 + 456*t^3 + 3456*t^4 + O(t^5)
    sage: C.frobenius_polynomial().reverse()(t)/((1-t)*(1-p*t)) + O(t^5)
    1 + 6*t + 81*t^2 + 456*t^3 + 3456*t^4 + O(t^5)


    sage: p = 49999;
    sage: x = PolynomialRing(GF(p),"x").gen();
    sage: CyclicCover(5, x^5 + x ).frobenius_polynomial()
    x^12 + 299994*x^10 + 37498500015*x^8 + 2499850002999980*x^6 + 93742500224997000015*x^4 + 1874812507499850001499994*x^2 + 15623125093747500037499700001
    sage: CyclicCover(5, 2*x^5 + x).frobenius_polynomial()
    x^12 + 299994*x^10 + 37498500015*x^8 + 2499850002999980*x^6 + 93742500224997000015*x^4 + 1874812507499850001499994*x^2 + 15623125093747500037499700001

    sage: p = 107;
    sage: x = PolynomialRing(GF(p),"x").gen();
    sage: CyclicCover(2, x^5 + x).frobenius_matrix()
    [              O(107^2)      89*107 + O(107^2)               O(107^2)               O(107^2)]
    [     89*107 + O(107^2)               O(107^2)               O(107^2)               O(107^2)]
    [              O(107^2)               O(107^2)               O(107^2) 105 + 5*107 + O(107^2)]
    [              O(107^2)               O(107^2) 89 + 53*107 + O(107^2)               O(107^2)]
    sage: CyclicCover(2, 3*x^5 + x).frobenius_matrix()
    [              O(107^2)      14*107 + O(107^2)               O(107^2)               O(107^2)]
    [     69*107 + O(107^2)               O(107^2)               O(107^2)               O(107^2)]
    [              O(107^2)               O(107^2)               O(107^2) 61 + 58*107 + O(107^2)]
    [              O(107^2)               O(107^2) 69 + 53*107 + O(107^2)               O(107^2)]
    sage: CyclicCover(3, x^3 + x).frobenius_matrix()
    [          0           0      O(107)      O(107)]
    [          0           0 52 + O(107)      O(107)]
    [     O(107) 35 + O(107)           0           0]
    [44 + O(107)      O(107)           0           0]
    sage: CyclicCover(3, 3*x^3 + x).frobenius_matrix()
    [          0           0      O(107)      O(107)]
    [          0           0 79 + O(107)      O(107)]
    [     O(107) 42 + O(107)           0           0]
    [30 + O(107)      O(107)           0           0]


"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2018   Vishal Arul <varul@mit.edu>,
#                       Alex Best <alex.j.best@gmail.com>,
#                       Edgar Costa <edgarcosta@math.dartmouth.edu>,
#                       Richard Magner <rmagner@bu.edu>,
#                       Nicholas Triantafillou <ngtriant@mit.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.arith.misc import euler_phi
from sage.functions.other import ceil, binomial, floor
from sage.functions.log import log
from sage.rings.all import PolynomialRing, PowerSeriesRing
from sage.rings.padics.factory import Zp, Zq, Qq
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.matrix.constructor import matrix, zero_matrix
from sage.modules.free_module_element import vector
from sage.schemes.hyperelliptic_curves.hypellfrob import interval_products
from sage.misc.cachefunc import cached_method

from .charpoly_frobenius import charpoly_frobenius
from . import cycliccover_generic

class CyclicCover_finite_field(cycliccover_generic.CyclicCover_generic):
    def __init__(self, AA, r, f, names=None, verbose = 0):
        cycliccover_generic.CyclicCover_generic.__init__(self, AA, r, f, names=names);
        self._verbose = verbose;
        self._init_frobQ = False;
        self._N0 = None;


    def _init_frob(self, desired_prec = None):
       if not self._init_frobQ or self._N0 != desired_prec:
            if self._r < 2 or self._d < 2:
                raise NotImplementedError("Only implemented for r, f.degree() >= 2");

            self._init_frobQ = True;

            self._Fq = self._f.base_ring()
            self._p = self._Fq.characteristic();
            self._q = self._Fq.cardinality();
            self._n = self._Fq.vector_space().dimension();
            self._epsilon = 0 if self._delta == 1 else 1;


            # our basis choice doesn't always give an integral matrix
            if self._epsilon == 0:
                self._extraprec = floor( log(self._r, self._p) + log( (2*self._genus + (self._delta - 2))/self._delta, self._p));
            else:
                self._extraprec = floor( log(self._r * 2 - 1,  self._p));

            self._nodenominators = self._extraprec == 0;

            if desired_prec is None:
                self._N0 = self._find_N0();
            else:
                self._N0 = desired_prec;

            self._plarge = self._p > self._d * self._r * (self._N0 + self._epsilon);

            # working prec
            if self._plarge:
                self._N = self._N0 + 1;
            else:
                self._N = self._find_N_43();

            # we will use the sqrt(p) version?
            self._sqrtp = self._plarge and self._p == self._q
            self._extraworkingprec = self._extraprec;
            if not self._plarge:
                # we might have some denominators showing up during horizontal and vertical reductions
                self._extraworkingprec =+ 2*ceil(log(self._d * self._r * (self._N0 + self._epsilon), self._p))


            # Rings
            if self._plarge and self._nodenominators:
                if self._n == 1:
                    #IntegerModRing is significantly faster than Zq
                    self._Zq = IntegerModRing(self._p ** (self._N));
                    if self._sqrtp:
                        self._Zq0 = IntegerModRing(self._p ** (self._N  - 1))
                    self._Qq = Qq(self._p, prec=self._N, type = 'capped-rel');
                    self._w = 1;
                else:
                    self._Zq = Zq(self._q, names='w', modulus=self._Fq.polynomial(),
                            prec = self._N, type = "capped-abs");
                    self._w = self._Zq.gen()
                    self._Qq = self._Zq.fraction_field()
            else:
                self._Zq = Qq(self._q, names='w', modulus=self._Fq.polynomial(),
                        prec = self._N + self._extraworkingprec);
                self._w = self._Zq.gen()
                self._Qq = self._Zq;
            self._Zp = Zp(self._p, prec = self._N + self._extraworkingprec);


            self._Zqx = PolynomialRing(self._Zq, "x");

            # Want to take a lift of f from Fq to Zq
            if self._n == 1: # When n = 1, can lift from Fp[x] to Z[x] and then to Zp[x]
                self._flift = self._Zqx( [elt.lift() for elt in  self._f.list()] ) ;
                self._frobf = self._Zqx( self._flift.list() )
            else: # When n > 1, need to be more careful with the lift
                self._flift = self._Zqx( [elt.polynomial().change_ring(ZZ)(self._Zq.gen()) for elt in self._f.list()]);

                self._frobf = self._Zqx( [elt.frobenius() for elt in self._flift.list()] )

            self._dflift = self._flift.derivative()

            # Set up local cache for Frob(f)^s
            frobpow = [None] * (self._N0  + 2) # This variable will store the powers of frob(f)
            frobpow[0] = self._Zqx(1)
            for k in range(0, self._N0 + 1):
                frobpow[k + 1] = self._frobf * frobpow[k]
            # We don't make it a polynomials as we need to keep track that the ith coefficient represents  (i*p)-th
            self._frobpow_list = [ elt.list() for elt in frobpow ]


            if self._sqrtp:
                ## precision of self._Zq0
                N = self._N - 1
                vandermonde = matrix(self._Zq0, N, N );
                for i in range(N):
                    vandermonde[i, 0] = 1;
                    for j in range(1, N):
                        vandermonde[i, j] = vandermonde[i, j-1] * (i + 1);
                self._vandermonde = vandermonde.inverse();

                self._horizontal_fat_s = {};
                self._vertical_fat_s = {};



    def _divide_vector(self, D, vect, R ):
        DQq = self._Qq(D).lift_to_precision( self._Qq.precision_cap() );
        m = 1/DQq;
        if not R.is_field():
            vectQq = vector(self._Qq, [m * self._Qq(elt).lift_to_precision( self._Qq.precision_cap()) for elt in vect]);
            return vector(R, [R(elt) for elt in vectQq]);
        else:
            return vector(R, [(m * elt).lift_to_precision() for elt in vect]);

    def _extend_frobpow(self, power):
        if power < len(self._frobpow_list):
            pass;
        else:
            frobpow = self._Zqx( self._frobpow_list[-1] );
            for k in range(len(self._frobpow_list), power + 1):
                frobpow *= self._frobf;
                self._frobpow_list.extend([ frobpow.list() ]);

        assert power < len(self._frobpow_list);



    def _N0_nodenominators(self):
        return max([ ceil( log( 2 * (2 * self._genus)/i, self._p) + (self._n * i) / ZZ(2)) for i in range(1, self._genus +1 )]);
    def _N0_RH(self):
        return ceil( log( 2 * binomial(2*self._genus, self._genus), self._p) + self._genus * self._n/ZZ(2));

    def _find_N0(self):
        if self._nodenominators:
            return self._N0_nodenominators();
        else:
            return self._N0_RH() + self._extraprec;

    def _find_N(self):
        # for p >> 0
        # N = N0 + 1
        # Theorem 7.3
        p = self._p;
        r = self._r;
        N0 = self._N0;
        def left_side_log(n):
            return floor(log(p*r*(n+1) - 3*r)/log(p));
        n = N0 + 1;
        while n  < N0 + left_side_log(n):
            n += 1
        return n


    def _find_N_43(self):
        # Goncalves - Theorem 4.3
        # for p >> 0, N = N0 + 2
        p = self._p;
        r = self._r;
        d = self._d;
        delta = self._delta;
        N0 = self._N0;
        n = N0;
        left_side = N0 +  floor(log((d * p * (r - 1) + r)/delta)/log(p))
        def right_side_log(n):
            return floor(log(p*(r*n - 1) -r)/log(p))
        n = left_side
        while n <= left_side + right_side_log(n):
            n += 1
        return n

    def _frob_sparse(self, i, j, N0):
        r"""
        Computes `Frob(x^i y^(-j) dx ) / dx` for y^r = f(x) with N0 terms

        INPUT:

        -   ``i`` - The power of x in the expression `Frob(x^i dx/y^j) / dx`

        -   ``j`` - The (negative) power of y in the expression `Frob(x^i dx/y^j) /
            dx`

        OUTPUT:

        ``frobij`` - a Matrix of size  (d * (N0 - 1) + ) x (N0)
                     that represents the Frobenius expansion of x^i dx/y^j modulo p^(N0 + 1)

                    the entry (l, s) corresponds to the coefficient associated to the monomial x ** (p * (i + 1 + l) -1) * y ** (p * -(j + r*s))
                    (l, s) --> (p * (i + 1 + l) -1, p * -(j + r*s))

        ALGORITHM:

        Computes::

        Frob(x^i dx/y^j) / dx = p * x ** (p * (i+1) - 1) * y ** (-j*p) * Sigma

        where::

        Sigma = \sum_{k = 0} ^{N0-1}
                        \sum_{s = 0} ^k
                            (-1) ** (k-s) * binomial(k, s)
                            * binomial(-j/r, k)
                            * self._frobpow[s]
                            * self._y ** (-self._r * self._p * s)
                = \sum_{s = 0} ^{N0 - 1}
                    \sum_{k = s} ^N0
                        (-1) ** (k-s) * binomial(k, s)
                        * binomial(-j/self._r, k)
                        * self._frobpow[s]
                        * self._y ** (-self._r*self._p*s)
                = \sum_{s = 0} ^{N0-1}
                        D_{j, s}
                        * self._frobpow[s]
                        * self._y ** (-self._r * self._p * s)
                = \sum_{s = 0} ^N0
                    \sum_{l = 0} ^(d*s)
                        D_{j, s} * self._frobpow[s][l]
                        * x ** (self._p ** l)
                        * y ** (-self._r * self._p ** s)

        and::

        D_{j, s} = \sum_{k = s} ^N0 (-1) ** (k-s) * binomial(k, s) * binomial(-j/self._r, k) )


        """
        r = self._r;
        Dj = [ self._Zq(sum( [(-1) ** (k-l) * binomial(k, l) * binomial(-ZZ(j)/r, k)  for k in range(l, N0)] )) for l in range(N0)];
        self._extend_frobpow(N0);
        frobij = matrix(self._Zq,  self._d * (N0 - 1) + 1, N0);
        for s in range(N0):
            for l in range(self._d * s + 1):
                frobij[ l,  s  ] = self._p * Dj[s] * self._frobpow_list[s][l];
        return frobij

    def _horizonal_matrix_reduction(self, s):
        r"""

        Return the tuple of tuples that represents the horizonal matrix reduction at poler order ``s``.

        INPUT:
            - ``s`` -- the integer s

        OUTPUT:

            A tuple of tuples ``( (D0, D1), (M0, M1) )``
            where MH_{e, s}  = M0 + e * M1 and DH_{e,s} = D0 + e * D1

        ALGORITHM:

        Let W_{e, s} to be the Qq-vector space of differential forms of the form
            G x^e y^{-s} dx
        where \deg G \leq d - 1.

        Let v = [G_0, ..., G_{d-1}] represent G

        There is a map
            MH_{e, s} : W_{e, s} --> W_{e-1, s}
        and a function to
            DH: \NN x \NN -->  Qq
        such that

        G x^e y^{-s} dx \cong H x^{e - 1} y^{-s} dx

        where H = DH(e, s)^{-1} * MH_{e,s} ( G )

        The matrix MH_{e, s} can be written as
             MH_{e, s}  = M0_{s} + e * M1_{s}
        similarly
            DH_{e,s} = D0_{s} + e * D1_{s}
        """

        f_co = self._flift.list()
        # DH_{e,s} = ((s -r)*d - e * r) * f_d
        m1 = -1 * self._r * f_co[-1] # r is a, g_co[-1] is lambda
        m0 = (s - self._r) * self._d * f_co[-1] # j is a*t + beta

        M1 = matrix(self._Zq, self._d, lambda m,n: m1 if m == n + 1 else self._r * f_co[m] if n == self._d - 1 else 0)
        M0 = matrix(self._Zq, self._d, lambda m,n: m0 if m == n + 1 else (self._r - s) * m * f_co[m] if n == self._d - 1 else 0)

        return ((m0, m1), (M0, M1))

    def _vertical_matrix_reduction(self, s_0):
        r"""
        Return the tuple of tuples that represents the vertical matrix reduction.

        OUTPUT:

            A tuple of tuples ``( (D0, D1), (M0, M1) )``
            where MV_t  = M0 + t * M1 and DV_t = D0 + t * D1

        ALGORITHM:

        """

        d = self._d;
        f_co = [0 for i in range(d-2)] + self._flift.list() + [0 for i in range(d-1)]
        fd_co = [0 for i in range(d-1)] + self._dflift.list() + [0 for i in range(d-0)]

        rows = [f_co[d-2-i:-i-1] for i in range(d-1)]
        rows += [fd_co[d-1-i:-i-1] for i in range(d)]

        m = matrix(rows).transpose().inverse()

        a_foo = m[0:d,0:d]
        b_foo = m[d-1:2*d-1,0:d]
        a_foo = matrix(d, d, lambda i, j: 1 if i == j and i != d-1 else 0) * a_foo
        foo = matrix(d, d , lambda i, j: j if i == j-1 else 0)
        bp_foo = foo * b_foo
        A_vert = a_foo.submatrix(0,0,d-1,d-1);
        Bd_vert = bp_foo.submatrix(0,0,d-1,d-1);
        M1 = (s_0-self._r) * A_vert + self._r * Bd_vert
        M2 = self._r * A_vert
        m1 = (s_0-self._r)
        m2 = self._r
        return ( (m1, m2), (M1, M2) )

    def _reduce_vector_horizontal(self, G, e, s, k = 1, method = 0):
        if self._sqrtp and k == self._p:
            vect =  self._reduce_vector_horizontal_BSGS(G, e, s);
        else:
            vect =  self._reduce_vector_horizontal_plain(G, e, s, k);
        return vect

    def _reduce_vector_horizontal_BSGS(self, G, e, s):
        r"""
        Input:
            - a vector -- G \in W_{e, s}
        Output:
            - a vector -- H \in W_{e - p, s} such that
                G x^e y^{-s} dx \cong H x^{e - p} y^{-s} dx
        """
        if G == 0:
            return G
        if self._verbose > 2:
            print("_reduce_vector_horizontal_BSGS(self, %s, %s, %s)" % (vector(self._Qq,G), e, s,));
        assert (e + 1) % self._p == 0
        (m0, m1), (M0, M1) = self._horizonal_matrix_reduction(s);
        vect = vector(self._Zq, G);
        # we do the first d reductions carefully
        D = 1;
        for i in reversed(range(e - self._d + 1, e + 1)):
            Mi = M0 + i*M1
            Di = m0 + i*m1

            vect = Mi*vect
            D *= Di
        assert Di % self._p == 0
        iD = 1/self._Zq0(D.lift()/self._p);
        vect = vector(self._Zq0, [iD * ZZ(elt.lift()/self._p) for elt in vect])
        #use BSGS
        iDH, MH = self._get_fat_horizontal_matrix(e, s);
        vect = iDH* (MH * vect.change_ring(self._Zq0));

        # last reduction
        i = e - self._p + 1;
        Mi = M0 + i*M1
        Di = 1/( m0 + i*m1);

        vect = Di * (Mi * vect.change_ring(self._Zq));

        if self._verbose > 2:
            print("done _reduce_vector_horizontal_BSGS(self, %s, %s, %s)" % (vector(self._Qq,G), e, s,));
            print("return %s\n" % (vector(self._Qq, vect),));
        return vect;

    def _initialize_fat_horizontal(self, s, L):
        assert self._sqrtp
        if s not in self._horizontal_fat_s:
            N = self._N - 1; # padic precision of self._Zq0
            d = self._d
            L0 = min(L, N);
            targets = [0] * (2*L0);
            for l in range(L0):
                targets[2*l] = self._p*l;
                targets[2*l + 1] = self._p*(l+1) - d - 1;
            (m0, m1), (M0, M1) = self._horizonal_matrix_reduction(s);
            M0, M1 = [ elt.change_ring(self._Zq0) for elt in [M0, M1] ]
            D0, D1 = [ matrix(self._Zq0, [elt]) for elt in [m0, m1]];
            MH = interval_products(M0, M1, targets);
            DH = [elt[0,0] for elt in interval_products(D0, D1, targets)];
            if L > N:
                #  f^{(r)}(0) p^r / r! for r = 0, ..., N-1,
                MT = [None] * N;
                DT = [0] * N;
                for r in range(N):
                    MT[r] = matrix(self._Zq0, d, d);
                    for h in range(N):
                        v = self._vandermonde[r, h];
                        for i in range(d):
                            for j in range(d):
                                MT[r][i, j] += v*MH[h][i, j];
                        DT[r] += v * DH[h];
                for k in range(N, L):
                    M = matrix(self._Zq0, d, d);
                    D = 0;
                    k1pow = self._Zq0(1); # power of k + 1
                    for h in range(N):
                        for i in range(d):
                            for j in range(d):
                                M[i,j] += k1pow * MT[h][i,j];
                        D += k1pow * DT[h];
                        k1pow *= (k + 1);
                    MH.append(M);
                    DH.append(D);

            iDH = [ 1/elt for elt in DH ];
            self._horizontal_fat_s[s] = [ (iDH[i], MH[i]) for i in range(0, L) ];
        else:
            assert len(self._horizontal_fat_s[s]) <= L





    def _get_fat_horizontal_matrix(self, e, s):
        return self._horizontal_fat_s[s][ (e + 1)/self._p - 1];



    def _reduce_vector_horizontal_plain(self, G, e, s, k = 1):
        r"""
        Input:
            - a vector -- G \in W_{e, s}
        Output:
            - a vector -- H \in W_{e - k, s} such that
                G x^e y^{-s} dx \cong H x^{e - k} y^{-s} dx
        """
        if self._verbose > 2:
            print("_reduce_vector_horizontal_plain(self, G = %s, e = %s, s = %s, k = %s)" % (vector(self._Qq,G), e, s, k,))
        if G == 0:
            return G
        (m0, m1), (M0, M1) = self._horizonal_matrix_reduction(s);
        vect = vector(self._Zq, G);
        D = self._Zq(1);
        assert k <= self._p, "more than p reductions at a time should be avoided!"
        assert e - k + 1 >= 0
        for i in reversed(range(e - k + 1, e + 1)):
            Mi = M0 + i*M1
            Di = m0 + i*m1

            vect = Mi*vect
            D *= Di;
            if self._plarge and Di % self._p == 0:
                assert (i  + self._d) % self._p == 0
                vect = self._divide_vector( D, vect, self._Zq);
                D = self._Zq(1);

        vect = self._divide_vector( D, vect, self._Zq);

        if self._verbose > 2:
            print("done _reduce_vector_horizontal_plain(self, %s, %s, %s, %s)" % (vector(self._Qq,G), e, s, k,));
            print("return %s\n" % (vector(self._Qq, vect),));
        return vect


    def _reduce_vector_vertical(self, G, s0, s, k = 1,method = 0):
        if self._sqrtp:
            return self._reduce_vector_vertical_BSGS(G, s0, s, k)
        else:
            return self._reduce_vector_vertical_plain(G, s0, s, k);

    def _reduce_vector_vertical_BSGS(self, G, s0, s, k):
        MV = self._get_fat_vertical_matrix(s0, s, k);
        return MV * G;

    def _initialize_fat_vertical(self, s0, max_upper_target ):
        L = floor((max_upper_target - self._epsilon)/self._p) + 1
        if s0 not in self._vertical_fat_s:
            (m0, m1), (M0, M1) = self._vertical_matrix_reduction(s0);
            D0, D1 = map( lambda y: matrix(self._Zq, [y]), [m0, m1]);
            targets = [0] * (2*L);
            for l in reversed(range(L)):
                targets[2*l] = max_upper_target - self._p*(L - l);
                targets[2*l + 1] = max_upper_target -  self._p*(L - 1 - l);
            if targets[0] < 0:
                targets[0] = self._epsilon
            MV = interval_products(M0, M1, targets)
            DV = interval_products(D0, D1, targets)
            for l in range(L):
                D = DV[l][0,0]
                if D % self._p == 0:
                    iD = 1/self._Zq(D.lift()/self._p)
                    MV[l] = matrix(self._Zq, [[ iD * ZZ(elt.lift()/self._p) for elt in row] for row in MV[l].rows()]);
                else:
                    MV[l] *= (1/D)
            self._vertical_fat_s[s0] = MV;
        else:
            assert len(self._vertical_fat_s[s0]) <= L

    def _get_fat_vertical_matrix(self, s0, s, k):
        if s0 in self._vertical_fat_s:
            if k < self._p:
                assert s - k == self._epsilon
                return self._vertical_fat_s[s0][0]
            elif k == self._p:
                return self._vertical_fat_s[s0][s//self._p]

        (m0, m1), (M0, M1) = self._vertical_matrix_reduction(s0);
        target =  [ s - k, s]
        MV = interval_products(M0, M1, target)[0]
        D0, D1 = map( lambda y: matrix(self._Zq, [y]), [m0, m1]);
        DV = interval_products(D0, D1, target)[0][0,0]
        if DV % self._p == 0:
            # we need to do the reduction here, as both terms might be divisible by p
            iDV = 1/self._Zq(DV.lift()/self._p)
            MV = matrix(self._Zq, [[ iDV * ZZ(elt.lift()/self._p) for elt in row] for row in MV.rows()]);
            return MV;
        else:
            return MV/DV


    def _reduce_vector_vertical_plain(self, G, s0, s, k = 1):
        r"""
        Input:
            - a vector -- G \in W_{-1, r*s + s0}
        Output:
            - a vector -- H \in W_{-1, r*(s - k) + s0} such that 
                G y^{-(r*s + s0)} dx \cong H y^{-(r*(s -k) + s0)} dx
        """
        if self._verbose > 2:
            print("_reduce_vector_vertical(self, G = %s, s0 = %s, s = %s, k = %s)" % (vector(self._Qq, G), s0,  s, k,))

        (m0, m1), (M0, M1) = self._vertical_matrix_reduction(s0);
        vect = vector(self._Zq, G);
        D = 1;
        assert k <= self._p, "more than p reductions at a time should be avoided!"
        assert s - k >= 0
        for i in reversed(range(s - k + 1,s + 1)):
            Mi = M0 + i*M1
            Di = m0 + i*m1
            vect = Mi*vect
            if self._plarge and Di % self._p != 0:
                D *= Di;
            else:
                vect = self._divide_vector( Di, vect, self._Zq);

        vect = self._divide_vector( D, vect, self._Zq);

        if self._verbose > 2:
            print("done _reduce_vector_vertical(self,  %s, %s, %s)" % (vector(self._Qq, G), s, k,))
            print("return %s\n" % (vector(self._Qq, vect),))

        return vect;

    def _frob(self, i, j, N0):
        r"""
        Compute Frob(x^i dx/y^j) / dx in terms of the cohomology basis,
        whose x and y exponents are constrained to be in a particular range.
        """
        # a Matrix that represents the Frobenius expansion of
        # x^i dx/y^j modulo p^(N0 + 1)
        # the entry (l, s) corresponds to the coefficient associated to the monomial x ** (p * (i + 1 + l) -1) * y ** (p * (j + r*s))
        assert N0 <= self._N0;
        frobij = self._frob_sparse(i, j, N0);
        # recall
        # Frob(x^i dx/y^j) / dx
        #       = p * x ** (p * (i+1) - 1) * y ** (j*p)
        #         * = \sum_{s = 0} ^{N0-1}
        #               \sum_{l = 0} ^(d*s)
        #                   D_{j, s} * Frobpow[s][l] * x ** (p ** l) y ** (r * p ** s)
        # H represents H(x) * y^(-p**s) s /dx
        # the entry (l, s) of frobij
        # corresponds to the monomial (p * (i + 1 + l) -1, p * -(j + r*s))
        d = self._d
        r = self._r
        p = self._p
        H = vector(self._Zq, d - 1);
        k = (p * j) // r;
        s0 = (p * j) % r

        if self._sqrtp:
            self._initialize_fat_vertical(s0,  k  +  p*(N0 - 1));

        for s in reversed(range(0, N0)):
            if self._sqrtp:
                # (i + 1) <= d
                self._initialize_fat_horizontal(p*j + p*r*s, d* s + (d-2) + 1); #d * (s + 1) )
            # G represents G(x) * x^(p ** l - 1) y^(-p(j + r*s)) /dx
            G = vector(self._Zq, d);
            for l in reversed(range(1, d * s + (i + 1) + 1)):
                if l >= (i + 1):
                    G[0] += frobij[l - (i+1), s];
                G = self._reduce_vector_horizontal(G, p * l - 1 , p*j + p*r*s, p);
            assert G[0] == 0;
            H += vector(G.list()[1:]);
            if s > 0:
                # y^(-p(j + r*s))  --- >  y^(-p(j + r*(s-1)))
                H = self._reduce_vector_vertical(H, s0, k + p*s, p);
        # H represents
        # H(x) y^-p*j = H(x) y^-(k*r + s0)
        # now we reduce the pole order to s0 + r*epsilon, where s0 = p *j % r
        while k > self._epsilon:
            steps = p if k - self._epsilon > p else k - self._epsilon;
            H = self._reduce_vector_vertical(H, s0, k, steps);
            k -= steps;
        assert k == self._epsilon
        H = [self._Qq(elt).add_bigoh(N0) for elt in H]
        return  matrix(H).transpose()

    def _frobenius_matrix_p(self, N0):
        r"""
        computes the matrix that represents the p-power Frobenius
        """
        assert self._init_frobQ;
        assert N0 <= self._N0;

        m = matrix(self._Qq, (self._d - 1)*(self._r - 1))

        # Want to build m, "slice by slice" using the output of _frob
        for j in range(1,self._r):
            s0 = (j*self._p) % self._r
            for i in range(0,self._d-1):
                m[(s0 - 1)*(self._d - 1):s0*(self._d - 1),i+(j-1)*(self._d-1)] = self._frob(i, j + self._epsilon*self._r, N0 )
        return m

    def frobenius_matrix(self, N = None):
        """
            Compute p-adic frobenius matrix to precision p^N. If N not
            supplied, a default value is selected, which is the minimum needed
            to recover the charpoly unambiguously.

        EXAMPLES::

            sage: p = 107;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: CyclicCover(2, x^5 + x).frobenius_matrix()
            [              O(107^2)      89*107 + O(107^2)               O(107^2)               O(107^2)]
            [     89*107 + O(107^2)               O(107^2)               O(107^2)               O(107^2)]
            [              O(107^2)               O(107^2)               O(107^2) 105 + 5*107 + O(107^2)]
            [              O(107^2)               O(107^2) 89 + 53*107 + O(107^2)               O(107^2)]
            sage: CyclicCover(2, 3*x^5 + x).frobenius_matrix()
            [              O(107^2)      14*107 + O(107^2)               O(107^2)               O(107^2)]
            [     69*107 + O(107^2)               O(107^2)               O(107^2)               O(107^2)]
            [              O(107^2)               O(107^2)               O(107^2) 61 + 58*107 + O(107^2)]
            [              O(107^2)               O(107^2) 69 + 53*107 + O(107^2)               O(107^2)]
            sage: CyclicCover(3, x^3 + x).frobenius_matrix()
            [          0           0      O(107)      O(107)]
            [          0           0 52 + O(107)      O(107)]
            [     O(107) 35 + O(107)           0           0]
            [44 + O(107)      O(107)           0           0]
            sage: CyclicCover(3, 3*x^3 + x).frobenius_matrix()
            [          0           0      O(107)      O(107)]
            [          0           0 79 + O(107)      O(107)]
            [     O(107) 42 + O(107)           0           0]
            [30 + O(107)      O(107)           0           0]
        """
        self._init_frob(N);
        FrobP = self._frobenius_matrix_p(self._N0);
        if self._n == 1:
            return FrobP;
        else:
            current = FrobP
            total = FrobP
            for i in range(0, self._n - 1):
                current = matrix([[entry.frobenius() for entry in row] for row in current])
                total = total * current
            return total;

    def _denominator(self):
        R = PolynomialRing(ZZ, "T");
        T = R.gen();
        denom = R(1)
        lc = self._f.list()[-1];
        if lc == 1: # MONIC
            for i in range(2, self._delta + 1):
                if self._delta % i == 0:
                    phi = euler_phi(i)
                    G = IntegerModRing(i);
                    ki = G(self._q).multiplicative_order()
                    denom = denom * (T**ki - 1)**(phi//ki)
            return denom
        else: # NON-MONIC
            x = PolynomialRing(self._Fq, "x").gen();
            f = x**self._delta - lc;
            L = f.splitting_field("a")
            roots = [r for r,_ in f.change_ring(L).roots()]
            roots_dict = dict([(r, i) for i,r in  enumerate(roots)])
            rootsfrob = [L.frobenius_endomorphism()(r) for r in roots]
            m = zero_matrix(len(roots))
            for i, r in enumerate(roots):
                m[i, roots_dict[rootsfrob[i]]] = 1
            return R(R(m.characteristic_polynomial())//(T - 1))


    @cached_method
    def frobenius_polynomial(self):
        r"""
        Returns the characteristic polynomial of Frobenius.

        EXAMPLES:


        Hyperelliptic curves::

            sage: p = 11;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: f = x^9 + 4*x^2 + 10*x + 4;
            sage: CyclicCover(2, f).frobenius_polynomial() == \
            ....: HyperellipticCurve(f).frobenius_polynomial()
            True
            sage: f = 2*x^5 + 4*x^3 + x^2 + 2*x + 1;
            sage: CyclicCover(2, f).frobenius_polynomial() == \
            ....: HyperellipticCurve(f).frobenius_polynomial()
            True
            sage: f = 2*x^6 + 4*x^4 + x^3 + 2*x^2 + x;
            sage: CyclicCover(2, f).frobenius_polynomial() == \
            ....: HyperellipticCurve(f).frobenius_polynomial()
            True
            sage: p = 1117;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: f = x^9 + 4*x^2 + 10*x + 4;
            sage: CyclicCover(2, f).frobenius_polynomial() == \
            ....: HyperellipticCurve(f).frobenius_polynomial()
            True
            sage: f = 2*x^5 + 4*x^3 + x^2 + 2*x + 1;
            sage: CyclicCover(2, f).frobenius_polynomial() == \
            ....: HyperellipticCurve(f).frobenius_polynomial()
            True

        Super Elliptic curves::

            sage: p = 11;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: CyclicCover(3, x^4 + 4*x^3 + 9*x^2 + 3*x + 1).frobenius_polynomial()
            x^6 + 21*x^4 + 231*x^2 + 1331
            sage: CyclicCover(4,x^5 + x + 1).frobenius_polynomial()
            x^12 - 4*x^11 + 38*x^10 - 140*x^9 + 743*x^8 - 2200*x^7 + 9812*x^6 - 24200*x^5 + 89903*x^4 - 186340*x^3 + 556358*x^2 - 644204*x + 1771561
            sage: p = 4999;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: CyclicCover(3, x^4 + 4*x^3 + 9*x^2 + 3*x + 1).frobenius_polynomial()
            x^6 + 180*x^5 + 20988*x^4 + 1854349*x^3 + 104919012*x^2 + 4498200180*x + 124925014999
            sage: CyclicCover(4,x^5 + x + 1).frobenius_polynomial()
            x^12 - 64*x^11 + 5018*x^10 - 488640*x^9 + 28119583*x^8 - 641791616*x^7 + 124245485932*x^6 - 3208316288384*x^5 + 702708407289583*x^4 - 61043359329111360*x^3 + 3133741752599645018*x^2 - 199800079984001599936*x + 15606259372500374970001


        Non Super Elliptic curves:

            sage: p = 5;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: C = CyclicCover(6, x^6 + 1);
            sage: C.frobenius_polynomial()
            x^20 + 50*x^18 + 1125*x^16 + 15000*x^14 + 131250*x^12 + 787500*x^10 + 3281250*x^8 + 9375000*x^6 + 17578125*x^4 + 19531250*x^2 + 9765625
            sage: R.<t> = PowerSeriesRing(Integers());
            sage: C.projective_closure().zeta_series(4,t)
            1 + 6*t + 81*t^2 + 456*t^3 + 3456*t^4 + O(t^5)
            sage: C.frobenius_polynomial().reverse()(t)/((1-t)*(1-p*t)) + O(t^5)
            1 + 6*t + 81*t^2 + 456*t^3 + 3456*t^4 + O(t^5)

            sage: x = PolynomialRing(GF(11),"x").gen();
            sage: CyclicCover(4, x^6 - 11*x^3 + 70*x^2 - x + 961).frobenius_polynomial()
            x^14 + 14*x^12 + 287*x^10 + 3025*x^8 + 33275*x^6 + 381997*x^4 + 2254714*x^2 + 19487171
            sage: x = PolynomialRing(GF(4999),"x").gen();
            sage: CyclicCover(4, x^6 - 11*x^3 + 70*x^2 - x + 961).frobenius_polynomial()
            x^14 - 4*x^13 - 2822*x^12 - 30032*x^11 + 37164411*x^10 - 152369520*x^9 + 54217349361*x^8 - 1021791160888*x^7 + 271032529455639*x^6 - 3807714457169520*x^5 + 4642764601604000589*x^4 - 18754988504199390032*x^3 - 8809934776794570547178*x^2 - 62425037490001499880004*x + 78015690603129374475034999

            sage: p = 11;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: CyclicCover(3, 5*x^3 - 5*x + 13).frobenius_polynomial()
            x^2 + 11
            sage: CyclicCover(3, x^6 + x^4 - x^3 + 2*x^2 - x - 1).frobenius_polynomial()
            x^8 + 32*x^6 + 462*x^4 + 3872*x^2 + 14641
            sage: p = 4999;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: CyclicCover(3, 5*x^3 - 5*x + 13).frobenius_polynomial()
            x^2 - 47*x + 4999
            sage: CyclicCover(3, x^6 + x^4 - x^3 + 2*x^2 - x - 1).frobenius_polynomial()
            x^8 + 122*x^7 + 4594*x^6 - 639110*x^5 - 82959649*x^4 - 3194910890*x^3 + 114804064594*x^2 + 15240851829878*x + 624500149980001

            sage: p = 11;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: CyclicCover(5, x^5 + x).frobenius_polynomial()
            x^12 + 4*x^11 + 22*x^10 + 108*x^9 + 503*x^8 + 1848*x^7 + 5588*x^6 + 20328*x^5 + 60863*x^4 + 143748*x^3 + 322102*x^2 + 644204*x + 1771561
            sage: CyclicCover(5, 2*x^5 + x).frobenius_polynomial()
            x^12 - 9*x^11 + 42*x^10 - 108*x^9 - 47*x^8 + 1782*x^7 - 8327*x^6 + 19602*x^5 - 5687*x^4 - 143748*x^3 + 614922*x^2 - 1449459*x + 1771561
            sage: p = 49999;
            sage: x = PolynomialRing(GF(p),"x").gen();
            sage: CyclicCover(5, x^5 + x ).frobenius_polynomial()
            x^12 + 299994*x^10 + 37498500015*x^8 + 2499850002999980*x^6 + 93742500224997000015*x^4 + 1874812507499850001499994*x^2 + 15623125093747500037499700001
            sage: CyclicCover(5, 2*x^5 + x).frobenius_polynomial()
            x^12 + 299994*x^10 + 37498500015*x^8 + 2499850002999980*x^6 + 93742500224997000015*x^4 + 1874812507499850001499994*x^2 + 15623125093747500037499700001


        """
        self._init_frob();
        F = self.frobenius_matrix(self._N0);

        denom = self._denominator();
        R = PolynomialRing(ZZ, "x");

        if self._nodenominators:
            min_val = 0;
        else:
            # are there any denominators in F?
            min_val = min( [ min( [self._Qq(elt).valuation() for elt in row] ) for row in F.rows()]);

        if min_val >= 0:
            prec = self._N0_nodenominators()
            charpoly_prec = [prec + i  for i in reversed(range(1,self._genus + 1))] + [prec ]*(self._genus + 1);
            cp = charpoly_frobenius(F, charpoly_prec, self._p, 1, self._n,  denom.list());
            return R(cp);
        else:
            cp = F.charpoly().reverse();
            denom = denom.reverse();
            PS = PowerSeriesRing(self._Zp, "T")
            cp = PS(cp)/PS(denom);
            cp = cp.padded_list(self._genus + 1);
            cpZZ = [None for _ in range(2*self._genus + 1)];
            cpZZ[0] = 1;
            cpZZ[-1] = self._p**self._genus;
            for i in range(1, self._genus + 1):
                cmod = cp[i]
                bound = binomial(2*self._genus, i)*self._p**(i*self._n * 0.5);
                localmod = self._p**( ceil(log(bound, self._p)));
                c = cmod.lift() % localmod;
                if c > bound:
                    c = -(-cmod.lift() % localmod)
                cpZZ[i] = c;
                if i != self._genus + 1:
                    cpZZ[2*self._genus - i] = c*self._p**(self._genus - i)
            cpZZ.reverse()
            return R(cpZZ)
