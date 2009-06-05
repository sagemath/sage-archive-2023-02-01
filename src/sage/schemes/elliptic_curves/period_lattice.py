r"""
Period lattices of elliptic curves and related functions.

Let `E` be an elliptic curve defined over a number field `K`
(including `\QQ`).  We attach a period lattice (a discrete rank 2
subgroup of `\CC`) to each embedding of `K` into `\CC`.

In the case of real embeddings, the lattice is stable under complex
conjugation and is called a real lattice.  These have two types:
rectangular, (the real curve has two connected components and positive
discriminant) or non-rectangular (one connected component, negative
discriminant).

The periods are computed to arbitrary precision using the AGM (Gauss's
Arithmetic-Geometric Mean).

EXAMPLES::

    sage: K.<a> = NumberField(x^3-2)
    sage: E = EllipticCurve([0,1,0,a,a])

First we try a real embedding::

    sage: emb = K.embeddings(RealField())[0]
    sage: L = E.period_lattice(emb); L
    Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a over Number Field in a with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
    From: Number Field in a with defining polynomial x^3 - 2
    To:   Algebraic Real Field
    Defn: a |--> 1.259921049894873?

The first basis period is real::

    sage: L.basis()
    (3.81452977217855, 1.90726488608927 + 1.34047785962440*I)
    sage: L.is_real()
    True

For a basis `\omega_1,\omega_2` normalised so that `\omega_1/\omega_2`
is in the fundamental region of the upper half-plane, use the function
``normalised_basis()`` instead::

    sage: L.normalised_basis()
    (1.90726488608927 - 1.34047785962440*I, -1.90726488608927 - 1.34047785962440*I)

Next a complex embedding::

    sage: emb = K.embeddings(ComplexField())[0]
    sage: L = E.period_lattice(emb); L
    Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a over Number Field in a with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
    From: Number Field in a with defining polynomial x^3 - 2
    To:   Algebraic Field
    Defn: a |--> -0.6299605249474365? - 1.091123635971722?*I

In this case, the basis `\omega_1`, `\omega_2` is always normalised so
that `\tau = \omega_1/\omega_2` is in the fundamental region in the
upper half plane::

    sage: w1,w2 = L.basis(); w1,w2
    (-1.37588604166076 - 2.58560946624443*I, -2.10339907847356 + 0.428378776460622*I)
    sage: L.is_real()
    False
    sage: tau = w1/w2; tau
    0.387694505032876 + 1.30821088214407*I
    sage: L.normalised_basis()
    (-1.37588604166076 - 2.58560946624443*I, -2.10339907847356 + 0.428378776460622*I)

AUTHORS:

- ?: initial version.

- John Cremona:

  - Adapted to handle real embeddings of number fields, September 2008.

  - Added basis_matrix function, November 2008

  - Added support for complex embeddings, May 2009.

"""

from sage.modules.free_module import FreeModule_generic_pid
from sage.rings.all import ZZ, QQ, RealField, ComplexField, is_RealField, PolynomialRing, QQbar, AA
from sage.rings.real_mpfr import RealNumber as RealNumber
from sage.rings.complex_number import ComplexNumber as ComplexNumber
from sage.rings.number_field.number_field import refine_embedding
from sage.rings.infinity import Infinity
from sage.functions.other import imag

class PeriodLattice(FreeModule_generic_pid):
    """
    The class for the period lattice of an algebraic variety.
    """
    pass

class PeriodLattice_ell(PeriodLattice):
    r"""
    The class for the period lattice of an elliptic curve.

    Currently supported are elliptic curves defined over `\QQ`, and
    elliptic curves defined over a number field with a real or complex
    embedding, where the lattice constructed depends on that
    embedding.
    """

    def __init__(self, E, embedding=None):
        r"""
        Initialises the period lattice by storing the elliptic curve and the embedding.

        INPUT:

        - ``E`` -- an elliptic curve

        - ``embedding`` (defult: ``None``) -- an embedding of the base
          field `K` of ``E`` into a real or complex field.  If
          ``None``:

          - use the built-in coercion to `\RR` for `K=\QQ`;

          - use the first embedding into `\RR` given by
          ``K.embeddings(RealField())``, if there are any;

          - use the first embedding into `\CC` given by
          ``K.embeddings(ComplexField())``, if `K` is totally complex.

        .. note::

           No periods are computed on creation of the lattice; see the
           functions ``basis()``, ``normalised_basis()`` and
           ``real_period()`` for precision setting.

        EXAMPLES:

        This function is not normally called directly, but will be
        called by the period_lattice() function of classes
        ell_number_field and ell_rational_field::

            sage: from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
            sage: E = EllipticCurve('37a')
            sage: PeriodLattice_ell(E)
            Period lattice associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = PeriodLattice_ell(E,emb); L
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a over Number Field in a with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
            From: Number Field in a with defining polynomial x^3 - 2
            To:   Algebraic Real Field
            Defn: a |--> 1.259921049894873?

            sage: emb = K.embeddings(ComplexField())[0]
            sage: L = PeriodLattice_ell(E,emb); L
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a over Number Field in a with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
            From: Number Field in a with defining polynomial x^3 - 2
            To:   Algebraic Field
            Defn: a |--> -0.6299605249474365? - 1.091123635971722?*I

        TESTS::

            sage: from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = PeriodLattice_ell(E,emb)
            sage: L == loads(dumps(L))
            True
        """
        # First we cache the elliptic curve with this period lattice:

        self.E = E

        # Next we cache the embedding into QQbar or AA which extends
        # the given embedding:

        K = E.base_field()
        if embedding is None:
            embs = K.embeddings(AA)
            real = len(embs)>0
            if not real:
                embs = K.embeddings(QQbar)
            embedding = embs[0]
        else:
            embedding = refine_embedding(embedding,Infinity)
            real = embedding(K.gen()).imag().is_zero()

        self.embedding = embedding

        # Next we compute and cache (in self.real_flag) the type of
        # the lattice: +1 for real rectangular, -1 for real
        # non-rectangular, 0 for non-real:

        self.real_flag = 0
        if real:
            self.real_flag = +1
            if embedding(E.discriminant())<0:
                self.real_flag = -1

        # The following algebraic data associated to E and the
        # embedding is cached:
        #
        # Ebar: the curve E base-changed to QQbar (or AA)
        # f2: the 2-division polynomial of Ebar
        # ei: the roots e1, e2, e3 of f2, as elements of QQbar (or AA)
        #
        # The ei are used both for period computation and elliptic
        # logarithms.

        self.Ebar = self.E.change_ring(self.embedding)
        self.f2 = self.Ebar.two_division_polynomial()
        if self.real_flag == 1: # positive discriminant
            self._ei = self.f2.roots(AA,multiplicities=False)
            self._ei.sort()  # e1 < e2 < e3
            e1, e2, e3 = self._ei
        elif self.real_flag == -1: # negative discriminant
            self._ei = self.f2.roots(QQbar,multiplicities=False)
            self._ei = list(sorted(self._ei,key=lambda z: z.imag()))
            e1, e3, e2 = self._ei # so e3 is real
            e3 = AA(e3)
            self._ei = [e1, e2, e3]
        else:
            self._ei = self.f2.roots(QQbar,multiplicities=False)
            e1, e2, e3 = self._ei

        # The quantities sqrt(e_i-e_j) are cached (as elements of
        # QQbar) to be used in period computations:

        self._abc = (e3-e1).sqrt(), (e3-e2).sqrt(), (e2-e1).sqrt()

        PeriodLattice.__init__(self, base_ring=ZZ, rank=2, degree=1, sparse=False)

    def __cmp__(self, other):
        r"""
        Comparison function for period lattices

        TESTS::

            sage: from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
            sage: K.<a> = NumberField(x^3-2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: embs = K.embeddings(ComplexField())
            sage: L1,L2,L3 = [PeriodLattice_ell(E,e) for e in embs]
            sage: L1 < L2 < L3
            True

        """
        if not isinstance(other, PeriodLattice_ell): return -1
        t = cmp(self.E, other.E)
        if t: return t
        a = self.E.base_field().gen()
        return cmp(self.embedding(a), other.embedding(a))

    def __repr__(self):
        """
        Returns the string representation of this period lattice.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb); L
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 + a*x + a over Number Field in a with defining polynomial x^3 - 2 with respect to the embedding Ring morphism:
            From: Number Field in a with defining polynomial x^3 - 2
            To:   Algebraic Real Field
            Defn: a |--> 1.259921049894873?
        """
        if self.E.base_field() is QQ:
            return "Period lattice associated to %s"%(self.E)
        else:
            return "Period lattice associated to %s with respect to the embedding %s"%(self.E, self.embedding)

    def __call__(self, P, prec=None):
        r"""
        Return The elliptic logarithm of a point `P`.

        INPUT:

        - ``P`` (point) -- a point on the elliptic curve associated
          with this period lattice.

        - ``prec`` (default: ``None``) -- precision in bits (default
          precision if ``None``).

        OUTPUT:

        (complex number) The elliptic logarithm of the point `P` with
        respect to this period lattice.  If `E` is the elliptic curve
        and `\sigma:K\to\CC` the embedding, the the returned value `z`
        is such that `z\pmod{L}` maps to `\signa(P)` under the
        standard Weierstrass isomorphism from `\CC/L` to `\sigma(E)`.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.period_lattice()
            sage: E.discriminant() > 0
            True
            sage: L.real_flag
            1
            sage: P = E([-1,1])
            sage: P.is_on_identity_component ()
            False
            sage: L(P, prec=96)
            0.4793482501902193161295330101 + 0.985868850775824102211203849...*I
            sage: Q=E([3,5])
            sage: Q.is_on_identity_component()
            True
            sage: L(Q, prec=96)
            1.931128271542559442488585220

        An example with negative discriminant, and a torsion point::

            sage: E = EllipticCurve('11a1')
            sage: L = E.period_lattice()
            sage: E.discriminant() < 0
            True
            sage: L.real_flag
            -1
            sage: P = E([16,-61])
            sage: L(P)
            0.253841860855911
            sage: L.real_period() / L(P)
            5.00000000000000
        """
        return self.elliptic_logarithm(P,prec)

    def basis(self, prec=None, algorithm='sage'):
        r"""
        Return a basis for this period lattice as a 2-tuple.

        INPUT:

        - ``prec`` (default: ``None``) -- precision in bits (default
          precision if ``None``).

        - ``algorithm`` (string, default 'sage') -- choice of
          implementation (for real ambeddings only) between 'sage'
          (native Sage implementation) or 'pari' (use the pari
          library: only available for real embeddings).

        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice is
        `\ZZ\omega_1 + \ZZ\omega_2`.  If the lattice is real then
        `\omega_1` is real and positive, `\Im(\omega_2)>0` and
        `\Re(\omega_1/\omega_2)` is either `0` (for rectangular
        lattices) or `\frac{1}{2}` (for non-rectangular lattices).
        Otherwise, `\omega_1/\omega_2` is in the fundamental region of
        the upper half-plane.  If the latter normalisation is required
        for real lattices, use the function ``normalised_basis()``
        instead.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().basis()
            (2.99345864623196, 2.45138938198679*I)

        This shows that the issue reported at trac \#3954 is fixed::

            sage: E = EllipticCurve('37a')
            sage: b1 = E.period_lattice().basis(prec=30)
            sage: b2 = E.period_lattice().basis(prec=30)
            sage: b1 == b2
            True

        This shows that the issue reported at trac \#4064 is fixed::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().basis(prec=30)[0].parent()
            Real Field with 30 bits of precision
            sage: E.period_lattice().basis(prec=100)[0].parent()
            Real Field with 100 bits of precision

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.basis(64)
            (3.81452977217854509, 1.90726488608927255 + 1.34047785962440202*I)

            sage: emb = K.embeddings(ComplexField())[0]
            sage: L = E.period_lattice(emb)
            sage: w1,w2 = L.basis(); w1,w2
            (-1.37588604166076 - 2.58560946624443*I, -2.10339907847356 + 0.428378776460622*I)
            sage: L.is_real()
            False
            sage: tau = w1/w2; tau
            0.387694505032876 + 1.30821088214407*I
        """
        # We divide into two cases: (1) Q, or a number field with a
        # real embedding; (2) a number field with a complex embedding.
        # In each case the periods are computed by a different
        # internal function.

        if self.is_real():
            return self._compute_periods_real(prec=prec, algorithm=algorithm)
        else:
            return self._compute_periods_complex(prec=prec)

    def normalised_basis(self, prec=None, algorithm='sage'):
        r"""
        Return a normalised basis for this period lattice as a 2-tuple.

        INPUT:

        - ``prec`` (default: ``None``) -- precision in bits (default
          precision if ``None``).

        - ``algorithm`` (string, default 'sage') -- choice of
          implementation (for real ambeddings only) between 'sage'
          (native Sage implementation) or 'pari' (use the pari
          library: only available for real embeddings).

        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice has
        the form `\ZZ\omega_1 + \ZZ\omega_2`.  The basis is normalised
        so that `\omega_1/\omega_2` is in the fundamental region of
        the upper half-plane.  For an alternative normalisation for
        real lattices (with the first period real), use the function
        basis() instead.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().normalised_basis()
            (2.99345864623196, -2.45138938198679*I)

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.normalised_basis(64)
            (1.90726488608927255 - 1.34047785962440202*I, -1.90726488608927255 - 1.34047785962440202*I)

            sage: emb = K.embeddings(ComplexField())[0]
            sage: L = E.period_lattice(emb)
            sage: w1,w2 = L.normalised_basis(); w1,w2
            (-1.37588604166076 - 2.58560946624443*I, -2.10339907847356 + 0.428378776460622*I)
            sage: L.is_real()
            False
            sage: tau = w1/w2; tau
            0.387694505032876 + 1.30821088214407*I
        """
        w1, w2 = periods = self.basis(prec=prec, algorithm=algorithm)
        if self.real_flag:
            periods, mat = normalise_periods(w1,w2)
        return periods


    def _compute_periods_real(self, prec=None, algorithm='sage'):
        r"""
        Internal function to compute the periods (real embedding case).

        INPUT:


        - `prec` (int or ``None`` (default)) -- floating point
          precision (in bits); if None, use the default precision.

        - `algorithm` (string, default 'sage') -- choice of implementation between
          - `pari`: use the pari library

          - `sage`: use a native Sage implementation (with the same underlying algorithm).


        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice has
        the form `\ZZ\omega_1 + \ZZ\omega_2`, `\omega_1` is real and
        `\omega_1/\omega_2` has real part either `0` or `frac{1}{2}`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: embs = K.embeddings(CC)
            sage: Ls = [E.period_lattice(e) for e in embs]
            sage: [L.is_real() for L in Ls]
            [False, False, True]
            sage: Ls[2]._compute_periods_real(100)
            (3.8145297721785450936365098936,
            1.9072648860892725468182549468 + 1.3404778596244020196600112394*I)
            sage: Ls[2]._compute_periods_real(100, algorithm='pari')
            (3.8145297721785450936365098936,
            -1.9072648860892725468182549468 + 1.3404778596244020196600112394*I)
        """
        if prec is None:
            prec = RealField().precision()
        R = RealField(prec)
        C = ComplexField(prec)

        if algorithm=='pari':
            if self.E.base_field() is QQ:
                periods = self.E.pari_curve(prec).omega().python()
                return (R(periods[0]), C(periods[1]))

            from sage.libs.pari.all import pari
            E_pari = pari([R(self.embedding(ai).real()) for ai in self.E.a_invariants()]).ellinit(precision=prec)
            periods = E_pari.omega().python()
            return (R(periods[0]), C(periods[1]))

        if algorithm!='sage':
            raise ValueError, "invalid value of 'algorithm' parameter"

        pi = R.pi()
        if self.real_flag == 1: # positive discriminant
            # Up to now everything has been exact in AA, but now we
            # must go transcendental.  Only now is the desired
            # precision used!
            a, b, c = (R(x) for x in self._abc)
            w1 = R(pi/a.agm(b))   # real
            w2 = C(0,pi/a.agm(c)) # pure imaginary
        else:
            # Up to now everything has been exact in AA, but now we
            # must go transcendental.  Only now is the desired
            # precision used!
            a, b, c = (C(x) for x in self._abc)
            w1 = R(pi/a.abs().agm(a.real().abs())) # real & positive
            w2 = C(pi*C.gen()/a.agm(c))
            w2 = C(w1/2,w2.imag().abs())

        return (w1,w2)

    def _compute_periods_complex(self, prec=None):
        r"""
        Internal function to compute the periods (complex embedding case).

        INPUT:

        - `prec` (int or ``None`` (default)) -- floating point precision (in bits); if None,
          use the default precision.

        OUTPUT:

        (tuple of Complex) `(\omega_1,\omega_2)` where the lattice has
        the form `\ZZ\omega_1 + \ZZ\omega_2` and `(\omega_1/\omega_2)`
        is in te fundamental region of the upper half plane.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: embs = K.embeddings(CC)
            sage: Ls = [E.period_lattice(e) for e in embs]
            sage: [L.is_real() for L in Ls]
            [False, False, True]
            sage: Ls[2]._compute_periods_complex(100)
            (1.9072648860892725468182549468 - 1.3404778596244020196600112394*I, -1.9072648860892725468182549468 - 1.3404778596244020196600112394*I)
        """
        if prec is None:
            prec = RealField().precision()
        C = ComplexField(prec)

        # Up to now everything has been exact in AA, but now we
        # must go transcendental.  Only now is the desired
        # precision used!
        pi = C.pi()
        a, b, c = (C(x) for x in self._abc)
        w1 = pi/a.agm(b)
        w2 = pi*C.gen()/a.agm(c)
        w1w2, mat = normalise_periods(w1,w2)
        return w1w2

    def is_real(self):
        r"""
        Return True if this period lattice is real.

        EXAMPLES::

            sage: f = EllipticCurve('11a')
            sage: f.period_lattice().is_real()
            True

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K,[0,0,0,i,2*i])
            sage: emb = K.embeddings(ComplexField())[0]
            sage: L = E.period_lattice(emb)
            sage: L.is_real()
            False

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: [E.period_lattice(emb).is_real() for emb in K.embeddings(CC)]
            [False, False, True]


        ALGORITHM:

        The lattice is real if it is associated to a real embedding;
        such lattices are stable under conjugation.
        """
        return self.real_flag!=0

    def is_rectangular(self):
        r"""
        Return True if this period lattice is rectangular.

        .. note::

           Only defined for real lattices; a RuntimeError is raised for
           non-real lattices.

        EXAMPLES::

            sage: f = EllipticCurve('11a')
            sage: f.period_lattice().basis()
            (1.26920930427955, 0.634604652139777 + 1.45881661693850*I)
            sage: f.period_lattice().is_rectangular()
            False

        ::

            sage: f = EllipticCurve('37b')
            sage: f.period_lattice().basis()
            (1.08852159290423, 1.76761067023379*I)
            sage: f.period_lattice().is_rectangular()
            True

        ALGORITHM:

        The period lattice is rectangular precisely if the
        discriminant of the Weierstrass equation is positive, or
        equivalently if the number of real components is 2.
        """
        if self.is_real():
            return self.real_flag == +1
        raise RuntimeError, "Not defined for non-real lattices."

    def real_period(self, prec = None, algorithm='sage'):
        """
        Returns the real period of this period lattice.

        INPUT:

        - ``prec`` (int or ``None`` (default)) -- real precision in
          bits (default real precision if ``None``)

        - ``algorithm`` (string, default 'sage') -- choice of
          implementation (for real ambeddings only) between 'sage'
          (native Sage implementation) or 'pari' (use the pari
          library: only available for real embeddings).

        .. note::

           Only defined for real lattices; a RuntimeError is raised for
           non-real lattices.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().real_period()
            2.99345864623196

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.real_period(64)
            3.81452977217854509
        """
        if self.is_real():
            return self.basis(prec,algorithm)[0]
        raise RuntimeError, "Not defined for non-real lattices."

    def omega(self, prec = None):
        r"""
        Returns the real or complex volume of this period lattice.

        INPUT:

        - ``prec`` (int or ``None``(default)) -- real precision in
          bits (default real precision if ``None``)

        OUTPUT:

        (real) For real lattices, this is the real period times the
        number of connected components.  For non-real lattices it is
        the complex area.

        .. note::

           If the curve is defined over `\QQ` and is given by a
           *minimal* Weierstrass equation, then this is the correct
           period in the BSD conjecture, i.e., it is the least real
           period * 2 when the period lattice is rectangular.  More
           generally the product of this quantity over all embeddings
           appears in the generalised BSD formula.


        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().omega()
            5.98691729246392

        This is not a minimal model::

            sage: E = EllipticCurve([0,-432*6^2])
            sage: E.period_lattice().omega()
            0.486109385710056

        If you were to plug the above omega into the BSD conjecture, you
        would get nonsense.   The following works though::

            sage: F = E.minimal_model()
            sage: F.period_lattice().omega()
            0.972218771420113

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.omega(64)
            3.81452977217854509

        A complex example (taken from J.E.Cremona and E.Whitley,
        *Periods of cusp forms and elliptic curves over imaginary
        quadratic fields*, Mathematics of Computation 62 No. 205
        (1994), 407-429)::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0,1-i,i,-i,0])
            sage: L = E.period_lattice(K.embeddings(CC)[0])
            sage: L.omega()
            8.80694160502647
        """
        if self.is_real():
            n_components = (self.real_flag+3)//2
            return self.real_period(prec) * n_components
        else:
            return self.complex_area()

    def basis_matrix(self, prec=None, normalised=False):
        r"""
        Return the basis matrix of this period lattice.

        INPUT:

        - ``prec`` (int or ``None``(default)) -- real precision in
          bits (default real precision if ``None``).

        - ``normalised`` (bool, default None) -- if True and the
          embedding is real, use the normalised basis (see
          ``normalised_basis()``) instead of the default.

        OUTPUT:

        A 2x2 real matrix whose rows are the lattice basis vectors,
        after identifying `\CC` with `\RR^2`.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().basis_matrix()
            [ 2.99345864623196 0.000000000000000]
            [0.000000000000000  2.45138938198679]

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.basis_matrix(64)
            [ 3.81452977217854509 0.000000000000000000]
            [ 1.90726488608927255  1.34047785962440202]

        See \#4388::

            sage: L = EllipticCurve('11a1').period_lattice()
            sage: L.basis_matrix()
            [ 1.26920930427955 0.000000000000000]
            [0.634604652139777  1.45881661693850]
            sage: L.basis_matrix(normalised=True)
            [0.634604652139777 -1.45881661693850]
            [-1.26920930427955 0.000000000000000]

        ::

            sage: L = EllipticCurve('389a1').period_lattice()
            sage: L.basis_matrix()
            [ 2.49021256085505 0.000000000000000]
            [0.000000000000000  1.97173770155165]
            sage: L.basis_matrix(normalised=True)
            [ 2.49021256085505 0.000000000000000]
            [0.000000000000000 -1.97173770155165]
        """
        from sage.matrix.all import Matrix

        if normalised:
            return Matrix([list(w) for w in self.normalised_basis(prec)])

        w1,w2 = self.basis(prec)
        if self.is_real():
            return Matrix([[w1,0],list(w2)])
        else:
            return Matrix([list(w) for w in w1,w2])

    def complex_area(self, prec=None):
        """
        Return the area of a fundamental domain for the period lattice
        of the elliptic curve.

        INPUT:

        - ``prec`` (int or ``None``(default)) -- real precision in
          bits (default real precision if ``None``).

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().complex_area()
            7.33813274078958

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: embs = K.embeddings(ComplexField())
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: [E.period_lattice(emb).is_real() for emb in K.embeddings(CC)]
            [False, False, True]
            sage: [E.period_lattice(emb).complex_area() for emb in embs]
            [6.02796894766694, 6.02796894766694, 5.11329270448345]
        """
        w1,w2 = self.basis(prec)
        return (w1*w2.conjugate()).imag().abs()

    def sigma(self, z, prec = None, flag=0):
        r"""
        Returns the value of the Weierstrass sigma function for this elliptic curve  period lattice.

        INPUT:

        - ``z`` -- a complex number

        - ``prec`` (default: ``None``) -- real precision in bits
            (default real precision if None).

        - ``flag`` --

            0: (default) ???;

            1: computes an arbitrary determination of log(sigma(z))

            2, 3: same using the product expansion instead of theta series. ???

        .. note::

           The reason for the ???'s above, is that the PARI
           documentation for ellsigma is very vague.  Also this is
           only implemented for curves defined over `\QQ`.

        TODO:

        This function does not use any of the PeriodLattice functions
        and so should be moved to ell_rational_field.

        EXAMPLES::

            sage: EllipticCurve('389a1').period_lattice().sigma(CC(2,1))
            2.60912163570108 - 0.200865080824587*I
        """
        if prec is None:
            prec = RealField().precision()
        try:
            return self.E.pari_curve(prec).ellsigma(z, flag)
        except AttributeError:
            raise NotImplementedError, "sigma function not yet implemented for period lattices of curves not defined over Q."

    def curve(self):
        r"""
        Return the elliptic curve associated with this period lattice.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.period_lattice()
            sage: L.curve() is E
            True

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(K.embeddings(RealField())[0])
            sage: L.curve() is E
            True

            sage: L = E.period_lattice(K.embeddings(ComplexField())[0])
            sage: L.curve() is E
            True
        """
        return self.E

    def ei(self):
        r"""
        Return the x-coordinates of the 2-division points of the elliptic curve associated with this period lattice, as elements of QQbar.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.period_lattice()
            sage: L.ei()
            [-1.107159871688768?, 0.2695944364054446?, 0.8375654352833230?]

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(K.embeddings(RealField())[0])
            sage: L.ei()
            [0.?e-19 - 1.122462048309373?*I, 0.?e-19 + 1.122462048309373?*I, -1]

        sage: L = E.period_lattice(K.embeddings(ComplexField())[0])
        sage: L.ei()
        [-1.000000000000000? + 0.?e-1...*I,
        -0.9720806486198328? - 0.561231024154687?*I,
        0.9720806486198328? + 0.561231024154687?*I]
        """
        return self._ei

    def elliptic_logarithm(self, P, prec=None):
        r"""
        Return the elliptic logarithm of a point.

        INPUT:

        - ``P`` (point) -- A point on the elliptic curve associated with this period lattice.

        - ``prec`` (default: ``None``) -- real precision in bits
            (default real precision if None).

        OUTPUT:

        (complex number) The elliptic logarithm of the point `P` with
        respect to this period lattice.  If `E` is the elliptic curve
        and `\sigma:K\to\CC` the embedding, the the returned value `z`
        is such that `z\pmod{L}` maps to `\signa(P)` under the
        standard Weierstrass isomorphism from `\CC/L` to `\sigma(E)`.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.period_lattice()
            sage: E.discriminant() > 0
            True
            sage: L.real_flag
            1
            sage: P = E([-1,1])
            sage: P.is_on_identity_component ()
            False
            sage: L.elliptic_logarithm(P, prec=96)
            0.4793482501902193161295330101 + 0.985868850775824102211203849...*I
            sage: Q=E([3,5])
            sage: Q.is_on_identity_component()
            True
            sage: L.elliptic_logarithm(Q, prec=96)
            1.931128271542559442488585220

        Note that this is actually the inverse of the Weierstrass isomorphism::

            sage: pari(E).ellztopoint(_)
            [3.00000000000000 + 0.E-18*I, 5.00000000000000 + 0.E-18*I]

        An example with negative discriminant, and a torsion point::

            sage: E = EllipticCurve('11a1')
            sage: L = E.period_lattice()
            sage: E.discriminant() < 0
            True
            sage: L.real_flag
            -1
            sage: P = E([16,-61])
            sage: L.elliptic_logarithm(P)
            0.253841860855911
            sage: L.real_period() / L.elliptic_logarithm(P)
            5.00000000000000

        An example where precision is problematic:

            sage: E = EllipticCurve([1, 0, 1, -85357462, 303528987048]) #18074g1
            sage: P = E([4458713781401/835903744, -64466909836503771/24167649046528, 1])
            sage: L = E.period_lattice()
            sage: L.ei()
            [5334.003952567705? - 1.964393150436?e-6*I, 5334.003952567705? + 1.964393150436?e-6*I, -10668.25790513541?]
            sage: L.elliptic_logarithm(P,prec=100)
            0.27656204014107061464076203097

        """
        if not P.curve() is self.E:
            raise ValueError, "Point is on the wrong curve"
        if self.real_flag==0:
            raise NotImplementedError, "elliptic log not implemented for complex embeddings"
        if prec is None:
            prec = RealField().precision()
        prec2 = prec*2
        R = RealField(prec2)
        C = ComplexField(prec2)
        if P.is_zero():
            return ComplexField(prec)(0)
        pi = R.pi()
        a1,a2,a3 = [self.embedding(a) for a in self.E.ainvs()[:3]]
        xP, yP = P.xy()
        wP = 2*yP+a1*xP+a3
        e1,e2,e3 = self._ei
        if self.real_flag==-1:
            z = self._abc[1] # sqrt(e3-e2)
            beta = z.norm()
            alpha = 2*(e3-e2).real()
            a = beta.sqrt()*2
            b = (alpha+2*beta).sqrt()
            #a = 2*z.abs()
            #b = 2*z.real()
            c = (xP-e3+beta)/(xP-e3).sqrt()
        else:
            on_egg = (xP<e3)
            if on_egg:
                y3 = -(a1*e1+a3)/2
                lam = (yP-y3)/(xP-e1)
                x3 = lam*(lam+a1)-a2-xP-e1
                yP = lam*(xP-x3)-yP-a1*x3-a3
                xP = x3
                wP = 2*yP+a1*xP+a3
            a = self._abc[0] # sqrt(e3-e1)
            b = self._abc[1] # sqrt(e3-e2)
            c = (xP-e1).sqrt()
        # So far the values of a,b,c are algebraic (in AA)
        a = R(a)
        b = R(b)
        c = R(c)
        a,b,c = extended_agm_iteration(a,b,c)
        w1, w2 = self.basis(prec)
        if self.real_flag==-1:
            z = (a/c).arcsin()
            if wP*((xP-e3)*(xP-e3)-beta*beta) >= 0:
                z = pi - z
            if wP > 0:
                z += pi
            z /= a
        elif self.real_flag==+1:
            z = (a/c).arcsin()/a
            if wP > 0:
                z = w1 - z
            if on_egg:
                z += w2/2
        return ComplexField(prec)(z)


def reduce_tau(tau):
    r"""
    Transform a point in the upper half plane to the fundamental region.

    INPUT:

    - ``tau`` (complex) -- a complex number with positive imaginary part

    OUTPUT:

    (tuple) `(\tau',[a,b,c,d])` where `a,b,c,d` are integers such that

      - `ad-bc=1`;
      - `\tau`=(a\tau+b)/(c\tau+d)`;
      - `|\tau'|\ge1`;
      - `|\Re(\tau')|\le\frac{1}{2}`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.period_lattice import reduce_tau
        sage: reduce_tau(CC(1.23,3.45))
        (0.230000000000000 + 3.45000000000000*I, [1, -1, 0, 1])
        sage: reduce_tau(CC(1.23,0.0345))
        (-0.463960069171512 + 1.35591888067914*I, [-5, 6, 4, -5])
        sage: reduce_tau(CC(1.23,0.0000345))
        (0.130000000001761 + 2.89855072463768*I, [13, -16, 100, -123])
    """
    assert tau.imag()>0
    tau_orig = tau
    a, b = ZZ(1), ZZ(0)
    c, d = b, a
    k = tau.real().round()
    tau -= k
    a -= k*c
    b -= k*d
    while tau.abs()<0.999:
        tau = -1/tau
        a, b, c, d = c, d, -a, -b
        k = tau.real().round()
        tau -= k
        a -= k*c
        b -= k*d
    assert a*d-b*c==1
    assert tau.abs()>=0.999 and tau.real().abs() <= 0.5
    return tau, [a,b,c,d]

def normalise_periods(w1,w2):
    r"""
    Normalise the period basis `(w_1,w_2)` so that `w_1/w_2` is in the fundamental region.

    INPUT:

    - ``w1,w2`` (complex) -- two complex numbers with non-real ratio

    OUTPUT:

    (tuple) `((\omega_1',\omega_2'),[a,b,c,d])` where `a,b,c,d` are
    integers such that

      - `ad-bc=\pm1`;
      - `(\omega_1',\omega_2') = (a\omega_1+b\omega_2,c\omega_1+d\omega_2)`;
      - `\tau=\omega_1'/\omega_2'` is in the upper half plane;
      - `|\tau|\ge1` and `|\Re(\tau)|\le\frac{1}{2}`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.period_lattice import reduce_tau, normalise_periods
        sage: w1 = CC(1.234, 3.456)
        sage: w2 = CC(1.234, 3.456000001)
        sage: w1/w2    # in lower half plane!
        0.999999999743367 - 9.16334785827644e-11*I
        sage: w1w2, abcd = normalise_periods(w1,w2)
        sage: a,b,c,d = abcd
        sage: w1w2 == (a*w1+b*w2, c*w1+d*w2)
        True
        sage: w1w2[0]/w1w2[1]
        1.23400010389203e9*I
        sage: a*d-b*c # note change of orientation
        -1

    """
    tau = w1/w2
    s = +1
    if tau.imag()<0:
        w2 = -w2
        tau = -tau
        s = -1
    tau, abcd = reduce_tau(tau)
    a, b, c, d = abcd
    if s<0:
        abcd = (a,-b,c,-d)
    return (a*w1+b*w2,c*w1+d*w2), abcd


def extended_agm_iteration(a,b,c):
    r"""
    Internal function for the extended AGM used in elliptic logarithm computation.
    INPUT:

    - ``a``, ``b``, ``c`` (real or complex) -- three real or complex numbers.

    OUTPUT:

    (3-tuple) `(a_0,b_0,c_0)`, the limit of the iteration `(a,b,c) \mapsto ((a+b)/2,\sqrt{ab},(c+\sqrt(c^2+b^2-a^2))/2)`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.period_lattice import extended_agm_iteration
        sage: extended_agm_iteration(RR(1),RR(2),RR(3))
        (1.45679103104691, 1.45679103104691, 3.21245294970054)
        sage: extended_agm_iteration(CC(1,2),CC(2,3),CC(3,4))
        (1.46242448156430 + 2.47791311676267*I,
        1.46242448156430 + 2.47791311676267*I,
        3.22202144343535 + 4.28383734262540*I)

    TESTS::

        sage: extended_agm_iteration(1,2,3)
        Traceback (most recent call last):
        ...
        ValueError: values must be real or complex numbers

    """
    if not isinstance(a, (RealNumber,ComplexNumber)):
        raise ValueError, "values must be real or complex numbers"
    while True:
        a1 = (a + b)/2
        b1 = (a*b).sqrt()
        delta = (b**2 - a**2)/c**2
        f = (1 + (1 + delta).sqrt())/2
        if f.abs()==1:
            return a,b,c
        c*=f
        a,b = a1,b1
