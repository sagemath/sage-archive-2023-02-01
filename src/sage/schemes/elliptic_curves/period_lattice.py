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
from sage.rings.number_field.number_field import refine_embedding
from sage.rings.infinity import Infinity

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
        """
        self.E = E
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
        self.real_flag = 0
        if real:
            self.real_flag = +1
            if embedding(E.discriminant())<0:
                self.real_flag = -1

        PeriodLattice.__init__(self, base_ring=ZZ, rank=2, degree=1, sparse=False)

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
            (3.81452977217854509, 1.90726488608927254 + 1.34047785962440202*I)

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
            (1.90726488608927254 - 1.34047785962440202*I, -1.90726488608927254 - 1.34047785962440202*I)

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

        f = self.E.two_division_polynomial()

        if self.real_flag == 1: # positive discriminant
            # We cannot coerce directly from K to AA:
            ei = PolynomialRing(AA,'x')([self.embedding(c) for c in list(f)]).roots(multiplicities=False)
            ei.sort()
            e1, e2, e3 = ei # e1 < e2 < e3
            a = (e3-e1).sqrt() #
            b = (e3-e2).sqrt() # all real
            c = (e2-e1).sqrt()

            # Up to now everything has been exact in AA, but now we
            # must go transcendental.  Only now is the desired
            # precision used!
            pi = R.pi()
            a = R(a)
            b = R(b)
            c = R(c)
            w1 = R(pi/a.agm(b))   # real
            w2 = C(0,pi/a.agm(c)) # pure imaginary
        else:
            # We cannot coerce directly from K to QQbar:
            ei = PolynomialRing(QQbar,'x')([self.embedding(c) for c in list(f)]).roots(multiplicities=False)
            ei.sort(cmp=lambda z1,z2:cmp(z1.imag(),z2.imag()))
            e1, e3, e2 = ei # so e3 is real
            a = (e3-e1).sqrt()
            c = (e2-e1).sqrt()
            # Up to now everything has been exact in AA, but now we
            # must go transcendental.  Only now is the desired
            # precision used!

            pi = C.pi()
            a = C(a)
            c = C(c)
            w1 = R(pi/a.abs().agm(a.real().abs())) # real & positive
            w2 = C(pi*C.gen()/a.agm(c))
            w2 = C(w1/2,w2.imag().abs())
            # we subtract w1 from w2 to make the real part of w2/w1
            # -1/2 instead of +1/2, for compatibility with pari.

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

        f = self.E.two_division_polynomial()
        e1, e2, e3 = PolynomialRing(QQbar,'x')([self.embedding(c) for c in list(f)]).roots(multiplicities=False)
        a = (e3-e1).sqrt()
        b = (e3-e2).sqrt()
        c = (e2-e1).sqrt()
        # Up to now everything has been exact in AA, but now we
        # must go transcendental.  Only now is the desired
        # precision used!
        pi = C.pi()
        a = C(a)
        b = C(b)
        c = C(c)
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
            (1.26920930427955, 0.634604652139775 + 1.45881661693850*I)
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

    def real_period(self, prec = None):
        """
        Returns the real period of this period lattice.

        INPUT:

        - ``prec`` (int or ``None`` (default)) -- real precision in
          bits (default real precision if ``None``)

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
            return self.basis(prec)[0]
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
            [ 1.90726488608927254  1.34047785962440202]

        See \#4388::

            sage: L = EllipticCurve('11a1').period_lattice()
            sage: L.basis_matrix()
            [ 1.26920930427955 0.000000000000000]
            [0.634604652139775  1.45881661693850]
            sage: L.basis_matrix(normalised=True)
            [0.634604652139775 -1.45881661693850]
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
            [6.02796894766694, 6.02796894766694, 5.11329270448346]
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


