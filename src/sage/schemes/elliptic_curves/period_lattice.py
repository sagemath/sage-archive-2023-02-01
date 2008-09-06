from sage.modules.free_module import FreeModule_generic_pid
from sage.rings.all import ZZ, QQ, RealField
from sage.misc.misc import prec_bits_to_words

class PeriodLattice(FreeModule_generic_pid):
    """
    The class for the period lattice of an algebraic variety.
    """
    pass

class PeriodLattice_ell(PeriodLattice):
    """
    The class for the period lattice of an elliptic curve.

    Currently supported are elliptic curves defined over Q, and
    elliptic curves defined over a number field with a real embedding,
    where the lattice constructed depends on that embedding.
    Extending this support to complex embeddings is not possible using
    pari library functions, but could be implemented using the complex
    AGM.
    """

    def __init__(self, E, real_embedding=None):
        """
        Initializes the period lattice by storing the elliptic curve.

        INPUT:
            E -- an elliptic curve
            real_embedding -- an embedding of the base field K of E
                              into a real field.  If None (the
                              default) uses the first embedding given
                              by K.embeddings(RealField()), except
                              that when K=Q, the built-in coercion is
                              used instead.

        NOTE: No periods are computed on creation of the lattice; see
        the functions basis() and real_period() for precision matters.

        EXAMPLES:
            This function is not normally called directly, but will be
            called by the period_lattice() function of classes
            ell_number_field and ell_rational_field.

            sage: from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
            sage: E = EllipticCurve('37a')
            sage: PeriodLattice_ell(E)
            Period lattice associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = PeriodLattice_ell(E,emb); L
            Period lattice associated to Elliptic Curve defined by y^2  = x^3 + x^2 + a*x + a over Number Field in a with defining polynomial x^3 - 2 with respect to the real embedding Ring morphism:
            From: Number Field in a with defining polynomial x^3 - 2
            To:   Real Field with 53 bits of precision
            Defn: a |--> 1.25992104989487
        """
        self.E = E
        if real_embedding is None and not E.base_field() is QQ:
            try:
                self.real_embedding = E.base_field().embeddings(RealField())[0]
            except IndexError:
                raise IndexError, "Initialization of PeriodLattice not implemened over fields with no real embedding"
        else:
            self.real_embedding = real_embedding

        # From now on, self.real_embedding is None iff field is Q
        if self.real_embedding is None:
            self.n_components = 2 if E.discriminant()>0 else 1
        else:
            self.n_components = 2 if self.real_embedding(E.discriminant())>0 else 1

        PeriodLattice.__init__(self, base_ring=ZZ, rank=2, degree=1, sparse=False)

    def __repr__(self):
        """
        Returns the string representation of this period lattice.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb); L
            Period lattice associated to Elliptic Curve defined by y^2  = x^3 + x^2 + a*x + a over Number Field in a with defining polynomial x^3 - 2 with respect to the real embedding Ring morphism:
            From: Number Field in a with defining polynomial x^3 - 2
            To:   Real Field with 53 bits of precision
            Defn: a |--> 1.25992104989487
        """
        if self.real_embedding is None:
            return "Period lattice associated to %s"%(self.E)
        else:
            return "Period lattice associated to %s with respect to the real embedding %s"%(self.E, self.real_embedding)

    def basis(self, prec=None):
        r"""
        Return a basis for this period lattice as a 2-tuple.

        The basis has the form $[\omega_1, \omega_2]$, where
        $\Im(\omega_1/\omega_2) > 0$ and $\omega_1$ is real.

        INPUT:
            prec -- real precision in bits (default real precision if None)

        OUTPUT:
            omega_1 -- complex number
            omega_2 -- complex number

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().basis()                 # machine-independent
            (2.99345864623195963, 2.45138938198679006*I)

        This shows that the issue reported at trac \#3954 is fixed:
            sage: E = EllipticCurve('37a')
            sage: b1 = E.period_lattice().basis(prec=30)
            sage: b2 = E.period_lattice().basis(prec=30)
            sage: b1 == b2
            True

        This shows that the issue reported at trac \#4064 is fixed:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().basis(prec=30)[0].parent()
            Real Field with 32 bits of precision # 32-bit
            Real Field with 64 bits of precision # 64-bit
            sage: E.period_lattice().basis(prec=100)[0].parent()
            Real Field with 128 bits of precision

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.basis(64)
            (3.81452977217854509, -1.90726488608927255 + 1.34047785962440202*I)

        """
        if prec is None:
            prec = RealField().precision()
        # When the field is Q we get the basis straight from pari:
        if self.real_embedding is None:
            return tuple(self.E.pari_curve(prec).omega().python())

        # Otherwise we refine the precision of the real embedding if
        # necessary, and then call pari:
        from sage.libs.pari.all import pari
        from sage.rings.number_field.number_field import refine_embedding
        self.real_embedding = refine_embedding(self.real_embedding,prec)
        prec_words = prec_bits_to_words(prec)
        E_pari = pari([self.real_embedding(ai) for ai in self.E.a_invariants()]).ellinit(precision=prec_words)
        return tuple(E_pari.omega().python())

    def is_rectangular(self):
        r"""
        Return True precisely if the period lattice of self
        is rectangular.

        EXAMPLES:
            sage: f = EllipticCurve('11a')
            sage: f.period_lattice().basis()              # machine-independent
            (1.26920930427955342, 0.634604652139776711 + 1.45881661693849523*I)

            sage: f.period_lattice().is_rectangular()
            False
            sage: f = EllipticCurve('37b')
            sage: f.period_lattice().basis()              # machine-independent
            (1.08852159290422917, 1.76761067023378948*I)
            sage: f.period_lattice().is_rectangular()
            True

        ALGORITHM:
            The period lattice is rectangular precisely if the
            discriminant of the Weierstrass equation is positive, or
            equivalently if the number of real components is 2.
        """
        return self.n_components == 2

    def real_period(self, prec = None):
        """
        Returns the real period of this period lattice.

        INPUT:
            prec -- real precision in bits (default real precision if None)

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().real_period()  # machine-independent
            2.99345864623195963

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.real_period(64)
            3.81452977217854509
        """
        return self.basis(prec)[0]

    def omega(self, prec = None):
        """
        Returns the real period of this period lattice, scaled by the number of real compnents.

        If self is given by a \emph{minimal Weierstrass equation} then
        this is the correct period in the BSD conjecture, i.e., it is
        the least real period * 2 when the period lattice is
        rectangular.

        INPUT:
            prec -- real precision in bits (default real precision if None)

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().omega()  # machine-independent
            5.98691729246391926

        This is not a minimal model.
            sage: E = EllipticCurve([0,-432*6^2])
            sage: E.period_lattice().omega()  # machine-independent
            0.486109385710056430

        If you were to plug the above omega into the BSD conjecture, you
        would get nonsense.   The following works though:
            sage: F = E.minimal_model()
            sage: F.period_lattice().omega()  # machine-independent
            0.972218771420112860

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.omega(64)
            3.81452977217854509
        """
        return self.real_period(prec) * self.n_components

    def complex_area(self, prec=None):
        """
        Return the area of a fundamental domain for the period lattice
        of the elliptic curve.

        INPUT:
            prec -- real precision in bits (default real precision if None)

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().complex_area()  # machine-independent
            7.33813274078957674

            sage: K.<a> = NumberField(x^3-2)
            sage: emb = K.embeddings(RealField())[0]
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: L = E.period_lattice(emb)
            sage: L.complex_area(64)
            5.11329270448345399
        """
        w1,w2 = self.basis(prec)
        return (w1*w2.imag()).real()

    def sigma(self, z, prec = None, flag=0):
        """
        Returns the value of the Weierstrass sigma function of the lattice
        associated to this elliptic curve E.

        INPUT:
            z -- a complex number
            prec -- real precision in bits (default real precision if None)
            flag -- 0 - default ???
                    1 - computes an arbitrary determination of log(sigma(z))
                    2, 3 - same using the product expansion instead of theta series.
                           ???
        OUTPUT:
            a complex number

        NOTE: The reason for the ???'s above, is that the PARI documentation for
              ellsigma is very vague.  Also this is only implemented for curves
              defined over Q.

        EXAMPLES:
            sage: EllipticCurve('389a1').period_lattice().sigma(CC(2,1))
            2.609121635701083769 - 0.20086508082458695134*I
        """
        if prec is None:
            prec = RealField().precision()
        try:
            return self.E.pari_curve(prec).ellsigma(z, flag)
        except AttributeError:
            raise NotImplementedError, "sigma function not yet implemented for period lattices of curves not defined over Q."
