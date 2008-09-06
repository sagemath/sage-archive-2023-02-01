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
            into a real field.  If None (the default) uses the first
            embedding given by K.embeddings(RealField()), except when
            K=Q and the built-in coercion is used instead.

        EXAMPLES:
            See the documentation for
            EllipticCurve_number_field.period_lattice() and
            EllipticCurve_rational_field.period_lattice().

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
        self -- an elliptic curve
        prec -- real precision in bits (default real precision if None)

        OUTPUT:
        omega_1 -- complex number
        omega_2 -- complex number

        EXAMPLES:
        sage: E = EllipticCurve('37a')
        sage: E.period_lattice().basis()
        (2.99345864623195963, 2.45138938198679006*I)  # 32-bit

        This shows that the issue reported at trac \#3954 is fixed:
        sage: E = EllipticCurve('37a')
        sage: b1 = E.period_lattice().basis(prec=30)
        sage: b2 = E.period_lattice().basis(prec=30)
        sage: b1 == b2
        True

        This shows that the issue reported at trac \#4064 is fixed:
        sage: E = EllipticCurve('37a')
        sage: E.period_lattice().basis(prec=30)[0].parent()
        Real Field with 32 bits of precision
        sage: E.period_lattice().basis(prec=100)[0].parent()
        Real Field with 128 bits of precision
        """
        # When the field is Q we get the basis straight from pari:
        if self.real_embedding is None:
            return tuple(self.E.pari_curve(prec).omega().python())

        # Otherwise we refine the precision of the real embedding if
        # necessary, and then call pari:
        from sage.libs.pari.all import pari
        if not prec is None:
            from sage.rings.number_field.number_field import refine_embedding
            self.real_embedding = refine_embedding(self.real_embedding,prec)
        else:
            prec = RealField().precision()
        prec_words = prec_bits_to_words(prec)
        E_pari = pari([self.real_embedding(ai) for ai in self.E.a_invariants()]).ellinit(precision=prec_words)
        return tuple(E_pari.omega().python())

    def is_rectangular(self):
        r"""
        Return True precisely if the period lattice of self
        is rectangular.

        EXAMPLES:
            sage: f = EllipticCurve('11a')
            sage: f.period_lattice().basis()
            (1.26920930427955342, 0.634604652139776711 + 1.45881661693849523*I)  # 32-bit

            sage: f.period_lattice().is_rectangular()
            False
            sage: f = EllipticCurve('37b')
            sage: f.period_lattice().basis()
            (1.08852159290422917, 1.76761067023378948*I)  # 32-bit
            sage: f.period_lattice().is_rectangular()
            True

        ALGORITHM:
            The period lattice is rectangular precisely if the
            discriminant of the Weierstrass equation is positive, or
            equivalently if the number of real components is 2.
        """
        return self.n_components == 2

    def omega(self, prec = None):
        """
        Returns the real period.

        If self is given by a \emph{minimal Weierstrass equation} then
        this is the correct period in the BSD conjecture, i.e., it is
        the least real period * 2 when the period lattice is
        rectangular.

        INPUT:
        prec -- real precision in bits (default real precision if None)

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().omega()
            5.98691729246391926  # 32-bit

        This is not a minimal model.
            sage: E = EllipticCurve([0,-432*6^2])
            sage: E.period_lattice().omega()
            0.486109385710056430  # 32-bit

        If you were to plug the above omega into the BSD conjecture, you
        would get nonsense.   The following works though:
            sage: F = E.minimal_model()
            sage: F.period_lattice().omega()
            0.972218771420112860  # 32-bit
        """
        return self.basis(prec)[0] * self.n_components

    def complex_area(self):
        """
        Return the area of a fundamental domain for the period lattice
        of the elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().complex_area()
            7.33813274078957674  # 32-bit
        """
        w1,w2 = self.basis()
        return (w1*w2.imag()).real()

    def sigma(self, z, prec = None, flag=0):
        """
        Returns the value of the Weierstrass sigma function of the lattice
        associated to this elliptic curve E.

        INPUT:
            z -- a complex number
            prec -- the real precision desired (default real precision
                    if None)
            flag -- 0 - default ???
                    1 - computes an arbitrary determination of log(sigma(z))
                    2, 3 - same using the product expansion instead of theta series.
                           ???
        OUTPUT:
            a complex number

        NOTE: The reason for the ???'s above, is that the PARI documentation for
              ellsigma is very vague.  Also this is only implemented for curves
              defined over Q.
        """
        try:
            return self.E.pari_curve(prec).ellsigma(z, flag)
        except AttributeError:
            raise NotImplementedError, "sigma function not yet implemented for period lattices of curves not defined over Q."
