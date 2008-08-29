from sage.modules.free_module import FreeModule_generic_pid
from sage.rings.all import ZZ, QQ, RealField

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
        -- an elliptic curve
        OUTPUT:
        omega_1 -- complex number
        omega_2 -- complex number

        EXAMPLES:
        sage: E = EllipticCurve('37a')
        sage: E.period_lattice().basis()
        (2.993458646231959629832009979452508177797583791370132985340523378563250356987, 2.451389381986790060854224831866525225349617289144796614656471406129152899999*I)    # 32-bit
        (2.99345864623195962983200997945250817779758379137013298534052337856325035698668290412039406739705147343584052710494728819414438723737202525437537667109326, 2.45138938198679006085422483186652522534961728914479661465647140612915289999925689289113212802918108871268421886966184797547519986661675580167893816478303*I)   # 64-bit

        This shows that the issue reported at trac \#3954 is fixed:
        sage: E = EllipticCurve('37a')
        sage: b1 = E.period_lattice().basis(prec=30)
        sage: b2 = E.period_lattice().basis(prec=30)
        sage: b1 == b2
        True
        """
        # When the field is Q we get the basis straight from pari:
        if self.real_embedding is None:
            return tuple(self.E.pari_curve(prec).omega().python(precision=prec))

        # Otherwise we refine the precision of the real embedding if
        # necessary, and then call pari:
        from sage.libs.pari.all import pari
        if not prec is None:
            from sage.rings.number_field.number_field import refine_embedding
            self.real_embedding = refine_embedding(self.real_embedding,prec)
        E_pari = pari([self.real_embedding(ai) for ai in self.E.a_invariants()]).ellinit(precision=prec)
        return tuple(E_pari.omega().python(precision=prec))

    def is_rectangular(self):
        r"""
        Return True precisely if the period lattice of self
        is rectangular.

        EXAMPLES:
            sage: f = EllipticCurve('11a')
            sage: f.period_lattice().basis()
            (1.269209304279553421688794616754547305219492241830608667967136921230408338613, 0.6346046521397767108443973083772736526097461209153043339835684606152041693064 + 1.458816616938495229330889612903675257159243428952665161469618762450537896609*I)                                     # 32-bit
            (1.26920930427955342168879461675454730521949224183060866796713692123040833861277772269036230592151260731164529627832128743728170032847684397649271401057075, 0.634604652139776710844397308377273652609746120915304333983568460615204169306388861345181152960756303655822648139160643718640850164238421988246357005285375 + 1.45881661693849522933088961290367525715924342895266516146961876245053789660902872639765673368315820172095257526042401249237362183079269125313009041993832*I)             # 64-bit

            sage: f.period_lattice().is_rectangular()
            False
            sage: f = EllipticCurve('37b')
            sage: f.period_lattice().basis()
            (1.088521592904229173504308311539594823105140504301377799086597419750048367573, 1.767610670233789475881323144497815233734289378984139837146363810096739201810*I) # 32-bit
            (1.08852159290422917350430831153959482310514050430137779908659741975004836757281196724615591294512604175793056433512324867543024134734839104934760089947025, 1.76761067023378947588132314449781523373428937898413983714636381009673920180953691706599273805495417094215579677634900614786897226142483706622542207437740*I) # 64-bit
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

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().omega()
            5.986917292463919259664019958905016355595167582740265970681046757126500713973     # 32-bit
            5.98691729246391925966401995890501635559516758274026597068104675712650071397336580824078813479410294687168105420989457638828877447474405050875075334218652        # 64-bit

        This is not a minimal model.
            sage: E = EllipticCurve([0,-432*6^2])
            sage: E.period_lattice().omega()
            0.4861093857100564298972304561738255425526098601923921971195408561181781048715    # 32-bit
            0.486109385710056429897230456173825542552609860192392197119540856118178104871498709353307487730842084963451161261340032305532890657753313985258848453458110       # 64-bit

        If you were to plug the above omega into the BSD conjecture, you
        would get nonsense.   The following works though:
            sage: F = E.minimal_model()
            sage: F.period_lattice().omega()
            0.9722187714201128597944609123476510851052197203847843942390817122363562097430    # 32-bit
            0.972218771420112859794460912347651085105219720384784394239081712236356209742997418706614975461684169926902322522680064611065781315506627970517696906916220      # 64-bit
        """
        return self.basis(prec)[0] * self.n_components

    def complex_area(self):
        """
        Return the area of a fundamental domain for the period lattice
        of the elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice().complex_area()
            7.338132740789576739070721003332305588006176586123733717543180156079096606979     # 32-bit
            7.33813274078957673907072100333230558800617658612373371754318015607909660697945809438214607592923817142863798604406024044503049908233884534256274529672707        # 64-bit
        """
        w1,w2 = self.basis()
        return (w1*w2.imag()).real()

    def sigma(self, z, prec = 50, flag=0):
        """
        Returns the value of the Weierstrass sigma function of the lattice
        associated to this elliptic curve E.

        INPUT:
            z -- a complex number
            prec -- the precision desired
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
