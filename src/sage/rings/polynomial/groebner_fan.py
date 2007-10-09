r"""
Groebner Fans

SAGE provides much of the functionality of gfan, which is a software
package whose main function is to enumerate all reduced Gr\"obner
bases of a polynomial ideal. The reduced Gr\"obner bases yield the
maximal cones in the Gr\"obner fan of the ideal. Several
subcomputations can be issued and additional tools are included. Among
these the highlights are:

\begin{itemize}

\item Commands for computing tropical varieties.

\item Interactive walks in the Gr\"obner fan of an ideal.

\item Commands for graphical renderings of Gr\"obner fans
      and monomial ideals.

\end{itemize}


AUTHORS:

   -- Anders Nedergaard Jensen: Wrote the gfan C++ program, which
      implements algorithms many of which were invented by Jensen,
      Komei Fukuda, and Rekha Thomas.   All the underlying hard work
      of the Groebner fans functionality of \sage depends on
      this C++ program.

   -- William Stein (2006-04-20): Wrote first version of the \sage
      code for working with Groebner fans.

   -- Tristram Bogart (bogart@math): the design of the \sage interface
      to gfan is joint work with Tristram Bogart, who also supplied
      numerous examples.

EXAMPLES:
    sage: x,y = QQ['x,y'].gens()
    sage: i = ideal(x^2 - y^2 + 1)
    sage: g = i.groebner_fan()
    sage: g.reduced_groebner_bases()
    [[x^2 - y^2 + 1], [-x^2 + y^2 - 1]]

TESTS:
    sage: x,y = QQ['x,y'].gens()
    sage: i = ideal(x^2 - y^2 + 1)
    sage: g = i.groebner_fan()
    sage: g == loads(dumps(g))
    True

"""

__doc_exclude = ['to_intvec', 'multiple_replace', 'forall', \
                 'SageObject', 'gfan', 'is_MPolynomialIdeal', \
                 'is_MPolynomialRing', 'MPolynomialRing', \
                 'QQ', 'ZZ']

import os

import string

import pexpect

from sage.misc.multireplace import multiple_replace
from sage.misc.misc import forall

from sage.structure.sage_object import SageObject
from sage.interfaces.gfan import gfan
from multi_polynomial_ideal import is_MPolynomialIdeal
from multi_polynomial_ring import is_MPolynomialRing, MPolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ

def to_intvec(w):
    try:
        w = eval(w.replace(' ', ','))
    except SyntaxError, msg:
        raise SyntaxError, "%s\n%s"%(w,msg)
    # now w is a tuple of Python ints
    M = ZZ**len(w)    # make a free module
    return M(w)


class GroebnerFan(SageObject):
    def _repr_(self):
        return "Groebner fan of the ideal:\n%s"%self.__ideal

    def __init__(self, I, is_groebner_basis=False, symmetry=None, verbose=False):
        """
        INPUT:
            I -- ideal in a multivariate polynomial ring
            is_groebner_basis -- bool (default False).  if True, then
                         I.gens() must be a Groebner basis with
                         respect to the standard degree lexicographic
                         term order.
            symmetry -- default: None; if not None, describes symmetries
                        of the ideal
            verbose  -- default: False; if True, printout useful info
                        during computations

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: I = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y])
            sage: G = I.groebner_fan(); G
            Groebner fan of the ideal:
            Ideal (x^2*y - z, -x + y^2*z, x*z^2 - y) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        self.__is_groebner_basis = is_groebner_basis
        self.__symmetry = symmetry
        if symmetry:
            print "WARNING! Symmetry option not yet implemented!!"
        self.__verbose = verbose
        if not is_MPolynomialIdeal(I):
            raise TypeError, "I must be a multivariate polynomial ideal"
        S = I.ring()
        R = S.base_ring()
        if not R.is_field():
            raise NotImplementedError, "Groebner fan computation only implemented over fields"
        if not (R is QQ or (R.is_finite() and R.is_prime() and R.order() <= 32749)):
            # sage: previous_prime (2^15)
            # 32749
            raise NotImplementedError, "Groebner fan computation only implemented over Q or GF(p) for p <= 32749."
        if S.ngens() > 52:
            raise NotImplementedError, "Groebner fan computation only implemented for rings in at most 52 variables."

        self.__ideal = I
        self.__ring = S

    def __eq__(self,right):
        return type(self) == type(right) and self.ideal() == right.ideal()

    def ideal(self):
        """
        Return the ideal the was used to define this Groebner fan.
        """
        return self.__ideal

    def _gfan_maps(self):
        """
        INPUT:
            none
        OUTPUT:
            -- map from SAGE ring to gfan ring
            -- map from gfan ring to SAGE ring

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: G._gfan_maps()
            (Ring morphism:
              From: Multivariate Polynomial Ring in x, y, z over Rational Field
              To:   Multivariate Polynomial Ring in a, b, c over Rational Field
              Defn: x |--> a
                    y |--> b
                    z |--> c,
             Ring morphism:
              From: Multivariate Polynomial Ring in a, b, c over Rational Field
              To:   Multivariate Polynomial Ring in x, y, z over Rational Field
              Defn: a |--> x
                    b |--> y
                    c |--> z)
        """
        try:
            return self.__gfan_maps
        except AttributeError:
            S = self.__ring
            n = S.ngens()

            # Define a polynomial ring in n variables
            # that are named a,b,c,d, ..., z, A, B, C, ...
            R = S.base_ring()
            T = MPolynomialRing(R, n, string.ascii_letters[:n])

            # Define the homomorphism that sends the
            # generators of S to the generators of T.
            phi = S.hom(T.gens())

            # Define the homomorphism that sends the
            # generators of T to the generators of S.
            psi = T.hom(S.gens())
            self.__gfan_maps = (phi, psi)
            return self.__gfan_maps

    def _gfan_vardict(self):
        try:
            return self.__gfan_vardict
        except AttributeError:
            S = self.__ring
            n = S.ngens()
            L = string.ascii_letters[:n]
            d = {}
            v = S.variable_names()
            for i in range(n):
                d[L[i]] = v[i]
            self.__gfan_vardict = d
            return d

    def _gfan_ideal(self):
        """
        Return the ideal in gfan's notation (with variables
        mapped to a,b,c, etc.)

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: G._gfan_ideal()
            '{a^2*b - c, b^2*c - a, a*c^2 - b}'
        """
        try:
            return self.__gfan_ideal
        except AttributeError:
            to_gfan, _ = self._gfan_maps()
            J = to_gfan(self.__ideal)
            s = str(J.gens())
            s = s.replace('(','{').replace(')','}').replace(',}','}')
            self.__gfan_ideal = s
            return s

    def ring(self):
        """
        Return the multivariate polynomial ring.
        """
        return self.__ring

    def _gfan_reduced_groebner_bases(self):
        try:
            return self.__gfan_reduced_groebner_bases
        except AttributeError:
            B = self.gfan()
            self.__gfan_reduced_groebner_bases = B
            return B

    def characteristic(self):
        """
        Return the characteristic of the base ring.
        """
        return self.__ring.characteristic()

    def universal_groebner_basis(self):
        try:
            return self.__universal_groebner_basis
        except AttributeError:
            U0 = self.gfan(cmd='polynomialsetunion',
                          I=self._gfan_reduced_groebner_bases().replace(' ',','))
            U = multiple_replace(self._gfan_vardict(), U0).split(' ')[1:-1]
            self.__universal_groebner_basis = U
            return U

    def reduced_groebner_bases(self):
        """
        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: X = G.reduced_groebner_bases()
            sage: len(X)
            33
            sage: X[0]
            [z^15 - z, y - z^11, x - z^9]
            sage: X[0].ideal()
            Ideal (z^15 - z, y - z^11, x - z^9) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: X[:5]
            [[z^15 - z, y - z^11, x - z^9],
            [-y + z^11, y*z^4 - z, y^2 - z^8, x - z^9],
            [-y^2 + z^8, y*z^4 - z, y^2*z^3 - y, y^3 - z^5, x - y^2*z],
            [-y^3 + z^5, y*z^4 - z, y^2*z^3 - y, y^4 - z^2, x - y^2*z],
            [-y^4 + z^2, y^6*z - y, y^9 - z, x - y^2*z]]
        """
        try:
            return self.__reduced_groebner_bases
        except AttributeError:
            G0 = self._gfan_reduced_groebner_bases()

            # change the variable names back using one big substitution?
            G = multiple_replace(self._gfan_vardict(), G0)
            G = G.replace('{{','').replace('}}','').split('} {')
            G0 = G0.replace('{{','').replace('}}','').split('} {')
            S = self.__ring
            X = [ReducedGroebnerBasis(self, [S(f) for f in G[i].split()], G0[i]) for i in range(len(G))]
            self.__reduced_groebner_bases = X
            return X

    def _gfan_mod(self):
        """
        Return the extra options to the gfan command that are used by
        this object to account for working modulo a prime or in the
        presence of extra symmetries.
        """
        try:
            return self.__gfan_mod
        except AttributeError:
            mod = ''
            p = self.characteristic()
            if p != 0:
                mod += ' --mod %s'%p
            else:
                mod += ''

            if self.__is_groebner_basis:
                mod += ' -g'

            if self.__symmetry:
                mod += ' --symmetry'

            self.__gfan_mod = mod
            return self.__gfan_mod

    def gfan(self, cmd='', I=None, format=True):
        if I is None:
            I = self._gfan_ideal()
        # todo -- put something in here (?) when self.__symmetry isn't None...
        cmd += self._gfan_mod()
        s = gfan(I, cmd, verbose=self.__verbose, format=format)
        if s.strip() == '{':
            raise RuntimeError, "Error running gfan command %s on %s"%(cmd, self)
        return s

    def __iter__(self):
        for x in self.reduced_groebner_bases():
            yield x

    def __getitem__(self, i):
        return self.reduced_groebner_bases()[i]

    def buchberger(self):
        """
        Computes and returns a lexicographic reduced Groebner basis for the ideal.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x - z^3, y^2 - x + x^2 - z^3*x]).groebner_fan()
            sage: G.buchberger ()
            [-z^3 + y^2, -z^3 + x]
        """
        try:
            return self.__buchberger
        except AttributeError:
            B = self.gfan(cmd='buchberger')
            B = multiple_replace(self._gfan_vardict(), B)
            S = self.__ring
            B = [S(f) for f in B[1:-1].split()]
            self.__buchberger = B
            return B

    def fvector(self):
        """
        Return the f-vector for the Grobner fan.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.fvector ()
            (1, 4, 3)
        """
        try:
            return self.__fvector
        except AttributeError:
            f = self.gfan(cmd='fvector', I=self._gfan_reduced_groebner_bases().replace(' ',','))
            f = to_intvec(f)
            self.__fvector = f
            return f


    def homogeneity_space(self):
        """
        Return the homogeneity space of a the list of polynomials that
        define this Groebner fan.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: H = G.homogeneity_space()
        """
        try:
            return self.__homogeneity_space
        except AttributeError:
            h = self.gfan(cmd='homogeneityspace', format=False)
            self.__homogeneity_space = h
            return h

    def render(self, file, larger=False, shift=0, show=False):
        """
        Render a Groebner fan as an xfig file.

        More precisely, the output is a drawing of the Groebner
        fan intersected with a traingle.  The corners of the
        triangle are (1,0,0) to the right, (0,1,0) to the left
        and (0,0,1) at the top.  If there are more than three
        variables in the ring we extend these coordinates
        with zeros.

        INPUT:
            shift -- shift the positions of the variables in
                     the drawing.  For example, with shift=1,
                     the corners will be b (right), c (left),
                     and d (top).  The shifting is done modulo
                     the number of variables in the polynomial
                     ring.   The default is 0.
            larger -- bool (default: False); if True, make
                     the triangle larger so that the shape of
                     of the Groebner region appears.
            show  -- bool (default: True); if True, pop up
                     xfig displaying the resulting file (this
                     requires xfig).

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.render('a.fig')

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage.: G.render('a.fig', show=True, larger=True)
        """
        cmd = 'render'
        if shift:
            cmd += ' --shiftVariables %s'%shift
        if larger:
            cmd += ' -L'
        s = self.gfan(cmd, I=self._gfan_reduced_groebner_bases().replace(' ',','), format=False)
        open(file,'w').write(s)
        if show:
            os.system('xfig %s 2>/dev/null 1>/dev/null&'%file)

    def _gfan_stats(self):
        """
        Return various statistics about this Groebner fan.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G._gfan_stats()
            {'Number of reduced Groebner bases': 3, 'Maximal total degree of a Groebner basis': 4, 'Dimension of homogeneity space': 0, 'Number of variables': 2, 'Minimal total degree of a Groebner basis': 2}
        """
        try:
            return self.__stats
        except AttributeError:
            s = self.gfan(cmd='stats', I=self._gfan_reduced_groebner_bases().replace(' ',','), format=False)
            d = {}
            for v in s.split('\n'):
                if len(v) > 0:
                    a,b = v.split(':')
                    d[a] = ZZ(b)
            self.__stats = d
            return d

    def dimension_of_homogeneity_space(self):
        """
        Return the dimension of the homogeneity space.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.dimension_of_homogeneity_space()
            0
        """
        return self._gfan_stats()['Dimension of homogeneity space']

    def maximal_total_degree_of_a_groebner_basis(self):
        """
        Return the maximal total degree of any Groebner basis.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.maximal_total_degree_of_a_groebner_basis()
            4
        """
        return self._gfan_stats()['Maximal total degree of a Groebner basis']

    def minimal_total_degree_of_a_groebner_basis(self):
        """
        Return the minimal total degree of any Groebner basis.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.minimal_total_degree_of_a_groebner_basis()
            2
        """
        return self._gfan_stats()['Minimal total degree of a Groebner basis']

    def number_of_reduced_groebner_bases(self):
        """
        Return the number of reduced Groebner bases.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.number_of_reduced_groebner_bases()
            3
        """
        return self._gfan_stats()['Number of reduced Groebner bases']

    def number_of_variables(self):
        """
        Return the number of variables.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.number_of_variables()
            2

            sage: R = PolynomialRing(QQ,'x',10)
            sage: R.inject_variables(globals())
            Defining x0, x1, x2, x3, x4, x5, x6, x7, x8, x9
            sage: G = ideal([x0 - x9, sum(R.gens())]).groebner_fan()
            sage: G.number_of_variables()
            10
        """
        return self.__ring.ngens()

    def tropical_basis(self, check=True):
        """
        Return a tropical basis for the tropical curve associated to
        this ideal.

        INPUT:
            check -- bool (default: True); if True raises a ValueError
                     exception if this ideal does not define a tropical curve
                     (i.e., the condition that R/I has dimension equal
                     to 1 + the dimension of the homogeneity space is
                     not satisfied).
        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3, order='lex')
            sage: G = R.ideal([y^3-3*x^2, z^3-x-y-2*y^3+2*x^2]).groebner_fan()
            sage: G
            Groebner fan of the ideal:
            Ideal (-3*x^2 + y^3, 2*x^2 - x - 2*y^3 - y + z^3) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: G.tropical_basis ()
            [-4*x^2 - x - y + z^3, -3*x^2 + y^3]
        """
        try:
            return self.__tropical_basis
        except AttributeError:
            pass
        cmd = 'tropicalbasis'

        I = self.ideal()
        hom, _ = forall(I.gens(), lambda x : x.is_homogeneous())
        if not hom:
            cmd += ' -h'
        if check:
            if I.dimension() != 1 + self.dimension_of_homogeneity_space():
                raise ValueError, "The ideal does not define a tropical curve."

        B = self.gfan(cmd)
        B = multiple_replace(self._gfan_vardict(), B)[1:-1]
        S = self.__ring
        B = [S(f) for f in B.split()]
        self.__tropical_basis = B
        return B

    def interactive(self, *args, **kwds):
        """
        See the documentation for self[0].interative()
        """
        self[0].interactive(*args, **kwds)


class ReducedGroebnerBasis(SageObject, list):
    def __init__(self, groebner_fan, gens, gfan_gens):
        self.__groebner_fan = groebner_fan
        list.__init__(self, gens)
        self.__gfan_gens = '{' + gfan_gens.replace(' ',',') + '}'

    def _repr_(self):
        return list.__repr__(self)

    def _gfan_gens(self):
        return self.__gfan_gens

    def _gfan(self):
        return self.__groebner_fan

    def interactive(self, latex=False, flippable=False, wall=False,
                    inequalities=False, weight=False):
        """
        Do an interactive walk of the Groebner fan starting at this
        reduced Groebner basis.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage.: G[0].interactive()
            Initializing gfan interactive mode
            *********************************************
            *     Press control-C to return to SAGE     *
            *********************************************
            ....
        """
        cmd = 'gfan_interactive'
        if latex:
            cmd += ' -L'
        if flippable:
            cmd += ' -f'
        if wall:
            cmd += ' -w'
        if inequalities:
            cmd += ' -i'
        if weight:
            cmd += ' -W'
        cmd += self.__groebner_fan._gfan_mod()
        E = pexpect.spawn(cmd)
        print "Initializing gfan interactive mode"
        #E.sendline(self._gfan_ideal())
        E.sendline(self.__gfan_gens)
        print "*"*45
        print "*     Press control-C to return to SAGE     *"
        print "*"*45
        try:
            E.interact()
        except OSError:
            print "Returning to SAGE."

    def groebner_cone(self, restrict=False):
        """
        Return defining inequalities for the full-dimensional
        Groebner cone associated to this marked minimal reduced
        Groebner basis.

        INPUT:
            restrict -- bool (default: False); if True, add an inequality for
                        each coordinate, so that the cone is restricted to
                        the positive orthant.

        OUTPUT:
            tuple of integer vectors

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G[1].groebner_cone()
            ((-1, 2), (1, -1))
            sage: [g.groebner_cone() for g in G]
            [((0, 1), (1, -2)), ((-1, 2), (1, -1)), ((-1, 2), (-1, 1), (1, 0))]
            sage: G[1].groebner_cone(restrict=True)
            ((-1, 2), (1, -1), (1, 0), (0, 1))
        """
        try:
            return self.__groebner_cone[restrict]
        except AttributeError:
            self.__groebner_cone = {}
        except KeyError:
            pass
        cmd = 'groebnercone'
        if restrict:
            cmd += ' --restrict'
        gf = self.__groebner_fan
        c = gf.gfan(cmd=cmd, I=self.__gfan_gens)
        c = c.replace(') (','),(')
        v = c[1:-1].split(',')
        v = tuple([to_intvec(x) for x in v])
        self.__groebner_cone[restrict] = v
        return v


    def ideal(self):
        """
        Return the ideal generated by this basis.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x - z^3, y^2 - 13*x]).groebner_fan()
            sage: G[0].ideal()
            Ideal (-13*z^3 + y^2, -z^3 + x) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.__groebner_fan.ring().ideal(self)

    def weight_vector(self):
        """
        Return the weight vector of this reduced Groebner basis.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x - z^3, y^2 - 3*x]).groebner_fan()
            sage: [g.weight_vector() for g in G]
            [(4, 2, 1), (3, 1, 1), (4, 3, 2)]

            sage: R.<x,y,z> = PolynomialRing(GF(3),3)
            sage: G = R.ideal([x - z^3, y^2 - 3*x]).groebner_fan()
            sage: [g.weight_vector() for g in G]
            [(4, 1, 1), (2, 1, 1)]
        """
        try:
            return self.__weight_vector
        except AttributeError:
            cmd = 'weightvector'
            gf = self.__groebner_fan
            w = to_intvec(gf.gfan(cmd='weightvector', I=self.__gfan_gens))
            self.__weight_vector = w
            return w


