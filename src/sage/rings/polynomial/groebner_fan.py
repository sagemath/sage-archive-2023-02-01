r"""
Groebner Fans

\SAGE provides much of the functionality of gfan, which is a software
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
      Komei Fukuda, and Rekha Thomas. All the underlying hard work of
      the Gr\"obner fans functionality of \sage depends on this C++
      program.

   -- William Stein (2006-04-20): Wrote first version of the \SAGE
      code for working with Groebner fans.

   -- Tristram Bogart (bogart@math): the design of the \SAGE interface
      to gfan is joint work with Tristram Bogart, who also supplied
      numerous examples.

   -- Marshall Hampton (2008-03-25): Rewrote various functions to use
      gfan-0.3.
      This is still a work in progress, comments are appreciated on
      sage-devel@googlegroups.com (or personally at hamptonio@gmail.com).

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

REFERENCES:
     Anders N. Jensen; Gfan, a software system for Gr\"obner fans;
     available at
     \url{http://www.math.tu-berlin.de/~jensen/software/gfan/gfan.html}
"""

import os

from string import ascii_letters

import pexpect

from sage.misc.multireplace import multiple_replace
from sage.misc.misc import forall, tmp_filename

from sage.structure.sage_object import SageObject
from sage.interfaces.gfan import gfan
from sage.rings.polynomial.multi_polynomial_ideal import is_MPolynomialIdeal
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing, MPolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.plot.plot import Graphics, line, polygon

def prefix_check(str_list):
    """
    Checks if any strings in a list are prefixes of another string in
    the list.

    EXAMPLES:
        sage: from sage.rings.polynomial.groebner_fan import prefix_check
        sage: prefix_check(['z1','z1z1'])
        False
        sage: prefix_check(['z1','zz1'])
        True
    """
    for index1 in range(len(str_list)):
        for index2 in range(len(str_list)):
            string1 = str_list[index1]
            string2 = str_list[index2]
            if index1 != index2 and str_list[index1][0:len(string2)].find(string2) != -1:
                return False
    return True


def max_degree(list_of_polys):
    """
    Computes the maximum degree of a list of polynomials

    EXAMPLES:
        sage: from sage.rings.polynomial.groebner_fan import max_degree
        sage: R.<x,y> = PolynomialRing(QQ,2)
        sage: p_list = [x^2-y,x*y^10-x]
        sage: max_degree(p_list)
        11.0
    """
    return max([float(qf.degree()) for qf in list_of_polys])

class PolyhedralCone(SageObject):

    def __init__(self, gfan_polyhedral_cone, ring = QQ):
        """
        Converts polymake/gfan data on a polyhedral cone into a sage
        class.  Currently (18-03-2008) needs a lot of work.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.facets()
            [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
        """
        cone_keys = ['AMBIENT_DIM','DIM','IMPLIED_EQUATIONS', 'LINEALITY_DIM', 'LINEALITY_SPACE','FACETS', 'RELATIVE_INTERIOR_POINT']
        poly_lines = gfan_polyhedral_cone.split('\n')
        self.cone_dict = {}
        key_ind = 0
        cur_key = None
        for ting in poly_lines:
            if cone_keys.count(ting) > 0:
                cur_key = ting
                self.cone_dict[cur_key] = []
            elif cur_key and ting != '':
                self.cone_dict[cur_key].append(ting)
        self._facets = []
        for facet in self.cone_dict['FACETS']:
            temp_facet = facet.split('\t')[0]
            temp_facet = temp_facet.split(' ')
            temp_facet = [int(x) for x in temp_facet]
            self._facets.append(temp_facet)
        self._ambient_dim = int(self.cone_dict['AMBIENT_DIM'][0])
        self._dim = int(self.cone_dict['DIM'][0])
        self._lineality_dim = int(self.cone_dict['LINEALITY_DIM'][0])
        rel_int_pt_str = self.cone_dict['RELATIVE_INTERIOR_POINT'][0]
        self._relative_interior_point = [int(q) for q in rel_int_pt_str.split(' ')]

    def _repr_(self):
        """
        Returns a basic description of the polyhedral cone.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a # indirect doctests
            Polyhedral cone in 3 dimensions of dimension 3
        """
        return "Polyhedral cone in %s dimensions of dimension %s"%(str(self.ambient_dim()), str(self.dim()))

    def facets(self):
        """
        Returns the inward facet normals of the Groebner cone.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.facets()
            [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
        """
        return self._facets

    def ambient_dim(self):
        """
        Returns the ambient dimension of the Groebner cone.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.ambient_dim()
            3
        """
        return self._ambient_dim

    def dim(self):
        """
        Returns the dimension of the Groebner cone.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.dim()
            3
        """
        return self._dim

    def lineality_dim(self):
        """
        Returns the lineality dimension of the Groebner cone.  This is
        the just the difference between the ambient dimension and the
        dimension of the cone.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.lineality_dim()
            0
        """
        return self._lineality_dim

    def relative_interior_point(self):
        """
        Returns a point in the relative interior of the Groebner cone.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf[0].groebner_cone()
            sage: a.relative_interior_point()
            [1, 1, 1]
        """
        return self._relative_interior_point

class PolyhedralFan(SageObject):
    def __init__(self, gfan_polyhedral_fan):
        """
        Converts polymake/gfan data on a polyhedral fan into a sage
        class.  Currently (18-03-2008) needs a lot of work.

        INPUT:
            gfan_polyhedral_fan -- output from gfan of a polyhedral fan.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: i2 = ideal(x*z + 6*y*z - z^2, x*y + 6*x*z + y*z - z^2, y^2 + x*z + y*z)
            sage: gf2 = i2.groebner_fan(verbose = False)
            sage: pf = gf2.polyhedralfan()
            sage: pf.rays()
            [[1, 0, 0], [-2, -1, 0], [1, 1, 0], [0, -1, 0], [-1, 1, 0]]
        """
        fan_keys = ['AMBIENT_DIM','DIM','LINEALITY_DIM','RAYS','N_RAYS',
                    'LINEALITY_SPACE','ORTH_LINEALITY_SPACE','F_VECTOR',
                    'CONES','MAXIMAL_CONES','PURE']
        poly_lines = gfan_polyhedral_fan.split('\n')
        self.fan_dict = {}
        key_ind = 0
        cur_key = None
        for ting in poly_lines:
            if fan_keys.count(ting) > 0:
                cur_key = ting
                self.fan_dict[cur_key] = []
            elif cur_key and ting != '':
                self.fan_dict[cur_key].append(ting)
        self._ambient_dim = int(self.fan_dict['AMBIENT_DIM'][0])
        self._dim = int(self.fan_dict['DIM'][0])
        self._lineality_dim = int(self.fan_dict['LINEALITY_DIM'][0])
        self._rays = []
        for ray in self.fan_dict['RAYS']:
            temp_ray = ray.split('\t')[0]
            temp_ray = temp_ray.split(' ')
            temp_ray = [int(x) for x in temp_ray]
            self._rays.append(temp_ray)
        self._str = gfan_polyhedral_fan

    def _repr_(self):
        """
        Returns a basic description of the polyhedral fan.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: pf = gf.polyhedralfan()
            sage: pf # indirect doctest
            Polyhedral fan in 3 dimensions of dimension 3
        """
        return "Polyhedral fan in %s dimensions of dimension %s"%(str(self.ambient_dim()), str(self.dim()))

    def _str_(self):
        r"""
        Returns the raw output of gfan as a string.  This should only
        be needed internally as all relevant output is converted to
        sage objects.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: pf = gf.polyhedralfan()
            sage: pf._str_()
            '_application PolyhedralFan\n_version 2.2\n_type PolyhedralFan\n\nAMBIENT_DIM\n3\n\nDIM\n3\n\nLINEALITY_DIM\n0\n\nRAYS\n1 0 0\t# 0\n0 1 0\t# 1\n0 0 1\t# 2\n\nN_RAYS\n3\n\nLINEALITY_SPACE\n\nORTH_LINEALITY_SPACE\n0 0 1\n0 1 0\n1 0 0\n\nF_VECTOR\n1 3 3 1\n\nCONES\n{}\t# Dimension 0\n{0}\t# Dimension 1\n{1}\n{2}\n{0 1}\t# Dimension 2\n{0 2}\n{1 2}\n{0 1 2}\t# Dimension 3\n\nMAXIMAL_CONES\n{0 1 2}\t# Dimension 3\n\nPURE\n1\n'
        """
        return self._str

    def ambient_dim(self):
        """
        Returns the ambient dimension of the Groebner fan.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf.polyhedralfan()
            sage: a.ambient_dim()
            3
        """
        return self._ambient_dim

    def dim(self):
        """
        Returns the dimension of the Groebner fan.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf.polyhedralfan()
            sage: a.dim()
            3
        """
        return self._dim

    def lineality_dim(self):
        """
        Returns the lineality dimension of the Groebner fan.
        This is the just the difference between the ambient dimension
        and the dimension of the cone.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-2]).groebner_fan()
            sage: a = gf.polyhedralfan()
            sage: a.lineality_dim()
            0
        """
        return self._lineality_dim

    def rays(self):
        """
        Returns a list of rays of the polyhedral fan.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: i2 = ideal(x*z + 6*y*z - z^2, x*y + 6*x*z + y*z - z^2, y^2 + x*z + y*z)
            sage: gf2 = i2.groebner_fan(verbose = False)
            sage: pf = gf2.polyhedralfan()
            sage: pf.rays()
            [[1, 0, 0], [-2, -1, 0], [1, 1, 0], [0, -1, 0], [-1, 1, 0]]
        """
        return self._rays


class GroebnerFan(SageObject):

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
        if prefix_check([str(R_gen) for R_gen in I.ring().gens()]) != True:
            raise RuntimeError, "Ring variables cannot contain each other as prefixes"
        S = I.ring()
        R = S.base_ring()
        if not R.is_field():
            raise NotImplementedError, "Groebner fan computation only implemented over fields"
        if not (R is QQ or (R.is_finite() and R.is_prime_field() and R.order() <= 32749)):
            # sage: previous_prime (2^15)
            # 32749
            raise NotImplementedError, "Groebner fan computation only implemented over Q or GF(p) for p <= 32749."
        if S.ngens() > 52:
            raise NotImplementedError, "Groebner fan computation only implemented for rings in at most 52 variables."

        self.__ideal = I
        self.__ring = S

    def _repr_(self):
        """
        Describes the Groebner fan and its corresponding ideal.

        EXAMPLES:
            sage: R.<q,u> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([q-u,u^2-1]).groebner_fan()
            sage: gf # indirect doctest
            Groebner fan of the ideal:
            Ideal (q - u, u^2 - 1) of Multivariate Polynomial Ring in q, u over Rational Field

        """
        return "Groebner fan of the ideal:\n%s"%self.__ideal

    def __eq__(self,right):
        """
        Tests equality of Groeber fan objects.

        EXAMPLES:
            sage: R.<q,u> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([q^2-u,u^2-q]).groebner_fan()
            sage: gf2 = R.ideal([u^2-q,q^2-u]).groebner_fan()
            sage: gf.__eq__(gf2)
            True
        """
        return type(self) == type(right) and self.ideal() == right.ideal()

    def ideal(self):
        """
        Return the ideal the was used to define this Groebner fan.

        EXAMPLES:
            sage: R.<x1,x2> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x1^3-x2,x2^3-2*x1-2]).groebner_fan()
            sage: gf.ideal()
            Ideal (x1^3 - x2, x2^3 - 2*x1 - 2) of Multivariate Polynomial Ring in x1, x2 over Rational Field
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
            T = MPolynomialRing(R, n, ascii_letters[:n])

            # Define the homomorphism that sends the
            # generators of S to the generators of T.
            phi = S.hom(T.gens())

            # Define the homomorphism that sends the
            # generators of T to the generators of S.
            psi = T.hom(S.gens())
            self.__gfan_maps = (phi, psi)
            return self.__gfan_maps

    def _gfan_ring(self):
        """
        Return the ring in gfan's notation

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: G._gfan_ring()
            'Q[x, y, z]'
        """
        if self.__ideal.base_ring() == QQ:
            ring_str = 'Q' + str(self.__ideal.ring().gens()).replace('(','[').replace(')',']')
        else:
            ring_str = 'Z/' + str(self.__ideal.base_ring().characteristic()) + 'Z'+ str(self.__ideal.ring().gens()).replace('(','[').replace(')',']')
        return ring_str

    def _gfan_ideal(self):
        """
        Return the ideal in gfan's notation.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: G._gfan_ideal()
            'Q[x, y, z]{x^2*y-z,y^2*z-x,x*z^2-y}'
        """
        try:
            return self.__gfan_ideal
        except AttributeError:
            ideal_gen_str = '{' + (str(self.__ideal.gens()).replace(' ', '').replace("'",""))[1:-1] + '}'
            ring_str = self._gfan_ring()
            self.__gfan_ideal = ring_str + ideal_gen_str
            return self.__gfan_ideal

    def ring(self):
        """
        Return the multivariate polynomial ring.

        EXAMPLES:
            sage: R.<x1,x2> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x1^3-x2,x2^3-x1-2]).groebner_fan()
            sage: gf.ring()
            Multivariate Polynomial Ring in x1, x2 over Rational Field
        """
        return self.__ring

    def _gfan_reduced_groebner_bases(self):
        """
        A string of the reduced Groebner bases of the ideal as output by gfan.

        EXAMPLES:
            sage: R.<a,b> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([a^3-b^2,b^2-a-1]).groebner_fan()
            sage: gf._gfan_reduced_groebner_bases()
            'Q[a,b]{{b^6-1+2*b^2-3*b^4,a+1-b^2},{b^2-1-a,a^3-1-a}}'
        """
        try:
            return self.__gfan_reduced_groebner_bases
        except AttributeError:
            B = self.gfan()
            B = B.replace('\n','')
            self.__gfan_reduced_groebner_bases = B
            return B

    def characteristic(self):
        """
        Return the characteristic of the base ring.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: i1 = ideal(x*z + 6*y*z - z^2, x*y + 6*x*z + y*z - z^2, y^2 + x*z + y*z)
            sage: gf = i1.groebner_fan()
            sage: gf.characteristic()
            0
        """
        return self.__ring.characteristic()

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
            sage: R3.<x,y,z> = PolynomialRing(GF(2477),3)
            sage: gf = R3.ideal([300*x^3-y,y^2-z,z^2-12]).groebner_fan()
            sage: gf.reduced_groebner_bases()
            [[z^2 - 12, y^2 - z, x^3 + 933*y],
            [-y^2 + z, y^4 - 12, x^3 + 933*y],
            [z^2 - 12, -300*x^3 + y, x^6 - 1062*z],
            [-828*x^6 + z, -300*x^3 + y, x^12 + 200]]
        """
        try:
            return self.__reduced_groebner_bases
        except AttributeError:
            G = self._gfan_reduced_groebner_bases()
            if G.find(']') != -1:
                G = G.split(']')[1]
            G = G.replace('{{','').replace('}}','').split('},{')
            S = self.__ring
            #print G
            #print [([f for f in G[i].split()], G[i]) for i in range(len(G))]
            X = [ReducedGroebnerBasis(self, [S(f) for f in G[i].split(',')], G[i]) for i in range(len(G))]
            self.__reduced_groebner_bases = X
            return X

    def _gfan_mod(self):
        """
        Return the extra options to the gfan command that are used by
        this object to account for working modulo a prime or in the
        presence of extra symmetries.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x^3-y,y^3-x-1]).groebner_fan()
            sage: gf._gfan_mod()
            ''
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
        """
        Returns the gfan output as a string given an input cmd; the
        default is to produce the list of reduced Groebner bases in
        gfan format.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x^3-y,y^3-x-1]).groebner_fan()
            sage: gf.gfan()
            'Q[x,y]\n{{\ny^9-1-y+3*y^3-3*y^6,\nx+1-y^3}\n,\n{\ny^3-1-x,\nx^3-y}\n,\n{\ny-x^3,\nx^9-1-x}\n}\n'
        """
        if I is None:
            I = self._gfan_ideal()
        # todo -- put something in here (?) when self.__symmetry isn't None...
        cmd += self._gfan_mod()
        s = gfan(I, cmd, verbose=self.__verbose, format=format)
        if s.strip() == '{':
            raise RuntimeError, "Error running gfan command %s on %s"%(cmd, self)
        return s

    def __iter__(self):
        """
        Returns an iterator for the reduced Groebner bases.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([x^3-y,y^3-x-1]).groebner_fan()
            sage: a = gf.__iter__()
            sage: a.next()
            [y^9 - 3*y^6 + 3*y^3 - y - 1, -y^3 + x + 1]
        """
        for x in self.reduced_groebner_bases():
            yield x

    def __getitem__(self, i):
        """
        Gets a reduced groebner basis

        EXAMPLES;
            sage: R4.<w1,w2,w3,w4> = PolynomialRing(QQ,4)
            sage: gf = R4.ideal([w1^2-w2,w2^3-1,2*w3-w4^2,w4^2-w1]).groebner_fan()
            sage: gf[0]
            [w4^12 - 1, -1/2*w4^2 + w3, -w4^4 + w2, -w4^2 + w1]
        """
        return self.reduced_groebner_bases()[i]

    def buchberger(self):
        """
        Computes and returns a lexicographic reduced Groebner basis
        for the ideal.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x - z^3, y^2 - x + x^2 - z^3*x]).groebner_fan()
            sage: G.buchberger()
            [-z^3 + y^2, -z^3 + x]
        """
        try:
            return self.__buchberger
        except AttributeError:
            B = self.gfan(cmd='buchberger')
            if B.find(']') != -1:
                B = B.split(']')[1]
            B = B.replace('}','').replace('{','')
            S = self.__ring
            B = [S(f) for f in B.split(',')]
            self.__buchberger = B
            return B

    def polyhedralfan(self):
        """
        Returns a polyhedral fan object corresponding to the reduced
        Groebner bases.

        EXAMPLES:
            sage: R3.<x,y,z> = PolynomialRing(QQ,3)
            sage: gf = R3.ideal([x^8-y^4,y^4-z^2,z^2-1]).groebner_fan()
            sage: pf = gf.polyhedralfan()
            sage: pf.rays()
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        """
        try:
            return self.__polyhedralfan
        except AttributeError:
            f = self.gfan(cmd='topolyhedralfan', I=self._gfan_reduced_groebner_bases())
            self.__polyhedralfan = PolyhedralFan(f)
            return PolyhedralFan(f)

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

    def render(self, file = None, larger=False, shift=0, rgbcolor = (0,0,0), polyfill = max_degree):
        """
        Render a Groebner fan as sage graphics or save as an xfig
        file.

        More precisely, the output is a drawing of the Groebner fan
        intersected with a triangle.  The corners of the triangle are
        (1,0,0) to the right, (0,1,0) to the left and (0,0,1) at the
        top.  If there are more than three variables in the ring we
        extend these coordinates with zeros.

        INPUT:
	    file  -- a filename if you prefer the output saved to a file.
	             This will be in xfig format.
            shift -- shift the positions of the variables in
                     the drawing.  For example, with shift=1,
                     the corners will be b (right), c (left),
                     and d (top).  The shifting is done modulo
                     the number of variables in the polynomial
                     ring.   The default is 0.
            larger -- bool (default: False); if True, make
                     the triangle larger so that the shape of
                     of the Groebner region appears.  Affects the xfig file
		     but probably not the sage graphics (?)
            rgbcolor -- This will not affect the saved xfig file, only the sage graphics
                     produced.
            polyfill -- Whether or not to fill the cones with a color determined by the highest degree in each reduced Groebner basis for that cone.


        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G.render()

            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: G = R.ideal([x^2*y - z, y^2*z - x, z^2*x - y]).groebner_fan()
            sage: G.render(larger=True)
        """
        cmd = 'render'
        if shift:
            cmd += ' --shiftVariables %s'%shift
        if larger:
            cmd += ' -L'
        s = self.gfan(cmd, I=self._gfan_reduced_groebner_bases().replace(' ',','), format=False)
        if file != None:
            open(file,'w').write(s)
        sp = s.split('\n')
        sp2 = []
        for x in sp[9:]:
            xs = x.split(' ')
            y = []
            if x[0:3] != '2 3' and len(xs) > 1:
                for q in xs:
                    if q != '':
                        y.append(q)
                sp2.append(y)
        sp3 = []
        for j in range(len(sp2)):
            temp = []
            for i in range(0,len(sp2[j])-1,2):
                temp.append([float(sp2[j][i])/1200.0, float(sp2[j][i+1])/1200.0])
            sp3.append(temp)
        r_lines = Graphics()
        for x in sp3:
            r_lines = r_lines + line(x, rgbcolor = rgbcolor)
        if polyfill:
            vals = [polyfill(q) for q in self.reduced_groebner_bases()]
            vmin = min(vals)
            vmax = max(vals)
            for index in range(len(sp3)):
                r_lines = r_lines + polygon(sp3[index], hue = .3 + .4*(vals[index]-vmin)/(vmax-vmin))
        return r_lines


    def _gfan_stats(self):
        """
        Return various statistics about this Groebner fan.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G._gfan_stats()
            {'Number of reduced Groebner bases': 3,
             'Maximal total degree of a Groebner basis': 4,
             'Dimension of homogeneity space': 0,
             'Number of variables': 2,
             'Minimal total degree of a Groebner basis': 2}
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

    def tropical_basis(self, check=True, verbose = False):
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
            sage: G.tropical_basis()
            [-3*x^2 + y^3, 2*x^2 - x - 2*y^3 - y + z^3, 3/4*x + y^3 + 3/4*y - 3/4*z^3]
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
        if B.find(']') != -1:
                B = B.split(']')[1]
        S = self.__ring
        B = B.replace('\n','')
        B = B.replace('{','').replace('}','').split(',')
        if verbose: print S, B
        X = [S(f) for f in B]
        self.__tropical_basis = X
        return X

    def interactive(self, *args, **kwds):
        """
        See the documentation for self[0].interative()
        This does not work with the notebook.

        EXAMPLES:
            sage: print "This is not easily doc-testable; please write a good one!"
            This is not easily doc-testable; please write a good one!
        """
        self[0].interactive(*args, **kwds)

    def tropical_intersection(self, ideal_arg = False, *args, **kwds):
	"""
	Returns information about the tropical intersection of the
        polynomials defining the ideal.

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: i1 = ideal(x*z + 6*y*z - z^2, x*y + 6*x*z + y*z - z^2, y^2 + x*z + y*z)
            sage: gf = i1.groebner_fan()
            sage: pf = gf.tropical_intersection()
            sage: pf.rays()
            [[-1, 0, 0]]
	"""
        try:
            return self.__tropical_intersection
        except AttributeError:
            f = self.gfan(cmd='tropicalintersection', I = self._gfan_ideal())
            pf = PolyhedralFan(f)
            self.__tropical_intersection = pf
            return pf



class ReducedGroebnerBasis(SageObject, list):
    def __init__(self, groebner_fan, gens, gfan_gens):
        """
        A class for representing reduced Groebner bases as produced by gfan.

        INPUT:
            groebner_fan -- a GroebnerFan object from an ideal
            gens -- the generators of the ideal
            gfan_gens -- the generators as a gfan string

        EXAMPLES:
            sage: R.<a,b> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([a^2-b^2,b-a-1]).groebner_fan()
            sage: from sage.rings.polynomial.groebner_fan import ReducedGroebnerBasis
            sage: ReducedGroebnerBasis(gf,gf[0],gf[0]._gfan_gens())
            [b - 1/2, a + 1/2]
        """
        self.__groebner_fan = groebner_fan
        list.__init__(self, gens)
        self.__gfan_gens = '{' + gfan_gens.replace(' ',',') + '}'
        self.__ring = groebner_fan._gfan_ring()

    def _repr_(self):
        """
        Returns the reduced Groebner basis as a string.

        EXAMPLES:
            sage: R.<z1,zz1> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([z1^2*zz1-1,zz1-2]).groebner_fan()
            sage: rgb1 = gf.reduced_groebner_bases()[0]
            sage: rgb1 # indirect doctest
            [zz1 - 2, z1^2 - 1/2]
        """
        return list.__repr__(self)

    def _gfan_gens(self):
        """
        Returns the reduced Groebner basis as a string in gfan format.

        EXAMPLES:
            sage: R.<z1,zz1> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([z1^2*zz1-1,zz1-2]).groebner_fan()
            sage: rgb1 = gf.reduced_groebner_bases()[0]
            sage: rgb1._gfan_gens()
            '{zz1-2,z1^2-1/2}'
        """
        return self.__gfan_gens

    def _gfan(self):
        """
        Returns a description of the Groebner fan this basis was derived from.

        EXAMPLES:
            sage: R.<z1,zz1> = PolynomialRing(QQ,2)
            sage: gf = R.ideal([z1^2*zz1-1,zz1-2]).groebner_fan()
            sage: rgb1 = gf.reduced_groebner_bases()[0]
            sage: rgb1._gfan()
            Groebner fan of the ideal:
            Ideal (z1^2*zz1 - 1, zz1 - 2) of Multivariate Polynomial Ring in z1, zz1 over Rational Field
        """
        return self.__groebner_fan

    def interactive(self, latex=False, flippable=False, wall=False,
                    inequalities=False, weight=False):
        """
        Do an interactive walk of the Groebner fan starting at this
        reduced Groebner basis.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: G = R.ideal([y^3 - x^2, y^2 - 13*x]).groebner_fan()
            sage: G[0].interactive()      # not tested
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
            sage: poly_cone = G[1].groebner_cone()
            sage: poly_cone.facets()
            [[-1, 2], [1, -1]]
            sage: [g.groebner_cone().facets() for g in G]
            [[[0, 1], [1, -2]], [[-1, 2], [1, -1]], [[-1, 1], [1, 0]]]
            sage: G[1].groebner_cone(restrict=True).facets()
            [[-1, 2], [1, -1]]
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
        c = gf.gfan(cmd=cmd, I=self.__ring + self.__gfan_gens)
        return PolyhedralCone(c)


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



