"""
This file defines a class for an isogeny class of an elliptic curve.

AUTHORS:

David Roe (2012-03-29) -- initial version.
"""

##############################################################################
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.structure.sage_object import SageObject
from sage.misc.lazy_attribute import lazy_attribute
import constructor
import sage.databases.cremona
from sage.rings.all import ZZ
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field

class IsogenyClass_EC(SageObject):
    """
    Isogeny class of an elliptic curve.

    The current implementation chooses a curve from each isomorphism
    class in the isogeny class, since over `\QQ` there is a unique
    reduced minimal model in each isomorphism class.

    EXAMPLES::


    """
    def __init__(self, E, label=None, empty=False):
        """
        Over `\QQ` we use curves since minimal models exist and there is a canonical choice of one.

        INPUT:

        - ``label`` -- string or None, a Cremona or LMFDB label, used
          in printing

        EXAMPLES::

            sage: cls = EllipticCurve('1011b1').isogeny_class()
            sage: print "\n".join([repr(E) for E in cls.curves])
            Elliptic Curve defined by y^2 + x*y = x^3 - 8*x - 9 over Rational Field
            Elliptic Curve defined by y^2 + x*y = x^3 - 23*x + 30 over Rational Field
        """
        self.E = E
        self._label = label
        if not empty:
            self._compute()

    def __len__(self):
        """
        The length is just the number of curves in the class.

        EXAMPLES::

            sage: E = EllipticCurve('15a')
            sage: len(E.isogeny_class()) # indirect doctest
            8
        """
        return len(self.curves)

    def __iter__(self):
        """
        Iterator over curves in the class.

        EXAMPLES::

            sage: E = EllipticCurve('15a')
            sage: all(C.conductor() == 15 for C in E.isogeny_class()) # indirect doctest
            True
        """
        return iter(self.curves)

    def __getitem__(self, i):
        """
        Gets the `i`th curve in the class.

        EXAMPLES::

            sage: E = EllipticCurve('990j1')
            sage: iso = E.isogeny_class(order="lmfdb") # orders lexicographically on a-invariants
            sage: iso[2] == E # indirect doctest
            True
        """
        return self.curves[i]

    def index(self, C):
        """
        Returns the index of a curve in this class.

        INPUT:

        - ``C`` -- an elliptic curve in this isogeny class.

        OUTPUT:

        - ``i`` -- an integer so that the ``i``th curve in the class
          is isomorphic to ``C``

        EXAMPLES::

            sage: E = EllipticCurve('990j1')
            sage: iso = E.isogeny_class(order="lmfdb") # orders lexicographically on a-invariants
            sage: iso.index(E.short_weierstrass_model())
            2
        """
        # This will need updating once we start talking about curves over more general number fields
        if not isinstance(C, EllipticCurve_rational_field):
            raise ValueError("x not in isogeny class")
        return self.curves.index(C.minimal_model())

    def __cmp__(self, other):
        """
        Returns 0 if self and other are the same isogeny class.

        If they are different, compares the sorted underlying lists of
        curves.

        Note that two isogeny classes with different orderings will
        compare as the same.  If you want to include the ordering,
        just compare the list of curves.

        EXAMPLES::

            sage: E = EllipticCurve('990j1')
            sage: EE = EllipticCurve('990j4')
            sage: E.isogeny_class() == EE.isogeny_class() # indirect doctest
            True
        """
        # This will need updating once we start talking about curves over more general number fields
        if isinstance(other, IsogenyClass_EC):
            return cmp(sorted(self.curves), sorted(other.curves))
        return cmp(type(self), type(other))

    def __hash__(self):
        """
        Hash is based on the a-invariants of the sorted list of
        minimal models.

        EXAMPLES::

            sage: E = EllipticCurve('990j1')
            sage: C = E.isogeny_class()
            sage: hash(C) == hash(tuple(sorted([curve.a_invariants() for curve in C.curves]))) # indirect doctest
            True
        """
        try:
            return self._hash
        except AttributeError:
            self._hash = hash(tuple(sorted([E.a_invariants() for E in self.curves])))
            return self._hash

    def _repr_(self):
        """
        The string representation depends on whether an LMFDB or
        Cremona label for the curve is known when this isogeny class
        is constructed.

        EXAMPLES:

        If the curve is constructed from an LMFDB label then that
        label is used::

            sage: E = EllipticCurve('462.f3')
            sage: E.isogeny_class() # indirect doctest
            Elliptic curve isogeny class 462.f

        If the curve is constructed from a Cremona label then that
        label is used::

            sage: E = EllipticCurve('990j1')
            sage: E.isogeny_class()
            Elliptic curve isogeny class 990j

        Otherwise the representation is determined from the first
        curve in the class::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.isogeny_class()
            Isogeny class of Elliptic Curve defined by y^2 + x*y = x^3 - x^2 + 4*x + 3 over Rational Field
        """
        if self._label:
            return "Elliptic curve isogeny class %s"%(self._label)
        else:
            return "Isogeny class of %r"%(self.curves[0])

    def __contains__(self, x):
        """
        INPUT:

        - ``x`` -- a Python object.

        OUTPUT:

        - boolean -- True iff ``x`` is an elliptic curve in this isogeny class.

        EXAMPLES::

            sage: cls = EllipticCurve('15a3').isogeny_class()
            sage: E = EllipticCurve('15a7'); E in cls
            True
            sage: E.short_weierstrass_model() in cls
            True
        """
        if not isinstance(x, EllipticCurve_rational_field):
            return False
        return x.minimal_model() in self.curves

    @cached_method
    def matrix(self, fill=True):
        """
        Returns the matrix whose entries give the minimal degrees of
        isogenies between curves in this class.

        INPUT:

        - ``fill`` -- boolean (default True).  If False then the
          matrix will contain only zeros and prime entries; if True it
          will fill in the other degrees.

        EXAMPLES::

            sage: isocls = EllipticCurve('15a3').isogeny_class()
            sage: isocls.matrix()
            [ 1  2  2  2  4  4  8  8]
            [ 2  1  4  4  8  8 16 16]
            [ 2  4  1  4  8  8 16 16]
            [ 2  4  4  1  2  2  4  4]
            [ 4  8  8  2  1  4  8  8]
            [ 4  8  8  2  4  1  2  2]
            [ 8 16 16  4  8  2  1  4]
            [ 8 16 16  4  8  2  4  1]
            sage: isocls.matrix(fill=False)
            [0 2 2 2 0 0 0 0]
            [2 0 0 0 0 0 0 0]
            [2 0 0 0 0 0 0 0]
            [2 0 0 0 2 2 0 0]
            [0 0 0 2 0 0 0 0]
            [0 0 0 2 0 0 2 2]
            [0 0 0 0 0 2 0 0]
            [0 0 0 0 0 2 0 0]
        """
        if self._mat is None:
            self._compute_matrix()
        mat = self._mat
        if fill and mat[0,0] == 0:
            from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix
            mat = fill_isogeny_matrix(mat)
        if not fill and mat[0,0] == 1:
            from sage.schemes.elliptic_curves.ell_curve_isogeny import unfill_isogeny_matrix
            mat = unfill_isogeny_matrix(mat)
        return mat

    @cached_method
    def isogenies(self, fill=False):
        """
        Returns a list of lists of isogenies and 0s, corresponding to the entries of :meth:`matrix`

        INPUT:

        - ``fill`` -- boolean (default False).  Whether to only return
          prime degree isogenies.  Currently only implemented for
          ``fill=False``.

        OUTPUT:

        - a list of lists, where the ``j``th entry of the ``i``th list
          is either zero or a prime degree isogeny from ``i``th curve
          in this class to the ``j``th curve.

        WARNING:

        - The domains and codmains of the isogenies will have the same
          Weierstrass equation as the curves in this class, but they
          may not be identical python objects in the current
          implementation.

        EXAMPLES::

            sage: isocls = EllipticCurve('15a3').isogeny_class()
            sage: f = isocls.isogenies()[0][1]; f
            Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 5*x + 2 over Rational Field to Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 80*x + 242 over Rational Field
            sage: f.domain() == isocls.curves[0] and f.codomain() == isocls.curves[1]
            True
        """
        if fill:
            raise NotImplementedError
        isogenies = self._maps
        if isogenies is None:
            self._compute_isogenies()
            isogenies = self._maps
        return isogenies

    @cached_method
    def graph(self):
        """
        Returns a graph whose vertices correspond to curves in this class, and whose edges correspond to prime degree isogenies.

        .. note:

            There are only finitely many possible isogeny graphs for
            curves over `\QQ` [M78].  This function tries to lay out the graph
            nicely by special casing each isogeny graph.

        .. note:

            The vertices are labeled 1 to n rather than 0 to n-1 to
            correspond to LMFDB and Cremona labels.

        EXAMPLES::

            sage: isocls = EllipticCurve('15a3').isogeny_class()
            sage: G = isocls.graph()
            sage: sorted(G._pos.items())
            [(1, [-0.8660254, 0.5]), (2, [-0.8660254, 1.5]), (3, [-1.7320508, 0]), (4, [0, 0]), (5, [0, -1]), (6, [0.8660254, 0.5]), (7, [0.8660254, 1.5]), (8, [1.7320508, 0])]

        REFERENCES:

        .. [M78] B. Mazur.  Rational Isogenies of Prime Degree.
          *Inventiones mathematicae* 44,129-162 (1978).
        """
        from sage.graphs.graph import Graph
        M = self.matrix(fill = False)
        n = M.nrows() # = M.ncols()
        G = Graph(M, format='weighted_adjacency_matrix')
        N = self.matrix(fill = True)
        D = dict([(v,self.curves[v]) for v in G.vertices()])
        # The maximum degree classifies the shape of the isogeny
        # graph, though the number of vertices is often enough.
        # This only holds over Q, so this code will need to change
        # once other isogeny classes are implemented.
        if n == 1:
            # one vertex
            pass
        elif n == 2:
            # one edge, two vertices.  We align horizontally and put
            # the lower number on the left vertex.
            G.set_pos(pos={0:[-0.5,0],1:[0.5,0]})
        else:
            maxdegree = max(max(N))
            if n == 3:
                # o--o--o
                centervert = [i for i in range(3) if max(N.row(i)) < maxdegree][0]
                other = [i for i in range(3) if i != centervert]
                G.set_pos(pos={centervert:[0,0],other[0]:[-1,0],other[1]:[1,0]})
            elif maxdegree == 4:
                # o--o<8
                centervert = [i for i in range(4) if max(N.row(i)) < maxdegree][0]
                other = [i for i in range(4) if i != centervert]
                G.set_pos(pos={centervert:[0,0],other[0]:[0,1],other[1]:[-0.8660254,-0.5],other[2]:[0.8660254,-0.5]})
            elif maxdegree == 27:
                # o--o--o--o
                centers = [i for i in range(4) if list(N.row(i)).count(3) == 2]
                left = [j for j in range(4) if N[centers[0],j] == 3 and j not in centers][0]
                right = [j for j in range(4) if N[centers[1],j] == 3 and j not in centers][0]
                G.set_pos(pos={left:[-1.5,0],centers[0]:[-0.5,0],centers[1]:[0.5,0],right:[1.5,0]})
            elif n == 4:
                # square
                opp = [i for i in range(1,4) if not N[0,i].is_prime()][0]
                other = [i for i in range(1,4) if i != opp]
                G.set_pos(pos={0:[1,1],other[0]:[-1,1],opp:[-1,-1],other[1]:[1,-1]})
            elif maxdegree == 8:
                # 8>o--o<8
                centers = [i for i in range(6) if list(N.row(i)).count(2) == 3]
                left = [j for j in range(6) if N[centers[0],j] == 2 and j not in centers]
                right = [j for j in range(6) if N[centers[1],j] == 2 and j not in centers]
                G.set_pos(pos={centers[0]:[-0.5,0],left[0]:[-1,0.8660254],left[1]:[-1,-0.8660254],centers[1]:[0.5,0],right[0]:[1,0.8660254],right[1]:[1,-0.8660254]})
            elif maxdegree == 18:
                # two squares joined on an edge
                centers = [i for i in range(6) if list(N.row(i)).count(3) == 2]
                top = [j for j in range(6) if N[centers[0],j] == 3]
                bl = [j for j in range(6) if N[top[0],j] == 2][0]
                br = [j for j in range(6) if N[top[1],j] == 2][0]
                G.set_pos(pos={centers[0]:[0,0.5],centers[1]:[0,-0.5],top[0]:[-1,0.5],top[1]:[1,0.5],bl:[-1,-0.5],br:[1,-0.5]})
            elif maxdegree == 16:
                # tree from bottom, 3 regular except for the leaves.
                centers = [i for i in range(8) if list(N.row(i)).count(2) == 3]
                center = [i for i in centers if len([j for j in centers if N[i,j] == 2]) == 2][0]
                centers.remove(center)
                bottom = [j for j in range(8) if N[center,j] == 2 and j not in centers][0]
                left = [j for j in range(8) if N[centers[0],j] == 2 and j != center]
                right = [j for j in range(8) if N[centers[1],j] == 2 and j != center]
                G.set_pos(pos={center:[0,0],bottom:[0,-1],centers[0]:[-0.8660254,0.5],centers[1]:[0.8660254,0.5],left[0]:[-0.8660254,1.5],right[0]:[0.8660254,1.5],left[1]:[-1.7320508,0],right[1]:[1.7320508,0]})
            elif maxdegree == 12:
                # tent
                centers = [i for i in range(8) if list(N.row(i)).count(2) == 3]
                left = [j for j in range(8) if N[centers[0],j] == 2]
                right = []
                for i in range(3):
                    right.append([j for j in range(8) if N[centers[1],j] == 2 and N[left[i],j] == 3][0])
                G.set_pos(pos={centers[0]:[-0.75,0],centers[1]:[0.75,0],left[0]:[-0.75,1],right[0]:[0.75,1],left[1]:[-1.25,-0.75],right[1]:[0.25,-0.75],left[2]:[-0.25,-0.25],right[2]:[1.25,-0.25]})
        G.set_vertices(D)
        G.relabel(range(1,n+1))
        return G

    @cached_method
    def reorder(self, order):
        """
        Return a new isogeny class with the curves reordered.

        INPUT:

        - ``order`` -- None, a string or an iterable over all curves
          in this class.  See
          :meth:`sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.isogeny_class`
          for more details.

        OUTPUT:

        - Another :class:`IsogenyClass_EC` with the curves reordered
          (and matrices and maps changed as appropriate)

        EXAMPLES::

            sage: isocls = EllipticCurve('15a1').isogeny_class()
            sage: print "\n".join([repr(C) for C in isocls.curves])
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 5*x + 2 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + 35*x - 28 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 135*x - 660 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 80*x + 242 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 110*x - 880 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 2160*x - 39540 over Rational Field
            sage: isocls2 = isocls.reorder('lmfdb')
            sage: print "\n".join([repr(C) for C in isocls2.curves])
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 2160*x - 39540 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 135*x - 660 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 110*x - 880 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 80*x + 242 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 5*x + 2 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 over Rational Field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + 35*x - 28 over Rational Field
        """
        if order is None or isinstance(order, basestring) and order == self._algorithm:
            return self
        if isinstance(order, basestring):
            if order == "lmfdb":
                reordered_curves = sorted(self.curves, key = lambda E: E.a_invariants())
            else:
                reordered_curves = list(self.E.isogeny_class(algorithm=order))
        elif isinstance(order, (list, tuple, IsogenyClass_EC)):
            reordered_curves = list(order)
            if len(reordered_curves) != len(self.curves):
                raise ValueError("Incorrect length")
        else:
            raise TypeError("order parameter should be a string, list of curves or isogeny class")
        need_perm = self._mat is not None
        cpy = self.copy()
        curves = []
        perm = []
        for E in reordered_curves:
            try:
                j = self.curves.index(E)
            except ValueError:
                try:
                    j = self.curves.index(E.minimal_model())
                except ValueError:
                    raise ValueError("order does not yield a permutation of curves")
            curves.append(self.curves[j])
            if need_perm: perm.append(j+1)
        cpy.curves = tuple(curves)
        if need_perm:
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            perm = SymmetricGroup(len(self.curves))(perm)
            cpy._mat = perm.matrix() * self._mat * (~perm).matrix()
            if self._maps is not None:
                n = len(self._maps)
                cpy._maps = [self._maps[perm(i+1)-1] for i in range(n)]
                for i in range(n):
                    cpy._maps[i] = [cpy._maps[i][perm(j+1)-1] for j in range(n)]
        else:
            cpy._mat = None
            cpy._maps = None
        return cpy

    @abstract_method
    def _compute(self):
        """
        This function is called during initialization and should fill
        in ``self.curves``, and possibly ``self._mat`` and
        ``self._maps`` as well, depending on the algorithm.

        .. note:

            This is currently redundant; it used to be required when mwrank was used.
        """
        pass

    @abstract_method
    def _compute_matrix(self):
        """
        For algorithms that don't compute the isogeny matrix at first,
        this function should fill in ``self._mat``.

        EXAMPLES::

            sage: EllipticCurve('11a1').isogeny_class('database').matrix() # indirect doctest
            [ 1  5  5]
            [ 5  1 25]
            [ 5 25  1]
        """
        pass

    @abstract_method
    def _compute_isogenies(self):
        """
        For algorithms that don't compute the isogenies at first, this
        function should fill in ``self._maps``.

        .. note:

            This is currently redundant; it used to be required when mwrank was used.
        """
        pass

class IsogenyClass_EC_Rational(IsogenyClass_EC):
    """
    Isogeny classes for elliptic curves over `\QQ`.
    """
    def __init__(self, E, algorithm="sage", label=None, empty=False):
        """
        INPUT:

        - ``E`` -- an elliptic curve over Q.

        - ``algorithm`` -- a string (default "sage").  One of the
          following:

          - "sage" -- Use sage's implementation to compute the curves,
            matrix and isogenies

          - "database" -- Use the Cremona database (only works if the
            curve is in the database)

        - ``label`` -- a string, the label of this isogeny class
          (e.g. '15a' or '37.b').  Used in printing.

        - ``empty`` -- don't compute the curves right now (used when reordering)

        EXAMPLES::

            sage: isocls = EllipticCurve('389a1').isogeny_class(); isocls
            Elliptic curve isogeny class 389a
            sage: E = EllipticCurve([0, 0, 0, 0, 1001]) # conductor 108216108
            sage: E.isogeny_class(order='database')
            Traceback (most recent call last):
            ...
            RuntimeError: unable to to find Elliptic Curve defined by y^2 = x^3 + 1001 over Rational Field in the database
            sage: TestSuite(isocls).run()
        """
        self._algorithm = algorithm
        IsogenyClass_EC.__init__(self, E, label=label, empty=empty)

    def copy(self):
        """
        Returns a copy (mostly used in reordering).

        EXAMPLES::

            sage: isocls = EllipticCurve('389a1').isogeny_class()
            sage: isocls2 = isocls.copy()
            sage: isocls is isocls2
            False
            sage: isocls == isocls2
            True
        """
        ans = IsogenyClass_EC_Rational(self.E, self._algorithm, self._label, empty=True)
        # The following isn't needed internally, but it will keep
        # things from breaking if this is used for something other
        # than reordering.
        ans.curves = self.curves
        ans._mat = None
        ans._maps = None
        return ans

    def _compute(self):
        """
        Computes the list of curves, and possibly the matrix and
        prime-degree isogenies (depending on the algorithm selected).

        EXAMPLES::

            sage: isocls = EllipticCurve('48a1').isogeny_class('sage').copy()
            sage: isocls._mat
            sage: isocls._compute(); isocls._mat
            [0 2 2 2 0 0]
            [2 0 0 0 2 2]
            [2 0 0 0 0 0]
            [2 0 0 0 0 0]
            [0 2 0 0 0 0]
            [0 2 0 0 0 0]
        """
        algorithm = self._algorithm
        from sage.schemes.elliptic_curves.ell_curve_isogeny import fill_isogeny_matrix, unfill_isogeny_matrix
        from sage.matrix.all import MatrixSpace
        self._maps = None
        if algorithm == "database":
            try:
                label = self.E.cremona_label(space=False)
            except RuntimeError:
                raise RuntimeError("unable to to find %s in the database"%self.E)
            db = sage.databases.cremona.CremonaDatabase()
            curves = db.isogeny_class(label)
            if len(curves) == 0:
                raise RuntimeError("unable to to find %s in the database"%self.E)
            # All curves will have the same conductor and isogeny class,
            # and there are are most 8 of them, so lexicographic sorting is okay.
            self.curves = tuple(sorted(curves, key = lambda E: E.cremona_label()))
            self._mat = None
        elif algorithm == "sage":
            curves = [self.E.minimal_model()]
            ijl_triples = []
            l_list = None
            i = 0
            while i<len(curves):
                E = curves[i]
                isogs = E.isogenies_prime_degree(l_list)
                for phi in isogs:
                    Edash = phi.codomain()
                    l = phi.degree()
                    # look to see if Edash is new.  Note that the
                    # curves returned by isogenies_prime_degree() are
                    # standard minimal models, so it suffices to check
                    # equality rather than isomorphism here.
                    try:
                        j = curves.index(Edash)
                    except ValueError:
                        j = len(curves)
                        curves.append(Edash)
                    ijl_triples.append((i,j,l,phi))
                if l_list is None:
                    l_list = [l for l in set([ZZ(f.degree()) for f in isogs])]
                i = i+1
            self.curves = tuple(curves)
            ncurves = len(curves)
            self._mat = MatrixSpace(ZZ,ncurves)(0)
            self._maps = [[0]*ncurves for i in range(ncurves)]
            for i,j,l,phi in ijl_triples:
                self._mat[i,j] = l
                self._maps[i][j]=phi
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm

    def _compute_matrix(self):
        """
        Computes the matrix, assuming that the list of curves is computed.

        EXAMPLES::

            sage: isocls = EllipticCurve('1225h1').isogeny_class('database')
            sage: isocls._mat
            sage: isocls._compute_matrix(); isocls._mat
            [ 0 37]
            [37  0]
        """
        self._mat = self.E.isogeny_class(order=self.curves)._mat

    def _compute_isogenies(self):
        """
        EXAMPLES::

            sage: E = EllipticCurve('15a1')
            sage: isocls = E.isogeny_class()
            sage: maps = isocls.isogenies() # indirect doctest
            sage: f = maps[0][1]
            sage: f.domain() == isocls[0] and f.codomain() == isocls[1]
            True
        """
        recomputed = self.E.isogeny_class(order=self.curves)
        self._mat = recomputed._mat
        # The domains and codomains here will be equal, but not the same Python object.
        self._maps = recomputed._maps
