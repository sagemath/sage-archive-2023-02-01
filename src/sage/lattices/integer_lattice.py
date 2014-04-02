# -*- coding: utf-8 -*-
"""
Discrete Subgroups of `\ZZ^n`.

AUTHORS:

- Martin Albrecht (2014-03): initial version

- Jan Pöschko (2012-08): some code in this module was taken from Jan Pöschko's
  2012 GSoC project

TESTS::

    sage: L = RealLattice(random_matrix(ZZ, 10, 10))
    sage: TestSuite(L).run()
"""
 
from copy import copy
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.lattices.real_lattice import RealLattice
from sage.libs.pari.pari_instance import pari
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.method_decorator import MethodDecorator
from sage.modules.free_module_element import vector
from sage.modules.vector_integer_dense import Vector_integer_dense
from sage.structure.parent import Parent

class IntegerLattice(RealLattice):
    """
    This class represents lattices over the integers.

    Lattices are discrete subgroups of $\\RR^n$. However, here we restrict out
    attention to subgroups of $\\ZZ^n$, which are common in many applications.

    Lattices are represented by a basis which is not unique.  In particular, the
    basis representing our lattice might change during the life-time of an
    object as the basis is improved when new information such as short vectors
    become available.

    EXAMPLE::

         sage: L = RealLattice(sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True)); L
         Lattice of degree 10 and rank 10 over Integer Ring
         Basis matrix:
         [ 0  1  2  0  1  2 -1  0 -1 -1]
         [ 0  1  0 -3  0  0  0  0  3 -1]
         [ 1  1 -1  0 -3  0  0  1  2 -2]
         [-1  2 -1 -1  2 -2  1 -1  0 -1]
         [ 1  0 -4  2  0  1 -2 -1  0  0]
         [ 2  3  0  1  1  0 -2  3  0  0]
         [-2 -3 -2  0  0  1 -1  1  3 -2]
         [-3  0 -1  0 -2 -1 -2  1 -1  1]
         [ 1  4 -1  1  2  2  1  0  3  1]
         [-1 -1  0 -3 -1  2  2  3 -1  0]
         sage: L.shortest_vector()
         (0, 1, 2, 0, 1, 2, -1, 0, -1, -1)
    """
    def __init__(self, basis, lll_reduce=True):
        """
        Construct a new lattice.

        INPUT:

        - ``basis``
          - a list of vectors or
          - a matrix over the integers
          - an element of an absolute order

        - ``lll_reduce`` -- (default: ``True``) run LLL reduction on the basis
          on construction

        EXAMPLES::

            sage: RealLattice([[1,0,-2],[0,2,5], [0,0,7]])
            Lattice of degree 3 and rank 3 over Integer Ring
            Basis matrix:
            [ 1  0 -2]
            [ 1 -2  0]
            [ 2  2  1]

            sage: RealLattice(random_matrix(ZZ, 5, 5, x=-2^20, y=2^20))
            Lattice of degree 5 and rank 5 over Integer Ring
            Basis matrix:
            [  -7945 -381123   85872 -225065   12924]
            [-158254  120252  189195 -262144 -345323]
            [ 232388  -49556  306585  -31340  401528]
            [-353460  213748  310673  158140  172810]
            [-287787  333937 -145713 -482137  186529]

            sage: K.<a> = NumberField(x^8+1)
            sage: O = K.ring_of_integers()
            sage: f = O.random_element(); f
            a^7 - a^6 + 4*a^5 - a^4 + a^3 + 1

            sage: RealLattice(f)
            Lattice of degree 8 and rank 8 over Integer Ring
            Basis matrix:
            [ 0  1  0  1  0  3  3  0]
            [ 1  0  0  1 -1  4 -1  1]
            [ 0  0  1  0  1  0  3  3]
            [-4  1 -1  1  0  0  1 -1]
            [ 1 -3  0  0  0  3  0 -2]
            [ 0 -1  1 -4  1 -1  1  0]
            [ 2  0 -3 -1  0 -3  0  0]
            [-1  0 -1  0 -3 -3  0  0]

        .. note::

            If your input is already LLL reduced and you do not want to pay the
            cost for running LLL again, set ``lll_reduce=False`` and set
            :attr:`IntegerLattice._is_LLL_reduced` to ``True``.
        """
        self._basis = matrix(ZZ, basis)

        if not lll_reduce and self._basis.nrows() > self._basis.ncols():
            raise ValueError("Basis nrows (%d) > ncols (%d), so not a basis."%(self._basis.nrows(), self._basis.ncols()))

        if not lll_reduce and self.rank() != self._basis.nrows():
            raise ValueError("Input basis has rank %d but expecting rank %d."%(self.rank(), self._basis.nrows()))

        self._is_LLL_reduced = False

        Parent.__init__(self, facade=ZZ**self.degree(), category=CommutativeAdditiveGroups())

        if lll_reduce:
            self.LLL()

    @cached_method
    def _inverse_of_basis(self):
        """
        We cache ``~self.basis`` as it is used to test if a vector is in this
        lattice.

        EXAMPLE::

            sage: L = RealLattice(Matrix(ZZ, 4, 4, [[2,1,0,-1],[1,7,3,1],[0,3,54,3],[-1,1,3,673]])); L
            Lattice of degree 4 and rank 4 over Integer Ring
            Basis matrix:
            [   2    1    0   -1]
            [  -1    6    3    2]
            [   6  -20   42   -6]
            [ 291 -107   -9  473]

            sage: L._inverse_of_basis()
            [  8807/25538   -841/25538     65/25538     23/25538]
            [  4843/38307   4325/38307   -623/76614     -4/12769]
            [-3523/229842 14983/229842  4399/229842     -5/76614]
            [ -7043/38307   1802/38307   -233/76614     19/12769]

        TESTS::

            sage: L = RealLattice(random_matrix(ZZ, 5, 10))
            sage: L._inverse_of_basis()
            Traceback (most recent call last):
            ...
            ValueError: Basis does not have full rank.

        """
        if self.rank() == self.degree():
            return ~self.basis
        else:
            raise ValueError("Basis does not have full rank.")

    def _element_constructor_(self, x):
        """
        Make sure x is a valid member of self, and return the constructed
        element.

        EXAMPLE::

            sage: L = RealLattice(Matrix(ZZ, 4, 4, [[2,1,0,-1],[1,7,3,1],[0,3,54,3],[-1,1,3,673]])); L
            Lattice of degree 4 and rank 4 over Integer Ring
            Basis matrix:
            [   2    1    0   -1]
            [  -1    6    3    2]
            [   6  -20   42   -6]
            [ 291 -107   -9  473]

            sage: L(0) # indirect doctest
            (0, 0, 0, 0)

            sage: vector(ZZ,4,(1,0,0,1)) in L
            False

            sage: L.shortest_vector() in L
            True

        TESTS::

            sage: L = RealLattice(random_matrix(ZZ, 5, 10))
            sage: L.shortest_vector() in L
            True

        """
        x = self.facade_for()[0](x)

        if self.degree() == self.rank():
            y = x * self._inverse_of_basis()
        else:
            y = self._basis.solve_left(x)

        try:
            y = y.change_ring(ZZ)
        except TypeError:
            raise ValueError("x is not in this lattice.")

        x.set_immutable()
        return x

    def _an_element_(self):
        return self.facade_for()[0].zero()

    def LLL(self, *args, **kwds):
        """
        Run LLL on the current basis.

        A lattice basis `(b_1, b_2, ..., b_d)` is `(δ,η)`-LLL-reduced if the two
        following conditions hold:

        -  For any `i>j`, we have `|μ_{i, j}| ≦ η`,

        -  For any `i<d`, we have `δ|b_i^*|^2 ≦ |b_{i+1}^* + μ_{i+1, i} b_i^*|^2`,

        where `μ_{i,j} = <b_i, b_j^*>/<b_j^*,b_j^*>` and `b_i^*` is the `i`-th vector
        of the Gram-Schmidt orthogonalisation of `(b_1, b_2, …, b_d)`.

        The default reduction parameters are `δ=3/4` and `η=0.501`.

        The parameters `δ` and `η` must satisfy: `0.25 < δ ≦ 1.0` and
        `0.5 ≦ η < √δ`. Polynomial time complexity is only guaranteed
        for `δ < 1`.

        INPUT:

        - ``*args`` - passed through to :func:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL`
        - ``*kwds`` - passed through to :func:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL`

        EXAMPLE::

            sage: A = random_matrix(ZZ, 10, 10, x=-2000, y=2000)
            sage: L = RealLattice(A, lll_reduce=False); L
            Lattice of degree 10 and rank 10 over Integer Ring
            Basis matrix:
            [ -645 -1037 -1775 -1619  1721 -1434  1766  1701  1669  1534]
            [ 1303   960  1998 -1838  1683 -1332   149   327  -849 -1562]
            [-1113 -1366  1379   669    54  1214 -1750  -605 -1566  1626]
            [-1367  1651   926  1731  -913   627   669 -1437  -132  1712]
            [ -549  1327 -1353    68  1479 -1803  -456  1090  -606  -317]
            [ -221 -1920 -1361  1695  1139   111 -1792  1925  -656  1992]
            [-1934   -29    88   890  1859  1820 -1912 -1614 -1724  1606]
            [ -590 -1380  1768   774   656   760  -746  -849  1977 -1576]
            [  312  -242 -1732  1594  -439 -1069   458 -1195  1715    35]
            [  391  1229 -1815   607  -413  -860  1408  1656  1651  -628]
            sage: min(v.norm().n() for v in L.basis)
            3346.57...

            sage: L.LLL()
            sage: L
            Lattice of degree 10 and rank 10 over Integer Ring
            Basis matrix:
            [ -888    53  -274   243   -19   431   710   -83   928   347]
            [  448  -330   370  -511   242  -584    -8  1220   502   183]
            [ -524  -460   402  1338  -247  -279 -1038   -28  -159  -794]
            [  166  -190  -162  1033  -340   -77 -1052  1134  -843   651]
            [  -47 -1394  1076  -132   854  -151   297  -396  -580  -220]
            [-1064   373  -706   601  -587 -1394   424   796   -22  -133]
            [-1126   398   565 -1418  -446  -890  -237  -378   252   247]
            [ -339   799   295   800   425  -605  -730 -1160   808   666]
            [  755 -1206  -918  -192 -1063   -37  -525   -75   338   400]
            [  382  -199 -1839  -482   984   -15  -695   136   682   563]
            sage: L.basis[0].norm().n()
            1613.74...
        """
        self._basis = [v for v in self._basis.LLL(*args, **kwds) if v]
        self._basis = matrix(ZZ, len(self._basis), len(self._basis[0]), self._basis)
        self.gram_matrix.clear_cache()
        self._is_LLL_reduced = True

    def BKZ(self, *args, **kwds):
        """
        Run Block Korkine-Zolotareff reduction on the current basis.

        .. note::

            If ``block_size == L.rank()`` where ``L`` is this latice, then this
            function performs Hermite-Korkine-Zolotareff (HKZ) reduction.

        INPUT:

        - ``*args`` - passed through to :func:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.BKZ`
        - ``*kwds`` - passed through to :func:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.BKZ`

        EXAMPLE::

            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=100, q=2^60, seed=42)
            sage: L = RealLattice(A, lll_reduce=False)
            sage: min(v.norm().n() for v in L.basis)
            4.17330740711759e15

            sage: L.LLL()
            sage: min(v.norm().n() for v in L.basis)
            5.09901951359278

            sage: L.BKZ(block_size=10)
            sage: min(v.norm().n() for v in L.basis)
            4.12310562561766

        """
        self._basis = [v for v in self._basis.BKZ(*args, **kwds) if v]
        self._basis = matrix(ZZ, len(self._basis), len(self._basis[0]), self._basis)
        self.gram_matrix.clear_cache()
        self._is_LLL_reduced = True

    def HKZ(self, *args, **kwds):
        """
        Hermite-Korkine-Zolotarev (HKZ) reduce the basis.

        A basis `B` of a lattice `L`, with orthogonalized basis `B^*` such
        that `B = M \cdot B^*` is HKZ reduced, if and only if, the following properties are
        satisfied:
        
        #. The basis `B` is size-reduced, i.e., all off-diagonal coefficients of
           `M` satisfy `|μ_{i,j}| ≦ 1/2`

        #. The vector `b_1` realizes the first minimum `λ_1(L)`.

        #. The projection of the vectors `b_2,…,b_r` orthogonally to `b_1` form
           an HKZ reduced basis.
       
        .. note::

            This is realised by calling
            :func:`sage.lattices.integer_lattice.IntegerLattice.BKZ` with
            ``block_size == self.rank()``.
        
        INPUT:

        - ``*args`` - passed through to :func:`sage.lattices.integer_lattice.IntegerLattice.BKZ`
        - ``*kwds`` - passed through to :func:`sage.lattices.integer_lattice.IntegerLattice.BKZ`


        EXAMPLE::

            sage: L = sage.crypto.gen_lattice(type='random', n=1, m=40, q=2^60, seed=42, lattice=True)
            sage: L.HKZ()
            sage: L.basis[0]
            (-1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, -1, 1, 0, 1, 1, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0)        
        """
        self.BKZ(block_size=self.rank())
        
        
    @cached_method
    def gram_matrix(self):
        """
        Return $B \cdot B^T$ where $B$ is the current basis.

        EXAMPLE::

            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=100, q=2^60, seed=42)
            sage: L = RealLattice(A, lll_reduce=False)
            sage: L.basis * L.basis.T == L.gram_matrix()
            True

            sage: L.LLL()
            sage: L.basis * L.basis.T == L.gram_matrix()
            True
        """
        return self._basis * self._basis.T

    def degree(self):
        """
        Return the number of columns f the basis matrix.

        EXAMPLE::

            sage: K.<a> = NumberField(x^8+1)
            sage: O = K.ring_of_integers()
            sage: RealLattice(O.random_element()).degree()
            8
        """
        return self._basis.ncols()

    @cached_method
    def rank(self):
        """
        Return the rank of the basis matrix.

        EXAMPLE::

            sage: K.<a> = NumberField(x^8+1)
            sage: O = K.ring_of_integers()
            sage: RealLattice(O.random_element()).rank()
            8
        """
        return self._basis.rank()

    @cached_method
    def volume(self):
        """
        Return $vol(L)$ which is $\\sqrt{\det(B \cdot B^T)}$ for any basis $B$.

        EXAMPLE::

            sage: L = sage.crypto.gen_lattice(m=10, seed=42, lattice=True)
            sage: L.volume()
            14641
        """
        if self.rank() == self.degree():
            return abs(self._basis.determinant())
        else:
            return self.gram_matrix().determinant().sqrt()

    @cached_method
    def discriminant(self):
        """
        Return $|\det(G)|$, i.e. the absolute value of the determinant of the
        gram matrix $B \cdot B^T$ for any basis $B$.

        EXAMPLE::

            sage: L = sage.crypto.gen_lattice(m=10, seed=42, lattice=True)
            sage: L.discriminant()
            214358881
        """
        return abs(self.gram_matrix().determinant())

    @cached_method
    def is_unimodular(self):
        """
        Return True if this lattice is unimodular.

        EXAMPLES::

            sage: L = RealLattice([[1, 0], [0, 1]])
            sage: L.is_unimodular()
            True
            sage: RealLattice([[2, 0], [0, 3]]).is_unimodular()
            False
        """
        return self.volume() == 1

    def base_ring(self):
        """
        Return $\\ZZ$

        TESTS::

             sage: RealLattice(random_matrix(ZZ, 2, 2)).base_ring()
             Integer Ring

        """
        return ZZ

    def _repr_(self):
        """

        TESTS::

            sage: RealLattice([[1,0,-2],[0,2,5], [0,0,7]]) #indirect doctest
            Lattice of degree 3 and rank 3 over Integer Ring
            Basis matrix:
            [ 1  0 -2]
            [ 1 -2  0]
            [ 2  2  1]
        """
        return "Lattice of degree %s and rank %s over %s\nBasis matrix:\n%s"%(
            self.degree(), self.rank(), self.base_ring(), self.basis)

    @property
    def basis(self):
        """
        Return the basis $B$.

        .. note::

            During the life-time of this lattice object the basis $B$ returned
            might change, as improved lattice bases are computed.

        TESTS::

            sage: RealLattice([[1,0,-2],[0,2,5], [0,0,7]]).basis
            [ 1  0 -2]
            [ 1 -2  0]
            [ 2  2  1]
        """
        return copy(self._basis)

    @cached_method
    def shortest_vector(self, update_basis=True, algorithm="fplll", *args, **kwds):
        """
        Return a shortest vector.

        INPUT:

        - ``update_basis`` -- (default: True) set this flag if the found vector
          should be used to improve the basis

        - ``algorithm`` -- (default: "fplll") either "fplll" or "pari"

        - ``*args`` - passed through to underlying implementation

        - ``*kwds`` - passed through to underlying implementation

        EXAMPLE::

            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=30, q=2^40, seed=42)
            sage: L = RealLattice(A, lll_reduce=False)
            sage: min(v.norm().n() for v in L.basis)
            6.03890756700000e10

            sage: L.shortest_vector().norm().n()
            3.74165738677394

            sage: L = RealLattice(A, lll_reduce=False)
            sage: min(v.norm().n() for v in L.basis)
            6.03890756700000e10

            sage: L.shortest_vector(algorithm="pari").norm().n()
            3.74165738677394

        """
        if algorithm == "pari":
            if not self._is_LLL_reduced:
                self.LLL()
            qf = self.gram_matrix()

            count, length, vectors = pari(qf).qfminim(0, None)
            v = vectors.python().columns()[0]
            w = v*self._basis
        elif algorithm == "fplll":
            from sage.libs.fplll.fplll import FP_LLL
            L = FP_LLL(self._basis)
            w = L.shortest_vector(*args, **kwds)
        else:
            raise ValueError("Algorithm '%s' unknown."%algorithm)

        if update_basis:
            self.update_basis(w)
        return w

    def update_basis(self, w):
        """
        Inject the vector w and run LLL to update the basis.

        EXAMPLE::

            sage: A = sage.crypto.gen_lattice(type='random', n=1, m=30, q=2^40, seed=42)
            sage: L = RealLattice(A)
            sage: B = L.basis
            sage: v = L.shortest_vector(update_basis=False)
            sage: L.update_basis(v)
            sage: bool(L.basis[0].norm() < B[0].norm())
            True

        """
        w = matrix(ZZ, w)

        L = w.stack(self._basis).LLL()
        assert(L[0] == 0)
        self._basis = L.matrix_from_rows(range(1,L.nrows()))
        self.gram_matrix.clear_cache()
        self._is_LLL_reduced = True

    @cached_method
    def hermite_form_basis(self):
        """
        Return Hermite Normal form of basis, i.e. a unique representation for this lattice.

        EXAMPLES::

            sage: RealLattice(random_matrix(ZZ, 10, 11)).hermite_form_basis()
            [           1            0            0            0            0            0            0            0            0 107070676021    440848398]
            [           0            1            0            0            0            0            0            0            0  96760618325    398398190]
            [           0            0            1            0            0            0            0            1            0  12932839669     53249143]
            [           0            0            0            1            0            0            0            0            0 105468971463    434253605]
            [           0            0            0            0            1            0            0            0            0  69813072820    287445476]
            [           0            0            0            0            0            1            0            0            0 100089000492    412102334]
            [           0            0            0            0            0            0            1            2            0  73712648156    303501425]
            [           0            0            0            0            0            0            0            3            0  54135613288    222895745]
            [           0            0            0            0            0            0            0            0            1  37920244961    156131257]
            [           0            0            0            0            0            0            0            0            0 140081789582    576766999]

        """
        return self.basis.hermite_form()

    def __eq__(self, other):
        """
        TESTS::

            sage: A = random_matrix(ZZ, 20, 10)
            sage: L0 = RealLattice(A)
            sage: L1 = RealLattice(L0.basis)
            sage: L0 == L1
            True
            sage: L0 == RealLattice(random_matrix(ZZ, 10, 10))
            False
        """
        if not isinstance(other, IntegerLattice):
            return False
        return self.hermite_form_basis() == other.hermite_form_basis()

    def random_element(self, distribution="default", *args, **kwds):
        """
        Return a random element from the lattice.

        Let $B$ be the current basis of this lattice. Then the following
        distributions are supported.

        DISTRIBUTIONS:

        - ``default`` - return $v \cdot B$ where $v$ is sampled by calling
          ``FreeModule(ZZ, B.rank()).random_element(*args, **kwds)``.

        INPUT:

        - ``distribution`` - see above.
        - ``*args`` passed through
        - ``*kwds`` passed through


        EXAMPLE::

            sage: K.<a> = NumberField(x^8+1)
            sage: O = K.ring_of_integers()
            sage: f = O.random_element(); f
            a^7 + 2*a^6 - a^5 + a^4 + 2*a - 8
            sage: L = RealLattice(f)
            sage: v = L.random_element(); v
            (-7, 5, -8, -9, 19, -1, 12, -38)
            sage: g = sum(v[i]*a^i for i in range(8))
            sage: g in ideal(f)
            True

        """
        if distribution == "default":
            return (ZZ**self.rank()).random_element() * self._basis
        else:
            raise NotImplementedError("Distribution '%s' not implemented.")

    @cached_method
    def voronoi_cell(self, radius=None):
        """
        Compute the Voronoi cell of a lattice, returning a Polyhedron.

        INPUT:

        - ``radius`` -- radius of ball containing considered vertices
          (default: automatic determination).

        OUTPUT:

        The Voronoi cell as a Polyhedron instance.

        The result is cached so that subsequent calls to this function
        return instantly.

        EXAMPLES::

            sage: L = RealLattice([[1, 0], [0, 1]])
            sage: V = L.voronoi_cell()
            sage: V.Vrepresentation()
            (A vertex at (1/2, -1/2), A vertex at (1/2, 1/2), A vertex at (-1/2, 1/2), A vertex at (-1/2, -1/2))

        The volume of the Voronoi cell is the square root of the discriminant of the lattice::

            sage: L = RealLattice(Matrix(ZZ, 4, 4, [[0,0,1,-1],[1,-1,2,1],[-6,0,3,3,],[-6,-24,-6,-5]])); L
            Lattice of degree 4 and rank 4 over Integer Ring
            Basis matrix:
            [  0   0   1  -1]
            [  1  -1   2   1]
            [ -6   0   3   3]
            [ -6 -24  -6  -5]
            sage: V = L.voronoi_cell()
            sage: V.volume()
            678
            sage: sqrt(L.discriminant())
            678

        Lattices not having full dimension are handled as well::

            sage: L = RealLattice([[2, 0, 0], [0, 2, 0]])
            sage: V = L.voronoi_cell()
            sage: V.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, -1, 0) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0)

        ALGORITHM:

        Uses parts of the algorithm from [Vit1996].

        REFERENCES:

        .. [Vit1996] E. Viterbo, E. Biglieri. Computing the Voronoi Cell
          of a Lattice: The Diamond-Cutting Algorithm.
          IEEE Transactions on Information Theory, 1996.
        """
        if not self._is_LLL_reduced:
            self.LLL()

        from diamond_cutting import calculate_voronoi_cell
        return calculate_voronoi_cell(self._basis, radius=radius)

    def voronoi_relevant_vectors(self):
        """
        Compute the embedded vectors inducing the Voronoi cell.

        OUTPUT:

        The list of Voronoi relevant vectors.

        EXAMPLES::

            sage: L = RealLattice([[3, 0], [4, 0]])
            sage: L.voronoi_relevant_vectors()
            [(-1, 0), (1, 0)]
        """
        V = self.voronoi_cell()

        def defining_point(ieq):
            """
            Compute the point defining an inequality.

            INPUT:

            - ``ieq`` - an inequality in the form [c, a1, a2, ...]
              meaning a1 * x1 + a2 * x2 + ... ≦ c

            OUTPUT:

            The point orthogonal to the hyperplane defined by ``ieq``
            in twice the distance from the origin.
            """
            c = ieq[0]
            a = ieq[1:]
            n = sum(y ** 2 for y in a)
            return vector([2 * y * c / n for y in a])

        return [defining_point(ieq) for ieq in V.inequality_generator()]

    def closest_vector(self, t):
        """
        Compute the closest vector in the embedded lattice to a given vector.

        INPUT:

        - ``t`` -- the target vector to compute the closest vector to.

        OUTPUT:

        The vector in the lattice closest to ``t``.

        EXAMPLES::

            sage: L = RealLattice([[1, 0], [0, 1]])
            sage: L.closest_vector((-6, 5/3))
            (-6, 2)

        ALGORITHM:

        Uses the algorithm from [Mic2010].

        REFERENCES:

        .. [Mic2010] D. Micciancio, P. Voulgaris. A Deterministic Single
          Exponential Time Algorithm for Most Lattice Problems based on
          Voronoi Cell Computations.
          Proceedings of the 42nd ACM Symposium Theory of Computation, 2010.
        """
        voronoi_cell = self.voronoi_cell()

        def projection(M, v):
            Mt = M.transpose()
            P = Mt * (M * Mt) ** (-1) * M
            return P * v

        t = projection(matrix(self.basis), vector(t))

        def CVPP_2V(t, V, voronoi_cell):
            t_new = t
            while not voronoi_cell.contains(t_new.list()):
                v = max(V, key=lambda v: t_new * v / v.norm() ** 2)
                t_new = t_new - v
            return t - t_new

        V = self.voronoi_relevant_vectors()
        t = vector(t)
        p = 0
        while not (ZZ(2 ** p) * voronoi_cell).contains(t):
            p += 1
        t_new = t
        i = p
        while i >= 1:
            V_scaled = [v * (2 ** (i - 1)) for v in V]
            t_new = t_new - CVPP_2V(t_new, V_scaled, ZZ(2 ** (i - 1)) * voronoi_cell)
            i -= 1
        return t - t_new
