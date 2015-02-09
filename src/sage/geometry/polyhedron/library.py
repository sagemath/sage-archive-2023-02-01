r"""
Library of commonly used, famous, or interesting polytopes

REFERENCES:

..  [Fetter2012]
    Hans L. Fetter,
    "A Polyhedron Full of Surprises",
    Mathematics Magazine 85 (2012), no. 5, 334-342.
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.rings.all import Integer, RR, QQ, ZZ, RDF
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.combinat.permutation import Permutations
from sage.groups.perm_gps.permgroup_named import AlternatingGroup
from sage.misc.functional import norm
from sage.functions.other import sqrt, floor, ceil
from sage.functions.trig import sin, cos
from sage.misc.decorators import rename_keyword

from constructor import Polyhedron



#########################################################################
class Polytopes():
    """
    A class of constructors for commonly used, famous, or interesting
    polytopes.

    TESTS::

        sage: TestSuite(polytopes).run(skip='_test_pickling')
    """

    @staticmethod
    def orthonormal_1(dim_n=5):
        """
        A matrix of rational approximations to orthonormal vectors to
        ``(1,...,1)``.

        INPUT:

        - ``dim_n`` - the dimension of the vectors

        OUTPUT:

        A matrix over ``QQ`` whose rows are close to an orthonormal
        basis to the subspace normal to ``(1,...,1)``.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.library import Polytopes
            sage: m = Polytopes.orthonormal_1(5)
            sage: m
            [ 70711/100000   -7071/10000             0             0             0]
            [    1633/4000     1633/4000 -81649/100000             0             0]
            [   7217/25000    7217/25000    7217/25000  -43301/50000             0]
            [ 22361/100000  22361/100000  22361/100000  22361/100000  -44721/50000]
        """
        pb = []
        for i in range(0,dim_n-1):
            pb.append([1.0/(i+1)]*(i+1) + [-1] + [0]*(dim_n-i-2))
        m = matrix(RDF,pb)
        new_m = []
        for i in range(0,dim_n-1):
            new_m.append([RDF(100000*q/norm(m[i])).ceil()/100000 for q in m[i]])
        return matrix(QQ,new_m)

    @staticmethod
    def project_1(fpoint):
        """
        Take a ndim-dimensional point and projects it onto the plane
        perpendicular to (1,1,...,1).

        INPUT:

          - ``fpoint`` - a list of ndim numbers

        EXAMPLES::

            sage: from sage.geometry.polyhedron.library import Polytopes
            sage: Polytopes.project_1([1,1,1,1,2])
            [1/100000, 1/100000, 1/50000, -559/625]
        """
        dim_n = len(fpoint)
        p_basis = [list(q) for q in Polytopes.orthonormal_1(dim_n)]
        out_v = []
        for v in p_basis:
            out_v.append(sum([fpoint[ind]*v[ind] for ind in range(dim_n)]))
        return out_v

    @staticmethod
    def _pfunc(i,j,perm):
        """
        An internal utility function for constructing the Birkhoff polytopes.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.library import Polytopes
            sage: Polytopes._pfunc(1,2,Permutations(3)[0])
            0
        """
        if perm[i-1] == j:
            return 1
        else:
            return 0


    @rename_keyword(deprecation=11634, field='base_ring')
    def regular_polygon(self, n, base_ring=QQ):
        """
        Return a regular polygon with `n` vertices.  Over the rational
        field the vertices may not be exact.

        INPUT:

        - ``n`` -- a positive integer, the number of vertices.

        - ``base_ring`` -- a ring in which the coordinates will lie.

        EXAMPLES::

            sage: octagon = polytopes.regular_polygon(8)
            sage: len(octagon.vertices())
            8
            sage: polytopes.regular_polygon(3).vertices()
            (A vertex at (-125283617/144665060, -500399958596723/1000799917193445),
             A vertex at (0, 1),
             A vertex at (94875313/109552575, -1000799917193444/2001599834386889))
            sage: polytopes.regular_polygon(3, base_ring=RealField(100)).vertices()
            (A vertex at (0.00000000000000000000000000000, 1.0000000000000000000000000000),
             A vertex at (0.86602540378443864676372317075, -0.50000000000000000000000000000),
             A vertex at (-0.86602540378443864676372317076, -0.50000000000000000000000000000))
            sage: polytopes.regular_polygon(3, base_ring=RealField(10)).vertices()
            (A vertex at (0.00, 1.0),
             A vertex at (0.87, -0.50),
             A vertex at (-0.86, -0.50))
        """
        try:
            omega = 2*base_ring.pi()/n
        except AttributeError:
            omega = 2*RR.pi()/n
        verts = []
        for i in range(n):
            t = omega*i
            verts.append([base_ring(t.sin()), base_ring(t.cos())])
        return Polyhedron(vertices=verts, base_ring=base_ring)


    def Birkhoff_polytope(self, n):
        """
        Return the Birkhoff polytope with n! vertices.  Each vertex
        is a (flattened) n by n permutation matrix.

        INPUT:

        - ``n`` -- a positive integer giving the size of the permutation matrices.

        EXAMPLES::

            sage: b3 = polytopes.Birkhoff_polytope(3)
            sage: b3.n_vertices()
            6
        """
        verts = []
        for p in Permutations(range(1,n+1)):
            verts += [ [Polytopes._pfunc(i,j,p) for j in range(1,n+1)
                        for i in range(1,n+1) ] ]
        return Polyhedron(vertices=verts)


    def n_simplex(self, dim_n=3, project = True):
        """
        Return a rational approximation to a regular simplex in
        dimension ``dim_n``.

        INPUT:

        - ``dim_n`` -- The dimension of the simplex, a positive
          integer.

        - ``project`` -- Optional argument, whether to project
          orthogonally.  Default is True.

        OUTPUT:

        A Polyhedron object of the ``dim_n``-dimensional simplex.

        EXAMPLES::

            sage: s5 = polytopes.n_simplex(5)
            sage: s5.dim()
            5
        """
        verts = Permutations([0 for i in range(dim_n)] + [1]).list()
        if project: verts = [Polytopes.project_1(x) for x in verts]
        return Polyhedron(vertices=verts)


    @rename_keyword(deprecation=11634, field='base_ring')
    def icosahedron(self, base_ring=QQ):
        """
        Return an icosahedron with edge length 1.

        INPUT:

        - ``base_ring`` -- Either ``QQ`` or ``RDF``.

        OUTPUT:

        A Polyhedron object of a floating point or rational
        approximation to the regular 3d icosahedron.

        If ``base_ring=QQ``, a rational approximation is used and the
        points are not exactly the vertices of the icosahedron. The
        icosahedron's coordinates contain the golden ratio, so there
        is no exact representation possible.

        EXAMPLES::

            sage: ico = polytopes.icosahedron()
            sage: sum(sum( ico.vertex_adjacency_matrix() ))/2
            30
        """
        if base_ring == QQ:
            g = QQ(1618033)/1000000 # Golden ratio approximation
            r12 = QQ(1)/2
        elif base_ring == RDF:
            g = RDF( (1 + sqrt(5))/2 )
            r12 = RDF( QQ(1)/2 )
        else:
            raise ValueError("field must be QQ or RDF.")
        verts = [i([0,r12,g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,r12,-g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,-r12,g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,-r12,-g/2]) for i in AlternatingGroup(3)]
        return Polyhedron(vertices=verts, base_ring=base_ring)


    @rename_keyword(deprecation=11634, field='base_ring')
    def dodecahedron(self, base_ring=QQ):
        """
        Return a dodecahedron.

        INPUT:

        - ``base_ring`` -- Either ``QQ`` (in which case a rational
          approximation to the golden ratio is used) or ``RDF``.

        EXAMPLES::

            sage: d12 = polytopes.dodecahedron()
            sage: d12.n_inequalities()
            12
        """
        return self.icosahedron(base_ring=base_ring).polar()


    def small_rhombicuboctahedron(self):
        """
        Return an Archimedean solid with 24 vertices and 26 faces.

        EXAMPLES::

            sage: sr = polytopes.small_rhombicuboctahedron()
            sage: sr.n_vertices()
            24
            sage: sr.n_inequalities()
            26
        """
        verts = [ [-3/2, -1/2, -1/2], [-3/2, -1/2, 1/2], [-3/2, 1/2, -1/2],
                  [-3/2, 1/2, 1/2], [-1/2, -3/2, -1/2], [-1/2, -3/2, 1/2],
                  [-1/2, -1/2, -3/2], [-1/2,-1/2, 3/2], [-1/2, 1/2, -3/2],
                  [-1/2, 1/2, 3/2], [-1/2, 3/2, -1/2], [-1/2, 3/2, 1/2],
                  [1/2, -3/2, -1/2], [1/2, -3/2, 1/2], [1/2, -1/2,-3/2],
                  [1/2, -1/2, 3/2], [1/2, 1/2, -3/2], [1/2, 1/2, 3/2],
                  [1/2, 3/2,-1/2], [1/2, 3/2, 1/2], [3/2, -1/2, -1/2],
                  [3/2, -1/2, 1/2], [3/2, 1/2,-1/2], [3/2, 1/2, 1/2] ]
        return Polyhedron(vertices=verts)


    @rename_keyword(deprecation=11634, field='base_ring')
    def great_rhombicuboctahedron(self, base_ring=QQ):
        """
        Return an Archimedean solid with 48 vertices and 26 faces.

        EXAMPLES::

            sage: gr = polytopes.great_rhombicuboctahedron()
            sage: gr.n_vertices()
            48
            sage: gr.n_inequalities()
            26
        """
        v1 = QQ(131739771357/54568400000)
        v2 = QQ(104455571357/27284200000)
        verts = [ [1, v1, v2],
                  [1, v2, v1],
                  [v1, 1, v2],
                  [v1, v2, 1],
                  [v2, 1, v1],
                  [v2, v1, 1] ]
        verts = verts + [[x[0],x[1],-x[2]] for x in verts]
        verts = verts + [[x[0],-x[1],x[2]] for x in verts]
        verts = verts + [[-x[0],x[1],x[2]] for x in verts]
        if base_ring!=QQ:
            verts = [base_ring(v) for v in verts]
        return Polyhedron(vertices=verts, base_ring=base_ring)


    def rhombic_dodecahedron(self):
        """
        This face-regular, vertex-uniform polytope is dual to the
        cuboctahedron. It has 14 vertices and 12 faces.

        EXAMPLES::

            sage: rd = polytopes.rhombic_dodecahedron()
            sage: rd.n_vertices()
            14
            sage: rd.n_inequalities()
            12
        """
        v = [ [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1], [-1, 1, 1],
              [-1, 1, -1], [-1, -1, 1], [-1, -1, -1], [0, 0, 2], [0, 2, 0],
              [2, 0, 0], [0, 0, -2], [0, -2, 0], [-2, 0, 0] ]
        return Polyhedron(vertices=v)


    def cuboctahedron(self):
        """
        An Archimedean solid with 12 vertices and 14 faces.  Dual to
        the rhombic dodecahedron.

        EXAMPLES::

            sage: co = polytopes.cuboctahedron()
            sage: co.n_vertices()
            12
            sage: co.n_inequalities()
            14
        """
        one = Integer(1)
        v = [ [0, -one/2, -one/2], [0, one/2, -one/2], [one/2, -one/2, 0],
              [one/2, one/2, 0], [one/2, 0, one/2], [one/2, 0, -one/2],
              [0, one/2, one/2], [0, -one/2, one/2], [-one/2, 0, one/2],
              [-one/2, one/2, 0], [-one/2, 0, -one/2], [-one/2, -one/2, 0] ]
        return Polyhedron(vertices=v)


    @rename_keyword(deprecation=11634, field='base_ring')
    def buckyball(self, base_ring=QQ):
        """
        Also known as the truncated icosahedron, an Archimedean solid.
        It has 32 faces and 60 vertices.  Rational coordinates are not
        exact.

        EXAMPLES::

            sage: bb = polytopes.buckyball()
            sage: bb.n_vertices()
            60
            sage: bb.n_inequalities()   # number of facets
            32
            sage: bb.base_ring()
            Rational Field
        """
        # Note: QQ would give some incorrecty subdivided facets
        p = self.icosahedron(base_ring=RDF).edge_truncation()
        if base_ring==RDF:
            return p
        # Converting with low precision to save time.
        new_ieqs = [[int(1000*x)/QQ(1000) for x in y] for y in p.inequalities()]
        return Polyhedron(ieqs=new_ieqs)


    def pentakis_dodecahedron(self):
        """
        This face-regular, vertex-uniform polytope is dual to the
        truncated icosahedron.  It has 60 faces and 32 vertices.

        EXAMPLES::

            sage: pd = polytopes.pentakis_dodecahedron()
            sage: pd.n_vertices()
            32
            sage: pd.n_inequalities()   # number of facets
            60
        """
        return self.buckyball().polar()

    def Kirkman_icosahedron(self):
        """
        A non-uniform icosahedron with interesting properties.

        See [Fetter2012]_ for details.

        OUTPUT:

        The Kirkman icosahedron, a 3-dimensional polyhedron
        with 20 vertices, 20 faces, and 38 edges.

        EXAMPLES::

            sage: KI = polytopes.Kirkman_icosahedron()
            sage: KI.f_vector()
            (1, 20, 38, 20, 1)
            sage: vertices = KI.vertices()
            sage: edges = [[vector(edge[0]),vector(edge[1])] for edge in KI.bounded_edges()]
            sage: edge_lengths = [norm(edge[0]-edge[1]) for edge in edges]
            sage: union(edge_lengths)
            [7, 8, 9, 11, 12, 14, 16]
        """
        vertices = [[-12, -4, 0], [-12, 4, 0], [-9, -6, -6],
                    [-9, -6, 6], [-9, 6, -6], [-9, 6, 6], [-6, 0, -12],
                    [-6, 0, 12], [0, -12, -8], [0, -12, 8], [0, 12, -8],
                    [0, 12, 8], [6, 0, -12], [6, 0, 12], [9, -6, -6],
                    [9, -6, 6], [9, 6, -6], [9, 6, 6], [12, -4, 0],
                    [12, 4, 0]]
        return Polyhedron(vertices=vertices)


    def twenty_four_cell(self):
        """
        Return the standard 24-cell polytope.

        OUTPUT:

        A Polyhedron object of the 4-dimensional 24-cell, a regular
        polytope. The coordinates of this polytope are exact.

        EXAMPLES::

            sage: p24 = polytopes.twenty_four_cell()
            sage: v = p24.vertex_generator().next()
            sage: for adj in v.neighbors(): print adj
            A vertex at (-1/2, -1/2, -1/2, 1/2)
            A vertex at (-1/2, -1/2, 1/2, -1/2)
            A vertex at (-1, 0, 0, 0)
            A vertex at (-1/2, 1/2, -1/2, -1/2)
            A vertex at (0, -1, 0, 0)
            A vertex at (0, 0, -1, 0)
            A vertex at (0, 0, 0, -1)
            A vertex at (1/2, -1/2, -1/2, -1/2)
        """
        verts = []
        q12 = QQ(1)/2
        base = [q12,q12,q12,q12]
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        verts.append([x for x in base])
                        base[3] = base[3]*(-1)
                    base[2] = base[2]*(-1)
                base[1] = base[1]*(-1)
            base[0] = base[0]*(-1)
        verts = verts + Permutations([0,0,0,1]).list()
        verts = verts + Permutations([0,0,0,-1]).list()
        return Polyhedron(vertices=verts)


    def six_hundred_cell(self):
        """
        Return the standard 600-cell polytope.

        OUTPUT:

        A Polyhedron object of the 4-dimensional 600-cell, a regular
        polytope.  In many ways this is an analogue of the
        icosahedron.  The coordinates of this polytope are rational
        approximations of the true coordinates of the 600-cell, some
        of which involve the (irrational) golden ratio.

        EXAMPLES::

            sage: p600 = polytopes.six_hundred_cell() # not tested - very long time
            sage: len(list(p600.bounded_edges())) # not tested - very long time
            120
        """
        verts = []
        q12 = QQ(1)/2
        base = [q12,q12,q12,q12]
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        verts.append([x for x in base])
                        base[3] = base[3]*(-1)
                    base[2] = base[2]*(-1)
                base[1] = base[1]*(-1)
            base[0] = base[0]*(-1)
        verts += Permutations([0,0,0,1]).list()
        verts += Permutations([0,0,0,-1]).list()
        g = QQ(1618033)/1000000 # Golden ratio approximation
        verts = verts + [i([q12,g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,-g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,-g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,-g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,-g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        return Polyhedron(vertices=verts)


    @rename_keyword(deprecation=11634, field='base_ring')
    def cyclic_polytope(self, dim_n, points_n, base_ring=QQ):
        """
        Return a cyclic polytope.

        INPUT:

        - ``dim_n`` -- positive integer. the dimension of the polytope.

        - ``points_n`` -- positive integer. the number of vertices.

        - ``base_ring`` -- either ``QQ`` (default) or ``RDF``.

        OUTPUT:

        A cyclic polytope of dim_n with points_n vertices on the
        moment curve ``(t,t^2,...,t^n)``, as Polyhedron object.

        EXAMPLES::

            sage: c = polytopes.cyclic_polytope(4,10)
            sage: c.n_inequalities()
            35
        """
        verts = [[t**i for i in range(1,dim_n+1)] for t in range(points_n)]
        return Polyhedron(vertices=verts, base_ring=base_ring)


    def hypersimplex(self, dim_n, k, project = True):
        """
        The hypersimplex in dimension dim_n with d choose k vertices,
        projected into (dim_n - 1) dimensions.

        INPUT:

        - ``n`` -- the numbers ``(1,...,n)`` are permuted

        - ``project`` -- If ``False``, the polyhedron is left in
          dimension ``n``.

        OUTPUT:

        A Polyhedron object representing the hypersimplex.

        EXAMPLES::

            sage: h_4_2 = polytopes.hypersimplex(4,2) # combinatorially equivalent to octahedron
            sage: h_4_2.n_vertices()
            6
            sage: h_4_2.n_inequalities()
            8
        """
        vert0 = [0]*(dim_n-k) + [1]*k
        verts = Permutations(vert0).list()
        if project:
            verts = [Polytopes.project_1(x) for x in verts]
        return Polyhedron(vertices=verts)


    def permutahedron(self, n, project = True):
        """
        The standard permutahedron of (1,...,n) projected into n-1
        dimensions.

        INPUT:

        - ``n`` -- the numbers ``(1,...,n)`` are permuted

        - ``project`` -- If ``False`` the polyhedron is left in dimension ``n``.

        OUTPUT:

        A Polyhedron object representing the permutahedron.

        EXAMPLES::

            sage: perm4 = polytopes.permutahedron(4)
            sage: perm4
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices
            sage: polytopes.permutahedron(5).show()    # long time
            Graphics3d Object
        """
        verts = range(1,n+1)
        verts = Permutations(verts).list()
        if project:
            verts = [Polytopes.project_1(x) for x in verts]
        p = Polyhedron(vertices=verts)
        return p


    def n_cube(self, dim_n):
        """
        Return a cube in the given dimension

        INPUT:

        - ``dim_n`` -- integer. The dimension of the cube.

        OUTPUT:

        A Polyhedron object of the ``dim_n``-dimensional cube, with
        exact coordinates.

        EXAMPLES::

            sage: four_cube = polytopes.n_cube(4)
            sage: four_cube.is_simple()
            True
        """
        if dim_n == 1:
            return Polyhedron(vertices = [[1],[-1]])

        pre_cube = polytopes.n_cube(dim_n-1)
        vertices = [];
        for pre_v in pre_cube.vertex_generator():
            vertices.append( [ 1] + [v for v in pre_v] );
            vertices.append( [-1] + [v for v in pre_v] );
        return Polyhedron(vertices = vertices)


    def cross_polytope(self, dim_n):
        """
        Return a cross-polytope in dimension ``dim_n``. These are
        the generalization of the octahedron.

        INPUT:

        - ``dim_n`` -- integer. The dimension of the cross-polytope.

        OUTPUT:

        A Polyhedron object of the ``dim_n``-dimensional cross-polytope,
        with exact coordinates.

        EXAMPLES::

            sage: four_cross = polytopes.cross_polytope(4)
            sage: four_cross.is_simple()
            False
            sage: four_cross.n_vertices()
            8
        """
        verts = Permutations([0 for i in range(dim_n-1)] + [1]).list()
        verts += Permutations([0 for i in range(dim_n-1)] + [-1]).list()
        return Polyhedron(vertices=verts)


    def parallelotope(self, generators):
        r"""
        Return the parallelotope spanned by the generators.

        INPUT:

        - ``generators`` -- an iterable of anything convertible to vector
          (for example, a list of vectors) such that the vectors all
          have the same dimension.

        OUTPUT:

        The parallelotope. This is the multi-dimensional
        generalization of a parallelogram (2 generators) and a
        parallelepiped (3 generators).

        EXAMPLES::

            sage: polytopes.parallelotope([ (1,0), (0,1) ])
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: polytopes.parallelotope([[1,2,3,4],[0,1,0,7],[3,1,0,2],[0,0,1,0]])
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 16 vertices
        """
        try:
            generators = [ vector(QQ,v) for v in generators ]
            base_ring = QQ
        except TypeError:
            generators = [ vector(RDF,v) for v in generators ]
            base_ring = RDF

        from sage.combinat.combination import Combinations
        par =  [ 0*generators[0] ]
        par += [ sum(c) for c in Combinations(generators) if c!=[] ]
        return Polyhedron(vertices=par, base_ring=base_ring)



polytopes = Polytopes()
