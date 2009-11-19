from sage.misc.misc import prod
from sage.combinat.family import Family
from root_lattice_realization import RootLatticeRealization

class WeightLatticeRealization(RootLatticeRealization):
    def _test_weight_lattice_realization(self, **options):
        """
        Runs sanity checks on this weight lattice realization:
         - scalar products between the fundamental weights and simple coroots
         - rho, highest_root, ...

        See also: :class:`Test``

        EXAMPLES::
            sage: RootSystem(['A',3]).root_lattice()._test_root_lattice_realization()
        """
        tester = self._tester(**options)
        Lambda     = self.fundamental_weights()
        alphacheck = self.simple_coroots()
        tester.assertEqual(Lambda.keys(), self.index_set())

        for i in self.index_set():
            assert(Lambda[i].is_dominant())
            for j in self.index_set():
                tester.assertEqual(Lambda[j].scalar(alphacheck[i]), (1 if i==j else 0))

        tester.assert_(self.rho().is_dominant())
        if self.root_system.is_finite() and self.root_system.is_irreducible():
            tester.assert_(self.highest_root().is_dominant())

    # Should this be a method or an attribute?
    # same question for the roots, ...
    # Should this use rename to set a nice name for this family?
    def fundamental_weights(self):
        r"""
        Returns the family `(\Lambda_i)_{i\in I}` of the fundamental weights.

        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: f = e.fundamental_weights()
            sage: [f[i] for i in [1,2,3]]
            [(1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0)]
        """
        if not hasattr(self,"_fundamental_weights"):
            self._fundamental_weights = Family(self.index_set(),
                                               self.fundamental_weight)
            # self._fundamental_weights.rename("Lambda")
            # break some doctests.
        return self._fundamental_weights

    def rho(self):
        """
        EXAMPLES::

            sage: RootSystem(['A',3]).ambient_lattice().rho()
            (3, 2, 1, 0)
        """
        return sum(self.fundamental_weights())

    def signs_of_alcovewalk(self, walk):
        r"""
        Let walk = `[i_1,\dots,i_n]` denote an alcove walk starting
        from the fundamental alcove `y_0`, crossing at step 1 the
        wall `i_1`, and so on.

        For each `k`, set `w_k = s_{i_1} \circ s_{i_k}`, and denote
        by `y_k = w_k(y_0)` the alcove reached after `k` steps. Then,
        `y_k` is obtained recursively from `y_{k-1}` by applying the
        following reflection:

        .. math::

              y_k = s_{w_{k-1} \alpha_{i_k}} y_{k-1}

        The step is said positive if `w_{k-1} \alpha_{i_k}` is a
        negative root (considering `w_{k-1}` as element of the classical
        Weyl group and `\alpha_{i_k}` as a classical root) and
        negative otherwise.

        This function returns a list of the form `[+1,+1,-1,...]`,
        where the `k^{th}` entry denotes whether the `k^{th}` step was
        positive or negative

        See equation 3.4, of Ram: Alcove walks ..., arxiv:math/0601343v1 [math.RT]

        EXAMPLES::

            sage: L = RootSystem(['C',2,1]).weight_lattice()
            sage: L.signs_of_alcovewalk([1,2,0,1,2,1,2,0,1,2])
            [-1, -1, 1, -1, 1, 1, 1, 1, 1, 1]
            sage: L = RootSystem(['A',2,1]).weight_lattice()
            sage: L.signs_of_alcovewalk([0,1,2,1,2,0,1,2,0,1,2,0])
            [1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1]
        """
        lattice_classical = self.root_system.cartan_type().classical().root_system().ambient_space()
        W = lattice_classical.weyl_group()
        simple_reflections = W.simple_reflections()
        alphacheck = lattice_classical.alphacheck()
        rho = lattice_classical.rho()
        word = W.unit()
        signs = []
        for s in walk:
            if ((alphacheck[s]).scalar((word).action(rho)) > 0):
                signs.append(-1)
            else:
                signs.append(1)
            word = simple_reflections[s]*word
        return signs

    def rho_classical(self):
        """
        For an affine type in a weight space, rho_classical is the analog of
        rho in the classical parabolic subgroup. it lives in the level 0.

        EXAMPLES::

            sage: RootSystem(['C',4,1]).weight_space().rho_classical()
            -4*Lambda[0] + Lambda[1] + Lambda[2] + Lambda[3] + Lambda[4]
            sage: WS = RootSystem(['D',4,1]).weight_space()
            sage: WS.rho_classical().scalar(WS.null_coroot())
            0
        """
        rho = self.rho()
        Lambda = self.fundamental_weights()
        return rho - (rho.level()/Lambda[0].level()) * Lambda[0]



    # Should it be a method of highest_weight?
    def weyl_dimension(self, highest_weight):
        """
        EXAMPLES::

            sage: RootSystem(['A',3]).ambient_lattice().weyl_dimension([2,1,0,0])
            20
        """
        highest_weight = self(highest_weight)
        assert(highest_weight.is_dominant())
        rho = self.rho()
        n = prod([(rho+highest_weight).dot_product(x) for x in self.positive_roots()])
        d = prod([ rho.dot_product(x) for x in self.positive_roots()])
        return n/d

    def plot(self, size=[[0],[0]], projection='usual', simple_roots=True, fundamental_weights=True, alcovewalks=[]):
        r"""
        Return a graphics object built from a space of weight(space/lattice).
        There is a different technic to plot if the Cartan type is affine or not.
        The graphics returned is a Graphics object.

        This function is experimental, and is subject to short term evolutions.

        EXAMPLES::

          By default, the plot returned has no axes and the ratio between axes is 1.
            sage: G = RootSystem(['C',2]).weight_lattice().plot()
            sage: G.axes(True)
            sage: G.set_aspect_ratio(2)

          For a non affine Cartan type, the plot method work for type with 2 generators,
          it will draw the hyperlane(line for this dimension) accrow the fundamentals weights.
            sage: G = RootSystem(['A',2]).weight_lattice().plot()
            sage: G = RootSystem(['B',2]).weight_lattice().plot()
            sage: G = RootSystem(['G',2]).weight_lattice().plot()

          The plot returned has a size of one fundamental polygon by default. We can
          ask plot to give a bigger plot by using the argument size
            sage: G = RootSystem(['G',2,1]).weight_space().plot(size = [[0..1],[-1..1]])
            sage: G = RootSystem(['A',2,1]).weight_space().plot(size = [[-1..1],[-1..1]])

          A very important argument is the projection which will draw the plot. There are
          some usual projections is this method. If you want to draw in the plane a very
          special Cartan type, Sage will ask you to specify the projection. The projection
          is a matrix over a ring. In practice, calcul over float is a good way to draw.
            sage: L = RootSystem(['A',2,1]).weight_space()
            sage: G = L.plot(projection=matrix(RR, [[0,0.5,-0.5],[0,0.866,0.866]]))
            sage: G = RootSystem(['C',2,1]).weight_space().plot()

          By default, the plot method draw the simple roots, this can be disabled by setting
          the argument simple_roots=False
            sage: G = RootSystem(['A',2]).weight_space().plot(simple_roots=False)

          By default, the plot method draw the fundamental weights,this can be disabled by
          setting the argument fundamental_weights=False
            sage: G = RootSystem(['A',2]).weight_space().plot(fundamental_weights=False, simple_roots=False)

          There is in a plot an argument to draw alcoves walks. The good way to do this is
          to use the crystals theory. the plot method contains only the drawing part...
            sage: L = RootSystem(['A',2,1]).weight_space()
            sage: G = L.plot(size=[[-1..1],[-1..1]],alcovewalks=[[0,2,0,1,2,1,2,0,2,1]])
        """

        from sage.plot.plot import Graphics
        from sage.plot.line import line
        from cartan_type import CartanType
        from sage.matrix.constructor import matrix
        from sage.rings.all import QQ, RR
        from sage.plot.arrow import arrow
        from sage.plot.point import point

        # We begin with an empty plot G
        G = Graphics()

        ct = self.cartan_type()
        n = ct.n

        # Define a set of colors
        # TODO : Colors in option ?
        colors=[(0,1,0),(1,0,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1)]

        # plot the affine types:
        if ct.is_affine():

            # Check the projection
            # TODO : try to have usual_projection for main plotable types
            if projection == 'usual':
                if ct == CartanType(['A',2,1]):
                    projection = matrix(RR, [[0,0.5,-0.5],[0,0.866,0.866]])
                elif ct == CartanType(['C',2,1]):
                    projection = matrix(QQ, [[0,1,1],[0,0,1]])
                elif ct == CartanType(['G',2,1]):
                    projection = matrix(RR, [[0,0.5,0],[0,0.866,1.732]])
                else:
                    raise 'There is no usual projection for this Cartan type, you have to give one in argument'

            assert(n + 1 == projection.ncols())
            assert(2 == projection.nrows())

            # Check the size is correct with the lattice
            assert(len(size) == n)

            # Select the center of the translated fundamental polygon to plot
            translation_factors = ct.translation_factors()
            simple_roots = self.simple_roots()
            translation_vectors = [translation_factors[i]*simple_roots[i] for i in ct.classical().index_set()]

            initial = [[]]
            for i in range(n):
                prod_list = []
                for elem in size[i]:
                    for partial_list in initial:
                        prod_list.append( [elem]+partial_list );
                initial = prod_list;

            part_lattice = []
            for combinaison in prod_list:
                elem_lattice = self.zero()
                for i in range(n):
                    elem_lattice = elem_lattice + combinaison[i]*translation_vectors[i]
                part_lattice.append(elem_lattice)

            # Get the vertices of the fundamental alcove
            fundamental_weights = self.fundamental_weights()
            vertices = map(lambda x: (1/x.level())*x, fundamental_weights.list())

            # Recup the group which act on the fundamental polygon
            classical = self.weyl_group().classical()

            for center in part_lattice:
                for w in classical:
                    # for each center of polygon and each element of classical
                    # parabolic subgroup, we have to draw an alcove.

                    #first, iterate over pairs of fundamental weights, drawing lines border of polygons:
                    for i in range(1,n+1):
                        for j in range(i+1,n+1):
                            p1=projection*((w.action(vertices[i])).to_vector() + center.to_vector())
                            p2=projection*((w.action(vertices[j])).to_vector() + center.to_vector())
                            G+=line([p1,p2],rgbcolor=(0,0,0),thickness=2)

                    #next, get all lines from point to a fundamental weight, that separe different
                    #chanber in a same polygon (important: associate a color with a fundamental weight)
                    pcenter = projection*(center.to_vector())
                    for i in range(1,n+1):
                        p3=projection*((w.action(vertices[i])).to_vector() + center.to_vector())
                        G+=line([p3,pcenter], rgbcolor=colors[n-i+1])

            #Draw alcovewalks
            #FIXME : The good way to draw this is to use the alcoves walks works made in Cristals
            #The code here just draw like example and import the good things.
            rho = (1/self.rho().level())*self.rho()
            W = self.weyl_group()
            for walk in alcovewalks:
                target = W.from_reduced_word(walk).action(rho)
                for i in range(len(walk)):
                    walk.pop()
                    origin = W.from_reduced_word(walk).action(rho)
                    G+=arrow(projection*(origin.to_vector()),projection*(target.to_vector()), rgbcolor=(0.6,0,0.6), width=1, arrowsize=5)
                    target = origin

        else:
            # non affine plot

            # Check the projection
            # TODO : try to have usual_projection for main plotable types
            if projection == 'usual':
                if ct == CartanType(['A',2]):
                    projection = matrix(RR, [[0.5,-0.5],[0.866,0.866]])
                elif ct == CartanType(['B',2]):
                    projection = matrix(QQ, [[1,0],[1,1]])
                elif ct == CartanType(['C',2]):
                    projection = matrix(QQ, [[1,1],[0,1]])
                elif ct == CartanType(['G',2]):
                    projection = matrix(RR, [[0.5,0],[0.866,1.732]])
                else:
                    raise 'There is no usual projection for this Cartan type, you have to give one in argument'

            # Get the fundamental weights
            fundamental_weights = self.fundamental_weights()
            WeylGroup = self.weyl_group()

            #Draw not the alcove but the cones delimited by the hyperplanes
            #The size of the line depend of the fundamental weights.
            pcenter = projection*(self.zero().to_vector())
            for w in WeylGroup:
                for i in range(1,n+1):
                    p3=3*projection*((w.action(fundamental_weights[i])).to_vector())
                    G+=line([p3,pcenter], rgbcolor=colors[n-i+1])

        #Draw the simple roots
        if simple_roots:
            SimpleRoots = self.simple_roots()
            if ct.is_affine():
                G+=arrow((0,0), projection*(SimpleRoots[0].to_vector()), rgbcolor=(0,0,0))
            for j in range(1,n+1):
                G+=arrow((0,0),projection*(SimpleRoots[j].to_vector()), rgbcolor=colors[j])

        #Draw the fundamental weights
        if fundamental_weights:
            FundWeight = self.fundamental_weights()
            for j in range(1,n+1):
                G+=point(projection*(FundWeight[j].to_vector()), rgbcolor=colors[j], pointsize=60)

        G.set_aspect_ratio(1)
        G.axes(False)
        return G

