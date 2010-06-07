r"""
Crystals
"""
#*****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.cachefunc import CachedFunction
from sage.misc.abstract_method import abstract_method
from sage.categories.category import Category
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.latex import latex
from sage.combinat import ranker

class Crystals(Category):
    """
    The category of crystals

    See :mod:`sage.combinat.crystals` for an introduction to crystals.

    EXAMPLES::

        sage: C = Crystals()
        sage: C
        Category of crystals
        sage: C.super_categories()
        [Category of... enumerated sets]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    Parents in this category should implement the following methods:

     - either a method ``cartan_type`` or an attribute ``_cartan_type``

     - ``module_generators``: a list (or container) of distinct elements
       which generate the crystal using `f_i`

    Furthermore, their elements should implement the following methods:

     - x.e(i) (returning `e_i(x)`)

     - x.f(i) (returning `f_i(x)`)

    EXAMPLES::

        sage: from sage.misc.abstract_method import abstract_methods_of_class
        sage: abstract_methods_of_class(Crystals().element_class)
        {'required': ['e', 'f'], 'optional': []}

    TESTS::

        sage: TestSuite(C).run()
        sage: B = Crystals().example()
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: Crystals().super_categories()
            [Category of enumerated sets]
        """
        return [EnumeratedSets()]

    class ParentMethods:

        #TODO: implement tests for the TestSuite to check whether the crystal satisfies the Stembridge rules

        def _an_element_(self):
            """
            Returns an element of ``self``

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.an_element()
                1
            """
            return self.first()

        def weight_lattice_realization(self):
            """
            Returns the weight lattice realization for the root system
            associated to ``self``. This default implementation uses
            the ambient space of the root system.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.weight_lattice_realization()
                Ambient space of the Root system of type ['A', 5]
                sage: K = KirillovReshetikhinCrystal(['A',2,1], 1, 1)
                sage: K.weight_lattice_realization()
                Weight lattice of the Root system of type ['A', 2, 1]
            """
            F = self.cartan_type().root_system()
            if F.ambient_space() is None:
                return F.weight_lattice()
            else:
                return F.ambient_space()

        def cartan_type(self):
            """
            Returns the Cartan type of the crystal

            EXAMPLES::
                sage: C = CrystalOfLetters(['A',2])
                sage: C.cartan_type()
                ['A', 2]
            """
            return self._cartan_type

        def index_set(self):
            """
            Returns the index set of the Dynkin diagram underlying the crystal

            EXAMPLES:
                sage: C = CrystalOfLetters(['A', 5])
                sage: C.index_set()
                [1, 2, 3, 4, 5]
            """
            return self.cartan_type().index_set()

        def Lambda(self):
            """
            Returns the fundamentals weights in the weight lattice
            realization for the root system associated the crystal

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.Lambda()
                Finite family {1: (1, 0, 0, 0, 0, 0), 2: (1, 1, 0, 0, 0, 0), 3: (1, 1, 1, 0, 0, 0), 4: (1, 1, 1, 1, 0, 0), 5: (1, 1, 1, 1, 1, 0)}
            """
            return self.weight_lattice_realization().fundamental_weights()

        def crystal_morphism(self, g, index_set = None, automorphism = lambda i : i, direction = 'down', direction_image = 'down',
                             similarity_factor = None, similarity_factor_domain = None, cached = False, acyclic = True):
            """
            Constructs a morphism from the crystal ``self`` to another crystal.
            The input `g` can either be a function of a (sub)set of elements of self to
            element in another crystal or a dictionary between certain elements.
            Usually one would map highest weight elements or crystal generators to each
            other using g.
            Specifying index_set gives the opportunity to define the morphism as `I`-crystals
            where `I =` index_set. If index_set is not specified, the index set of self is used.
            It is also possible to define twisted-morphisms by specifying an automorphism on the
            nodes in te Dynkin diagram (or the index_set).
            The option direction and direction_image indicate whether to use `f_i` or `e_i` in
            self or the image crystal to construct the morphism, depending on whether the direction
            is set to 'down' or 'up'.
            It is also possible to set a similarity_factor. This should be a dictionary between
            the elements in the index set and positive integers. The crystal operator `f_i` then gets
            mapped to `f_i^{m_i}` where `m_i =` similarity_factor[i].
            Setting similarity_factor_domain to a dictionary between the index set and positive integers
            has the effect that `f_i^{m_i}` gets mapped to `f_i` where `m_i =` similarity_factor_domain[i].
            Finally, it is possible to set the option `acyclic = False`. This calculates an isomorphism
            for cyclic crystals (for example finite affine crystals). In this case the input function `g`
            is supposed to be given as a dictionary.

            EXAMPLES::

                sage: C2 = CrystalOfLetters(['A',2])
                sage: C3 = CrystalOfLetters(['A',3])
                sage: g = {C2.module_generators[0] : C3.module_generators[0]}
                sage: g_full = C2.crystal_morphism(g)
                sage: g_full(C2(1))
                1
                sage: g_full(C2(2))
                2
                sage: g = {C2(1) : C2(3)}
                sage: g_full = C2.crystal_morphism(g, automorphism = lambda i : 3-i, direction_image = 'up')
                sage: [g_full(b) for b in C2]
                [3, 2, 1]
                sage: T = CrystalOfTableaux(['A',2], shape = [2])
                sage: g = {C2(1) : T(rows=[[1,1]])}
                sage: g_full = C2.crystal_morphism(g, similarity_factor = {1:2, 2:2})
                sage: [g_full(b) for b in C2]
                [[[1, 1]], [[2, 2]], [[3, 3]]]
                sage: g = {T(rows=[[1,1]]) : C2(1)}
                sage: g_full = T.crystal_morphism(g, similarity_factor_domain = {1:2, 2:2})
                sage: g_full(T(rows=[[2,2]]))
                2

                sage: B1=KirillovReshetikhinCrystal(['A',2,1],1,1)
                sage: B2=KirillovReshetikhinCrystal(['A',2,1],1,2)
                sage: T=TensorProductOfCrystals(B1,B2)
                sage: T1=TensorProductOfCrystals(B2,B1)
                sage: La = T.weight_lattice_realization().fundamental_weights()
                sage: t = [b for b in T if b.weight() == -3*La[0] + 3*La[1]][0]
                sage: t1 = [b for b in T1 if b.weight() == -3*La[0] + 3*La[1]][0]
                sage: g={t:t1}
                sage: f=T.crystal_morphism(g,acyclic = False)
                sage: [[b,f(b)] for b in T]
                [[[[[1]], [[1, 1]]], [[[1, 1]], [[1]]]],
                [[[[1]], [[1, 2]]], [[[1, 1]], [[2]]]],
                [[[[1]], [[2, 2]]], [[[1, 2]], [[2]]]],
                [[[[1]], [[1, 3]]], [[[1, 1]], [[3]]]],
                [[[[1]], [[2, 3]]], [[[1, 2]], [[3]]]],
                [[[[1]], [[3, 3]]], [[[1, 3]], [[3]]]],
                [[[[2]], [[1, 1]]], [[[1, 2]], [[1]]]],
                [[[[2]], [[1, 2]]], [[[2, 2]], [[1]]]],
                [[[[2]], [[2, 2]]], [[[2, 2]], [[2]]]],
                [[[[2]], [[1, 3]]], [[[2, 3]], [[1]]]],
                [[[[2]], [[2, 3]]], [[[2, 2]], [[3]]]],
                [[[[2]], [[3, 3]]], [[[2, 3]], [[3]]]],
                [[[[3]], [[1, 1]]], [[[1, 3]], [[1]]]],
                [[[[3]], [[1, 2]]], [[[1, 3]], [[2]]]],
                [[[[3]], [[2, 2]]], [[[2, 3]], [[2]]]],
                [[[[3]], [[1, 3]]], [[[3, 3]], [[1]]]],
                [[[[3]], [[2, 3]]], [[[3, 3]], [[2]]]],
                [[[[3]], [[3, 3]]], [[[3, 3]], [[3]]]]]
            """
            if index_set is None:
                index_set = self.index_set()
            if similarity_factor is None:
                similarity_factor = dict( (i,1) for i in index_set )
            if similarity_factor_domain is None:
                similarity_factor_domain = dict( (i,1) for i in index_set )
            if direction == 'down':
                e_string = 'e_string'
            else:
                e_string = 'f_string'
            if direction_image == 'down':
                f_string = 'f_string'
            else:
                f_string = 'e_string'

            if acyclic:
                if type(g) == dict:
                    g = g.__getitem__

                def morphism(b):
                    for i in index_set:
                        c = getattr(b, e_string)([i for k in range(similarity_factor_domain[i])])
                        if c is not None:
                            d = getattr(morphism(c), f_string)([automorphism(i) for k in range(similarity_factor[i])])
                            if d is not None:
                                return d
                            else:
                                raise ValueError, "This is not a morphism!"
                            #now we know that b is hw
                    return g(b)

            else:
                import copy
                morphism = copy.copy(g)
                known = set( g.keys() )
                todo = copy.copy(known)
                images = set( [g[x] for x in known] )
                # Invariants:
                # - images contains all known morphism(x)
                # - known contains all elements x for which we know morphism(x)
                # - todo  contains all elements x for which we haven't propagated to each child
                while todo <> set( [] ):
                    x = todo.pop()
                    for i in index_set:
                        eix  = getattr(x, f_string)([i for k in range(similarity_factor_domain[i])])
                        eigx = getattr(morphism[x], f_string)([automorphism(i) for k in range(similarity_factor[i])])
                        if bool(eix is None) <> bool(eigx is None):
                            # This is not a crystal morphism!
                            raise ValueError, "This is not a morphism!" #, print("x="x,"g(x)="g(x),"i="i)
                        if (eix is not None) and (eix not in known):
                            todo.add(eix)
                            known.add(eix)
                            morphism[eix] = eigx
                            images.add(eigx)
                # Check that the module generators are indeed module generators
                assert(len(known) == self.cardinality())
                # We may not want to systematically run those tests,
                # to allow for non bijective crystal morphism
                # Add an option CheckBijective?
                if not ( len(known) == len(images) and len(images) == images.pop().parent().cardinality() ):
                    return(None)
                return morphism.__getitem__

            if cached:
                return morphism
            else:
                return CachedFunction(morphism)

        def digraph(self):
            """
            Returns the DiGraph associated to self.

            EXAMPLES::

                sage: C = Crystals().example(5)
                sage: C.digraph()
                Digraph on 6 vertices

            TODO: add more tests
            """
            from sage.graphs.all import DiGraph
            d = {}
            for x in self:
                d[x] = {}
                for i in self.index_set():
                    child = x.f(i)
                    if child is None:
                        continue
                    d[x][child]=i
            return DiGraph(d)

        def latex_file(self, filename):
            r"""
            Exports a file, suitable for pdflatex, to 'filename'. This requires
            a proper installation of ``dot2tex`` in sage-python. For more
            information see the documentation for ``self.latex()``.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.latex_file('/tmp/test.tex') #optional requires dot2tex
            """
            header = r"""\documentclass{article}
            \usepackage[x11names, rgb]{xcolor}
            \usepackage[utf8]{inputenc}
            \usepackage{tikz}
            \usetikzlibrary{snakes,arrows,shapes}
            \usepackage{amsmath}
            \usepackage[active,tightpage]{preview}
            \newenvironment{bla}{}{}
            \PreviewEnvironment{bla}

            \begin{document}
            \begin{bla}"""

            footer = r"""\end{bla}
            \end{document}"""

            f = open(filename, 'w')
            f.write(header + self.latex() + footer)
            f.close()

        def latex(self):
            r"""
            Returns the crystal graph as a latex string. This can be exported
            to a file with self.latex_file('filename').

            See :meth:`Graph.layout_graphviz` for further documentation on which
            packages are required.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.latex()         #optional requires dot2tex
                ...
                sage: view(C, pdflatex = True, tightpage = True) # optional
            """
            try:
                from dot2tex.dot2tex import Dot2TikZConv
            except ImportError:
                print "dot2tex not available.  Install after running \'sage -sh\'"
                return

            # In MuPAD, 'autosize' is an option, but this doesn't seem to work here.
            options = {'format':'tikz', 'crop':True, 'usepdflatex':True, 'figonly':True}

            content = (Dot2TikZConv(options)).convert(self.dot_tex())

            return content

        _latex_ = latex

        def metapost(self, filename, thicklines=False, labels=True, scaling_factor=1.0, tallness=1.0):
            r"""
            Use C.metapost("filename.mp",[options]), where options can be:

            thicklines = True (for thicker edges) labels = False (to suppress
            labeling of the vertices) scaling_factor=value, where value is a
            floating point number, 1.0 by default. Increasing or decreasing the
            scaling factor changes the size of the image. tallness=1.0.
            Increasing makes the image taller without increasing the width.

            Root operators e(1) or f(1) move along red lines, e(2) or f(2)
            along green. The highest weight is in the lower left. Vertices with
            the same weight are kept close together. The concise labels on the
            nodes are strings introduced by Berenstein and Zelevinsky and
            Littelmann; see Littelmann's paper Cones, Crystals, Patterns,
            sections 5 and 6.

            For Cartan types B2 or C2, the pattern has the form

            a2 a3 a4 a1

            where c\*a2 = a3 = 2\*a4 =0 and a1=0, with c=2 for B2, c=1 for C2.
            Applying e(2) a1 times, e(1) a2 times, e(2) a3 times, e(1) a4 times
            returns to the highest weight. (Observe that Littelmann writes the
            roots in opposite of the usual order, so our e(1) is his e(2) for
            these Cartan types.) For type A2, the pattern has the form

            a3 a2 a1

            where applying e(1) a1 times, e(2) a2 times then e(3) a1 times
            returns to the highest weight. These data determine the vertex and
            may be translated into a Gelfand-Tsetlin pattern or tableau.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 2])
                sage: C.metapost('/tmp/test.mp') #optional

            ::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.metapost('/tmp/test.mp')
                Traceback (most recent call last):
                ...
                NotImplementedError
            """
            # FIXME: those tests are not robust
            # Should use instead self.cartan_type() == CartanType(['B',2])
            if self.cartan_type()[0] == 'B' and self.cartan_type()[1] == 2:
                word = [2,1,2,1]
            elif self.cartan_type()[0] == 'C' and self.cartan_type()[1] == 2:
                word = [2,1,2,1]
            elif self.cartan_type()[0] == 'A' and self.cartan_type()[1] == 2:
                word = [1,2,1]
            else:
                raise NotImplementedError
            size = self.cardinality()
            string_data = []
            for i in range(size):
                turtle = self.list()[i]
                string_datum = []
                for j in word:
                    turtlewalk = 0
                    while not turtle.e(j) == None:
                        turtle = turtle.e(j)
                        turtlewalk += 1
                    string_datum.append(turtlewalk)
                string_data.append(string_datum)

            if self.cartan_type()[0] == 'A':
                if labels:
                    c0 = int(55*scaling_factor)
                    c1 = int(-25*scaling_factor)
                    c2 = int(45*tallness*scaling_factor)
                    c3 = int(-12*scaling_factor)
                    c4 = int(-12*scaling_factor)
                else:
                    c0 = int(45*scaling_factor)
                    c1 = int(-20*scaling_factor)
                    c2 = int(35*tallness*scaling_factor)
                    c3 = int(12*scaling_factor)
                    c4 = int(-12*scaling_factor)
                outstring = "verbatimtex\n\\magnification=600\netex\n\nbeginfig(-1);\nsx:=35; sy:=30;\n\nz1000=(%d,0);\nz1001=(%d,%d);\nz1002=(%d,%d);\nz2001=(-3,3);\nz2002=(3,3);\nz2003=(0,-3);\nz2004=(7,0);\nz2005=(0,7);\nz2006=(-7,0);\nz2007=(0,7);\n\n"%(c0,c1,c2,c3,c4)
            else:
                if labels:
                    outstring = "verbatimtex\n\\magnification=600\netex\n\nbeginfig(-1);\n\nsx := %d;\nsy=%d;\n\nz1000=(2*sx,0);\nz1001=(-sx,sy);\nz1002=(-16,-10);\n\nz2001=(0,-3);\nz2002=(-5,3);\nz2003=(0,3);\nz2004=(5,3);\nz2005=(10,1);\nz2006=(0,10);\nz2007=(-10,1);\nz2008=(0,-8);\n\n"%(int(scaling_factor*40),int(tallness*scaling_factor*40))
                else:
                    outstring = "beginfig(-1);\n\nsx := %d;\nsy := %d;\n\nz1000=(2*sx,0);\nz1001=(-sx,sy);\nz1002=(-5,-5);\n\nz1003=(10,10);\n\n"%(int(scaling_factor*35),int(tallness*scaling_factor*35))
            for i in range(size):
                if self.cartan_type()[0] == 'A':
                    [a1,a2,a3] = string_data[i]
                else:
                    [a1,a2,a3,a4] = string_data[i]
                shift = 0
                for j in range(i):
                    if self.cartan_type()[0] == 'A':
                        [b1,b2,b3] = string_data[j]
                        if b1+b3 == a1+a3 and b2 == a2:
                            shift += 1
                    else:
                        [b1,b2,b3,b4] = string_data[j]
                        if b1+b3 == a1+a3 and b2+b4 == a2+a4:
                            shift += 1
                if self.cartan_type()[0] == 'A':
                    outstring = outstring +"z%d=%d*z1000+%d*z1001+%d*z1002;\n"%(i,a1+a3,a2,shift)
                else:
                    outstring = outstring +"z%d=%d*z1000+%d*z1001+%d*z1002;\n"%(i,a1+a3,a2+a4,shift)
            outstring = outstring + "\n"
            if thicklines:
                outstring = outstring +"pickup pencircle scaled 2\n\n"
            for i in range(size):
                for j in range(1,3):
                    dest = self.list()[i].f(j)
                    if not dest == None:
                        dest = self.list().index(dest)
                        if j == 1:
                            col = "red;"
                        else:
                            col = "green;  "
                        if self.cartan_type()[0] == 'A':
                            [a1,a2,a3] = string_data[i] # included to facilitate hand editing of the .mp file
                            outstring = outstring+"draw z%d--z%d withcolor %s   %% %d %d %d\n"%(i,dest,col,a1,a2,a3)
                        else:
                            [a1,a2,a3,a4] = string_data[i]
                            outstring = outstring+"draw z%d--z%d withcolor %s   %% %d %d %d %d\n"%(i,dest,col,a1,a2,a3,a4)
            outstring += "\npickup pencircle scaled 3;\n\n"
            for i in range(self.cardinality()):
                if labels:
                    if self.cartan_type()[0] == 'A':
                        outstring = outstring+"pickup pencircle scaled 15;\nfill z%d+z2004..z%d+z2006..z%d+z2006..z%d+z2007..cycle withcolor white;\nlabel(btex %d etex, z%d+z2001);\nlabel(btex %d etex, z%d+z2002);\nlabel(btex %d etex, z%d+z2003);\npickup pencircle scaled .5;\ndraw z%d+z2004..z%d+z2006..z%d+z2006..z%d+z2007..cycle;\n"%(i,i,i,i,string_data[i][2],i,string_data[i][1],i,string_data[i][0],i,i,i,i,i)
                    else:
                        outstring = outstring+"%%%d %d %d %d\npickup pencircle scaled 1;\nfill z%d+z2005..z%d+z2006..z%d+z2007..z%d+z2008..cycle withcolor white;\nlabel(btex %d etex, z%d+z2001);\nlabel(btex %d etex, z%d+z2002);\nlabel(btex %d etex, z%d+z2003);\nlabel(btex %d etex, z%d+z2004);\npickup pencircle scaled .5;\ndraw z%d+z2005..z%d+z2006..z%d+z2007..z%d+z2008..cycle;\n\n"%(string_data[i][0],string_data[i][1],string_data[i][2],string_data[i][3],i,i,i,i,string_data[i][0],i,string_data[i][1],i,string_data[i][2],i,string_data[i][3],i,i,i,i,i)
                else:
                    outstring += "drawdot z%d;\n"%i
            outstring += "\nendfig;\n\nend;\n\n"

            f = open(filename, 'w')
            f.write(outstring)
            f.close()

        def dot_tex(self):
            r"""
            Returns a dot_tex string representation of ``self``.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',2])
                sage: C.dot_tex()
                'digraph G { \n  node [ shape=plaintext ];\n  N_0 [ label = " ", texlbl = "$1$" ];\n  N_1 [ label = " ", texlbl = "$2$" ];\n  N_2 [ label = " ", texlbl = "$3$" ];\n  N_0 -> N_1 [ label = " ", texlbl = "1" ];\n  N_1 -> N_2 [ label = " ", texlbl = "2" ];\n}'
            """
            import re
            rank = ranker.from_list(self.list())[0]
            vertex_key = lambda x: "N_"+str(rank(x))

            # To do: check the regular expression
            # Removing %-style comments, newlines, quotes
            # This should probably be moved to sage.misc.latex
            quoted_latex = lambda x: re.sub("\"|\r|(%[^\n]*)?\n","", latex(x))

            result = "digraph G { \n  node [ shape=plaintext ];\n"

            for x in self:
                result += "  " + vertex_key(x) + " [ label = \" \", texlbl = \"$"+quoted_latex(x)+"$\" ];\n"
            for x in self:
                for i in self.index_set():
                    child = x.f(i)
                    if child is None:
                        continue
    #                result += "  " + vertex_key(x) + " -> "+vertex_key(child)+ " [ label = \" \", texlbl = \""+quoted_latex(i)+"\" ];\n"
                    if i == 0:
                        option = "dir = back, "
                        (source, target) = (child, x)
                    else:
                        option = ""
                        (source, target) = (x, child)
                    result += "  " + vertex_key(source) + " -> "+vertex_key(target)+ " [ "+option+"label = \" \", texlbl = \""+quoted_latex(i)+"\" ];\n"
            result+="}"
            return result

        def plot(self, **options):
            """
            Returns the plot of self as a directed graph.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: show_default(False) #do not show the plot by default
                sage: C.plot()
                Graphics object consisting of 17 graphics primitives
            """
            return self.digraph().plot(edge_labels=True,vertex_size=0,**options)

    class ElementMethods:

        def index_set(self):
            """
            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).index_set()
                [1, 2, 3, 4, 5]
            """
            return self.parent().index_set()

        def cartan_type(self):
            """
            Returns the cartan type associated to ``self``

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C(1).cartan_type()
                ['A', 5]
            """
            return self.parent().cartan_type()

        def weight(self):
            """
            Returns the weight of this crystal element

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).weight()
                (1, 0, 0, 0, 0, 0)
            """
            return self.Phi() - self.Epsilon()

        @abstract_method
        def e(self, i):
            r"""
            Returns `e_i(x)` if it exists or ``None`` otherwise.

            This method should be implemented by the element class of
            the crystal.

            EXAMPLES::

                sage: C = Crystals().example(5)
                sage: x = C[2]; x
                3
                sage: x.e(1), x.e(2), x.e(3)
                (None, 2, None)
            """

        @abstract_method
        def f(self, i):
            r"""
            Returns `f_i(x)` if it exists or ``None`` otherwise.

            This method should be implemented by the element class of
            the crystal.

            EXAMPLES::

                sage: C = Crystals().example(5)
                sage: x = C[1]; x
                2
                sage: x.f(1), x.f(2), x.f(3)
                (None, 3, None)
            """

        def epsilon(self, i):
            r"""
            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).epsilon(1)
                0
                sage: C(2).epsilon(1)
                1
            """
            assert i in self.index_set()
            x = self
            eps = 0
            while True:
                x = x.e(i)
                if x is None:
                    break
                eps = eps+1
            return eps

        def phi(self, i):
            r"""
            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).phi(1)
                1
                sage: C(2).phi(1)
                0
            """
            assert i in self.index_set()
            x = self
            phi = 0
            while True:
                x = x.f(i)
                if x is None:
                    break
                phi = phi+1
            return phi

        def phi_minus_epsilon(self, i):
            """
            Returns `\phi_i - \epsilon_i` of self. There are sometimes
            better implementations using the weight for this. It is used
            for reflections along a string.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).phi_minus_epsilon(1)
                1
            """
            return self.phi(i) - self.epsilon(i)

        def Epsilon(self):
            """
            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(0).Epsilon()
                (0, 0, 0, 0, 0, 0)
                sage: C(1).Epsilon()
                (0, 0, 0, 0, 0, 0)
                sage: C(2).Epsilon()
                (1, 0, 0, 0, 0, 0)
            """
            Lambda = self.parent().Lambda()
            return sum(self.epsilon(i) * Lambda[i] for i in self.index_set())

        def Phi(self):
            """
            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(0).Phi()
                (0, 0, 0, 0, 0, 0)
                sage: C(1).Phi()
                (1, 0, 0, 0, 0, 0)
                sage: C(2).Phi()
                (1, 1, 0, 0, 0, 0)
            """
            Lambda = self.parent().Lambda()
            return sum(self.phi(i) * Lambda[i] for i in self.index_set())

        def f_string(self, list):
            r"""
            Applies `f_{i_r} ... f_{i_1}` to self for `list = [i_1, ..., i_r]`

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',3])
                sage: b = C(1)
                sage: b.f_string([1,2])
                3
                sage: b.f_string([2,1])

            """
            b = self
            for i in list:
                b = b.f(i)
                if b is None:
                    return None
            return b

        def e_string(self, list):
            r"""
            Applies `e_{i_r} ... e_{i_1}` to self for `list = [i_1, ..., i_r]`

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',3])
                sage: b = C(3)
                sage: b.e_string([2,1])
                1
                sage: b.e_string([1,2])

            """
            b = self
            for i in list:
                b = b.e(i)
                if b is None:
                    return None
            return b

        def s(self, i):
            r"""
            Returns the reflection of ``self`` along its `i`-string

            EXAMPLES::

                sage: C = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: b=C(rows=[[1,1],[3]])
                sage: b.s(1)
                [[2, 2], [3]]
                sage: b=C(rows=[[1,2],[3]])
                sage: b.s(2)
                [[1, 2], [3]]
                sage: T=CrystalOfTableaux(['A',2],shape=[4])
                sage: t=T(rows=[[1,2,2,2]])
                sage: t.s(1)
                [[1, 1, 1, 2]]
            """
            d = self.phi_minus_epsilon(i)
            b = self
            if d > 0:
                for j in range(d):
                    b = b.f(i)
            else:
                for j in range(-d):
                    b = b.e(i)
            return b

        def demazure_operator(self, i, truncated = False):
            r"""
            Returns the list of the elements one can obtain from
            ``self`` by application of `f_i`.  If the option
            "truncated" is set to True, then ``self`` is not included
            in the list.

            REFERENCES:

                .. [L1995] Peter Littelmann, Crystal graphs and Young tableaux,
                   J. Algebra 175 (1995), no. 1, 65--87.

                .. [K1993] Masaki Kashiwara, The crystal base and Littelmann's refined Demazure character formula,
                   Duke Math. J. 71 (1993), no. 3, 839--858.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: t=T(rows=[[1,2],[2]])
                sage: t.demazure_operator(2)
                [[[1, 2], [2]], [[1, 3], [2]], [[1, 3], [3]]]
                sage: t.demazure_operator(2, truncated = True)
                [[[1, 3], [2]], [[1, 3], [3]]]
                sage: t.demazure_operator(1, truncated = True)
                []
                sage: t.demazure_operator(1)
                [[[1, 2], [2]]]

                sage: K = KirillovReshetikhinCrystal(['A',2,1],2,1)
                sage: t = K(rows=[[3],[2]])
                sage: t.demazure_operator(0)
                [[[2, 3]], [[1, 2]]]
            """
            if truncated:
                l = []
            else:
                l = [self]
            for k in range(self.phi(i)):
                l.append(self.f_string([i for j in range(k+1)]))
            return(l)

        def is_highest_weight(self, index_set = None):
            r"""
            Returns ``True`` if ``self`` is a highest weight.
            Specifying the option ``index_set`` to be a subset `I` of the
            index set of the underlying crystal, finds all highest
            weight vectors for arrows in `I`.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).is_highest_weight()
                True
                sage: C(2).is_highest_weight()
                False
                sage: C(2).is_highest_weight(index_set = [2,3,4,5])
                True
            """
            if index_set is None:
                index_set = self.index_set()
            return all(self.e(i) is None for i in index_set)
