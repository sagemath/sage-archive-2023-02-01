r"""
Crystals

Let `T` be a CartanType with index set `I`, and
`W` be a realization of the type `T` weight
lattice.

A type `T` crystal `C` is a colored oriented graph
equipped with a weight function from the nodes to some realization
of the type `T` weight lattice such that:


-  Each edge is colored with a label in `i \in I`.

-  For each `i\in I`, each node `x` has:


   -  at most one `i`-successor `f_i(x)`;

   -  at most one `i`-predecessor `e_i(x)`.


   Furthermore, when they exist,


   -  `f_i(x)`.weight() = x.weight() - `\alpha_i`;

   -  `e_i(x)`.weight() = x.weight() +
      `\alpha_i`.



This crystal actually models a representation of a Lie algebra if
it satisfies some further local conditions due to Stembridge, see
J. Stembridge, *A local characterization of simply-laced crystals*,
Trans. Amer. Math. Soc. 355 (2003), no. 12, 4807-4823.

EXAMPLES:

We construct the type `A_5` crystal on letters (or in
representation theoretic terms, the highest weight crystal of type
`A_5` corresponding to the highest weight
`\Lambda_1`)

::

    sage: C = CrystalOfLetters(['A',5]); C
    The crystal of letters for type ['A', 5]

It has a single highest weight element::

    sage: C.highest_weight_vectors()
    [1]

A crystal is a CombinatorialClass; and we can count and list its
elements in the usual way::

    sage: C.cardinality()
    6
    sage: C.list()
    [1, 2, 3, 4, 5, 6]

as well as use it in for loops

::

    sage: [x for x in C]
    [1, 2, 3, 4, 5, 6]

Here are some more elaborate crystals (see their respective
documentations)::

    sage: Tens = TensorProductOfCrystals(C, C)
    sage: Spin = CrystalOfSpins(['B', 3])
    sage: Tab  = CrystalOfTableaux(['A', 3], shape = [2,1,1])
    sage: Fast  = FastCrystal(['B', 2], shape = [3/2, 1/2])

One can get (currently) crude plotting via::

    sage: Tab.plot()

For rank two crystals, there is an alternative method of getting
metapost pictures. For more information see C.metapost?

Caveat: this crystal library, although relatively featureful for
classical crystals, is still in an early development stage, and the
syntax details may be subject to changes.

TODO:


-  Vocabulary and conventions:


   -  elements or vectors of a crystal?

   -  For a classical crystal: connected / highest weight /
      irreducible

   -  ...


-  More introductory doc explaining the mathematical background

-  Layout instructions for plot() for rank 2 types

-  Streamlining the latex output

-  Littelmann paths and/or alcove paths (this would give us the
   exceptional types)

-  RestrictionOfCrystal / DirectSumOfCrystals

-  Crystal.crystal_morphism

-  Affine crystals

-  Kirillov-Reshetikhin crystals


Most of the above features (except Littelmann/alcove paths) are in
MuPAD-Combinat (see lib/COMBINAT/crystals.mu), which could provide
inspiration.
"""

#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
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
#****************************************************************************
# Acknowledgment: most of the design and implementation of this
# library is heavily inspired from MuPAD-Combinat.
#****************************************************************************

from sage.misc.latex import latex
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.combinat.combinat import CombinatorialClass
from sage.graphs.graph import DiGraph
from sage.combinat import ranker
from sage.combinat.tools import transitive_ideal
from sage.combinat.root_system.weyl_characters import WeylCharacterRing, WeylCharacter
from sage.combinat.backtrack import GenericBacktracker

## MuPAD-Combinat's Cat::Crystal
# FIXME: crystals, like most parent should have unique data representation
class Crystal(CombinatorialClass, Parent):
    r"""
    The abstract class of crystals

    instances of this class should have the following attributes:


    -  cartan_type

    -  index_set the index set of the cartan type

    -  module_generators a list (or container) of distinct elements
       which generate the crystal using `f_i`

    -  weight_lattice_realization
    """

    def weight_lattice_realization(self):
        """
        Returns the weight lattice realization for the root system
        associated to self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A', 5])
            sage: C.weight_lattice_realization()
            Ambient space of the Root system of type ['A', 5]
        """
        return self.cartan_type.root_system().ambient_space()

    def Lambda(self):
        """
        Returns the fundamentals weights in the weight lattice realization
        for the root system associated to self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A', 5])
            sage: C.Lambda()
            Finite family {1: (1, 0, 0, 0, 0, 0), 2: (1, 1, 0, 0, 0, 0), 3: (1, 1, 1, 0, 0, 0), 4: (1, 1, 1, 1, 0, 0), 5: (1, 1, 1, 1, 1, 0)}
        """
        return self.weight_lattice_realization().fundamental_weights()

    def check(self):
        r"""
        Runs sanity checks on the crystal:


        -  Checks that count, list, and __iter__ are consistent. For a
           ClassicalCrystal, this in particular checks that the number of
           elements returned by the brute force listing and the iterator
           __iter__ are consistent with the Weyl dimension formula.

        -  Should check Stembridge's rules, etc.


        EXAMPLES::

            sage: C = CrystalOfLetters(['A', 5])
            sage: C.check()
            True
        """
        # Those tests could be lifted up to CombinatorialClass
        list1 = self.list()
        set1 = set(list1)
        list2 = [c for c in self]
        set2 = set(list2)
        list3 = Crystal.list(self)
        set3 = set(list3)
        return len(set1) == len(list1) \
               and len(set2) == len(list2) \
               and len(set3) == len(list3) \
               and len(set1) == self.cardinality() \
               and set1 == set2 \
               and set2 == set3

    def list(self):
        """
        Returns a list of the elements of self obtained by continually
        apply the `f_i` operators to the module generators of
        self.

        EXAMPLES::

            sage: from sage.combinat.crystals.crystals import Crystal
            sage: C = CrystalOfLetters(['A', 5])
            sage: l = Crystal.list(C)
            sage: l.sort(); l
            [1, 2, 3, 4, 5, 6]
        """
        # To be generalized to some transitiveIdeal
        # To be moved in a super category CombinatorialModule
        result = set(self.module_generators)
        todo = result.copy()
        while len(todo) > 0:
            x = todo.pop()
            for i in self.index_set:
                y = x.f(i)
                if y == None or y in result:
                    continue
                todo.add(y)
                result.add(y)
        return list(result)

    def digraph(self):
        """
        Returns the DiGraph associated to self.

        EXAMPLES::

            sage: from sage.combinat.crystals.crystals import Crystal
            sage: C = CrystalOfLetters(['A', 5])
            sage: Crystal.digraph(C)
            Digraph on 6 vertices
        """
        d = {}
        for x in self:
            d[x] = {}
            for i in self.index_set:
                child = x.f(i)
                if child is None:
                    continue
                d[x][child]=i
        return DiGraph(d)

    def character(self, R):
        """
        INPUT: R, a WeylCharacterRing. Produces the character of the
        crystal.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',2])
            sage: T = TensorProductOfCrystals(C, C)
            sage: A2 = WeylCharacterRing(C.cartan_type); A2
            The Weyl Character Ring of Type [A,2] with Integer Ring coefficients
            sage: chi = T.character(A2); chi
            A2(1,1,0) + A2(2,0,0)
            sage: chi.check(verbose = true)
            [9, 9]
        """
        if not R.cartan_type == self.cartan_type:
            raise ValueError, "ring does not have the right Cartan type"
        hlist = {}
        mlist = {}

	for x in self.highest_weight_vectors():
            k = x.weight()
            if k in hlist:
                hlist[k] += 1
            else:
                hlist[k] = 1
        for x in self.list():
            k = x.weight()
            if k in mlist:
                mlist[k] += 1
            else:
                mlist[k] = 1
        return WeylCharacter(R, hlist, mlist)

    def latex_file(self, filename):
        r"""
        Exports a file, suitable for pdflatex, to 'filename'. This requires
        a proper installation of dot2tex in sage-python. For more
        information see the documentation for self.latex().

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
        Returns the crystal graph as a bit of latex. This can be exported
        to a file with self.latex_file('filename').

        This requires dot2tex to be installed in sage-python.

        Here some tips for installation:


        -  Install graphviz = 2.14

        -  Download pyparsing-1.4.11.tar.gz pydot-0.9.10.tar.gz
           dot2tex-2.7.0.tar.gz (see the dot2tex web page for download links)
           (note that the most recent version of pydot may not work. Be sure
           to install the 0.9.10 version.) Install each of them using the
           standard python install, but using sage-python:

           ::

                           # FIX ACCORDING TO YOUR Sage INSTALL
                           export sagedir=/opt/sage/
                           export sagepython=$sagedir/local/bin/sage-python

                           # Use downloaded version nums
                           for package in pyparsing-1.4.11 pydot-0.9.10 dot2tex-2.7.0; do\
                                   tar zxvf $package.tar.gz;\
                                   cd $package;\
                                   sudo $sagepython setup.py install;\
                                   cd ..;\
                               done


        -  Install pgf-2.00 inside your latex tree In short:


           -  untaring in /usr/share/texmf/tex/generic

           -  clean out remaining pgf files from older version

           -  run texhash



        You should be done! To test, go to the dot2tex-2.7.0/examples
        directory, and type::

                    $sagedir//local/bin/dot2tex balls.dot > balls.tex
                    pdflatex balls.tex
                    open balls.pdf \#your favorite viewer here


        EXAMPLES::

            sage: C = CrystalOfLetters(['A', 5])
            sage: C.latex() #optional requires dot2tex
            ...
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

    def metapost(self, filename, thicklines=False, labels=True, scaling_factor=1.0, tallness=1.0):
        """Use C.metapost("filename.mp",[options])
        where options can be:

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
        if self.cartan_type[0] == 'B' and self.cartan_type[1] == 2:
            word = [2,1,2,1]
        elif self.cartan_type[0] == 'C' and self.cartan_type[1] == 2:
            word = [2,1,2,1]
        elif self.cartan_type[0] == 'A' and self.cartan_type[1] == 2:
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

        if self.cartan_type[0] == 'A':
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
            if self.cartan_type[0] == 'A':
                [a1,a2,a3] = string_data[i]
            else:
                [a1,a2,a3,a4] = string_data[i]
            shift = 0
            for j in range(i):
                if self.cartan_type[0] == 'A':
                    [b1,b2,b3] = string_data[j]
                    if b1+b3 == a1+a3 and b2 == a2:
                        shift += 1
                else:
                    [b1,b2,b3,b4] = string_data[j]
                    if b1+b3 == a1+a3 and b2+b4 == a2+a4:
                        shift += 1
            if self.cartan_type[0] == 'A':
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
                    if self.cartan_type[0] == 'A':
                        [a1,a2,a3] = string_data[i] # included to facilitate hand editing of the .mp file
                        outstring = outstring+"draw z%d--z%d withcolor %s   %% %d %d %d\n"%(i,dest,col,a1,a2,a3)
                    else:
                        [a1,a2,a3,a4] = string_data[i]
                        outstring = outstring+"draw z%d--z%d withcolor %s   %% %d %d %d %d\n"%(i,dest,col,a1,a2,a3,a4)
        outstring += "\npickup pencircle scaled 3;\n\n"
        for i in range(self.cardinality()):
            if labels:
                if self.cartan_type[0] == 'A':
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
        Returns a dot_tex version of self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',2])
            sage: C.dot_tex()
            'digraph G { \n  node [ shape=plaintext ];\n  N_0 [ label = " ", texlbl = "$\\text{1}$" ];\n  N_1 [ label = " ", texlbl = "$\\text{2}$" ];\n  N_2 [ label = " ", texlbl = "$\\text{3}$" ];\n  N_0 -> N_1 [ label = " ", texlbl = "1" ];\n  N_1 -> N_2 [ label = " ", texlbl = "2" ];\n}'
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
            for i in self.index_set:
                child = x.f(i)
                if child is None:
                    continue
                result += "  " + vertex_key(x) + " -> "+vertex_key(child)+ " [ label = \" \", texlbl = \""+quoted_latex(i)+"\" ];\n"
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

class CrystalElement(Element):
    r"""
    The abstract class of crystal elements

    Sub classes should implement:


    -  x.e(i) (returning `e_i(x)`)

    -  x.f(i) (returning `f_i(x)`)

    -  x.weight()
    """

    def index_set(self):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: C(1).index_set()
            [1, 2, 3, 4, 5]
        """
        return self._parent.index_set

    def weight(self):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: C(1).weight()
            (1, 0, 0, 0, 0, 0)
        """
	return self.Phi() - self.Epsilon()

    def e(self, i):
        r"""
        Returns `e_i(x)` if it exists or None otherwise. This is
        to be implemented by subclasses of CrystalElement.

        TESTS::

            sage: from sage.combinat.crystals.crystals import CrystalElement
            sage: C = CrystalOfLetters(['A',5])
            sage: CrystalElement.e(C(1), 1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def f(self, i):
        r"""
        Returns `f_i(x)` if it exists or None otherwise. This is
        to be implemented by subclasses of CrystalElement.

        TESTS::

            sage: from sage.combinat.crystals.crystals import CrystalElement
            sage: C = CrystalOfLetters(['A',5])
            sage: CrystalElement.f(C(1), 1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

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
	return sum(self.epsilon(i) * self._parent.Lambda()[i] for i in self.index_set())

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
	return sum(self.phi(i) * self._parent.Lambda()[i] for i in self.index_set())

    def s(self, i):
	r"""
	Returns the reflection of self along its `i`-string

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
        d = self.phi(i)-self.epsilon(i)
	b = self
	if d > 0:
	    for j in range(d):
		b = b.f(i)
    	else:
	    for j in range(-d):
		b = b.e(i)
	return b

    def is_highest_weight(self):
	r"""
	Returns True if self is a highest weight.

	EXAMPLES::

	    sage: C = CrystalOfLetters(['A',5])
	    sage: C(1).is_highest_weight()
	    True
	    sage: C(2).is_highest_weight()
	    False
	"""
	return all(self.e(i) == None for i in self.index_set())

class CrystalBacktracker(GenericBacktracker):
    def __init__(self, crystal):
        """
        Time complexity: `O(nf)` amortized for each produced
        element, where `n` is the size of the index set, and f is
        the cost of computing `e` and `f` operators.

        Memory complexity: O(depth of the crystal)

        Principle of the algorithm:

        Let C be a classical crystal. It's an acyclic graph where all
        connected component has a unique element without predecessors (the
        highest weight element for this component). Let's assume for
        simplicity that C is irreducible (i.e. connected) with highest
        weight element u.

        One can define a natural spanning tree of `C` by taking
        `u` as rot of the tree, and for any other element
        `y` taking as ancestor the element `x` such that
        there is an `i`-arrow from `x` to `y` with
        `i` minimal. Then, a path from `u` to `y`
        describes the lexicographically smallest sequence
        `i_1,\dots,i_k` such that
        `(f_{i_k} \circ f_{i_1})(u)=y`.

        Morally, the iterator implemented below just does a depth first
        search walk through this spanning tree. In practice, this can be
        achieved recursively as follow: take an element `x`, and
        consider in turn each successor `y = f_i(x)`, ignoring
        those such that `y = f_j(x')` for some `x'` and
        `j<i` (this can be tested by computing `e_j(y)`
        for `j<i`).

        EXAMPLES::

            sage: from sage.combinat.crystals.crystals import CrystalBacktracker
            sage: C = CrystalOfTableaux(['B',3],shape=[3,2,1])
            sage: CB = CrystalBacktracker(C)
            sage: len(list(CB))
            1617
        """
        GenericBacktracker.__init__(self, None, None)
        self._crystal = crystal

    def _rec(self, x, state):
        """
        Returns an iterator for the (immediate) children of x in the search
        tree.

        EXAMPLES::

            sage: from sage.combinat.crystals.crystals import CrystalBacktracker
            sage: C = CrystalOfLetters(['A', 5])
            sage: CB = CrystalBacktracker(C)
            sage: list(CB._rec(C(1), 'n/a'))
            [(2, 'n/a', True)]
        """
        #We will signal the initial case by having a object and state
        #of None and consider it separately.
        if x is None and state is None:
            for gen in self._crystal.highest_weight_vectors():
                yield gen, "n/a", True
            return

        # Run through the children y of x
        for i in self._crystal.index_set:
            y = x.f(i)
            if y is None:
                continue
            # Ignore those which can be reached by an arrow with smaller label
            hasParent = False
            for j in x.index_set():
                if j == i:
                    break
                if not y.e(j) is None:
                    hasParent = True
                    break
            if hasParent:
                continue

            # yield y and all elements further below
            yield y, "n/a", True



class ClassicalCrystal(Crystal):
    r"""
    The abstract class of classical crystals
    """
    list  = CombinatorialClass.list#__list_from_iterator
    def __iter__(self):
        r"""
        Returns an iterator over the elements of the crystal.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: [x for x in C]
            [1, 2, 3, 4, 5, 6]

        TESTS::

            sage: C = CrystalOfLetters(['D',4])
            sage: D = CrystalOfSpinsPlus(['D',4])
            sage: E = CrystalOfSpinsMinus(['D',4])
            sage: T=TensorProductOfCrystals(D,E,generators=[[D.list()[0],E.list()[0]]])
            sage: U=TensorProductOfCrystals(C,E,generators=[[C(1),E.list()[0]]])
            sage: T.cardinality()
            56
            sage: T.check()
            True
            sage: U.check()
            True

        Bump's systematic tests::

            sage: fa3 = lambda a,b,c: CrystalOfTableaux(['A',3],shape=[a+b+c,b+c,c])
            sage: fb3 = lambda a,b,c: CrystalOfTableaux(['B',3],shape=[a+b+c,b+c,c])
            sage: fc3 = lambda a,b,c: CrystalOfTableaux(['C',3],shape=[a+b+c,b+c,c])
            sage: fb4 = lambda a,b,c,d: CrystalOfTableaux(['B',4],shape=[a+b+c+d,b+c+d,c+d,d])
            sage: fd4 = lambda a,b,c,d: CrystalOfTableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d])
            sage: fd5 = lambda a,b,c,d,e: CrystalOfTableaux(['D',5],shape=[a+b+c+d+e,b+c+d+e,c+d+e,d+e,e])
            sage: def fd4spinplus(a,b,c,d):\
                 C = CrystalOfTableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d]);\
                 D = CrystalOfSpinsPlus(['D',4]);\
                 return TensorProductOfCrystals(C,D,generators=[[C[0],D[0]]])
            sage: def fb3spin(a,b,c):\
                 C = CrystalOfTableaux(['B',3],shape=[a+b+c,b+c,c]);\
                 D = CrystalOfSpins(['B',3]);\
                 return TensorProductOfCrystals(C,D,generators=[[C[0],D[0]]])

        TODO: choose a good panel of values for a,b,c ... both for basic
        systematic tests and for conditionally run more computationally
        involved tests

        ::

            sage: fb4(1,0,1,0).check()
            True

        ::

            #sage: fb4(1,1,1,1).check() # expensive: the crystal is of size 297297
            #True
        """
        return iter(CrystalBacktracker(self))


    def highest_weight_vectors(self):
        r"""
        Returns a list of the highest weight vectors of self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: C.highest_weight_vectors()
            [1]

        ::

            sage: C = CrystalOfLetters(['A',2])
            sage: T = TensorProductOfCrystals(C,C,C,generators=[[C(2),C(1),C(1)],[C(1),C(2),C(1)]])
            sage: T.highest_weight_vectors()
            [[2, 1, 1], [1, 2, 1]]
        """
        # Implementation: selects among the module generators those that are
        # highest weight and cache the result
        try:
            return self._highest_weight_vectors
        except AttributeError:
            pass

        self._highest_weight_vectors = [g for g in self.module_generators if g.is_highest_weight()]
        return self._highest_weight_vectors

    def highest_weight_vector(self):
        r"""
        Returns the highest weight vector if there is a single one;
        otherwise, raise an error.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: C.highest_weight_vector()
            1
        """
        hw = self.highest_weight_vectors();
        if len(hw) == 1:
            return hw[0]
        else:
            raise RuntimeError("The crystal does not have exactly one highest weight vector")

    def cardinality(self):
        r"""
        Returns the number of elements of the crystal, using Weyl's
        dimension formula on each connected component

        EXAMPLES::

            sage: from sage.combinat.crystals.crystals import ClassicalCrystal
            sage: C = CrystalOfLetters(['A', 5])
            sage: ClassicalCrystal.cardinality(C)
            6
        """
        return sum(self.weight_lattice_realization().weyl_dimension(x.weight())
                   for x in self.highest_weight_vectors())


class AffineCrystal(Crystal):
    r"""
    The abstract class of affine crystals
    """
    pass
