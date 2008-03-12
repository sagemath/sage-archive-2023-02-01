r"""
Crystals

Let $T$ be a CartanType with index set $I$, and $W$ be a realization of the
type $T$ weight lattice.

A type $T$ crystal $C$ is an oriented graph equipped with a weight
function the nodes to some realization of the type $T$ weight lattice
such that:
\begin{itemize}
\item each edge has a label in $I$
\item for each $i$ in $I$, each node $x$ has:
    - at most one $i$-successor $f_i(x)$
    - at most one $i$-predecessor $e_i(x)$
   Furthermore, when the exists,
    - $f_i(x)$.weight() = x.weight() - $\alpha_i$
    - $e_i(x)$.weight() = x.weight() + $\alpha_i$

This crystal actually models a representation of a Lie algebra if it
satisfies some further local conditions due to Stembridge.

EXAMPLES:

We construct the type $A_5$ crystal on letters

    sage: C = CrystalOfLetters(['A',5]); C
    The crystal of letters for type ['A', 5]

It has a single highest weight element:
    sage: C.module_generators
    [1]

A crystal is a CombinatorialClass; and we can count and list its elements
in the usual way:
    sage: C.count()    # todo: not implemented
    6
    sage: C.list()
    [1, 2, 3, 4, 5, 6]

as well as use it in for loops
    sage: [x for x in C]
    [1, 2, 3, 4, 5, 6]

Here are some more elaborate crystals (see their respective documentations):
    sage: Tens = TensorProductOfCrystals(C, C)
    sage: Spin = CrystalOfSpins(['B', 3])
    sage: Tab  = CrystalOfTableaux(['A', 3], shape = [2,1,1])
    sage: Fast  = FastCrystal(['B', 2], shape = [3/2, 1/2])

One can get (currently) crude ploting via:

    sage: Tab.plot()

Thanks to graphviz and dot2tex (which need to be installed), one can produce
high quality LaTeX output of the crystal graph:

    f=open('/tmp/crystal.tex', 'w')
    f.write(C.latex())
    f.close()


Here some tips for installing dot2tex within sage:
 - Install graphviz >= 2.14
 - Make sure sage-python is >= 2.4
 - Download pyparsing-1.4.10.tgz pydot-0.9.10.tgz dot2tex-2.7.0.tgz
   (see the dot2tex web page for download links)
   Install each of them using the standard python install, but using sage-python:

        sagedir=/opt/sage-2.10-opteron-ubuntu64-x86_64-Linux/  # FIX ACCORDING TO YOUR SAGE INSTALL
        sagepython=$sagedir/local/bin/sage-python
	for package in pyparsing-1.4.10 pydot-0.9.10 dot2tex-2.7.0; do\
            tar zxvf $package.tgz;\
            cd $package;\
            sudo $sagepython setup.py install;\
            cd ..;\
        done

 - Install pgf-2.00 inside your latex tree
   In short:
    - untaring in in /usr/share/texmf/tex/generic
    - clean out remaining pgf files from older version
    - run texhash

You should be done!
To test, go to the dot2tex-2.7.0/examples directory, and type:

	$sagedir//local/bin/dot2tex balls.dot > balls.tex
        pdflatex balls.tex
        open balls.pdf # your favorite viewer here



For rank two crystals, there is an alternative method of getting
metapost pictures. For more information see C.metapost?

Caveat: this crystal library, although relatively featureful for
classical crystals, is still in an early development stage, and the
syntax details may be subject to changes.

TODO:
 - Vocabulary and conventions:
   - elements or vectors of a crystal?
   - For a classical crystal: connected / highest weight / irreducible
   - ...
 - More introductory doc explaining the mathematical background
 - Layout instructions for plot() for rank 2 types
 - Streamlining the latex output
 - Littelmann paths and/or alcove paths (this would give us the exceptional types)
 - RestrictionOfCrystal / DirectSumOfCrystals
 - Crystal.crystal_morphism
 - Affine crystals

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

from sage.misc.latex           import latex
from sage.structure.parent     import Parent
from sage.structure.element    import Element
from sage.combinat.combinat    import CombinatorialClass
from sage.combinat.cartan_type import CartanType
from sage.graphs.graph         import DiGraph
from sage.combinat             import ranker

## MuPAD-Combinat's Cat::Crystal
# FIXME: crystals, like most parent should have unique data representation
class Crystal(CombinatorialClass, Parent):
    r"""
    The abstract class of crystals

    instances of this class should have the following attributes:
    \begin{itemize}
    \item cartan_type
    \item index_set
        the index set of the cartan type
    \item module_generators
        a list (or container) of distinct elements which generate the crystal
    \item weight_lattice_realization
    \end{itemize}
    """

    def weight_lattice_realization(self):
	return self.cartanType.root_system().ambient_lattice()

    def Lambda(self):
	return self.weight_lattice_realization().fundamental_weights()

    def check(self):
        r"""
        Runs sanity checks on the crystal:
        \begin{itemize}
        \item Checks that count, list, and __iter__ are
        consistent. For a ClassicalCrystal, this in particular checks
        that the number of elements returned by the brute force
        listing and the iterator __iter__ are consistent with the Weyl
        dimension formula.
        \item Should check Stembridge's rules, etc.
        \end{itemize}
        """
        # Those tests could be lifted up to CombinatorialClass
        list1 = self.list();        set1 = set(list1)
        list2 = [ c for c in self]; set2 = set(list2)
        list3 = Crystal.list(self); set3 = set(list3)
        return len(set1) == len(list1)   and \
               len(set2) == len(list2)   and \
               len(set3) == len(list3)   and \
               len(set1) == self.count() and \
               set1 == set2              and \
               set2 == set3

    def list(self):
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
        dict = {}
        for x in self:
            dict[x] = {}
            for i in self.index_set:
                child = x.f(i)
                if child is None:
                    continue
                dict[x][child]=i
        return DiGraph(dict)

    def latex(self):
        try:
            from dot2tex.dot2tex import Dot2TikZConv
        except ImportError:
            print "dot2tex not available.  Install after running \'sage -sh\'"
            return "Game over, man!"
        conv = Dot2TikZConv()
        return conv.convert(self.dot_tex())

    def metapost(self, filename, thicklines=False, labels=True, scaling_factor=1.0, tallness=1.0):
        """Use C.metapost("filename.mp",[options])
        where options can be:

            thicklines = True (for thicker edges)
            labels = False (to suppress labeling of the vertices)
            scaling_factor=value, where value is a floating point
              number, 1.0 by default. Increasing or decreasing the
              scaling factor changes the size of the image.
            tallness=1.0. Increasing makes the image taller without
              increasing the width.

        Root operators e(1) or f(1) move along red lines, e(2) or
        f(2) along green. The highest weight is in the lower left.
        Vertices with the same weight are kept close together.
        The concise labels on the nodes are strings introduced by
        Berenstein and Zelevinsky and Littelmann; see Littelmann's
        paper Cones, Crystals, Patterns, sections 5 and 6.

        For Cartan types B2 or C2, the pattern has the form

              a2  a3  a4
                  a1

        where c*a2 >= a3 >= 2*a4 >=0 and a1>=0, with c=2 for B2, c=1 for C2.
        Applying e(2) a1 times, e(1) a2 times, e(2) a3 times, e(1) a4 times
        returns to the highest weight. (Observe that Littelmann writes the
        roots in opposite of the usual order, so our e(1) is his e(2) for
        these Cartan types.) For type A2, the pattern has the form

              a3   a2
                 a1

        where applying e(1) a1 times, e(2) a2 times then e(3) a1 times
        returns to the highest weight. These data determine the vertex
        and may be translated into a Gelfand-Tsetlin pattern or tableau.
        """
        if self.cartanType[0] == 'B' and self.cartanType[1] == 2:
            word = [2,1,2,1]
        elif self.cartanType[0] == 'C' and self.cartanType[1] == 2:
            word = [2,1,2,1]
        elif self.cartanType[0] == 'A' and self.cartanType[1] == 2:
            word = [1,2,1]
        else:
            raise NotImplementedError
        size = self.count()
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

        if self.cartanType[0] == 'A':
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
            if self.cartanType[0] == 'A':
                [a1,a2,a3] = string_data[i]
            else:
                [a1,a2,a3,a4] = string_data[i]
            shift = 0
            for j in range(i):
                if self.cartanType[0] == 'A':
                    [b1,b2,b3] = string_data[j]
                    if b1+b3 == a1+a3 and b2 == a2:
                        shift += 1
                else:
                    [b1,b2,b3,b4] = string_data[j]
                    if b1+b3 == a1+a3 and b2+b4 == a2+a4:
                        shift += 1
            if self.cartanType[0] == 'A':
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
                    if self.cartanType[0] == 'A':
                        [a1,a2,a3] = string_data[i] # included to facilitate hand editing of the .mp file
                        outstring = outstring+"draw z%d--z%d withcolor %s   %% %d %d %d\n"%(i,dest,col,a1,a2,a3)
                    else:
                        [a1,a2,a3,a4] = string_data[i]
                        outstring = outstring+"draw z%d--z%d withcolor %s   %% %d %d %d %d\n"%(i,dest,col,a1,a2,a3,a4)
        outstring += "\npickup pencircle scaled 3;\n\n"
        for i in range(self.count()):
            if labels:
                if self.cartanType[0] == 'A':
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
        rank = ranker.from_list(self.list())[0]
        def vertex_key(x):
            return "N_"+str(rank(x))
        def quoted_latex(x):
            import re
            # To do: check the regular expression
            # Removing %-style comments, newlines, quotes
            return re.sub("\"|\r|(%[^\n]*)?\n","", latex(x))

        result = "digraph G {\n"
        for x in self:
            result += vertex_key(x) + " [ label = \" \", texlbl = \"$"+quoted_latex(x)+"$\" ];\n"
        for x in self:
            for i in self.index_set:
                child = x.f(i)
                if child is None:
                    continue
                result += vertex_key(x)+ " -> "+vertex_key(child)+ " [ label = \" \", texlbl = \""+quoted_latex(i)+"\" ];\n"
        result+="}"
        return result

    def plot(self, **options):
        return self.digraph().plot(edge_labels=True,vertex_size=0,**options)

class CrystalElement(Element):
    r"""
    The abstract class of crystal elements

    Sub classes should implement:
    \begin{itemize}
    \item x.e(i)        (returning $e_i(x)$)
    \item x.f(i)        (returning $f_i(x)$)
    \item x.weight()
    \end{itemize}
    """

    def index_set(self):
        return self._parent.index_set

    def weight(self):
	return self.Phi() - self.Epsilon()

    def e(self, i):
        r"""
        Returns $e_i(x)$ if it exists or None otherwise
        """
        raise NotImplementedError

    def f(self, i):
        r"""
        Returns $f_i(x)$ if it exists or None otherwise
        """
        raise NotImplementedError

    def epsilon(self, i):
        r"""
        TESTS:
            # rather minimal tests
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
        TESTS:
            # rather minimal tests
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
	return sum(self.epsilon(i) * self._parent.Lambda()[i-1] for i in self.index_set())

    def Phi(self):
	return sum(self.phi(i) * self._parent.Lambda()[i-1] for i in self.index_set())

    def is_highest_weight(self):
	r"""
        TEST:
	    sage: C = CrystalOfLetters(['A',5])
	    sage: C(1).is_highest_weight()
	    True
	    sage: C(2).is_highest_weight()
	    False
	"""
	return all(self.e(i) == None for i in self.index_set())

class ClassicalCrystal(Crystal):
    r"""
    The abstract class of classical crystals
    """

    list  = CombinatorialClass.list#__list_from_iterator
    def __iter__(self):
        r"""
        Returns an iterator over the elements of the crystal.

        Time complexity: $O(nf)$ amortized for each produced element,
        where $n$ is the size of the index set, and f is the cost of
        computing $e$ and $f$ operators.

        Memory complexity: O(depth of the crystal)

        Principle of the algorithm:

        Let C be a classical crystal. It's an acyclic graph where all
        connected componnent has a unique element without predecessors
        (the highest weight element for this component). Let's assume
        for simplicity that C is irreducible (i.e. connected) with
        highest weigth element u.

        One can define a natural spanning tree of $C$ by taking $u$ as
        rot of the tree, and for any other element $y$ taking as
        ancestor the element $x$ such that there is an $i$-arrow from
        $x$ to $y$ with $i$ minimal. Then, a path from $u$ to $y$
        describes the lexicographically smallest sequence
        $i_1,\dots,i_k$ such that $(f_{i_k} \circ f_{i_1})(u)=y$.

        Morally, the iterator implemented below just does a depth
        first search walk through this spanning tree. In practice,
        this can be achieved recursively as follow: take an element
        $x$, and consider in turn each successor $y = f_i(x)$,
        ignoring those such that $y = f_j(x')$ for some $x'$ and $j<i$
        (this can be tested by computing $e_j(y)$ for $j<i$).

        It probably would be more efficient (and about as readable) to
        unroll the recursion by managing the stack by hand. This would
        avoid yieding down through the whole call stack for each
        element of the crystal.

        EXAMPLES:
            sage: C = CrystalOfLetters(['A',5])
            sage: [x for x in C]
            [1, 2, 3, 4, 5, 6]


        TESTS:
            sage: C = CrystalOfLetters(['D',4])
            sage: D = CrystalOfSpinsPlus(['D',4])
            sage: E = CrystalOfSpinsMinus(['D',4])
            sage: T=TensorProductOfCrystals(D,E,generators=[[D.list()[0],E.list()[0]]])
            sage: U=TensorProductOfCrystals(C,E,generators=[[C(1),E.list()[0]]])
            sage: len(T)    # triggered bug reported by Daniel Bump Sun, 24 Feb 2008 17:06:37 -0800
            56
            sage: T.check()
            True
            sage: U.check()
            True

            Bump's systematic tests:

            sage: def fa3(a,b,c):\
                 return CrystalOfTableaux(['A',3],shape=[a+b+c,b+c,c])
            sage: def fb3(a,b,c):\
                 return CrystalOfTableaux(['B',3],shape=[a+b+c,b+c,c])
            sage: def fb3spin(a,b,c):\
                 C = CrystalOfTableaux(['B',3],shape=[a+b+c,b+c,c]);\
                 D = CrystalOfSpins(['B',3]);\
                 return TensorProductOfCrystals(C,D,generators=[[C.list()[0],D.list()[0]]])
            sage: def fb4(a,b,c,d):\
                 return CrystalOfTableaux(['B',4],shape=[a+b+c+d,b+c+d,c+d,d])
            sage: def fc3(a,b,c):\
                 return CrystalOfTableaux(['C',3],shape=[a+b+c,b+c,c])
            sage: def fd4(a,b,c,d):\
                 return CrystalOfTableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d])
            sage: def fd4spinplus(a,b,c,d):\
                 C = CrystalOfTableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d]);\
                 D = CrystalOfSpinsPlus(['D',4]);\
                 return TensorProductOfCrystals(C,D,generators=[[C.list()[0],D.list()[0]]])
            sage: def fd5(a,b,c,d,e):\
                 return CrystalOfTableaux(['D',5],shape=[a+b+c+d+e,b+c+d+e,c+d+e,d+e,e])


            TODO: choose a good panel of values for a,b,c ... both for
            basic systematic tests and for conditionally run more
            computatinaly involved tests

            sage: fb4(1,0,1,0).check()
            True

            #sage: fb4(1,1,1,1).check() # expensive: the crystal is of size 297297
            #True


        """
        def rec(x):
            for i in self.index_set: # Run through the children y of x
                y = x.f(i)
                if y is None:
                    continue
                # Ignore those which can be reached by an arrow with smaller label
                hasParent = False
                for j in x.index_set():
                    if j == i:
                        break
                    if not y.e(j) == None:
                        hasParent = True
                        break
                if hasParent:
                    continue
                # yield y and all elements further below
                yield y
                for z in rec(y):
                    yield z

        for generator in self.highest_weight_vectors():
            yield generator
            for x in rec(generator):
                yield x

    iterator = __iter__

    def highest_weight_vectors(self):
        r"""
        Returns the highest weight vectors
        """
        # Implementation: selects among the module generators those that are highest weight
        # and cache the result
        if not self.__dict__.has_key('_highest_weight_vectors'): # What's the right idiom for testing the existence of an attribute
            self._highest_weight_vectors = []
            for x in self.module_generators: # What's the right idiom for 'select'
                if x.is_highest_weight():
                    self._highest_weight_vectors.append(x)
        return self._highest_weight_vectors

    def highest_weight_vector(self):
        r"""
        Returns the higest weight vector if there is a single one
        Raise an error otherwise
        """
        hw = self.highest_weight_vectors();
        if len(hw) == 1:
            return hw[0]
        else:
            raise RuntimeError("The crystal does not have exactly one highest weight vector")

    def count(self):
        r"""
        Returns the number of elements of the crystal, using Weyl's dimension formula on each
        connected component
        """
        return sum(self.weight_lattice_realization().weyl_dimension(x.weight())
                   for x in self.highest_weight_vectors())

class AffineCrystal(Crystal):
    r"""
    The abstract class of affine crystals
    """
    pass
