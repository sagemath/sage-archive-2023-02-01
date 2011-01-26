r"""
Substitutions over unit cube faces (Rauzy fractals)

This module implements the `E_1^*(\sigma)` substitution
associated with a one-dimensional substitution `\sigma`,
that acts on unit faces of dimension `(d-1)` in `\RR^d`.

This module defines the following classes and functions:

- ``Face`` - a class to model a face

- ``Patch`` - a class to model a finite set of faces

- ``E1Star`` - a class to model the `E_1^*(\sigma)` application
  defined by the substitution sigma

See the documentation of these objects for more information.

The convention for the choice of the unit faces and the
definition of `E_1^*(\sigma)` varies from article to article.
Here, unit faces are defined by

.. MATH::

    \begin{array}{ccc}
    \,[x, 1]^* & = & \{x + \lambda e_2 + \mu e_3 : \lambda, \mu \in [0,1]\} \\
    \,[x, 2]^* & = & \{x + \lambda e_1 + \mu e_3 : \lambda, \mu \in [0,1]\} \\
    \,[x, 3]^* & = & \{x + \lambda e_1 + \mu e_2 : \lambda, \mu \in [0,1]\}
    \end{array}

and the generalized substitution `E_1^*(\sigma)` is defined by

.. MATH::

    E_1^*(\sigma)([x,i]^*) =
    \bigcup_{k = 1,2,3} \; \bigcup_{s | \sigma(k) = pis}
    [M^{-1}(x + \ell(s)), k]^*,

where `\ell(s)` is the abelianized of `s`, and `M` is the matrix of `\sigma`.

AUTHORS:

- Franco Saliola (2009): initial version
- Vincent Delecroix, Timo Jolivet, Stepan Starosta (2010-05): redesign
- Timo Jolivet (2010-08, 2010-09): redesign

REFERENCES:

.. [AI] P. Arnoux, S. Ito,
   Pisot substitutions and Rauzy fractals,
   Bull. Belg. Math. Soc. 8 (2), 2001, pp. 181--207

.. [AIS] Y. Sano, P. Arnoux, S. Ito,
   Higher dimensional extensions of substitutions and their dual maps,
   J. Anal. Math. 83, 2001, pp. 183--206

EXAMPLES:

We start by drawing a simple three-face patch::

    sage: from sage.combinat.e_one_star import E1Star, Face, Patch
    sage: x = [Face((0,0,0),1), Face((0,0,0),2), Face((0,0,0),3)]
    sage: P = Patch(x)
    sage: P
    Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]
    sage: P.plot()                   #not tested

We apply a substitution to this patch, and draw the result::

    sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
    sage: E = E1Star(sigma)
    sage: E(P)
    Patch: [[(1, 0, -1), 1]*, [(0, 1, -1), 2]*, [(0, 0, 0), 3]*, [(0, 0, 0), 1]*, [(0, 0, 0), 2]*]
    sage: E(P).plot()                #not tested

.. NOTE::

    - The type of a face is given by an integer in ``[1, ..., d]``
      where ``d`` is the length of the vector of the face.

    - The alphabet of the domain and the codomain of `\sigma` must be
      equal, and they must be of the form ``[1, ..., d]``, where ``d``
      is a positive integer corresponding to the length of the vectors
      of the faces on which `E_1^*(\sigma)` will act.

::

    sage: x = [Face((0,0,0),1), Face((0,0,0),2), Face((0,0,0),3)]
    sage: P = Patch(x)
    sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
    sage: E = E1Star(sigma)
    sage: E(P)
    Patch: [[(1, 0, -1), 1]*, [(0, 1, -1), 2]*, [(0, 0, 0), 3]*, [(0, 0, 0), 1]*, [(0, 0, 0), 2]*]

We can also use the Patch method ``apply_substitution``, which changes the patch::

    sage: sigma = WordMorphism({1:[1,3,1], 2:[1], 3:[1,1,3,2]})
    sage: E = E1Star(sigma)
    sage: p = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: p.apply_substitution(E, 5)
    sage: p
    Patch of 723 faces

The application of an ``E1Star`` substitution assigns to each new face the color of its preimage.
The ``repaint`` method allows us to repaint the faces of a patch.
A single color can also be assigned to every face, by specifying a list of a single color::

    sage: p = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: p.apply_substitution(E, 5)
    sage: p.repaint(['green'])
    sage: p.plot()                   #not tested

A list of colors allows us to color the faces sequentially::

    sage: p = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: p.apply_substitution(E)
    sage: p.repaint(['red', 'yellow', 'green', 'blue', 'black'])
    sage: p.apply_substitution(E, 3)
    sage: p.plot()                   #not tested

All the color schemes from ``matplotlib.cm.datad.keys()`` can be used::

    sage: p = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: p.repaint(cmap='summer')
    sage: p.apply_substitution(E, 3)
    sage: p.plot()                   #not tested
    sage: p.repaint(cmap='hsv')
    sage: p.apply_substitution(E, 2)
    sage: p.plot()                   #not tested

It is also possible to specify a dictionary to color the faces according to their type::

    sage: p = Patch([Face((0,0,0),t) for t in [1,2,3]])
    sage: p.apply_substitution(E, 5)
    sage: p.repaint({1:(0.7, 0.7, 0.7), 2:(0.5,0.5,0.5), 3:(0.3,0.3,0.3)})
    sage: p.plot()                   #not tested
    sage: p.repaint({1:'red', 2:'yellow', 3:'green'})
    sage: p.plot()                   #not tested

Let us look at a nice big patch in 3D::

    sage: sigma = WordMorphism({1:[1,2], 2:[3], 3:[1]})
    sage: E = E1Star(sigma)
    sage: x = [Face((0,0,0),t) for t in [1,2,3]]
    sage: P = Patch(x)
    sage: P.add(P.translate([-1,1,0]))
    sage: P.apply_substitution(E, 11)
    sage: P.plot3d()                 #not tested

Plotting with TikZ pictures is possible::

    sage: sigma = WordMorphism({1:[1,2], 2:[3], 3:[1]})
    sage: E = E1Star(sigma)
    sage: x = [Face((0,0,0),t) for t in [1,2,3]]
    sage: P = Patch(x)
    sage: P.apply_substitution(E, 3)
    sage: s = P.plot_tikz()
    sage: print s
    \begin{tikzpicture}[x={(-0.216506cm,-0.125000cm)}, y={(0.216506cm,-0.125000cm)}, z={(0.000000cm,0.250000cm)}]
    \def\loza#1#2#3#4#5#6{
      \definecolor{facecolor}{rgb}{#4,#5,#6}
      \fill[fill=facecolor, draw=black, shift={(#1, #2, #3)}]
      (0, 0, 0) -- (0, 1, 0) -- (0, 1, 1) -- (0, 0, 1) -- cycle;
    }
    \def\lozb#1#2#3#4#5#6{
      \definecolor{facecolor}{rgb}{#4,#5,#6}
      \fill[fill=facecolor, draw=black, shift={(#1, #2, #3)}]
      (0, 0, 0) -- (0, 0, 1) -- (1, 0, 1) -- (1, 0, 0) -- cycle;
    }
    \def\lozc#1#2#3#4#5#6{
      \definecolor{facecolor}{rgb}{#4,#5,#6}
      \fill[fill=facecolor, draw=black, shift={(#1, #2, #3)}]
      (0, 0, 0) -- (1, 0, 0) -- (1, 1, 0) -- (0, 1, 0) -- cycle;
    }
    \loza{0}{0}{1}{1.00}{0.00}{0.00}
    \lozc{-1}{0}{2}{1.00}{0.00}{0.00}
    \lozb{-1}{1}{1}{1.00}{0.00}{0.00}
    \loza{0}{0}{0}{1.00}{0.00}{0.00}
    \loza{1}{-1}{0}{0.00}{1.00}{0.00}
    \lozc{0}{-1}{1}{0.00}{1.00}{0.00}
    \lozb{0}{0}{0}{0.00}{1.00}{0.00}
    \loza{1}{0}{-1}{0.00}{0.00}{1.00}
    \lozc{0}{0}{0}{0.00}{0.00}{1.00}
    \end{tikzpicture}

Note that everything can be done in any dimension, but that the plotting
features only work in dimension three (that is, with three-letter alphabets)::

    sage: x = [Face((0,0,0,0),1), Face((0,0,0,0),4)]
    sage: P = Patch(x)
    sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1,4], 4:[1]})
    sage: E = E1Star(sigma)
    sage: E(P)
    Patch: [[(1, 0, 0, -1), 1]*, [(0, 1, 0, -1), 2]*, [(0, 0, 1, -1), 3]*, [(0, 0, 0, 0), 4]*, [(0, 0, 0, 0), 3]*]

::

    sage: sigma = WordMorphism({1:[1,2],2:[1,3],3:[1,4],4:[1,5],5:[1,6],6:[1,7],7:[1,8],8:[1,9],9:[1,10],10:[1,11],11:[1,12],12:[1]})
    sage: E = E1Star(sigma)
    sage: E
    E_1^*(WordMorphism: 1->12, 10->1,11, 11->1,12, 12->1, 2->13, 3->14, 4->15, 5->16, 6->17, 7->18, 8->19, 9->1,10)
    sage: P = Patch([Face((0,0,0,0,0,0,0,0,0,0,0,0),t) for t in [1,2,3]])
    sage: for x in E(P): print x
    [(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1), 1]*
    [(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1), 2]*
    [(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1), 3]*
    [(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1), 4]*
    [(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1), 5]*
    [(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1), 6]*
    [(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1), 7]*
    [(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1), 8]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1), 9]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1), 10]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1), 11]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 12]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 1]*
    [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 2]*
"""
#*****************************************************************************
#       Copyright (C) 2010 Franco Saliola <saliola@gmail.com>
#                          Vincent Delecroix <20100.delecroix@gmail.com>
#                          Timo Jolivet <timo.jolivet@gmail.com>
#                          Stepan Starosta <stepan.starosta@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.integer import Integer
from sage.misc.functional import det
from sage.structure.sage_object import SageObject
from sage.combinat.words.morphism import WordMorphism
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.plot.plot import Graphics
from sage.plot.plot import rainbow
from sage.plot.colors import Color
from sage.plot.polygon import polygon
from sage.rings.integer_ring import ZZ
from sage.misc.latex import LatexExpr
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

# matplotlib color maps, loaded on-demand
cm = None

class Face(SageObject):
    r"""
    A class to model a unit face of arbitrary dimension.

    A unit face in dimension `d` is represented by
    a `d`-dimensional vector ``v`` and a type ``t`` in `\{1, \ldots, d\}`.
    The type of the face corresponds to the canonical unit vector
    to which the face is orthogonal.
    The optional ``color`` argument is used in plotting functions.

    INPUT:

    - ``v`` - tuple of integers
    - ``t`` - integer in ``[1, ..., len(v)]``, type of the face. The face of type `i`
      is orthogonal to the canonical vector `e_i`.
    - ``color`` - color (optional, default: ``None``) color of the face,
      used for plotting only. If ``None``, its value is guessed from the
      face type.

    EXAMPLES::

        sage: from sage.combinat.e_one_star import E1Star, Face, Patch
        sage: f = Face((0,2,0), 3)
        sage: f.vector()
        (0, 2, 0)
        sage: f.type()
        3

    ::

        sage: f = Face((0,2,0), 3, color=(0.5, 0.5, 0.5))
        sage: f.color()
        RGB color (0.5, 0.5, 0.5)
    """
    def __init__(self, v, t, color=None):
        r"""
        Face constructor. See class doc for more information.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,2,0), 3)
            sage: f.vector()
            (0, 2, 0)
            sage: f.type()
            3
        """
        self._vector = (ZZ**len(v))(v)
        self._vector.set_immutable()

        if not(isinstance(t, Integer) and 1 <= t <= len(v)):
            raise ValueError, 'The type must be an integer between 1 and len(v)'
        self._type = t

        if color is None:
            if self._type == 1:
                color = Color((1,0,0))
            elif self._type == 2:
                color = Color((0,1,0))
            elif self._type == 3:
                color = Color((0,0,1))
            else:
                color = Color()
        self._color = Color(color)

    def __repr__(self):
        r"""
        String representation of a face.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,0,0,3), 3)
            sage: f
            [(0, 0, 0, 3), 3]*

        ::

            sage: f = Face((0,0,0,3), 3)
            sage: f
            [(0, 0, 0, 3), 3]*
        """
        return "[%s, %s]*"%(self.vector(), self.type())

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,0,0,3), 3)
            sage: g = Face((0,0,0,3), 3)
            sage: f == g
            True
        """
        return (isinstance(other, Face) and
                self.vector() == other.vector() and
                self.type() == other.type() )


    def __hash__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,0,0,3), 3)
            sage: g = Face((0,0,0,3), 3)
            sage: hash(f) == hash(g)
            True
        """
        return hash((self.vector(), self.type()))

    def vector(self):
        r"""
        Return the vector of the face.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,2,0), 3)
            sage: f.vector()
            (0, 2, 0)
        """
        return self._vector

    def type(self):
        r"""
        Returns the type of the face.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,2,0), 3)
            sage: f.type()
            3

        ::

            sage: f = Face((0,2,0), 3)
            sage: f.type()
            3
        """
        return self._type

    def color(self, color=None):
        r"""
        Returns or change the color of the face.

        INPUT:

        - ``color`` - string, rgb tuple, color (optional, default: ``None``)
          the new color to assign to the face. If ``None``, it returns the
          color of the face.

        OUTPUT:

            color

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,2,0), 3)
            sage: f.color()
            RGB color (0.0, 0.0, 1.0)
            sage: f.color('red')
            sage: f.color()
            RGB color (1.0, 0.0, 0.0)

        """
        if color is None:
            return self._color
        else:
            self._color = Color(color)

    def _plot(self, projmat, face_contour, opacity):
        r"""
        Returns a 2D graphic object representing the face.

        INPUT:

        - ``projmat`` - 2*3 projection matrix
        - ``face_contour`` - dict, maps the face type to vectors describing
          the contour of unit faces

        OUTPUT:

            2D graphic object

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,0,3), 3)
            sage: projmat = matrix(2, [-1.7320508075688772*0.5, 1.7320508075688772*0.5, 0, -0.5, -0.5, 1])
            sage: face_contour = {1: map(vector, [(0,0,0),(0,1,0),(0,1,1),(0,0,1)]), 2: map(vector, [(0,0,0),(0,0,1),(1,0,1),(1,0,0)]), 3: map(vector, [(0,0,0),(1,0,0),(1,1,0),(0,1,0)])}
            sage: G = f._plot(projmat, face_contour, 0.75)
        """
        v = self.vector()
        t = self.type()
        G = polygon([projmat*(u+v) for u in face_contour[t]], alpha=opacity,
               thickness=1, rgbcolor=self.color())
        return G

    def _plot3d(self, face_contour):
        r"""
        3D reprensentation of a unit face (Jmol).

        INPUT:

        - ``face_contour`` - dict, maps the face type to vectors describing
          the contour of unit faces

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,0,3), 3)
            sage: face_contour = {1: map(vector, [(0,0,0),(0,1,0),(0,1,1),(0,0,1)]), 2: map(vector, [(0,0,0),(0,0,1),(1,0,1),(1,0,0)]), 3: map(vector, [(0,0,0),(1,0,0),(1,1,0),(0,1,0)])}
            sage: G = f._plot3d(face_contour)      #not tested
        """
        v = self.vector()
        t = self.type()
        c = self.color()
        G = polygon([u+v for u in face_contour[t]], rgbcolor=c)
        return G

class Patch(SageObject):
    r"""
    A class to model a collection of faces.

    INPUT:

    - ``faces`` - finite iterable of faces
    - ``face_contour`` - dict (optional, default:``None``) maps the face
      type to vectors describing the contour of unit faces. If None,
      defaults contour are assumed for faces of type 1, 2, 3 or 1, 2, 3.
      Used in plotting methods only.

    EXAMPLES::

        sage: from sage.combinat.e_one_star import E1Star, Face, Patch
        sage: x = [Face((0,0,0),t) for t in [1,2,3]]
        sage: P = Patch(x)
        sage: P
        Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]

    ::

        sage: face_contour = {}
        sage: face_contour[1] =  map(vector, [(0,0,0),(0,1,0),(0,1,1),(0,0,1)])
        sage: face_contour[2] = map(vector, [(0,0,0),(0,0,1),(1,0,1),(1,0,0)])
        sage: face_contour[3] = map(vector, [(0,0,0),(1,0,0),(1,1,0),(0,1,0)])
        sage: Patch([Face((0,0,0),t) for t in [1,2,3]], face_contour=face_contour)
        Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]
    """
    def __init__(self, faces, face_contour=None):
        r"""
        Constructor of a patch (set of faces). See class doc for more information.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]
            sage: y = Face((0,1,0), 3)
            sage: P.add(y)
            sage: P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*, [(0, 1, 0), 3]*]
        """
        self._faces = list(faces)

        if not face_contour is None:
            self._face_contour = face_contour

        else:
            self._face_contour = {
                    1: map(vector, [(0,0,0),(0,1,0),(0,1,1),(0,0,1)]),
                    2: map(vector, [(0,0,0),(0,0,1),(1,0,1),(1,0,0)]),
                    3: map(vector, [(0,0,0),(1,0,0),(1,1,0),(0,1,0)])
            }

    def __eq__(self, other):
        r"""
        Equality test for Patch.

        INPUT:

        - ``other`` - an object

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: x = [Face((0,0,0),1), Face((0,0,0),2), Face((0,0,0),3)]
            sage: y = [Face((0,1,0),1), Face((0,0,0),3)]
            sage: P = Patch(x)
            sage: Q = Patch(y)
            sage: R = Patch(x)
            sage: P == Q
            False
            sage: P == R
            True
            sage: P == 4
            False

        ::

            sage: s = WordMorphism({1:[1,3], 2:[1,2,3], 3:[3]})
            sage: t = WordMorphism({1:[1,2,3], 2:[2,3], 3:[3]})
            sage: P = Patch([Face((0,0,0), 1), Face((0,0,0), 2), Face((0,0,0), 3)])
            sage: E1Star(s)(P) == E1Star(t)(P)
            False
            sage: E1Star(s * t)(P) == E1Star(t)(E1Star(s)(P))
            True
        """
        return (isinstance(other, Patch) and set(self) == set(other))

    def __len__(self):
        r"""
        Returns the number of faces contained in the patch.

        OUPUT:

            integer

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: len(P)       #indirect doctest
            3
        """
        return len(self._faces)

    def __repr__(self):
        r"""
        String representation of a patch.

        Displays all the faces if there less than 20,
        otherwise displays the number of faces.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]

        ::

            sage: x = [Face((0,0,a),1) for a in range(25)]
            sage: P = Patch(x)
            sage: P
            Patch of 25 faces
        """
        if len(self) <= 20:
            return "Patch: %s"%self._faces
        else:
            return "Patch of %s faces"%len(self)

    def __iter__(self):
        r"""
        Return an iterator over the faces of the patch.

        OUTPUT:

            iterator

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: list(P)
            [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*]
        """
        return iter(self._faces)

    def __getitem__(self, key):
        r"""
        INPUT:

        - ``key`` - integer between 0 and len(self) - 1

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2,3]]
            sage: P = Patch(x)
            sage: P[1]
            [(0, 0, 0), 2]*
        """
        return self._faces[key]

    def __add__(self, other):
        r"""
        Addition of patches (union). A face already in both patches will
        appear twice in the new patch. No verification of unicity is done.

        INPUT:

        - ``other`` - an object

        EXAMPLES::

            sage: from sage.combinat.e_one_star import Face, Patch, E1Star
            sage: U = Patch([Face([0,0,0], 1), Face([0,0,0], 2), Face([0,0,0], 3)])
            sage: V = U.translate([1,-1,0])
            sage: U + V
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(0, 0, 0), 3]*, [(1, -1, 0), 1]*, [(1, -1, 0), 2]*, [(1, -1, 0), 3]*]
        """
        return Patch(self._faces + other._faces)

    def add(self, faces):
        r"""
        Add a face or many faces to the patch. This changes the patch.

        .. NOTE::

            The faces to be added are assumed to be distinct from those
            already in the patch. A face already in the patch will be
            added twice. No verification of unicity is done.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: x = [Face((0,0,0),t) for t in [1,2]]
            sage: P = Patch(x)
            sage: P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*]
            sage: P.add(Face((1,2,3), 3))
            sage: P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(1, 2, 3), 3]*]
            sage: P.add([Face((1,2,3), 3), Face((2,3,3), 2)])
            sage: P
            Patch: [[(0, 0, 0), 1]*, [(0, 0, 0), 2]*, [(1, 2, 3), 3]*, [(1, 2, 3), 3]*, [(2, 3, 3), 2]*]
        """
        if isinstance(faces, Face):
            self._faces.append(faces)
        else:
            self._faces.extend(faces)

    def faces_of_type(self, t):
        r"""
        Returns a list of the faces that have type ``t``.

        INPUT:

        - ``t`` - integer or any other type

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),1), Face((1,2,0),3), Face((1,2,0),1)])
            sage: P.faces_of_type(1)
            [[(0, 0, 0), 1]*, [(1, 2, 0), 1]*]
        """
        return [face for face in self if face.type() == t]

    def translate(self, v):
        r"""
        Translates all the faces by the vector ``v``.

        This does not change self. It returns a copy.

        INPUT:

        - ``v`` - vector or tuple

        OUTPUT:

            a patch

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),1), Face((1,2,0),3), Face((1,2,0),1)])
            sage: P.translate([-1,-2,0])
            Patch: [[(-1, -2, 0), 1]*, [(0, 0, 0), 3]*, [(0, 0, 0), 1]*]
        """
        v = vector(v)
        return Patch(Face(f.vector()+v, f.type(), f.color()) for f in self)

    def repaint(self, cmap='prism'):
        r"""
        Repaints all the faces of self from the given color map.

        This changes self.

        INPUT:

        -  ``cmap`` - color map (default: ``'prism'``). It can be one of the
           following :

           - string - A coloring map. For available coloring map names type:
             ``sorted(colormaps)``
           - list - a list of colors to assign cyclically to the faces.
             A list of a single color colors all the faces with the same color.
           - dict - a dict of face types mapped to colors, to color the
             faces according to their type.
           - ``{}``, the empty dict - shorcut for
             ``{1:'red', 2:'green', 3:'blue'}``.

        EXAMPLES:

        Using a color map::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: color = (0, 0, 0)
            sage: p = Patch([Face((0,0,0),t,color) for t in [1,2,3]])
            sage: for f in p: f.color()
            RGB color (0.0, 0.0, 0.0)
            RGB color (0.0, 0.0, 0.0)
            RGB color (0.0, 0.0, 0.0)
            sage: p.repaint()
            sage: p[1].color()
            RGB color (0.198..., 0.912..., 0.0)

        Using a list of colors::

            sage: p = Patch([Face((0,0,0),t,color) for t in [1,2,3]])
            sage: p.repaint([(0.9, 0.9, 0.9), (0.65,0.65,0.65), (0.4,0.4,0.4)])
            sage: for f in p: f.color()
            RGB color (0.9000..., 0.9000..., 0.9000...)
            RGB color (0.6500..., 0.6500..., 0.6500...)
            RGB color (0.4000..., 0.4000..., 0.4000...)

        Using a dictionary to color faces according to their type::

            sage: p = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: p.repaint({1:'black', 2:'yellow', 3:'green'})
            sage: p.plot()                   #not tested
            sage: p.repaint({})
            sage: p.plot()                   #not tested
        """
        if cmap == {}:
            cmap = {1: 'red', 2:'green', 3:'blue'}

        if isinstance(cmap, dict):
            for f in self:
                f.color(cmap[f.type()])

        elif isinstance(cmap, list):
            L = len(cmap)
            for i, f in enumerate(self):
                f.color(cmap[i % L])

        elif isinstance(cmap, str):
            # matplotlib color maps
            global cm
            if not cm:
                from matplotlib import cm

            if not cmap in cm.datad.keys():
                raise RuntimeError("Color map %s not known (type sorted(colors) for valid names)" % cmap)
            cmap = cm.__dict__[cmap]
            dim = float(len(self))
            for i,f in enumerate(self):
                f.color(cmap(i/dim)[:3])

        else:
            raise TypeError, "Type of cmap (=%s) must be dict, list or str" %cmap

    def apply_substitution(self, E, iterations=1):
        r"""
        Applies the substitution ``E`` on the patch.

        The color of every new face in the image is given the same color as its preimage.

        This changes self.

        INPUT:

        - ``E`` - an instance of the ``E1Star`` class. Its domain alphabet must
          be of the same size as the dimension of the ambient space of the
          faces.
        - ``iterations`` - integer (optional, default: 1)
          number of times the substitution E is applied

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: p = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: p.apply_substitution(E)
            sage: p.plot()                   #not tested

        ::

            sage: x = (0,0,0)
            sage: p = Patch([Face(x, 1, 'red'), Face(x, 2, 'yellow'), Face(x, 3, 'green')])
            sage: p.apply_substitution(E, 4)
            sage: p
            Patch of 31 faces
            sage: p.plot()                   #not tested
        """
        if not isinstance(E, E1Star):
            raise TypeError, "E must be an instance of E1Star"

        self._faces = E(self, iterations=iterations)._faces

    def plot(self, projmat=None, opacity=0.75):
        r"""
        Returns a 2D graphic object depicting the patch.

        INPUT:

        - ``projmat`` - matrix (optional, default: ``None``) the projection
          matrix. Its number of lines must be two. Its number of columns
          must equal the dimension of the ambient space of the faces. If
          ``None``, the isometric projection is used by default.

        - ``opacity`` - float between ``0`` and ``1`` (optional, default: ``0.75``)
          opacity of the the face (edges cannot be seen with opacity equal to ``1``)

        .. WARNING::

            Only three-dimensional Faces and Patches can be plotted.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.plot()                   #not tested

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.apply_substitution(E, 5)
            sage: P.plot()

        Plot with a different projection matrix::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: M = matrix(2, 3, [1,0,-1,0.3,1,-3])
            sage: P.apply_substitution(E, 3)
            sage: P.plot(projmat=M)
        """
        if len(self[0].vector()) != 3:
            raise NotImplementedError, "Plotting only works for three-dimensional Patches and faces."

        #projmat = matrix(2, 3, [1,0,-1/sqrt(2),0,1,-1/sqrt(2)])
        if projmat == None:
            projmat = matrix(2, [-1.7320508075688772*0.5, 1.7320508075688772*0.5, 0, -0.5, -0.5, 1])

        G = Graphics()
        for face in self:
            G += face._plot(projmat, self._face_contour, opacity)
        G.set_aspect_ratio(1)
        return G

    def plot3d(self):
        r"""
        Returns a 3D graphics object depicting the patch.

        .. WARNING::

            Only three-dimensional Faces and Patches can be plotted.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.plot3d()                #not tested

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.apply_substitution(E, 5)
            sage: P.repaint()
            sage: P.plot3d()                #not tested
        """
        if len(self[0].vector()) != 3:
            raise NotImplementedError, "Plotting only works for three-dimensional Patches and faces."

        from sage.all import pi
        face_list = [face._plot3d(self._face_contour) for face in self]
        P = sum(face_list)
        return P

    def plot_tikz(self, projmat=None, print_macros=True,
            print_tikz_env=True, edgecolor='black', scale=0.25,
            drawzero=False, extra_code_before='', extra_code_after=''):
        r"""
        Returns a string containing some TikZ code to be included into
        a LaTeX document.

        .. WARNING::

            Only three-dimensional Faces and Patches can be plotted.

        INPUT:

        - ``projmat`` - matrix (optional, default: ``None``) the projection
          matrix. Its number of lines must be two. Its number of columns
          must equal the dimension of the ambient space of the faces. If
          ``None``, the isometric projection is used by default.
        - ``print_macros`` - bool (optional, default: ``True``) if ``True``,
          the three lozenge macros are printed
        - ``print_tikz_env`` - bool (optional, default: ``True``) if ``True``,
          the tikzpicture environment are printed
        - ``edgecolor`` - string (optional, default: ``'black'``) either
          ``'black'`` or ``'facecolor'`` (color of unit face edges)
        - ``scale`` - real number (optional, default: ``0.25``) scaling
          constant for the whole figure
        - ``drawzero`` - bool (optional, default: ``False``) if ``True``,
          mark the origin by a black dot
        - ``extra_code_before`` - string (optional, default: ``''``) extra code to
          include in the tikz picture
        - ``extra_code_after`` - string (optional, default: ``''``) extra code to
          include in the tikz picture

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: s = P.plot_tikz()
            sage: print s
            \begin{tikzpicture}[x={(-0.216506cm,-0.125000cm)}, y={(0.216506cm,-0.125000cm)}, z={(0.000000cm,0.250000cm)}]
            \def\loza#1#2#3#4#5#6{
              \definecolor{facecolor}{rgb}{#4,#5,#6}
              \fill[fill=facecolor, draw=black, shift={(#1, #2, #3)}]
              (0, 0, 0) -- (0, 1, 0) -- (0, 1, 1) -- (0, 0, 1) -- cycle;
            }
            \def\lozb#1#2#3#4#5#6{
              \definecolor{facecolor}{rgb}{#4,#5,#6}
              \fill[fill=facecolor, draw=black, shift={(#1, #2, #3)}]
              (0, 0, 0) -- (0, 0, 1) -- (1, 0, 1) -- (1, 0, 0) -- cycle;
            }
            \def\lozc#1#2#3#4#5#6{
              \definecolor{facecolor}{rgb}{#4,#5,#6}
              \fill[fill=facecolor, draw=black, shift={(#1, #2, #3)}]
              (0, 0, 0) -- (1, 0, 0) -- (1, 1, 0) -- (0, 1, 0) -- cycle;
            }
            \loza{0}{0}{0}{1.00}{0.00}{0.00}
            \lozb{0}{0}{0}{0.00}{1.00}{0.00}
            \lozc{0}{0}{0}{0.00}{0.00}{1.00}
            \end{tikzpicture}

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.apply_substitution(E, 2)
            sage: s = P.plot_tikz()

        Plot using shades of gray (useful for article figures)::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: P.repaint([(0.9, 0.9, 0.9), (0.65,0.65,0.65), (0.4,0.4,0.4)])
            sage: P.apply_substitution(E, 4)
            sage: P
            Patch of 31 faces
            sage: s = P.plot_tikz()

        Plotting with various options::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: M = matrix(2, 3, map(float, [1,0,-0.7071,0,1,-0.7071]))
            sage: P.apply_substitution(E, 3)
            sage: s = P.plot_tikz(projmat=M, edgecolor='facecolor', scale=0.6, drawzero=True)

        Adding X, Y, Z axes using the extra code feature::

            sage: length = 1.5
            sage: space = 0.3
            sage: axes = ''
            sage: axes += "\\draw[->, thick, black] (0,0,0) -- (%s, 0, 0);\n" % length
            sage: axes += "\\draw[->, thick, black] (0,0,0) -- (0, %s, 0);\n" % length
            sage: axes += "\\node at (%s,0,0) {$x$};\n" % (length + space)
            sage: axes += "\\node at (0,%s,0) {$y$};\n" % (length + space)
            sage: axes += "\\node at (0,0,%s) {$z$};\n" % (length + space)
            sage: axes += "\\draw[->, thick, black] (0,0,0) -- (0, 0, %s);\n" % length
            sage: cube = Patch([Face((0,0,0),1), Face((0,0,0),2), Face((0,0,0),3)])
            sage: options = dict(print_macros=False,scale=0.5,drawzero=True,extra_code_before=axes)
            sage: r = cube.plot_tikz(**options)
            sage: len(r)
            611
            sage: print r
            \begin{tikzpicture}[x={(-0.433013cm,-0.250000cm)}, y={(0.433013cm,-0.250000cm)}, z={(0.000000cm,0.500000cm)}]
            \draw[->, thick, black] (0,0,0) -- (1.50000000000000, 0, 0);
            \draw[->, thick, black] (0,0,0) -- (0, 1.50000000000000, 0);
            \node at (1.80000000000000,0,0) {$x$};
            \node at (0,1.80000000000000,0) {$y$};
            \node at (0,0,1.80000000000000) {$z$};
            \draw[->, thick, black] (0,0,0) -- (0, 0, 1.50000000000000);
            \loza{0}{0}{0}{1.00}{0.00}{0.00}
            \lozb{0}{0}{0}{0.00}{1.00}{0.00}
            \lozc{0}{0}{0}{0.00}{0.00}{1.00}
            \node[circle,fill=black,draw=black,minimum size=1.5mm,inner sep=0pt] at (0,0,0) {};
            \end{tikzpicture}
            <BLANKLINE>
        """
        if len(self[0].vector()) != 3:
            raise NotImplementedError, "Plotting only works for three-dimensional Patches and faces."

        if projmat == None:
            projmat = matrix(2, [-1.7320508075688772*0.5, 1.7320508075688772*0.5, 0, -0.5, -0.5, 1])*scale

        e1 = projmat*vector([1,0,0])
        e2 = projmat*vector([0,1,0])
        e3 = projmat*vector([0,0,1])

        # string s contains the TiKZ code of the patch
        s = ''

        if print_tikz_env:
            s += '\\begin{tikzpicture}[x={(%fcm,%fcm)}, y={(%fcm,%fcm)}, z={(%fcm,%fcm)}]\n'%(e1[0], e1[1], e2[0], e2[1], e3[0], e3[1])

        if print_macros:
            for type in self._face_contour:
                # chr(97) is 'a'
                s += '\\def\\loz%s#1#2#3#4#5#6{\n' % chr(int(type) + 96)
                s += '  \\definecolor{facecolor}{rgb}{#4,#5,#6}\n'
                s += '  \\fill[fill=facecolor, draw=%s, shift={(#1, #2, #3)}]\n'%edgecolor
                s += ' -- '.join(map(str, self._face_contour[type])) + ' -- cycle;\n'
                s += '}\n'

        s += extra_code_before

        for face in self:
            t = chr(int(face.type()) + 96)
            x, y, z = face.vector()
            r, g, b = face.color()
            s += '\\loz%s{%d}{%d}{%d}{%.2f}{%.2f}{%.2f}\n' %(t, x, y, z, r, g, b)

        s += extra_code_after

        if drawzero:
            s += '\\node[circle,fill=black,draw=black,minimum size=1.5mm,inner sep=0pt] at (0,0,0) {};\n'

        if print_tikz_env:
            s += '\\end{tikzpicture}\n'

        return LatexExpr(s)

    _latex_ = plot_tikz

class E1Star(SageObject):
    r"""
    A class to model the `E_1^*(\sigma)` map associated with
    a unimodular substitution `\sigma`.

    INPUT:

    - ``sigma`` - unimodular ``WordMorphism``, i.e. such that its incidence
      matrix has determinant `\pm 1`.

    .. NOTE::

        The alphabet of the domain and the codomain of `\sigma` must be
        equal, and they must be of the form ``[1, ..., d]``, where ``d``
        is a positive integer corresponding to the length of the vectors
        of the faces on which `E_1^*(\sigma)` will act.

    EXAMPLES::

        sage: from sage.combinat.e_one_star import E1Star, Face, Patch
        sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
        sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
        sage: E = E1Star(sigma)
        sage: E(P)
        Patch: [[(1, 0, -1), 1]*, [(0, 1, -1), 2]*, [(0, 0, 0), 3]*, [(0, 0, 0), 1]*, [(0, 0, 0), 2]*]

    ::

        sage: x = [Face((0,0,0,0),1), Face((0,0,0,0),4)]
        sage: P = Patch(x)
        sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1,4], 4:[1]})
        sage: E = E1Star(sigma)
        sage: E(P)
        Patch: [[(1, 0, 0, -1), 1]*, [(0, 1, 0, -1), 2]*, [(0, 0, 1, -1), 3]*, [(0, 0, 0, 0), 4]*, [(0, 0, 0, 0), 3]*]
    """
    def __init__(self, sigma):
        r"""
        E1Star constructor. See class doc for more information.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E
            E_1^*(WordMorphism: 1->12, 2->13, 3->1)
        """
        if not isinstance(sigma, WordMorphism):
            raise TypeError, "sigma (=%s) must be an instance of WordMorphism"%sigma

        if sigma.domain().alphabet() != sigma.codomain().alphabet():
            raise ValueError, "The domain and codomain of (%s) must be the same."%sigma

        if abs(det(matrix(sigma))) != 1:
            raise ValueError, "The substitution (%s) must be unimodular."%sigma

        first_letter = sigma.codomain().alphabet()[0]
        if not isinstance(first_letter, Integer) or (first_letter < 1):
            raise ValueError, "The substitution (%s) must be defined on positive integers."%sigma

        self._sigma = WordMorphism(sigma)
        self._d = self._sigma.domain().size_of_alphabet()

    @cached_method
    def matrix(self):
        r"""
        Returns the matrix associated with self.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E.matrix()
            [1 1 1]
            [1 0 0]
            [0 1 0]
        """
        return self._sigma.incidence_matrix()

    @cached_method
    def inverse_matrix(self):
        r"""
        Returns the inverse of the matrix associated with self.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E.inverse_matrix()
            [ 0  1  0]
            [ 0  0  1]
            [ 1 -1 -1]

        """
        return self.matrix().inverse()

    @lazy_attribute
    def _base_iter(self):
        r"""
        Returns a base for the iteration of the application of self on set
        of faces. (Exploits the linearity of `E_1^*(\sigma)` to optimize computation.)

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E._base_iter.keys()
            [1, 2, 3]
            sage: E._base_iter[1]
            [((1, 0, -1), 1), ((0, 1, -1), 2), ((0, 0, 0), 3)]
            sage: E._base_iter[2]
            [((0, 0, 0), 1)]
            sage: E._base_iter[3]
            [((0, 0, 0), 2)]

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: sorted(E._base_iter.keys())
            [1, 2, 3]
            sage: E._base_iter[1]
            [((1, 0, -1), 1), ((0, 1, -1), 2), ((0, 0, 0), 3)]
            sage: E._base_iter[2]
            [((0, 0, 0), 1)]
            sage: E._base_iter[3]
            [((0, 0, 0), 2)]

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: sorted(E._base_iter.keys())
            [1, 2, 3]
            sage: E._base_iter[1]
            [((1, 0, -1), 1), ((0, 1, -1), 2), ((0, 0, 0), 3)]
            sage: E._base_iter[1]
            [((1, 0, -1), 1), ((0, 1, -1), 2), ((0, 0, 0), 3)]
            sage: E._base_iter[2]
            [((0, 0, 0), 1)]
            sage: E._base_iter[3]
            [((0, 0, 0), 2)]

        ::

            sage: sigma = WordMorphism({1:[1,3,1], 2:[1], 3:[1,1,3,2]})
            sage: E = E1Star(sigma)
            sage: sorted(E._base_iter.keys())
            [1, 2, 3]
            sage: E._base_iter[1]
            [((1, -1, 0), 1), ((0, 0, 0), 1), ((0, 0, 0), 2), ((0, -1, 1), 3), ((0, -2, 1), 3)]
            sage: E._base_iter[2]
            [((0, 0, 0), 3)]
            sage: E._base_iter[3]
            [((0, 1, 0), 1), ((-1, 0, 1), 3)]

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: sorted(E._base_iter.keys())
            [1, 2, 3]
            sage: E._base_iter[1]
            [((1, 0, -1), 1), ((0, 0, 0), 3)]
            sage: E._base_iter[2]
            [((0, 0, 0), 1)]
            sage: E._base_iter[3]
            [((0, 0, 0), 2)]

        ::

            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1,4], 4:[1]})
            sage: E = E1Star(sigma)
            sage: sorted(E._base_iter.keys())
            [1, 2, 3, 4]
            sage: E._base_iter[1]
            [((1, 0, 0, -1), 1), ((0, 1, 0, -1), 2), ((0, 0, 1, -1), 3), ((0, 0, 0, 0), 4)]
            sage: E._base_iter[2]
            [((0, 0, 0, 0), 1)]
            sage: E._base_iter[3]
            [((0, 0, 0, 0), 2)]
            sage: E._base_iter[4]
            [((0, 0, 0, 0), 3)]
        """
        alphabet = self._sigma.domain().alphabet()
        X = {}
        for k in alphabet:
            subst_im = self._sigma.image(k)
            for n, letter in enumerate(subst_im):
                suffix = subst_im[n+1:]
                if not letter in X:
                    X[letter] = []
                v = self.inverse_matrix()*vector(suffix.parikh_vector())
                X[letter].append((v, k))
        return X

    def __call__(self, patch, iterations=1):
        r"""
        Apply a generalized substitution to a Patch; this returns a new
        object.

        The color of every new face in the image is given the same color as its preimage.

        INPUT:

        - ``patch`` - a patch
        - ``iterations`` - integer (optional, default: 1) number of iterations

        OUTPUT:

            a patch

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: P = Patch([Face((0,0,0),t) for t in [1,2,3]])
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: R = E(P)
            sage: len(R)
            5

        ::

            sage: x = (0,0,0)
            sage: p = Patch([Face(x, 1, 'red'), Face(x, 2, 'yellow'), Face(x, 3, 'green')])
            sage: p = E(p, iterations=4)
            sage: p
            Patch of 31 faces
        """
        old_faces = patch
        for i in xrange(iterations):
            new_faces = []
            for f in old_faces:
                new_faces.extend(self._call_on_face(f, color=f.color()))
            old_faces = new_faces

        return Patch(new_faces)

    def __mul__(self, other):
        r"""
        Return the product of self and other.

        The product satisfies the following rule :
        `E_1^*(\sigma\circ\sigma') = E_1^*(\sigma')` \circ  E_1^*(\sigma)`

        INPUT:

        - ``other`` - an instance of E1Star

        OUTPUT:

            an instance of E1Star

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: s = WordMorphism({1:[2],2:[3],3:[1,2]})
            sage: t = WordMorphism({1:[1,3,1],2:[1],3:[1,1,3,2]})
            sage: E1Star(s) * E1Star(t)
            E_1^*(WordMorphism: 1->1, 2->1132, 3->1311)
            sage: E1Star(t * s)
            E_1^*(WordMorphism: 1->1, 2->1132, 3->1311)
        """
        if not isinstance(other, E1Star):
            raise TypeError, "other (=%s) must be an instance of E1Star" % other
        return E1Star(other.sigma() * self.sigma())

    def _call_on_face(self, face, color=None):
        r"""
        Returns an iterator of faces obtained by applying self on the face.

        INPUT:

        - ``face`` - a face
        - ``color`` - string, RGB tuple or color, (optional, default: None)
          RGB color

        OUTPUT:

            iterator of faces

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: f = Face((0,2,0), 1)
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: list(E._call_on_face(f))
            [[(3, 0, -3), 1]*, [(2, 1, -3), 2]*, [(2, 0, -2), 3]*]
        """
        if len(face.vector()) != self._d:
            raise ValueError, "The dimension of the faces must be equal to the size of the alphabet of the substitution."
        x_new = self.inverse_matrix() * face.vector()
        t = face.type()
        return (Face(x_new + v, k, color=color) for v, k in self._base_iter[t])

    def __repr__(self):
        r"""
        String representation of a patch.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: E
            E_1^*(WordMorphism: 1->12, 2->13, 3->1)
        """
        return "E_1^*(%s)" % str(self._sigma)

    def sigma(self):
        r"""
        Returns the ``WordMorphism`` associated with self.

        EXAMPLES::

            sage: from sage.combinat.e_one_star import E1Star, Face, Patch
            sage: sigma = WordMorphism({1:[1,2], 2:[1,3], 3:[1]})
            sage: E = E1Star(sigma)
            sage: print E.sigma()
            WordMorphism: 1->12, 2->13, 3->1
        """
        return self._sigma


