"""
Polytopes

This module provides access to polymake, which 'has been developed
since 1997 in the Discrete Geometry group at the Institute of
Mathematics of Technische Universitat Berlin. Since 2004 the
development is shared with Fachbereich Mathematik, Technische
Universitat Darmstadt.  The system offers access to a wide variety of
algorithms and packages within a common framework. polymake is
flexible and continuously expanding. The software supplies C++ and
perl interfaces which make it highly adaptable to individual needs.'

AUTHOR:
    -- Ewgenij Gawrilow, Michael Joswig: main authors of polymake
    -- William Stein: SAGE interface
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.misc.all import SAGE_TMP, tmp_filename
from sage.rings.all import Integer, QQ
from sage.structure.all import Sequence
from sage.modules.all import VectorSpace
from sage.structure.sage_object import SageObject

import os

path = '%s/local/polymake/bin/'%os.environ['SAGE_ROOT']
polymake_command = path + 'polymake'

if os.path.exists(path):
    os.environ['PATH'] = '%s:'%path + os.environ['PATH']

tmp_file = '%s/tmp.poly'%SAGE_TMP

class Polytope(SageObject):
    def __init__(self, datafile, desc):
        self.__data = datafile
        self.__desc = desc

    def _repr_(self):
        s = self.__desc
        # vertices, facets, points, inequalities
        try:
            s += '\nVertices:\n%s'%self.__vertices
        except AttributeError:
            pass
        try:
            s += '\nFacets:\n%s'%self.__facets
        except AttributeError:
            pass
        return s



    def __add__(self, other):
        """

        """
        if not isinstance(other, Polytope):
            raise TypeError, "other (=%s) must be a polytope"%other
        output_file = tmp_filename()
        infile1 = tmp_filename()
        open(infile1,'w').write(self.__data)
        infile2 = tmp_filename()
        open(infile2,'w').write(other.__data)
        cmd = "minkowski_sum %s 1 %s 1 %s"%(output_file, infile1,
                                            infile2)
        stdin, stdout, stderr = os.popen3(cmd)
        stdin.close()
        err = stderr.read()
        if len(err) > 0:
            raise RuntimeError, err
        print stdout.read(), err
        S = polymake.from_data(open(output_file).read())
        os.unlink(infile1)
        os.unlink(infile2)
        os.unlink(output_file)
        return S


    def data(self):
        return self.__data

    def write(self, filename):
        open(filename,'w').write(self.__data)

    def cmd(self, cmd):
        cmd = cmd.upper()
        # First check if the value of the command
        # is already known.
        D = self.data()
        cmd2='\n%s\n'%cmd
        i = D.find(cmd2)
        if i != -1:
            j = D[i:].find('\n\n')
            if j == -1:
                j = len(D)
            else:
                j += i
            return D[i+len(cmd2)-1:j]

        F = tmp_file
        open(F,'w').write(self.__data)
        c = '%s %s %s'%(polymake_command, F, cmd)
        stdin, stdout, stderr = os.popen3(c)
        stdin.close()
        err = stderr.read()
        if len(err) > 0:
            raise RuntimeError, err
        ans = stdout.read()
        if len(ans) == 0:
            raise ValueError, "%s\nError executing polymake command %s"%(
                err,cmd)
        self.__data = open(F).read()
        return ans


    def facets(self):
        """
        EXAMPLES:
            sage: P = polymake.convex_hull([[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])   # optional: needs polymake
            sage: P.facets()                            # optional
            [(0, 0, 0, 1), (0, 1, 0, 0), (0, 0, 1, 0), (1, 0, 0, -1), (1, 0, -1, 0), (1, -1, 0, 0)]
        """
        try:
            return self.__facets
        except AttributeError:
            pass
        s = self.cmd('FACETS')
        s = s.rstrip().split('\n')[1:]
        if len(s) == 0:
            ans = Sequence([], immutable=True)
        else:
            n = len(s[0].split())
            V = VectorSpace(QQ, n)
            ans = Sequence((V(x.split()) for x in s), immutable=True)
        self.__facets = ans
        return ans

    def vertices(self):
        """
        EXAMPLES:
            sage: P = polymake.convex_hull([[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])     # optional: needs polymake
            sage: P.vertices()                            # optional
            [(1, 0, 0, 0), (1, 0, 0, 1), (1, 0, 1, 0), (1, 0, 1, 1), (1, 1, 0, 0), (1, 1, 0, 1), (1, 1, 1, 0), (1, 1, 1, 1)]
        """
        try:
            return self.__vertices
        except AttributeError:
            pass
        s = self.cmd('VERTICES')
        s = s.rstrip().split('\n')[1:]
        if len(s) == 0:
            ans = Sequence([], immutable=True)
        else:
            n = len(s[0].split())
            V = VectorSpace(QQ, n)
            ans = Sequence((V(x.split()) for x in s), immutable=True)
        self.__vertices = ans
        return ans

    def visual(self):
        try:
            self.cmd('visual')
        except ValueError:
            pass

    def graph(self):
        try:
            return self.__graph
        except AttributeError:
            pass
        g = self.cmd('GRAPH')
        return g

    def is_simple(self):
        r"""
        Return True if this polytope is simple.

        A polytope is \emph{simple} if the degree of each
        vertex equals the dimension of the polytope.

        EXAMPLES:
            sage: P = polymake.convex_hull([[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])        # optional: needs polymake
            sage: P.is_simple()                              # optional
            True

        AUTHORS:
            -- Edwin O'Shea (2006-05-02): Definition of simple.
        """
        try:
            return self.__is_simple
        except AttributeError:
            pass
        s = self.cmd('SIMPLE')
        i = s.find('\n')
        self.__is_simple = bool(int(s[i:]))
        return self.__is_simple



    def n_facets(self):
        """
        EXAMPLES:
            sage: P = polymake.convex_hull([[1,0,0,0], [1,0,0,1], [1,0,1,0], [1,0,1,1],  [1,1,0,0], [1,1,0,1], [1,1,1,0], [1,1,1,1]])     # optional: needs polymake
            sage: P.n_facets()                            # optional
            6
        """
        try:
            return self.__n_facets
        except AttributeError:
            pass
        s = self.cmd('N_FACETS')
        i = s.find('\n')
        self.__n_facets = Integer(s[i:])
        return self.__n_facets

class Polymake:
    def __repr__(self):
        return "Object that makes polytopes."

    def __make(self, cmd, name):
        os.system(cmd)
        try:
            d = open(tmp_file).read()
        except IOError:
            raise RuntimeError, "You may need to install the polymake package"
        return Polytope(d, name)

    def reconfigure(self):
        """
        Reconfigure polymake.

        Remember to run polymake.reconfigure() as soon as you have
        changed the customization file and/or installed missing
        software!
        """
        os.system("polymake --reconfigure")

    def associahedron(self, dimension):
        return self.__make('associahedron %s %s'%(tmp_file, dimension),
                           '%s-dimensional associahedron'%dimension)

    def birkhoff(self, n):
        return self.__make('birkhoff %s %s'%(tmp_file, n),
                           'Birkhoff %s'%n)


    def cell24(self):
        """
        EXAMPLES:
            sage: polymake.cell24()            # optional: needs polymake
            The 24-cell
        """
        return self.__make('24-cell %s'%tmp_file,
                           'The 24-cell')

    def convex_hull(self, points=[]):
        """

        EXAMPLES:
            sage: R.<x,y,z> = PolynomialRing(QQ,3)
            sage: f = x^3 + y^3 + z^3 + x*y*z
            sage: e = f.exponents()
            sage: a = [[1] + list(v) for v in e]
            sage: a
            [[1, 3, 0, 0], [1, 0, 3, 0], [1, 1, 1, 1], [1, 0, 0, 3]]
            sage: n = polymake.convex_hull(a)       # optional: needs polymake
            sage: n                                 # optional
            Convex hull of points [[1, 0, 0, 3], [1, 0, 3, 0], [1, 1, 1, 1], [1, 3, 0, 0]]
            sage: n.facets()                        # optional
            [(0, 1, 0, 0), (3, -1, -1, 0), (0, 0, 1, 0)]
            sage: n.is_simple()                     # optional
            True
            sage: n.graph()                         # optional
            'GRAPH\n{1 2}\n{0 2}\n{0 1}\n\n'
        """
        f = 'POINTS\n'
        points.sort()
        for p in points:
            f += ' '.join(str(x) for x in p) + '\n'
        f += '\n'
        return Polytope(f, 'Convex hull of points %s'%points)

    def cube(self, dimension, scale=0):
        return self.__make('cube %s %s %s'%(tmp_file, dimension, scale),
                           'Cube of dimension %s (scale %s)'%(dimension, scale))

    def from_data(self, data):
        return Polytope(data, 'A Polytope')

    def rand01(self, d, n, seed=None):
        cmd = 'rand01 %s %s %s'%(tmp_file, d, n)
        if not seed is None:
            cmd += ' -seed %s'%seed
        return self.__make(cmd,
              '%s-dimensional 0/1-polytope with %s random vertices (uniform distribution)'%(d, n))



polymake = Polymake()
