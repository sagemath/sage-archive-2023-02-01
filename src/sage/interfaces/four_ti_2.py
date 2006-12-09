"""
Interface to 4ti2.

4ti2 is a software package for algebraic, geometric and combinatorial
problems on linear spaces.



"""

import os

import sage.misc.misc
from sage.rings.all import ZZ

seq = 0
tmp = None

def tmp_dir():
    global tmp
    if tmp is None:
        tmp = sage.misc.misc.tmp_dir('4ti2')
    return tmp

bin='%s/local/lib/4ti2/'%os.environ['SAGE_ROOT']

## class FourTi2:
##     r"""
##     A 4ti2 Interface Object.

##     4ti2 is a software package for algebraic, geometric and
##     combinatorial problems on linear spaces.

##     If you use any of these function in your work, cite SAGE and the following
##     reference:
##     \begin{verbatim}
##         @Misc{4ti2,
##         author =   {4ti2 team},
##         title =    {4ti2 -- A software package for algebraic, geometric and
##                    combinatorial problems on linear spaces},
##         howpublished = {Available at www.4ti2.de}
##         }
##     \end{verbatim}
##     """
##     def __init__(self, mat='', lat=''):
##         global seq
##         self._seq = seq
##         seq = seq + 1
##         self._dir = tmp_dir()
##         self._mat = mat
##         self._lat = lat
##         self._dat = {}
##         open('%s/%s'%(self._dir, self._seq), 'w').write(self._mat)

##     def __repr__(self):
##         s = "4ti2 Linear Combinatorial Object Defined By matrix:\n%s"%self._mat
##         if self._lat:
##             s += '\nand lattice:\n%s'%self._lat
##         return s

##     def __call__(self, command, extension, options='', verbose=False):
##         try:
##             return self._dat[command]
##         except KeyError:
##             pass
##         ans = call_4ti2(self._dir, self._seq, command, extension, options, verbose)
##         self._dat[command] = ans
##         return ans

##     ########################################################################
##     # Commands:
##     ########################################################################
##     def circuits(self, options='', verbose=False):
##         """
##         """
##         return self('circuits', 'cir', options='', verbose=False)

##     def genmodel(self):
##         """
##         """
##         return self('genmodel', 'x', options='', verbose=False)

##     def gensymm(self):
##         """
##         """
##         return self('gensymm', 'x', options='', verbose=False)

##     def graver(self):
##         """
##         """
##         return self('graver', 'gra')

##     def groebner(self):
##         """
##         """
##         return self('groebner', 'gro')

##     def hilbert(self):
##         """
##         """
##         return self('hilbert', 'hil')

##     def markov(self):
##         """
##         """
##         return self('markov', 'mar')

##     def minimize(self):
##         """
##         """
##         return self('minimize', 'x')

##     def normalform(self):
##         """
##         """
##         self.groebner()  #??
##         return self('normalform', 'x')

##     def output(self):
##         """
##         """
##         return self('output', 'x')

##     def qsolve(self):
##         """
##         """
##         return self('qsolve', 'qhom')

##     def rays(self):
##         """
##         """
##         return self('rays', 'ray')

##     def walk(self):
##         """
##         """
##         return self('walk', 'x')

##     def zbasis(self):
##         """
##         """
##         return self('zbasis', 'lat')

## four_ti_2 = FourTi2


class Cone:
    def integer_points(self):
        pass

class LatticePointsInCone:
    def __init__(self, cone, lattice=None):
        self._cone = cone
        self._lattice = lattice

class Lattice:
    pass

class Lattice_as_kernel(Lattice):
    # every single function should work on this.
    def __init__(self, matrix):
        self.__matrix = matrix

    def intersection(self, cone):
        return LatticePointsInCone(self, cone)

    def groebner(self):
        pass

    def markov(self):
        pass

    def walk(self):
        pass

    def normal_form(self):
        pass

    def zbasis(self):  # Z-basis  -- lattice generators
        pass

class Lattice_with_gens(Lattice):
    def __init__(self, matrix):
        self.__matrix = matrix

def matrix_to_4ti2(x):
    e = str(x).replace(',',' ').replace('[','').replace(']','')
    return '%s %s %s'%(x.nrows(), x.ncols(), e)

def call_4ti2(dir, project, command, options='', verbose=True):

    s = 'cd "%s"; "%s"/%s %s %s >out 2>err'%(dir, bin, command,  options, project)
    if verbose:
        print s
    if os.system(s):
        raise RuntimeError, "Error running 4ti2 command: %s"%s

    out = open('%s/out'%dir).read().strip()
    err = open('%s/err'%dir).read().strip()

    if verbose:
        print out
        print err

    return (out, err)

def tex_rel(x):
    if x == '<=':
        return '\\leq{}'
    elif x == '>=':
        return '\\geq{}'
    elif x == '==':
        return '='
    raise NotImplementedError

class ExtremalRays(list):
    def __init__(self, s, linear_system):
        """
        INPUT:
            s -- string of the form n i
        """
        self._s = s
        self._linear_system = linear_system
        v = s.split('\n')
        w = v[0].strip().split()
        self._number_of_vectors = int(w[0])
        self._degree_of_vectors = int(w[1])
        # ASSUMPTION: carriage return between each vector
        self._rows = v[1:]
        self._M = ZZ ** self._degree_of_vectors

    def __repr__(self):
        return "List of %s Extremal Rays"%self._number_of_vectors

    def linear_system(self):
        return self._linear_system

    def __len__(self):
        return self._number_of_vectors

    def __getitem__(self, i):
        return self._M(self._rows[i].strip().split())


class LinearSystem:
    r"""
    A linear system of relations $A x \eps b$, where $\eps$ is a list
    of relations (<=, =, >=), $A$ is a matrix of integers, and $b$ is a
    column vector.  Also the possible orthants of $x$ can be constrained.

    LinearSystem(matrix, rel, rhs, sign)

    INPUT:
         matrix -- a matrix with integer entries

         rel -- a list of signs (as strings, e.g., '<=', '==', '>='), which
                 has length the number of rows of the matrix.

         rhs --  (default: 0 vector) a list (or vector) of integers of
                 length the number of rows.

         sign -- (default: 0 vector) a list (or vector) of sign
                 constraints, where each entry is

                       -1 -- nonpositive: this particular entry
                        where you put -1 is nonpositive

                        0 -- unconstrained: this particular entry
                        where you put 0 is not constrained

                        1 -- nonnegative: this particular entry
                        where you put 1 is nonnegative

                        2 -- nonpositive union nonnegative: the
                        entry where you put 0

                 If you put a 2 in any position, what you can
                 think of happening behind the scenes is that the
                 computation is run with both a -1 and a 1 in that
                 position, and the output is the union of the
                 resuls.  E.g., if you specify exactly 3 '2's,
                 then you have a union of 8 cases.  (Internally
                 the comptutation is not actually done by
                 splitting up like this.)

    EXAMPLES:
        sage: from sage.interfaces.fortytwo import LinearSystem
        sage: m = matrix(ZZ, 2, 4, [1,1,1,1, 1,2,3,4])
        sage: rel = ['<=', '<=']
        sage: rhs = 0
        sage: signs = [1,2,2,0]
        sage: L = LinearSystem(m, rel, rhs, signs)
        sage: L
        Linear System defined by
        Matrix:
        [1 1 1 1]
        [1 2 3 4]
        Relations: ['<=', '<=']
        Right Hand Side: [0, 0]
        Signs: [1, 2, 2, 0]
    """
    def __init__(self, matrix, rel, rhs=0, signs=0):
        matrix.set_immutable()
        self._matrix = matrix
        if len(rel) != matrix.nrows():
            raise ValueError, "number of rel must equal numbers rows of matrix"

        self._rel = rel

        if not isinstance(rel, list):
            raise TypeError, 'rel must be a list'
        for x in rel:
            if not x in ['<=', '==', '>=']:
                raise ValueError, "each rel must be <=, ==, or >="

        if signs == 0:
            signs = [0]*matrix.ncols()
        self._signs = signs

        if not isinstance(signs, list):
            raise TypeError, 'signs must be a list'
        for x in signs:
            if not x in [-1,0,1,2]:
                raise ValueError, "each sign must be -1,0,1,2"
        if len(signs) != matrix.ncols():
            raise ValueError, "number of signs must equal number of columns of matrix"

        if rhs == 0:
            rhs = [0]*matrix.nrows()
        if not isinstance(rhs, list):
            raise TypeError, 'rhs must be a list'

        # check that rhs are integers later.
        if len(rhs) != matrix.nrows():
            raise ValueError, "length of rhs must equal number of rows of matrix"

        self._rhs = rhs

        self._dat = {}

        self._dir = tmp_dir()
        global seq
        seq = seq + 1
        self._project = seq

        self.export_to_4ti2(self._dir, self._project)

    def __repr__(self):
        s = "Linear System defined by\nMatrix:\n%s\nRelations: %s\nRight Hand Side: %s\nSigns: %s"%(
            self._matrix, self._rel, self._rhs, self._signs)
        return s

    def _latex_(self):
        rel = '\\begin{array}{c} %s \\end{array}'%('\\\\\n'.join([tex_rel(x) for x in self._rel]))
        rhs = self._matrix.new_matrix(len(self._rhs), 1, self._rhs)
        signs = self._signs
        s = '\\begin{array}{cccc} %s & \\mathbf{x} & %s & %s\\\\ \n %s & & & & \\end{array}'%(
            self._matrix._latex_(),rel, rhs._latex_(), signs)
        return s


    def call_4ti2(self, command, options='', verbose=False):
        return call_4ti2(self._dir, self._project, command, options, verbose)

    def _read_file(self, extension):
        try:
            return open('%s/%s.%s'%(self._dir, self._project, extension)).read().strip()
        except IOError:
            return ''

    def matrix(self):
        return self._matrix

    def export_to_4ti2(self, dir, project):
        if not os.path.exists(dir):
            os.makedirs(dir)

        # put matrix in foo
        mat = open('%s/%s'%(dir, project), 'w')
        mat.write(matrix_to_4ti2(self._matrix))

        # put rel in foo.rel
        rel = open('%s/%s.rel'%(dir, project), 'w')
        rel.write('1 %s %s'%(len(self._rel), ' '.join([str(x[0]) for x in self._rel])))

        # put sign in foo.sign
        signs = open('%s/%s.sign'%(dir, project), 'w')
        signs.write('1 %s %s'%(len(self._signs), ' '.join([str(s) for s in self._signs])))

        # put rhs in foo.rhs
        rhs = open('%s/%s.rhs'%(dir, project), 'w')
        rhs.write('1 %s %s'%(len(self._rhs), ' '.join([str(s) for s in self._rhs])))

    def qsolve(self, options=''):
        """
        This function returns the extreme rays in the cone determined
        by this system of linear constraints.

        Every cone can be written as a sum of a pointed cone and a
        vector space.

        OUTPUT:
            -- extreme rays of the pointed cone:
                    these are specified by giving
                       number of vectors        number of variables
                       each vector v; the ray is {a*v for a nonnegative reals}.
            -- generators (defined over ZZ) for the RR-vector space.

        EXAMPLES:
            sage: from sage.interfaces.fortytwo import LinearSystem
            sage: m = matrix(ZZ, 2, 4, [1,1,1,1, 1,2,3,4])
            sage: rel = ['<=', '<=']
            sage: rhs = 0
            sage: signs = [1,2,2,0]
            sage: L = LinearSystem(m, rel, rhs, signs)
            sage: a = L.qsolve();
            sage: print a[0]
            10 4
             0 -2  0  1
             0  0 -4  3
             0  0  0 -1
             0  0  1 -1
             0  1  0 -1
             1  0 -3  2
             1  0  0 -1
             2 -3  0  1
             0  1 -2  1
             0 -1  2 -1
            sage: print a[1]
            0 6
        """
        for x in self._rhs:
            if x != 0:
                raise NotImplementedError, "qsolve is currently only implemented for homogeneous systems (i.e., with rhs=0)"
        out, err = self.call_4ti2('qsolve', options=options)
        qhom = ExtremalRays(self._read_file('qhom'), self)
        qfree = self._read_file('qfree')
        return (qhom, qfree)

    def circuits(self):
        r"""
        The support minimal elements of $Ax \eps 0$, where $\eps$ are
        the list of relations.

        \note{Currently rhs must be 0.}

        ALGORITHM: Calls qsolve with signs all 2.
        """
        if rhs != 0:
            raise NotImplementedError
        pass

    def rays(self):
        """
        ALGORITHM: Calls qsolve with signs all 1.
        """
        pass


    def graver(self):
        """
        Calls zsolve with
        """
        pass

    def hilbert(self):
        """
        """
        pass



class LinearSystemOverZZ(LinearSystem):
    def solve(self, bound=None):
        if not bound is None:
            # the extensions will be .ub and .lb
            raise NotImplementedError

        out, err = self.call_4ti2('zsolve')
        zhom = self._read_file('zhom')
        zinhom = self._read_file('zinhom')
        zfree = self._read_file('zfree')

        return (zhom, zinhom, zfree)

class LinearSystemOverQQ(LinearSystem):
    pass


linear_system_over_ZZ = LinearSystemOverZZ

linear_system_over_QQ = LinearSystemOverQQ


class ToricIdeal:
    pass

