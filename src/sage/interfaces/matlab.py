r"""
Interface to MATLAB

According to their website, MATLAB is "a high-level language and
interactive environment that enables you to perform computationally
intensive tasks faster than with traditional programming languages
such as C, C++, and Fortran."

The commands in this section only work if you have the
"matlab" interpreter installed and available in
your PATH.  It's not necessary to install any special
\sage packages.

EXAMPLES:
    sage: matlab.eval('2+2')
    'ans = 4'

    sage: a = matlab(10)
    sage: a**10
    1e+10

AUTHORS:
   -- William Stein (2006-10-11)

\subsection{Tutorial}
EXAMPLES:
    sage: matlab('4+10')
    14
    sage: matlab('date')
    11-Oct-2006
    sage: matlab('5*10 + 6')
    56
    sage: matlab('(6+6)/3')
    4
    sage: matlab('9')^2
    81
    sage: a = matlab(10); b = matlab(20); c = matlab(30)
    sage: avg = (a+b+c)/3
    sage: avg
    20
    sage: parent(avg)
    Matlab

    sage: my_scalar = matlab('3.1415')
    sage: my_scalar
    3.1415
    sage: my_vector1 = matlab('[1,5,7]')
    sage: my_vector1
    1     5     7
    sage: my_vector2 = matlab('[1;5;7]')
    sage: my_vector2
    1
    5
    7
    sage: my_vector1 * my_vector2
    75

    sage: row_vector1 = matlab('[1 2 3]')
    sage: row_vector2 = matlab('[3 2 1]')
    sage: matrix_from_row_vec = matlab('[%s; %s]'%(row_vector1.name(), row_vector2.name()))
    sage: matrix_from_row_vec
    1     2     3
    3     2     1

    sage: column_vector1 = matlab('[1;3]')
    sage: column_vector2 = matlab('[2;8]')
    sage: matrix_from_col_vec = matlab('[%s %s]'%(column_vector1.name(), column_vector2.name()))
    sage: matrix_from_col_vec
    1     2
    3     8

    sage: my_matrix = matlab('[8, 12, 19; 7, 3, 2; 12, 4, 23; 8, 1, 1]')
    sage: my_matrix

     8    12    19
     7     3     2
    12     4    23
     8     1     1

    sage: combined_matrix = matlab('[%s, %s]'%(my_matrix.name(), my_matrix.name()))
    sage: combined_matrix
     8    12    19     8    12    19
     7     3     2     7     3     2
    12     4    23    12     4    23
     8     1     1     8     1     1

    sage: tm = matlab('0.5:2:10')
    sage: tm
    0.5000    2.5000    4.5000    6.5000    8.5000

    sage: my_vector1 = matlab('[1,5,7]')
    sage: my_vector1(1)
    1
    sage: my_vector1(2)
    5
    sage: my_vector1(3)
    7

Matrix indexing works as follows:
    sage: my_matrix = matlab('[8, 12, 19; 7, 3, 2; 12, 4, 23; 8, 1, 1]')
    sage: my_matrix(3,2)
    4


Setting using paranthesis cannot work (because of how the Python language
works).  Use square brackets or the set function:

sage: my_matrix = matlab('[8, 12, 19; 7, 3, 2; 12, 4, 23; 8, 1, 1]')
sage: my_matrix.set(2,3, 1999)
sage: my_matrix
           8          12          19
           7           3        1999
          12           4          23
           8           1           1


"""

##############################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL), Version 2.
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

import os

from expect import Expect, ExpectElement

from sage.misc.misc import verbose

#import sage.matrix.matrix_space

class Matlab(Expect):
    """
    Interface to the Matlab interpreter.

    EXAMPLES:
        sage: a = matlab('[ 1, 1, 2; 3, 5, 8; 13, 21, 33 ]')    # optional
        sage: b = matlab('[ 1; 3; 13]')                         # optional
        sage: c = a * b                                         # optional
        sage: print c                                           # optional
    """
    def __init__(self, maxread=100, script_subdirectory="",
                 logfile=None, server=None):
        Expect.__init__(self,
                        name = 'matlab',
                        prompt = '>> ',
                        command = "matlab -nodisplay",
                        maxread = maxread,
                        server = server,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=100)

    def __reduce__(self):
        return reduce_load_Matlab, tuple([])

    def _read_in_file_command(self, filename):
        return 'source("%s");'%filename

    def _quit_string(self):
        return 'quit;'

    def _install_hints(self):
        return """
        You must obtain the program MATLAB in order to use MATLAB
        from SAGE.   You can read all about MATLAB at
                  http://www.mathworks.com/

        You might have to buy MATLAB (list price: $1900).
        """

    def _start(self):
        Expect._start(self)

    def whos(self):
        return self.eval('whos')

    def get_via_file(self, var_name):
        t = self._temp_file(var_name)
        self.eval('save -text "%s" %s'%(t,var_name))
        r = open(t).read()
        os.unlink(t)
        return r.strip('\n')

    def set_via_file(self, var_name, x):
        t = self._temp_file(var_name)
        open(t,'w').write(x)
        print 'load "%s" %s'%(t, var_name)
        self.eval('load "%s" %s'%(t, var_name))
        #os.unlink(t)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s;'%(var,value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError, "Error executing code in Matlab\nCODE:\n\t%s\nMatlab ERROR:\n\t%s"%(cmd, out)

    def get(self, var):
        """
        Get the value of the variable var.
        """
        s = self.eval('%s'%var)
        i = s.find('=')
        return s[i+1:].strip('\n')

    def console(self):
        matlab_console()

    def version(self):
        return matlab_version()[1:]

    def sage2matlab_matrix_string(self, A):
        """
        Return an matlab matrix from a SAGE matrix.

        INPUT:
            A SAGE matrix with entries in the rationals or reals.

        OUTPUT:
            A string that evaluates to an Matlab matrix.

        EXAMPLES:
            sage: M33 = MatrixSpace(QQ,3,3)
            sage: A = M33([1,2,3,4,5,6,7,8,0])
            sage: matlab.sage2matlab_matrix_string(A)   # requires optional matlab
            '[1, 2, 3; 4, 5, 6; 7, 8, 0]'

        AUTHOR: David Joyner and William Stein
        """
        return str(A.rows()).replace('), (', '; ').replace('(', '').replace(')','')

    def _object_class(self):
        return MatlabElement


class MatlabElement(ExpectElement):
    def __getitem__(self, n):
        raise RuntimeError, "Use parenthesis for MATLAB matrices instead."

    def _matrix_(self, R):
        r"""
        Return \sage matrix from this matlab element.

        EXAMPLES:
            sage: A = matlab('[1,2;3,4]')       # optional matlab package
            sage: matrix(ZZ, A)
            [1 2]
            [3 4]
            sage: A = matlab('[1,2;3,4.5]')     # optional matlab package
            sage: matrix(RR, A)
            [1.0000000000000000 2.0000000000000000]
            [3.0000000000000000 4.5000000000000000]
        """
        from sage.matrix.all import MatrixSpace
        s = str(self).strip()
        v = s.split('\n ')
        nrows = len(v)
        if nrows == 0:
            return MatrixSpace(R,0,0)(0)
        ncols = len(v[0].split())
        M = MatrixSpace(R, nrows, ncols)
        v = sum([[x for x in w.split()] for w in v], [])
        return M(v)

    def set(self, i, j, x):
        P = self._check_valid()
        z = P(x)
        P.eval('%s(%s,%s) = %s'%(self.name(), i, j, z.name()))

# An instance
matlab = Matlab(script_subdirectory='user')

def reduce_load_Matlab():
    return matlab


import os
def matlab_console():
    """
    This requires that the optional matlab program be installed and in
    your PATH, but no optional \sage packages need be installed.

    EXAMPLES:
        sage.: matlab_console()
                                       < M A T L A B >
                           Copyright 1984-2006 The MathWorks, Inc.
        ...
        >> 2+3

        ans =

             5

        >> quit

    Typing quit exits the matlab console and returns you to SAGE.
    matlab, like SAGE, remembers its history from one session to
    another.
    """
    os.system('matlab -nodisplay')


def matlab_version():
    """
    Return the version of Matlab installed.

    EXAMPLES:
        sage: matlab_version()    # optional matlab package
        '2.1.73'
    """
    return str(matlab('version')).strip()

