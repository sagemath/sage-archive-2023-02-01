r"""
Interface to Macaulay2

\note{You must have \code{Macaulay2} installed on your computer
for this interface to work. Macaulay2 is not included with \sage,
but you can obtain it from \url{http://www.math.uiuc.edu/Macaulay2/}.}

SAGE provides an interface to the Macaulay2 computational algebra
system. This system provides extensive functionality for commutative
algebra. You do not have to install any optional packages.

The Macaulay2 interface offers three pieces of functionality:
\begin{enumerate}

\item \code{Macaulay2_console()} -- A function that dumps you
into an interactive command-line Macaulay2 session.

\item \code{Macaulay2(expr)} -- Evaluation of arbitrary Macaulay2
expressions, with the result returned as a string.

\item \code{Macaulay2.new(expr)} -- Creation of a SAGE object that wraps a
Macaulay2 object.  This provides a Pythonic interface to Macaulay2.  For
example, if \code{f=Macaulay2.new(10)}, then \code{f.gcd(25)} returns the
GCD of $10$ and $25$ computed using Macaulay2.

\end{enumerate}

EXAMPLES:
    sage: macaulay2('3/5 + 7/11')
    68/55
    sage: f = macaulay2('f = i -> i^3')
    sage: f
    f
    sage: f(5)
    125

    sage: R = macaulay2('ZZ/5[x,y,z]')
    sage: R
    ZZ/5 [x, y, z]
    sage: x = macaulay2('x')
    sage: y = macaulay2('y')
    sage: (x+y)^5
    x^5+y^5
    sage: parent((x+y)^5)
    Macaulay2

    sage: R = macaulay2('QQ[x,y,z,w]')
    sage: f = macaulay2('x^4 + 2*x*y^3 + x*y^2*w + x*y*z*w + x*y*w^2 + 2*x*z*w^2 + y^4 + y^3*w + 2*y^2*z*w + z^4 + w^4')
    sage: f
    x^4+2*x*y^3+y^4+z^4+x*y^2*w+y^3*w+x*y*z*w+2*y^2*z*w+x*y*w^2+2*x*z*w^2+w^4

    sage: g = f * macaulay2('x+y^5')
    sage: g.factor()
      new Product from {new Power from
      {x^4+2*x*y^3+y^4+z^4+x*y^2*w+y^3*w+x*y*z*w+2*y^2*z*w+x*y*w^2+2*x*z*w^2+w^4, 1},
      new Power from {y^5+x, 1}}


AUTHORS:
   -- Kiran Kedlaya and David Roe (2006-02-05, during SAGE coding sprint)
   -- William Stein (2006-02-09): inclusion in SAGE; prompt uses regexp,
             calling of Macaulay2 functions via __call__.
   -- William Stein (2006-02-09): fixed bug in reading from file and
             improved output cleaning.

"""

#*****************************************************************************
#       Copyright (C) 2006 Kiran S. Kedlaya <kedlaya@mit.edu>
#                          David Roe <roed@mit.edu>
#                          William Stein <wstein@ucsd.edu>
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
#*****************************************************************************

#
# NOTE TO SELVES: it would be nice if this interface implemented a bit
# of the Macaulay2 class hierarchy, in a manner similar to the Singular
# interface.
#

import os

from expect import Expect, ExpectElement

from sage.misc.misc import verbose

from re import search

def clean_output(s):
    i = s.find('=')
    if i == -1:
        return ''
    return s[i+1:]

    #rex1 = 'o[0-9]+ = (.*)\\n\\no[0-9]+ : .*\\n'
    #rex2 = 'o[0-9]+ = (.*)\\n'
    #m = search(rex1, s)
    #if not m is None:
    #    return m.group(1)
    #m = search(rex2, s)
    #if not m is None:
    #    return m.group(1)
    #return ""

class Macaulay2(Expect):
    """
    Interface to the Macaulay2 interpreter.
    """
    def __init__(self, maxread=10000, script_subdirectory="", logfile=None, server=None):
        Expect.__init__(self,
                        name = 'Macaulay2',
                        prompt = 'i[0-9]* : ',
                        command = "M2",
                        maxread = maxread,
                        server = server,
                        script_subdirectory = script_subdirectory,
                        verbose_start = False,
                        logfile=None,
                        eval_using_file_cutoff=500)

    # no appropriate clear command in Macaulay2

    def _read_in_file_command(self, filename):
        return 'value get "%s"'%filename

    def eval(self, code):
        """
        Send the code x to the Macaulay2 interpreter and return the output
        as a string suitable for input back into Macaulay2, if possible.
        """
        ans = Expect.eval(self, 'toExternalString(%s)'%code)
        if ans.find("stdio:") != -1:
            raise RuntimeError, "Error evaluating Macaulay2 code.\nIN:%s\nOUT:%s"%(code, ans)
        return clean_output(ans)

    def get(self, var):
        """
        Get the value of the variable var.
        """
        return self.eval("%s"%var)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s;'%(var,value)
        ans = self.eval(cmd)
        #ans = self._eval_line(cmd, allow_use_file=True)
        if ans.find("stdio:") != -1:
            raise RuntimeError, "Error evaluating Macaulay2 code.\nIN:%s\nOUT:%s"%(cmd, ans)

    def _object_class(self):
        return Macaulay2Element

    def console(self, readline=True):
        Macaulay2_console(readline=readline)

class Macaulay2Element(ExpectElement):
    def __getitem__(self, n):
        return self.parent().new('%s#%s'%(self._name, n))

    def __call__(self, x):
        self._check_valid()
        P = self.parent()
        r = P(x)
        return P('%s %s'%(self.name(), r.name()))



# An instance
macaulay2 = Macaulay2(script_subdirectory='user')

# Cleverly run Macaulay2 with the benefit of readline, which
# is something the usual Macaulay2 doesn't provide!
# See
#    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/363500

import os, sys
def generic_console(command, readline=True):
    if not readline:
        os.system(command)
        return
    f1 = os.popen(command + ' ', 'w')
    f1.flush()
    try:
        while True:
            sys.stdout.write('')
            try:
                line = raw_input('')
                f1.writelines(line+'\n')
                f1.flush()
            except KeyboardInterrupt:
                f1.close()
                break
    except EOFError:
        pass
    sys.stdout.write('\n')

def macaulay2_console(readline=True):
    generic_console('M2', readline)
