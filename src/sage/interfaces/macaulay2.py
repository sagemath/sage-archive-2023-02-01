r"""
Interface to Macaulay2

\note{You must have \code{Macaulay2} installed on your computer
for this interface to work. Macaulay2 is not included with \sage,
but you can obtain it from \url{http://www.math.uiuc.edu/Macaulay2/}.
Note additional optional \sage packages are required.}

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
   -- Kiran Kedlaya (2006-02-12): added ring and ideal constructors,
             list delimiters, is_Macaulay2Element, sage_polystring,
             __floordiv__, __mod__, __iter__, __len__; stripped extra
             leading space and trailing newline from output.

TODO:
   -- get rid of all numbers in output, e.g., in ideal function below.
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

import os

from expect import Expect, ExpectElement

from sage.misc.misc import verbose

# from re import search

def clean_output(s):
    i = s.find('= ')
    if i == -1:
        return ''
    return s[i+2:-1]

class Macaulay2(Expect):
    """
    Interface to the Macaulay2 interpreter.
    """
    def __init__(self, maxread=10000, script_subdirectory="",
                 logfile=None, server=None):
        Expect.__init__(self,
                        name = 'macaulay2',
                        prompt = 'i[0-9]* : ',
                        command = "M2 --no-debug --no-readline --silent ",
                        maxread = maxread,
                        server = server,
                        script_subdirectory = script_subdirectory,
                        verbose_start = False,
                        logfile = None,
                        eval_using_file_cutoff=500)

    # Macaulay2 provides no "clear" function. However, Macaulay2 does provide
    # garbage collection; since expect automatically reuses variable names,
    # garbage collection in SAGE properly sets up garbage collection in
    # Macaulay2.

    def _read_in_file_command(self, filename):
        return 'value get "%s"'%filename

    def eval(self, code):
        """
        Send the code x to the Macaulay2 interpreter and return the output
        as a string suitable for input back into Macaulay2, if possible.
        """
        # TODO: in some cases change toExternalString to toString??
        ans = Expect.eval(self, 'toExternalString(%s)'%code)
        if ans.find("stdio:") != -1:
            if 'cannot be converted to external string' in ans:
                return clean_output(Expect.eval(self, '%s'%code))
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
        ans = Expect.eval(self, cmd)
        if ans.find("stdio:") != -1:
            raise RuntimeError, "Error evaluating Macaulay2 code.\nIN:%s\nOUT:%s"%(cmd, ans)

    def _object_class(self):
        return Macaulay2Element

    def console(self):
        macaulay2_console()

    def _left_list_delim(self):
        return '{'

    def _right_list_delim(self):
        return '}'

    def _true_symbol(self):
        return 'true'

    def _false_symbol(self):
        return 'false'

    def _equality_symbol(self):
        return '=='


### Constructors

    def ideal(self, *gens):
        """
        Return the ideal generated by gens.

        INPUT:
            gens -- list or tuple of Macaulay2 objects (or objects that can be
                    made into Macaulay2 objects via evaluation)
        OUTPUT:
            the Macaulay2 ideal generated by the given list of gens

        EXAMPLES:
            sage: R2 = macaulay2.ring('QQ', '[x, y]'); R2            # optional
            QQ [x, y, MonomialOrder => Lex, MonomialSize => 16]
            sage: I = macaulay2.ideal( ('y^2 - x^3', 'x - y') ); I   # optional
            ideal (- x  + y , x - y)
            <BLANKLINE>
            o17 : Ideal of sage2
            sage: J = I^3; J.gb()                                    # optional
            GroebnerBasis[status: done; S-pairs encountered up to degree 9]
            <BLANKLINE>
            o21 : GroebnerBasis
            sage: J.gb().generators()                                # optional
            map(sage2^{{0}}, sage2^{{-9}, {-7}, {-5}, {-3}}, {{y^9-3*y^8+3*y^7-y^6, x*y^6-2*x*y^5+x*y^4-y^7+2*y^6-y^5, x^2*y^3-x^2*y^2-2*x*y^4+2*x*y^3+y^5-y^4, x^3-3*x^2*y+3*x*y^2-y^3}})
        """
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        gens2 = []
        for g in gens:
            if not isinstance(g, Macaulay2Element):
                gens2.append(self(g))
            else:
                gens2.append(g)
        return self('ideal {%s}'%(",".join([g.name() for g in gens2])))

    def ring(self, base_ring='ZZ', vars='[x]', order='Lex'):
        r"""
        Create a Macaulay2 ring.

        INPUT:
            base_ring -- base ring (see examples below)
            vars -- a tuple or string that defines the variable names
            order -- string -- the monomial order (default: 'Lex')

        OUTPUT:
            a Macaulay2 ring (with base ring ZZ)

        EXAMPLES:
        This is a ring in variables named a through d over the finite field
        of order 7, with graded reverse lex ordering:
            sage: R1 = macaulay2.ring('ZZ/7', '[a..d]', 'GRevLex'); R1  # optional
            ZZ/7 [a, b, c, d, MonomialSize => 16]
            sage: R1.char()                                             # optional
            7

        This is a polynomial ring over the rational numbers:
            sage: R2 = macaulay2.ring('QQ', '[x, y]'); R2               # optional
            QQ [x, y, MonomialOrder => Lex, MonomialSize => 16]
        """
        varstr = str(vars)[1:-1]
        return self.new('%s[%s, MonomialSize=>16, MonomialOrder=>%s]'%(
            base_ring, varstr, order))

class Macaulay2Element(ExpectElement):

    def __iter__(self):
        for i in range(len(self)):  # zero-indexed!
            yield self[i]

    def __len__(self):
        self._check_valid()
        return int(self.parent()("#%s"%self.name()))

    def __getitem__(self, n):
        self._check_valid()
        return self.parent().new('%s#%s'%(self.name(), n))

    def __call__(self, x):
        self._check_valid()
        P = self.parent()
        r = P(x)
        return P('%s %s'%(self.name(), r.name()))

    def gen(self, n):
        self._check_valid()
        return self.parent().new('%s_%s'%(self._name, int(n)))

    def __floordiv__(self, x):
        """
        Quotient of division of self by other.  This is denoted //.

        EXAMPLE:
            sage: R = PolynomialRing(GF(7), 2, 'xy', macaulay2=True)   # optional
            sage: x, y = R.gens()                                      # optional
            sage: f = (x^3 + 2*y^2*x)^7; f                             # optional
            2*x^7*y^14 + x^21
            sage: h = macaulay2(f); h                                  # optional
            x^21+2*x^7*y^14
            sage: f1 = (x^2 + 2*y*x)                                   # optional
            sage: h1 = macaulay2(f1)                                   # optional
            sage: f2 = (x^3 + 2*y*x)                                   # optional
            sage: h2 = macaulay2(f2)                                   # optional
            sage: h % [h1,h2]                                          # optional
            -3*x*y
            sage: u = h // [h1,h2]; u                                  # optional
            [x^19-2*x^18*y-3*x^17*y^2-x^16*y^3+2*x^15*y^4+3*x^14*y^5+x^13*y^6-2*x^12*y^7-3*x^11*y^8-x^10*y^9+2*x^9*y^10+3*x^8*y^11+x^7*y^12-2*x^6*y^13-x^5*y^14+2*x^4*y^15+3*x^3*y^16+x^2*y^17-x*y^17-3*x*y^16-2*x*y^15+x*y^14+3*x*y^13+2*x*y^12-x*y^11-3*x*y^10-2*x*y^9+x*y^8+3*x*y^7+2*x*y^6-x*y^5-3*x*y^4-2*x*y^3+x*y^2+3*x*y+2*x+2*y^18-y^17-3*y^16-2*y^15+y^14+3*y^13+2*y^12-y^11-3*y^10-2*y^9+y^8+3*y^7+2*y^6-y^5-3*y^4-2*y^3+y^2+3*y, -2*y^18+y^17+3*y^16+2*y^15-y^14-3*y^13-2*y^12+y^11+3*y^10+2*y^9-y^8-3*y^7-2*y^6+y^5+3*y^4+2*y^3-y^2-3*y-2]
            sage: h == u[0]*h1 + u[1]*h2 + (h % [h1,h2])               # optional
            True
        """
        if isinstance(x, (list, tuple)):
            y = self.parent(x)
            z = self.parent().new('%s // matrix{%s}'%(self.name(), y.name()))
            return list(z.entries().flatten())
        else:
            return self.parent().new('%s // %s'%(self.name(), x.name()))

    def __mod__(self, x):
        """
        Remainder of division of self by other.  This is denoted %.

        EXAMPLE:
            sage: R = PolynomialRing(GF(7), 2, ['x','y'], macaulay2=True)   # optional
            sage: x, y = R.gens()                                           # optional
            sage: f = (x^3 + 2*y^2*x)^7; f                                  # optional
            2*x^7*y^14 + x^21
            sage: h= f._macaulay2_(); h                                     # optional
            x^21+2*x^7*y^14
            sage: f1 = (x^2 + 2*y*x)                                        # optional
            sage: h1 = f1._macaulay2_()                                     # optional
            sage: f2 = (x^3 + 2*y*x)                                        # optional
            sage: h2 = f2._macaulay2_()                                     # optional
            sage: h % [h1,h2]                                               # optional
            -3*x*y
            sage: u = h // [h1,h2]                                          # optional
            sage: h == u[0]*h1 + u[1]*h2 + (h % [h1,h2])                    # optional
            True
        """
        if isinstance(x, (list, tuple)):
            y = self.parent(x)
            return self.parent().new('%s %% matrix{%s}'%(self.name(), y.name()))
        else:
            return self.parent().new('%s %% %s'%(self.name(), x.name()))

    def is_zero(self):
        P = self.parent()
        return P.eval('%s == 0'%self.name()) == 'true'

    def sage_polystring(self):
	"""
	If this Macaulay2 element is a polynomial, return a string
	representation of this polynomial that is suitable for
	evaluation in Python.  Thus * is used for multiplication
	and ** for exponentiation.   This function is primarily
	used internally.

	EXAMPLES:
            sage: R = macaulay2.ring('QQ','(x,y)')               # optional
            sage: f = macaulay2('x^3 + 3*y^11 + 5')              # optional
            sage: f                                              # optional
            x^3+3*y^11+5
            sage: f.sage_polystring()                            # optional
            'x**3+3*y**11+5'
	"""
        return str(self).replace('^','**')

    def structure_sheaf(self):
        """
        EXAMPLES:
            sage: S = macaulay2('QQ[a..d]')                     # optional
            sage: R = S/macaulay2('a^3+b^3+c^3+d^3')            # optional
            sage: X = R.Proj()                                  # optional
            sage: X.structure_sheaf()                           # optional
            new SheafOfRings from {symbol ring => (QQ [a, b, c, d])/(a^3+b^3+c^3+d^3),
                 symbol variety => sage0}
        """
        return self.parent()('OO_%s'%self.name())

def is_Macaulay2Element(x):
    return isinstance(x, Macaulay2Element)

# An instance
macaulay2 = Macaulay2(script_subdirectory='user')

import os, sys

def macaulay2_console():
    os.system('M2')




