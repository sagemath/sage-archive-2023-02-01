r"""
Interface to Axiom

TODO:
   * Evaluation using a file is not done.   Any input
     line with more than a few thousand characters would
     hang the system, so currently it automatically raises
     an exception.
   * All completions of a given command.
   * Interactive help.

Axiom is a free GPL-compatible (modified BSD license)  general purpose
computer algebra system whose development started in 1973 at IBM.  It
contains symbolic manipulation algorithms, as well as implementations
of special functions, including elliptic functions and generalized
hypergeometric functions. Moreover, Axiom has implementations of many
functions relating to the invariant theory of the symmetric group $S_n$.
For many links to Axiom documentation see
         \url{http://wiki.axiom-developer.org}.

AUTHORS OF THIS MODULE:
    -- Bill Page (2006-10): Created this (based on maxima interfac)
       NOTE: Bill Page put a huge amount of effort into the SAGE Axiom interface
             over several days during the SAGE Days 2 coding sprint.  This is
             contribution is greatly appreciated.
    -- William Stein (2006-10): misc touchup.

If the string "error" (case insensitive) occurs in the output of
anything from axiom, a RuntimeError exception is raised.

EXAMPLES:
We evaluate a very simple expression in axiom.
    sage: axiom('3 * 5')                     # optional
    15
    sage: a = axiom(3) * axiom(5); a         # optional
    15

The type of a is AxiomElement, i.e., an element of the axiom interpreter.
    sage: type(a)                            # optional
    <class 'sage.interfaces.axiom.AxiomElement'>
    sage: parent(a)                          # optional
    Axiom

The underlying Axiom type of a is also available, via the type method:
    sage: a.type()                           # optional
    PositiveInteger


We factor $x^5 - y^5$ in Axiom in several different ways.
The first way yields a Axiom object.
    sage: F = axiom.factor('x^5 - y^5'); F      # optional
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    sage: type(F)                               # optional
    <class 'sage.interfaces.axiom.AxiomElement'>
    sage: F.type()                              # optional
    Factored Polynomial Integer


Note that Axiom objects are normally displayed using ``ASCII art''.
In some cases you can see a normal linear representation of any Axiom
object x, using \code{str(x)}.  This can be useful for moving axiom
data to other systems.
    sage: a = axiom('2/3'); a          # optional
    2
    -
    3
    sage: str(a)                       # optional
    '2/3'
    sage: a = axiom('x^2 + 3/7')       # optional
    sage: str(a)                       # optional
    'x*x+3/7'

The \code{axiom.eval} command evaluates an expression in axiom and
returns the result as a string.  This is exact as if we typed in the
given line of code to axiom; the return value is what Axiom would
print out.

    sage: print axiom.eval('factor(x^5 - y^5)')   # optional
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    <BLANKLINE>
    Type: Factored Polynomial Integer

We can create the polynomial $f$ as a Axiom polynomial, then call
the factor method on it.  Notice that the notation \code{f.factor()}
is consistent with how the rest of \sage works.
    sage: f = axiom('x^5 - y^5')                  # optional
    sage: f^2                                     # optional
     10     5 5    10
    y   - 2x y  + x
    sage: f.factor()                              # optional
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

Control-C interruption works well with the axiom interface,
because of the excellent implementation of axiom.  For example,
try the following sum but with a much bigger range, and hit
control-C.
    sage.:  f = axiom('(x^5 - y^5)^10000')       # optional
    Interrupting Axiom...
    ...
    <type 'exceptions.TypeError'>: Ctrl-c pressed while running Axiom

Symbolic constants:
    sage: axiom(e + pi)                      # optional
    %e + %pi
"""

#\subsection{Tutorial}
#We follow the tutorial at
#\url{http://wiki.axiom-developer.org/AxiomTutorial}.
#\subsection{Interactivity}
#Axiom has a non-interactive mode that is initiated via the command
#")read file.input".
#\subsection{Long Input}
#The Axiom interface reads in even very long input (using files) in a
#robust manner, as long as you are creating a new object.
#\note{Using \code{axiom.eval} for long input
#is much less robust, and is not recommended.}
#
#    sage: t = '"%s"'%10^10000   # ten thousand character string.    (optional)
#    sage: a = axiom(t)                                              # optional


###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
###########################################################################

import os, re

from expect import Expect, ExpectElement, FunctionElement, ExpectFunction, tmp
from pexpect import EOF

from sage.misc.misc import verbose, DOT_SAGE, SAGE_ROOT

from sage.misc.multireplace import multiple_replace

cnt = 0
seq = 0

COMMANDS_CACHE = '%s/axiom_commandlist_cache.sobj'%DOT_SAGE

# The Axiom commands ")what thing det" ")show Matrix" and ")display
# op det" commands, gives a list of all identifiers that begin in
# a certain way.  This could maybe be useful somehow... (?)  Also
# axiom has a lot a lot of ways for getting documentation from the
# system -- this could also be useful.

class Axiom(Expect):
    """
    Interface to the Axiom interpreter.
    """
    def __init__(self, script_subdirectory=None, logfile=None, server=None):
        """
        Create an instance of the Axiom interpreter.
        """
        eval_using_file_cutoff = 200
        self.__eval_using_file_cutoff = eval_using_file_cutoff
        Expect.__init__(self,
                        name = 'axiom',
                        prompt = '\([0-9]+\) -> ',
                        command = "axiom -nox -noclef",
                        maxread = 10,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        init_code = [')lisp (si::readline-off)'],
                        logfile = logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return AxiomExpectFunction(self, attrname)

    def _start(self):
        Expect._start(self)
        self._expect.expect(self._prompt)
        out = self._eval_line(')set functions compile on', reformat=False)
        out = self._eval_line(')set output length 245', reformat=False)
        out = self._eval_line(')set message autoload off', reformat=False)
        self._expect.expect(self._prompt)

    def _eval_line_using_file(self, line, tmp):
        F = open(tmp, 'w')
        F.write(line)
        F.close()
        if self._expect is None:
            self._start()
        # For some reason this trivial comp
        # keeps certain random freezes from occuring.  Do not remove this.
        # The space before the \n is also important.
        self._expect.sendline(')read "%s"\n'%tmp)
        self._expect.expect(self._prompt)
        return ''

    def __reduce__(self):
        return reduce_load_Axiom, tuple([])

    def _quit_string(self):
        return ')lisp (quit)'

    def _eval_line(self, line, reformat=True, allow_use_file=False,
                   wait_for_prompt=True):
        if not wait_for_prompt:
            return Expect._eval_line(self, line)
        line = line.rstrip().rstrip(';')
        if line == '':
            return ''
        if len(line) > 3000:
            raise NotImplementedError, "evaluation of long input lines (>3000 characters) in Axiom not yet implemented."
        global seq
        seq += 1
        if self._expect is None:
            self._start()
        if allow_use_file and self.__eval_using_file_cutoff and \
                            len(line) > self.__eval_using_file_cutoff:
            return self._eval_line_using_file(line, tmp)
        try:
            E = self._expect
            # debug
            verbose("in = '%s'"%line,level=3)
            E.sendline(line)
            self._expect.expect(self._prompt)
            out = self._expect.before
            # debug
            verbose("out = '%s'"%out,level=3)
        except EOF:
          if self._quit_string() in line:
             return ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()

        if '>> Error detected within library code:' in out or \
           'Cannot find a definition or applicable library operation named' in out:
            raise RuntimeError, out

        if not reformat:
            return out
        if 'error' in out:
            return out
        #out = out.lstrip()
        i = out.find('\n')
        out = out[i+1:]
        outs = out.split("\n")
        i = 0
        outline = ''
        for line in outs:
            line = line.rstrip()
            # print "'%s'"%line
            if line[:4] == '   (':
                i = line.find('(')
                i += line[i:].find(')')
                if line[i+1:] == "":
                    i = 0
                    outs = outs[1:]
                break;
        out = "\n".join(line[i+1:] for line in outs[1:])
        return out

    ###########################################
    # Interactive help
    ###########################################

    def help(self, s):
        import sage.server.support
        if sage.server.support.EMBEDDED_MODE:
            e = os.system('asq -op "%s"< /dev/null'%s)
        else:
            e = os.system('asq -op "%s"'%s)
        if e:
            print "Help system not available."

    def example(self, s):
        import sage.server.support
        if sage.server.support.EMBEDDED_MODE:
            e = os.system('asq -doc "%s" < /dev/null'%s)
        else:
            e = os.system('asq -doc "%s"'%s)
        if e:
            print "Help system not available."

    describe = help

    def demo(self):
        import sage.server.support
        if sage.server.support.EMBEDDED_MODE:
            os.system('axiom -ht < /dev/null')
        else:
            os.system('axiom -ht')

    def completions(self, s):
        """
        Return all commands that complete the command starting with the
        string s.   This is like typing s[tab] in the maple interpreter.
        """
        s = self.eval('apropos(%s)'%s).replace('\\ - ','-')
        return [x for x in s[1:-1].split(',') if x[0] != '?']

    def _commands(self):
        """
        Return list of all commands defined in Axiom.
        """
        try:
            return self.__commands
        except AttributeError:
            self.__commands = sum([self.completions(chr(97+n)) for n in range(26)], [])
        return self.__commands

    def _object_class(self):
        return AxiomElement

    def _true_symbol(self):
        return 'true'

    def _false_symbol(self):
        return 'false'

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s := %s'%(var, value)
        out = self._eval_line(cmd, reformat=False)

        if out.find("error") != -1:
            raise TypeError, "Error executing code in Axiom\nCODE:\n\t%s\nAxiom ERROR:\n\t%s"%(cmd, out)


    def get(self, var):
        """
        Get the string value of the Axiom variable var.
        """
        s = self._eval_line(str(var))
        i = s.rfind('Type:')
        return s[:i].rstrip()

    def console(self):
        axiom_console()


class AxiomElement(ExpectElement):
    def __call__(self, x):
        self._check_valid()
        P = self.parent()
        return P('%s(%s)'%(self.name(), x))

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: a = axiom(1); b = axiom(2)                  # optional
            sage: a == b                                      # optional
            False
            sage: a < b                                       # optional
            True
            sage: a > b                                       # optional
            False
            sage: b < a                                       # optional
            False
            sage: b > a                                       # optional
            True

        We can also compare more complicated object such as functions:
            sage: f = axiom('sin(x)'); g = axiom('cos(x)')    # optional
            sage: f == g                                      # optional
            False
        """
        if not isinstance(other, AxiomElement):
            return -1
        P = self.parent()
        t = P._true_symbol()
        if P('%s = %s'%(self.name(), other.name())).__repr__().strip() == t:
            return 0
        elif P('%s > %s'%(self.name(), other.name())).__repr__().strip() == t:
            return 1
        else:
            return -1  # everything is supposed to be comparable in Python, so we define
                       # the comparison thus when no comparable in interfaced system.
    def numer(self):
        P = self.parent()
        return P('numeric(%s)'%self._name)

    def real(self):
        return self.realpart()

    def imag(self):
        return self.imagpart()

    def str(self):
        """
        Get the linear string representation of this object, if possible (often it isn't).
        """
        P = self._check_valid()
        s = P.eval('unparse(%s::InputForm)'%self._name)
        if 'translation error' in s:
            raise RuntimeError, s
        s = multiple_replace({'\r\n':'', # fix stupid Fortran-ish
                              'DSIN(':'sin(',
                              'DCOS(':'cos(',
                              'DTAN(':'tan(',
                              'DSINH(':'sinh('}, s)
        return re.search(r'"(.*)"',s).groups(0)[0]

    def __repr__(self):
        P = self._check_valid()
        return P.get(self._name)

    def __str__(self):
        return self.str()

    def type(self):
        P = self._check_valid()
        s = P._eval_line(self.name())
        i = s.rfind('Type:')
        return AxiomType(s[i+5:].strip())

    def __float__(self):
        return float(str(self.numer()))

    def __len__(self):
        """
        Return the length of a list.

        EXAMPLES:
            sage: v = axiom('[x^i for i in 0..5]')            # optional
            sage: len(v)                                      # optional
            6
        """
        P = self._check_valid()
        s = P.eval('# %s '%self.name())
        i = s.rfind('Type')
        return int(s[:i-1])

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return AxiomFunctionElement(self, attrname)

    def __getitem__(self, n):
        r"""
        Return the n-th element of this list.

        \note{Lists are 1-based.}

        EXAMPLES:
            sage: v = axiom('[i*x^i for i in 0..5]'); v          # optional
                   2   3   4   5
            [0,x,2x ,3x ,4x ,5x ]
            sage: v[4]                                           # optional
              3
            3x
            sage: v[1]                                           # optional
            0
            sage: v[10]                                          # optional
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        n = int(n)
        if n <= 0 or n > len(self):
            raise IndexError, "index out of range"
        P = self._check_valid()
        if not isinstance(n, tuple):
            return P.new('%s(%s)'%(self._name, n))
        else:
            return P.new('%s(%s)'%(self._name, str(n)[1:-1]))


    def subst(self, val):
        P = self.parent()
        return P('subst(%s, %s)'%(self.name(), val))

    def comma(self, args):
        self._check_valid()
        P = self.parent()
        return P('%s, %s'%(self.name(), args))

    def _latex_(self):
        self._check_valid()
        P = self.parent()
        s = axiom._eval_line('outputAsTex(%s)'%self.name(), reformat=False)
        if not '$$' in s:
            raise RuntimeError, "Error texing axiom object."
        i = s.find('$$')
        j = s.rfind('$$')
        s = s[i+2:j]
        s = multiple_replace({'\r':'', '\n':' ',
                              ' \\sp ':'^',
                              '\\arcsin ':'\\sin^{-1} ',
                              '\\arccos ':'\\cos^{-1} ',
                              '\\arctan ':'\\tan^{-1} '},
            re.sub(r'\\leqno\(.*?\)','',s)) # no eq number!
        return s


class AxiomFunctionElement(FunctionElement):
    def _sage_doc_(self):
        return self._obj.parent().help(self._name)

class AxiomExpectFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)

class AxiomType:
    def __init__(self, x):
        self.__x = str(x)
    def __repr__(self):
        return self.__x

def is_AxiomElement(x):
    return isinstance(x, AxiomElement)

# An instance
axiom = Axiom(script_subdirectory=None)

def reduce_load_Axiom():
    return axiom

import os
def axiom_console():
    os.system('axiom -nox')

def __doctest_cleanup():
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
