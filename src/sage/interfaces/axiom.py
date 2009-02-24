r"""
Interface to Axiom

TODO:

- Evaluation using a file is not done. Any input line with more than a
  few thousand characters would hang the system, so currently it
  automatically raises an exception.

- All completions of a given command.

- Interactive help.

Axiom is a free GPL-compatible (modified BSD license) general
purpose computer algebra system whose development started in 1973
at IBM. It contains symbolic manipulation algorithms, as well as
implementations of special functions, including elliptic functions
and generalized hypergeometric functions. Moreover, Axiom has
implementations of many functions relating to the invariant theory
of the symmetric group `S_n`. For many links to Axiom
documentation see http://wiki.axiom-developer.org.

AUTHORS:

- Bill Page (2006-10): Created this (based on maxima interface)


  .. note::

     Bill Page put a huge amount of effort into the Sage Axiom
     interface over several days during the Sage Days 2 coding
     sprint. This is contribution is greatly appreciated.

- William Stein (2006-10): misc touchup.

- Bill Page (2007-08): Minor modifications to support axiom4sage-0.3

.. note::

   The axiom4sage-0.3.spkg is based on an experimental version of the
   FriCAS fork of the Axiom project by Waldek Hebisch that uses
   pre-compiled cached Lisp code to build Axiom very quickly with
   clisp.

If the string "error" (case insensitive) occurs in the output of
anything from axiom, a RuntimeError exception is raised.

EXAMPLES: We evaluate a very simple expression in axiom.

::

    sage: axiom('3 * 5')                     # optional
    15
    sage: a = axiom(3) * axiom(5); a         # optional
    15

The type of a is AxiomElement, i.e., an element of the axiom
interpreter.

::

    sage: type(a)                            # optional
    <class 'sage.interfaces.axiom.AxiomElement'>
    sage: parent(a)                          # optional
    Axiom

The underlying Axiom type of a is also available, via the type
method::

    sage: a.type()                           # optional
    PositiveInteger

We factor `x^5 - y^5` in Axiom in several different ways.
The first way yields a Axiom object.

::

    sage: F = axiom.factor('x^5 - y^5'); F      # optional
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    sage: type(F)                               # optional
    <class 'sage.interfaces.axiom.AxiomElement'>
    sage: F.type()                              # optional
    Factored Polynomial Integer

Note that Axiom objects are normally displayed using "ASCII art".

::

    sage: a = axiom(2/3); a          # optional
      2
      -
      3
    sage: a = axiom('x^2 + 3/7'); a      # optional
       2   3
      x  + -
           7

The ``axiom.eval`` command evaluates an expression in
axiom and returns the result as a string. This is exact as if we
typed in the given line of code to axiom; the return value is what
Axiom would print out.

::

    sage: print axiom.eval('factor(x^5 - y^5)')   # optional
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

    Type: Factored Polynomial Integer

We can create the polynomial `f` as a Axiom polynomial,
then call the factor method on it. Notice that the notation
``f.factor()`` is consistent with how the rest of Sage
works.

::

    sage: f = axiom('x^5 - y^5')                  # optional
    sage: f^2                                     # optional
       10     5 5    10
      y   - 2x y  + x
    sage: f.factor()                              # optional
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

Control-C interruption works well with the axiom interface, because
of the excellent implementation of axiom. For example, try the
following sum but with a much bigger range, and hit control-C.

::

    sage:  f = axiom('(x^5 - y^5)^10000')       # not tested
    Interrupting Axiom...
    ...
    <type 'exceptions.TypeError'>: Ctrl-c pressed while running Axiom

::

    sage: axiom('1/100 + 1/101')                  # optional
       201
      -----
      10100
    sage: a = axiom('(1 + sqrt(2))^5'); a         # optional
         +-+
      29\|2  + 41

TESTS: We check to make sure the subst method works with keyword
arguments.

::

    sage: a = axiom(x+2); a  #optional
    x + 2
    sage: a.subst(x=3)       #optional
    5

We verify that Axiom floating point numbers can be converted to
Python floats.

::

    sage: float(axiom(2))     #optional
    2.0
"""

###########################################################################
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#                     2007 Bill Page
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
###########################################################################

import os, re

from expect import Expect, ExpectElement, FunctionElement, ExpectFunction
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
    def __init__(self, script_subdirectory=None, logfile=None, server=None, server_tmpdir=None):
        """
        Create an instance of the Axiom interpreter.

        TESTS::

            sage: axiom == loads(dumps(axiom))
            True
        """
        eval_using_file_cutoff = 200
        self.__eval_using_file_cutoff = eval_using_file_cutoff
        Expect.__init__(self,
                        name = 'axiom',
                        prompt = '\([0-9]+\) -> ',
                        command = "axiom -nox -noclef",
                        maxread = 10,
                        script_subdirectory = script_subdirectory,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        #init_code = [')lisp (si::readline-off)'],
                        logfile = logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)

    def _function_class(self):
        """
        Return the AxiomExpectFunction class.

        EXAMPLES::

            sage: axiom._function_class()
            <class 'sage.interfaces.axiom.AxiomExpectFunction'>
            sage: type(axiom.gcd)
            <class 'sage.interfaces.axiom.AxiomExpectFunction'>
        """
        return AxiomExpectFunction

    def _object_class(self):
        """
        EXAMPLES::

            sage: axiom._object_class()
            <class 'sage.interfaces.axiom.AxiomElement'>
            sage: type(axiom(2)) #optional -- requires Axiom
            <class 'sage.interfaces.axiom.AxiomElement'>
        """
        return AxiomElement

    def _function_element_class(self):
        """
        Returns the Axiom function element class.

        EXAMPLES::

            sage: axiom._function_element_class()
            <class 'sage.interfaces.axiom.AxiomFunctionElement'>
            sage: type(axiom(2).gcd) #optional -- requires Axiom
            <class 'sage.interfaces.axiom.AxiomFunctionElement'>
        """
        return AxiomFunctionElement

    def _start(self):
        """
        Start the Axiom interpreter.

        EXAMPLES::

            sage: a = Axiom()
            sage: a.is_running()
            False
            sage: a._start()     #optional -- requires axiom
            sage: a.is_running() #optional
            True
            sage: a.quit()       #optional
        """
        Expect._start(self)
        out = self._eval_line(')set functions compile on', reformat=False)
        out = self._eval_line(')set output length 245', reformat=False)
        out = self._eval_line(')set message autoload off', reformat=False)

    def _read_in_file_command(self, filename):
        r"""
        EXAMPLES::

            sage: axiom._read_in_file_command('test.input')
            ')read test.input \n'
            sage: axiom._read_in_file_command('test')
            Traceback (most recent call last):
            ...
            ValueError: the filename must end with .input

        ::

            sage: filename = tmp_filename()+'.input'
            sage: f = open(filename, 'w')
            sage: f.write('xx := 22;\n')
            sage: f.close()
            sage: axiom.read(filename)    #optional -- requires Axiom
            sage: axiom.get('xx')         #optional
            '22'
        """
        if not filename.endswith('.input'):
            raise ValueError, "the filename must end with .input"

        # For some reason this trivial comp
        # keeps certain random freezes from occuring.  Do not remove this.
        # The space before the \n is also important.
        return ')read %s \n'%filename

    def __reduce__(self):
        """
        EXAMPLES::

            sage: axiom.__reduce__()
            (<function reduce_load_Axiom at 0x...>, ())
            sage: f, args = _
            sage: f(*args)
            Axiom
        """
        return reduce_load_Axiom, tuple([])

    def _quit_string(self):
        """
        Returns the string used to quit Axiom.

        EXAMPLES::

            sage: axiom._quit_string()
            ')lisp (quit)'

        ::

            sage: a = Axiom()
            sage: a.is_running()
            False
            sage: a._start()     #optional -- requires axiom
            sage: a.is_running() #optional
            True
            sage: a.quit()       #optional
            sage: a.is_running() #optional
            False
        """
        return ')lisp (quit)'

    def _eval_line(self, line, reformat=True, allow_use_file=False,
                   wait_for_prompt=True):
        """
        EXAMPLES::

            sage: print axiom._eval_line('2+2')  #optional -- requires Axiom
              4
                                                       Type: PositiveInteger
        """
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
            return self._eval_line_using_file(line)
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

    def _commands(self):
        """
        Returns a list of commands available. This is done by parsing the
        result of the first section of the output of ')what things'.

        EXAMPLES::

            sage: cmds = axiom._commands() #optional -- requires Axiom
            sage: len(cmds) > 100  #optional
            True
            sage: '<' in cmds      #optional
            True
            sage: 'factor' in cmds #optional
            True
        """
        s = self.eval(")what things")
        start = '\r\n\r\n#'
        i = s.find(start)
        end = "To get more information about"
        j = s.find(end)
        s = s[i+len(start):j].split()
        return s


    def trait_names(self, verbose=True, use_disk_cache=True):
        """
        Returns a list of all the commands defined in Axiom and optionally
        (per default) store them to disk.

        EXAMPLES::

            sage: c = axiom.trait_names(use_disk_cache=False, verbose=False) #optional
            sage: len(c) > 100  #optional
            True
            sage: 'factor' in c  #optional
            True
            sage: '**' in c     #optional
            False
            sage: 'upperCase?' in c  #optional
            False
            sage: 'upperCase_q' in c #optional
            True
            sage: 'upperCase_e' in c #optional
            True
        """
        try:
            return self.__trait_names
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__trait_names = sage.misc.persist.load(COMMANDS_CACHE)
                    return self.__trait_names
                except IOError:
                    pass
            if verbose:
                print "\nBuilding Axiom command completion list (this takes"
                print "a few seconds only the first time you do it)."
                print "To force rebuild later, delete %s."%COMMANDS_CACHE
            v = self._commands()

            #Process we now need process the commands to strip out things which
            #are not valid Python identifiers.
            import re
            valid = re.compile('[^a-zA-Z0-9_]+')
            names = [x for x in v if valid.search(x) is None]

            #Change everything that ends with ? to _q and
            #everything that ends with ! to _e
            names += [x[:-1]+"_q" for x in v if x.endswith("?")]
            names += [x[:-1]+"_e" for x in v if x.endswith("!")]

            self.__trait_names = names
            if len(v) > 200:
                # Axiom is actually installed.
                sage.misc.persist.save(v, COMMANDS_CACHE)
            return names

    def set(self, var, value):
        """
        Set the variable var to the given value.

        EXAMPLES::

            sage: axiom.set('xx', '2')    #optional -- requires Axiom
            sage: axiom.get('xx')         #optional
            '2'
        """
        cmd = '%s := %s'%(var, value)
        out = self._eval_line(cmd, reformat=False)

        if out.find("error") != -1:
            raise TypeError, "Error executing code in Axiom\nCODE:\n\t%s\nAxiom ERROR:\n\t%s"%(cmd, out)


    def get(self, var):
        r"""
        Get the string value of the Axiom variable var.

        EXAMPLES::

            sage: axiom.set('xx', '2')    #optional -- requires Axiom
            sage: axiom.get('xx')         #optional
            '2'
            sage: a = axiom('(1 + sqrt(2))^5') #optional
            sage: axiom.get(a.name())          #optional
            '     +-+\r\n  29\\|2  + 41'
        """
        s = self._eval_line(str(var))
        i = s.rfind('Type:')
        s = s[:i].rstrip().lstrip("\n")
        if '\n' not in s:
            s = s.strip()
        return s

    def console(self):
        """
        Spawn a new Axiom (FriCAS) command-line session.

        EXAMPLES::

            sage: axiom.console() #not tested
                             FriCAS (AXIOM fork) Computer Algebra System
                                     Version: FriCAS 2007-07-19
                          Timestamp: Saturday October 20, 2007 at 20:08:37
            -----------------------------------------------------------------------------
               Issue )copyright to view copyright notices.
               Issue )summary for a summary of useful system commands.
               Issue )quit to leave AXIOM and return to shell.
            -----------------------------------------------------------------------------
        """
        axiom_console()


class AxiomElement(ExpectElement):
    def __call__(self, x):
        """
        EXAMPLES::

            sage: f = axiom(x+2) #optional -- requires Axiom
            sage: f(2)           #optional
            4
        """
        self._check_valid()
        P = self.parent()
        return P('%s(%s)'%(self.name(), x))

    def __cmp__(self, other):
        """
        EXAMPLES::

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

        We can also compare more complicated object such as functions::

            sage: f = axiom('sin(x)'); g = axiom('cos(x)')    # optional
            sage: f == g                                      # optional
            False
        """

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: two = axiom(2)  #optional -- requires Axiom
            sage: two == 2        #optional
            True
            sage: two == 3        #optional
            False
            sage: two < 3         #optional
            True
            sage: two > 1         #optional
            True
        """
        P = self.parent()
        if 'true' in P.eval("(%s = %s) :: Boolean"%(self.name(),other.name())):
            return 0
        elif 'true' in P.eval("(%s < %s) :: Boolean"%(self.name(), other.name())):
            return -1
        elif 'true' in P.eval("(%s > %s) :: Boolean"%(self.name(),other.name())):
            return 1

        # everything is supposed to be comparable in Python, so we define
        # the comparison thus when no comparable in interfaced system.
        if (hash(self) < hash(other)):
            return -1
        else:
            return 1

    def type(self):
        """
        Returns the type of an AxiomElement.

        EXAMPLES::

            sage: axiom(x+2).type()  #optional -- requires Axiom
            Polynomial Integer
        """
        P = self._check_valid()
        s = P._eval_line(self.name())
        i = s.rfind('Type:')
        return P(s[i+5:].strip())

    def __len__(self):
        """
        Return the length of a list.

        EXAMPLES::

            sage: v = axiom('[x^i for i in 0..5]')            # optional
            sage: len(v)                                      # optional
            6
        """
        P = self._check_valid()
        s = P.eval('# %s '%self.name())
        i = s.rfind('Type')
        return int(s[:i-1])

    def __getitem__(self, n):
        r"""
        Return the n-th element of this list.

        .. note::

           Lists are 1-based.

        EXAMPLES::

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

    def comma(self, *args):
        """
        Returns a Axiom tuple from self and args.

        EXAMPLES::

            sage: two = axiom(2)  #optional -- requires Axiom
            sage: two.comma(3)    #optional
            [2,3]
            sage: two.comma(3,4)  #optional
            [2,3,4]
            sage: _.type()        #optional
            Tuple PositiveInteger
        """
        P = self._check_valid()
        args = list(args)
        for i, arg in enumerate(args):
            if not isinstance(arg, AxiomElement) or arg.parent() is not P:
                args[i] = P(arg)
        cmd = "(" + ",".join([x.name() for x in [self]+args]) + ")"
        return P(cmd)

    def _latex_(self):
        """
        EXAMPLES::

            sage: a = axiom(1/2) #optional -- requires Axiom
            sage: latex(a)       #optional
            1 \over 2
        """
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
    def __init__(self, object, name):
        """
        TESTS::

            sage: a = axiom('"Hello"') #optional -- requires Axiom
            sage: a.upperCase_q        #optional
            upperCase?
            sage: a.upperCase_e        #optional
            upperCase!
            sage: a.upperCase_e()      #optional
            "HELLO"
        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        FunctionElement.__init__(self, object, name)


class AxiomExpectFunction(ExpectFunction):
    def __init__(self, parent, name):
        """
        TESTS::

            sage: axiom.upperCase_q
            upperCase?
            sage: axiom.upperCase_e
            upperCase!
        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        ExpectFunction.__init__(self, parent, name)


def is_AxiomElement(x):
    """
    Returns True of x is of type AxiomElement.

    EXAMPLES::

        sage: from sage.interfaces.axiom import is_AxiomElement
        sage: is_AxiomElement(axiom(2)) #optional -- requires Axiom
        True
        sage: is_AxiomElement(2)
        False
    """
    return isinstance(x, AxiomElement)

# An instance
axiom = Axiom(script_subdirectory=None)

def reduce_load_Axiom():
    """
    Returns the Axiom interface object defined in
    sage.interfaces.axiom.

    EXAMPLES::

        sage: from sage.interfaces.axiom import reduce_load_Axiom
        sage: reduce_load_Axiom()
        Axiom
    """
    return axiom

import os
def axiom_console():
    """
    Spawn a new Axiom (FriCAS) command-line session.

    EXAMPLES::

        sage: axiom_console() #not tested
                         FriCAS (AXIOM fork) Computer Algebra System
                                 Version: FriCAS 2007-07-19
                      Timestamp: Saturday October 20, 2007 at 20:08:37
        -----------------------------------------------------------------------------
           Issue )copyright to view copyright notices.
           Issue )summary for a summary of useful system commands.
           Issue )quit to leave AXIOM and return to shell.
        -----------------------------------------------------------------------------
    """
    os.system('axiom -nox')

def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.axiom import __doctest_cleanup
        sage: a = Axiom()
        sage: two = a(2)     #optional -- requires Axiom
        sage: a.is_running() #optional
        True
        sage: __doctest_cleanup()
        sage: a.is_running()
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
