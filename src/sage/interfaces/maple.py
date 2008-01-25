r"""
Interface to Maple

AUTHOR:
    -- William Stein (2005): maple interface
    -- Gregg Musiker (2006-02-02): tutorial
    -- William Stein (2006-03-05): added tab completion, e.g., maple.[tab],
       and help, e.g, maple.sin?.

You must have the optional commercial Maple interpreter installed and
available as the command \code{maple} in your PATH in order to use
this interface.  You do not have to install any optional \sage packages.


    Type \code{maple.[tab]} for a list of all the functions available
    from your Maple install.  Type \code{maple.[tab]?} for Maple's
    help about a given function.  Type \code{maple(...)} to create
    a new Maple object, and \code{maple.eval(...)} to run a string
    using Maple (and get the result back as a string).

EXAMPLES:
    sage: maple('3 * 5')
    15
    sage: maple.eval('ifactor(2005)')
    '``(5)*``(401)'
    sage: maple.ifactor(2005)
    ``(5)*``(401)
    sage: maple.fsolve('x^2=cos(x)+4', 'x=0..5')
    1.914020619
    sage: maple.factor('x^5 - y^5')
    (x-y)*(x^4+x^3*y+x^2*y^2+x*y^3+y^4)

If the string "error" (case insensitive) occurs in the
output of anything from Maple, a RuntimeError exception is raised.

\subsection{Tutorial}

AUTHOR:
    -- Gregg Musiker (2006-02-02): initial version.


This tutorial is based on the Maple Tutorial for number theory
from  \url{http://www.math.mun.ca/~drideout/m3370/numtheory.html}.

There are several ways to use the Maple Interface in \SAGE.  We will
discuss two of those ways in this tutorial.

\begin{enumerate}
\item If you have a maple expression such as
\begin{verbatim}
factor( (x^5-1));
\end{verbatim}
We can write that in sage as

    sage: maple('factor(x^5-1)')
    (x-1)*(x^4+x^3+x^2+x+1)

Notice, there is no need to use a semicolon.

\item Since \SAGE is written in Python, we can also import maple
commands and write our scripts in a pythonic way.
For example, \code{factor()} is a maple command, so we can also factor
in \sage using

    sage: maple('(x^5-1)').factor()
    (x-1)*(x^4+x^3+x^2+x+1)

where \code{expression.command()} means the same thing as
\code{command(expression)} in Maple.  We will use this second type of
syntax whenever possible, resorting to the first when needed.

    sage: maple('(x^12-1)/(x-1)').simplify()
    x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1

\end{enumerate}

The normal command will always reduce a rational function to the
lowest terms. The factor command will factor a polynomial with
rational coefficients into irreducible factors over the ring of
integers. So for example,

    sage: maple('(x^12-1)').factor( )
    (x-1)*(x+1)*(x^2+x+1)*(x^2-x+1)*(x^2+1)*(x^4-x^2+1)

    sage: maple('(x^28-1)').factor( )
    (x-1)*(x^6+x^5+x^4+x^3+x^2+x+1)*(x+1)*(1-x+x^2-x^3+x^4-x^5+x^6)*(x^2+1)*(x^12-x^10+x^8-x^6+x^4-x^2+1)


Another important feature of maple is its online help.  We can access
this through sage as well.  After reading the description of the
command, you can press q to immediately get back to your original
prompt.

% NOTE: DOESN'T BRING UP NEW SCREEN IN SSH

Incidentally you can always get into a maple console by the command

    sage: maple.console()          # not tested
    sage: !maple                   # not tested

Note that the above two commands are slightly different, and the first
is preferred.

For example, for help on the maple command fibonacci, we type

    sage: maple.help('fibonacci')  # not tested, since it uses a pager

We see there are two choices.  Type

    sage: maple.help('combinat, fibonacci')   # not tested, since it uses a pager

We now see how the Maple command fibonacci works under the
combinatorics package.  Try typing in

    sage: maple.fibonacci(10)
    fibonacci(10)

You will get fibonacci(10) as output since Maple has not loaded the
combinatorics package yet.  To rectify this type

    sage: maple('combinat[fibonacci]')(10)
    55

instead.

If you want to load the combinatorics package for future calculations,
in \sage this can be done as

    sage: maple.with_package('combinat')

or

    sage: maple.load('combinat')

Now if we type \code{maple.fibonacci(10)}, we get the correct output:

    sage: maple.fibonacci(10)
    55

Some common maple packages include \code{combinat}, \code{linalg}, and
\code{numtheory}.  To produce the first 19 Fibonacci
numbers, use the sequence command.

    sage: maple('seq(fibonacci(i),i=1..19)')
    1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584,
    4181

Two other useful Maple commands are ifactor and isprime. For example

    sage: maple.isprime(maple.fibonacci(27))
    false
    sage: maple.ifactor(maple.fibonacci(27))
    ``(2)*``(17)*``(53)*``(109)

Note that the isprime function that is included with \sage (which uses
PARI) is better than the Maple one (it is faster and gives a provably
correct answer, whereas Maple is sometimes wrong).

    sage: alpha = maple('(1+sqrt(5))/2')
    sage: beta = maple('(1-sqrt(5))/2')
    sage: f19  = alpha^19 - beta^19/maple('sqrt(5)')
    sage: f19
    (1/2+1/2*5^(1/2))^19-1/5*(1/2-1/2*5^(1/2))^19*5^(1/2)
    sage: f19.simplify()                # somewhat randomly ordered output...
    6765+5778/5*5^(1/2)


Let's say we want to write a maple program now that squares a number
if it is positive and cubes it if it is negative.  In maple, that
would look like

\begin{verbatim}
mysqcu := proc(x)
if x > 0 then x^2;
else x^3; fi;
end;
\end{verbatim}
In SAGE, we write

   sage: mysqcu = maple('proc(x) if x > 0 then x^2 else x^3 fi end')
   sage: mysqcu(5)
   25
   sage: mysqcu(-5)
   -125

More complicated programs should be put in a separate file and
loaded.
"""

#############################################################################
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#############################################################################

from __future__ import with_statement

import os

from expect import Expect, ExpectElement, ExpectFunction, FunctionElement, gc_disabled

import pexpect

from sage.misc.misc import verbose, DOT_SAGE
from sage.misc.pager import pager

COMMANDS_CACHE = '%s/maple_commandlist_cache.sobj'%DOT_SAGE

class Maple(Expect):
    """
    Interface to the Maple interpreter.

    Type \code{maple.[tab]} for a list of all the functions available
    from your Maple install.  Type \code{maple.[tab]?} for Maple's
    help about a given function.  Type \code{maple(...)} to create
    a new Maple object, and \code{maple.eval(...)} to run a string
    using Maple (and get the result back as a string).
    """
    def __init__(self, maxread=100, script_subdirectory="", server=None, server_tmpdir=None, logfile=None):
        """
        Create an instance of the Maple interpreter.
        """
        Expect.__init__(self,
                        name = 'maple',
                        prompt = '#-->',
                        command = "maple -t",
                        maxread = maxread,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        server = server,
                        server_tmpdir = server_tmpdir,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=1)  # very important that this is 1
        # It's very important to use file i/o for everything,
        # since maple stupid command line interface always
        # dumps you into the editor when an error occurs,
        # and I can find no way to turn it off!!

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MapleFunction(self, attrname)

    def _keyboard_interrupt(self):
        print "Interrupting %s..."%self
        self._expect.sendline(chr(3))  # send ctrl-c
        self._expect.expect(self._prompt)
        self._expect.expect(self._prompt)
        raise RuntimeError, "Ctrl-c pressed while running %s"%self

    def __reduce__(self):
        return reduce_load_Maple, tuple([])

    def _read_in_file_command(self, filename):
        return 'read "%s"'%filename

    def _quit_string(self):
        return 'quit'

    def _install_hints(self):
        """
        Hints for installing Maple on your computer.

        AUTHOR:
            -- William Stein and Justin Walker (2006-02-12).
        """
        return """

In order to use the Maple interface you need to have Maple installed
and have a script in your PATH called "maple" that runs the
command-line version of Maple.  Alternatively, you could use a remote
connection to a server running Maple; for hints, type
    print maple._install_hints_ssh()

  (1) You might have to buy Maple (http://webstore.maplesoft.com/).

  (2) * LINUX: The maple script comes standard with your Maple install.

      * APPLE OS X:
          (a) create a file called maple (in your PATH), with the following contents:
             #!/bin/sh
             /Library/Frameworks/Maple.framework/Versions/Current/bin/maple $@
          (b) Save the file.
          (c) Make the file executable.
                chmod +x maple

      * WINDOWS:
        You must install Maple-for-Linux into the VMware machine (sorry, that's
        the only way at present).
"""

    def expect(self):
        return self._expect

    def console(self):
        maple_console()

##     def killall(self):
##         """
##         Kill all running instances of the maple interpreter
##         on this system.

##         TODO: When SAGE exists it doesn't correctly by default kill
##         all running Maple interpreters, for some strange reason.
##         Calling this function uses the kill and pidof operating system
##         programs to find all instances of cmaple and kill them.
##         """
##         import os
##         self._expect = None
##         while True:
##             pid = os.popen("pidof cmaple").read()[:-1]
##             if len(pid) > 0:
##                 os.system('kill -9 %s'%pid)
##             else:
##                 break

    def completions(self, s):
        """
        Return all commands that complete the command starting with the
        string s.   This is like typing s[Ctrl-T] in the maple interpreter.
        """
        bs = chr(8)*len(s)
        if self._expect is None:
            self._start()
        E = self._expect
        E.sendline('%s%s%s'%(s,chr(20),bs))
        t = E.timeout
        E.timeout=0.3  # since some things have no completion
        try:
            E.expect('----')
        except pexpect.TIMEOUT:
            E.timeout = t
            return []
        E.timeout = t
        v = E.before
        E.expect(self._prompt)
        E.expect(self._prompt)
        return v.split()[2:]

    def _commands(self):
        """
        Return list of all commands defined in maple.
        """
        try:
            v = sum([self.completions(chr(65+n)) for n in range(26)], []) + \
                sum([self.completions(chr(97+n)) for n in range(26)], [])
        except RuntimeError:
            print "\n"*3
            print "*"*70
            print "WARNING: You do not have a working version of Maple installed!"
            print "*"*70
            v = []
        v.sort()
        return v

    def trait_names(self, verbose=True, use_disk_cache=True):
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
                print "\nBuilding Maple command completion list (this takes"
                print "a few seconds only the first time you do it)."
                print "To force rebuild later, delete %s."%COMMANDS_CACHE
            v = self._commands()
            self.__trait_names = v
            if len(v) > 200:
                # Maple is actually installed.
                sage.misc.persist.save(v, COMMANDS_CACHE)
            return v

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True):
        line += ';'
        with gc_disabled():
            z = Expect._eval_line(self, line, allow_use_file=allow_use_file,
                    wait_for_prompt=wait_for_prompt).replace('\\\n','').strip()
            if z.lower().find("error") != -1:
                # The following was very tricky to figure out.
                # When an error occurs using Maple, unfortunately,
                # Maple also dumps one into the line where the
                # error occured with that line copied in.  This
                # totally messes up the pexpect interface.  However,
                # I think the following few lines successfully
                # "clear things out", i.e., delete the text from
                # the edit buffer and get a clean prompt.
                e = self.expect()
                e.sendline('%s__sage__;'%(chr(8)*len(line)))
                e.expect('__sage__;')
                e.expect(self._prompt)
                raise RuntimeError, "An error occured running a Maple command:\nINPUT:\n%s\nOUTPUT:\n%s"%(line, z)
        return z

    def cputime(self, t=None):
        if t is None:
            return float(self('time()'))
        else:
            return float(self('time() - %s'%float(t)))

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s:=%s:'%(var,value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError, "Error executing code in Maple\nCODE:\n\t%s\nMaple ERROR:\n\t%s"%(cmd, out)

    def get(self, var):
        """
        Get the value of the variable var.
        """
        s = self.eval('printf("%%q",%s)'%var)
        return s

    def get_using_file(self, var):
        """
        Get the value of the variable var using a file.

        (I would make this the default for values that are bigger than
        a few thousand characters.  However, it's not at all obvious
        how to figure out if the string representation of an object is
        big ahead of time!  We assume it is for now, if the string
        used to create the object was big.)
        """
        # pdehaye/ 20070819: This might be obsolete
        s = self.eval('save %s, "%s"'%(var, self._remote_tmpfile()))
        s = open(self._remote_tmpfile()).read().replace('\\\n','')
        i = s.find('=')
        return s[i+2:-2]

    def _object_class(self):
        return MapleElement

    def _equality_symbol(self):
        return '='

    def _true_symbol(self):
        return 'true'

    def _assign_symbol(self):
        return ":="

    def _source(self, s):
        """
        Trys to return the source code of a Maple function str
        as a string.

        EXAMPLES:
            sage: print maple._source('curry').strip()  #optional requires maple
            p -> subs('_X' = args[2 .. nargs], () -> p(_X, args))
            sage: maple._source('ZZZ')                  #optional requires maple
            Traceback (most recent call last):
            ...
            Exception: no source code could be found
        """
        cmd = 'echo "interface(verboseproc=2): print(%s);" | maple -q'%s
        src = os.popen(cmd).read()
        if src.strip() == s:
            raise Exception, "no source code could be found"
        else:
            return src

    def source(self, s):
        """
        Display the Maple source (if possible) about s.  This is the same as
        returning the output produced by the following Maple commands:

        interface(verboseproc=2): print(s)

        INPUT:
            s -- a string representing the function whose source code you
                 want
        """
        try:
            pager()(self._source(s))
        except Exception:
            pager()('No source code could be found.')

    def _help(self, str):
        return os.popen('echo "?%s" | maple -q'%str).read()

    def help(self, str):
        """
        Display Maple help about str.  This is the same as typing "?str" in
        the Maple console.

        INPUT:
            str -- a string to search for in the maple help system
        """
        pager()(self._help(str))

    def with_package(self, package):
        """
        Make a package of Maple procedures available in the
        interpreter.

        INPUT:
            package -- string

        EXAMPLES:
        Some functions are unknown to Maple until you use with to include
        the appropriate package.

            sage: maple.quit()   # optional -- to reset maple.
            sage: maple('partition(10)')              # optional
            partition(10)
            sage: maple('bell(10)')                   # optional
            bell(10)
            sage: maple.with_package('combinat')               # optional
            sage: maple('partition(10)')               # optional
            [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 2], [1, 1, 1, 1, 1, 1, 2, 2], [1, 1, 1, 1, 2, 2, 2], [1, 1, 2, 2, 2, 2], [2, 2, 2, 2, 2], [1, 1, 1, 1, 1, 1, 1, 3], [1, 1, 1, 1, 1, 2, 3], [1, 1, 1, 2, 2, 3], [1, 2, 2, 2, 3], [1, 1, 1, 1, 3, 3], [1, 1, 2, 3, 3], [2, 2, 3, 3], [1, 3, 3, 3], [1, 1, 1, 1, 1, 1, 4], [1, 1, 1, 1, 2, 4], [1, 1, 2, 2, 4], [2, 2, 2, 4], [1, 1, 1, 3, 4], [1, 2, 3, 4], [3, 3, 4], [1, 1, 4, 4], [2, 4, 4], [1, 1, 1, 1, 1, 5], [1, 1, 1, 2, 5], [1, 2, 2, 5], [1, 1, 3, 5], [2, 3, 5], [1, 4, 5], [5, 5], [1, 1, 1, 1, 6], [1, 1, 2, 6], [2, 2, 6], [1, 3, 6], [4, 6], [1, 1, 1, 7], [1, 2, 7], [3, 7], [1, 1, 8], [2, 8], [1, 9], [10]]
            sage: maple('bell(10)')                   # optional
            115975
            sage: maple('fibonacci(10)')              # optional
            55
        """
        self.eval('with(%s)'%package)

    load = with_package

    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
        # Unfortunately MAPLE does not have a clear command.
        # The next best thing is to set equal to the constant
        # 0, so that memory will be freed.
    #    self.eval("%s=0;"%var)

class MapleFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M._help(self._name)

    def _sage_src_(self):
        """
        Returns the source code of self.  This is the function that eventually
        gets called when doing maple.gcd?? for example.

        EXAMPLES:
            sage: print maple.curry._sage_src_().strip() #optional requires maple
            p -> subs('_X' = args[2 .. nargs], () -> p(_X, args))
            sage: maple.ZZZ._sage_src_()                 #optional requires maple
            Traceback (most recent call last):
            ...
            Exception: no source code could be found

        """
        M = self._parent
        return M._source(self._name)

class MapleFunctionElement(FunctionElement):
    def _sage_doc_(self):
        return self._obj.parent()._help(self._name)
    def _sage_src_(self):
        """
        Returns the source code of self.

        EXAMPLES:
            sage: g = maple('gcd')                   #optional requires maple
            sage: print g.curry._sage_src_().strip() #optional requires maple
            p -> subs('_X' = args[2 .. nargs], () -> p(_X, args))
            sage: m = maple('2')                     #optional requires maple
            sage: m.ZZZ._sage_src_()                 #optional requires maple
            Traceback (most recent call last):
            ...
            Exception: no source code could be found

        """
        return self._obj.parent()._source(self._name)

class MapleElement(ExpectElement):
    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MapleFunctionElement(self, attrname)

    def __float__(self):
        M = self.parent()
        return float(maple.eval('evalf(%s)'%self.name()))

    def __cmp__(self, other):
        """
        Compare equality between self and other, using maple.

        EXAMPLES:
            sage: a = maple(5)
            sage: b = maple(5)
            sage: a == b
            True
            sage: a == 5
            True

            sage: c = maple(3)
            sage: a == c
            False
            sage: a < c
            False
            sage: a < 6
            True
            sage: c <= a
            True

            sage: M = matrix(ZZ, 2, range(1,5))
            sage: Mm = maple(M)
            sage: Mm == Mm
            True
            sage: Mm < 5
            True
            sage: (Mm < 5) == (M < 5)
            True
            sage: 5 < Mm
            False

        TESTS:
            sage: x = var('x')
            sage: t = maple((x+1)^2)
            sage: u = maple(x^2+2*x+1)
            sage: u == t # todo: not implemented
            True         # returns False, should use 'testeq' in maple
            sage: maple.eval('testeq(%s = %s)'%(t.name(),u.name()))
            'true'
        """
        P = self.parent()
        if P.eval("evalb(%s %s %s)"%(self.name(), P._equality_symbol(),
                                 other.name())) == P._true_symbol():
            return 0
        # Maple does not allow comparing objects of different types and
        # it raises an error in this case.
        # We catch the error, and return True for <
        try:
            if P.eval("evalb(%s %s %s)"%(self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
                return -1
        except RuntimeError, e:
            msg = str(e)
            if 'is not valid' in msg and 'to < or <=' in msg:
                if (hash(str(self)) < hash(str(other))):
                    return -1
                else:
                    return 1
            else:
                raise RuntimeError, e
        if P.eval("evalb(%s %s %s)"%(self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
            return 1
        # everything is supposed to be comparable in Python, so we define
        # the comparison thus when no comparable in interfaced system.
        if (hash(str(self)) < hash(str(other))):
            return -1
        else:
            return 1

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: t = maple(5); u = maple(3)
            sage: t*u
            15
            sage: M = matrix(ZZ,2,range(4))
            sage: Mm = maple(M)
            sage: Mm*Mm
            Matrix(2, 2, [[2,3],[6,11]])

            sage: v = vector(ZZ,2,[2,3])
            sage: vm = maple(v)
            sage: vm*Mm
            Vector[row](2, [6,11])

            sage: t*Mm
            Matrix(2, 2, [[0,5],[10,15]])
        """
        P = self._check_valid()
        try:
            return P.new('%s . %s'%(self._name, right._name))
        except Exception, msg:
            raise TypeError,msg

    def trait_names(self):
        return self.parent().trait_names()

    def __repr__(self):
        """
        Return a string representation of self.

        EXAMPLES:
            sage: x = var('x')
            sage: maple(x)
            x
            sage: maple(5)
            5
            sage: M = matrix(QQ,2,range(4))
            sage: maple(M)
            Matrix(2, 2, [[0,1],[2,3]])
        """
        self._check_valid()
        return self.parent().get(self._name)

    def _latex_(self):
        r"""
        You can output Maple expressions in latex.

        EXAMPLES:
            sage: print latex(maple('(x^4 - y)/(y^2-3*x)'))      # optional
            {\frac {{x}^{4}-y}{{y}^{2}-3\,x}}
            sage: print latex(maple(pi - e^3))                   # optional
            \pi - \left( {e^{1}} \right) ^{3}

        \note{Some expressions might require the Maple style file
        \code{maple2e.sty} in order to latex correctly.}
        """
        return self.parent().eval('latex(%s)'%self.name())

# An instance
maple = Maple(script_subdirectory='user')

def reduce_load_Maple():
    return maple


import os
def maple_console():
    os.system('maple')


def __doctest_cleanup():
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()



"""
The following only works in Maple >= 9, I guess, but could
be useful.

From Jaap Spies:  In addition Maple has a nice feature the function

 > FunctionAdvisor();

 > FunctionAdvisor(topics, quiet);
      [DE, analytic_extension, asymptotic_expansion, branch_cuts,
      branch_points, calling_sequence, class_members,
      classify_function, definition, describe, differentiation_rule,
      function_classes, identities, integral_form,
      known_functions, relate, series, singularities, special_values,
      specialize, sum_form, synonyms]

 > FunctionAdvisor(syntax, hypergeom);
                                            hypergeom([a, b], [c], z)

Eventually this could be used to do an intelligent command
completion.
"""
