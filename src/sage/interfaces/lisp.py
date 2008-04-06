r"""nodoctest
Lisp Interface

EXAMPLES:
    sage: lisp.eval('(* 4 5)')
    '20'
    sage: a = lisp(3); b = lisp(5)
    sage: a + b
    8
    sage: a * b
    15
    sage: a / b
    3/5
    sage: a - b
    -2
    sage: a.sin()
    0.14112
    sage: b.cos()
    0.2836622
    sage: a.exp()
    20.085537
    sage: lisp.eval('(+ %s %s)'%(a.name(), b.name()))
    '8'

One can define functions and the interface supports object-oriented
notation for calling them:
    sage: lisp.eval('(defun factorial (n) (if (= n 1) 1 (* n (factorial (- n 1)))))')
    'FACTORIAL'
    sage: lisp('(factorial 10)')
    3628800
    sage: lisp(10).factorial()
    3628800
    sage: a = lisp(17)
    sage: a.factorial()
    355687428096000

AUTHORS:
    -- William Stein (first version)
    -- William Stein (2007-06-20): significant improvements.
"""

##########################################################################
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#
##########################################################################

from __future__ import with_statement

import random

from expect import Expect, ExpectElement, ExpectFunction, FunctionElement, gc_disabled
from sage.misc.misc import verbose, UNAME, is_64bit
from sage.structure.element import RingElement

is_64bit_linux = is_64bit() and UNAME == "Linux"

class Lisp(Expect):
    def __init__(self,
                 maxread=100000, script_subdirectory=None,
                 logfile=None,
                 server=None,
                 server_tmpdir=None):
        Expect.__init__(self,

                        # The capitalized version of this is used for printing.
                        name = 'lisp',

                        # This is regexp of the input prompt.  If you can change
                        # it to be very obfuscated that would be better.   Even
                        # better is to use sequence numbers.
                        prompt = '\[[0-9]+\]> ',

                        # This is the command that starts up your program
                        command = "clisp --silent -on-error abort",

                        maxread = maxread,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory = script_subdirectory,

                        # If this is true, then whenever the user presses Control-C to
                        # interrupt a calculation, the whole interface is restarted.
                        restart_on_ctrlc = False,

                        # If true, print out a message when starting
                        # up the command when you first send a command
                        # to this interface.
                        verbose_start = False,

                        logfile=logfile,

                        # If an input is longer than this number of characters, then
                        # try to switch to outputing to a file.
                        eval_using_file_cutoff=1024)

        self.__seq = 0
        self.__in_seq = 1

    def eval(self, code, strip=True):
        with gc_disabled():
            self._synchronize()
            code = str(code)
            code = code.strip()
            code = code.replace('\n',' ')
            x = []
            for L in code.split('\n'):
                if L != '':
                    try:
                        s = self.__in_seq + 1
                        pr = '\[%s\]>'%s
                        M = self._eval_line(L, wait_for_prompt=self._prompt)
                        if is_64bit_linux:
                            phrase = '[C\x1b[C\n'
                        else:
                            phrase = L
                        i = M.rfind(phrase)
                        if i > 1:
                            M = M[i+len(phrase):]
                        x.append(M.strip())
                        self.__in_seq = s
                    except KeyboardInterrupt:
                        # DO NOT CATCH KeyboardInterrupt, as it is being caught
                        # by _eval_line
                        # In particular, do NOT call self._keyboard_interrupt()
                        raise
                    except TypeError, s:
                        return 'error evaluating "%s":\n%s'%(code,s)
            return '\n'.join(x)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '(setq %s %s)'%(var, value)
        out = self.eval(cmd)
        if '***' in out:
            raise TypeError, "Error executing code in SAGE\nCODE:\n\t%s\nSAGE ERROR:\n\t%s"%(cmd, out)

    def get(self, var):
        out = self.eval(var).lstrip().lstrip(var).lstrip()
        return out

    def _start(self, *args, **kwds):
        Expect._start(self, *args, **kwds)
        self.__in_seq = 1

    def _synchronize(self):
        E = self._expect
        if E is None:
            self._start()
            E = self._expect
        r = random.randrange(2147483647)
        s = str(r+1)
        cmd = "(+ 1 %s)\n"%r
        E.sendline(cmd)
        E.expect(s)
        E.expect(self._prompt)

    def _repr_(self):
        return 'Lisp Interpreter'

    def __reduce__(self):
        return reduce_load_lisp, tuple([])

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return LispFunction(self, attrname)

    def _quit_string(self):
        raise NotImplementedError

    def _read_in_file_command(self, filename):
        raise NotImplementedError

    def trait_names(self):
        ## [[implement giving a list of all functions and identifiers in the system]]
        raise NotImplementedError

    def read(self, filename):
        # [[implement loading of the contents of filename into the system]]
        raise NotImplementedError


    def kill(self, var):
        # [[send code that kills the variable with given name in the system.]]
        pass

    def console(self):
        # run the console command (defined below).
        pass

    def version(self):
        # run the version command (defined below)
        pass

    def _object_class(self):
        return LispElement

    def _true_symbol(self):
        # return the string rep of truth, i.e., what the system outputs
        # when you type 1==1.
        raise NotImplementedError

    def _false_symbol(self):
        # return the string rep of truth, i.e., what the system outputs
        # when you type 1==2.
        raise NotImplementedError

    def _equality_symbol(self):
        # return the symbol for checking equality, e.g., == or eq.
        raise NotImplementedError

    def help(self, command):
        # return help on a given command.
        raise NotImplementedError

    def function_call(self, function, args=[]):
        if function == '':
            raise ValueError, "function name must be nonempty"
        if function[:2] == "__":
            raise AttributeError
        if not isinstance(args, list):
            args = [args]
        for i in range(len(args)):
            if not isinstance(args[i], ExpectElement):
                args[i] = self.new(args[i])
        return self.new("(%s %s)"%(function, ",".join([s.name() for s in args])))

class LispElement(ExpectElement):
    """
    Describe elements of your system here.
    """
    def _add_(self, right):
        P = self._check_valid()
        return P.new('(+ %s %s)'%(self._name, right._name))

    def _sub_(self, right):
        P = self._check_valid()
        return P.new('(- %s %s)'%(self._name, right._name))

    def _mul_(self, right):
        P = self._check_valid()
        return P.new('(* %s %s)'%(self._name, right._name))

    def _div_(self, right):
        P = self._check_valid()
        return P.new('(/ %s %s)'%(self._name, right._name))

    def __pow__(self, n):
        return RingElement.__pow__(self, n)

class LispFunctionElement(FunctionElement):
    def _sage_doc_(self):
        M = self._obj.parent()
        return M.help(self._name)


class LispFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)



def is_LispElement(x):
    return isinstance(x, LispElement)

# An instance
lisp = Lisp()

def reduce_load_Lisp():
    return lisp

import os
def lisp_console():
    os.system('clisp')


