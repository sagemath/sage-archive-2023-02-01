r"""
Interface to MuPAD

AUTHOR:
    -- William Stein

You must have the optional commercial MuPAD interpreter installed and
available as the command \code{mupkern} in your PATH in order to use
this interface.  You do not have to install any optional \sage
packages.
"""

#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#############################################################################

seq = 0
PROMPT = '___SAGE___'

import os

from expect import (Expect, ExpectElement, ExpectFunction,
                    FunctionElement, AsciiArtString)

import pexpect

from sage.misc.misc import verbose, DOT_SAGE
from sage.misc.pager import pager

class Mupad(Expect):
    """
    Interface to the MuPAD interpreter.
    """
    def __init__(self, maxread=1000, script_subdirectory="", logfile=None):
        """
        Create an instance of the MuPAD interpreter.
        """
        Expect.__init__(self,
                        name = 'MuPAD',
                        prompt = '>>',
                        command = "mupkern -P e",
                        maxread = maxread,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        logfile = logfile)

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MupadFunction(self, attrname)

    def _keyboard_interrupt(self):
        print "Interrupting %s..."%self
        self._expect.sendline(chr(3))  # send ctrl-c
        self._expect.expect(PROMPT)
        self._expect.expect(PROMPT)
        raise RuntimeError, "Ctrl-c pressed while running %s"%self

    def __reduce__(self):
        return reduce_load_mupad, tuple([])

    def _read_in_file_command(self, filename):
        return 'read "%s"'%filename

    def _start(self, alt_message=None, block_during_init=True):
        Expect._start(self, alt_message=alt_message, block_during_init=block_during_init)
        self._expect.sendline('Pref::promptString("%s");'%PROMPT)
        self._expect.expect(PROMPT)
        self._expect.send('1;')
        self._expect.expect(PROMPT)

    def _quit_string(self):
        return 'quit'

    def _install_hints(self):
        """
        Hints for installing MuPAD on your computer.
        """
        return """
In order to use the MuPAD interface you need to have MuPAD installed
and have a script in your PATH called "mupkern" that runs the
command-line version of MuPAD.

  (1) You might have to buy MuPAD.

  (2) * LINUX: The mupkern script comes standard with your Mupad install.

      * APPLE OS X:
         ???
"""

    def expect(self):
        return self._expect

    def console(self):
        mupad_console()

    def eval(self, code, strip=True):
        s = Expect.eval(self, code)
        return AsciiArtString(s)

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True,
                   need_output=True):
        if self._expect is None:
            self._start()
        if not need_output:
            E = self._expect
            E.sendline(line)
            return

        global seq
        seq += 1
        START = '__start__(%s+1)'%seq
        END = '__end__(%s+1)'%seq
        line = '%s; %s; %s;'%(START, line, END)
        START = '__start__(%s)'%(seq+1)
        END = '__end__(%s)'%(seq+1)

        E = self._expect
        E.sendline(line)
        E.expect(PROMPT)
        z = E.before
        i = z.find(START)
        if i == -1:
            raise RuntimeError, "%s\nError evaluating code in MuPAD"%z
        z = z[i+len(START)+2:]
        z = z.rstrip().rstrip(END).rstrip('"').rstrip().strip('\n').strip('\r').strip('\n').replace('\\\r\n','')
        i = z.find('Error: ')
        if i != -1:
            raise RuntimeError, z[i + 7:]
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
        i = out.find('Error: ')
        if i != -1:
            raise RuntimeError, out[i + 7:]

    def get(self, var):
        """
        Get the value of the variable var.
        """
        s = self.eval('%s'%var)
        i = s.find('=')
        return s[i+1:]

    def get_using_file(self, var):
        """
        Get the value of the variable var using a file.

        (I would make this the default for values that are bigger than
        a few thousand characters.  However, it's not at all obvious
        how to figure out if the string representation of an object is
        big ahead of time!  We assume it is for now, if the string
        used to create the object was big.)
        """
        s = self.eval('save %s, "%s"'%(var, tmp))
        s = open(tmp).read().replace('\\\n','')
        i = s.find('=')
        return s[i+2:-2]

    def _object_class(self):
        return MupadElement

    def _equality_symbol(self):
        return '=='

    def _assign_symbol(self):
        return ":="

    def _help(self, str):
        return os.popen('echo "?%s" | mupad -q'%str).read()

class MupadFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M._help(self._name)

class MupadFunctionElement(FunctionElement):
    def _sage_doc_(self):
        return self._obj.parent()._help(self._name)


class MupadElement(ExpectElement):
    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MupadFunctionElement(self, attrname)

    def trait_names(self):
        return self.parent().trait_names()

    def __repr__(self):
        self._check_valid()
        return self.parent().get(self._name)


# An instance
mupad = Mupad(script_subdirectory='user')

def reduce_load_mupad():
    return mupad

import os
def mupad_console():
    os.system('mupkern')


def __doctest_cleanup():
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()

