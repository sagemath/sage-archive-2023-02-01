"""
Interface to Maple

You must have the maple interpreter installed and available as
the command "maple" in your PATH in order to use this interface.

EXAMPLES:
    sage: maple('3 * 5')                       # needs optional maple
    15
    sage: maple.eval('ifactor(2005)')          # needs optional maple
    '``(5)*``(401)'
    sage: maple.ifactor(2005)                  # needs optional maple
    ``(5)*``(401)
    sage: maple.fsolve('x^2=cos(x)+4', 'x=0..5')   # needs optional maple
    1.914020619
    sage: maple.factor('x^5 - y^5')            # needs optional maple
    (x-y)*(x^4+x^3*y+x^2*y^2+x*y^3+y^4)

If the string "error" (case insensitive) occurs in the
output of anything from Maple, a RuntimeError exception is raised.

"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

from expect import Expect, ExpectElement, tmp

from sage.misc.misc import verbose

class Maple(Expect):
    """
    Interface to the Maple interpreter.
    """
    def __init__(self, maxread=100, script_subdirectory="", logfile=None):
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
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=1)  # very important that this is 1
        # It's very important to use file i/o for everything,
        # since maple stupid command line interface always
        # dumps you into the editor when an error occurs,
        # and I can find no way to turn it off!!

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

    def _eval_line(self, line, allow_use_file=True):
        line += ';'
        z = Expect._eval_line(self, line, allow_use_file=allow_use_file).\
                                  replace('\\\n','').strip()
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

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s:=%s;'%(var,value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError, "Error executing code in Maple\nCODE:\n\t%s\nMaple ERROR:\n\t%s"%(cmd, out)

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
        return MapleElement

    def help(self, str):
        """
        Display Maple help about str.  This is the same as typing "?str" in
        the Maple console.

        INPUT:
            str -- a string to search for in the maple help system
        """
        os.system('echo "?%s" | maple -q |less'%str)

    def with(self, package):
        """
        Make a package of Maple procedures available in the
        interpreter.

        INPUT:
            package -- string

        EXAMPLES:
        Some functions are unknown to Maple until you use with to include
        the appropriate package.

            sage: maple('partition(10)')               # optional
            partition(10)
            sage: maple('bell(10)')                    # optional
            bell(10)
            sage: maple.with('combinat')               # optional
            sage: maple('partition(10)')               # optional
            [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 2], [1, 1, 1, 1, 1,
            1, 2, 2], [1, 1, 1, 1, 2, 2, 2], [1, 1, 2, 2, 2, 2], [2, 2, 2, 2, 2], [1, 1, 1
            , 1, 1, 1, 1, 3], [1, 1, 1, 1, 1, 2, 3], [1, 1, 1, 2, 2, 3], [1, 2, 2, 2, 3],
            [1, 1, 1, 1, 3, 3], [1, 1, 2, 3, 3], [2, 2, 3, 3], [1, 3, 3, 3], [1, 1, 1, 1,
            1, 1, 4], [1, 1, 1, 1, 2, 4], [1, 1, 2, 2, 4], [2, 2, 2, 4], [1, 1, 1, 3, 4],
            [1, 2, 3, 4], [3, 3, 4], [1, 1, 4, 4], [2, 4, 4], [1, 1, 1, 1, 1, 5], [1, 1, 1
            , 2, 5], [1, 2, 2, 5], [1, 1, 3, 5], [2, 3, 5], [1, 4, 5], [5, 5], [1, 1, 1, 1
            , 6], [1, 1, 2, 6], [2, 2, 6], [1, 3, 6], [4, 6], [1, 1, 1, 7], [1, 2, 7], [3,
            7], [1, 1, 8], [2, 8], [1, 9], [10]]
            sage: maple('bell(10)')                   # optional
            115975
            sage: maple('fibonacci(10)')              # optional
            55
        """
        self.eval('with(%s)'%package)

    load = with

    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
        # Unfortunately MAPLE does not have a clear command.
        # The next best thing is to set equal to the constant
        # 0, so that memory will be freed.
    #    self.eval("%s=0;"%var)

class MapleElement(ExpectElement):
    def __float__(self):
        M = self.parent()
        return float(maple.eval('evalf(%s)'%self.name()))

    def __repr__(self):
        self._check_valid()
        return self.parent().get(self._name)

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
