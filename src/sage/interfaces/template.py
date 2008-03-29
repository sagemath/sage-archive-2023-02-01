r"""nodoctest [[remove the nodoctest, so doctests will run]]
Interface Template

[[Describe the math software you are interfacing with here.]]

[[Replace this by something relevant to your system.]]
    Type \code{gp.[tab]} for a list of all the functions available
    from your Gp install.  Type \code{gp.[tab]?} for Gp's
    help about a given function.  Type \code{gp(...)} to create
    a new Gp object, and \code{gp.eval(...)} to run a string
    using Gp (and get the result back as a string).


EXAMPLES:

[[Go through a standard tutorial for your software package
and do it via your SAGE interface.]]

AUTHORS:
    -- William Stein (template)
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

from expect import Expect, ExpectElement, ExpectFunction, FunctionElement
from sage.misc.misc import verbose

class MySystem(Expect):
    """
    [[Some basic help about your system.  This is what
      will be displayed when somebody write mysystem?.]]
    """
    def __init__(self,
                 maxread=100000, script_subdirectory=None,
                 logfile=None,
                 server=None,
                 server_tmpdir=None):
        Expect.__init__(self,

                        # The capitalized version of this is used for printing.
                        name = 'mysystem',

                        # This is regexp of the input prompt.  If you can change
                        # it to be very obfuscated that would be better.   Even
                        # better is to use sequence numbers.
                        prompt = '>> ',

                        # This is the command that starts up your program
                        command = "mysystem",

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

    def _repr_(self):
        return 'MySystem Interpreter'

    def __reduce__(self):
        return reduce_load_mysystem, tuple([])

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return MySystemFunction(self, attrname)

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
        return MySystemElement

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

class MySystemElement(ExpectElement):
    """
    Describe elements of your system here.
    """
    def trait_names(self):
        # This is if your system doesn't really have types.  If you have types
        # this function should only return the relevant methods that take self
        # as their first argument.
        return self.parent().trait_names()


class MySystemFunctionElement(FunctionElement):
    def _sage_doc_(self):
        M = self._obj.parent()
        return M.help(self._name)


class MySystemFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)



def is_MySystemElement(x):
    return isinstance(x, MySystemElement)

# An instance
mysystem = MySystem()

def reduce_load_MySystem():
    return mysystem

import os
def mysystem_console():
    # This will only spawn local processes
    os.system('mysystem')


def mysystem_version():
    """
    EXAMPLES:
        sage: mysystem.version()
        ???
    """
    raise NotImplementedError
