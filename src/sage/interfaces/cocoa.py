r"""nodoctest
Interface to Cocoa 4

[[Describe the math software you are interfacing with here.]]

[[Replace this by something relevant to your system.]]
    Type \code{cocoa.[tab]} for a list of all the functions available
    from your Cocoa install.  Type \code{cocoa.[tab]?} for Cocoa's
    help about a given function.  Type \code{cocoa(...)} to create
    a new Cocoa object, and \code{cocoa.eval(...)} to run a string
    using Cocoa (and get the result back as a string).


EXAMPLES:

[[Go through a standard tutorial for your software package
and do it via your SAGE interface.]]

AUTHORS:
    -- William Stein (2006-10)
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

class Cocoa(Expect):
    """
    [[Some basic help about your system.  This is what
      will be displayed when somebody write cocoa?.]]
    """
    def __init__(self,
                 maxread=100000, script_subdirectory=None,
                 logfile=None,
                 server=None):
        Expect.__init__(self,

                        # The capitalized versionof this is used for printing.
                        name = 'cocoa',

                        # This is regexp of the input prompt.  If you can change
                        # it to be very obfuscated that would be better.   Even
                        # better is to use sequence numbers.
                        prompt = '-------------------------------',

                        # This is the command that starts up your program
                        command = "cocoa",

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
        return 'A Cocoa 4 Interpreter'

    def __reduce__(self):
        return reduce_load_cocoa, tuple([])

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return CocoaFunction(self, attrname)

    def _quit_string(self):
        return "\\q"

    def _read_in_file_command(self, filename):
        return 'read("%s")'%filename

    def trait_names(self):
        ## [[implement giving a list of all functions and identifiers in the system]]
        pass

    def read(self, filename):
        # [[implement loading of the contents of filename into the system]]
        pass


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
        return CocoaElement

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
        pass

class CocoaElement(ExpectElement):
    """
    Describe elements of your system here.
    """
    def trait_names(self):
        # This is if your system doesn't really have types.  If you have types
        # this function should only return the relevant methods that thake self
        # as their firce argument.
        return self.parent().trait_names()


class GpFunctionElement(FunctionElement):
    def _sage_doc_(self):
        M = self._obj.parent()
        return M.help(self._name)


class GpFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)



def is_CocoaElement(x):
    return isinstance(x, CocoaElement)

# An instance
cocoa = Cocoa()

def reduce_load_Cocoa():
    return cocoa

import os
def cocoa_console():
    os.system('cocoa')

def cocoa_version():
    """
    EXAMPLES:
        sage: cocoa.version()            # optional

    """
