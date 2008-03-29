r"""
Interface to LiE

LiE is a software package under development at CWI since
January 1988.  Its purpose is to enable mathematicians and
physicists to obtain on-line information as well as to
interactively perform computations of a Lie group theoretic
nature.  It focuses on the representation theory of complex
semisimple (reductive) Lie groups and algebras, and on the
structure of their Weyl groups and root systems.

Type \code{lie.[tab]} for a list of all the functions available
from your LiE install.  Type \code{lie.[tab]?} for LiE's
help about a given function.  Type \code{lie(...)} to create
a new LiE object, and \code{lie.eval(...)} to run a string
using LiE (and get the result back as a string).

To access the LiE interpreter directly, run lie_console().


EXAMPLES:
    sage: a4 = lie('A4')  # optional -- requires lie package
    sage: lie.diagram('A4')          # optional
    O---O---O---O
    1   2   3   4
    A4

    sage: lie.diagram(a4)            # optional
    O---O---O---O
    1   2   3   4
    A4

    sage: a4.diagram()               # optional
    O---O---O---O
    1   2   3   4
    A4

    sage: a4.Cartan()                # optional
         [[ 2,-1, 0, 0]
         ,[-1, 2,-1, 0]
         ,[ 0,-1, 2,-1]
         ,[ 0, 0,-1, 2]
         ]
    sage: lie.LR_tensor([3,1],[2,2]) # optional
    1X[5,3]


\subsection{Tutorial}

The following examples are taken from Section 2.1 of the LiE manual.

You can perform basic arithmetic operations in LiE.

    sage: lie.eval('19+68') # optional
    '87'
    sage: a = lie('1111111111*1111111111') # optional
    sage: a # optional
    1234567900987654321
    sage: a/1111111111 # optional
    1111111111
    sage: a = lie('345') # optional
    sage: a^2+3*a-5 # optional
    120055
    sage: _ / 7*a # optional
    5916750

Vectors in LiE are created using square brackets.  Notice that
the indexing in LiE is 1-based, unlike Python/SAGE which is
0-based.

    sage: v = lie('[3,2,6873,-38]') # optional
    sage: v # optional
    [3,2,6873,-38]
    sage: v[3] # optional
    6873
    sage: v+v # optional
    [6,4,13746,-76]
    sage: v*v # optional
    47239586
    sage: v+234786 # optional
    [3,2,6873,-38,234786]
    sage: v-3 # optional
    [3,2,-38]
    sage: v^v # optional
    [3,2,6873,-38,3,2,6873,-38]

You can also work with matrices in LiE.

    sage: m = lie('[[1,0,3,3],[12,4,-4,7],[-1,9,8,0],[3,-5,-2,9]]') # optional
    sage: m # optional
         [[ 1, 0, 3,3]
         ,[12, 4,-4,7]
         ,[-1, 9, 8,0]
         ,[ 3,-5,-2,9]
         ]
    sage: print lie.eval('*'+m._name) # optional
         [[1,12,-1, 3]
         ,[0, 4, 9,-5]
         ,[3,-4, 8,-2]
         ,[3, 7, 0, 9]
         ]

    sage: m^3 # optional
         [[ 220,   87, 81, 375]
         ,[-168,-1089, 13,1013]
         ,[1550,  357,-55,1593]
         ,[-854, -652, 98,-170]
         ]
    sage: v*m # optional
    [-6960,62055,55061,-319]
    sage: m*v # optional
    [20508,-27714,54999,-14089]
    sage: v*m*v # optional
    378549605
    sage: m+v # optional
         [[ 1, 0,   3,  3]
         ,[12, 4,  -4,  7]
         ,[-1, 9,   8,  0]
         ,[ 3,-5,  -2,  9]
         ,[ 3, 2,6873,-38]
         ]

    sage: m-2 # optional
         [[ 1, 0, 3,3]
         ,[-1, 9, 8,0]
         ,[ 3,-5,-2,9]
         ]


LiE handles multivariate (Laurent) polynomials.

    sage: lie('X[1,2]') # optional
    1X[1,2]
    sage: -3*_ # optional
    -3X[1,2]
    sage: _ + lie('4X[-1,4]') # optional
    4X[-1,4] - 3X[ 1,2]
    sage: _^2 # optional
    16X[-2,8] - 24X[ 0,6] +  9X[ 2,4]
    sage: lie('(4X[-1,4]-3X[1,2])*(X[2,0]-X[0,-4])') # optional
    -4X[-1, 0] + 3X[ 1,-2] + 4X[ 1, 4] - 3X[ 3, 2]
    sage: _ - _ # optional
    0X[0,0]


You can call LiE's built-in functions using lie.functionname .

    sage: lie.partitions(6) # optional
         [[6,0,0,0,0,0]
         ,[5,1,0,0,0,0]
         ,[4,2,0,0,0,0]
         ,[4,1,1,0,0,0]
         ,[3,3,0,0,0,0]
         ,[3,2,1,0,0,0]
         ,[3,1,1,1,0,0]
         ,[2,2,2,0,0,0]
         ,[2,2,1,1,0,0]
         ,[2,1,1,1,1,0]
         ,[1,1,1,1,1,1]
         ]
    sage: lie.diagram('E8') # optional
            O 2
            |
            |
    O---O---O---O---O---O---O
    1   3   4   5   6   7   8
    E8


You can define your own functions in LiE using lie.eval .  Once you've defined
a function (say f), you can call it using lie.f ; however, user-defined functions
do not show up when using tab-completion.

    sage: lie.eval('f(int x) = 2*x') # optional
    ''
    sage: lie.f(984) # optional
    1968
    sage: lie.eval('f(int n) = a=3*n-7; if a < 0 then a = -a fi; 7^a+a^3-4*a-57') # optional
    ''
    sage: lie.f(2) # optional
    -53
    sage: lie.f(5) # optional
    5765224



LiE's help can be accessed through lie.help('functionname') where
functionname is the function you want to receive help for.

   sage: print lie.help('diagram') # optional
   diagram(g).   Prints the Dynkin diagram of g, also indicating
      the type of each simple component printed, and labeling the nodes as
      done by Bourbaki (for the second and further simple components the
      labels are given an offset so as to make them disjoint from earlier
      labels). The labeling of the vertices of the Dynkin diagram prescribes
      the order of the coordinates of root- and weight vectors used in LiE.

This can also be accessed with lie.functionname? .



With the exception of groups, all LiE data types can be converted into
native SAGE data types by calling the .sage() method.

Integers:

    sage: a = lie('1234') # optional
    sage: b = a.sage(); b # optional
    1234
    sage: type(b) # optional
    <type 'sage.rings.integer.Integer'>

Vectors:

    sage: a = lie('[1,2,3]')# optional
    sage: b = a.sage(); b # optional
    [1, 2, 3]
    sage: type(b) # optional
    <type 'list'>

Matrices:

    sage: a = lie('[[1,2],[3,4]]') # optional
    sage: b = a.sage(); b # optional
    [1 2]
    [3 4]
    sage: type(b) # optional
    <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>


Polynomials:

    sage: a = lie('X[1,2] - 2*X[2,1]') # optional
    sage: b = a.sage(); b # optional
    -2*x0^2*x1 + x0*x1^2
    sage: is_MPolynomial(b) # optional
    True


Text:

    sage: a = lie('"text"') # optional
    sage: b = a.sage(); b # optional
    'text'
    sage: type(b) # optional
    <type 'str'>


LiE can be programmed using the SAGE interface as well. Section 5.1.5
of the manual gives an example of a function written in LiE's language
which evaluates a polynomial at a point.  Below is a (roughly) direct
translation of that program into Python / SAGE.

    sage: def eval_pol(p, pt): # optional
    ...       s = 0
    ...       for i in range(1,p.length().sage()+1):
    ...           m = 1
    ...           for j in range(1,pt.size().sage()+1):
    ...               m *= pt[j]^p.expon(i)[j]
    ...           s += p.coef(i)*m
    ...       return s
    sage: a = lie('X[1,2]') # optional
    sage: b1 = lie('[1,2]') # optional
    sage: b2 = lie('[2,3]') # optional
    sage: eval_pol(a, b1) # optional
    4
    sage: eval_pol(a, b2) # optional
    18



AUTHORS:
    -- Mike Hansen 2007-08-27
    -- William Stein (template)
"""

##########################################################################
#
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#
##########################################################################

from expect import Expect, ExpectElement, ExpectFunction, FunctionElement, AsciiArtString
from sage.misc.misc import verbose, DOT_SAGE, SAGE_LOCAL


COMMANDS_CACHE = '%s/lie_commandlist_cache.sobj'%DOT_SAGE
HELP_CACHE = '%s/lie_helpdict_cache.sobj'%DOT_SAGE

class LiE(Expect):
    r"""
    Interface to the LiE interpreter.

    Type \code{lie.[tab]} for a list of all the functions available
    from your LiE install.  Type \code{lie.[tab]?} for LiE's
    help about a given function.  Type \code{lie(...)} to create
    a new LiE object, and \code{lie.eval(...)} to run a string
    using LiE (and get the result back as a string).

    """
    def __init__(self,
                 maxread=100000, script_subdirectory=None,
                 logfile=None,
                 server=None):
        Expect.__init__(self,

                        # The capitalized versionof this is used for printing.
                        name = 'LiE',

                        # This is regexp of the input prompt.  If you can change
                        # it to be very obfuscated that would be better.   Even
                        # better is to use sequence numbers.
                        prompt = '> ',

                        # This is the command that starts up your program
                        command = "bash "+ SAGE_LOCAL + "/bin/lie",

                        maxread = maxread,
                        server=server,
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

        self.__trait_names_dict = None
        self.__trait_names_list = None
        self.__help_dict = None

    def __read_info_files(self):
        use_disk_cache = True

        import sage.misc.persist
        if use_disk_cache:
            try:
                trait_dict = sage.misc.persist.load(COMMANDS_CACHE)
                help_dict = sage.misc.persist.load(HELP_CACHE)
                v = []
                for key in trait_dict:
                    v += trait_dict[key]
                self.__trait_names_list = v
                self.__trait_names_dict = trait_dict
                self.__help_dict = help_dict
                return
            except IOError:
                pass


        #Go through INFO.3 and get the necessary information
        filenames = ['INFO.3', 'INFO.0']
        commands = {}
        commands['vid'] = []
        help = {}


        for f in filenames:
            filename = SAGE_LOCAL + "/lib/lie/" + f
            info = open(filename)
            prev_command = ""
            help_text = ""
            for line in info:
                #If the line doesn't start with an "@", then
                #it is part of the help text for the previous
                #command
                if len(line) == 0 or line[0] != "@":
                    if prev_command != "":
                        help_text += line
                    continue


                #Do not add not completions that do not start with an
                #alphabetical character or that contain 'silence'
                if len(line) > 1 and (not line[1].isalpha() or line.find('silence') != -1):
                    help[prev_command] = help.get(prev_command, "") + help_text
                    help_text = ""
                    prev_command = ""
                    continue


                #At this point we should be at the start of a new
                #command definition


                #Get the type of the first argument of the command
                i = line.find('(')
                if line[i+1] == ")":
                    t = 'vid'
                else:
                    t = line[i+1:i+4]

                #Save the help text for the command
                help[prev_command] = help.get(prev_command, "") + help_text
                help_text = ""
                prev_command = line[1:i]

                #Add the commad
                if t in commands:
                    commands[t].append(line[1:i])
                else:
                    commands[t] = [ line[1:i] ]

            #Take care of the last help text which doesn't get processed
            #since there's no following @ symbol
            help[prev_command] = help.get(prev_command, "") + help_text

        info.close()


        #Build the list of all possible command completions
        l = []
        for key in commands:
            l += commands[key]

        #Save the data
        self.__trait_names_dict = commands
        self.__trait_names_list = l
        self.__help_dict = help

        #Write them to file
        if use_disk_cache:
            sage.misc.persist.save(commands, COMMANDS_CACHE)
            sage.misc.persist.save(help, HELP_CACHE)

    def _repr_(self):
        return 'LiE Interpreter'

    def __reduce__(self):
        return reduce_load_lie, tuple([])

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return LiEFunction(self, attrname)

    def _quit_string(self):
        return 'quit'

    def _read_in_file_command(self, filename):
        raise NotImplementedError


    def trait_names(self, type=None, verbose=False, use_disk_cache=True):
        if self.__trait_names_dict is None:
            self.__read_info_files()
        if type:
            return self.__trait_names_dict[type]
        else:
            return self.__trait_names_list

    def _an_element_impl(self):
        return self(0)

    def read(self, filename):
        # [[implement loading of the contents of filename into the system]]
        self.eval('read %s'%filename)


    def kill(self, var):
        # [[send code that kills the variable with given name in the system.]]
        pass

    def console(self):
        # run the console command (defined below).
        lie_console()

    def version(self):
        # run the version command (defined below)
        pass

    def _object_class(self):
        return LiEElement

    def _true_symbol(self):
        # return the string rep of truth, i.e., what the system outputs
        # when you type 1==1.
        return '1'

    def _false_symbol(self):
        # return the string rep of false, i.e., what the system outputs
        # when you type 1==2.
        return '0'

    def _equality_symbol(self):
        # return the symbol for checking equality, e.g., == or eq.
        return '=='

    def help(self, command):
        # return help on a given command.
        if self.__help_dict is None:
            self.__read_info_files()
        try:
            return self.__help_dict[command]
        except KeyError:
            return "Could not find help for " + command

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True):
        #if line.find('\n') != -1:
        #    raise ValueError, "line must not contain any newlines"
        if allow_use_file and self._eval_using_file_cutoff and len(line) > self._eval_using_file_cutoff:
            return self._eval_line_using_file(line, tmp)
        try:
            if self._expect is None:
                self._start()
            E = self._expect
            try:
                if len(line) >= 4096:
                    raise RuntimeError, "Sending more than 4096 characters with %s on a line may cause a hang and you're sending %s characters"%(self, len(line))
                E.sendline(line)
                if wait_for_prompt == False:
                    return ''

            except OSError, msg:
                raise RuntimeError, "%s\nError evaluating %s in %s"%(msg, line, self)

            if len(line)>0:
                try:
                    if isinstance(wait_for_prompt, basestring):
                        E.expect(wait_for_prompt)
                    else:
                        E.expect(self._prompt)
                except pexpect.EOF, msg:
                    try:
                        if self._read_in_file_command(tmp) in line:
                            raise pexpect.EOF, msg
                    except NotImplementedError:
                        pass
                    if self._quit_string() in line:
                        # we expect to get an EOF if we're quitting.
                        return ''
                    raise RuntimeError, "%s\n%s crashed executing %s"%(msg,
                                                   self, line)
                out = E.before
            else:
                out = '\n\r'
        except KeyboardInterrupt:
            self._keyboard_interrupt()
            raise KeyboardInterrupt, "Ctrl-c pressed while running %s"%self
        i = out.find("\n")
        j = out.rfind("\r")

        err = max( out.find("\n(in"), out.find('not defined'), out.find('Argument types')  )
        if err != -1:
            raise RuntimeError, "An error occured running a LiE command:\n%s"%(out[i+1:j].replace('\r\n','\n'))

        return out[i+1:j].replace('\r\n','\n')


    def eval(self, code, strip=True, **kwds):
        s = Expect.eval(self,code, strip=True, **kwds)
        #return s.strip()
        if len(s) > 0 and s.find("\n") != -1:
            return s
        else:
            return s.strip()

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s'%(var,value)
        out = self.eval(cmd)
        i = min( out.find('not defined'), out.find('\(in'), out.find('Argument types') )
        if i != -1:
            raise RuntimeError, out

    def get(self, var):
        """
        Get the value of the variable var.
        """
        s = self.eval('%s'%var)
        return s

    def get_using_file(self, var):
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
        #Handle the functions that do not return a value that can
        #be assigned to a variable
        if function in ['diagram', 'setdefault', 'print_tab', 'type', 'factor', 'void', 'gcol']:
            return AsciiArtString(self.eval("%s(%s)"%(function, ",".join([s.name() for s in args]))))
        #Otherwise, create a new object
        else:
            return self.new("%s(%s)"%(function, ",".join([s.name() for s in args])))

class LiEElement(ExpectElement):
    """
    Describe elements of your system here.
    """
    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return LiEFunctionElement(self, attrname)

    def trait_names(self):
        # This is if your system doesn't really have types.  If you have types
        # this function should only return the relevant methods that take self
        # as their first argument.
        return self.parent().trait_names(type=self.type())

    def type(self):
        t = self.parent().eval('type(%s)'%self._name)
        i = t.find(':')
        return t[i+1:].strip()


    def _matrix_(self, R):
        if self.type() == 'mat':
            return R( eval( str(self).replace('\n','').strip() ) )
        else:
            raise ValueError, "not a matrix"

    def _vector_(self, R):
        if self.type() == 'vec':
            return R(eval(str(self)))
        else:
            raise ValueError, "not a vector"

    def __cmp__(self, other):
        P = self.parent()
        if P.eval("%s %s %s"%(self.name(), P._equality_symbol(),
                                 other.name())) == P._true_symbol():
            return 0
        elif P.eval("%s %s %s"%(self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
            return -1
        elif P.eval("%s %s %s"%(self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
            return 1
        else:
            return -1  # everything is supposed to be comparable in Python, so we define
                       # the comparison thus when no comparable in interfaced system.

    def sage(self):
        t = self.type()
        if t == 'grp':
            raise ValueError, "cannot convert Lie groups to native SAGE objects"
        elif t == 'mat':
            import sage.matrix.constructor
            return  sage.matrix.constructor.matrix( eval( str(self).replace('\n','').strip())  )
        elif t == 'pol':
            import sage.misc.misc
            from sage.rings.all import PolynomialRing, QQ

            #Figure out the number of variables
            s = str(self)
            open_bracket = s.find('[')
            close_bracket = s.find(']')
            nvars = len(s[open_bracket:close_bracket].split(','))

            #create the polynomial ring
            R = PolynomialRing(QQ, nvars, 'x')
            x = R.gens()
            pol = R(0)

            #Split up the polynomials into terms
            terms = []
            for termgrp in s.split(' - '):
                #The first entry in termgrp has
                #a negative coefficient
                termgrp = "-"+termgrp.strip()
                terms += termgrp.split('+')
            #Make sure we don't accidently add a negative
            #sign to the first monomial
            if s[0] != "-":
                terms[0] = terms[0][1:]

            #go through all the terms in s
            for term in terms:
                xpos = term.find('X')
                coef = eval(term[:xpos].strip())
                exps = eval(term[xpos+1:].strip())
                monomial = sage.misc.misc.prod(map(lambda i: x[i]**exps[i] , range(nvars)))
                pol += coef * monomial

            return pol
        elif t == 'tex':
            return repr(self)
        elif t == 'vid':
            return None
        else:
            return ExpectElement.sage(self)


class LiEFunctionElement(FunctionElement):
    def _sage_doc_(self):
        M = self._obj.parent()
        return M.help(self._name)


class LiEFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)



def is_LiEElement(x):
    return isinstance(x, LiEElement)

# An instance
lie = LiE()

def reduce_load_lie():
    return lie

import os
def lie_console():
    os.system('bash `which lie`')


def lie_version():
    raise NotImplementedError
