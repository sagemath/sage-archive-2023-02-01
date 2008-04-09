r"""
R Interface

The following examples try to follow "An Introduction to R" which can
be found at http://cran.r-project.org/doc/manuals/R-intro.html .

EXAMPLES:

2 Simple manipulations; numbers and vectors

The simplest data structure in R is the numeric vector which
consists of an ordered collection of numbers.  To create a
vector named x using the R interface in SAGE, you pass the
R interpreter object a list or tuple of numbers.
    sage: x = r([10.4,5.6,3.1,6.4,21.7]); x
    [1] 10.4  5.6  3.1  6.4 21.7

You can invert elements of a vector x in R by using the
invert operator or by doing 1/x.
    sage: ~x
    [1] 0.09615385 0.17857143 0.32258065 0.15625000 0.04608295
    sage: 1/x
    [1] 0.09615385 0.17857143 0.32258065 0.15625000 0.04608295

The following assignemnt creates a vector y with 11 entries which
consists of two copies of x with a 0 in between.
    sage: y = r([x,0,x]); y
    [1] 10.4  5.6  3.1  6.4 21.7  0.0 10.4  5.6  3.1  6.4 21.7

=Vector Arithmetic=

The following command generates a new vector v of length 11 constructed
by adding together (element by element) 2*x repeated 2.2 times, y
repeated just once, and 1 repeated 11 times.
    sage: v = 2*x+y+1; v
    [1] 32.2 17.8 10.3 20.2 66.1 21.8 22.6 12.8 16.9 50.8 43.5

One can compute the sum of the elements of an R vector in the following
two ways:
    sage: sum(x)
    [1] 47.2
    sage: x.sum()
    [1] 47.2

One can calculate the sample variance of a list of numbers:
    sage: ((x-x.mean())^2/(x.length()-1)).sum()
    [1] 53.853
    sage: x.var()
    [1] 53.853


    sage: x.sort()
    [1] 3.1  5.6  6.4 10.4 21.7
    sage: x.min()
    [1] 3.1
    sage: x.max()
    [1] 21.7
    sage: x
    [1] 10.4  5.6  3.1  6.4 21.7



    sage: r(-17).sqrt()
    [1] NaN
    sage: r('-17+0i').sqrt()
    [1] 0+4.123106i

=Generating Regular Sequences=

    sage: r('1:10')
    [1] 1  2  3  4  5  6  7  8  9 10

Because 'from' is a keyword in Python, it can't be used
as a keyword argument.  Instead, 'from_' can be passed, and
R will recognize it as the correct thing.
    sage: r.seq(length=10, from_=-1, by=.2)
    [1] -1.0 -0.8 -0.6 -0.4 -0.2  0.0  0.2  0.4  0.6  0.8

    sage: x = r([10.4,5.6,3.1,6.4,21.7]);
    sage: x.rep(2)
    [1] 10.4  5.6  3.1  6.4 21.7 10.4  5.6  3.1  6.4 21.7
    sage: x.rep(times=2)
    [1] 10.4  5.6  3.1  6.4 21.7 10.4  5.6  3.1  6.4 21.7
    sage: x.rep(each=2)
    [1] 10.4 10.4  5.6  5.6  3.1  3.1  6.4  6.4 21.7 21.7


=Logical Vectors=


=Missing Values=
    sage: na = r('NA')
    sage: z = r([1,2,3,na])
    sage: z
    [1]  1  2  3 NA
    sage: ind = r.is_na(z)
    sage: ind
    [1] FALSE FALSE FALSE  TRUE
    sage: zero = r(0)
    sage: zero / zero
    [1] NaN
    sage: inf = r('Inf')
    sage: inf-inf
    [1] NaN
    sage: r.is_na(inf)
    [1] FALSE
    sage: r.is_na(inf-inf)
    [1] TRUE
    sage: r.is_na(zero/zero)
    [1] TRUE
    sage: r.is_na(na)
    [1] TRUE
    sage: r.is_nan(inf-inf)
    [1] TRUE
    sage: r.is_nan(zero/zero)
    [1] TRUE
    sage: r.is_nan(na)
    [1] FALSE


=Character Vectors=

    sage: labs = r.paste('c("X","Y")', '1:10', sep='""'); labs
    [1] "X1"  "Y2"  "X3"  "Y4"  "X5"  "Y6"  "X7"  "Y8"  "X9"  "Y10"


=Index vectors; selecting and modifying subsets of a data set=

    sage: na = r('NA')
    sage: x = r([10.4,5.6,3.1,6.4,21.7,na]); x
    [1] 10.4  5.6  3.1  6.4 21.7   NA
    sage: x['!is.na(self)']
    [1] 10.4  5.6  3.1  6.4 21.7


    sage: x = r([10.4,5.6,3.1,6.4,21.7,na]); x
    [1] 10.4  5.6  3.1  6.4 21.7   NA
    sage: (x+1)['(!is.na(self)) & self>0']
    [1] 11.4  6.6  4.1  7.4 22.7
    sage: x = r([10.4,-2,3.1,-0.5,21.7,na]); x
    [1] 10.4 -2.0  3.1 -0.5 21.7   NA
    sage: (x+1)['(!is.na(self)) & self>0']
    [1] 11.4  4.1  0.5 22.7

=Distributions=

sage: r.options(width="60");
sage: rr = r.dnorm(r.seq(-3,3,0.1))
sage: rr
 [1] 0.004431848 0.005952532 0.007915452 0.010420935
 [5] 0.013582969 0.017528300 0.022394530 0.028327038
 [9] 0.035474593 0.043983596 0.053990967 0.065615815
[13] 0.078950158 0.094049077 0.110920835 0.129517596
[17] 0.149727466 0.171368592 0.194186055 0.217852177
[21] 0.241970725 0.266085250 0.289691553 0.312253933
[25] 0.333224603 0.352065327 0.368270140 0.381387815
[29] 0.391042694 0.396952547 0.398942280 0.396952547
[33] 0.391042694 0.381387815 0.368270140 0.352065327
[37] 0.333224603 0.312253933 0.289691553 0.266085250
[41] 0.241970725 0.217852177 0.194186055 0.171368592
[45] 0.149727466 0.129517596 0.110920835 0.094049077
[49] 0.078950158 0.065615815 0.053990967 0.043983596
[53] 0.035474593 0.028327038 0.022394530 0.017528300
[57] 0.013582969 0.010420935 0.007915452 0.005952532
[61] 0.004431848

=Convert R datastructures to Python/Sage=

If possible native objects like lists, and then work directly
on them.

sage: rr = r.dnorm(r.seq(-3,3,0.1))
sage: sum(rr._sage_())
9.9772125168981081

Or you get a dictionary to be able to access all the information.

sage: rs = r.summary(r.c(1,4,3,4,3,2,5,1))
sage: rs
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  1.000   1.750   3.000   2.875   4.000   5.000
sage: rs._sage_() == {'DATA': [1, 1.75, 3, 2.875, 4, 5], '_Names': ['Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.'], '_r_class': 'table'}
True

AUTHORS:
    -- Mike Hansen (2007-11-01)
    -- William Stein
    -- Harald Schilly (2008-03-20)
"""

##########################################################################
#
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#                     2007 Mike Hansen   <mhansen@gmail.com>
#                     2008 Harald Schilly <harald.schilly@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#
##########################################################################

from keyword import iskeyword
from expect import Expect, ExpectElement, ExpectFunction, FunctionElement
from sage.misc.misc import verbose, DOT_SAGE
from sage.misc.sage_eval import sage_eval
import pexpect
from random import randrange
import re

COMMANDS_CACHE = '%s/r_commandlist.cache'%DOT_SAGE
PROMPT = 'ObFuSc4t3dpR0mp7> '
#there is a mirror network, but lets take #1 for now
RRepositoryURL = "http://cran.r-project.org/"
RFilteredPackages = [ '.GlobalEnv' ]
# crosscheck with https://svn.r-project.org/R/trunk/src/main/names.c
# but package:base should cover this. i think.
RBaseCommands = [ 'c',
"NULL", "NA", "True", "False", "Inf", "NaN"
]
# patterns for _sage_()
rel_re_param   = re.compile('\s([\w\.]+)\s=')
rel_re_xrange  = re.compile('([\d]+):([\d]+)')
rel_re_integer = re.compile('([^\d])([\d]+)L')

class R(Expect):
    """
    Hello! This is the R Engine from inside Sage speaking to you!

    R is a comprehensive collection of methods for statistics,
    modelling, bioinformatics, data analysis and much more.
    Coninues here: http://www.r-project.org/about.html

    Resources:
     * http://r-project.org/ provides more information about R.
     * http://rseek.org/ R's own search engine.

    Examples:
         sage: r.summary(r.c(1,2,3,111,2,3,2,3,2,5,4))
         Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
         1.00    2.00    3.00   12.55    3.50  111.00

    """
    def __init__(self, stacksize=10000000,   # 10MB
                 maxread=100000, script_subdirectory=None,
                 server_tmpdir = None,
                 logfile=None,
                 server=None,
                 init_list_length=1024):
        Expect.__init__(self,

                  # The capitalized version of this is used for printing.
                  name = 'r',

                  # This is regexp of the input prompt.  If you can change
                  # it to be very obfuscated that would be better.   Even
                  # better is to use sequence numbers.
                  # options(prompt=\"<prompt> \")
                  prompt = '> ', #default, later comes the change

                  # This is the command that starts up your program
                  command = "R --vanilla --quiet",

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
        self.__var_store_len = 0
        self.__init_list_length = init_list_length
        self._prompt_wait = [self._prompt]
        #only raw output, easier to read an R users are used to it
        #in general the output is just some sort of "rendering" and this
        #could mess up things
        #for standard view use the _sage_ object tranlation
        #self._remove_indices_re = re.compile(r"^\s*\[\d*\]", re.M)


    def _start(self):
        Expect._start(self)
        # width is line with, what's a good value? maximum is 10000!
        # pager needed to replace help view from less to printout
        # option device= is for plotting, is set to x11, NULL would be better?
        self._change_prompt(PROMPT)
        self.eval('options(prompt=\"%s\",width=100,pager="cat",device="png")'%PROMPT)
        self.expect().expect(PROMPT)
        self.eval('options(repos="%s")'%RRepositoryURL)
        # don't abort on errors, just raise them!
        # necessary for non-interactive execution
        self.eval('options(error = expression(NULL))')
        #if verbose: print "Initializing R, this takes some secs ..."

    def convRlistPyL(self,l):
        pl = []
        #l.replace('\n','')
        l = "".join(l.split("\n"))
        l = l[2:len(l)-1]
        for l in l.split(","):
            pl.append(l.strip().strip('"'))

        #if verbose: print "convRlistPyL: %s"%l
        return pl

    def install_packages(self,*a,**k):
        print """
Please do this outside the interactive environment!
Exit Sage and call R: $ ./sage -R
or use the console inside Sage: r.console()
"""

    # do this via R command
    #def getAllInstalledLibraries(self):
    #   """
    #   returns the list of all currently installed and recognized
    #   R packages.
    #   """
    #   if self._expect is None:
    #       self._start()
    #   # scanning all libraries
    #   # .Library is a path, list.files =^= bash's ls()
    #   # .libPaths() gives a vector of all library paths (could be more than one)
    #   # pack at the end puts it on the stack for return
    #   #wrong: self.execute('for(l in .libPaths()) { for(p in list.files(l)) { pack=c(pack,p) } }')
    #   RPackages = self.convRlistPyL(self.eval('dput(.packages(all=TRUE))'))
    #   #if verbose: print RPackages
    #   return RPackages


    def _repr_(self):
        return 'R Interpreter'

    def __reduce__(self):
        return reduce_load_r, tuple([])

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return RFunction(self, attrname)

    def _quit_string(self):
        return 'quit(save="no")'

    def _read_in_file_command(self, filename):
        raise NotImplementedError

    def read(self, filename):
        """
        Calls R's source function on a read-only file connection
        to @filename, which evaluates the given file in the
        current - possibly global - environment.
        """
        self.eval('source(file=file("%s",open="r"))'%filename)

    def _install_hints(self):
        return r"""
To install R point your browser to http://cran.r-project.org/mirrors.html
http://r-project.org/ provides more information about R.
"""

    ##########################
    def _source(self, s):
        """
        Tries to return the source code of a R function str
        as a string.

        EXAMPLES:
            print r._source("print.anova")
        """
        if s[-2:] == "()":
            s = s[-2:]
        return self.eval('%s'%s) #that's right, no quotes or brackets!

    def source(self, s):
        """
        Display the R source (if possible) about s.

        EXAMPLES:
            print r.source("print.anova")

        INPUT:
            s -- a string representing the function whose source code you
                 want
        """
        return self._source(s)
    ######################

    def kill(self, var):
        # [[send code that kills the variable with given name in the system.]]
        return "rm %s"%var

    def console(self):
        # run the console command (defined below).
        return "R --vanilla --quiet"

    #def version(self):
    #    # run the version command (defined below)
    #    s = self.eval('version')
    #    version_string = 'version.string '
    #    i = s.find(version_string)
    #    return s[len(version_string)+i:]

    def library(self,l): #overwrites library
        """
        load a library for R
        get a list with r.available_packages()
        """
        success = self.eval('require("%s")'%l)
        if not success:
           print "Could not load library %s"%l
        else:
           try:
               #we need to rebuild keywords!
               del self.__trait_names
           except AttributeError:
               pass
           self.trait_names()
           return success

    require = library #verwrites require

    def available_packages(self): #overwrites available_packages
        """
        list of ["package name", "version"] of available packages
        NOTE: needs internet connection, CRAN server is already specified
        """
        p = self.new('available_packages()')
        if p:
            p = p._sage_()
            s = p['_Dim'][0]
            print 'There are %s packages available.\n'%s
            l = [[a['DATA'][i],a['DATA'][s+1+i]] for i in xrange(0,s)]
            return l
        else:
            print 'nothing found, no return: %s'%p

    def _object_class(self):
        return RElement

    def _true_symbol(self):
        # return the string rep of truth, i.e., what the system outputs
        # when you type 1==1.
        return "TRUE"

    def _false_symbol(self):
        # return the string rep of truth, i.e., what the system outputs
        # when you type 1==2.
        return "FALSE"

    def _equality_symbol(self):
        # return the symbol for checking equality, e.g., == or eq.
        return "=="

    def help(self, command):
        """return help on a given command."""
        return self.eval('help("%s")'%command) #?cmd is only an unsafe shortcut

    def _assign_symbol(self):
        return " <- "

    def _left_list_delim(self):
        return "c("

    def _right_list_delim(self):
        return ")"

    def console(self):
        r_console()

    def function_call(self, function, args=[], kwargs={}):
        if function == '':
            raise ValueError, "function name must be nonempty"
        if function[:2] == "__":
            raise AttributeError
        if not isinstance(args, list):
            args = [args]
        for i in range(len(args)):
            if not isinstance(args[i], ExpectElement):
                args[i] = self.new(args[i])

        for key in kwargs:
            if not isinstance(kwargs[key], RElement):
                kwargs[key] = self.new(kwargs[key])

        return self.new("%s(%s)"%(function, ",".join([s.name() for s in args] + [self._sage_to_r_name(key)+'='+kwargs[key].name() for key in kwargs ] )))

    def call(self, function_name, *args, **kwargs):
        return self.function_call(function_name, args=args, kwargs=kwargs)

    def _an_element_impl(self):
        return self(0)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s <- %s'%(var,value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError, "Error executing code in R\nCODE:\n\t%s\nR ERROR:\n\t%s"%(cmd, out)

    def get(self, var):
        """
        Get the value of the variable var.
        """
        s = self.eval('%s'%var)
        #return self._remove_indices_re.sub("", s).strip()
        return s

    def na(self):
        """
        Returns the NA in R.

        EXAMPLES:
            sage: r.na()
            [1] NA
        """
        return self('NA')

    #def completions(self, s):
    #    """
    #    Return all commands that complete the command starting with the
    #    string s.   This is like typing s[Ctrl-T] in the R interpreter.
    #    """
    #    return []

    def _commands(self):
        """
        Return list of all commands defined in R.
        """
        # jump to R's source: https://svn.r-project.org/R/trunk/src/main/gram.c .y
        v = RBaseCommands
        try:
            ll = self.eval('dput(search())') # loaded libs
            ll = self.convRlistPyL(ll)
            #if verbose: print ll
            for lib in ll:
                if lib in RFilteredPackages: continue  #not the bad ones
                if lib.find("package:") != 0: continue #only packages
                #if verbose: print lib
                raw = self.eval('dput(objects("%s"))'%lib)
                #if verbose: print raw
                raw = self.convRlistPyL(raw)
                # "." in R are really like "_" in python
                raw = map(lambda x: x.replace(".","_"), raw)

                #TODO are there others? many of them are shortcuts or
                #should be done on another level, like selections in lists
                #instead of calling obj.[[( fun-args) or other crazy stuff like that

                #TODO further filtering, check if strings are now
                #really functions, in R: exists(s, mode = "function"))
                # (apply to vector with sapply(vec,func))

                #filter only python compatible identifiers
                from re import compile as re_compile
                valid = re_compile('[^a-zA-Z0-9_]+')
                raw = [x for x in raw if valid.search(x) == None]
                #if verbose: print raw[:min(len(raw),10)]
                print '.', #progress
                v += raw

        except RuntimeError:
            print "\n"*3
            print "*"*70
            print "WARNING: You do not have a working version of R installed!"
            print "*"*70
            v = []
        v.sort()
        return v

    def trait_names(self, verbose=True, use_disk_cache=False):
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
                print "\nBuilding R command completion list (this takes"
                print "a few seconds only the first time you do it)."
                print "To force rebuild later, delete %s."%COMMANDS_CACHE
            print "progress: ",
            v = self._commands()
            self.__trait_names = v
            if len(v) > 200:
                pass
                # R is actually installed.
                # TODO persist save doesn't work
            #    sage.misc.persist.save(v, COMMANDS_CACHE)
            return v

    def _sendstr(self, str):
        if self._expect is None:
            self._start()
        try:
            os.write(self._expect.child_fd, str)
        except OSError:
            self._crash_msg()
            self.quit()
            self._sendstr(str)

    def _synchronize(self):
        """
        Synchronize pexpect interface.

        This put a random integer (plus one!) into the output stream,
        then waits for it, thus resynchronizing the stream.  If the
        random integer doesn't appear within 1 second, R is sent
        interrupt signals.

        This way, even if you somehow left R in a busy state
        computing, calling _synchronize gets everything fixed.
        """
        if self._expect is None: return
        r = randrange(2147483647)
        s = str(r+1)
        cmd = "1+%s;\n"%r
        self._sendstr(cmd)
        try:
            self._expect_expr(timeout=0.5)
            if not s in self._expect.before:
                self._expect_expr(s,timeout=0.5)
                self._expect_expr(timeout=0.5)
        except pexpect.TIMEOUT, msg:
            self._interrupt()
        except pexpect.EOF:
            self._crash_msg()
            self.quit()

    def _expect_expr(self, expr=None, timeout=None):
        if expr is None:
            expr = self._prompt_wait
        if self._expect is None:
            self._start()
        try:
            if timeout:
                i = self._expect.expect(expr,timeout=timeout)
            else:
                i = self._expect.expect(expr)
            if i > 0:
                v = self._expect.before
                self.quit()
                raise ValueError, "%s\nComputation failed due to a bug in R -- NOTE: R had to be restarted."%v
        except KeyboardInterrupt, msg:
            #print self._expect.before
            i = 0
            while True:
                try:
                    print "Control-C pressed.  Interrupting R. Please wait a few seconds..."
                    self._sendstr('quit;\n'+chr(3))
                    self._sendstr('quit;\n'+chr(3))
                    self.interrupt()
                    self.interrupt()
                except KeyboardInterrupt:
                    i += 1
                    if i > 10:
                        break
                    pass
                else:
                    break
            raise KeyboardInterrupt, msg


    def plot(self, *args, **kwargs):
        #return self.function_call('plot', args=args, kwargs=kwargs)
        rp = RPlot()
        return rp.plot(self, *args, **kwargs)


    def eval(self, *args, **kwargs):
        # TODO split code at ";" outside of quotes and send them as individual
        #      lines withouth ";".
        if len(args) == 0:
            #called empty
            return "You called R with no arguments! e.g. try eval('1+1') or eval('matrix(seq(1,9),3)')."
        else:
            code = args[0]
            ret = Expect.eval(self, code, *args,**kwargs)
            return ret


    #####################
    def _r_to_sage_name(self, s):
        """
        EXAMPLES:
            sage: f = r._r_to_sage_name
            sage: f('t.test')
            't_test'
            sage: f('attr<-')
            'attr__'
            sage: f('parent.env<-')
            'parent_env__'
            sage: f('class')
            'class_'
        """
        s = s.replace('.', '_')
        s = s.replace('<-', '__')
        if iskeyword(s):
            s += '_'
        return s

    def _sage_to_r_name(self, s):
        """
        EXAMPLES:
            sage: f = r._sage_to_r_name
            sage: f('t_test')
            't.test'
            sage: f('attr__')
            'attr<-'
            sage: f('parent_env__')
            'parent.env<-'
            sage: f('class_')
            'class'
        """
        if len(s) > 1 and s[-2:] == "__":
            s = s[:-2] + '<-'
        if len(s) > 0 and s[-1] == '_':
            s = s[:-1]
        s = s.replace('_', '.')
        return s

    # TODO r.treerings is not a function, but a dataset
    #      tried to do this here but no luck
    #def __getattr__(self, s):
    #    """
    #    use R's 'get' function here
    #    """

    #    print '__getattr__(self,%s)'%s
    #    if '_' == s[:1] or 'get' == s:
    #        return self.s
    #    else:
    #        return self.get(s)

    def __getitem__(self, s):
        """
        Returns the RFunction s.

        EXAMPLES:
            sage: r['as.data.frame']
            as.data.frame
            sage: r['print']
            print
        """
        return RFunction(self, s, r_name=True)

# TODO how works plotting? ? ?
# sage: p = r.plot(r.c(1,2,3,4,3,4,3))
# sage: p.show()
# this just hangs/crashed
#####################
# png(file="myplot.png", bg="transparent")
# plot(1:10)
# rect(1, 5, 3, 7, col="white")
# dev.off()

from sage.structure.sage_object import SageObject
class RPlot(SageObject):
    r = None
    plotcmd = 'plot' #there is much more, ggplot2 package, lattice, mosaicplot, ...
    args = ""
    kwargs = ""

    #def __init__(self,Rint):
    #    self.r = Rint

    def save(self, filename=None, device = 'png', verbose=0, block=True, extra_opts=''):
        """
        save R plot

        format must be a suiteable R engine, use png!
        """
        self.r.execute('%s(file="%s", bg="transparent"'%(device,filename))
        #TODO it must be possible to issue multiple plot commands onto the same image
        self.r.execute('%s(%s, %s)'%(self.plotcmd, self.args, self.kwargs))
        self.r.execute('dev.off()')

    def show(self, verbose=0, extra_opts=''):
        import sage.plot.plot
        if sage.plot.plot.DOCTEST_MODE:
            filename = sage.misc.misc.graphics_filename()
            self.save(SAGE_TMP + '/test.png', verbose=verbose, extra_opts=extra_opts)
            return
        if sage.plot.plot.EMBEDDED_MODE:
            filename = sage.misc.misc.graphics_filename()
            self.save(filename, verbose=verbose, extra_opts=extra_opts)
            return
        filename = sage.misc.misc.tmp_filename() + '.png'
        self.save(filename, verbose=verbose, extra_opts=extra_opts)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(sage.misc.viewer.browser(), filename))

    def plot(self, rint, plotcmd = 'plot', *args, **kwargs):
        self.r         = rint
        self.plotcmd   = plotcmd
        self.args      = ", ".join(args)
        self.kwargs    = ", ".join(kwargs)
        return self




class RElement(ExpectElement):
    """
    Describe elements of your system here.
    """
    def trait_names(self):
        # This is if your system doesn't really have types.  If you have types
        # this function should only return the relevant methods that take self
        # as their first argument.
        # TODO: rewrite it, just take methods(class=class(self))
        return self.parent().trait_names()

    def __len__(self):
        """
        EXAMPLES:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: len(x)
            5
        """
        return sage_eval(self.parent().eval('dput(length(%s))'%self.name())[:-1] )

    def __getattr__(self, attrname):
        self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        return RFunctionElement(self, attrname)

    def __getitem__(self, n):
        """
        EXAMPLES:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x[0]
            numeric(0)
            sage: x[1]
            [1] 10.4
            sage: x[-1]
            [1] 5.6  3.1  6.4 21.7
            sage: x[-2]
            [1] 10.4  3.1  6.4 21.7
            sage: x[-3]
            [1] 10.4  5.6  6.4 21.7
            sage: x[4]
            [1] 6.4

        """
        P = self._check_valid()
        if isinstance(n, basestring):
            n = n.replace('self', self._name)
            return P.new('%s[%s]'%(self._name, n))
        elif not isinstance(n, tuple):
            return P.new('%s[%s]'%(self._name, n))
        else:
            return P.new('%s[%s]'%(self._name, str(n)[1:-1]))

    def __nonzero__(self):
        """
        bool(self) will only return True if self == 0 contains
        a FALSE.

        EXAMPLES:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: bool(x)
            True
            sage: y = r([0,0,0,0])
            sage: bool(y)
            False
            sage: bool(r(0))
            False
            sage: bool(r(1))
            True

        """
        return "FALSE" in repr(self==0)

    def _comparison(self, other, symbol):
        """
        TESTS:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x._comparison(10.4, "==")
            [1] TRUE FALSE FALSE FALSE FALSE

        """
        P = self.parent()
        other = P(other)
        return P('%s %s %s'%(self.name(), symbol, other.name()))

    def __eq__(self, other):
        """
        TESTS:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x == 10.4
            [1] TRUE FALSE FALSE FALSE FALSE
        """
        return self._comparison(other, "==")

    def __lt__(self, other):
        """
        TESTS:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x < 7
            [1] FALSE  TRUE  TRUE  TRUE FALSE

        """
        return self._comparison(other, "<")

    def __gt__(self, other):
        """
        TESTS:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x > 8
            [1] TRUE FALSE FALSE FALSE  TRUE

        """
        return self._comparison(other, ">")

    def __le__(self, other):
        """
        TESTS:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x <= 10.4
            [1] TRUE  TRUE  TRUE  TRUE FALSE

        """
        return self._comparison(other, "<=")

    def __ge__(self, other):
        """
        TESTS:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x >= 10.4
            [1] TRUE FALSE FALSE FALSE  TRUE

        """
        return self._comparison(other, ">=")

    def __ne__(self, other):
        """
        TESTS:
            sage: x = r([10.4,5.6,3.1,6.4,21.7])
            sage: x != 10.4
            [1] FALSE  TRUE  TRUE  TRUE  TRUE

        """
        return self._comparison(other, "!=")

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
    def dot_product(self, other):
        """
        Implements the notation self . other.

        EXAMPLES:
            sage: c = r.c(1,2,3,4)
            sage: c.dot_product(c.t())
                 [,1] [,2] [,3] [,4]
            [1,]    1    2    3    4
            [2,]    2    4    6    8
            [3,]    3    6    9   12
            [4,]    4    8   12   16

            sage: v = r([3,-1,8])
            sage: v.dot_product(v)
                 [,1]
            [1,]   74
        """
        P = self._check_valid()
        Q = P(other)
        #the R operator is %*% for matrix multiplication
        return P('%s %%*%% %s'%(self.name(), Q.name()))

        # TODO implement R's %o% operator for 'outer products'
        # outer is already robj.outer(robj) or
        # sage: v = r([3,-1,8])
        # sage: r.outer(v,v) || v.outer(v)
        #      [,1] [,2] [,3]
        # [1,]    9   -3   24
        # [2,]   -3    1   -8
        # [3,]   24   -8   64


    def _sage_(self):
        """
        Returns Sage representation of the R object.
        R objects are basic C structures, of different kind,
        that can be stacked together.
        This is similar to Python lists with variable objects,
        including lists in lists.
        If R lists have names, they are translated to a Python
        dictionary, anonymous list entries are called #{number}

        sage: rs = r.summary(r.c(1,4,3,4,3,2,5,1))
        sage: rs._sage_() == {'DATA': [1, 1.75, 3, 2.875, 4, 5], '_Names': ['Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.'], '_r_class': 'table'}
        True
        """
        self._check_valid()
        P = self.parent()

        #helper functions and patterns for _sage_
        #patterns are global in __init__
        #rel_re_param matches .Dim =, data.name = or .Names ...
        #rel_re_xrange translates \d:\d to xrange(d,d)
        #rel_re_integer digitsL to digits (L for raw long?)

        def ___subst_dots(x): return x.group().replace('.','_')
        def ___subst_xrange(x):
            g = x.groups()
            g1 = int(g[1]) + 1
            return 'xrange(%s,%s)'%(g[0],g1)
        def ___subst_integer(x): return '%sint(%s)'%x.groups()

        #thats the core of the trick: using dput
        # dput prints out the internal structure of R's data elements
        # options via .deparseOpts(control=...)
        # TODO dput also works with a file, if things get huge!
        exp = P.eval('dput(%s)'%self.name())

        #if verbose: print '_sage_ step 0: %s'%exp

        # preprocess expression
        # example what this could be:
        # structure(list(statistic = structure(0.233549683248457, .Names = "t"),
        # parameter = structure(5.58461538461538, .Names = "df"), p.value = 0.823657802106985,
        # conf.int = structure(c(-2.41722062247400, 2.91722062247400
        # ), conf.level = 0.95), estimate = structure(c(2.75, 2.5), .Names = c("mean of x",
        # "mean of y")), null.value = structure(0, .Names = "difference in means"),
        # alternative = "two.sided", method = "Welch Two Sample t-test",
        # data.name = "c(1, 2, 3, 5) and c(1, 2, 3, 4)"), .Names = c("statistic",
        # "parameter", "p.value", "conf.int", "estimate", "null.value",
        # "alternative", "method", "data.name"), class = "htest")

        # R's structure (from help):
        # structure(.Data, ...)
        #    .Data: an object which will have various attributes attached to it.
        #    ...: attributes, specified in 'tag=value' form, which will be
        #         attached to data.
        #For historical reasons (these names are used when deparsing),
        # attributes '".Dim"', '".Dimnames"', '".Names"', '".Tsp"' and
        # '".Label"' are renamed to '"dim"', '"dimnames"', '"names"',
        # '"tsp"' and '"levels"'.



        # we want this in a single line
        exp.replace('\n','')
        exp = "".join(exp.split("\n"))
        # python compatible parameters
        exp = re.sub(rel_re_param, ___subst_dots, exp)

        # we just call R's class attribute _r_class, is a string, no real meaning
        exp = re.sub(' class = "', ' _r_class = "',exp)
        # also for list and structure
        # TODO: check that we are outside of quotes ""
        exp = re.sub(' structure\(', ' _r_structure(', exp)
        exp = re.sub('^structure\(', '_r_structure(', exp) #special case

        exp = re.sub(' list\(', ' _r_list(', exp)
        exp = re.sub('\(list\(', '(_r_list(', exp)

        #if verbose: print '_sage_ step 1: %s'%exp

        # R's xrange is number:number with last one INCLUDING!
        exp = re.sub(rel_re_xrange, ___subst_xrange, exp)

        # If R want's to keep integers, it appends Ls to numbers
        # e.g. r.eval('dput(matrix(1:9,3))')
        #      'structure(1:9, .Dim = c(3L, 3L))'
        exp = re.sub(rel_re_integer, ___subst_integer, exp)

        # speciality, the call = <function name>
        # this is really bad
        # hack: replace the function name by the r.function name
        #       and wrap into quotes. (last one is TODO ?)
        exp = re.sub(', call = ', ' , call = r.', exp)
        exp = re.sub(' (sage\d+)', r' r.\1', exp)
        #callexp = re.match(', _r_call = "([^,]+),', exp)
        #print callexp
        # seems to work for
        # rr = r.summary(r.princomp(r.matrix(r.c(1,2,3,4,3,4,1,2,2),4)))
        # rr._sage_()
        # but the call expression get's evaluated. why?!? TODO


        #if verbose: print '_sage_ step 2: %s'%exp

        #class c(list):
        #   def __init__(self,*args):
        #      self[:] = list(args)
        #      self.__analyzeMe()
        #   def __analyzeMe(self):
        #      return self

        ################

        # translation:
        # c is an ordered list
        # list is a dictionary (where _Names give the entries names.
        #    map entries in names to (value, name) in each entry?
        # structure is .. see above .. strucuture(DATA,**kw)
        # TODO: thinking of just replacing c( with ( to get a long tuple?

        # this would work, but there are only 255 args max allowed in python -_-
        #def c(*args):
        #   return list(args)

        def __convert__RList__to__PyList(exp):
             from re import compile as re_compile
             from re import split   as re_split
             splt = re_compile('(c\(|\(|\))') # c( or ( or )
             lvl = 0
             ret = []
             for token in re_split(splt, exp):
                if token == 'c(':
                   ret.append('[ ')
                   lvl += 1
                elif token == '(':
                   ret.append(token)
                   if lvl > 0 : lvl += 1
                elif token == ')':
                   if lvl == 1:
                      ret.append(']')
                      lvl -= 1
                   else:
                      ret.append(token)
                      if lvl > 0: lvl -= 1
                else:
                    ret.append(token)
                #print "TOKEN ", token, " LEVEL ", lvl
             return ''.join(ret)

        exp = __convert__RList__to__PyList(exp)

        if verbose: print '_sage_ step 4: %s'%exp

        def _r_list(*ar, **kw):
           ret = dict(kw)
           i = 0
           for k in ar:
              ret['#%s'%i] = k
           return ret

        def _r_structure(__DATA__, **kw):
           if '_Dim' in kw: #we have a matrix
               # TODO what about more than 2 dimensions?
               #      additional checks!!
               try:
                  from sage.matrix.constructor import matrix
                  d = kw.get('_Dim')
                  # TODO: higher dimensions? happens often in statistics
                  if len(d) != 2:
                      raise TypeError
                  #since R does it the other way round, we assign
                  #transposed and then transpose the matrix :)
                  m = matrix(d[1], d[0], [i for i in __DATA__])
                  return m.transpose()
               except TypeError:
                  pass
           d = dict(DATA=__DATA__)
           d.update(kw)
           return d

        # R's NA and NULL stands for not available in a list of values
        NA   = None
        NULL = None
        FALSE = False
        TRUE  = True

        # the evaluation starts, like in misc/sage_eval.py
        return eval(exp)


    def _latex_(self):
        """
        Return Latex representation of this R object.

        This calls the tex command in R

        EXAMPLES:
            sage: r('sqrt(2) + 1/3 + asin(5)')._latex_()
            '\\sin^{-1}\\cdot5+\\sqrt{2}+{{1}\\over{3}}'
        """
        self._check_valid()
        P = self.parent()
        # latex is in Hmisc, this is currently not part of Sage's R!!!
        if P.loadLibrary('Hmisc'):
            #s = P._eval_line('latex(%s, file="");'%self.name(), reformat=False)
            s = P.eval('latex(%s, file="");'%self.name())
            return s
        else:
            return "Sorry, R package 'Hmisc' is not available. This one is necessary for\
                    R to latex translations :("

class RFunctionElement(FunctionElement):
    def _sage_doc_(self):
        M = self._obj.parent()
        return M.help(self._name)

    def _sage_src_(self):
        M = self._obj.parent()
        return M.source(self._name)

    def __call__(self, *args, **kwargs):
        return self._obj.parent().function_call(self._name, args=[self._obj] + list(args), kwargs=kwargs)


class RFunction(ExpectFunction):
    def __init__(self, parent, name, r_name=False):
        self._parent = parent
        if r_name:
            self._name = name
        else:
            self._name = parent._sage_to_r_name(name)

    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)

    def _sage_src_(self):
        M = self._parent
        return M.source(self._name)

    def __call__(self, *args, **kwargs):
        return self._parent.function_call(self._name, args=list(args), kwargs=kwargs)

def is_RElement(x):
    return isinstance(x, RElement)

# An instance
r = R()

def reduce_load_R():
    return r

import os
def r_console():
    # This will only spawn local processes
    os.system('R --vanilla')

## def r_version():
##     """
##     EXAMPLES:
##         sage: r.version()
##         ???
##     """
##     raise NotImplementedError
