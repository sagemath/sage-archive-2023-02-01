r"""
Interface to GP/Pari

    Type \code{gp.[tab]} for a list of all the functions available
    from your Gp install.  Type \code{gp.[tab]?} for Gp's
    help about a given function.  Type \code{gp(...)} to create
    a new Gp object, and \code{gp.eval(...)} to run a string
    using Gp (and get the result back as a string).


EXAMPLES:
We illustrate objects that wrap GP objects (gp is the PARI interpreter):

    sage: M = gp('[1,2;3,4]')
    sage: M
    [1, 2; 3, 4]
    sage: M * M
    [7, 10; 15, 22]
    sage: M + M
    [2, 4; 6, 8]
    sage: M.matdet()
    -2


    sage: E = gp.ellinit([1,2,3,4,5])
    sage: E.ellglobalred()
    [10351, [1, -1, 0, -1], 1]
    sage: E.ellan(20)
    [1, 1, 0, -1, -3, 0, -1, -3, -3, -3, -1, 0, 1, -1, 0, -1, 5, -3, 4, 3]

    sage: primitive_root(7)
    3
    sage: x = gp("znlog( Mod(2,7), Mod(3,7))")
    sage: 3^x % 7
    2

    sage: print gp("taylor(sin(x),x)")
    x - 1/6*x^3 + 1/120*x^5 - 1/5040*x^7 + 1/362880*x^9 - 1/39916800*x^11 + 1/6227020800*x^13 - 1/1307674368000*x^15 + O(x^16)

GP has a powerful very efficient algorithm for numerical computation
of integrals.

    sage: gp("a = intnum(x=0,6,sin(x))")
    0.03982971334963397945434770208               # 32-bit
    0.039829713349633979454347702077075594548     # 64-bit
    sage: gp("a")
    0.03982971334963397945434770208               # 32-bit
    0.039829713349633979454347702077075594548     # 64-bit
    sage: gp.kill("a")
    sage: gp("a")
    a

Note that gp ASCII plots \emph{do} work in SAGE, as follows:

    sage.: print gp.eval("plot(x=0,6,sin(x))")
    [ plot of sin(x) ]

The GP interface reads in even very long input (using files) in a
robust manner, as long as you are creating a new object.

    sage: t = '"%s"'%10^10000   # ten thousand character string.
    sage: a = gp.eval(t)
    sage: a = gp(t)

In \sage, the PARI large galois groups datafiles should be installed
by default:
    sage: f = gp('x^9 - x - 2')
    sage: f.polgalois()
    [362880, -1, 34, "S9"]

AUTHORS:
    -- William Stein
    -- David Joyner (some examples)
    -- William Stein (2006-03-01): added tab completion for methods:
              gp.[tab] and x = gp(blah); x.[tab]
    -- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta
    -- William Stein (2006-05-17): updated to work with PARI 2.2.13-beta
"""

##########################################################################
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#
##########################################################################

from expect import Expect, ExpectElement, ExpectFunction, FunctionElement
from sage.misc.misc import verbose
from sage.libs.pari.all import pari
import sage.rings.all

class Gp(Expect):
    """
    Interface to the PARI gp interpreter.

    Type \code{gp.[tab]} for a list of all the functions available
    from your Gp install.  Type \code{gp.[tab]?} for Gp's
    help about a given function.  Type \code{gp(...)} to create
    a new Gp object, and \code{gp.eval(...)} to run a string
    using Gp (and get the result back as a string).
    """
    def __init__(self, stacksize=10000000,   # 10MB
                 maxread=100000, script_subdirectory=None,
                 logfile=None,
                 server=None,
                 init_list_length=1024):
        Expect.__init__(self,
                        name = 'pari',
                        prompt = '\\? ',
                        command = "gp --emacs --fast --quiet --stacksize %s"%stacksize,
                        maxread = maxread,
                        server=server,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        logfile=logfile,
                        eval_using_file_cutoff=50)
        self.__seq = 0
        self.__var_store_len = 0
        self.__init_list_length = init_list_length

    def _repr_(self):
        return 'GP/PARI interpreter'

    def __reduce__(self):
        return reduce_load_GP, tuple([])

    def __getattr__(self, attrname):
        if attrname[:1] == "_":
            raise AttributeError
        return GpFunction(self, attrname)

    def _quit_string(self):
        return "\\q"

    def _read_in_file_command(self, filename):
        return 'read("%s")'%filename

    def trait_names(self):
        try:
            b = self.__builtin
        except AttributeError:
            b = self.eval('?*').split()
            self.__builtin = b
        return b + self.eval('?0').split()

    def get_precision(self):
        """
        Return the current PARI precision for real number computations.
        """
        a = self.eval('\\p')
        i = a.find('=')
        j = a.find('significant')
        return int(a[i+1:j])

    get_real_precision = get_precision

    def set_precision(self, prec=None):
        """
        Sets the current PARI precision (in decimal digits) for real number computations,
        and returns the old one.
        """
        old = self.get_precision()
        self.eval('\\p %s'%int(prec))
        return old

    set_real_precision = set_precision

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True):
        line = line.strip()
        if len(line) == 0:
            return ''
        a = Expect._eval_line(self, line,
                              allow_use_file=allow_use_file,
                              wait_for_prompt=wait_for_prompt)
        if a.find("the PARI stack overflows") != -1:
            verbose("automatically doubling the PARI stack and re-executing current input line")
            self.eval("allocatemem()")
            return self._eval_line(line)
        else:
            return a

    def read(self, filename):
        s = self.eval('read("%s")'%filename)
        if 'error' in s:
            raise IOError, s


    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s;'%(var,value)
        out = self.eval(cmd)
        if out.find('***') != -1:
            raise TypeError, "Error executing code in GP/PARI:\nCODE:\n\t%s\nGP/PARI ERROR:\n%s"%(cmd, out)


    def get(self, var):
        """
        Get the value of the variable var.
        """
        return self.eval('print(%s)'%var)

    def kill(self, var):
        self.eval('kill(%s)'%var)

    #def xclear(self, var):
        #"""
        #Clear the variable named var.
        #"""
        #for varname based memory -- only 65000 variables and then dead.
        #self.eval('kill(%s)'%var)
        # for array-based memory this is best:
        #self.eval('%s=0'%var)
        # However, I've commented it out, since PARI doesn't seem
        # to ever free any memory on its stack anyways.
        # Killing variables as above takes a lot of time in some
        # cases, also.

    def _next_var_name(self):
        self.__seq += 1
        if self.__seq >= self.__var_store_len:
            if self.__var_store_len == 0:
                self.eval('sage=vector(%s,k,0);'%self.__init_list_length)
                self.__var_store_len = self.__init_list_length
            else:
                self.eval('sage=concat(sage, vector(%s,k,0));'%self.__var_store_len)
                self.__var_store_len *= 2
                verbose("doubling PARI/sage object vector: %s"%self.__var_store_len)
        return 'sage[%s]'%self.__seq

    def console(self):
        gp_console()

    def version(self):
        return gp_version()

    def _object_class(self):
        return GpElement

    def _true_symbol(self):
        return '1'

    def _false_symbol(self):
        return '0'

    def _equality_symbol(self):
        return '=='

    def help(self, command):
        return self.eval('?%s'%command).strip()

    def new_with_bits_prec(self, s, precision = 0):
        old_prec = self.get_real_precision()
        prec = int(precision/3.4) - 1
        if not precision:
            self.set_real_precision(prec)
            x = self(s)
            self.set_real_precision(old_prec)
        else:
            x = self(s)
        return x


class GpElement(ExpectElement):
    """
    EXAMPLES:
    This example illustrates dumping and loading GP elements to compressed strings.

        sage: a = gp(39393)
        sage: loads(a.dumps()) == a
        True

    Since dumping and loading uses the string representation of the object,
    it need not result in an identical object from the point of view of PARI:

        sage: E = gp('ellinit([1,2,3,4,5])')
        sage: loads(E.dumps()) == E
        False
        sage: loads(E.dumps())
        [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351, [-1.618909932267371342378000940, -0.3155450338663143288109995302 - 2.092547096911958607981689447*I, -0.3155450338663143288109995302 + 2.092547096911958607981689447*I]~, 2.780740013766729771063197627, -1.390370006883364885531598814 + 1.068749776356193066159263547*I, -1.554824121162190164275074561 + 3.415713103000000000000000000 E-29*I, 0.7774120605810950821375372806 - 1.727349756386839866714149879*I, 2.971915267817909670771647951] # 32-bit
        [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351, [-1.6189099322673713423780009396072169750, -0.31554503386631432881099953019639151248 - 2.0925470969119586079816894466366945829*I, -0.31554503386631432881099953019639151248 + 2.0925470969119586079816894466366945829*I]~, 2.7807400137667297710631976271813584994, -1.3903700068833648855315988135906792497 + 1.0687497763561930661592635474375038788*I, -1.5548241211621901642750745610982915039 + 7.9528267991764473360000000000000000000 E-39*I, 0.77741206058109508213753728054914575197 - 1.7273497563868398667141498789110695181*I, 2.9719152678179096707716479509361896060]  # 64-bit
        sage: E
        [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351, [-1.618909932267371342378000940, -0.3155450338663143288109995302 - 2.092547096911958607981689447*I, -0.3155450338663143288109995302 + 2.092547096911958607981689447*I]~, 2.780740013766729771063197627, -1.390370006883364885531598814 + 1.068749776356193066159263547*I, -1.554824121162190164275074561 + 3.415713103 E-29*I, 0.7774120605810950821375372806 - 1.727349756386839866714149879*I, 2.971915267817909670771647951]   # 32-bit
        [1, 2, 3, 4, 5, 9, 11, 29, 35, -183, -3429, -10351, 6128487/10351, [-1.6189099322673713423780009396072169750, -0.31554503386631432881099953019639151248 - 2.0925470969119586079816894466366945829*I, -0.31554503386631432881099953019639151248 + 2.0925470969119586079816894466366945829*I]~, 2.7807400137667297710631976271813584994, -1.3903700068833648855315988135906792497 + 1.0687497763561930661592635474375038788*I, -1.5548241211621901642750745610982915039 + 7.952826799176447336 E-39*I, 0.77741206058109508213753728054914575197 - 1.7273497563868398667141498789110695181*I, 2.9719152678179096707716479509361896060]  # 64-bit

    The two elliptic curves look the same, but internally the
    floating point numbers are slightly different.
    """

    def __long__(self):
        """
        Return Python long.
        """
        return long(str(self))

    def __float__(self):
        """
        Return Python float.
        """
        return float(pari(str(self)))

    def __bool__(self):
        P = self._check_valid()
        return P.eval('%s == 0'%(self.name())) == '1'

    def __getattr__(self, attrname):
        self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        return GpFunctionElement(self, attrname)

    def _complex_mpfr_field_(self, CC):
        """
        EXAMPLES:
            sage: CC(gp(1+15*I))
             1.00000000000000 + 15.0000000000000*I
            sage: CC(gp(11243.9812+15*I))
             11243.9812000000 + 15.0000000000000*I
            sage: ComplexField(10)(gp(11243.9812+15*I))
             1.1e4 + 15*I
        """
        GP = self.parent()
        orig = GP.get_real_precision()
        GP.set_real_precision(int(CC.prec()*3.33)+1)
        real = str(self.real()).replace(' E','e')
        imag = str(self.imag()).replace(' E','e')
        GP.set_real_precision(orig)
        return sage.rings.all.ComplexNumber( CC, real, imag )

    def __len__(self):
        return int(self.length())

    def __del__(self):
        return  # clearing object is pointless, since it wastes time, and PARI/GP
                # doesn't really free used memory anyways!

    # This is tempting -- but the (la)tex output is very very
    # out of date, e.g., for matrices it uses \pmatrix (which
    # causes an error if amsmath is loaded) and for rationals
    # it does nothing, etc.
    #def _latex_(self):
    #    P = self._check_valid()
    #    return P.eval('printtex(%s)'%self.name())

    def trait_names(self):
        return self.parent().trait_names()


class GpFunctionElement(FunctionElement):
    def _sage_doc_(self):
        M = self._obj.parent()
        return M.help(self._name)


class GpFunction(ExpectFunction):
    def _sage_doc_(self):
        M = self._parent
        return M.help(self._name)




def is_GpElement(x):
    return isinstance(x, GpElement)

# An instance
gp = Gp()

def reduce_load_GP():
    return gp

import os
def gp_console():
    os.system('gp')


def gp_version():
    """
    EXAMPLES:
        sage: gp.version()    # random output
        ((2, 3, 1), 'GP/PARI CALCULATOR Version 2.3.1 (0)')
    """
    v = gp.eval(r'\v')
    i = v.find("Version ")
    w = v[i+len("Version "):]
    i = w.find(' ')
    w = w[:i]
    t = tuple([int(n) for n in w.split('.')])
    k = v.find('\n')
    return t, v[:k].strip()
