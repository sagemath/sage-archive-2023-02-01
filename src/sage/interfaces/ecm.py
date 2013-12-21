r"""
The Elliptic Curve Factorization Method

\sage includes GMP-ECM, which is a highly optimized implementation
of Lenstra's elliptic curve factorization method.  See
\url{http://ecm.gforge.inria.fr/}
for more about GMP-ECM.
"""

import os, pexpect
from math import ceil, floor

from sage.rings.integer import Integer
from sage.misc.misc import verbose, tmp_filename
from sage.misc.decorators import rename_keyword

import cleaner

import sage.misc.package

def nothing():
    pass

class ECM:
    def __init__(self, B1=10, B2=None, **kwds):
        r"""
        Create an interface to the GMP-ECM elliptic curve method
        factorization program.

        See \url{http://ecm.gforge.inria.fr/}.

        AUTHORS:
            These people wrote GMP-ECM:
              Pierrick Gaudry, Jim Fougeron,
              Laurent Fousse, Alexander Kruppa,
              Dave Newman, Paul Zimmermann

        William Stein and Robert Bradshaw -- wrote the Sage interface to GMP-ECM

        INPUT:
            B1 -- stage 1 bound
            B2 -- stage 2 bound (or interval B2min-B2max)

            x0 -- x        use x as initial point
            sigma -- s     use s as curve generator [ecm]
            A -- a         use a as curve parameter [ecm]
            k -- n         perform >= n steps in stage 2
            power -- n     use x^n for Brent-Suyama's extension
            dickson -- n   use n-th Dickson's polynomial for Brent-Suyama's extension
            c -- n         perform n runs for each input
            pm1 --         perform P-1 instead of ECM
            pp1 --         perform P+1 instead of ECM
            q --           quiet mode
            v --           verbose mode
            timestamp --   print a time stamp with each number
            mpzmod --      use GMP's mpz_mod for mod reduction
            modmuln --     use Montgomery's MODMULN for mod reduction
            redc --        use Montgomery's REDC for mod reduction
            nobase2 --     disable special base-2 code
            base2 -- n     force base 2 mode with 2^n+1 (n>0) or 2^n-1 (n<0)
            save -- file   save residues at end of stage 1 to file
            savea -- file  like -save, appends to existing files
            resume -- file resume residues from file, reads from
                           stdin if file is "-"
            primetest --   perform a primality test on input
            treefile -- f  store product tree of F in files f.0 f.1 ...
            i -- n         increment B1 by this constant on each run
            I -- f         auto-calculated increment for B1 multiplied by 'f' scale factor
            inp -- file    Use file as input (instead of redirecting stdin)
            b --           Use breadth-first mode of file processing
            d --           Use depth-first mode of file processing (default)
            one --         Stop processing a candidate if a factor is found (looping mode)
            n --           run ecm in 'nice' mode (below normal priority)
            nn --          run ecm in 'very nice' mode (idle priority)
            t -- n         Trial divide candidates before P-1, P+1 or ECM up to n
            ve -- n        Verbosely show short (< n character) expressions on each loop
            cofdec --      Force cofactor output in decimal (even if expressions are used)
            B2scale -- f   Multiplies the default B2 value by f
            go -- val      Preload with group order val, which can be a simple expression,
                           or can use N as a placeholder for the number being factored.
            prp -- cmd     use shell command cmd to do large primality tests
            prplen -- n    only candidates longer than this number of digits are 'large'
            prpval -- n    value>=0 which indicates the prp command foundnumber to be PRP.
            prptmp -- file outputs n value to temp file prior to running (NB. gets deleted)
            prplog -- file otherwise get PRP results from this file (NB. gets deleted)
            prpyes -- str  literal string found in prplog file when number is PRP
            prpno -- str   literal string found in prplog file when number is composite
        """
        self.__cmd = self.__startup_cmd(B1, B2, kwds)

    def __startup_cmd(self, B1, B2, kwds):
        options = ' '.join(['-%s %s'%(x,v) for x, v in kwds.iteritems()])
        s = 'ecm %s %s '%(options, B1)
        if not B2 is None:
            s += str(B2)
        return s

    def __call__(self, n, watch=False):
        n = Integer(n)
        self._validate(n)
        cmd = 'echo "%s" | %s'%(n, self.__cmd)
        if watch:
            t = tmp_filename()
            os.system('%s | tee %s'%(cmd, t))
            ou = open(t).read()
            os.unlink(t)
        else:
            from subprocess import Popen, PIPE
            x = Popen(cmd, shell=True,
            stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
            x.stdin.close()
            ou = x.stderr.read() + '\n' + x.stdout.read()
        if 'command not found' in ou:
            err = ou + '\n' + 'You must install GMP-ECM.\n'
            err += sage.misc.package.package_mesg('ecm-6.1.3')
            raise RuntimeError, err
        return ou

    def interact(self):
        """
        Interactively interact with the ECM program.
        """
        print "Enter numbers to run ECM on them."
        print "Press control-C to exit."
        os.system(self.__cmd)


    _recommended_B1_list = {15: 2000,
                             20: 11000,
                             25: 50000,
                             30: 250000,
                             35: 1000000,
                             40: 3000000,
                             45: 11000000,
                             50: 44000000,
                             55: 110000000,
                             60: 260000000,
                             65: 850000000,
                             70: 2900000000 }
    """Recommended settings from http://www.mersennewiki.org/index.php/Elliptic_Curve_Method."""

    def __B1_table_value(self, factor_digits, min=15, max=70):
        """Coerces to a key in _recommended_B1_list."""
        if factor_digits < min: factor_digits = min
        if factor_digits > max: raise ValueError('Too many digits to be factored via the elliptic curve method.')
        return 5*ceil(factor_digits/5)

    def recommended_B1(self, factor_digits):
        r"""
        Recommended settings from \url{http://www.mersennewiki.org/index.php/Elliptic_Curve_Method}.
        """
        return self._recommended_B1_list[self.__B1_table_value(factor_digits)]

    @rename_keyword(deprecation=6094, method="algorithm")
    def one_curve(self, n, factor_digits=None, B1=2000, algorithm="ECM", **kwds):
        """
        Run one single ECM (or P-1/P+1) curve on input n.

        INPUT:
            n -- a positive integer
            factor_digits -- decimal digits estimate of the wanted factor
            B1 -- stage 1 bound (default 2000)
            algorithm -- either "ECM" (default), "P-1" or "P+1"
        OUTPUT:
            a list [p,q] where p and q are integers and n = p * q.
            If no factor was found, then p = 1 and q = n.
            WARNING: neither p nor q is guaranteed to be prime.
        EXAMPLES:
            sage: f = ECM()
            sage: n = 508021860739623467191080372196682785441177798407961
            sage: f.one_curve(n, B1=10000, sigma=11)
            [1, 508021860739623467191080372196682785441177798407961]
            sage: f.one_curve(n, B1=10000, sigma=1022170541)
            [79792266297612017, 6366805760909027985741435139224233]
            sage: n = 432132887883903108009802143314445113500016816977037257
            sage: f.one_curve(n, B1=500000, algorithm="P-1")
            [67872792749091946529, 6366805760909027985741435139224233]
            sage: n = 2088352670731726262548647919416588631875815083
            sage: f.one_curve(n, B1=2000, algorithm="P+1", x0=5)
            [328006342451, 6366805760909027985741435139224233]
        """
        n = Integer(n)
        self._validate(n)
        if not factor_digits is None:
            B1 = self.recommended_B1(factor_digits)
        if algorithm == "P-1":
            kwds['pm1'] = ''
        elif algorithm == "P+1":
            kwds['pp1'] = ''
        else:
           if not algorithm == "ECM":
              err = "unexpected algorithm: " + algorithm
              raise ValueError, err
        self.__cmd = self._ECM__startup_cmd(B1, None, kwds)
        child = pexpect.spawn(self.__cmd)
        cleaner.cleaner(child.pid, self.__cmd)
        child.timeout = None
        child.__del__ = nothing   # work around stupid exception ignored error
        child.expect('[ECM]')
        child.sendline(str(n))
        child.sendline("bad") # child.sendeof()
        while True:
            try:
               child.expect('(Using B1=(\d+), B2=(\d+), polynomial ([^,]+), sigma=(\d+)\D)|(Factor found in step \d:\s+(\d+)\D)|(Error - invalid number)')
               info = child.match.groups()
               # B1 is info[1], B2 is info[2], poly is info[3], sigma is info[4],
               # step is info[5], factor is info[6], cofactor is info[7]
               if not info[0] is None:
                  # got Using B1=... line
                  self.last_params = { 'B1' : child.match.groups()[1],
                                       'B2' : child.match.groups()[2],
                                       'poly' : child.match.groups()[3],
                                       'sigma' : child.match.groups()[4] }
               elif info[7] != None:
                  # got Error - invalid number, which means the curve did
                  # end without finding any factor, and the next input 'bad'
                  # was given to GMP-ECM
                  child.kill(0)
                  return [1, n]
               else:
                  # got Factor found...
                  p = Integer(info[6])
                  child.kill(0)
                  return [p, n/p]

            except pexpect.EOF:
               child.kill(0)
               return [1, n]
            child.kill(0)


    def find_factor(self, n, factor_digits=None, B1=2000, **kwds):
        """
        Splits off a single factor of n.
        See ECM.factor()

        OUTPUT:
            list of integers whose product is n

        EXAMPLES:
            sage: f = ECM()
            sage: n = 508021860739623467191080372196682785441177798407961
            sage: f.find_factor(n)
            [79792266297612017, 6366805760909027985741435139224233]

        Note that the input number can't have more than 4095 digits:
            sage: f=2^2^14+1
            sage: ecm.find_factor(f)
            Traceback (most recent call last):
            ...
            ValueError: n must have at most 4095 digits
        """
        n = Integer(n)
        self._validate(n)
        if not 'c' in kwds: kwds['c'] = 1000000000
        if not 'I' in kwds: kwds['I'] = 1
        if not factor_digits is None:
            B1 = self.recommended_B1(factor_digits)
        kwds['one'] = ''
        kwds['cofdec'] = ''
        self.__cmd = self._ECM__startup_cmd(B1, None, kwds)
        self.last_params = { 'B1' : B1 }
        child = pexpect.spawn(self.__cmd)
        cleaner.cleaner(child.pid, self.__cmd)
        child.timeout = None
        child.__del__ = nothing   # program around stupid exception ignored error
        child.expect('[ECM]')
        child.sendline(str(n))
        child.sendline("bad") # child.sendeof()
        while True:

            try:
                child.expect('(Using B1=(\d+), B2=(\d+), polynomial ([^,]+), sigma=(\d+)\D)|(Factor found in step \d:\s+(\d+)\D)|(Error - invalid number)')
                info = child.match.groups()
                if not info[0] is None:
                    self.last_params = { 'B1' : child.match.groups()[1],
                                         'B2' : child.match.groups()[2],
                                         'poly' : child.match.groups()[3],
                                         'sigma' : child.match.groups()[4] }
                elif info[7] != None:
                    child.kill(0)
                    self.primality = [False]
                    return [n]
                else:
                    p = Integer(info[6])
                    child.expect('(input number)|(prime factor)|(composite factor)')
                    if not child.match.groups()[0] is None:
                        child.kill(0)
                        return self.find_factor(n, B1=4+floor(float(B1)/2), **kwds)
                    else:
                        # primality testing is cheap compared to factoring, but has already been done
                        # return [p, n/p]
                        self.primality = [not child.match.groups()[1] is None]
                        child.expect('((prime cofactor)|(Composite cofactor)) (\d+)\D')
                        q = Integer(child.match.groups()[3])
                        self.primality += [not child.match.groups()[1] is None]
                        child.kill(0)
                        return [p, q]


            except pexpect.EOF:
                child.kill(0)
                self.primality = [False]
                return [n]
        child.kill(0)


    def factor(self, n, factor_digits=None, B1=2000, **kwds):
        """
        Returns a list of integers whose product is n, computed using
        GMP-ECM, and PARI for small factors.

        ** WARNING: There is no guarantee that the factors returned are
        prime. **

        INPUT:
            n -- a positive integer
            factor_digits -- optional guess at how many digits are in the smallest factor.
            B1 -- initial lower bound, defaults to 2000 (15 digit factors)
            kwds -- arguments to pass to ecm-gmp. See help for ECM for more details.

        OUTPUT:
            a list of integers whose product is n

        NOTE:
            Trial division should typically be performed before using
            this method.  Also, if you suspect that n is the product
            of two similarly-sized primes, other methods (such as a
            quadratic sieve -- use the qsieve command) will usually be
            faster.

        EXAMPLES:
            sage: ecm.factor(602400691612422154516282778947806249229526581)
            [45949729863572179, 13109994191499930367061460439]

            sage: ecm.factor((2^197 + 1)/3)           # takes a long time
            [197002597249, 1348959352853811313, 251951573867253012259144010843]
        """
        if B1 < 2000 or len(str(n)) < 15:
            return sum([[p]*e for p, e in Integer(n).factor()], [])

        factors = self.find_factor(n, factor_digits, B1, **kwds)
        factors.sort()
        if len(factors) == 1:
            return factors
        assert len(factors) == 2
        _primality = [self.primality[0], self.primality[1]]
        try:
            last_B1 = self.last_params['B1']
        except AttributeError:
            self.last_params = {}
            self.last_params['B1'] = 10
            last_B1 = 10
        if not _primality[1]:
            factors[1:2] = self.factor(factors[1], B1=last_B1, **kwds)
            _primality[1:2] = self.primality
        if not _primality[0]:
            factors[0:1] = self.factor(factors[0], B1=last_B1, **kwds)
            _primality[0:1] = self.primality
        self.primality = _primality
        factors.sort()
        return factors


    def get_last_params(self):
        """
        Returns the parameters (including the curve) of the last ecm run.
        In the case that the number was factored successfully, this will return the parameters that yielded the factorization.

        INPUT:
            none

        OUTPUT:
            The parameters for the most recent factorization.

        EXAMPLES:
            sage: ecm.factor((2^197 + 1)/3)             # long time
            [197002597249, 1348959352853811313, 251951573867253012259144010843]
            sage: ecm.get_last_params()                 # random output
            {'poly': 'x^1', 'sigma': '1785694449', 'B1': '8885', 'B2': '1002846'}

        """
        return self.last_params



    def time(self, n, factor_digits, verbose=0):
        """
        Gives an approximation for the amount of time it will take to find a factor
        of size factor_digits in a single process on the current computer.
        This estimate is provided by GMP-ECM's verbose option on a single run of a curve.

        INPUT:
            n -- a positive integer
            factor_digits -- the (estimated) number of digits of the smallest factor

        EXAMPLES::

            sage: n = next_prime(11^23)*next_prime(11^37)

            sage: ecm.time(n, 20)                  # not tested
            Expected curves: 77     Expected time: 7.21s
            sage: ecm.time(n, 25)                  # not tested
            Expected curves: 206    Expected time: 1.56m
            sage: ecm.time(n, 30, verbose=1)       # not tested
            GMP-ECM 6.1.3 [powered by GMP 4.2.1] [ECM]

            Input number is 304481639541418099574459496544854621998616257489887231115912293 (63 digits)
            Using MODMULN
            Using B1=250000, B2=128992510, polynomial Dickson(3), sigma=2307628716
            dF=2048, k=3, d=19110, d2=11, i0=3
            Expected number of curves to find a factor of n digits:
            20      25      30      35      40      45      50      55      60      65
            8       50      430     4914    70293   1214949 2.5e+07 5.9e+08 1.6e+10 2.7e+13
            Step 1 took 6408ms
            Using 16 small primes for NTT
            Estimated memory usage: 3862K
            Initializing tables of differences for F took 16ms
            Computing roots of F took 128ms
            Building F from its roots took 408ms
            Computing 1/F took 608ms
            Initializing table of differences for G took 12ms
            Computing roots of G took 120ms
            Building G from its roots took 404ms
            Computing roots of G took 120ms
            Building G from its roots took 412ms
            Computing G * H took 328ms
            Reducing  G * H mod F took 348ms
            Computing roots of G took 120ms
            Building G from its roots took 408ms
            Computing G * H took 328ms
            Reducing  G * H mod F took 348ms
            Computing polyeval(F,G) took 1128ms
            Step 2 took 5260ms
            Expected time to find a factor of n digits:
            20      25      30      35      40      45      50      55      60      65
            1.58m   9.64m   1.39h   15.93h  9.49d   164.07d 9.16y   218.68y 5825y   1e+07y
            Expected curves: 4914   Expected time: 1.39h

        """
        self._validate(n)
        B1 = self.recommended_B1(factor_digits)
        self.__cmd = self._ECM__startup_cmd(B1, None, {'v': ' '})
        child = pexpect.spawn(self.__cmd)
        cleaner.cleaner(child.pid, self.__cmd)
        child.timeout = None
        child.expect('[ECM]')
        child.sendline(str(n))
        try:
            child.sendeof()
        except Exception:
            pass
        child.expect('20\s+25\s+30\s+35\s+40\s+45\s+50\s+55\s+60\s+65')
        if verbose:
            print child.before,
            print child.after,
        child.expect('(\d\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s', timeout=None)
        offset = (self.__B1_table_value(factor_digits, 20, 65)-20)/5
        curve_count = child.match.groups()[int(offset)]
        if verbose:
            print child.before,
            print child.after,
        child.expect('20\s+25\s+30\s+35\s+40\s+45\s+50\s+55\s+60\s+65', timeout=None)
        if verbose:
            print child.before,
            print child.after,
        child.expect('(\d\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s', timeout=None)
        if verbose:
            print child.before,
            print child.after
        time = child.match.groups()[int(offset)]
        child.kill(0)
        print "Expected curves:", curve_count, "\tExpected time:", time

    def _validate(self, n):
        """
        Verify that n is positive and has at most 4095 digits.

        INPUT:
            n

        This function raises a ValueError if the two conditions listed above
        are not both satisfied.  It is here because GMP-ECM silently ignores
        all digits of input after the 4095th!

        EXAMPLES:
            sage: ecm = ECM()
            sage: ecm._validate(0)
            Traceback (most recent call last):
            ...
            ValueError: n must be positive
            sage: ecm._validate(10^5000)
            Traceback (most recent call last):
            ...
            ValueError: n must have at most 4095 digits
        """
        n = Integer(n)
        if n <= 0:
            raise ValueError, "n must be positive"
        if n.ndigits() > 4095:
            raise ValueError, "n must have at most 4095 digits"


# unique instance
ecm = ECM()
