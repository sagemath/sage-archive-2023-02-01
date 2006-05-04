r"""
Interface to GMP-ECM

"""

import os, pexpect

import sage.rings.integer
from sage.misc.all import verbose, get_verbose, tmp_filename

import sage.misc.package


class ECM:
    def __init__(self, B1=10, B2=None, **kwds):
        r"""
        Create an interface to the GMP-ECM elliptic curve method
        factorization program.

        See \url{http://www.komite.net/laurent/soft/ecm/}.

        AUTHORS:
            These people wrote GMP-ECM:
              Jim Fougeron, Laurent Fousse,
              Alexander Kruppa, Dave Newman
              Paul Zimmermann

        William Stein -- wrote the SAGE interface to GMP-ECM

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
        n = sage.rings.integer.Integer(n)
        cmd = 'echo "%s" | %s'%(n, self.__cmd)
        if watch:
            t = tmp_filename()
            os.system('%s | tee %s'%(cmd, t))
            ou = open(t).read()
            os.unlink(t)
        else:
            i,o,e = os.popen3(cmd)
            i.close()
            ou = e.read() + '\n' + o.read()
        if 'command not found' in ou:
            err = ou + '\n' + 'You must install GMP-ECM.\n'
            err += sage.misc.package.package_mesg('ecm-6.0.1')
            raise RuntimeError, err
        return ou

    def interact(self):
        """
        Interactively interact with the ECM program.
        """
        print "Enter numbers to run ECM on them."
        print "Press control-C to exit."
        os.system(self.__cmd)



