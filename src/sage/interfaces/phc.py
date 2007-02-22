r"""nodoctest
Interface to PHC.

WARNING -- currently totally broken!! -- do not use.

PHC computes numerical information about systems of polynomials over
the complex numbers.

PHC implements polynomial homotopy methods to exploit structure in
order to better approximate all isolated solutions. The package also
includes extra tools to handle positive dimensional solution
components.

AUTHOR:
   -- PHC was written by J. Verschelde, R. Cools, and many others (?)
   -- William Stein and Kelly ?? -- first version of interface to PHC
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       TODO -- add yourself as copyright holder.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
########################################################################

import copy, os

import sage.misc.misc

class PHC:
    """
    ... say what it is...

    EXAMPLES:
        sage: from sage.interfaces.phc import phc
        sage: R.<x,y> = PolynomialRing(QQ,2)
        sage: v = phc([x^2 + 1, x*y - 1])          # optional -- you must have phc install
        sage.: print v
         2
        x^2+ 1;
        x*y-1;
        <BLANKLINE>
        THE SOLUTIONS :
        <BLANKLINE>
        2 2
        ===========================================================================
        solution 1 :    start residual :  1.225E-16   #iterations : 1   success
        t :  0.00000000000000E+00   0.00000000000000E+00
        m : 1
        the solution for t :
         x :  0.00000000000000E+00  -1.00000000000000E+00
         y :  0.00000000000000E+00   1.00000000000000E+00
        == err :  1.225E-16 = rco :  3.333E-01 = res :  0.000E+00 = complex regular ==
        solution 2 :    start residual :  1.225E-16   #iterations : 1   success
        t :  0.00000000000000E+00   0.00000000000000E+00
        m : 1
        the solution for t :
         x :  0.00000000000000E+00   1.00000000000000E+00
         y :  0.00000000000000E+00  -1.00000000000000E+00
        == err :  1.225E-16 = rco :  3.333E-01 = res :  0.000E+00 = complex regular ==
        ===========================================================================
        A list of 2 solutions has been refined :
        Number of regular solutions   : 2.
        Number of singular solutions  : 0.
        Number of real solutions      : 0.
        Number of complex solutions   : 2.
        Number of clustered solutions : 0.
        Number of failures            : 0.
        ===========================================================================
        Frequency tables for correction, residual, condition, and distances :
        FreqCorr :  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 : 2
        FreqResi :  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 : 2
        FreqCond :  2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 : 2
        FreqDist :  2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 : 2
        Small correction terms and residuals counted to the right.
        Well conditioned and distinct roots counted to the left.
        <BLANKLINE>
        TIMING INFORMATION for Solving the polynomial system
        The elapsed time in seconds was                  0.016529000 =  0h 0m 0s 17ms
        User time in seconds was                         0.014512000 =  0h 0m 0s 15ms
        System CPU time in seconds was                   0.002017000 =  0h 0m 0s  2ms
        Non-I/O page faults was                          0
        I/O page faults was                              0
        Signals delivered was                            0
        Swaps was                                        0
        Total context switches was                       0
        <BLANKLINE>
          ---------------------------------------------------------------------
          |                    TIMING INFORMATION SUMMARY                     |
          ---------------------------------------------------------------------
          |   root counts  |  start system  |  continuation  |   total time   |
          ---------------------------------------------------------------------
          |  0h 0m 0s  0ms |  0h 0m 0s  0ms |  0h 0m 0s  0ms |  0h 0m 0s 15ms |
          ---------------------------------------------------------------------
        PHC ran from 24 October 2006, 07:33:07 till 24 October 2006, 07:33:07.
        The total elapsed time is 0 seconds.
        <BLANKLINE>
    """
    def _input_file(self, polys):
        """
        This is used internally to implement the PHC interface.

        INPUT:
            polys -- a list of polynomials in a SAGE polynomial ring
                     over a field that embeds into the complex
                     numbers.

        OUTPUT:
            -- a PHC input file (as a text string) that
                describes these polynomials.
            -- a mapping between variable names.

        PHC has the following constraint on polynomials:
            An unknown can be denoted by at most 3 characters.  The
            first character must be a letter and the other two
            characters must be different from '+', '-', '*', '^', '/',
            ';', '(' and ')'.  The letter i means sqrt(-1), whence it
            does not represent an unknown.  The number of unknowns may
            not exceed the declared dimension.

        Thus just like for the Gfan interface, it is necessary to map
        the SAGE variables (which can be arbitrarily long) to
        variables that satisfy the above constraints.


        IMPORTANT NOTE -- I just tried a 4-character variable name, and it
        worked fine -- they must have improved their program.
        """
        if not isinstance(polys, (list, tuple)):
            raise TypeError, 'polys must be a list or tuple'
        s = '%s\n'%len(polys)
        if len(polys) == 0:
            ngens = 1
        else:
            ngens = polys[0].parent().ngens()
        s += '%s\n'%ngens

        # TODO: actually implement the variable mapping (copy code from gfan?)
        # For now, we'll assume all variables satisfy the restrictive
        # PHC constraints, i.e., length at most 3.
        for f in polys:
            s += '%s;\n'%f           # note the semicolon *terminators*

        # Current just return "none" for the mapping.
        return s, None

    def __call__(self, polys, mode='b',
                 cmd_list=['y', 'input', 'output'], verbose=False):
        """
        Returns as a string the result of running PHC with the given polynomials,
        mode, and command list.  The command list is a sequence of commands that
        would get input interactively to PHC if you ran it from the command line.

        INPUT:
            polys -- a list of multivariate polynomials (elements of a multivariate
                     polynomial ring).
            mode -- (default: 'b' for batch mode) one of the the letters
                    '0', 'a', b', 'c', 'd', 'e', 'f',
                    'k', 'l', 'm', 'p', 'q', 'r', 's', 'v', 'w', 'z'
            cmd_list -- a list of commands, which will get input to PHC.
                    Use the special commands 'input' and 'output' to represent
                    the tempory input filename and output filename, respectively.
                    default: ['y', 'input', 'output']
            verbose -- print lots of verbose information about what this function does.

        OUTPUT:
            a string.
        """

        # Get three temporary file names (these will be in SAGE_HOME/.sage/tmp/pid)
        input_filename = sage.misc.misc.tmp_filename()
        output_filename = sage.misc.misc.tmp_filename()
        cmd_filename = sage.misc.misc.tmp_filename()
        log_filename = sage.misc.misc.tmp_filename()

        # Get the input polynomial text
        input, mapping = self._input_file(polys)
        if verbose:
            print "Writing the input file to %s"%input_filename
        open(input_filename,'w').write(input)

        if verbose:
            print "The following file will be the input polynomial file to phc."
            print input

        # Replace 'input' and 'output' in the command list by the
        # temporary file names.
        # CRUCIAL: We *must* make a copy of cmd_list, otherwise we would
        # be modifying it in place, which is very bad.
        cmd_list = copy.copy(cmd_list)
        try:
            cmd_list[cmd_list.index('input')] =  input_filename
        except ValueError:
            pass
        try:
            cmd_list[cmd_list.index('output')] = output_filename
        except ValueError:
            pass

        # Create a single carriage-return separated string from the command list
        cmd_in = '\n'.join([str(x) for x in cmd_list])

        if verbose:
            print "The following command list will be fed to phc from stdin:"
            print cmd_in

        # Write the input command list to a temporary file.
        open(cmd_filename,'w').write(cmd_in + '\n\n')

        # Create the phc command line>
        cmd = 'phc -%s < %s > %s'%(mode, cmd_filename, log_filename)

        if verbose:
            print "The phc command line is:"
            print cmd

        # Do it -- make the system call.
        e = os.system(cmd)

        # Was there an error?
        if e:
            if os.system('which phc'):
                raise RuntimeError, "\n\n" + '-'*70 + '\n Install the *free* PHCpack somewhere on your computer\n' + \
                      " so that the command 'phc' is in your PATH.\n" + \
                      " Get it at http://www.cs.kuleuven.ac.be/~nines/research/phcpack/Download/index.html" + '\n' + '-'*70
            # todo -- why? etc.
            raise RuntimeError, open(log_filename).read() + "\nError running phc."

        if not os.path.exists(output_filename):
            raise RuntimeError, "The output file does not exist; something went wrong running phc."

        # Read the output produced by PHC
        out = open(output_filename).read()

        # Delete the temporary files
        os.unlink(input_filename)
        os.unlink(output_filename)
        os.unlink(cmd_filename)
        os.unlink(log_filename)

        # All done
        return out


################################

# The unique phc interface instance.
phc = PHC()

