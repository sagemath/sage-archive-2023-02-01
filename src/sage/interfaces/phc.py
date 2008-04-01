r"""
Interface to PHC.

PHC computes numerical information about systems of polynomials over
the complex numbers.

PHC implements polynomial homotopy methods to exploit structure in
order to better approximate all isolated solutions. The package also
includes extra tools to handle positive dimensional solution
components.

AUTHOR:
   -- PHC was written by J. Verschelde, R. Cools, and many others (?)
   -- William Stein and Kelly ?? -- first version of interface to PHC
   -- Marshall Hampton -- second version of interface to PHC

TODO:
    This either needs to be able to handle arbitrary PHCpack command
    options, or have enough special purpose functions added to get the
    full functionality.  One important missing special case is having
    a root-counting function that could do mixed volumes and
    structured Bezout counts.  Another nice feature to have would be a
    graphical output of the path-tracking, as has been done for Maple.
    Currently Marshall Hampton and Alex Jokela are working on this, please contact us if you are interested.
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2007 Marshall Hampton <hamptonio@gmail.com>
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
from sage.rings.complex_field import ComplexField
from sage.rings.all import CC
from sage.rings.integer import Integer

import pexpect

def get_solution_dicts(output_file_contents, input_ring, get_failures = True):
    '''
    Returns a list of dictionaries of variable:value (key:value)
    pairs.  Only used internally; see the solution_dict function in
    the PHC_Object class definition for details.
    '''
    output_list = output_file_contents.splitlines()
    test = 'False'
    solution_dicts = []
    for solution_line in range(len(output_list)-1,-1,-1):
        if output_list[solution_line].find('THE SOLUTIONS') == 0:
            break
    var_number = int(output_list[solution_line+2].split(' ')[1])
    sol_number = int(output_list[solution_line+2].split(' ')[0])
    for i in range(solution_line + 1,len(output_list)):
        if output_list[i].count('the solution for t') == 1:
            if output_list[i-3].count('success') > 0 or get_failures == True:
                temp_dict = {}
                for j in range(1,var_number+1):
                    rawsplit = output_list[i+j].split(': ')[1].split(' ')
                    for extras in range(rawsplit.count('')):
                        rawsplit.remove('')
                    temp_var = output_list[i+j].split(': ')[0].replace(' ','')
                    temp_dict[input_ring(temp_var)] = CC(rawsplit[0],rawsplit[1])
                solution_dicts.append(temp_dict)
    return solution_dicts

def get_variable_list(output_file_contents):
    '''
    Returns the variables, as strings, in the order in which PHCpack has processed them.
    '''
    output_list = output_file_contents.splitlines()
    for solution_line in range(len(output_list)-1,-1,-1):
        if output_list[solution_line].find('THE SOLUTIONS') == 0:
            break
    var_number = int(output_list[solution_line+2].split(' ')[1])
    varlist = []
    for var_ind in range(var_number):
        var = output_list[solution_line + 8 + var_ind].split(' ')[1]
        varlist.append(var)
    return varlist

class PHC_Object:
    """
    A container for data from the PHCpack program - lists of float
    solutions, etc.
    """
    def __init__(self, output_file_contents, input_ring):
	'''
	INPUT:
	    output_file_contents: the string output of PHCpack
	    input_ring: for coercion of the variables into the desired ring.
	'''
        self.output_file_contents = output_file_contents
	self.input_ring = input_ring

    def solution_dicts(self, get_failures = False):
        """
        Returns a list of solutions in dictionary form: variable:value.

	INPUT:
            self -- for access to self_out_file_contents, the string
	    of raw PHCpack output.

            get_failures (optional) -- a boolean.  The default (False)
	    is to not process failed homotopies.  These either lie on
	    positive-dimensional components or at infinity.

        OUTPUT:

	    solution_dicts: a list of dictionaries.  Each dictionary
	    element is of the form variable:value, where the variable
	    is an element of the input_ring, and the value is in
	    ComplexField.

        """
        try:
            return self.__solution_dicts
        except AttributeError:
            pass
        solution_dicts = get_solution_dicts(self.output_file_contents, self.input_ring, get_failures = get_failures)
        self.__solution_dicts = solution_dicts
        return solution_dicts

    def solutions(self, get_failures = False):
	"""
	Returns a list of solutions in the ComplexField.  Use the variable_list function to get the order of variables used by PHCpack, which is usually different than the term order of the input_ring.

	INPUT:
	    self -- for access to self_out_file_contents, the string
	    of raw PHCpack output.
	    get_failures (optional) -- a boolean.  The default (False)
	    is to not process failed homotopies.  These either lie on
	    positive-dimensional components or at infinity.

        OUTPUT:
	    solutions: a list of lists of ComplexField-valued solutions.
	"""
	try:
	    return self.__solutions
	except AttributeError:
            pass
        solution_dicts = get_solution_dicts(self.output_file_contents, self.input_ring, get_failures = get_failures)
        self.__solution_dicts = solution_dicts
	solutions = [sol_dict.values() for sol_dict in solution_dicts]
	self.__solutions = solutions
        return solutions

    def variable_list(self):
        """
	Returns the variables, as strings, in the order in which
	PHCpack has processed them.
        """
        try:
            return self.__var_list
        except AttributeError:
            pass
        var_list = get_variable_list(self.output_file_contents)
        self.__var_list = var_list
        return var_list

class PHC:
    """
    A class to interface with PHCpack, for computing numerical
    homotopies and root counts.

    EXAMPLES:
        sage: from sage.interfaces.phc import phc
        sage: R.<x,y> = PolynomialRing(CDF,2)
	sage: testsys = [x^2 + 1, x*y - 1]
        sage: phc.mixed_volume(testsys)        # optional -- you must have phc install
        2
        sage: v = phc.blackbox(testsys, R)     # optional
        sage: sols = v.solutions()             # optional
        sage: sols.sort()                      # optional
	sage: sols                             # optional
        [[-1.00000000000000*I, 1.00000000000000*I], [1.00000000000000*I, -1.00000000000000*I]]
	sage: sol_dict = v.solution_dicts()    # optional
        sage: sols_from_dict = [sorted(d.items()) for d in sol_dict] # optional
        sage: sorted(sols_from_dict)           # optional
        [[(y, -1.00000000000000*I), (x, 1.00000000000000*I)], [(y, 1.00000000000000*I), (x, -1.00000000000000*I)]]
	sage: residuals = [[test_equation.change_ring(CDF).subs(sol) for test_equation in testsys] for sol in v.solution_dicts()]      # optional
	sage: residuals                             # optional
	[[0, 0], [0, 0]]
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
        #if len(polys) == 0:
        #    ngens = 1
        #else:
        #    ngens = polys[0].parent().ngens()
        #s += '%s\n'%ngens

        # TODO: actually implement the variable mapping (copy code from gfan?)
        # For now, we'll assume all variables satisfy the restrictive
        # PHC constraints, i.e., length at most 3.
        for f in polys:
            s += '%s;\n'%f           # note the semicolon *terminators*

        # Current just return "none" for the mapping.
        return s, None


    def mixed_volume(self, polys, verbose=False):
        """
        Computes the mixed volume of the polynomial system given by the input polys.

        INPUT:
            polys -- a list of multivariate polynomials (elements of a multivariate
                     polynomial ring).
            verbose -- print lots of verbose information about what this function does.

        OUTPUT:
            The mixed volume.

        EXAMPLES:
            sage: from sage.interfaces.phc import phc
            sage: R.<x,y> = PolynomialRing(CDF,2)
            sage: testsys = [x^2 + 1, x*y - 1]
            sage: phc.mixed_volume(testsys)        # optional -- you must have phc install
            2

        """
        # Get three temporary file names (these will be in SAGE_HOME/.sage/tmp/pid)
        input_filename = sage.misc.misc.tmp_filename()
        output_filename = sage.misc.misc.tmp_filename()
        log_filename = sage.misc.misc.tmp_filename()

        # Get the input polynomial text
        input, mapping = self._input_file(polys)
        if verbose:
            print "Writing the input file to %s"%input_filename
        open(input_filename,'w').write(input)

        if verbose:
            print "The following file will be the input polynomial file to phc."
            print input

        # Create a phc process
        child_phc = pexpect.spawn('phc -m')
        child_phc.sendline('y')
        child_phc.sendline(input_filename)
        child_phc.sendline(output_filename)
        child_phc.sendline('4')
        child_phc.sendline('n')
        child_phc.sendline('n')
        child_phc.sendline('n')
        child_phc.expect('results')
        dots = child_phc.read()
        if verbose:
            print "should be ... : " + dots
        child_phc.close()
        if not os.path.exists(output_filename):
            raise RuntimeError, "The output file does not exist; something went wrong running phc."

        # Read the output produced by PHC
        out = open(output_filename).read()

        # Delete the temporary files
        os.unlink(input_filename)
        os.unlink(output_filename)

        # All done
        out_lines = out.split('\n')
        for a_line in out_lines:
            # the two conditions below are necessary because of changes in output format
            if a_line.find('The mixed volume equals :') == 0 or a_line.find('common mixed volume :') == 0:
                if verbose: print 'found line: ' +  a_line
                mixed_vol = Integer(a_line.split(':')[1])
                break

        try:
            return mixed_vol
        except NameError:
            raise RuntimeError, "Mixed volume not found in output; something went wrong running phc."


    def blackbox(self, polys, input_ring, verbose=False):
        """
        Returns as a string the result of running PHC with the given polynomials
        under blackbox mode (the '-b' option).

        INPUT:
            polys -- a list of multivariate polynomials (elements of a multivariate
                     polynomial ring).
            input_ring: for coercion of the variables into the desired ring.
            verbose -- print lots of verbose information about what this function does.

        OUTPUT:
            a PHC_Object object containing the phcpack output string.
        """

        # Get three temporary file names (these will be in SAGE_HOME/.sage/tmp/pid)
        input_filename = sage.misc.misc.tmp_filename()
        output_filename = sage.misc.misc.tmp_filename()
        log_filename = sage.misc.misc.tmp_filename()

        # Get the input polynomial text
        input, mapping = self._input_file(polys)
        if verbose:
            print "Writing the input file to %s"%input_filename
        open(input_filename,'w').write(input)

        if verbose:
            print "The following file will be the input polynomial file to phc."
            print input

        # Create the phc command line>
        cmd = 'phc -b %s %s'%(input_filename, output_filename)

        if verbose:
            print "The phc command line is:"
            print cmd

        # Do it -- make the system call.
        e = os.system(cmd)

        # Was there an error?
        if e:
            if os.system('which phc'):
                print os.system('which phc') + '  PHC needs to be installed and in your path'
                raise RuntimeError
            # todo -- why? etc.
            raise RuntimeError, open(log_filename).read() + "\nError running phc."

        if not os.path.exists(output_filename):
            raise RuntimeError, "The output file does not exist; something went wrong running phc."

        # Read the output produced by PHC
        out = open(output_filename).read()

        # Delete the temporary files
        os.unlink(input_filename)
        os.unlink(output_filename)

        # All done
        return PHC_Object(out, input_ring)


################################

# The unique phc interface instance.
phc = PHC()

