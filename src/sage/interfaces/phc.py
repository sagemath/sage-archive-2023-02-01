r"""
Interface to PHC, which solves systems of polynomials over the complex
numbers.

PHC implements polynomial homotopy methods to exploit structure in
order to better approximate all isolated solutions. The package also
includes extra tools to handle positive dimensional solution
components.

AUTHOR:
   -- PHC was written by [[todo]]
   -- Interface written by [[todo]]

TOD:

   ...

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

import os

import sage.misc.misc

class PHC:
    """
    ... say what it is...
    """
    def __call__(self, input, cmd_list):
        input_file = sage.misc.misc.tmp_filename()
        output_file = sage.misc.misc.tmp_filename()
        cmd_file = sage.misc.misc.tmp_filename()
        open(input_file,'w').write(input)
        # more work on cmd_list.
        cmd_list[2] = input_file
        cmd_list[3] = output_file

        s = '\n'.join([str(x) for x in cmd_list])
        print s
        open(cmd_file,'w').write(s + '\n\n')
        e = os.system('phc < %s'%cmd_file)
        if e:
            # todo -- why? etc.
            raise RuntimeError, "Error running phc..."

        out = open(output_file).read()
        os.unlink(input_file)
        os.unlink(output_file)
        os.unlink(cmd_file)
        return out


# The instance
phc = PHC()

