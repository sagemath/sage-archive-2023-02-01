r"""
Interface to Groebner Fan

AUTHOR:
   -- Anders Nedergaard Jensen: Write gfan C++ program, which implements
      algorithms many of which were invented by Jensen, Komei
      Fukuda, and Rekha Thomas.
   -- William Stein (2006-03-18): wrote gfan interface (first version)

TODO -- not done; this is just the first few minutes of work
        on an interface!

   * at most 52 variables:
       use gfan_substitute to make easier (?)

   * --symmetry is really useful
            - permutations are 0-based *not* cycle notation; a <---> 0
     output is broken up much more nicely.

   * -- can work in Z/pZ for p <= 32749

   * -- can compute individual GB's for lex and revlex (via buchberger)
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

class Gfan:
    """
    Interface to Anders Jensen's Groebner Fan program.
    """
    def __call__(self, I, cmd='', verbose=False, format=True):
        if cmd != '' and cmd.lstrip()[0] != '-':
            cmd = 'gfan_%s'%cmd
        else:
            cmd = 'gfan'

        if verbose:
            print "gfan command:\n%s"%cmd
            print "gfan input:\n%s"%I

        if not verbose:
            cmd += ' 2>/dev/null'
        i, o = os.popen2(cmd)
        s = '{' + (str(I).replace(' ', '').replace("'",""))[1:-1] + '}'
        i.write(s)
        i.close()
        t = o.read()
        if format:
            t = t.replace('\n','').replace(',',' ')
        return t



# The instance
gfan = Gfan()

