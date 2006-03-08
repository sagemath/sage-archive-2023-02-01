"""
Pager for showing strings
"""

#############################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#############################################################################

# Currently we just use the Ipython pager.  If we want to use
# something else, we can just change this function.  Any code
# in sage that uses a pager should use this pager.

import IPython.genutils

pager = IPython.genutils.page
