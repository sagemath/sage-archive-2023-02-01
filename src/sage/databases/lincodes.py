"""nodoctest
Linear codes
"""

#*****************************************************************************
#
#      Sage: Copyright (C) 2005 William Stein <was@math.harvard.edu>
#                          2005 Steven Sivek
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

import os, re, urllib
from sage.rings.all import Integer

def parse_bound_html(text, n, k):
    lower, upper = -1, -1
    bound = re.compile(r'(?P<type>[LU])b\('
                       r'(?P<n>[0-9]*),(?P<k>[0-9]*)'
                       r'\) = (?P<val>[0-9]*) ');

    for line in re.split(r'[\s,]*\n', text):
        m = bound.search(line)
        if m and int(m.group('n')) == n and int(m.group('k')) == k:
            if m.group('type') == 'L':
                lower = int(m.group('val'))
            elif m.group('type') == 'U':
                upper = int(m.group('val'))
    return lower, upper

## def linear_code_bound(q, n, k, verbose=False):
##     r"""
##     Find bounds on the minimum distance of linear codes over GF(q)
##     with length n and dimension k, courtesy of
##          \url{http://www.win.tue.nl/~aeb/voorlincod.html}.
##     If no bounds are in the database, returns lower and upper
##     bounds of -1.

##     INPUT:
##         q -- integer, the field order, which must be in
##                       [2, 3, 4, 5, 7, 8, 9]
##         n -- integer, the length of the code
##         k -- integer, the dimension of the code
##         verbose -- bool (default=False), print verbose message

##     OUTPUT:
##         integer -- lower bound
##         integer -- upper bound
##         str -- text about why the bounds are as given

##     EXAMPLES:
##     To find lower and upper bounds for values q=7, n=32, k=8, type
##         sage: lower, upper, text = linear_code_bound(7, 32, 8)     # optional -- needs internet
##         sage: lower                   # optional
##         19
##         sage: upper                   # optional
##         21
##         sage: text                    # optional
##         'Lb(32,8) = 19 DG4\n\nUb(32,8) = 21 follows by a one-step Griesmer bound from:\nUb(10,7) = 3 is found by considering shortening to:\nUb(9,6) = 3 is found by construction B:\n[consider deleting the (at most) 6 coordinates of a word in the dual]'

##     When bounds are not known the upper and lower returned bounds are -1:
##         sage: linear_code_bound(9, 32, 200)          # optional -- needs internet
##         (-1, -1, '(error executing -why-)')

##     This function raises an IOError if an error occurs downloading
##     data or parsing it.  It raises a ValueError if the q input is
##     invalid.

##     AUTHOR: 2005-11-14: Steven Sivek
##     """
##     q = Integer(q)
##     if not q in [2, 3, 4, 5, 7, 8, 9]:
##         raise ValueError, "q (=%s) must be in [2,3,4,5,7,8,9]"%q
##     n = Integer(n)
##     k = Integer(k)

##     param = ("%s/%s/%s"%(q,n,k)).replace('L','')

##     url = "http://homepages.cwi.nl/htbin/aeb/lincodbd/"+param
##     if verbose:
##         print "Looking up the bounds at %s"%url
##     f = urllib.urlopen(url)
##     s = f.read()
##     f.close()

##     i = s.find("<PRE>")
##     j = s.find("</PRE>")
##     if i == -1 or j == -1:
##         raise IOError, "Error parsing data (missing pre tags)."
##     text = s[i+5:j].strip()
##     lower, upper = parse_bound_html(text, n, k)

##     return lower, upper, text

