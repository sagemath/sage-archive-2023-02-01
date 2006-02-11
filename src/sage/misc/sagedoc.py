"""nodoctest
Format SAGE documentation for viewing with IPython
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

substitutes = [('\\_','_'),\
               ('\\to', '-->'), \
               ('<BLANKLINE>',''), \
               ('\\leq', '<='), ('\\geq', '>='), \
               ('\\sage', 'SAGE'), \
               ('\\SAGE', 'SAGE'), \
               ('\\rm', ''), \
               ('\\cdot', ' *'), \
               ('$',''), ('\\',''), ('sage.:', 'sage:'), ('backslash','\\'), \
               ('begin{enumerate}',''), ('end{enumerate}',''), \
               ('begin{verbatim}',''), ('end{verbatim}',''), \
               ('mapsto', ' |--> '), \
               ('ldots', '...'), ('note{','NOTE: ')]


def _rmcmd(s, cmd, left='', right=''):
    c = '\\%s{'%cmd
    while True:
        i = s.find(c)
        if i == -1:
            return s
        nesting = 1
        j = i+len(c)+1
        while j < len(s) and nesting > 0:
            if s[j] == '{':
                nesting += 1
            elif s[j] == '}':
                nesting -= 1
            j += 1
        j -= 1  # j is position of closing '}'
        if j < len(s):
            s = s[:i] + left + s[i+len(c):j] + right + s[j+1:]
        else:
            return s

# I wanted to be cool and use regexp's, but they aren't really
# useful, since really this is a parsing problem, because of
# nesting of commands, etc.   Since it doesn't have to be
# super super fast (it's a page of text scrolled to the user),
# the above works fine.

#
## import re
## def _rmcmd(s, cmd, left='', right=''):
##     c = '\\%s{.*}'%cmd
##     r = re.compile(c, re.DOTALL)
##     while True:
##         m = r.search(s)
##         if m is None: break
##         s = s[:m.start()] + left + s[m.start()+len(cmd)+1:m.end()-1] \
##             + right + s[m.end():]
##     return s

def format(s):
    """
    Format SAGE documentation for viewing with IPython.
    """
    if not isinstance(s, str):
        raise TypeError, "s must be a string"
    s = _rmcmd(s, 'url')
    s = _rmcmd(s, 'code')
    s = _rmcmd(s, 'mbox')
    s = _rmcmd(s, 'note', 'NOTE: ', '')
    s = _rmcmd(s, 'emph', '*', '*')
    for a,b in substitutes:
        s = s.replace(a,b)
    return s
