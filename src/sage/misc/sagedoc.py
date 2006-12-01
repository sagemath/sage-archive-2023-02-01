"""nodoctest
Format SAGE documentation for viewing with IPython
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

substitutes = [('\\item', '*'), \
               ('\\_','_'),\
               ('\\to', '-->'), \
               ('<BLANKLINE>',''), \
               ('\\leq', '<='), \
               ('\\geq', '>='), \
               ('\\le', '<='), \
               ('\\ge', '>='), \
               ('\\bf', ''),\
               ('\\sage', 'SAGE'), \
               ('\\SAGE', 'SAGE'), \
               ('\\rm', ''), \
               ('cdots', '...'), \
               ('\\cdot', ' *'), \
               ('$',''), ('\\',''), ('sage.:', 'sage:'), ('backslash','\\'), \
               ('begin{enumerate}',''), ('end{enumerate}',''), \
               ('begin{itemize}',''), ('end{itemize}',''), \
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
    import sage.all
    import sage.server.support
    docs = set([])
    while True:
        i = s.find("<<<")
        if i == -1: break
        j = s[i+3:].find('>>>')
        if j == -1: break
        obj = s[i+3:i+3+j]
        if obj in docs:
            t = ''
        else:
            x = eval('sage.all.%s'%obj, locals())
            t0 = sage.server.support.get_def(x, obj)
            t1 = my_getdoc(x)
            t = 'Definition: ' + t0 + '\n\n' + t1
            docs.add(obj)
        s = s[:i] + '\n' + t + s[i+6+j:]

    s = _rmcmd(s, 'url')
    s = _rmcmd(s, 'code')
    s = _rmcmd(s, 'mbox')
    s = _rmcmd(s, 'text')
    s = _rmcmd(s, 'section')
    s = _rmcmd(s, 'subsection')
    s = _rmcmd(s, 'subsubsection')
    s = _rmcmd(s, 'note', 'NOTE: ', '')
    s = _rmcmd(s, 'emph', '*', '*')
    for a,b in substitutes:
        s = s.replace(a,b)
    return s

def format_src(s):
    """
    Format SAGE documentation for viewing with IPython.
    """
    if not isinstance(s, str):
        raise TypeError, "s must be a string"
    docs = set([])
    import sage.all
    while True:
        i = s.find("<<<")
        if i == -1: break
        j = s[i+3:].find('>>>')
        if j == -1: break
        obj = s[i+3:i+3+j]
        if obj in docs:
            t = ''
        else:
            x = eval('sage.all.%s'%obj, locals())
            t = my_getsource(x, False)
            docs.add(obj)
        if t is None:
            print x
            t = ''
        s = s[:i] + '\n' + t + s[i+6+j:]

    return s


###############################

def search_sage(s, extra=''):
    """
    Search sage source code for lines containing s.  The search is not
    case sensitive.

    For this to work the "sage" command must be in your PATH,
    and you must have grep installed.
    """
    from sage.misc.all import pager
    pager()(os.popen('sage -grep "%s" | grep "%s"'%(s,extra)).read())



###############################

#######################################
## Add detex'ing of documentation
#######################################
import sagedoc
import inspect
import sagex_inspect


def my_getdoc(obj):
    try:
        ds = obj._sage_doc_()
    except (AttributeError, TypeError):  # TypeError for interfaces
        try:
            ds = inspect.getdoc(obj)
        except:
            return None
    if ds is None:
        return None
    return sagedoc.format(ds)

def my_getsource(obj, is_binary):
    try:
        s = sagex_inspect.getsource(obj, is_binary)
        return sagedoc.format_src(s)
    except Exception, msg:
        print 'Error getting source', msg
        return None


