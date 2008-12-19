r"""nodoctest
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
#('\\item', '*'), \
substitutes = [('\\_','_'),\
               ('\\item', '* '), \
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
               ('$',''), ('\\',''), ('backslash','\\'), \
               ('begin{enumerate}',''), ('end{enumerate}',''), \
               ('begin{description}',''), ('end{description}',''), \
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

import re
itempattern = re.compile(r"\\item\[?([^]]*)\]? *(.*)")
itemreplace = r"* \1 \2"

def detex(s):
    s = _rmcmd(s, 'url')
    s = _rmcmd(s, 'code')
    s = _rmcmd(s, 'class')
    s = _rmcmd(s, 'mbox')
    s = _rmcmd(s, 'text')
    s = _rmcmd(s, 'section')
    s = _rmcmd(s, 'subsection')
    s = _rmcmd(s, 'subsubsection')
    s = _rmcmd(s, 'note', 'NOTE: ', '')
    s = _rmcmd(s, 'emph', '*', '*')

    s = re.sub(itempattern, itemreplace, s)

    for a,b in substitutes:
        s = s.replace(a,b)
    return s

def format(s):
    """
    Format SAGE documentation for viewing with IPython.
    """
    if not isinstance(s, str):
        raise TypeError, "s must be a string"

    # parse directives at beginning of docstring
    # currently, only 'nodetex' is supported
    # eventually, 'nodoctest' might be supported
    first_newline = s.find('\n')
    if first_newline > -1:
        first_line = s[:first_newline]
    else:
        first_line = s
    directives = [ d.lower() for d in first_line.split(',') ]

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

    if 'nodetex' not in directives:
        s = detex(s)
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

def search_src(string, extra1='', extra2='', extra3='', extra4='', extra5='', interact=True):
    r"""
    Search \Sage library source code for lines containing \code{string}.
    The search is not case sensitive.

    INPUT:
        -- string: a string to find in the \Sage source code.
        -- extra1, ..., extra5: additional strings to require.

    NOTE:
        The extraN parameters are present only because
        search_src(string, *extras, interact=None)
        is not parsed correctly by Python 2.6; see http://bugs.python.org/issue1909.

    EXAMPLES:
        sage: print search_src(" fetch(", "def", interact=False) # random # long
        matrix/matrix0.pyx:    cdef fetch(self, key):
        matrix/matrix0.pxd:    cdef fetch(self, key)

        sage: print search_src(" fetch(", "def", interact=False) # random # long
        matrix/matrix0.pyx:    cdef fetch(self, key):
        matrix/matrix0.pxd:    cdef fetch(self, key)

        sage: print search_src(" fetch(", "def", "pyx", interact=False) # random # long
        matrix/matrix0.pyx:    cdef fetch(self, key):
    """
    cmd = 'sage -grep "%s"' % string
    for extra in [ extra1, extra2, extra3, extra4, extra5 ]:
        if not extra:
            continue
        cmd += '| grep "%s"' % extra

    old_banner = os.environ.has_key('SAGE_BANNER')
    if old_banner:
        old_banner = os.environ['SAGE_BANNER']
    os.environ['SAGE_BANNER'] = "no"
    r = os.popen(cmd).read()
    if old_banner == False:
        del os.environ['SAGE_BANNER']
    else:
        os.environ['SAGE_BANNER'] = old_banner

    if not interact:
        return r
    from sage.server.support import EMBEDDED_MODE
    if EMBEDDED_MODE:   # I.e., running from the notebook
        print format_search_as_html('Source Code', r, string + extra)
    else:
        from sage.misc.all import pager
        pager()(r)

def search_def(name, extra1='', extra2='', extra3='', extra4='', interact=True):
    r"""
    Search \Sage library source code for function names containing \code{name}.
    The search is not case sensitive.

    INPUT:
        -- name: a string to find in the names of functions in the \Sage source code.
        -- extra1, ..., extra4: see search_src.

    EXAMPLES:
        sage: print search_def("fetch", interact=False) # random # long
        matrix/matrix0.pyx:    cdef fetch(self, key):
        matrix/matrix0.pxd:    cdef fetch(self, key)

        sage: print search_def("fetch", "pyx", interact=False) # random # long
        matrix/matrix0.pyx:    cdef fetch(self, key):
    """
    return search_src(name, extra1='def ',
                      extra2=extra1, extra3=extra2,
                      extra4=extra3, extra5=extra4, interact=interact)

def format_search_as_html(what, r, search):
    s = '<html>'
    s += '<font color="black">'
    s += '<h2>Search %s: %s</h2>'%(what, search)
    s += '</font>'
    s += '<font color="darkpurple">'
    s += '<ol>'

    files = set([])
    for L in r.splitlines()[4:]:
        i = L.find(':')
        if i != -1:
            files.add(L[:i])
    files = list(files)
    files.sort()
    for F in files:
        if F.endswith('.html'):
            url = '/doc/live/' + F
        else:
            # source code
            url = '/src/' + F
        s += '<li><a href="%s"><tt>%s</tt></a>\n'%(url, F)
    s += '</ol>'
    s += '</font>'
    s += '</html>'
    return s

def search_doc(s, extra=''):
    """
    Full text search of the SAGE HTML documentation for lines
    containing s.  The search is not case sensitive.
    """
    cmd = 'sage -grepdoc "%s" | grep "%s"'%(s,extra)

    old_banner = os.environ.has_key('SAGE_BANNER')
    if old_banner:
        old_banner = os.environ['SAGE_BANNER']
    os.environ['SAGE_BANNER'] = "no"
    r = os.popen(cmd).read()
    if old_banner == False:
        del os.environ['SAGE_BANNER']
    else:
        os.environ['SAGE_BANNER'] = old_banner

    from sage.server.support import EMBEDDED_MODE
    if EMBEDDED_MODE:   # I.e., running from the notebook
        print format_search_as_html('Documentation', r, s + extra)
    else:
        from sage.misc.all import pager
        pager()(r)



###############################

#######################################
## Add detex'ing of documentation
#######################################
import sagedoc
import inspect
import sageinspect


def my_getdoc(obj):
    try:
        ds = obj._sage_doc_()
    except (AttributeError, TypeError):  # TypeError for interfaces
        try:
            ds = sageinspect.sage_getdoc(obj)
        except:
            return None
    if ds is None:
        return None
    return sagedoc.format(ds)

def my_getsource(obj, is_binary):
    try:
        s = sageinspect.sage_getsource(obj, is_binary)
        return sagedoc.format_src(s)
    except Exception, msg:
        print 'Error getting source:', msg
        return None


