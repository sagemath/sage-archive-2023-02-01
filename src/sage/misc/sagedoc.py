r"""
Format Sage documentation for viewing with IPython and the notebook

AUTHORS:

- William Stein (2005): initial version.
- Nick Alexander (2007): nodetex functions
- Nick Alexander (2008): search_src, search_def improvements
- Martin Albrecht (2008-03-21): parse LaTeX description environments in sagedoc
- J. H. Palmieri (2009-04-11): fix for #5754 plus doctests
- Dan Drake (2009-05-21): refactor search_* functions, use system `find' instead of sage -grep
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
from subprocess import Popen, PIPE

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
    """
    Remove the LaTeX command ``cmd`` from the string ``s``.  This
    function is used by ``detex``.

    INPUT:

    - ``s`` - (string) string from which to remove the command

    - ``cmd`` - (string) command to be removed.  This should be a
      command which takes a single argument, like 'emph' or 'url'; the
      command is removed, but its argument is not.

    - ``left``, ``right`` - (string, optional, default '') add these
      strings at the left and right ends of the command. See the
      examples.

    EXAMPLES::

        sage: from sage.misc.sagedoc import _rmcmd
        sage: _rmcmd('Check out \\url{http://www.sagemath.org}.', 'url')
        'Check out http://www.sagemath.org.'
        sage: _rmcmd('Text in \\emph{italics} looks like this.', 'emph', '*', '*')
        'Text in *italics* looks like this.'
        sage: _rmcmd('This is a \\very{silly} example.', 'very', right='!?')
        'This is a silly!? example.'
    """
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
    """nodetex
    This strips LaTeX commands from a string; it is used by the
    ``format`` function to process docstrings for display from the
    command line interface.

    INPUT: ``s``, a string.

    OUTPUT: string

    EXAMPLES::

        sage: from sage.misc.sagedoc import detex
        sage: detex(r'Some math: `n \geq k`.  A website: \url{sagemath.org}.')
        'Some math: `n >= k`.  A website: sagemath.org.'
        sage: detex(r'More math: `x \mapsto y`.  {\bf Bold face}.')
        'More math: `x  |-->  y`.  { Bold face}.'
        sage: detex(r'$a, b, c, \ldots, z$')
        'a, b, c, ..., z'
    """
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
    Format Sage documentation ``s`` for viewing with IPython.

    This calls ``detex`` on ``s`` to convert LaTeX commands to plain
    text, and if ``s`` contains a string of the form "<<<obj>>>",
    then it replaces it with the docstring for "obj".

    INPUT: ``s`` - string

    OUTPUT: string

    EXAMPLES::

        sage: from sage.misc.sagedoc import format
        sage: identity_matrix(2).rook_vector.__doc__[115:184]
        'Let `A` be a general `m` by `n`\n        (0,1)-matrix with `m \\le n`. '
        sage: format(identity_matrix(2).rook_vector.__doc__[115:184])
        'Let `A` be a general `m` by `n`\n        (0,1)-matrix with `m <= n`. '

    If the first line of the string is 'nodetex', remove 'nodetex' but
    don't modify any TeX commands::

        sage: format("nodetex\n`x \\geq y`")
        '\n`x \\geq y`'

    Testing a string enclosed in triple angle brackets::

        sage: format('<<<identity_matrix')
        '<<<identity_matrix'
        sage: format('identity_matrix>>>')
        'identity_matrix>>>'
        sage: format('<<<identity_matrix>>>')[:28]
        '\nDefinition: identity_matrix'
    """
    if not isinstance(s, str):
        raise TypeError, "s must be a string"

    # parse directives at beginning of docstring
    # currently, only 'nodetex' is supported.
    # 'no' + 'doctest' may be supported eventually (don't type that as
    # one word, or the whole file will not be doctested).
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
            t0 = sage.misc.sageinspect.sage_getdef(x, obj)
            t1 = my_getdoc(x)
            t = 'Definition: ' + t0 + '\n\n' + t1
            docs.add(obj)
        s = s[:i] + '\n' + t + s[i+6+j:]

    if 'nodetex' not in directives:
        s = detex(s)
    else:
        # strip the 'nodetex' directive from s
        s = s.replace('nodetex', '', 1)
    return s

def format_src(s):
    """
    Format Sage source code ``s`` for viewing with IPython.

    If ``s`` contains a string of the form "<<<obj>>>", then it
    replaces it with the source code for "obj".

    INPUT: ``s`` - string

    OUTPUT: string

    EXAMPLES::

        sage: from sage.misc.sagedoc import format_src
        sage: format_src('unladen swallow')
        'unladen swallow'
        sage: format_src('<<<Sq>>>')[5:15]
        'Sq(*nums):'
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

def _findcmd(path, exts, what):
    """
    Return a list of strings that represent a command suitable for input
    to Popen. The command uses the system ``find`` and ``grep`` commands
    to search the Sage library.

    INPUT:

    - ``path`` - path to search, relative to ``$SAGE_ROOT``. You don't
      need to start with a forward slash.
    - ``exts`` - list of file extensions to search.
    - ``what`` - search term passed to ``grep``.

    OUTPUT:

    A list of strings.

    EXAMPLES::

        sage: from sage.misc.sagedoc import _findcmd
        sage: _findcmd('devel/sage/sage', ['py', 'pyx'], 'foo bar')
        ['find',
         .../devel/sage/sage/',
         '(',
         '-name',
         '*.py',
         '-o',
         '-name',
         '*.pyx',
         ')',
         '-exec',
         'grep',
         '-i',
         '-H',
         'foo bar',
         '{}',
         '+']
       sage: _findcmd('devel/sage/doc/output', ['html'], 'xyz')
       ['find',
        .../devel/sage/doc/output/',
        '(',
        '-name',
        '*.html',
        ')',
        '-exec',
        'grep',
        '-i',
        '-H',
        'xyz',
        '{}',
        '+']
    """
    extstr = '( ' + ' -o '.join(['-name *.%s' % ext for ext in exts]) + ' )'
    # we add on "what" separately, since we need to respect spaces there.
    return (r"""find %s/ %s -exec grep -i -H""" %
         (os.environ['SAGE_ROOT'] + '/' + path,
         extstr)).split() + [what, '{}', '+']


def _search_src_or_doc(what, string, extra1='', extra2='', extra3='', extra4='', extra5='', interact=True):
    r"""
    Search the Sage library or documentation for lines containing
    ``string`` and possibly some other terms. This function is used by
    :func:`search_src`, :func:`search_doc`, and :func:`search_def`.

    INPUT:

    - ``what``: either ``'src'`` or ``'doc'``, according to whether you
      are searching the documentation or source code.
    - the rest of the input is the same as :func:`search_src`,
      :func:`search_doc`, and :func:`search_def`.

    OUTPUT:

    If ``interact`` is ``False``, a string containing the results;
    otherwise, the results are presented according to whether you are
    using the notebook or command-line interface.

    EXAMPLES::

        sage: from sage.misc.sagedoc import _search_src_or_doc
        sage: print _search_src_or_doc('src', 'matrix(', 'incidence_structures', 'self', '^combinat', interact=False) # random # long time
        misc/sagedoc.py:        sage: _search_src_or_doc('src', 'matrix(', 'incidence_structures', 'self', '^combinat', interact=False)
        combinat/designs/incidence_structures.py:        M1 = self.incidence_matrix()
        combinat/designs/incidence_structures.py:        A = self.incidence_matrix()
        combinat/designs/incidence_structures.py:        M = transpose(self.incidence_matrix())
        combinat/designs/incidence_structures.py:    def incidence_matrix(self):
        combinat/designs/incidence_structures.py:        A = self.incidence_matrix()
        combinat/designs/incidence_structures.py:        A = self.incidence_matrix()
        combinat/designs/incidence_structures.py:        #A = self.incidence_matrix()

    TESTS:

    The examples are nice, but marking them "random" means we're not
    really testing if the function works, just that it completes. These
    tests aren't perfect, but are reasonable.

    ::

        sage: len(_search_src_or_doc('src', 'matrix(', 'incidence_structures', 'self', 'combinat', interact=False).splitlines()) > 1
        True

        sage: 'abvar/homology' in _search_src_or_doc('doc', 'homology', 'variety', interact=False)
        True

        sage: 'divisors' in _search_src_or_doc('src', '^ *def prime', interact=False)
        True
    """
    if what == 'src':
        path = 'devel/sage/sage/'
        exts = ['py', 'pyx', 'pxd']
        title = 'Source Code'
    else:
        path = 'devel/sage/doc/output/'
        exts = ['html']
        title = 'Documentation'

    cmds = [Popen(_findcmd(path, exts, string), stdout=PIPE)]
    for extra in [extra1, extra2, extra3, extra4, extra5]:
        if extra:
            cmds.append(Popen(['grep', '-i', extra], stdin=cmds[-1].stdout,
                              stdout=PIPE))

    strip = len(os.environ['SAGE_ROOT'] + path) + 2
    r = '\n'.join([l[strip:].lstrip("/") for l in cmds[-1].communicate()[0].splitlines()])

    if not interact:
        return r
    from sage.server.support import EMBEDDED_MODE
    if EMBEDDED_MODE:   # I.e., running from the notebook
        # format the search terms nicely
        terms = ('"%s"' % string) + ' '.join(['"%s"' % ex for ex in
          [extra1, extra2, extra3, extra4, extra5] if ex])
        print format_search_as_html(title, r, terms)
    else:
        # hard-code a 25-line screen into the pager; this works around a
        # problem with doctests: see
        # http://trac.sagemath.org/sage_trac/ticket/5806#comment:11
        from IPython.genutils import page
        page(r, screen_lines = 25)


def search_src(string, extra1='', extra2='', extra3='', extra4='', extra5='', interact=True):
    r"""
    Search Sage library source code for lines containing ``string``.
    The search is not case sensitive.

    The file paths in the output are relative to
    ``$SAGE_ROOT/devel/sage/sage``.

    INPUT:

    - ``string`` - a string to find in the Sage source code.

    - ``extra1``, ..., ``extra5`` - additional strings to require,
      passed as arguments to grep

    - ``interact`` - if ``False``, will simply return a string with all
      the matches. Otherwise, results will be displayed appropriately,
      according to whether you are using the notebook or the
      command-line interface. You should not ordinarily need to use
      this.

    The ``string`` and ``extraN`` arguments are passed on to the system
    ``grep``, so you can use whatever regular expressions ``grep``
    supports. The matches will always be case-insensitive.

    .. note::

        The ``extraN`` parameters are present only because
        ``search_src(string, *extras, interact=None)``
        is not parsed correctly by Python 2.6; see http://bugs.python.org/issue1909.

    EXAMPLES::

        sage: print search_src(" fetch(", "def", interact=False) # random # long time
        matrix/matrix0.pyx:    cdef fetch(self, key):
        matrix/matrix0.pxd:    cdef fetch(self, key)

        sage: print search_src(" fetch(", "def", "pyx", interact=False) # random # long time
        matrix/matrix0.pyx:    cdef fetch(self, key):

    A little recursive narcissism: let's do a doctest that searches for
    this function's doctests. Note that you can't put "sage:" in the
    doctest string because it will get replaced by the Python ">>>"
    prompt.

    ::

        sage: print search_src('^ *sage[:] .*search_src(', interact=False) # long time
        misc/sagedoc.py:        ... print search_src(" fetch(", "def", interact=False) # random # long time
        misc/sagedoc.py:        ... print search_src(" fetch(", "def", "pyx", interact=False) # random # long time
        misc/sagedoc.py:        ... print search_src('^ *sage[:] .*search_src(', interact=False) # long time
        misc/sagedoc.py:        ... len(search_src("matrix", interact=False).splitlines()) > 10000 # long time
        misc/sagedoc.py:        ... print search_src('matrix', 'column', 'row', '0', 'sub', 'start', interact=False) # random # long time


    TESTS:

    As of this writing, there's about 10900 lines in the Sage library that
    contain "matrix"; it seems safe to assume we'll continue to have
    over 10000 such lines::

        sage: len(search_src("matrix", interact=False).splitlines()) > 10000 # long time
        True

    Check that you can pass 5 parameters::

        sage: print search_src('matrix', 'column', 'row', '0', 'sub', 'start', interact=False) # random # long time
        matrix/matrix0.pyx:598:        Get The 2 x 2 submatrix of M, starting at row index and column
        matrix/matrix0.pyx:607:        Get the 2 x 3 submatrix of M starting at row index and column index
        matrix/matrix0.pyx:924:        Set the 2 x 2 submatrix of M, starting at row index and column
        matrix/matrix0.pyx:933:        Set the 2 x 3 submatrix of M starting at row index and column

    """
    return _search_src_or_doc('src', string, extra1=extra1, extra2=extra2, extra3=extra3, extra4=extra4, extra5=extra5, interact=interact)


def search_doc(string, extra1='', extra2='', extra3='', extra4='', extra5='', interact=True):
    """
    Search Sage HTML documentation for lines containing ``string``. The
    search is not case sensitive.

    The file paths in the output are relative to
    ``$SAGE_ROOT/devel/sage/doc/output``.

    INPUT:

    - ``string`` - a string to find in the Sage source code.

    - ``extra1``, ..., ``extra5`` - additional strings to require,
      passed as arguments to grep

    - ``interact`` - if ``False``, will simply return a string with all
      the matches. Otherwise, results will be displayed appropriately,
      according to whether you are using the notebook or the
      command-line interface. You should not ordinarily need to use
      this.

    The ``string`` and ``extraN`` arguments are passed on to the system
    ``grep``, so you can use whatever regular expressions ``grep``
    supports. The matches will always be case-insensitive.

    .. note::

        The ``extraN`` parameters are present only because
        ``search_src(string, *extras, interact=None)``
        is not parsed correctly by Python 2.6; see http://bugs.python.org/issue1909.

    EXAMPLES::

        sage: search_doc('this creates a polynomial ring') # random # this function has no output: it just prints a string
        html/en/tutorial/tour_polynomial.html:<p>This creates a polynomial ring and tells Sage to use (the string)
    """
    return _search_src_or_doc('doc', string, extra1=extra1, extra2=extra2, extra3=extra3, extra4=extra4, extra5=extra5, interact=interact)


def search_def(name, extra1='', extra2='', extra3='', extra4='', extra5='', interact=True):
    r"""
    Search Sage library source code for function definitions containing
    ``name``. The search is not case sensitive.

    INPUT:

    - ``name`` - a string to find in the names of functions in the Sage
      source code.
    - ``extra1``, ..., ``extra4`` - additional strings to require, as
      in :func:`search_src`.
    - ``interact`` - if ``False``, will simply return a string with all
      the matches. Otherwise, results will be displayed appropriately,
      according to whether you are using the notebook or the
      command-line interface. You should not ordinarily need to use
      this.

    The ``string`` and ``extraN`` arguments are passed on to the system
    ``grep``, so you can use whatever regular expressions ``grep``
    supports. The matches will always be case-insensitive.

    .. note::

        The regular expression used by this function only finds function
        definitions that are preceded by spaces, so if you use tabs on a
        "def" line, this function will not find it. As tabs are not
        allowed in Sage library code, this should not be a problem.

    .. note::

        The ``extraN`` parameters are present only because
        ``search_src(string, *extras, interact=None)``
        is not parsed correctly by Python 2.6; see http://bugs.python.org/issue1909.

    EXAMPLES::

        sage: print search_def("fetch", interact=False) # random # long time
        matrix/matrix0.pyx:    cdef fetch(self, key):
        matrix/matrix0.pxd:    cdef fetch(self, key)

        sage: print search_def("fetch", "pyx", interact=False) # random # long time
        matrix/matrix0.pyx:    cdef fetch(self, key):
    """
    return _search_src_or_doc('src', '^ *def.*%s' % name, extra1=extra1,
                      extra2=extra2, extra3=extra3, extra4=extra4,
                      extra5=extra5, interact=interact)


def format_search_as_html(what, r, search):
    """
    Format the output from ``search_src``, ``search_def``, or
    ``search_doc`` as html, for use in the notebook.

    INPUT:

    - ``what`` - (string) what is being search
    - ``r`` - (string) the results of the search
    - ``search`` - (string) what was being searched for

    This function parses ``r``: it should have the form 'FILENAME:
    string' where FILENAME is the file in which the string, which
    matched the search was found.  If FILENAME ends in '.html', then
    this is part of the documentation; otherwise, it is in the source
    code.  In either case, an appropriate link is created.

    EXAMPLES::

        sage: from sage.misc.sagedoc import format_search_as_html
        sage: format_search_as_html('Source', 'algebras/steenrod_algebra_element.py:        an antihomomorphism: if we call the antipode `c`, then', 'antipode antihomomorphism')
        '<html><font color="black"><h2>Search Source: antipode antihomomorphism</h2></font><font color="darkpurple"><ol><li><a href="/src/algebras/steenrod_algebra_element.py"><tt>algebras/steenrod_algebra_element.py</tt></a>\n</ol></font></html>'
        sage: format_search_as_html('Other', 'html/en/reference/sage/algebras/steenrod_algebra_element.html:an antihomomorphism: if we call the antipode <span class="math">c</span>, then', 'antipode antihomomorphism')
        '<html><font color="black"><h2>Search Other: antipode antihomomorphism</h2></font><font color="darkpurple"><ol><li><a href="/doc/live/reference/sage/algebras/steenrod_algebra_element.html"><tt>reference/sage/algebras/steenrod_algebra_element.html</tt></a>\n</ol></font></html>'
    """
    s = '<html>'
    s += '<font color="black">'
    s += '<h2>Search %s: %s</h2>'%(what, search)
    s += '</font>'
    s += '<font color="darkpurple">'
    s += '<ol>'

    files = set([])
    for L in r.splitlines():
        i = L.find(':')
        if i != -1:
            files.add(L[:i])
    files = list(files)
    files.sort()
    for F in files:
        if F.endswith('.html'):
            F = F.split('/', 2)[2]
            url = '/doc/live/' + F
        else:
            # source code
            url = '/src/' + F
        s += '<li><a href="%s"><tt>%s</tt></a>\n'%(url, F)
    s += '</ol>'
    s += '</font>'
    s += '</html>'
    return s



#######################################
## Add detex'ing of documentation
#######################################
import inspect
import sageinspect

def my_getdoc(obj):
    """
    Retrieve the documentation for ``obj``.

    INPUT: ``obj`` - a Sage object, function, etc.

    OUTPUT: its documentation (string)

    EXAMPLES::

        sage: from sage.misc.sagedoc import my_getdoc
        sage: s = my_getdoc(identity_matrix)
        sage: type(s)
        <type 'str'>
    """
    try:
        ds = obj._sage_doc_()
    except (AttributeError, TypeError):  # TypeError for interfaces
        try:
            ds = sageinspect.sage_getdoc(obj)
        except:
            return None
    if ds is None:
        return None
    return format(ds)

def my_getsource(obj, is_binary):
    """
    Retrieve the source code for ``obj``.

    INPUT:

    - ``obj`` - a Sage object, function, etc.
    - ``is_binary`` - (boolean) ignored argument.

    OUTPUT: its documentation (string)

    EXAMPLES::

        sage: from sage.misc.sagedoc import my_getsource
        sage: s = my_getsource(identity_matrix, True)
        sage: s[:19]
        'def identity_matrix'
    """
    try:
        s = sageinspect.sage_getsource(obj, is_binary)
        return format_src(s)
    except Exception, msg:
        print 'Error getting source:', msg
        return None
