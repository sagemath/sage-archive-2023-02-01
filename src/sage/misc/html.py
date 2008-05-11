"""
HTML typesetting for the notebok
"""

########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
########################################################################

from sage.misc.latex import latex
from sage.misc.sage_eval import sage_eval

def math_parse(s):
    r"""
    Turn the HTML-ish string s that can have $$ and $'s in it into
    pure HTML.  See below for a precise definition of what this means.

    INPUT:
        s -- a string
    OUTPUT:
        a string.

    Do the following:
    \begin{verbatim}
       * Replace all $ text $'s by
          <span class='math'> text </span>
       * Replace all $$ text $$'s by
          <div class='math'> text </div>
       * Replace all \$'s by $'.s  Note that in
         the above two cases nothing is done if the $
         is preceeded by a backslash.
    \end{verbatim}

    EXAMPLES:
        sage: sage.misc.html.math_parse('This is $2+2$.')
        'This is <span class="math">2+2</span>.'
        sage: sage.misc.html.math_parse('This is $$2+2$$.')
        'This is <div class="math">2+2</div>.'

    TESTS:
        sage: sage.misc.html.math_parse(r'This \$\$is $2+2$.')
        'This $$is <span class="math">2+2</span>.'
    """
    # Below t always has the "parsed so far" version of s, and s is
    # just the part of the original input s that hasn't been parsed.
    t = ''
    while True:
        i = s.find('$')
        if i == -1:
            # No dollar signs -- definitely done.
            return t + s
        elif i > 0 and s[i-1] == '\\':
            # A dollar sign with a backslash right before it, so
            # we ignore it by sticking it in the parsed string t
            # and skip to the next iteration.
            t += s[:i-1] + '$'
            s = s[i+1:]
            continue
        elif i+1 < len(s) and s[i+1] == '$':
            # Found a math environment. Double dollar sign so div mode.
            typ = 'div'
        else:
            # Found math environment. Single dollar sign so span mode.
            typ = 'span'

        # Now find the matching $ sign and form the span or div.
        j = s[i+2:].find('$')
        if j == -1:
            j = len(s)
            s += '$'
            if typ == 'div':
                s += '$$'
        else:
            j += i + 2
        if typ == 'div':
            txt = s[i+2:j]
        else:
            txt = s[i+1:j]
        t += s[:i] + '<%s class="math">%s</%s>'%(typ,
                      ' '.join(txt.splitlines()), typ)
        s = s[j+1:]
        if typ == 'div':
            s = s[1:]
    return t

class HTMLExpr(str):
    def __init__(self, x):
        str.__init__(self, x)

    def __repr__(self):
        return str(self)

class HTML:
    def __init__(self):
        """
        Display the given html expression in the notebook.

        INPUT:
            s -- a string

        OUTPUT:
            prints a code that embeds html in the output.

        By default in the notebook an output cell has two parts, first a plain text
        preformat part, then second a general html part (not pre).   If you call
        html(s) at any point then that adds something that will be displayed
        in the preformated part in html.

        EXAMPLES:
            sage: html('<a href="http://sagemath.org">sagemath</a>')
            <html><font color='black'><a href="http://sagemath.org">sagemath</a></font></html>
        """

    def __call__(self, s, globals=None, locals=None):
        """
        EXAMPLES:
            sage: html('<hr>')
            <html><font color='black'><hr></font></html>
        """
        return HTMLExpr(self.eval(s, globals, locals))

    def eval(self, s, globals=None, locals=None):
        """
        EXAMPLES:
            sage: html.eval('<hr>')
            <html><font color='black'><hr></font></html>
            ''
        """
        if globals is None:
            globals = {}
        if locals is None:
            locals = {}
        s = str(s)
        s = math_parse(s)
        t = ''
        while len(s) > 0:
            i = s.find('<sage>')
            if i == -1:
                 t += s
                 break
            j = s.find('</sage>')
            if j == -1:
                 t += s
                 break
            t += s[:i] + '<span class="math">%s</span>'%\
                     latex(sage_eval(s[6+i:j], locals=locals))
            s = s[j+7:]
        print "<html><font color='black'>%s</font></html>"%t
        return ''

html = HTML()
