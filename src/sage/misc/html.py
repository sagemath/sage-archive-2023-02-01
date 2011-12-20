"""
HTML typesetting for the notebook
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

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
       * Replace all \$'s by $'s.  Note that in
         the above two cases nothing is done if the $
         is preceeded by a backslash.
       * Replace all \[ text \]'s by
         <div class='math'> text </div>
    \end{verbatim}

    EXAMPLES:
        sage: sage.misc.html.math_parse('This is $2+2$.')
        'This is <span class="math">2+2</span>.'
        sage: sage.misc.html.math_parse('This is $$2+2$$.')
        'This is <div class="math">2+2</div>.'
        sage: sage.misc.html.math_parse('This is \\[2+2\\].')
        'This is <div class="math">2+2</div>.'
        sage: sage.misc.html.math_parse(r'This is \[2+2\].')
        'This is <div class="math">2+2</div>.'

    TESTS:
        sage: sage.misc.html.math_parse(r'This \$\$is $2+2$.')
        'This $$is <span class="math">2+2</span>.'
    """
    # first replace \\[ and \\] by <div class="math"> and </div>, respectively.
    while True:
        i = s.find('\\[')
        if i == -1:
            break
        else:
            s = s[:i] + '<div class="math">' + s[i+2:]
            j = s.find('\\]')
            if j == -1:  # missing right-hand delimiter, so add one
                s = s + '</div>'
            else:
                s = s[:j] + '</div>' + s[j+2:]

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
    def __repr__(self):
        return str(self)

class HTML:
    def __init__(self):
        """
        Display the given HTML expression in the notebook.

        INPUT:
            s -- a string

        OUTPUT:
            prints a code that embeds HTML in the output.

        By default in the notebook an output cell has two parts, first a plain text
        preformat part, then second a general HTML part (not pre).   If you call
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
        r"""
        EXAMPLES::

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

    def table(self, x, header = False):
        r"""
        Print a nested list as a HTML table.  Strings of html
        will be parsed for math inside dollar and double-dollar signs.
        2D graphics will be displayed in the cells.  Expressions will
        be latexed.


        INPUT:

          - ``x`` -- a list of lists (i.e., a list of table rows)
          - ``header`` -- a row of headers.  If True, then the first
             row of the table is taken to be the header.

        EXAMPLES::

            sage: html.table([(i, j, i == j) for i in [0..1] for j in [0..1]])
            <html>
            <div class="notruncate">
            <table class="table_form">
            <tbody>
            <tr class ="row-a">
            <td><span class="math">0</span></td>
            <td><span class="math">0</span></td>
            <td><span class="math">\mathrm{True}</span></td>
            </tr>
            <tr class ="row-b">
            <td><span class="math">0</span></td>
            <td><span class="math">1</span></td>
            <td><span class="math">\mathrm{False}</span></td>
            </tr>
            <tr class ="row-a">
            <td><span class="math">1</span></td>
            <td><span class="math">0</span></td>
            <td><span class="math">\mathrm{False}</span></td>
            </tr>
            <tr class ="row-b">
            <td><span class="math">1</span></td>
            <td><span class="math">1</span></td>
            <td><span class="math">\mathrm{True}</span></td>
            </tr>
            </tbody>
            </table>
            </div>
            </html>

            sage: html.table(["Functions $f(x)$", sin(x), cos(x)], header = True)
            <html>
            <div class="notruncate">
            <table class="table_form">
            <tbody>
            <tr>
            <th>Functions <span class="math">f(x)</span></th>
            </tr>
            <tr class ="row-a">
            <td><span class="math">\sin\left(x\right)</span></td>
            </tr>
            <tr class ="row-b">
            <td><span class="math">\cos\left(x\right)</span></td>
            </tr>
            </tbody>
            </table>
            </div>
            </html>

            sage: html.table([(x,n(sin(x), digits=2)) for x in [0..3]], header = ["$x$", "$\sin(x)$"])
            <html>
            <div class="notruncate">
            <table class="table_form">
            <tbody>
            <tr>
            <th><span class="math">x</span></th>
            <th><span class="math">\sin(x)</span></th>
            </tr>
            <tr class ="row-a">
            <td><span class="math">0</span></td>
            <td><span class="math">0.00</span></td>
            </tr>
            <tr class ="row-b">
            <td><span class="math">1</span></td>
            <td><span class="math">0.84</span></td>
            </tr>
            <tr class ="row-a">
            <td><span class="math">2</span></td>
            <td><span class="math">0.91</span></td>
            </tr>
            <tr class ="row-b">
            <td><span class="math">3</span></td>
            <td><span class="math">0.14</span></td>
            </tr>
            </tbody>
            </table>
            </div>
            </html>

        """
        import types
        from sage.misc.all import latex
        from itertools import cycle
        if isinstance(x, types.GeneratorType):
            x = list(x)
        if isinstance(x, (list, tuple)):
            rows = len(x)
            if rows > 0:
                # if the table has less then 100 rows, don't truncate the output in the notebook
                if rows <= 100:
                    print "<html>\n<div class=\"notruncate\">\n<table class=\"table_form\">\n<tbody>"
                else:
                    print "<html>\n<div class=\"truncate\">\n<table class=\"table_form\">\n<tbody>"

                if header is True:
                    header=x[0]
                    x = list(x[1:])

                if header is not False:
                    print "<tr>"
                    self._table_columns(header, True)
                    print "</tr>"

                for row_class, row in zip(cycle(["row-a", "row-b"]), x):
                    print "<tr class =\"%s\">" % row_class
                    self._table_columns(row, False)
                    print "</tr>"
                print "</tbody>\n</table>\n</div>\n</html>"

    def _table_columns(self, row, header=False):
        r"""
        Print the items of a list as the columns of a HTML table.

        TESTS::

            sage: html._table_columns(["a $x^2$",1, sin(x)])
            <td>a <span class="math">x^2</span></td>
            <td><span class="math">1</span></td>
            <td><span class="math">\sin\left(x\right)</span></td>
            sage: html._table_columns("a", header=True)
            <th>a</th>
        """
        column_tag = "<th>%s</th>" if header else "<td>%s</td>"
        from sage.plot.plot import Graphics
        import types
        if isinstance(row, types.GeneratorType):
            row = list(row)
        elif not isinstance(row, (list, tuple)):
            row = [row]

        for column in xrange(len(row)):
            if isinstance(row[column], Graphics):
                print column_tag % row[column].show(linkmode = True)
            elif isinstance(row[column], str):
                print column_tag % math_parse(row[column])
            else:
                print column_tag % ('<span class="math">%s</span>' % latex(row[column]))

    def iframe(self, url, height=400, width=800):
        r"""
        Put an existing web page into a worksheet.

        INPUT:

            - ``url`` -- a url string, either with or without URI scheme
              (defaults to "http").
            - ``height`` -- the number of pixels for the page height.
              Defaults to 400.
            - ``width`` -- the number of pixels for the page width.
              Defaults to 800.

        OUTPUT:

            Opens the url in a worksheet. If the url is a regular web page it
            will appear in the worksheet. This was originally intended to bring
            GeoGebra worksheets into Sage, but it can be used for many other
            purposes.

        EXAMPLES::

            sage: html.iframe("sagemath.org")
            <html><font color='black'><iframe height="400" width="800"
            src="http://sagemath.org"></iframe></font></html>
            sage: html.iframe("http://sagemath.org",30,40)
            <html><font color='black'><iframe height="30" width="40"
            src="http://sagemath.org"></iframe></font></html>
            sage: html.iframe("https://sagemath.org",30)
            <html><font color='black'><iframe height="30" width="800"
            src="https://sagemath.org"></iframe></font></html>
            sage: html.iframe("/home/admin/0/data/filename")
            <html><font color='black'><iframe height="400" width="800"
            src="/home/admin/0/data/filename"></iframe></font></html>
            sage: html.iframe('data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA'
            ... 'AUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBA'
            ... 'AO9TXL0Y4OHwAAAABJRU5ErkJggg=="')
            <html><font color='black'><iframe height="400" width="800"
            src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBAAO9TXL0Y4OHwAAAABJRU5ErkJggg==""></iframe></font></html>

        AUTHOR:

        - Bruce Cohen (2011-06-14)
        """
        if ":" not in url and not url.startswith('/'):
            url = "http://" + url
        string = ( '<iframe height="%d" width="%d" src="%s"></iframe>' %
                    (height, width, url) )
        return html(string)

html = HTML()
