from sage.misc.latex import latex
from sage.misc.sage_eval import sage_eval

def math_parse(s):
    r"""
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
    """
    t = ''
    while True:
        i = s.find('$')
        if i == -1:
            return t + s
        elif i > 0 and s[i-1] == '\\':
            t += s[:i-1] + '$'
            s = s[i+1:]
        elif i-1 < len(s) and s[i+1] == '$':
            typ = 'div'
        else:
            typ = 'span'
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
    def __call__(self, *args, **kwds):
        return HTMLExpr(self.eval(*args, **kwds))

    def eval(self, s, globals={}, locals={}):
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
