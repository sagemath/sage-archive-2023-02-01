from sage.misc.latex import latex
from sage.misc.sage_eval import sage_eval


class HTML:
    def eval(self, s, globals={}, locals={}):
        s = str(s)
        s = s.replace('<$>','<span class="math">')
        s = s.replace('</$>','</span>')
        s = s.replace('<$$>','<div class="math">')
        s = s.replace('</$$>','</div>')
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
        print "<html>%s</html>"%t
        return ''

html = HTML()
