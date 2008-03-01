SAGE_CELL_ID=0

def _text_box(f, var):
    print """<html>
    <input type='text' value='<?%s>' onchange='dynamic(%s, "sage.server.notebook.manipulate.%s="+this.value+"\\n%s()");'></input>
    <table bgcolor=black cellpadding=3><tr><td bgcolor=white>
<?TEXT>
    <table border=0 width=800px>
    <tr><td align=center>  <?HTML>  </td></tr></table>
    </td></tr></table>
    </html>
    """%(var, SAGE_CELL_ID, var, f.__name__)

def manipulate(f):
    """
    Decorate a function f to make a manipulatable version of f.
    """
    import inspect
    vars = inspect.getargspec(f)[0]
    print _text_box(f, vars[0])
    def g():
        return f(eval(vars[0]))
    return g
