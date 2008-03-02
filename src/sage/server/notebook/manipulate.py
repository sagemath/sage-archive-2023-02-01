r"""
Manipulate Sage functions in the notebook

This module implements a manipulate decorator for function in the Sage
notebook.

The controls are:
\begin{itemize}
    \item TextBox -- a text box
    \item Slider -- a slider
\end{itemize}

AUTHOR:
    -- William Stein (2008-03-02): initial version at Sage/Enthought
       Days 8 in Texas
    -- Jason Grout (2008-03): collaborators substantially on the
       design and prototypes.

TODO/PLAN:
   [ ] sliders
   [ ] default value
   [ ] more widgets
   [ ] better widget layout controls
   [ ] 100% doctest coverage
"""

import inspect

SAGE_CELL_ID = 0
vars = {}

def html(s):
    """
    Render the input string s in a form that tells the notebook
    to display it in the HTML portion of the output.

    INPUT:
        s -- a string

    OUTPUT:
        string -- html format
    """
    print "<html>%s</html>"%s

class ManipulateControl:
    """
    Base class for manipulate controls.
    """
    def __init__(self, f, var):
        """
        Create a new manipulate control.

        INPUT:
             f -- a Python function (that's being decorated)
             var -- name of variable that this control manipulates
             SAGE_CELL_ID -- uses this global variable
        """
        self.__var = var
        self.__cell_id = SAGE_CELL_ID
        self.__f = f

    def __repr__(self):
        return "A ManipulateControl (abstract base class)"

    def manipulate(self):
        """
        Return a string that when evaluated in Javascript calls the
        javascript manipulate function with appropriate inputs for
        this control.

        OUTPUT:
            string -- that is meant to be evaluated in Javascript
        """
        return 'manipulate(%s, "sage.server.notebook.manipulate.vars[%s][\\"%s\\"]=sage_eval(r\\"\\"\\""+%s+"\\"\\"\\", globals())\\n%s()");'%(
            self.cell_id(), self.cell_id(), self.var(), self.value(), self.function_name())

    def function_name(self):
        """
        Returns the name of the function that this control manipulates.

        OUTPUT:
            string -- name of a function as a string
        """
        return self.__f.__name__

    def var(self):
        """
        Return the name of the variable that this control manipulates.

        OUTPUT:
            string -- name of a variable as a string.
        """
        return self.__var

    def cell_id(self):
        """
        Return the id of the cell that contains this manipulate control.

        OUTPUT:
            integer -- id of cell that this control manipulates
        """
        return self.__cell_id

class TextBox(ManipulateControl):
    """
    A text box manipulate control.
    """
    def __repr__(self):
        return "A TextBox manipulate control"

    def value(self):
        """
        Return javascript string that will give the
        value of this control element.

        OUTPUT:
             string -- javascript
        """
        return "this.value"

    def render(self):
        """
        Render this control as a string.

        OUTPUT:
             string -- html format
        """
        return """
        %s: <input type='text' value='<?%s>' onchange='%s'></input>
        """%(self.var(), self.var(), self.manipulate())

class Slider(ManipulateControl):
    """
    A slider manipulate control.
    """
    def __init__(self, f, var, values):
        """
        Create a slider manipulate control that takes on the given
        list of values.
        """
        ManipulateControl.__init__(self, f, var)
        self.__values = values

    def __repr__(self):
        return "A Slider manipulate control"

    def value(self):
        """
        Return javascript string that will give the
        value of this control element.

        OUTPUT:
             string -- javascript
        """
        return "this.value"

    def render(self):
        """
        Render this control as a string.

        OUTPUT:
             string -- html format
        """
        return """
        SLIDER %s: <input type='text' value='<?%s>' onchange='%s'></input>
        """%(self.var(), self.var(), self.manipulate())


class ManipulateCanvas:
    """
    Base class for manipulate canvases. This is where all the controls
    along with the output of the manipulated function are layed out
    and rendered.
    """
    def __init__(self, controls):
        """
        Create a manipulate canvas.

        INPUT:
            controls -- a list of ManipulateControl instances.
        """
        self.__controls = controls

    def render_output(self):
        """
        Render in text (html) form the output portion of the manipulate canvas.

        The output contains two special tags, <?TEXT> and <?HTML>,
        which get replaced at runtime by the text and html parts
        of the output of running the function.

        OUTPUT:
            string -- html
        """

        return """
        <table bgcolor=black cellpadding=3><tr><td bgcolor=white>
        <?TEXT>
           <table border=0 width=800px>
           <tr><td align=center>  <?HTML>  </td></tr>
           </table>
        </td></tr></table>
        """

    def render_controls(self):
        """
        Render in text (html) form all the input controls.

        OUTPUT:
            string -- html
        """
        # This will need some sophisticated layout querying of the c's, maybe.
        return ''.join([c.render() for c in self.__controls])

    def render(self):
        """
        Render in text (html) the entire manipulate canvas.

        OUTPUT:
            string -- html
        """
        return "%s%s"%(self.render_controls(), self.render_output())


def manipulate(f):
    """
    Decorate a function f to make a manipulate version of f.
        @manipulate
        def foo(n,m):
            ...
    """

    (args, varargs, varkw, defaults) = inspect.getargspec(f)

    if defaults is None:
        defaults = []

    n = len(args) - len(defaults)
    controls = [TextBox(f, v) for v in args[:n]] + \
               [defaults[i].render(f, args[i+n]) for i in range(len(defaults))]

    C = ManipulateCanvas(controls)

    vars[SAGE_CELL_ID] = {}
    d = vars[SAGE_CELL_ID]
    for v in args:
         d[v] = ''

    html(C.render())

    def g():
        return f(*[d[args[i]] for i in range(len(args))])
    return g


######################################################
# Actual control objects that the user passes in
######################################################
class control:
    pass

class text_box(control):
    def __init__(self, default):
        """
        INPUT:
            default -- string (the default value)
        """
        self.__default = default

    def __repr__(self):
        return "A manipulate text box control with default value '%s'"%self.__default

    def render(self, f, var):
        """
        INPUT:
            f -- a Python function
            var -- a string (variable; one of the variable names input to f)
        """
        return TextBox(f, var)

class slider(control):
    def __init__(self, vmin, vmax=None, steps=30):
        if isinstance(vmin, list):
            self.__values = vmin
        else:
            if vmax is None:
                vmax = vmin
                vmin = 0
            steps = int(steps)
            if steps <= 0:
                self.__values = [vmin, vmax]
            else:
                step = (vmax-vmin)/steps  # hard coded
                self.__values = [vmin + i*step for i in range(steps)] + [vmax]

    def __repr__(self):
        return "A manipulate slider control [%s - %s]."%(min(self.__values), max(self.__values))

    def render(self, f, var):
        return Slider(f, var, self.__values)

