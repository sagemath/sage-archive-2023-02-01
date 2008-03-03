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

TODO:
   [ ] get sliders to work
   [ ] default values; values after move slider
   [ ] implement a color object

PLANS and IDEAS:
   [ ] automagically determine the type of control from the default
       value of the variable.  Here is how this will work:
        * u                      blank input field
        * u = (umin,umax)        slider; umin must not be a string
        * u = (umin,umax,du)     discrete slider
        * u = [1,2,3,4]          setter bar: automatically when there are are most 5
                                 elements; otherwise a drop down menu
        * u = ((xmin,xmax),(ymin,ymax))  2d slider
        * u = Graphic            a locator in a 2d plot  (Graphic is a 2d graphics objects)
        * u = True or u = False  a checkbox
        * u = color(...)         a color slider
        * u = "string"           text input field
        * u = ('label', obj)     obj can be any of the above; control is labeled with
                                 the given label
        * u = element            if parent(element)._manipulate_control_(element) is
                                 defined, then it will be used.  Otherwise make a
                                 text input that coerces input to parent(element)
   [ ] tag_cell('foo') -- makes it so one can refer to the current cell
       from elsewhere using the tag foo instead of the cell id number
       This involves doing something with SAGE_CELL_ID and the cell_id()
       method.
   [ ] 100% doctest coverage

JQUERY:
   [ ] tab_view -- represents an object in which clicking the tab
                     with label lbl[i] displays expr[i]
   [ ] slide_view -- represents an object in which the a list of objects
                     are displayed on successive slides.
   [ ] framed -- put a frame around an object
   [ ] panel -- put an object in a panel
   [ ] flot (?)

ELEMENTS:
   [ ] control: this models the input and other tags in html
          align -- left, right, top, texttop, middle, absmiddle, baseline, bottom, absbottom
          background -- the color of the background for the cell
          frame -- draw a frame around
          disabled -- disables the input element when it first loads
                      so that the user can not write text in it, or
                      select it.
          editable -- bool
          font_size -- integer
          maxlength -- the maximum number of characters allowed in a text field.
          name -- defines a unique name for the input element
          size -- the size of the input element
          type -- button, checkbox, file, password, radio, slider, text, setter_bar, drop_down
   [ ] setter bar (buttons)
   [ ] checkbox
   [ ] color slider
   [ ] blank input field
   [ ] 2d slider
   [ ] locator in a graphic

IDEAS for code:

@manipulate
def foo(x=range(10), y=slider(1,10)):
    ...

@manipulate
def foo(x=random_matrix(ZZ,2))
    ...


@framed
def foo(x,y):
    ...
"""

import inspect

SAGE_CELL_ID = 0
vars = {}

_k = 0
def new_adapt_name():
    global _k
    _k += 1
    return 'adapt%s'%_k


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

def html_slider(id, callback, margin=0):
    s = """<div id='%s' class='ui-slider-1' style="margin:%spx;"><span class='ui-slider-handle'></span></div>"""%(
        id, int(margin))

    # We now generat javascript that gets run after the above div gets
    # inserted. This happens because of the setTimeout function below
    # which gets passed an anonymous function.
    s += """
    <script>
    setTimeout(function() {
        $('#%s').slider();
        $('#%s').bind('click', function () { var position = $('#%s').slider('value',0); %s; });
    }, 1)
    </script>
    """%(id, id, id, callback)
    return s


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

    def adapt(self):
        """
        Return string representation of function that is called to
        adapt the values of this control to Python.
        """
        name = new_adapt_name()
        vars[name] = self._adapt
        return 'sage.server.notebook.manipulate.vars[\\"%s\\"]'%name

    def _adapt(self, x):
        return x

    def manipulate(self):
        """
        Return a string that when evaluated in Javascript calls the
        javascript manipulate function with appropriate inputs for
        this control.

        OUTPUT:
            string -- that is meant to be evaluated in Javascript
        """
        s = 'manipulate(%s, "sage.server.notebook.manipulate.vars[%s][\\"%s\\"]=sage_eval(r\\"\\"\\"%s("+%s+")\\"\\"\\", globals())\\n%s()");'%(
            self.cell_id(), self.cell_id(), self.var(), self.adapt(), self.value(), self.function_name())
        return s

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
        return 'this.value'

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
        return "position"

    def _adapt(self, position):
        # Input position is a string that evals to a float between 0 and 100.
        # we translate it into an index into self.__values
        v = self.__values
        i = int(len(v) * (float(position)/100.0))
        if i < 0:
            i = 0
        elif i >= len(v):
            i = len(v) - 1
        return v[i]

    def render(self):
        """
        Render this control as a string.

        OUTPUT:
             string -- html format
        """
        return html_slider('slider-%s-%s'%(self.var(), self.cell_id()), self.manipulate())


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
           <table border=0 width=800px height=500px>
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
    controls = [automatic_control(f, args[i], defaults[i-n] if i >= n else None) for i in range(len(args))]
    #controls = [automatic_control(f, v) for v in args[:n]] + \
    #           [defaults[i].render(f, args[i+n]) for i in range(len(defaults))]

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

def automatic_control(f, v, default):
    if isinstance(default, control):
        C = default
    elif isinstance(default, list):
        C = slider(default)
    else:
        C = text_box(str(default))
    return C.render(f, v)
