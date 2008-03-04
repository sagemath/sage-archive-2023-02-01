r"""
Manipulate Sage functions in the notebook

This module implements a manipulate decorator for function in the Sage
notebook.

The controls are:
\begin{itemize}
    \item InputBox -- a input box
    \item Slider -- a slider
\end{itemize}

AUTHOR:
    -- William Stein (2008-03-02): initial version at Sage/Enthought
       Days 8 in Texas
    -- Jason Grout (2008-03): collaborators substantially on the
       design and prototypes.

TODO:
   [X] get sliders to work; values after move slider

   [x] default values
   [x] get everything in the current version to work 100% bug free (including some style). post bundle.
       BUGS:
          [x] have default values set from the get go
          [x] spacing around sliders; also need to have labels
          [x] when re-evaluate input, make sure to clear output so cell-manipulate-id div is gone.
          [x] two manipulates in one cell -- what to do?
          [x] draw initial state
          [x] make manipulate canvas resizable
          [x] if you  use a manipulate control after restarting, doesn't work.   Need to reset it.  How?
                (to finish -- fix the error message in js.py:
                    /* TODO: Make error message more distinct! */
                    if (new_manip_output.indexOf('KeyError: '+id) != -1) {
                        evaluate_cell(id, 0);
                        new_manip_output = "";
                    }
          [x] display html parts of output as html

   [x] NO -- autoswitch to 1-cell mode:
           put                  slide_mode(); jump_to_slide(%s);   in wrap_in_outside_frame
        but feals all wrong.

   [ ] completely get rid of left clicking to switch wrap mode for
           manipulate objects: always in word wrap mode!

   [ ] implement a color object
   [ ] cool looking sliders:
        http://jqueryfordesigners.com/demo/slider-gallery.html

PLANS and IDEAS:
   [ ] automagically determine the type of control from the default
       value of the variable.  Here is how this will work:
        * u                      blank input field
        * u = (umin,umax)        slider; umin must not be a sequence
        * u = (umin,umax,du)     discrete slider
        * u = [1,2,3,4]          setter bar: automatically when there are are most 5
                                 elements; otherwise a drop down menu
        * u = ((xmin,xmax),(ymin,ymax))  2d slider
        * u = Graphic            a locator in a 2d plot  (Graphic is a 2d graphics objects)
        * u = True or u = False  a checkbox
        * u = color(...)         a color slider
        * u = "string"           input field
        * u = ('label', obj)     obj can be any of the above; control is labeled with
                                 the given label
        * u = element            if parent(element)._manipulate_control_(element) is
                                 defined, then it will be used.  Otherwise make an
                                 input that coerces input to parent(element)
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
   [ ] color slider:
          http://interface.eyecon.ro/demos/slider_colorpicker.html
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

from sage.misc.all import srange, sage_eval

import inspect


# Module scope variable that is always set equal to
# the current cell id (of the executing cell).

SAGE_CELL_ID = 0

# Dictionary that stores the state of all evaluated
# manipulate cells.
state = {}

_k = 0
def new_adapt_number():
    global _k
    _k += 1
    return _k


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

def html_slider(label, id, callback, steps, default=0, margin=0):
    s = """<table style='margin:0px;padding:0px;'><tr><td>%s</td><td><div id='%s' class='ui-slider-1' style='padding:0px;margin:%spx;'><span class='ui-slider-handle'></span></div></div></td></tr></table>"""%(
        label, id, int(margin))

    # We now generat javascript that gets run after the above div gets
    # inserted. This happens because of the setTimeout function below
    # which gets passed an anonymous function.
    s += """
    <script>
    setTimeout(function() {
        $('#%s').slider({
               stepping: 1, minValue: 0, maxValue: %s, startValue: %s,
               change: function () { var position = Math.ceil($('#%s').slider('value')); %s; }
        });
    }, 1);      /* setTimeout might be a hack? This could lead to a bug?  */
    </script>
    """%(id, steps-1, default, id, callback)
    return s


class ManipulateControl:
    """
    Base class for manipulate controls.
    """
    def __init__(self, f, var, default_value):
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
        self.__default_value = default_value
        self.__adapt_number = new_adapt_number()

    def __repr__(self):
        return "A ManipulateControl (abstract base class)"

    def default_value(self):
        return self.__default_value

    def adapt_number(self):
        """
        Return string representation of function that is called to
        adapt the values of this control to Python.
        """
        return self.__adapt_number

    def _adaptor(self, value, globs):
        """
        Adapt a user input, which is a string, to be an element selected
        by this control.

        INPUT:
            value -- the string the user typed in
            globs -- the globals interpreter variables, e.g.,
                     globals(), which is useful for evaling value.

        OUTPUT:
            object
        """
        return sage_eval(value, globs)

    def manipulate(self):
        """
        Return a string that when evaluated in Javascript calls the
        javascript manipulate function with appropriate inputs for
        this control.

        OUTPUT:
            string -- that is meant to be evaluated in Javascript
        """
        # The following is a crazy line to read because of all the backslashes and try/except.
        # All it does is run the manipulate function once after setting exactly one
        # dynamic variable.    If setting the dynamic variable fails, due to a KeyError
        s = 'manipulate(%s, "sage.server.notebook.manipulate.update(%s, \\"%s\\", %s, \\""+%s+"\\", globals())")'%(
            self.cell_id(), self.cell_id(), self.var(), self.adapt_number(), self.value_js())
        return s

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

class InputBox(ManipulateControl):
    """
    An input box manipulate control.
    """
    def __repr__(self):
        return "A InputBox manipulate control"

    def value_js(self):
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
        %s: <input type='text' value='%r' width=200px onchange='%s'></input>
        """%(self.var(), self.default_value(),  self.manipulate())

class Slider(ManipulateControl):
    """
    A slider manipulate control.
    """
    def __init__(self, f, var, values, default_position):
        """
        Create a slider manipulate control that takes on the given
        list of values.
        """
        ManipulateControl.__init__(self, f, var, values[default_position])
        self.__values = values
        self.__default_position = default_position

    def __repr__(self):
        return "A Slider manipulate control"

    def default_position(self):
        """
        Return the default position (as an integer) of the slider.
        """
        return self.__default_position

    def value_js(self):
        """
        Return javascript string that will give the
        value of this control element.

        OUTPUT:
             string -- javascript
        """
        return "position"

    def _adaptor(self, position, globs):
        """
        Adapt a user input, which is the slider position, to be an
        element selected by this control.

        INPUT:
            position -- position of the slider
            globs -- the globals interpreter variables (not used here).

        OUTPUT:
            object
        """
        v = self.__values
        # We have to cast to int, since it comes back as a float that
        # is too big.
        return v[int(position)]

    def render(self):
        """
        Render this control as a string.

        OUTPUT:
             string -- html format
        """
        return html_slider('<font color=black>%s</font> '%self.var(), 'slider-%s-%s'%(self.var(), self.cell_id()),
                           self.manipulate(), steps=len(self.__values),
                           default=self.default_position())


class ManipulateCanvas:
    """
    Base class for manipulate canvases. This is where all the controls
    along with the output of the manipulated function are layed out
    and rendered.
    """
    def __init__(self, controls, id):
        """
        Create a manipulate canvas.

        INPUT:
            controls -- a list of ManipulateControl instances.
        """
        self.__controls = controls
        self.__cell_id = id

    def cell_id(self):
        return self.__cell_id

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
        <div id='cell-manipulate-%s'><?START>
        <table border=0 bgcolor='#white' width=100%% height=100%%>
        <tr><td bgcolor=white align=center valign=top>
          <?TEXT>
        </td></tr>
        <tr><td  align=center valign=top><?HTML></td></tr>
        </table>
        </td></tr></table><?END></div>
        """%self.cell_id()

    def render_controls(self):
        """
        Render in text (html) form all the input controls.

        OUTPUT:
            string -- html
        """
        # This will need some sophisticated layout querying of the c's, maybe.
        return ''.join([c.render() for c in self.__controls])

    def wrap_in_outside_frame(self, inside):
        return """<div padding=6 id='div-manipulate-%s'> <table bgcolor='#c5c5c5'
                 width=600px height=400px cellpadding=15><tr><td bgcolor='#f9f9f9' valign=top align=center>%s</td>
                 </tr></table></div>
                 """%(self.cell_id(), inside)
##                  <script>
##                  setTimeout(function() {
##                  $('#div-manipulate-%s').resizable();
##                  $('#div-manipulate-%s').draggable();
##                  }, 1);</script>

    def render(self):
        """
        Render in text (html) the entire manipulate canvas.

        OUTPUT:
            string -- html
        """
        s = "%s%s"%(self.render_controls(), self.render_output())
        s = self.wrap_in_outside_frame(s)
        return s


def manipulate(f):
    """
    Decorate a function f to make a manipulate version of f.
        @manipulate
        def foo(n,m):
            ...

    Note -- it is safe to make several distinct manipulate cells with functions that
    have the same name.
    """

    (args, varargs, varkw, defaults) = inspect.getargspec(f)

    if defaults is None:
        defaults = []

    n = len(args) - len(defaults)
    controls = [automatic_control(f, args[i], defaults[i-n] if i >= n else None) for i in range(len(args))]

    C = ManipulateCanvas(controls, SAGE_CELL_ID)

    d = {}
    ad = {}
    state[SAGE_CELL_ID] = {'variables':d, 'adapt':ad}

    for con in controls:
        d[con.var()] = con.default_value()
        ad[con.adapt_number()] = con._adaptor

    html(C.render())

    def g():
        return f(*[d[args[i]] for i in range(len(args))])

    state[SAGE_CELL_ID]['function'] = g

    return g


######################################################
# Actual control objects that the user passes in
######################################################
class control:
    pass

class input_box(control):
    def __init__(self, default):
        """
        INPUT:
            default -- string (the default value)
        """
        self.__default = default

    def __repr__(self):
        return "A manipulate input box control with default value '%r'"%self.__default

    def render(self, f, var):
        """
        INPUT:
            f -- a Python function
            var -- a string (variable; one of the variable names input to f)
        """
        return InputBox(f, var, self.__default)

class slider(control):
    def __init__(self, vmin, vmax=None, steps=30, default=None):
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
        if len(self.__values) == 0:
            self.__values = [0]
        if default is None:
            self.__default = 0
        else:
            try:
                i = self.__values.index(default)
            except ValueError:
                i = 0
            self.__default = i

    def __repr__(self):
        return "A manipulate slider control [%s - %s]."%(min(self.__values), max(self.__values))

    def default(self):
        """
        Return default index into the list of values.

        OUTPUT:
            int
        """
        return self.__default

    def render(self, f, var):
        return Slider(f, var, self.__values, self.__default)

def automatic_control(f, v, default):
    if isinstance(default, control):
        C = default
    elif isinstance(default, list):
        C = slider(default)
    elif isinstance(default, tuple):
        if len(default) == 2:
            if isinstance(default[0], list):
                C = slider(default[0], default=default[1])
            else:
                C = slider(range(default[0], default[1]))
        elif len(default) == 3:
            C = slider(srange(default[0], default[1], default[2]))
        else:
            C = slider(list(default))
    else:
        C = input_box(default)
    return C.render(f, v)


def update(cell_id, var, adapt, value, globs):
    """
    INPUT:
        cell_id -- the id of a manipulate cell
        var -- a variable associated to that cell
        adapt -- the number of the adapt function
    """
    try:
        S = state[cell_id]
    except KeyError:
        print "__SAGE_MANIPULATE_RESTART__"
    else:
        # Look up the function that adapts inputs to have the right type
        adapt_function = S["adapt"][adapt]
        # Apply that function and save the result in the appropriate variables dictionary.
        S["variables"][var] = adapt_function(value, globs)
        # Finally call the manipulatable function, which will use the above variables.
        S['function']()
