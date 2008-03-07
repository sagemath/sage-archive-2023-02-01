#############################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

r"""
Interact Sage functions in the notebook

This module implements a interact decorator for function in the Sage
notebook.

AUTHORS:
    -- William Stein (2008-03-02): version 1.0 at Sage/Enthought Days 8 in Texas
    -- Jason Grout (2008-03): discussion, and first few prototypes

"""

"""
 ** PLANNING **

BUGS:
   [x] have default values set from the get go
   [x] spacing around sliders; also need to have labels
   [x] when re-evaluate input, make sure to clear output so cell-interact-id div is gone.
   [x] two interacts in one cell -- what to do?
   [x] draw initial state
   [x] make interact canvas resizable
   [x] if you  use a interact control after restarting, doesn't work.   Need to reset it.  How?
   [x] display html parts of output as html
   [x] default slider pos doesn't work, eg. def cont(q1=(-1,(-3,3)), q2=(1,(-3,3))):
   [x] change from interact to interact everywhere.

-->  [ ] edit/save breaks interact mode <--
          * picking up images that shouldn't.
          * controls completely stop working.

   [ ] tab completion in interact broken formating
   [ ] error exception reporting broken
   [ ] problems with html/pre/text formating, e.g., in TEXT mode and in interact cells
   [ ] replace special %interact by something very obfuscated to keep from having
       really weird mistakes that are hard for people to debug.
   [ ] flicker resize during update (hard???)
   [ ] cross-platform testing.
   [ ] cell order corruption
   [ ] slider is too narrow -- need to expand to window width


VERSION 1:
   [X] get sliders to work; values after move slider
   [x] default values
   [x] NO -- autoswitch to 1-cell mode:
           put                  slide_mode(); jump_to_slide(%s);   in wrap_in_outside_frame
        but feals all wrong.
   [x] completely get rid of left clicking to switch wrap mode for
           interact objects: always in word wrap mode or hide.
   [ ] shortcut ('label', v)
   [ ] 100% documentation and doctest coverage
   [ ] test saving and loading whole notebook to a file
   [ ] setter bar (buttons)
   [ ] drop down menu
   [ ] checkbox
   [ ] write several paragraps at the top of this file that describe
       the philosophy and use of these interact things
   [ ] put the docs for this in the reference manul.
   [ ] collection of about 20 good examples of how to use interact (doctested)
   [ ] make it so slider has width 100%.
   [ ] line up all the control in a single table so all labels and all
       controls exactly match up

VERSION 2:
   [ ] make option for text input that correctly gives something of
       the same type as the default input.
   [ ] matrix input control (method of matrix space) -- a spreadsheet like thing
   [ ] a 2d slider:
          u = ((xmin,xmax),(ymin,ymax))  2d slider   -- NOT IMPLEMENTED
   [ ] locator in a 2d graphic
   [ ] tab_view -- represents an object in which clicking the tab
                   with label lbl[i] displays expr[i]
   [ ] controls: make them easy to customize as below --
          location -- where to put the slider (?)
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

VERSION 3:
   [ ] protocol for objects to have their own interact function; make
       it so for any object obj in sage, one can do
             {{{
             interact(obj)
             }}}
       and get something useful as a result.
   [ ] flot -- some pretty but very simple javascript only graphics  (maybe don't do... unclear)
   [ ] zorn -- similar (maybe don't do... unclear)
   [ ] color sliders (?):
          u = color(...)         a color slider
          http://interface.eyecon.ro/demos/slider_colorpicker.html
          http://jqueryfordesigners.com/demo/slider-gallery.html
   [ ] tag_cell('foo') -- makes it so one can refer to the current cell
       from elsewhere using the tag foo instead of the cell id number
       This involves doing something with SAGE_CELL_ID and the cell_id() method.
   [ ] framed -- put a frame around an object


"""

import inspect

from sage.misc.all import srange, sage_eval

# SAGE_CELL_ID is a module scope variable that is always set equal to
# the current cell id (of the executing cell).  Code that sets this is
# inserted by the notebook server each time a worksheet cell is
# evaluated.
SAGE_CELL_ID = 0

# Dictionary that stores the state of all active interact cells.
state = {}

_k = 0
def new_adapt_number():
    """
    Return an integer, always counting up, and starting with 0.  This
    is used for saving the adapt methods for controls.  An adapt
    method is just a function that coerces data into some object,
    e.g., makes sure the control always produces int's.

    OUTPUT:
        integer

    EXAMPLES:
        sage: sage.server.notebook.interact.new_adapt_number()   # random output -- depends on when called
        1
    """
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

    EXAMPLES:
        sage: sage.server.notebook.interact.html('hello')
        <html>hello</html>
    """
    print "<html>%s</html>"%s

def html_slider(label, id, callback, steps, default=0, margin=0):
    """
    Return the HTML representation of a jQuery slider.

    INPUT:
        label   -- string that is inserted just to the left of the slider
        id      -- string -- the DOM id of the slider (better be unique)
        callback-- javascript that is executed whenever the slider is done moving
        steps   -- number of steps from minimum to maximum value.
        default -- (default: 0) the default position of the slider
        margin  -- (default: 0) size of margin to insert around the slider

    EXAMPLES:
    We create a jQuery HTML slider.    If you do the following in the notebook
    you should obtain a slider that when moved pops up a window showing its
    current position.
        sage: from sage.server.notebook.interact import html_slider, html
        sage: html(html_slider('my slider', 'slider-007', 'alert(position)', steps=5, default=2, margin=5))
        <html><table style='margin:0px;padding:0px;'><tr><td>my slider</td><td><div id='slider-007' class='ui-slider-1' style='padding:0px;margin:5px;'><span class='ui-slider-handle'></span></div></div></td></tr></table><script>setTimeout(function() { $('#slider-007').slider({
        stepping: 1, minValue: 0, maxValue: 4, startValue: 2,
        change: function () { var position = Math.ceil($('#slider-007').slider('value')); alert(position); }
        });}, 1);</script></html>
    """
    s = """<table style='margin:0px;padding:0px;'><tr><td>%s</td><td><div id='%s' class='ui-slider-1' style='padding:0px;margin:%spx;'><span class='ui-slider-handle'></span></div></div></td></tr></table>"""%(
        label, id, int(margin))

    # We now generate javascript that gets run after the above div
    # gets inserted. This happens because of the setTimeout function
    # below which gets passed an anonymous function.
    s += """<script>setTimeout(function() { $('#%s').slider({
    stepping: 1, minValue: 0, maxValue: %s, startValue: %s,
    change: function () { var position = Math.ceil($('#%s').slider('value')); %s; }
});}, 1);</script>"""%(id, steps-1, default, id, callback)

    return s


class InteractControl:
    """
    Base class for interact controls.  These are controls associated
    with a specific interact instance.  They have internal state
    information about the specific function being interactd, etc.
    """
    def __init__(self, f, var, default_value, label=None):
        """
        Create a new interact control.

        INPUT:
             f -- a Python function (that's being decorated)
             var -- string; name of variable that this control interacts
             default_value -- the default value of the variable
                              corresponding to this control.
             label -- string (default: None) label of this control; if None
                      then defaults to var.


        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: def f(x): print x
            ...
            sage: InteractControl(f, 'x', 5)
            A InteractControl (abstract base class)
        """
        self.__var = var
        self.__cell_id = SAGE_CELL_ID
        self.__f = f
        self.__default_value = default_value
        self.__adapt_number = new_adapt_number()
        if label is None:
            self.__label = var
        else:
            self.__label = label

    def __repr__(self):
        """
        String representation of interact control.

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: def f(x): print x
            ...
            sage: InteractControl(f, 'x', 5).__repr__()
            'A InteractControl (abstract base class)'
        """
        return "A InteractControl (abstract base class)"

    def value_js(self):
        """
        Javascript that when evaluated gives the current value of this
        control.  This should be redefined in a derived class.

        OUTPUT:
            string -- defaults to NULL -- this should be redefined.
        """
        return 'NULL'

    def label(self):
        """
        Return the text label of this interact control.

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: def f(x): print x
            ...
            sage: InteractControl(f, 'x', 5, label='the x value').label()
            'the x value'
        """
        return self.__label

    def default_value(self):
        """
        Return the default value of the variable corresponding to this
        interact control.

        OUTPUT:
            object

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: def f(x): print x
            ...
            sage: InteractControl(f, 'x', 19/3).default_value()
            19/3
        """
        return self.__default_value

    def adapt_number(self):
        """
        Return integer index into adapt dictionary of function that is
        called to adapt the values of this control to Python.

        OUTPUT:
            an integer

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: def f(x): print x
            ...
            sage: InteractControl(f, 'x', 19/3).adapt_number()       # random -- depends on call order
            2
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

        EXAMPLES:
            sage: sage.server.notebook.interact.InteractControl(math.sin, 'x', 1)._adaptor('2/3', globals())
            2/3
        """
        return sage_eval(value, globs)

    def interact(self):
        """
        Return a string that when evaluated in Javascript calls the
        javascript interact function with appropriate inputs for
        this control.

        OUTPUT:
            string -- that is meant to be evaluated in Javascript

        EXAMPLES:
            sage: sage.server.notebook.interact.InteractControl(math.cos, 'x', 1).interact()
            'interact(0, "sage.server.notebook.interact.update(0, \\"x\\", 5, \\""+NULL+"\\", globals())")'
        """
        # The following is a crazy line to read because of all the backslashes and try/except.
        # All it does is run the interact function once after setting exactly one
        # dynamic variable.    If setting the dynamic variable fails, due to a KeyError
        s = 'interact(%s, "sage.server.notebook.interact.update(%s, \\"%s\\", %s, \\""+%s+"\\", globals())")'%(
            self.cell_id(), self.cell_id(), self.var(), self.adapt_number(), self.value_js())
        return s

    def var(self):
        """
        Return the name of the variable that this control interacts.

        OUTPUT:
            string -- name of a variable as a string.

        EXAMPLES:
            sage: sage.server.notebook.interact.InteractControl(math.cos, 'theta', 1).var()
            'theta'
        """
        return self.__var

    def cell_id(self):
        """
        Return the id of the cell that contains this interact control.

        OUTPUT:
            integer -- id of cell that this control interacts

        EXAMPLES:
        The output below should equal the ID of the current cell.
            sage: sage.server.notebook.interact.InteractControl(math.cos, 'theta', 1).cell_id()
            0
        """
        return self.__cell_id

class InputBox(InteractControl):
    """
    An input box interact control.
    EXAMPLES:
        sage: sage.server.notebook.interact.InputBox(math.cos, 'theta', 1)
        A InputBox interact control
    """
    def __repr__(self):
        """
        String representation of an InputBox interactive control.

        EXAMPLES:
            sage: sage.server.notebook.interact.InputBox(math.cos, 'theta', 1).__repr__()
            'A InputBox interact control'
        """
        return "A InputBox interact control"

    def value_js(self):
        """
        Return javascript string that will give the
        value of this control element.

        OUTPUT:
             string -- javascript

        EXAMPLES:
            sage: sage.server.notebook.interact.InputBox(math.cos, 'theta', 1).value_js()
            'this.value'
        """
        return 'this.value'

    def render(self):
        """
        Render this control as a string.

        OUTPUT:
             string -- html format

        EXAMPLES:
            sage: sage.server.notebook.interact.InputBox(math.cos, 'theta', 1).render()
            '\n        theta: <input type=\'text\' value=\'1\' width=200px onchange=\'interact(..., "sage.server.notebook.interact.update(..., \\"theta\\", ..., \\""+this.value+"\\", globals())")\'></input>\n        '
        """
        return """
        %s: <input type='text' value='%r' width=200px onchange='%s'></input>
        """%(self.label(), self.default_value(),  self.interact())

class Slider(InteractControl):
    """
    A slider interact control.
    """
    def __init__(self, f, var, values, default_position, label=None):
        """
        Create a slider interact control that takes on the given
        list of values.

        INPUT:
            f -- function
            var -- string; name of variable being interactd
            values -- list; a list of the values that the slider will take on
            default_position -- int; default location that the slider is set to.
            label -- alternative label to the left of the slider,
                     instead of the variable.

        EXAMPLES:
        """
        InteractControl.__init__(self, f, var, values[default_position], label=label)
        self.__values = values
        self.__default_position = default_position

    def __repr__(self):
        return "Slider Interact Control: %s [%s--|%s|---%s]."%(
            self.label(), min(self.__values),
            self.default_value(), max(self.__values))

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
        return html_slider('<font color=black>%s</font> '%self.label(), 'slider-%s-%s'%(self.var(), self.cell_id()),
                           self.interact(), steps=len(self.__values),
                           default=self.default_position())


class InteractCanvas:
    """
    Base class for interact canvases. This is where all the controls
    along with the output of the interactd function are layed out
    and rendered.
    """
    def __init__(self, controls, id):
        """
        Create a interact canvas.

        INPUT:
            controls -- a list of InteractControl instances.
        """
        self.__controls = controls
        self.__cell_id = id

    def cell_id(self):
        return self.__cell_id

    def render_output(self):
        """
        Render in text (html) form the output portion of the interact canvas.

        The output contains two special tags, <?TEXT> and <?HTML>,
        which get replaced at runtime by the text and html parts
        of the output of running the function.

        OUTPUT:
            string -- html
        """
        return """
        <div id='cell-interact-%s'><?START>
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
        return """<div padding=6 id='div-interact-%s'> <table width=800px height=400px bgcolor='#c5c5c5'
                 cellpadding=15><tr><td bgcolor='#f9f9f9' valign=top align=left>%s</td>
                 </tr></table></div>
                 """%(self.cell_id(), inside)
    # The following could be used to make the interact frame resizable and/or draggable.
    # Neither effect is as cool as it sounds!
##                  <script>
##                  setTimeout(function() {
##                  $('#div-interact-%s').resizable();
##                  $('#div-interact-%s').draggable();
##                  }, 1);</script>

    def render(self):
        """
        Render in text (html) the entire interact canvas.

        OUTPUT:
            string -- html
        """
        s = "%s%s"%(self.render_controls(), self.render_output())
        s = self.wrap_in_outside_frame(s)
        return s


def interact(f):
    r"""
    Make a function with controls to enable interactivity with the
    variables.

    SOME EXAMPLES:
    We create a interact control with two inputs, a text input for
    the variable $a$ and a $y$ slider that runs through the range of
    integers from $0$ to $19$.
        sage: @interact
        ... def f(a=5, y=range(20)): print a + y
        ...
        <html>...

    Draw a plot interacting with the ``continuous'' variable $a$.  By
    default continuous variables have exactly 50 possibilities.
        sage: @interact
        ... def example(a=(0,2)):
        ...     show(plot(sin(x*(1+a*x)), (x,0,6)), figsize=4)
        <html>...

    Interact a variable in steps of 1 (we also use an unnamed
    function):
        sage: @interact
        ... def _(n=(10,100,1)):
        ...     show(factor(x^n - 1))
        <html>...

    Interact two variables:
        sage: @interact
        ... def _(a=(1,4), b=(0,10)):
        ...     show(plot(sin(a*x+b), (x,0,6)), figsize=3)
        <html>...

    DEFAULTS:
    Defaults for the variables of the input function determine
    interactive controls according to the following rules:
    \begin{itemize}
        \item u -- blank input field
        \item u = element -- input eval field with default = element
        \item u = (umin,umax) -- continuous slider (really 50 steps)
        \item u = (umin,umax,du) -- slider with step size du
        \item u = [1,2,3] -- buttons if at most 5; otherwise, drop down
        \item u = a bool -- a checkbox
        \item u = "string" -- input field that allows arbitrary strings
        \item u = ('label', v) -- labeled control with v any of these
        \item u = (default, (umin, umax[, du]))
    \end{itemize}

    MORE EXAMPLES:
    We give defaults and name the variables:
        sage: @interact
        ... def _(a=('first', (1,4)), b=(0,10)):
        ...       show(plot(sin(a*x+sin(b*x)), (x,0,6)), figsize=3)
        <html>...

    Another example involving labels and defaults, and the
    slider command.
        sage: @interact
        ... def _(a = slider(1, 4, default=2, label='Multiplier'),
        ...       b = slider(0, 10, default=0, label='Phase Variable')):
        ...     show(plot(sin(a*x+b), (x,0,6)), figsize=4)
        <html>...

    You can rotate and zoom into three graphics while
    interacting with a variable.
        sage: @interact
        ... def _(a=(0,1)):
        ...     x,y = var('x,y')
        ...     show(plot3d(sin(x*cos(y*a)), (x,0,5), (y,0,5)), figsize=4)
        <html>...

    A random polygon:
        sage: pts = [(random(), random()) for _ in xrange(20)]
        sage: @interact
        ... def _(n = [4..len(pts)]):
        ...     G = points(pts[:n],pointsize=60) + polygon(pts[:n], rgbcolor='black')
        ...     show(G, figsize=5, xmin=0, ymin=0)
        <html>...

    Two "sinks" displayed simultaneously via a contour plot and a 3d
    interactive plot:
        sage: @interact
        ... def cont(q1=(-1,(-3,3)), q2=(-2,(-3,3))):
        ...     x,y = var('x,y')
        ...     f = q1/sqrt((x+1)^2 + y^2) + q2/sqrt((x-1)^2+(y+0.5)^2)
        ...     g = f._fast_float_('x','y')   # should not be needed soon
        ...     C = contour_plot(g, (-2,2), (-2,2), plot_points=30, contours=15, cmap='cool')
        ...     show(C, figsize=3, aspect_ratio=1)
        ...     show(plot3d(g, (-2,2), (-2,2)), figsize=4)
        <html>...

    In the following example, we only generate data for a given n
    once, so that as one varies p the data doesn't not randomly
    change.  We do this by simply caching the results for each n
    in a dictionary.
        sage: data = {}
        sage: @interact
        ... def _(n=(500,(100,5000,1)), p=(1,(0.1,10))):
        ...     n = int(n)
        ...     if not data.has_key(n):
        ...         data[n] = [(random(), random()) for _ in xrange(n)]
        ...     show(points([(x^p,y^p) for x,y in data[n]], rgbcolor='black'), xmin=0, ymin=0, axes=False)
        <html>...

    A conchoid:
        sage: @interact
        ... def _(k=(1.2,(1.1,2)), k_2=(1.2,(1.1,2)), a=(1.5,(1.1,2))):
        ...     u, v = var('u,v')
        ...     f = (k^u*(1+cos(v))*cos(u), k^u*(1+cos(v))*sin(u), k^u*sin(v)-a*k_2^u)
        ...     show(parametric_plot3d(f, (u,0,6*pi), (v,0,2*pi), plot_points=[40,40], texture=(0,0.5,0)))
        <html>...
    """

    (args, varargs, varkw, defaults) = inspect.getargspec(f)

    if defaults is None:
        defaults = []

    n = len(args) - len(defaults)
    controls = [automatic_control(defaults[i-n] if i >= n else None) for i in range(len(args))]

    # Convert the controls to InteractControl objects
    controls = [controls[i].render(f, args[i]) for i in range(len(args))]

    C = InteractCanvas(controls, SAGE_CELL_ID)

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
    def __init__(self, label=None):
        self.__label = label

    def label(self):
        return self.__label

    def set_label(self, label):
        self.__label = label

class input_box(control):
    def __init__(self, default, label=None):
        """
        INPUT:
            default -- string (the default value)
        """
        self.__default = default
        control.__init__(self, label)

    def __repr__(self):
        return "A interact input box control with default value '%r'"%self.__default

    def render(self, f, var):
        """
        INPUT:
            f -- a Python function
            var -- a string (variable; one of the variable names input to f)
        """
        return InputBox(f, var, self.__default, label=self.label())

class slider(control):
    def __init__(self, vmin, vmax=None, steps=30, default=None, label=None):
        control.__init__(self, label=label)
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

        # determine the best choice of index into the list of values
        # for the user-selected default.
        if default is None:
            self.__default = 0
        else:
            try:
                i = self.__values.index(default)
            except ValueError:
                # here no index matches -- which is best?
                try:
                    v = [(abs(default - self.__values[j]), j) for j in range(len(self.__values))]
                    m = min(v)
                    i = m[1]
                except TypeError: # abs not defined on everything, so give up
                    i = 0
            self.__default = i

    def __repr__(self):
        return "Slider: %s [%s--|%s|---%s]."%(self.label(),
                  min(self.__values),
             self.__values[self.default_index()], max(self.__values))

    def default_index(self):
        """
        Return default index into the list of values.

        OUTPUT:
            int
        """
        return self.__default

    def render(self, f, var):
        """
        Render the interact control for the given function and
        variable.

        INPUT:
            f -- a function
            var -- string; variable name

        EXAMPLES:
            sage: S = slider(0,10, default=3, label='theta'); S
            Slider: theta [0--|3|---10].
            sage: S.render(lambda x: x^2, 'x')
            Slider Interact Control: theta [0--|3|---10].
        """
        return Slider(f, var, self.__values, self.__default, label=self.label())

def automatic_control(default):
    """
    Automagically determine the type of control from the default
    value of the variable.

    INPUT:
        default -- the default value for v given by the function

    OUTPUT:
        a interact control
    """
    label = None
    default_value = None
    if isinstance(default, control):
        C = default
        label = C.label()
    elif isinstance(default, list):
        C = slider(default)
    elif isinstance(default, tuple):
        if len(default) == 2:
            if isinstance(default[1], tuple):
                default_value = default[0]
                default = default[1]
            if isinstance(default[0], str):
                # recursive call to get automatic control for second part
                C = automatic_control(default[1])
                label = default[0]
            elif isinstance(default[0], list):
                C = slider(default[0], default=default[1])
            else:
                # The default 49.0 below is a sort of "heuristic value" so there are 50 steps
                C = slider(srange(default[0], default[1], (default[1]-default[0])/49.0,
                                  include_endpoint=True), default = default_value)
        elif len(default) == 3:
            C = slider(srange(default[0], default[1], default[2], include_endpoint=True))
        else:
            C = slider(list(default))
    else:
        C = input_box(default)

    C.set_label(label)
    return C


def update(cell_id, var, adapt, value, globs):
    """
    INPUT:
        cell_id -- the id of a interact cell
        var -- a variable associated to that cell
        adapt -- the number of the adapt function
    """
    try:
        S = state[cell_id]
    except KeyError:
        print "__SAGE_INTERACT_RESTART__"
    else:
        # Look up the function that adapts inputs to have the right type
        adapt_function = S["adapt"][adapt]
        # Apply that function and save the result in the appropriate variables dictionary.
        S["variables"][var] = adapt_function(value, globs)
        # Finally call the interactive function, which will use the above variables.
        S['function']()
