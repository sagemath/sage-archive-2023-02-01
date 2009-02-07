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
    -- Jason Grout (2008-03): discussion and first few prototypes
    -- Jason Grout (2008-05): input_grid control
"""

"""
 ** PLANNING **

NOTES:
   * There is no testing of pickling anywhere in this file.  This is
     because there is no reason one would ever pickle anything in this
     file, since everything is associated with particular state
     information of a notebook.

BUGS:
   [x] have default values set from the get go
   [x] spacing around sliders; also need to have labels
   [x] when re-evaluate input, make sure to clear output so cell-interact-id div is gone.
   [x] two interacts in one cell -- what to do?
   [x] draw initial state
   [x] make interact canvas resizable
   [x] if you  use a interact control after restarting, doesn't work.   Need to reset it.  How?
   [x] display html parts of output as html
   [x] default slider pos doesn't work, eg. def _(q1=(-1,(-3,3)), q2=(1,(-3,3))):
   [x] change from interact to interact everywhere.
   [x] edit/save breaks interact mode
          * picking up images that shouldn't.
          * controls completely stop working.
   [x] problems with html/pre/text formating, e.g., in TEXT mode and in interact cells
   [x] tab completion in interact broken formating
   [x] error exception reporting broken
   [x] replace special %interact by something very obfuscated to keep from having
       really weird mistakes that are hard for people to debug.
   [x] cell order corruption
   [x] cross-platform testing (good enough -- it's jquery)
   [x] can't enter "foo" in input_box now because of how strings are
       passed back and forth using single quotes.
   [x] possible issue with page title being undefined; don't know why
       or if that is connected with interactives
   [x] autorunning interact cells on load is being injectected into the
       i/o pexpect stream way too early.
   [x] what do published worksheets do??

VERSION 1:
   [X] get sliders to work; values after move slider
   [x] default values
   [x] NO -- autoswitch to 1-cell mode:
           put                  slide_mode(); jump_to_slide(%s);   in wrap_in_outside_frame
        but feals all wrong.
   [x] completely get rid of left clicking to switch wrap mode for
           interact objects: always in word wrap mode or hide.
   [x] shortcut ('label', v)
   [x] test saving and loading whole notebook to a file
   [x] collection of about 20 good examples of how to use interact (doctested)
   [x] interact(f) should also work; i.e., no need to use decorators -- done; won't be advertised, but
       at least fixing this improves code quality.
   [x] obfuscate ?START and ?END much more.
   [x] type checked input box
   [x] line up all the control in a single table so all labels and all
       controls exactly match up
   [x] button bar
   [x] drop down menu
   [x] checkbox
   [x] color selector
   [x] something to avoid server flood, e.g., any %interact request removes all other
       such requests from the queue in worksheet.py

   DOCS:
   [x] 100% documentation and doctest coverage
   [ ] put the docs for this in the reference manual
   [ ] put summary doc in notebook help page

VERSION 2:
   [ ] vertical scroll bars (maybe very easy via jsquery)
   [ ] small version of color selector
   [ ] button -- block of code gets run when it is clicked
   [ ] when click a button in a button bar it is highlighted and
       other buttons are not (via some sort of javascript)
   [ ] much more fine control over style of all controls
   [ ] safe/secure evaluation mode
   [ ] slider is too narrow -- need to expand to window width?
   [ ] fix the flicker resize during update (hard???)
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
          type -- button, checkbox, file, password, radio, slider, text, setter_bar

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

# Standard system libraries
from base64 import standard_b64encode, standard_b64decode
import inspect
import math
import types

# Sage libraries
from sage.misc.all import srange, sage_eval
from sage.plot.misc import Color
from sage.structure.element import is_Matrix

# SAGE_CELL_ID is a module scope variable that is always set equal to
# the current cell id (of the executing cell).  Code that sets this is
# inserted by the notebook server each time a worksheet cell is
# evaluated.
SAGE_CELL_ID = 0

# Dictionary that stores the state of all active interact cells.
state = {}

def reset_state():
    """
    Reset the interact state of this sage process.

    EXAMPLES:
        sage: sage.server.notebook.interact.state  # random output
        {1: {'function': <function g at 0x72aaab0>, 'variables': {'m': 3, 'n': 5}, 'adapt': {1: <bound method Slider._adaptor of Slider Interact Control: n [1--|1|---10].>, 2: <bound method Slider._adaptor of Slider Interact Control: m [1--|1|---10].>}}}
        sage: from sage.server.notebook.interact import reset_state
        sage: reset_state()
        sage: sage.server.notebook.interact.state
        {}
    """
    global state
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

def html_slider(id, values, callback, steps, default=0, margin=0):
    """
    Return the HTML representation of a jQuery slider.

    INPUT:
        id      -- string -- the DOM id of the slider (better be unique)
        values  -- 'null' or javascript string containing array of values on slider
        callback-- javascript that is executed whenever the slider is done moving
        steps   -- number of steps from minimum to maximum value.
        default -- (default: 0) the default position of the slider
        margin  -- (default: 0) size of margin to insert around the slider

    EXAMPLES:
    We create a jQuery HTML slider.    If you do the following in the notebook
    you should obtain a slider that when moved pops up a window showing its
    current position.
        sage: from sage.server.notebook.interact import html_slider, html
        sage: html(html_slider('slider-007', 'null', 'alert(position)', steps=5, default=2, margin=5))
        <html>...</html>
    """
    s = """<table><tr><td>
    	<div id='%s' class='ui-slider ui-slider-3' style='margin:%spx;'><span class='ui-slider-handle'></span></div>
    	</td>"""%(id,int(margin))
    if values != "null":
        s += "<td><font color='black' id='%s-lbl'></font></td>"%id
    s += "</tr></table>"

    # We now generate javascript that gets run after the above div
    # gets inserted. This happens because of the setTimeout function
    # below which gets passed an anonymous function.
    s += """<script>(function(){ var values = %(values)s; setTimeout(function() {
    $('#%(id)s').slider({
    	stepping: 1, min: 0, max: %(maxvalue)s, startValue: %(startvalue)s,
    	change: function (e,ui) { var position = ui.value; if(values!=null) $('#%(id)s-lbl').text(values[position]); %(callback)s; },
    	slide: function(e,ui) { if(values!=null) $('#%(id)s-lbl').text(values[ui.value]); }
    });
    if(values != null) $('#%(id)s-lbl').text(values[$('#%(id)s').slider('value')]);
    }, 1); })();</script>"""%{'values': values, 'id': id, 'maxvalue': steps-1, 'startvalue': default, 'callback': callback}
    # change 'change' to 'slide' and it changes the slider every time it moves;
    # needs much more work to actually work, since server gets flooded by
    # requests.

    return s

def html_rangeslider(id, values, callback, steps, default_l=0, default_r=1, margin=0):
    """
    Return the HTML representation of a jQuery range slider.

    INPUT:
        id      -- string -- the DOM id of the slider (better be unique)
        values  -- 'null' or javascript string containing array of values on slider
        callback-- javascript that is executed whenever the slider is done moving
        steps   -- number of steps from minimum to maximum value.
        default_l -- (default: 0) the default position of the left edge of the slider
        default_r -- (default: 1) the default position of the right edge of the slider
        margin  -- (default: 0) size of margin to insert around the slider

    EXAMPLES:
    We create a jQuery range slider. If you do the following in the notebook
    you should obtain a slider that when moved pops up a window showing its
    current position.
        sage: from sage.server.notebook.interact import html_rangeslider, html
        sage: html(html_rangeslider('slider-007', 'null', 'alert(pos[0]+", "+pos[1])', steps=5, default_l=2, default_r=3, margin=5))
        <html>...</html>
    """
    s = """<table>
    <tr><td><div id='%s' class='ui-slider ui-slider-3' style='margin:%spx;'>
    <span class='ui-slider-handle'></span><span class='ui-slider-handle'></span>
    </div></td></tr>"""%(id,int(margin))
    if values != "null":
        s += "<tr><td><font color='black' id='%s-lbl'></font></td></tr>"%id
    s += "</table>"

    # We now generate javascript that gets run after the above div
    # gets inserted. This happens because of the setTimeout function
    # below which gets passed an anonymous function.
    s += """<script>(function()
    {
        var values = %s;
        var pos = [%s, %s];
        var sel = '#%s';
        var updatePos = function()
        {
            pos[0]=$(sel).slider('value', 0);
            pos[1]=$(sel).slider('value', 1);
            if(values!=null) $(sel+'-lbl').text("("+values[pos[0]]+", "+values[pos[1]]+")");
        };
        setTimeout(function()
        {
            $(sel).slider(
            {
                range: true,
                stepping: 1,
                min: 0,
                max: %s,
                handles: [{start: %s},{start:%s}],
                change: function(e,ui){ updatePos(); %s; },
                slide: updatePos
            });
            updatePos();
        }, 1);
    })();</script>"""%(values, default_l, default_r, id, steps-1, default_l, default_r, callback)


    # change 'change' to 'slide' and it changes the slider every time it moves;
    # needs much more work to actually work, since server gets fludded by
    # requests.

    return s


def html_color_selector(id, change, input_change, default='000000'):
    """
    Return HTML representation of a jQuery color selector.

    INPUT:
        id -- integer; the id of the html div element that this selector should have
        change -- javascript code to execute when the color selector changes.
        default -- string (default: '000000'); default color as a 6-character
                   HTML hex string.

    OUTPUT:
        string -- HTML that creates the slider.

    EXAMPLES:
        sage: sage.server.notebook.interact.html_color_selector(0, 'alert("changed")', '', default='0afcac')
        '<table>...'
    """
    s = """<table><tr><td><div id='%s-picker'></div></td><td>
<input type='text' id='%s' name='color' onchange='%s;$.farbtastic("#%s-picker").setColor(this.value);' value='%s'/></td></tr></table>"""%(
         id,id,input_change,id,default)
    # We now generate javascript that gets run after the above div
    # gets inserted. This happens because of the setTimeout function
    # below which gets passed an anonymous function.
    # Figuring out the below took understanding jQuery much better,
    # and took me surprisingly long, especially the part involving
    # linkTo which sets the callback.
    s += """<script>setTimeout(function() {
          $('#%s-picker').farbtastic('#%s');
          $.farbtastic('#%s-picker').linkTo(function(color) {
              var t = get_element('%s');
              if(color!=t.value) {
                  t.value = color;
                  t.style.backgroundColor = color;
                  %s;
              }
              return;
            })
       }, 1);</script>"""%(id,id,id,id,change)
    return s


class InteractElement(object):
    def label(self):
        """
        Returns an empty label for this element. This should be
        overridden for subclasses that need a label.

        EXAMPLES:
            sage: from sage.server.notebook.interact import UpdateButton, InteractElement
            sage: b = UpdateButton(1)
            sage: isinstance(b, InteractElement)
            True
            sage: b.label()
            ''
        """
        return ""

    def set_canvas(self, canvas):
        """
        Sets the InteractCanvas on which this element appears.  This
        method is primarily called in the constructor for
        InteractCanvas.

        EXAMPLES:
            sage: from sage.server.notebook.interact import InputBox, InteractCanvas
            sage: B = InputBox('x',2)
            sage: canvas1 = InteractCanvas([B], 3)
            sage: canvas2 = InteractCanvas([B], 3)
            sage: B.canvas() is canvas2
            True
            sage: B.set_canvas(canvas1)
            sage: B.canvas() is canvas1
            True

        """
        self._canvas = canvas

    def canvas(self):
        """
        Returns the InteractCanvas associated to this element.  If no
        canvas has been set (via the set_canvas method), then this
        will return a ValueError.

        EXAMPLES:
            sage: from sage.server.notebook.interact import InputBox, InteractCanvas
            sage: B = InputBox('x',2)
            sage: canvas1 = InteractCanvas([B], 3)
            sage: canvas2 = InteractCanvas([B], 3)
            sage: B.canvas() is canvas2
            True

        """
        if hasattr(self, '_canvas'):
            return self._canvas
        else:
            raise ValueError, "this element does not have a canvas associated with it"


class InteractControl(InteractElement):
    def __init__(self, var, default_value, label=None):
        """
        Abstract base class for interact controls.  These are controls
        that are used in a specific interact.  They have internal
        state information about the specific function being interactd,
        etc.

        INPUT:
             var -- string; name of variable that this control interacts
             default_value -- the default value of the variable
                              corresponding to this control.
             label -- string (default: None) label of this control; if None
                      then defaults to var.

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: InteractControl('x', default_value=5)
            A InteractControl (abstract base class)
        """
        self.__var = var
        self.__cell_id = SAGE_CELL_ID
        self.__default_value = default_value
        self.__adapt_number = new_adapt_number()
        if label is None:
            self.__label = var
        else:
            self.__label = label

        InteractElement.__init__(self)

    def __repr__(self):
        """
        String representation of interact control.

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: InteractControl('x', default_value=5).__repr__()
            'A InteractControl (abstract base class)'
        """
        return "A InteractControl (abstract base class)"

    def value_js(self):
        """
        Javascript that when evaluated gives the current value of this
        control.  This should be redefined in a derived class.

        OUTPUT:
            string -- defaults to NULL -- this should be redefined.

        EXAMPLES:
            sage: sage.server.notebook.interact.InteractControl('x', default_value=5).value_js()
            'NULL'
        """
        return 'NULL'

    def label(self):
        """
        Return the text label of this interact control.

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: InteractControl('x', default_value=5, label='the x value').label()
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
            sage: InteractControl('x', 19/3).default_value()
            19/3
        """
        return self.__default_value

    def html_escaped_default_value(self):
        """
        Returns the HTML escaped default value of the variable
        corresponding to this interact control.  Note that any
        HTML that uses quotes around this should use double
        quotes and not single quotes.

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: InteractControl('x', '"cool"').html_escaped_default_value()
            '&quot;cool&quot;'
            sage: InteractControl('x',"'cool'").html_escaped_default_value()
            "'cool'"
        """
        import cgi
        return cgi.escape(str(self.default_value()), quote=True)

    def adapt_number(self):
        """
        Return integer index into adapt dictionary of function that is
        called to adapt the values of this control to Python.

        OUTPUT:
            an integer

        EXAMPLES:
            sage: from sage.server.notebook.interact import InteractControl
            sage: InteractControl('x', 19/3).adapt_number()       # random -- depends on call order
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
            sage: sage.server.notebook.interact.InteractControl('x', 1)._adaptor('2/3', globals())
            2/3
        """
        return sage_eval(value, globs)

    def interact(self, *args):
        """
        Return a string that when evaluated in Javascript calls the
        javascript interact function with appropriate inputs for
        this control.

        This method will check to see if there is a canvas attached to
        this control and whether or not controls should automatically
        update the output when their values change.  If no canvas is
        associated with this control, then the control will
        automatically update.

        OUTPUT:
            string -- that is meant to be evaluated in Javascript

        EXAMPLES:
            sage: sage.server.notebook.interact.InteractControl('x', 1).interact()
            'interact(..., "sage.server.notebook.interact.update(..., \\"x\\", ..., sage.server.notebook.interact.standard_b64decode(\\""+encode64(NULL)+"\\"), globals());sage.server.notebook.interact.recompute(0)")'
        """
        #We have to do a try/except block here since the
        #control may not have a canvas associated with it.
        try:
            auto_update = self.canvas().is_auto_update()
        except ValueError:
            auto_update = True

        # The following is a crazy line to read because of all the backslashes and try/except.
        # All it does is run the interact function once after setting exactly one
        # dynamic variable.    If setting the dynamic variable fails, due to a KeyError
        python_string = 'sage.server.notebook.interact.update(%s, \\"%s\\", %s, sage.server.notebook.interact.standard_b64decode(\\""+encode64(%s)+"\\"), globals())'%(
            self.cell_id(), self.var(), self.adapt_number(), self.value_js(*args))

        if auto_update:
            python_string += ';sage.server.notebook.interact.recompute(%s)'%self.cell_id()

        s = 'interact(%s, "%s")'%(self.cell_id(), python_string)
        return s

    def var(self):
        """
        Return the name of the variable that this control interacts.

        OUTPUT:
            string -- name of a variable as a string.

        EXAMPLES:
            sage: sage.server.notebook.interact.InteractControl('theta', 1).var()
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
            sage: sage.server.notebook.interact.InteractControl('theta', 1).cell_id()
            0
        """
        return self.__cell_id

class InputBox(InteractControl):
    def __init__(self, var, default_value, label=None, type=None, width = 80):
        """
        An input box interact control.

        InputBox(var, default_value, label, type)

        EXAMPLES:
            sage: sage.server.notebook.interact.InputBox('theta', 1, 'theta')
            An InputBox interactive control with theta=1 and label 'theta'
            sage: sage.server.notebook.interact.InputBox('theta', 1, 'theta', int)
            An InputBox interactive control with theta=1 and label 'theta'
        """
        InteractControl.__init__(self, var, default_value, label)
        self.__type = type
	self.__width = width

    def __repr__(self):
        """
        String representation of an InputBox interactive control.

        EXAMPLES:
            sage: sage.server.notebook.interact.InputBox('theta', 1).__repr__()
            "An InputBox interactive control with theta=1 and label 'theta'"
        """
        return 'An InputBox interactive control with %s=%r and label %r'%(
            self.var(), self.default_value(), self.label())

    def _adaptor(self, value, globs):
        """
        Adapt a user input, which is the text they enter, to be an
        element selected by this control.

        INPUT:
            value -- text entered by user
            globs -- the globals interpreter variables (not used here).

        OUTPUT:
            object

        EXAMPLES:
            sage: sage.server.notebook.interact.InputBox('theta', Color('red'), type=Color)._adaptor('#aaaaaa',globals())
            RGB color (0.6640625, 0.6640625, 0.6640625)
        """
        if self.__type is None:
            return sage_eval(value, globs)
        elif self.__type is str:
            return value
        elif self.__type is Color:
            return Color(value)
        else:
            return self.__type(sage_eval(value,globs))

    def value_js(self):
        """
        Return javascript string that will give the value of this
        control element.

        OUTPUT:
             string -- javascript

        EXAMPLES:
            sage: sage.server.notebook.interact.InputBox('theta', 1).value_js()
            'this.value'
        """
        if self.__type is bool:
            return 'this.checked'
        else:
            return 'this.value'

    def render(self):
        """
        Render this control as a string.

        OUTPUT:
             string -- html format

        EXAMPLES:
            sage: sage.server.notebook.interact.InputBox('theta', 1).render()
            '<input type=\'text\' value="1" size=80 onchange=\'interact(0, "sage.server.notebook.interact.update(0, \\"theta\\", ..., sage.server.notebook.interact.standard_b64decode(\\""+encode64(this.value)+"\\"), globals());sage.server.notebook.interact.recompute(0)")\'></input>'
        """
        if self.__type is bool:
            return """<input type='checkbox' %s width=200px onchange='%s'></input>"""%(
                'checked' if self.default_value() else '',  self.interact())
        elif self.__type is str:
            return """<input type='text' value="%s" size=%s onchange='%s'></input>"""%(
                self.html_escaped_default_value(), self.__width, self.interact())
        else:
            return """<input type='text' value="%s" size=%s onchange='%s'></input>"""%(
                self.html_escaped_default_value(), self.__width,  self.interact())

class ColorInput(InputBox):
    def value_js(self, n):
        """
        Return javascript that evaluates to value of this control.

        INPUT:
            n -- integer, either 0 or 1.

        If n is 0 return code for evaluation by the actual color control.
        If n is 1, return code for the text area that displays the current color.

        EXAMPLES:
            sage: C = sage.server.notebook.interact.ColorInput('c', Color('red'))
            sage: C.value_js(0)
            'color'
            sage: C.value_js(1)
            'this.value'
        """
        if n == 0:
            return 'color'
        else:
            return 'this.value'

    def render(self):
        """
        Render this color input box to html.

        EXAMPLES:
            sage: sage.server.notebook.interact.ColorInput('c', Color('red')).render()
            '<table>...'
        """
        return html_color_selector('color-selector-%s-%s'%(self.var(), self.cell_id()),
                     change=self.interact(0), input_change=self.interact(1),
                     default=self.default_value().html_color())



class InputGrid(InteractControl):
    def __init__(self, var, rows, columns, default_value=None, label=None, to_value=lambda x: x, width=4):
        """
        A grid interact control.

        INPUT
            var -- the variable
            rows -- the number of rows
            columns -- the number of columns
            default_value -- if this is a scalar, it is put in every
                cell; if it is a list, it is filled into the cells row by
                row; if it is a nested list, then it is filled into the
                cells according to the nesting structure.
            label -- the label for the control
            to_value -- a function which is applied to the nested list
                from user input when assigning the variable
            width -- the width of the input boxes

        EXAMPLES:
            sage: sage.server.notebook.interact.InputGrid('M', 2,2, default_value = 0, label='M')
            A 2 x 2 InputGrid interactive control with M=[[0, 0], [0, 0]] and label 'M'
            sage: sage.server.notebook.interact.InputGrid('M', 2,2, default_value = [[1,2],[3,4]], label='M')
            A 2 x 2 InputGrid interactive control with M=[[1, 2], [3, 4]] and label 'M'
            sage: sage.server.notebook.interact.InputGrid('M', 2,2, default_value = [[1,2],[3,4]], label='M', to_value=MatrixSpace(ZZ,2,2))
            A 2 x 2 InputGrid interactive control with M=[1 2]
            [3 4] and label 'M'
            sage: sage.server.notebook.interact.InputGrid('M', 1, 3, default_value=[1,2,3], to_value=lambda x: vector(flatten(x)))
            A 1 x 3 InputGrid interactive control with M=(1, 2, 3) and label 'M'
        """

        self.__rows = rows
        self.__columns = columns
        self.__to_value = to_value
        self.__width = width

        if type(default_value) != list:
            default_value = [[default_value for _ in range(columns)] for _ in range(rows)]
        elif not all(type(elt)==list for elt in default_value):
            default_value = [[default_value[i*columns+j] for j in xrange(columns)] for i in xrange(rows)]

        self.__default_value_grid = default_value

        InteractControl.__init__(self, var, self.__to_value(default_value), label)

    def __repr__(self):
        """
        String representation of an InputGrid interactive control.

        EXAMPLES:
            sage: sage.server.notebook.interact.InputGrid('M', 2,2).__repr__()
            "A 2 x 2 InputGrid interactive control with M=[[None, None], [None, None]] and label 'M'"
        """

        return 'A %r x %r InputGrid interactive control with %s=%r and label %r'%( self.__rows,
                                self.__columns, self.var(),  self.default_value(), self.label())


    def _adaptor(self, value, globs):
        """
        Adapt a user input, which is the text they enter, to be an
        element selected by this control.

        INPUT:
            value -- text entered by user
            globs -- the globals interpreter variables (not used here).

        OUTPUT:
            object

        EXAMPLES:
            sage: sage.server.notebook.interact.InputGrid('M', 1,3, default_value=[[1,2,3]], to_value=lambda x: vector(flatten(x)))._adaptor("[[4,5,6]]", globals())
            (4, 5, 6)
        """

        return self.__to_value(sage_eval(value, globs))

    def value_js(self):
        """
        Return javascript string that will give the value of this
        control element.

        OUTPUT:
             string -- javascript

        EXAMPLES:
            sage: sage.server.notebook.interact.InputGrid('M', 2,2).value_js()
            ' "[["+jQuery(this).parents("table").eq(0).find("tr").map(function(){return jQuery(this).find("input").map(function() {return jQuery(this).val();}).get().join(",");}).get().join("],[")+"]]" '
        """
        # Basically, given an input element in a table, it constructs
        # a python string representation of a list of lists from the
        # rows in the table.

        return """ "[["+jQuery(this).parents("table").eq(0).find("tr").map(function(){return jQuery(this).find("input").map(function() {return jQuery(this).val();}).get().join(",");}).get().join("],[")+"]]" """

    def render(self):
        """
        Render this control as a string.

        OUTPUT:
             string -- html format

        EXAMPLES:
            sage: sage.server.notebook.interact.InputGrid('M', 1,2).render()
            '<table><tr><td><input type=\'text\' value=\'None\' ...

        """
        table = "<table>"
        for i in range(self.__rows):
            table += "<tr>"
            for j in range(self.__columns):
                table += "<td><input type='text' value='%r' size='%s' onchange='%s'></input></td>"%(self.__default_value_grid[i][j], self.__width, self.interact())
            table += "</tr>"
        table += "</table>"

        return table



class Selector(InteractControl):
    def __init__(self, var, values, label=None, default=0,
                 nrows=None, ncols=None, width=None, buttons=False):
        """
        A drop down menu or a button bar that when pressed sets a
        variable to a given value.

        Selector(var, values, label=None, nrows=None, ncols=None)

        INPUT:
            var   -- string; variable name
            values-- list; button values
            label -- string (default: None) label off to the left for this button group
            default -- integer (default: 0) position of default value in values list.
            nrows -- integer (default: None) number of rows
            ncols -- integer (default: None) number of columns
            width -- integer (default: None) width of all the buttons
            buttons -- bool (default: False) if True use buttons instead of dropdown

        EXAMPLES:
            sage: sage.server.notebook.interact.Selector('x', [1..5], 'alpha', default=2)
            Selector with 5 options for variable 'x'
            sage: sage.server.notebook.interact.Selector('x', [1..4], 'alpha', default=2, nrows=2, ncols=2, width=10, buttons=True)
            Selector with 4 options for variable 'x'
        """
        if len(values) > 0 and isinstance(values[0], tuple) and len(values[0]) == 2:
            vals = [z[0] for z in values]
            lbls = [str(z[1]) if z[1] is not None else None for z in values]
        else:
            vals = values
            lbls = [None]*len(vals)

        default = int(default)
        if default < 0 or default >= len(vals):
            default = 0

        InteractControl.__init__(self, var, vals[default], label)

        self.__default = default
        self.__buttons = buttons
        self.__values = vals
        self.__labels = lbls
        if nrows is None:
            if ncols is not None:
                nrows = len(values)/ncols
                if ncols * nrows < len(values):
                    nrows += 1
            else:
                nrows = 1 # temporary
        else:
            nrows = int(nrows)
            if nrows <= 0:
                nrows = 1
        if ncols is None:
            ncols = len(values)/nrows
            if ncols * nrows < len(values):
                ncols += 1

        self.__nrows = nrows
        self.__ncols = ncols

        if width is not None:
            self.__width = "width:%sex;"%width
        else:
            self.__width = ''

        self.__selected = 'background-color:orange;'

    def __repr__(self):
        """
        String representation of a Selector interactive control.

        EXAMPLES:
            sage: sage.server.notebook.interact.Selector('x', [1..5]).__repr__()
            "Selector with 5 options for variable 'x'"
        """
        return "Selector with %s options for variable '%s'"%(len(self.__values), self.var())

    def _adaptor(self, value, globs):
        """
        Adapt value of button or menu selection.

        The button value is just an integer, and this function adapts
        it to be the value that we associate with that button.

        INPUT:
            value -- value sent in via javascript
            globs -- the globals interpreter variables (not used here).

        OUTPUT:
            object

        EXAMPLES:
            sage: S = sage.server.notebook.interact.Selector('x', ['first',x^3+5])
            sage: S._adaptor(0,globals())
            'first'
            sage: S._adaptor(1,globals())
            x^3 + 5
        """
        return self.__values[int(value)]

    def use_buttons(self):
        """
        Whether or not to use buttons instead of a drop
        down menu for this select list.

        OUTPUT:
            bool

        EXAMPLES:
            sage: sage.server.notebook.interact.Selector('x', [1..5]).use_buttons()
            False
            sage: sage.server.notebook.interact.Selector('x', [1..5], buttons=True).use_buttons()
            True
        """
        return self.__buttons

    def value_js(self):
        """
        Return javascript string that will give the value of this
        control element.

        OUTPUT:
             string -- javascript

        EXAMPLES:
            sage: sage.server.notebook.interact.Selector('x', [1..5]).value_js()
            'this.options[this.selectedIndex].value'
            sage: sage.server.notebook.interact.Selector('x', [1..5], buttons=True).value_js()
            'this.value'
        """
        if self.use_buttons():
            return 'this.value'
        else:
            # Now we have to use a option selector.
            return 'this.options[this.selectedIndex].value'

    def render(self):
        """
        Render this control as a string.

        OUTPUT:
             string -- html format

        EXAMPLES:
            sage: sage.server.notebook.interact.Selector('x', [1..5]).render()
            '<select...</select>'
            sage: sage.server.notebook.interact.Selector('x', [1..5], buttons=True).render()
            '<table...</table>'
        """
        width = self.__width
        vals = self.__values
        lbls = self.__labels
        default = self.__default
        label = self.label()
        use_buttons = self.use_buttons()
        event = self.interact()
        if use_buttons:
        	#On selected buttons, border is set to inset, on unselected boxes - outset. This usually is default rendering.
        	if len(vals) > 1:
		       	event = '$("BUTTON", this.parentNode).css("border-style", "outset"); $(this).css("border-style", "inset"); %s'%event
        	s = '<table style="border:1px solid #dfdfdf;background-color:#efefef">'
        else:
            s = "<select onchange='%s;'>"%event
        i = 0
        for r in range(self.__nrows):
            if use_buttons:
                s += '\n<tr><td>'
            for c in range(self.__ncols):
                if i >= len(vals):
                    i += 1
                    continue
                style = width
                #if i == default:
                #    style += self.__selected
                if lbls[i] is None:
                    if isinstance(vals[i], str):
                        lbl = vals[i]
                    else:
                        lbl = repr(vals[i])
                else:
                    lbl = lbls[i]
                if use_buttons:
                    s += "<button style='%s%s' value='%s' onclick='%s'>%s</button>\n"%('border-style:inset;' if i==default and len(vals)>1 else 'border-style:outset;', style, i, event, lbl)
                else:
                    s += "<option value='%s' %s>%s</option>\n"%(i, 'selected' if i==default else '', lbl)
                i += 1
            if use_buttons:
                s += '</td></tr>'
        if use_buttons:
            s += '</table>'
        else:
            s += '</select>'
        return s


class SliderGeneric(InteractControl):
    def __init__(self, var, values, default_value, label=None, display_value=True):
        """
        An abstract slider interact control that takes on the given list of
        values.

        INPUT:
            var -- string; name of variable being interactd
            values -- list; a list of the values that the slider will take on
            default_value -- default valueoif slider.
            label -- alternative label to the left of the slider,
                     instead of the variable.
            display_value -- boolean, whether to display the current value
                             on the slider

        EXAMPLES:
            sage: sage.server.notebook.interact.SliderGeneric('x', [1..5], 2, 'alpha')
            Abstract Slider Interact Control: alpha [1--|2|---5]
        """
        InteractControl.__init__(self, var, default_value, label=label)
        self.__values = values
        self.__display_value = display_value

    def __repr__(self):
        """
        Return string representation of this slider control.

        EXAMPLES:
            sage: sage.server.notebook.interact.SliderGeneric('x', [1..5], 2, 'alpha').__repr__()
            'Abstract Slider Interact Control: alpha [1--|2|---5]'
        """
        return "Abstract Slider Interact Control: %s [%s--|%s|---%s]"%(
            self.label(), self.__values[0],
            self.default_value(), self.__values[-1])

    def values(self):
        """
        Return list of values the slider acts on.

        OUTPUT:
            list

        EXAMPLES:
            sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha').values()
            [1, 2, 3, 4, 5]
        """
        return self.__values

    def display_value(self):
        """
        Returns whether to display the value on the slider.

        OUTPUT:
            boolean

        EXAMPLES:
            sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha').display_value()
            True
        """
        return self.__display_value

    def values_js(self):
        """
        Returns Javascript array representation of values or null if display_value is False

        OUTPUT:
            string

        EXAMPLES:
            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha').values_js()
            '["1","2","3","4","5"]'
            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha', False).values_js()
            'null'
        """
        if self.__display_value == False:
            return "null"
        s = "["
        for i in self.__values:
            ie = str(i).replace("\\","\\\\").replace("\"","\\\"").replace("'","\\'")
            s += "\"%s\","%ie
        s = s[:-1] + ']'
        return s


class Slider(SliderGeneric):
    def __init__(self, var, values, default_position, label=None, display_value=True):
        """
        A slider interact control that takes on the given list of
        values.

        INPUT:
            var -- string; name of variable being interactd
            values -- list; a list of the values that the slider will take on
            default_position -- int; default location that the slider is set to.
            label -- alternative label to the left of the slider,
                     instead of the variable.
            display_value -- boolean, whether to display the current value right
                             of the slider

        EXAMPLES:
            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha')
            Slider Interact Control: alpha [1--|3|---5]
        """
        SliderGeneric.__init__(self, var, values, values[default_position], label=label, display_value=display_value)
        self.__default_position = default_position

    def __repr__(self):
        """
        Return string representation of this slider control.

        EXAMPLES:
            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha').__repr__()
            'Slider Interact Control: alpha [1--|3|---5]'
        """
        return "Slider Interact Control: %s [%s--|%s|---%s]"%(
            self.label(), self.values()[0],
            self.default_value(), self.values()[-1])

    def default_position(self):
        """
        Return the default position (as an integer) of the slider.

        EXAMPLES:
            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha').default_position()
            2
        """
        return self.__default_position

    def value_js(self):
        """
        Return javascript string that will give the
        value of this control element.

        OUTPUT:
             string -- javascript

        EXAMPLES:
            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha').value_js()
            'position'
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

        EXAMPLES:
            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha')._adaptor(2,globals())
            3
        """
        v = self.values()
        # We have to cast to int, since it comes back as a float that
        # is too big.
        return v[int(position)]

    def render(self):
        """
        Render this control as an HTML string.

        OUTPUT:
             string -- html format

        EXAMPLES:
            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha').render()
            '<table>...<div ...var values = ["1","2","3","4","5"];...'

            sage: sage.server.notebook.interact.Slider('x', [1..5], 2, 'alpha', display_value=False).render()
            '<table>...<div ...var values = null;...'
        """

        return html_slider('slider-%s-%s'%(self.var(), self.cell_id()),
                           self.values_js(), self.interact(), steps=len(self.values()),
                           default=self.default_position())


class RangeSlider(SliderGeneric):
    def __init__(self, var, values, default_position, label=None, display_value=True):
        """
        A range slider interact control that takes on the given list of
        values.

        INPUT:
            var -- string; name of variable being interactd
            values -- list; a list of the values that the slider will take on
            default_position -- (int,int); default location that the slider is set to.
            label -- alternative label to the left of the slider,
                     instead of the variable.
            display_value -- boolean, whether to display the current value below
                             the slider

        EXAMPLES:
            sage: sage.server.notebook.interact.RangeSlider('x', [1..5], (2,3), 'alpha')
            Range Slider Interact Control: alpha [1--|3==4|---5]
        """
        SliderGeneric.__init__(self, var, values, (values[default_position[0]], values[default_position[1]]), label=label, display_value=display_value)
        self.__default_position = default_position

    def __repr__(self):
        """
        Return string representation of this slider control.

        EXAMPLES:
            sage: sage.server.notebook.interact.RangeSlider('x', [1..5], (2,3), 'alpha').__repr__()
            'Range Slider Interact Control: alpha [1--|3==4|---5]'
        """
        return "Range Slider Interact Control: %s [%s--|%s==%s|---%s]"%(
            self.label(), self.values()[0],
            self.default_value()[0], self.default_value()[1], self.values()[-1])

    def default_position(self):
        """
        Return the default position (as an integer) of the slider.

        EXAMPLES:
            sage: sage.server.notebook.interact.RangeSlider('x', [1..5], (2,3), 'alpha').default_position()
            (2, 3)
        """
        return self.__default_position

    def value_js(self):
        """
        Return javascript string that will give the
        value of this control element.

        OUTPUT:
             string -- javascript

        EXAMPLES:
            sage: sage.server.notebook.interact.RangeSlider('x', [1..5], (2,3), 'alpha').value_js()
            "pos[0]+' '+pos[1]"
        """
        return "pos[0]+' '+pos[1]"

    def _adaptor(self, position, globs):
        """
        Adapt a user input, which is the slider position, to be an
        element selected by this control.

        INPUT:
            position -- position of the slider
            globs -- the globals interpreter variables (not used here).

        OUTPUT:
            object

        EXAMPLES:
            sage: sage.server.notebook.interact.RangeSlider('x', [1..5], (2,3), 'alpha')._adaptor("2 3",globals())
            (3, 4)
        """
        v = self.values()
        s = position.split(' ')
        # use of int() here matches it's use in Slider._adaptor
        return (v[int(s[0])], v[int(s[1])])

    def render(self):
        """
        Render this control as an HTML string.

        OUTPUT:
             string -- html format

        EXAMPLES:
            sage: sage.server.notebook.interact.RangeSlider('x', [1..5], (2,3), 'alpha').render()
            '<table>...<div ...var values = ["1","2","3","4","5"];...'

            sage: sage.server.notebook.interact.RangeSlider('x', [1..5], (2,3), 'alpha', display_value=False).render()
            '<table>...<div ...var values = null;...'
        """

        return html_rangeslider('slider-%s-%s'%(self.var(), self.cell_id()),
                           self.values_js(), self.interact(), steps=len(self.values()),
                           default_l=self.default_position()[0], default_r=self.default_position()[1])



class TextControl(InteractControl):
    def __init__(self, var, data):
        """
        A text field interact control

        INPUT:
            data -- the HTML value of the text field

        EXAMPLES:
            sage: sage.server.notebook.interact.TextControl('x', 'something')
            Text Interact Control: something
        """
        InteractControl.__init__(self, var, data, label='')
        self.__data = data

    def __repr__(self):
        """
        Return string representation of this control.

        EXAMPLES:
            sage: sage.server.notebook.interact.TextControl('x', 'something').__repr__()
            'Text Interact Control: something'
        """
        return 'Text Interact Control: %s'%self.default_value()

    def render(self):
        """
        Render this control as an HTML string.

        OUTPUT:
             string -- html format

        EXAMPLES:
            sage: sage.server.notebook.interact.TextControl('x', 'something').render()
            '<div ...>something</div>'
        """
        return '<div style="color:black; padding-bottom:5px">%s</div>'%self.default_value()


class InteractCanvas:
    def __init__(self, controls, id, **options):
        """
        Base class for interact canvases. This is where all the controls
        along with the output of the interactd function are layed out
        and rendered.

        INPUT:
            controls -- a list of InteractControl instances.
            id -- the id of the cell that contains this InteractCanvas.
            options -- any additional keyword arguments (for example, auto_update=False)

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: sage.server.notebook.interact.InteractCanvas([B], 3)
            Interactive canvas in cell 3 with 1 controls
        """
        for control in controls:
            control.set_canvas(self)

        self.__controls = controls
        self.__cell_id = id
        self.__options = options

    def __repr__(self):
        """
        Print representation of an interactive canvas.

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: sage.server.notebook.interact.InteractCanvas([B], 3).__repr__()
            'Interactive canvas in cell 3 with 1 controls'
        """
        return "Interactive canvas in cell %s with %s controls"%(
            self.__cell_id, len(self.__controls))

    def is_auto_update(self):
        """
        Returns True if any change of the values for the controls on
        this canvas should cause an update.  If auto_update=False was
        not specified in the constructor for this canvas, then this
        will default to True.

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: canvas = sage.server.notebook.interact.InteractCanvas([B], 3)
            sage: canvas.is_auto_update()
            True
            sage: canvas = sage.server.notebook.interact.InteractCanvas([B], 3, auto_update=False)
            sage: canvas.is_auto_update()
            False
        """
        return self.__options.get('auto_update', True)

    def cell_id(self):
        """
        Returns the cell id associated to this InteractCanvas.

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: canvas = sage.server.notebook.interact.InteractCanvas([B], 3)
            sage: canvas.cell_id()
            3
        """
        return self.__cell_id

    def controls(self):
        """
        Return list of controls in this canvas.

        WARNING: Returns a reference to a mutable list.

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: sage.server.notebook.interact.InteractCanvas([B], 3).controls()
            [An InputBox interactive control with x=2 and label 'x']
        """
        return self.__controls

    def cell_id(self):
        """
        Return the cell id that contains this interactive canvas.

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: sage.server.notebook.interact.InteractCanvas([B], 3).cell_id()
            3
        """
        return self.__cell_id

    def render_output(self):
        """
        Render in text (html) form the output portion of the interact canvas.

        The output contains two special tags, <?TEXT> and <?HTML>,
        which get replaced at runtime by the text and html parts
        of the output of running the function.

        OUTPUT:
            string -- html

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: sage.server.notebook.interact.InteractCanvas([B], 3).render_output()
            "<div ...</div>"
        """
        return """<div id='cell-interact-%s'><?__SAGE__START>
        <table border=0 bgcolor='#white' width=100%% height=100%%>
        <tr><td bgcolor=white align=left valign=top><pre><?__SAGE__TEXT></pre></td></tr>
        <tr><td  align=left valign=top><?__SAGE__HTML></td></tr>
        </table><?__SAGE__END></div>"""%self.cell_id()

    def render_controls(self):
        """
        Render in text (html) form all the input controls.

        OUTPUT:
            string -- html

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: sage.server.notebook.interact.InteractCanvas([B], 3).render_controls()
            '<table>...'
        """
        tbl_body = ''
        for c in self.__controls:
            if c.label() == '':
                tbl_body += '<tr><td colspan=2>%s</td></tr>\n'%c.render()
            else:
                tbl_body += '<tr><td align=right><font color="black">%s&nbsp;</font></td><td>%s</td></tr>\n'%(
                c.label(), c.render())
        return '<table>%s</table>'%tbl_body

    def wrap_in_outside_frame(self, inside):
        """
        Return the entire HTML for the interactive canvas, obtained by
        wrapping all the inside html of the canvas in a div and a
        table.

        INPUT:
            inside -- string (of HTML)

        OUTPUT:
            string of HTML

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: sage.server.notebook.interact.InteractCanvas([B], 3).wrap_in_outside_frame('<!--inside-->')
            "<!--notruncate--><div padding=6 id='div-interact-3'> ...</div>\n                 "
        """
        return """<!--notruncate--><div padding=6 id='div-interact-%s'> <table width=800px height=20px bgcolor='#c5c5c5'
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

        EXAMPLES:
            sage: B = sage.server.notebook.interact.InputBox('x',2)
            sage: sage.server.notebook.interact.InteractCanvas([B], 3).render()
            '<!--notruncate--><div padding=6 id=\'div-interact-3\'> ...</div>\n                 '
        """
        s = "%s%s"%(self.render_controls(), self.render_output())
        s = self.wrap_in_outside_frame(s)
        return s

class JavascriptCodeButton(InteractElement):
    def __init__(self, label, code):
        """
        This interact element displays a button which when clicked
        executes Javascript code in the notebook.

        EXAMPLES:
            sage: b = sage.server.notebook.interact.JavascriptCodeButton('Push me', 'alert("2")')
        """
        self.__label = label
        self.__code = code
        InteractElement.__init__(self)

    def render(self):
        r"""
        Returns the HTML to display this button.

        EXAMPLES:
            sage: b = sage.server.notebook.interact.JavascriptCodeButton('Push me', 'alert("2")')
            sage: b.render()
            '<input type="button" value="Push me" onclick=\'alert("2")\'>\n'

        """
        return '<input type="button" value="%s" onclick=\'%s\'>\n'%(self.__label, self.__code)

class UpdateButton(JavascriptCodeButton):
    def __init__(self, cell_id):
        r"""
        This interact element creates a button which when clicked
        causes the interact function in the cell cell_id to be
        recomputed with the current values of the variables.

        EXAMPLES:
            sage: b = sage.server.notebook.interact.UpdateButton(0)
            sage: b.render()
            '<input type="button" value="Update" onclick=\'interact(0, "sage.server.notebook.interact.recompute(0)")\'>\n'

        """
        s = 'interact(%s, "sage.server.notebook.interact.recompute(%s)")'%(cell_id, cell_id)
        JavascriptCodeButton.__init__(self, "Update", s)

def interact(f):
    r"""
    Use interact as a decorator to create interactive Sage notebook
    cells with sliders, text boxes, radio buttons, check boxes, and
    color selectors.  Simply put @interact on the line before a
    function definition in a cell by itself, and choose appropriate
    defaults for the variable names to determine the types of
    controls (see tables below).

    INPUT:
        f -- a Python function

    EXAMPLES:
    In each example below we use a single underscore for the function
    name.  You can use \emph{any} name you want; it does not have to
    be an underscore.

    We create an interact control with two inputs, a text input for
    the variable $a$ and a $y$ slider that runs through the range of
    integers from $0$ to $19$.
        sage: @interact
        ... def _(a=5, y=(0..20)): print a + y
        ...
        <html>...

    Draw a plot interacting with the ``continuous'' variable $a$.  By
    default continuous variables have exactly 50 possibilities.
        sage: @interact
        ... def _(a=(0,2)):
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

    Place a block of text among the controls:
        sage: @interact
        ... def _(t1=text_control("Factors an integer."), n="1"):
        ...     print factor(Integer(n))
        <html>...

    You do not have to use interact as a decorators; you can also
    simply write \code{interact(f)} where f is any Python function
    that you have defined, though this is frowned on.  E.g., f can
    also be a library function as long as it is written in Python:

        sage: interact(matrix)   # put ZZ, 2,2,[1..4] in boxes...
        <html>...

    If your the time to evaluate your function takes awhile, you may
    not want to have it reevaluated every time the inputs change.  In
    order to prevent this, you can add a keyword
    \code{auto_update=False} to your function to prevent it from
    updating whenever the values are changed.  This will cause a
    button labeled 'Update' to appear which you can click on to
    revaluate your function.

        sage: @interact
        ... def _(n=(10,100,1), auto_update=False):
        ...     show(factor(x^n - 1))
        <html>...

    DEFAULTS:
    Defaults for the variables of the input function determine
    interactive controls.  The standard controls are \code{input_box},
    \code{slider}, \code{range_slider}, \code{checkbox}, \code{selector},
    \code{input_grid}.  There is also a color selector and text control
    (see defaults below).

    \begin{itemize}
        \item u = input_box(default=None, label=None, type=None)
                         -- input box with given default; use type=str to
                            get input as an arbitrary string
        \item u = slider(vmin, vmax=None,step_size=1,default=None,label=None)
                         -- slider with given list of possible values; vmin can be a list
        \item u = range_slider(vmin, vmax=None,step_size=1,default=None,label=None)
                         -- range slider with given list of possible values;
                            vmin can be a list
        \item u = checkbox(default=True, label=None)
                         -- a checkbox
        \item u = selector(values, label=None, nrows=None, ncols=None, buttons=False)
                         -- a dropdown menu or buttons (get buttons if nrows,
                            ncols, or buttons is set, otherwise a dropdown menu)
        \item u = input_grid(nrows, ncols, default=None, label=None,
                             to_value=lambda x:x, width=4)
                         -- an editable grid of objects (a matrix or array)
        \item u = text_control(value='')
                         -- a block of text
    \end{itemize}

    You can create a color selector by setting the default value for a
    variable to Color(...).

    There are also some convenient defaults that allow you to make
    controls automatically without having to explicitly specify them.
    E.g., you can make $x$ a continuous slider of values between $u$
    and $v$ by just writing \code{x=(u,v)} in the argument list of
    your function.  These are all just convenient shortcuts for
    creating the controls listed above.

    \begin{itemize}
        \item u                 -- blank input_box field
        \item u = element       -- input_box with default=element, if element not below.
        \item u = (umin,umax)   -- continuous slider (really 100 steps)
        \item u = (umin,umax,du)-- slider with step size du
        \item u = list          -- buttons if len(list) at most 5; otherwise, drop down
        \item u = generator     -- a slider (up to 10000 steps)
        \item u = bool          -- a checkbox
        \item u = Color('blue') -- a 2d RGB color selector; returns Color object
        \item u = (default, v)  -- v as above, with given default value
        \item u = (label, v)    -- v as above, with given label (a string)
        \item u = matrix        -- an input_grid with to_value set to matrix.parent()
                                   and default values given by the matrix
    \end{itemize}

    WARNING: Suppose you would like to make a interactive with a
    default rgb color of (1,0,0), so the function would have signature
    \code{f(color=(1,0,0))}.  Unfortunately, the above shortcuts reinterpret
    the (1,0,0) as a discrete slider with step size 0 between 1 and 0.
    Instead you should do the following:
        sage: @interact
        ... def _(v = input_box((1,0,0))):
        ...       show(plot(sin,color=v))
        <html>...

    MORE EXAMPLES:
    We give an input box that allows one to enter completely arbitrary strings.
        sage: @interact
        ... def _(a=input_box('sage', label="Enter your name", type=str)):
        ...        print "Hello there %s"%a.capitalize()
        <html>...

    The scope of variables that you control via interact are local to
    the scope of the function being interacted with. However, by using
    the global Python keyword, you can still modify global variables
    as follows:
        sage: xyz = 10
        sage: @interact
        ... def _(a=('xyz',5)):
        ...       global xyz
        ...       xyz = a
        <html>...

    If you enter the above you obtain an interact canvas.  Entering
    values in the box, changes the global variable xyz.

        sage: @interact
        ... def _(title=["A Plot Demo", "Something silly", "something tricky"], a=input_box(sin(x*sin(x*sin(x))), 'function'),
        ...     clr = Color('red'), thickness=[1..30], zoom=(1,0.95,..,0.1), plot_points=(200..2000)):
        ...     html('<h1 align=center>%s</h1>'%title)
        ...     print plot_points
        ...     show(plot(a, -zoom*pi,zoom*pi, color=clr, thickness=thickness, plot_points=plot_points))
        <html>...

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

    An example where the range slider control is useful.
        sage: @interact
        ... def _(b = range_slider(-20, 20, 1, default=(-19,3), label='Range')):
        ...     plot(sin(x)/x, b[0], b[1]).show(xmin=b[0],xmax=b[1])
        <html>...

    An example using checkboxes, obtained by making the default values bools.
        sage: @interact
        ... def _(axes=('Show axes', True), square=False):
        ...       show(plot(sin, -5,5), axes=axes, aspect_ratio = (1 if square else None))
        <html>...

    An example generating a random walk that uses a checkbox control to determine
    whether points are placed at each step:
        sage: @interact
        ... def foo(pts = checkbox(True, "points"), n = (50,(10..100))):
        ...       s = 0; v = [(0,0)]
        ...       for i in range(n):
        ...            s += random() - 0.5
        ...            v.append((i, s))
        ...       L = line(v, rgbcolor='#4a8de2')
        ...       if pts: L += points(v, pointsize=20, rgbcolor='black')
        ...       show(L)
        <html>...

    You can rotate and zoom into 3D graphics while
    interacting with a variable.
        sage: @interact
        ... def _(a=(0,1)):
        ...     x,y = var('x,y')
        ...     show(plot3d(sin(x*cos(y*a)), (x,0,5), (y,0,5)), figsize=4)
        <html>...

    A random polygon:
        sage: pts = [(random(), random()) for _ in xrange(20)]
        sage: @interact
        ... def _(n = (4..len(pts)), c=Color('purple') ):
        ...     G = points(pts[:n],pointsize=60) + polygon(pts[:n], rgbcolor=c)
        ...     show(G, figsize=5, xmin=0, ymin=0)
        <html>...

    Two "sinks" displayed simultaneously via a contour plot and a 3d
    interactive plot:
        sage: @interact
        ... def _(q1=(-1,(-3,3)), q2=(-2,(-3,3))):
        ...     x,y = var('x,y')
        ...     f = q1/sqrt((x+1)^2 + y^2) + q2/sqrt((x-1)^2+(y+0.5)^2)
        ...     C = contour_plot(f, (-2,2), (-2,2), plot_points=30, contours=15, cmap='cool')
        ...     show(C, figsize=3, aspect_ratio=1)
        ...     show(plot3d(f, (x,-2,2), (y,-2,2)), figsize=4)
        <html>...

    This is similar to above, but you can select the color map from a dropdown menu:
        sage: @interact
        ... def _(q1=(-1,(-3,3)), q2=(-2,(-3,3)),
        ...    cmap=['autumn', 'bone', 'cool', 'copper', 'gray', 'hot', 'hsv',
        ...          'jet', 'pink', 'prism', 'spring', 'summer', 'winter']):
        ...     x,y = var('x,y')
        ...     f = q1/sqrt((x+1)^2 + y^2) + q2/sqrt((x-1)^2+(y+0.5)^2)
        ...     C = contour_plot(f, (x,-2,2), (y,-2,2), plot_points=30, contours=15, cmap=cmap)
        ...     show(C, figsize=3, aspect_ratio=1)
        <html>...

    A quadratic roots etch-a-sketch:
        sage: v = []
        sage: html('<h2>Quadratic Root Etch-a-sketch</h2>')
        <html><font color='black'><h2>Quadratic Root Etch-a-sketch</h2></font></html>
        sage: @interact
        ... def _(a=[-10..10], b=[-10..10], c=[-10..10]):
        ...       f = a*x^2 + b*x + c == 0; show(f)
        ...       soln = solve(a*x^2 + b*x + c == 0, x)[0].rhs()
        ...       show(soln)
        ...       P = tuple(CDF(soln))
        ...       v.append(P)
        ...       show(line(v, rgbcolor='purple') + point(P, pointsize=200))
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

    An input grid:
        sage: @interact
        ... def _(A=matrix(QQ,3,3,range(9)), v=matrix(QQ,3,1,range(3))):
        ...     try:
        ...         x = A\v
        ...         html('$$%s %s = %s$$'%(latex(A), latex(x), latex(v)))
        ...     except:
        ...         html('There is no solution to $$%s x=%s$$'%(latex(A), latex(v)))
        <html>...

    """

    (args, varargs, varkw, defaults) = inspect.getargspec(f)

    if defaults is None:
        defaults = []

    n = len(args) - len(defaults)
    controls = [automatic_control(defaults[i-n] if i >= n else None).render(arg)
                for i, arg in enumerate(args)]

    variables = {}
    adapt = {}
    state[SAGE_CELL_ID] = {'variables':variables, 'adapt':adapt}

    for control in controls:
        variables[control.var()] = control.default_value()
        adapt[control.adapt_number()] = control._adaptor

    #Replace the auto_update checkbox with a button that
    #will cause the cell to recompute itself.
    auto_update = variables.get('auto_update', True)
    if auto_update is False:
        i = args.index('auto_update')
        controls[i] = UpdateButton(SAGE_CELL_ID)

    C = InteractCanvas(controls, SAGE_CELL_ID, auto_update=auto_update)
    html(C.render())

    def _():
        z = f(*[variables[arg] for arg in args])
        if z: print z

    state[SAGE_CELL_ID]['function'] = _

    return f

######################################################
# Actual control objects that the user passes in
######################################################
class control:
    def __init__(self, label=None):
        """
        An interactive control object used with the interact command.
        This is the abstract base class.

        INPUTS:
            label -- a string

        EXAMPLES:
            sage: sage.server.notebook.interact.control('a control')
            Interative control 'a control' (abstract base class)
        """
        self.__label = label

    def __repr__(self):
        """
        Return string representation of this control.
        (It just mentions the label and that this is an abstract base class.)

        EXAMPLES:
            sage: sage.server.notebook.interact.control('a control').__repr__()
            "Interative control 'a control' (abstract base class)"
        """
        return "Interative control '%s' (abstract base class)"%self.__label

    def label(self):
        """
        Return the label of this control.

        OUTPUT:
            a string

        EXAMPLES:
            sage: sage.server.notebook.interact.control('a control').label()
            'a control'
            sage: selector([1,2,7], 'alpha').label()
            'alpha'
        """
        return self.__label

    def set_label(self, label):
        """
        Set the label of this control.

        INPUT:
            label -- a string

        EXAMPLES:
            sage: C = sage.server.notebook.interact.control('a control')
            sage: C.set_label('sage'); C
            Interative control 'sage' (abstract base class)
        """
        self.__label = label

class input_box(control):
    def __init__(self, default=None, label=None, type=None, width = 80):
        r"""
        An input box interactive control.  Use this in conjunction
        with the interact command.

        \code{input_box(default=None, label=None, type=None)}

        INPUT:
            default -- object; the default put in this input box
            label -- the label rendered to the left of the box.
            type -- coerce inputs to this; this doesn't have to be
                    an actual type, since anything callable will do.
            width -- width of text box in characters

        EXAMPLES:
            sage: input_box("2+2", 'expression')
            Interact input box labeled 'expression' with default value '2+2'
            sage: input_box('sage', label="Enter your name", type=str)
            Interact input box labeled 'Enter your name' with default value 'sage'
        """
        self.__default = default
        self.__type = type
        control.__init__(self, label)
	self.__width = width

    def __repr__(self):
        """
        Return print representation of this input box.

        EXAMPLES:
            sage: input_box("2+2", 'expression').__repr__()
            "Interact input box labeled 'expression' with default value '2+2'"
        """
        return "Interact input box labeled %r with default value %r"%(self.label(), self.__default)

    def default(self):
        """
        Return the default value of this input box.

        EXAMPLES:
            sage: input_box('2+2', 'Expression').default()
            '2+2'
            sage: input_box(x^2 + 1, 'Expression').default()
            x^2 + 1
            sage: checkbox(True, "Points").default()
            True
        """
        return self.__default

    def type(self):
        """
        Return the type that elements of this input box are coerced to
        or None if they are not coerced (they have whatever type they
        evaluate to).

        EXAMPLES:
            sage: input_box("2+2", 'expression', type=int).type()
            <type 'int'>
            sage: input_box("2+2", 'expression').type() is None
            True
        """
        return self.__type

    def render(self, var):
        r"""
        Return rendering of this input box as an InputBox to be used
        for an interact canvas.  Basically this specializes this
        input to be used for a specific function and variable.

        INPUT:
            var -- a string (variable; one of the variable names input to f)

        OUTPUT:
            InputBox -- an InputBox object.

        EXAMPLES:
            sage: input_box("2+2", 'Exp').render('x')
            An InputBox interactive control with x='2+2' and label 'Exp'
        """
        if self.__type is Color:
            return ColorInput(var, default_value=self.__default, label=self.label(), type=self.__type)
        else:
            return InputBox(var, default_value=self.__default, label=self.label(), type=self.__type, width = self.__width)


class input_grid(control):
    def __init__(self, nrows, ncols, default=None, label=None, to_value=lambda x: x, width=4):
        r"""
        An input grid interactive control.  Use this in conjunction
        with the interact command.

        INPUT:
            nrows -- integer
            ncols -- integer
            default -- object; the default put in this input box
            label -- the label rendered to the left of the box.
            to_value -- the grid output (list of rows) is sent through
                        this function.  This may reformat the data or
                        coerce the type.
            width -- size of each input box in characters

        NOTEBOOK EXAMPLE:
            @interact
            def _(m = input_grid(2,2, default = [[1,7],[3,4]],
                                 label='M=', to_value=matrix),
                  v = input_grid(2,1, default=[1,2],
                                 label='v=', to_value=matrix)):
                try:
                    x = m\v
                    html('$$%s %s = %s$$'%(latex(m), latex(x), latex(v)))
                except:
                    html('There is no solution to $$%s x=%s$$'%(latex(m), latex(v)))


        EXAMPLES:
            sage: input_grid(2,2, default = 0, label='M')
            Interact 2 x 2 input grid control labeled M with default value 0
            sage: input_grid(2,2, default = [[1,2],[3,4]], label='M')
            Interact 2 x 2 input grid control labeled M with default value [[1, 2], [3, 4]]
            sage: input_grid(2,2, default = [[1,2],[3,4]], label='M', to_value=MatrixSpace(ZZ,2,2))
            Interact 2 x 2 input grid control labeled M with default value [[1, 2], [3, 4]]
            sage: input_grid(1, 3, default=[[1,2,3]], to_value=lambda x: vector(flatten(x)))
            Interact 1 x 3 input grid control labeled None with default value [[1, 2, 3]]

        """
        self.__default = default
        self.__rows = nrows
        self.__columns = ncols
        self.__to_value = to_value
        self.__width = width
        control.__init__(self, label)

    def __repr__(self):
        """
        Return print representation of this input box.

        EXAMPLES:
            sage: input_grid(2,2, label='M').__repr__()
            'Interact 2 x 2 input grid control labeled M with default value None'

        """

        return 'Interact %r x %r input grid control labeled %s with default value %s'%( self.__rows,
                                self.__columns, self.label(),  self.default())


    def default(self):
        """
        Return the default value of this input grid.

        EXAMPLES:
            sage: input_grid(2,2, default=1).default()
            1
        """
        return self.__default


    def render(self, var):
        r"""
        Return rendering of this input grid as an InputGrid to be used
        for an interact canvas.  Basically this specializes this
        input to be used for a specific function and variable.

        INPUT:
            var -- a string (variable; one of the variable names input to f)

        OUTPUT:
            InputGrid -- an InputGrid object.

        EXAMPLES:
            sage: input_grid(2,2).render('x')
            A 2 x 2 InputGrid interactive control with x=[[None, None], [None, None]] and label 'x'

        """
        return InputGrid(var, rows=self.__rows, columns=self.__columns, default_value=self.__default, label=self.label(), to_value=self.__to_value, width=self.__width)



class checkbox(input_box):
    def __init__(self, default=True, label=None):
        """
        A checkbox interactive control.  Use this in conjecture
        with the interact command.

        INPUT:
            default -- bool (default: True); whether box should be checked or not
            label -- str or None (default: None) text label rendered to the left of the box

        EXAMPLES:
            sage: checkbox(False, "Points")
            Interact checkbox labeled 'Points' with default value False
            sage: checkbox(True, "Points")
            Interact checkbox labeled 'Points' with default value True
            sage: checkbox(True)
            Interact checkbox labeled None with default value True
            sage: checkbox()
            Interact checkbox labeled None with default value True
        """
        input_box.__init__(self, bool(default), label=label, type=bool)

    def __repr__(self):
        """
        Print representation of this checkbox.

        EXAMPLES:
            sage: checkbox(True, "Points").__repr__()
            "Interact checkbox labeled 'Points' with default value True"
        """
        return "Interact checkbox labeled %r with default value %r"%(self.label(), self.default())

class slider_generic(control):
    def __init__(self, vmin, vmax=None, step_size=None, label=None, display_value=True):
        control.__init__(self, label=label)
        self.__display_value = display_value
        if isinstance(vmin, list):
            vals = vmin
        else:
            if vmax is None:
                vmax = vmin
                vmin = 0

            #Compute step size; vmin and vmax are both defined here
            #500 is the length of the slider (in px)
            if step_size is None:
                step_size = (vmax-vmin)/499.0
            elif step_size <= 0:
                raise ValueError, "invalid negative step size -- step size must be positive"

            #Compute list of values
            num_steps = int(math.ceil((vmax-vmin)/float(step_size)))
            if num_steps <= 2:
                vals = [vmin, vmax]
            else:
                vals = srange(vmin, vmax, step_size, include_endpoint=True)
                if vals[-1] != vmax:
                    try:
                        if vals[-1] > vmax:
                            vals[-1] = vmax
                        else:
                            vals.append(vmax)
                    except (ValueError, TypeError):
                        pass

        #Is the list of values is small (len<=50), use the whole list.
        #Otherwise, use part of the list.
        if len(vals) == 0:
            self.__values = [0]
        elif(len(vals)<=500):
            self.__values = vals
        else:
	    vlen = (len(vals)-1)/499.0
	    self.__values = [vals[(int)(i*vlen)] for i in range(500)]

    def values(self):
        """
        Returns list of values that this slider takes on, in order.

        OUTPUT:
            list -- list of values

        WARNING: This is a reference to a mutable list.

        EXAMPLES:
            sage: sage.server.notebook.interact.slider(1,10,1/2).values()
            [1, 3/2, 2, 5/2, 3, 7/2, 4, 9/2, 5, 11/2, 6, 13/2, 7, 15/2, 8, 17/2, 9, 19/2, 10]
        """
        return self.__values

    def display_value(self):
        """
        Returns whether to display the value on the slider.

        OUTPUT:
            boolean

        EXAMPLES:
            sage.server.notebook.interact.slider_generic(1,10,1/2).display_value()
            True
        """
        return self.__display_value;


class slider(slider_generic):
    def __init__(self, vmin, vmax=None, step_size=None, default=None, label=None, display_value=True):
        r"""
        An interactive slider control, which can be used in conjunction
        with the interact command.

        \code{slider(vmin, vmax=None, step_size=1, default=None, label=None)}

        INPUT:
            vmin -- object or number
            vmax -- object or None; if None then vmin must be a list, and the slider
                    then varies over elements of the list.
            step_size -- integer (default: 1)
            default -- object or None; default value is ``closest'' in vmin or range
                       to this default.
            label -- string
            display_value -- boolean, whether to display the current value to the right
                             of the slider

        EXAMPLES:
        We specify both vmin and vmax.  We make the default 3, but
        since 3 isn't one of 3/17-th spaced values between 2 and 5,
        52/17 is instead chosen as the default (it is closest).
            sage: slider(2, 5, 3/17, 3, 'alpha')
            Slider: alpha [2--|52/17|---5]

        Here we give a list:
            sage: slider([1..10], None, None, 3, 'alpha')
            Slider: alpha [1--|3|---10]

        The elements of the list can be anything:
            sage: slider([1, 'x', 'abc', 2/3], None, None, 'x', 'alpha')
            Slider: alpha [1--|x|---2/3]
        """
        slider_generic.__init__(self, vmin, vmax, step_size, label, display_value)

        # determine the best choice of index into the list of values
        # for the user-selected default.
        if default is None:
            self.__default = 0
        else:
            try:
                i = self.values().index(default)
            except ValueError:
                # here no index matches -- which is best?
                try:
                    v = [(abs(default - self.values()[j]), j) for j in range(len(self.values()))]
                    m = min(v)
                    i = m[1]
                except TypeError: # abs not defined on everything, so give up
                    i = 0
            self.__default = i

    def __repr__(self):
        """
        Return string representation of this slider.

        EXAMPLES:
            sage: slider(2, 5, 1/5, 3, 'alpha').__repr__()
            'Slider: alpha [2--|3|---5]'
        """
        return "Slider: %s [%s--|%s|---%s]"%(self.label(),
                  self.values()[0],
             self.values()[self.default_index()], self.values()[-1])


    def default_index(self):
        """
        Return default index into the list of values.

        OUTPUT:
            int

        EXAMPLES:
            sage: slider(2, 5, 1/2, 3, 'alpha').default_index()
            2
        """
        return self.__default

    def render(self, var):
        """
        Render the interact control for the given function and
        variable.

        INPUT:
            var -- string; variable name

        EXAMPLES:
            sage: S = slider(0, 10, 1, default=3, label='theta'); S
            Slider: theta [0--|3|---10]
            sage: S.render('x')
            Slider Interact Control: theta [0--|3|---10]

            sage: slider(2, 5, 2/7, 3, 'alpha').render('x')
            Slider Interact Control: alpha [2--|20/7|---5]
        """
        return Slider(var, self.values(), self.__default, label=self.label(), display_value=self.display_value())


class range_slider(slider_generic):
    def __init__(self, vmin, vmax=None, step_size=None, default=None, label=None, display_value=True):
        r"""
        An interactive range slider control, which can be used in conjunction
        with the interact command.

        \code{range_slider(vmin, vmax=None, step_size=1, default=None, label=None)}

        INPUT:
            vmin -- object or number
            vmax -- object or None; if None then vmin must be a list, and the slider
                    then varies over elements of the list.
            step_size -- integer (default: 1)
            default -- (object, object) or None; default range is ``closest'' in vmin or range
                       to this default.
            label -- string
            display_value -- boolean, whether to display the current value below
                             the slider

        EXAMPLES:
        We specify both vmin and vmax.  We make the default (3,4) but
        since neither is one of 3/17-th spaced values between 2 and 5,
        the closest values: 52/17 and 67/17, are instead chosen as the
        default.
            sage: range_slider(2, 5, 3/17, (3,4), 'alpha')
            Range Slider: alpha [2--|52/17==67/17|---5]

        Here we give a list:
            sage: range_slider([1..10], None, None, (3,7), 'alpha')
            Range Slider: alpha [1--|3==7|---10]
        """
        slider_generic.__init__(self, vmin, vmax, step_size, label, display_value)

        # determine the best choice of index into the list of values
        # for the user-selected default.
        if default is None:
            self.__default = (0, 1)
        elif not isinstance(default,tuple) or len(default)!=2:
            raise TypeError("default value must be None or a 2-tuple.")
        else:
            dlist = []
            for i in [0, 1]:
                try:
                    d = self.values().index(default[i])
                except ValueError:
                    # here no index matches -- which is best?
                    try:
                        v = [(abs(default[i] - self.values()[j]), j) for j in range(len(self.values()))]
                        m = min(v)
                        d = m[1]
                    except TypeError: # abs not defined on everything, so give up
                        d = 0
                dlist.append(d)
            self.__default = (dlist[0], dlist[1])

    def __repr__(self):
        """
        Return string representation of this slider.

        EXAMPLES:
            sage: range_slider(2, 5, 1/5, (3,4), 'alpha').__repr__()
            'Range Slider: alpha [2--|3==4|---5]'
        """
        return "Range Slider: %s [%s--|%s==%s|---%s]"%(self.label(), self.values()[0],
             self.values()[self.default_index()[0]],
             self.values()[self.default_index()[1]], self.values()[-1])

    def default_index(self):
        """
        Return default index into the list of values.

        OUTPUT:
            (int, int)

        EXAMPLES:
            sage: range_slider(2, 5, 1/2, (3,4), 'alpha').default_index()
            (2, 4)
        """
        return self.__default

    def render(self, var):
        """
        Render the interact control for the given function and
        variable.

        INPUT:
            var -- string; variable name

        EXAMPLES:
            sage: S = range_slider(0, 10, 1, default=(3,7), label='theta'); S
            Range Slider: theta [0--|3==7|---10]
            sage: S.render('x')
            Range Slider Interact Control: theta [0--|3==7|---10]

            sage: range_slider(2, 5, 2/7, (3,4), 'alpha').render('x')
            Range Slider Interact Control: alpha [2--|20/7==4|---5]
        """
        return RangeSlider(var, self.values(), self.__default, label=self.label(), display_value=self.display_value())


class selector(control):
    def __init__(self, values, label=None, default=None,
                 nrows=None, ncols=None, width=None, buttons=False):
        r"""
        A drop down menu or a button bar that when pressed sets a
        variable to a given value.  Use this in conjunction with the
        interact command.

        \code{selector(values, label=None, nrows=None, ncols=None, buttons=False)}

        We use the same command to create either a drop down menu or
        selector bar of buttons, since conceptually the two controls
        do exactly the same thing -- they only look different.  If
        either nrows or ncols is given, then you get a buttons instead
        of a drop down menu.

        INPUT:
            values -- [val0, val1, val2, ...] or
                      [(val0, lbl0), (val1,lbl1), ...] where all labels must be
                                                       given or given as None.
            label -- (default: None); if given, this label is placed to
                                      the left of the entire button group
            default -- object (default: 0) default value in values list
            nrows -- (default: None); if given determines the number
                     of rows of buttons; if given buttons option below is set to True
            ncols -- (default: None); if given determines the number
                     of columns of buttons; if given buttons option below is set to True
            width -- (default: None); if given, all buttons are the same
                     width, equal to this in html ex units's.
            buttons -- (default: False); if True, use buttons

        EXAMPLES:
            sage: selector([1..5])
            Drop down menu with 5 options
            sage: selector([1,2,7], default=2)
            Drop down menu with 3 options
            sage: selector([1,2,7], nrows=2)
            Button bar with 3 buttons
            sage: selector([1,2,7], ncols=2)
            Button bar with 3 buttons
            sage: selector([1,2,7], width=10)
            Drop down menu with 3 options
            sage: selector([1,2,7], buttons=True)
            Button bar with 3 buttons

        We create an interactive that involves computing charpolys of matrices over various rings:
            sage: @interact
            ... def _(R=selector([ZZ,QQ,GF(17),RDF,RR]), n=(1..10)):
            ...      M = random_matrix(R, n)
            ...      show(M)
            ...      show(matrix_plot(M,cmap='Oranges'))
            ...      f = M.charpoly()
            ...      print f
            <html>...

        Here we create a drop-down
            sage: @interact
            ... def _(a=selector([(2,'second'), (3,'third')])):
            ...       print a
            <html>...
        """
        if nrows is not None or ncols is not None:
            buttons = True
        if default is None:
            default=0
        else:
            try:
                default = values.index(default)
            except IndexError:
                default = 0
        self.__values = values
        self.__nrows = nrows
        self.__ncols = ncols
        self.__width = width
        self.__default = default
        self.__buttons = buttons
        control.__init__(self, label)

    def __repr__(self):
        """
        Return print representation of this button.

        EXAMPLES:
            sage: selector([1,2,7], default=2).__repr__()
            'Drop down menu with 3 options'
        """
        if self.__buttons:
            return "Button bar with %s buttons"%len(self.__values)
        else:
            return "Drop down menu with %s options"%len(self.__values)

    def values(self):
        """
        Return the list of values or (val, lbl) pairs that this
        selector can take on.

        OUTPUT:
            list

        EXAMPLES:
            sage: selector([1..5]).values()
            [1, 2, 3, 4, 5]
            sage: selector([(5,'fifth'), (8,'eight')]).values()
            [(5, 'fifth'), (8, 'eight')]
        """
        return self.__values

    def default(self):
        """
        Return the default choice for this control.

        OUTPUT:
           int -- an integer, with 0 corresponding to the first choice.

        EXAMPLES:
            sage: selector([1,2,7], default=2).default()
            1
        """
        return self.__default

    def render(self, var):
        r"""
        Return rendering of this button as a Button instance to be
        used for an interact canvas.

        INPUT:
            var -- a string (variable; one of the variable names input to f)

        OUTPUT:
            Button -- a Button instance

        EXAMPLES:
            sage: selector([1..5]).render('alpha')
            Selector with 5 options for variable 'alpha'
        """
        return Selector(var, values=self.__values, label=self.label(),
                        default=self.__default,
                        nrows=self.__nrows, ncols=self.__ncols, width=self.__width,
                        buttons=self.__buttons)

class text_control(control):
    def __init__(self, value=''):
    	"""
    	Text that can be inserted among other interact controls.

        INPUT:
            value -- HTML for the control

        EXAMPLES:
            sage: text_control('something')
            Text field: something
    	"""
        self.__default = value
        control.__init__(self, '')

    def __repr__(self):
        """
        Return print representation of this control.

        EXAMPLES:
            sage: text_control('something')
            Text field: something
        """
        return "Text field: %s"%self.__default

    def render(self, var):
        """
        Return rendering of the text field

        INPUT:
            var -- a string (variable; one of the variable names input to f)

        OUTPUT:
            TextControl -- a TextControl instance
        """
        return TextControl(var, self.__default)


def automatic_control(default):
    """
    Automagically determine the type of control from the default
    value of the variable.

    INPUT:
        default -- the default value for v given by the function; see
                   the documentation to interact? for details.

    OUTPUT:
        a interact control

    EXAMPLES:
        sage: sage.server.notebook.interact.automatic_control('')
        Interact input box labeled None with default value ''
        sage: sage.server.notebook.interact.automatic_control(15)
        Interact input box labeled None with default value 15
        sage: sage.server.notebook.interact.automatic_control(('start', 15))
        Interact input box labeled 'start' with default value 15
        sage: sage.server.notebook.interact.automatic_control((1,250))
        Slider: None [1.0--|1.0|---250.0]
        sage: sage.server.notebook.interact.automatic_control(('alpha', (1,250)))
        Slider: alpha [1.0--|1.0|---250.0]
        sage: sage.server.notebook.interact.automatic_control((2,(0,250)))
        Slider: None [0.0--|2.00400801603|---250]
        sage: sage.server.notebook.interact.automatic_control(('alpha label', (2,(0,250))))
        Slider: alpha label [0.0--|2.00400801603|---250]
        sage: sage.server.notebook.interact.automatic_control((2, ('alpha label',(0,250))))
        Slider: alpha label [0.0--|2.00400801603|---250]
        sage: C = sage.server.notebook.interact.automatic_control((1,52, 5)); C
        Slider: None [1--|1|---52]
        sage: C.values()
        [1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 52]
        sage: sage.server.notebook.interact.automatic_control((17, (1,100,5)))
        Slider: None [1--|16|---100]
        sage: sage.server.notebook.interact.automatic_control([1..4])
        Button bar with 4 buttons
        sage: sage.server.notebook.interact.automatic_control([1..100])
        Drop down menu with 100 options
        sage: sage.server.notebook.interact.automatic_control((1..100))
        Slider: None [1--|1|---100]
        sage: sage.server.notebook.interact.automatic_control((5, (1..100)))
        Slider: None [1--|5|---100]
        sage: sage.server.notebook.interact.automatic_control(matrix(2,2))
        Interact 2 x 2 input grid control labeled None with default value [0, 0, 0, 0]
    """
    label = None
    default_value = None

    for _ in range(2):
        if isinstance(default, tuple) and len(default) == 2 and isinstance(default[0], str):
            label, default = default
        if isinstance(default, tuple) and len(default) == 2 and isinstance(default[1], (tuple, list, types.GeneratorType)):
            default_value, default = default

    if isinstance(default, control):
        C = default
        if label:
            C.set_label(label)
    elif isinstance(default, str):
        C = input_box(default, label=label, type=str)
    elif isinstance(default, bool):
        C = input_box(default, label=label, type=bool)
    elif isinstance(default, list):
        C = selector(default, default=default_value, label=label, buttons=len(default) <= 5)
    elif isinstance(default, types.GeneratorType):
        C = slider(list_of_first_n(default,10000), default=default_value, label=label)
    elif isinstance(default, Color):
        C = input_box(default, label=label, type=Color)
    elif isinstance(default, tuple):
        if len(default) == 2:
            C = slider(default[0], default[1], default=default_value, label=label)
        elif len(default) == 3:
            C = slider(default[0], default[1], default[2], default=default_value, label=label)
        else:
            C = slider(list(default), default=default_value, label=label)
    elif is_Matrix(default):
        C = input_grid(default.nrows(), default.ncols(), default=default.list(), to_value=default.parent())
    else:
        C = input_box(default, label=label)

    return C

def list_of_first_n(v,n):
    """
    Given an iterator v, return first n elements it produces as a list.

    INPUT:
        v -- an interator
        n -- an integer

    OUTPUT:
        list

    EXAMPLES:
        sage: sage.server.notebook.interact.list_of_first_n(Primes(), 10)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        sage: sage.server.notebook.interact.list_of_first_n((1..5), 10)
        [1, 2, 3, 4, 5]
        sage: sage.server.notebook.interact.list_of_first_n(QQ, 10)
        [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3]
    """
    if not hasattr(v, 'next'):
        v = v.__iter__()
    w = []
    while n > 0:
        try:
            w.append(v.next())
        except StopIteration:
            return w
        n -= 1
    return w


def update(cell_id, var, adapt, value, globs):
    """
    Called when updating the positions of an interactive control.
    Note that this just causes the values of the variables to be
    updated; it does not reevaluate the function with the new values.

    INPUT:
        cell_id -- the id of a interact cell
        var -- a variable associated to that cell
        adapt -- the number of the adapt function
        value -- new value of the control
        globs -- global variables.

    EXAMPLES:
    The following outputs __SAGE_INTERACT_RESTART__ to indicate that
    not all the state of the interrupt canvas has been setup yet (this
    setup happens when javascript calls certain functions).
        sage: sage.server.notebook.interact.update(0, 'a', 0, '5', globals())
        __SAGE_INTERACT_RESTART__
    """
    try:
        S = state[cell_id]
        # Look up the function that adapts inputs to have the right type
        adapt_function = S["adapt"][adapt]
        # Apply that function and save the result in the appropriate variables dictionary.
        S["variables"][var] = adapt_function(value, globs)
    except KeyError:
        # If you change this, make sure to change js.py as well.
        print "__SAGE_INTERACT_RESTART__"

def recompute(cell_id):
    """
    Evaluates the interact function associated to the cell
    cell_id. This typically gets called after a call to
    sage.server.notebook.interact.update.

    EXAMPLES:
    The following outputs __SAGE_INTERACT_RESTART__ to indicate that
    not all the state of the interrupt canvas has been setup yet (this
    setup happens when javascript calls certain functions).
        sage: sage.server.notebook.interact.recompute(10)
        __SAGE_INTERACT_RESTART__

    """
    try:
        S = state[cell_id]
        # Finally call the interactive function, which will use the above variables.
        S['function']()
    except KeyError:
        # If you change this, make sure to change js.py as well.
        print "__SAGE_INTERACT_RESTART__"

