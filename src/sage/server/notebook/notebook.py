r"""
SAGE Notebook Interface

AUTHORS:
    -- William Stein (2006-05-06): initial version
    -- Alex Clemesha
    -- Tom Boothby

SUPPORTED BROWSERS:
The SAGE notebook works with each of the following browsers:
\begin{itemize}
   \item Firefox
   \item Internet Explorer
\end{itemize}

  Note that Safari, Opera, and Konqueror are \emph{not} supported.

TUTORIAL:

Here are some things to try in the the notebook to get a feeling
for it (i.e., this is a very quick 3-minute tutorial):

Type "2+2" in the blank box and press "shift-enter" (like in Mathematica).
The line below "2+2" will turn reddish for a moment while SAGE fires
up and computes the answer.

Your cursor should now be in the next box down.   Type \code{a = 2^1000}
and press return, then "a" alone on the second line, then shift-return.
You'll see a big number.   Also, "a" will appear in the variable
browser in the upper left corner of the screen.  (Question -- what
should go in the variable browser.  I could make mousing over a variable
display all kinds of info about it... but what info?  It's a design
question, and feedback would be welcome.)    Next, click just to the
left of the big number in the blue-ish area.  The number will shrink
a little and change to occupy only one line.  You can see the whole
number using your browser's horizontal scroll bar.  Click again and
the number vanishes, to be replaced by a grey bar.  Click on the grey
bar and the number is back.  If you click "Hide Output" in the upper
right, all output disappears.

Next try graphics!  Type "show(plot(sin,0,10))" into an empty
box and hit shift-enter.   You'll get a graph of sin.   Try another
function, e.g.,
     \code{show(plot(lambda x: sin(x)^2 - cos(2*x)^3, -5,5))}
Click on the left side of the figure (twice) to make it disappear.

One really cool new feature of the SAGE notebook, is that you can
"queue up" a bunch of calculations in a row, *while* still editing the
notebook!  This is basically just like in Mathematica, and already
works pretty well.   (If you have a repeatable way to mess this up,
please let me know, since I don't know any at all right now.)
As an example, consider computing factorials, which takes a while
but not too long.  First, enter the following in a blank box and
press "shift-return":
         \code{ def f(n): return len(str(factorial(10^n)))}
This defines a function that takes a while to compute.   For example,
time the execution of "f(5)", by typing (in a new box), "f(5)", then
pressing "ctrl-return".  It should take a few seconds.   Next try
"f(6)", which takes quite a while (about 21 seconds on sage.math).
While f(6) is being computed, note 2 things:
   (a) the output line for f(6) is in light red, indicating that
       it is being computed
   (b) The Interrupt link in the upper right is not greyed out, i.e.,
       you can click and interrupt the computation.
While f(6) is computing (if it finishes first, restart it by
just hitting ctrl-enter in the box where "f(6)" is), try typing
"f(4)" in the next box.  You're allowed to give input, but the
result doesn't get computed immediately.  You can enter several more
lines as well, etc.; when the f(6) finally finishes, SAGE goes on
to compute "f(4)".   You can queue up dozens of calculations.  For
example, if you hit the "Evaluate All" link in the upper right, the
whole workbook is queued up for computation.  Try it.   When the
computation gets stuck on "f(6)", hit the interrupt button and
all the queued up calculations are cancelled.

Click "Hide Output" in the upper right.   You'll see just your
input and grey boxes; clicking on the grey boxes reveals output.

You can also embed nicely typeset math.  Try this:
\begin{verbatim}
   f = maxima('sin(x^2)')
   g = f.integrate('x')
   view(g)
\end{verbatim}

If this silently fails, type "view(g, debug=True)" instead.
You need latex and the "convert" and "gs" commands, which is just
an "apt-get install imagemagick gs" away...  Anyways, you get
a nicely typeset formula.  Try a matrix next:
   A = MatrixSpace(QQ, 5).random_element()
   view(A)

ADDING AND REMOVING CELLS:
How to add and remove input cells is an interesting design
decision.  I went with the following:
   (1) To add a new cell, click on a little black line that
       appears when you hover between any two cells, or above
       the top one.  This is what mathematica does.
   (2) To delete a cell delete all its contents, then
       hit backspace.  The cell vanishes forever.

You can also move back and forth between cells using
"tab" and "shift tab".  Right now, inserting new cells is
really the only thing in the interface that requires using
a mouse.  I could make it so there is some sort of keyboard
shortcut to insert a new cell... though there aren't many
options, since the web browser takes most of them already.

Anyway, please try the above out for adding and remove
cells, and let me know if it is usable.  There
is currently no support for moving and
reorganizing cells, and I'm not sure there should be, given
the somewhat linear nature of interaction with SAGE.
The underlying design of the system is such that adding
the ability to reorder cells wouldn't be too difficult.


INTROSPECTION:
The design for introspection, i.e., finding all completions,
and finding information information about an object is very
interesting.   Tom, Alex and I made a number of prototypes
that involve a little input line in the upper left corner
of the screen whose content varies to reflect all possible completions
of what you type into it.   In the end, I've settle on something
much different, since it seems very easy to use in practice, and
is more like IPython, so more familiar, and also better if you
don't like using a mouse.

  (1) to find all completions for the last identifier you're
   typing on a line, just hit "escape".  The completions are
   listed in the output box.  In particular, if you hit
   escape in a blank cell you get a list of all sage commands.
   If then type "Ab" followed by hitting escape, you get
   the abelian group commands.  Type one in and put a question
   mark afterwards and hit escape for help (see next step).

  (2) To find help for any object in a line, put a question
    mark after it and press "escape".  The object must exist
    in the current scope, e.g., so
      a = 5
      a?
    all in one box won't work.

  (3) To get source code about any object, put ?? after it
    and press escape.  The help appears below in the output
    box, but you can continue to type in your input box.

  (4) To get extensive help on an object, type "help(object)"
    and press return.  This works, since (a) I set the
    PAGER to "cat", and (b) I strip out control codes that
    appear in the output.  And this isn't annoying, since web
    browsers are very good for scrolling through long output.

  (5) To find all completions for a word anywhere in a line,
     type "?c" right after the word and press "escape".

  Try the above and give me feedback.  How does this feel to you?

OBJECTS:  When you start a notebook you give a name argument
to it, and it creates a directory.  Inside that directory there
will be many workbooks (which you can use all at once and easily
flip through -- not implemented yet), and an object store.
You can save and load objects (using save and load), and they'll
be listed in the box on the bottom let, e.g., try
   a = 5
   save(a, 'a')
and you'll see the "a" appear there.   You can load and save objects
from any workbook in any other one.

PASTING in EXAMPLES:

Code is evaluated by exec'ing (after preparsing). Only the output
of the last line of the cell is implicitly printed. If any line
starts with "sage:" or ">>>" the entire block is assumed to
contain text and examples, so only lines that begin with a
prompt are executed. Thus you can paste in *complete examples*
from the docs without any editing, and you can write input
cells that contains non-evaluated plain text mixed with
examples by starting the block with ">>>" or including an example.

SAVING:

The SAGE notebook is very persistent.  Every time you submit
a cell for computation, the state of the notebook is saved (a
few kb's file).  If you quit the notebook and reload, it will
have everything you typed from the previous session, along
with all output.   Todo: I'll make the saved state file backup
the last safe state, to avoid the potential for corruption
if the server Python process is killed while saving state.
Also, this could easily allow for a sophisticated undo function;
design ideas welcome...

ARCHITECTURE:

The SAGE Notebook is an AJAX application that can run either
entirely locally on your desktop machine, or partly on
a server and via a web browser that could be located somewhere
else.  There is currently no support for multiple connections, but
there will be as soon as I get multiple workbooks going (which
is already mostly there, since it was designed in from the start).
Anywhere, here are the components of the SAGE Notebook:

 (1) WEB SERVER:
     A Python process that uses the Python standard library's
     BaseHTTPServer.HTTPServer to create a web server.  This
     process also handles all requests from the web browser,
     e.g., organizing computation of cells, etc.  It doesn't
     do any actual mathematics, and only imports a very small
     subset of the SAGE library.  In particular, if you do
     "sage -notebook" at the command line, very little of
     SAGE is imported.  This could be separated out
     from SAGE to be its own Python module, in case people are
     generally interested in "a web version of IPython"...

 (2) SAGE SERVER:
     A Python process with all the SAGE libraries loaded; this
     is started by (1) when a web browser first requests that
     a cell be evaluated.  There's one of these Python processes
     for each workbook.

 (3) WEB BROWSER:
     The web browser runs a 500-line javascript (and 600 lines of css)
     program that Alex, Tom and I wrote, which implements a lot of the
     browser-side part of the SAGE notebook functionality.

When you use the SAGE Notebook, you are mainly interacting with a
javascript program.  When you do something serious, e.g., request
computation of some input, create a new cell, etc., a request
is made from your web browser to the web server telling it what
is going on.    If it's a calculation, the web server tells the
SAGE server to get started on the calculation, and tells the web
browser to check every delta (=0.2 seconds right now) whether
there is anything new with the calculation.   Your web browser
then checks back 5 times a second for updates, and when something
new appears it fills that in.  This continues until all calculations
are done. During this time, you can edit cells, submit more computations,
etc.   Note that output really is updated in real time.  For
example, try the following from the SAGE Notebook:

import time
for i in range(10):
    print i
    time.sleep(0.5)

You get to watch as the integers from 1 to 10 are "computed" in
real time.    Actually, getting it so output is reported in
real time is, I think, *crucial* to making a really usable
SAGE GUI -- users (i.e., me) want to run huge computations
and watch the output progress.

The architecture is also good from the point of view of being
able to interrupt running computations.  What happens when
you request an interrupt is that the web browser sends a message
to the web server, which in turn tells the SAGE server to
stop computing by sending it control-c's (up to about 10 seconds)
until it either stops, or if it's really frozen (due to a bug,
or calling into a C function that isn't properly wrapped
in signal handling, or maybe you run some crazy interactive
program like pine via "os.system('...')"), it'll just kill
that SAGE server 2 and restart it.  The result is that
the user doesn't get a frozen web browser or browser interface
at any point, and even if the whole SAGE process went down and
froze, at least all your input and output from your session is
still there in your browser.  The only thing you've lost is
the definition of all your variables.  Hit "shift-enter"
a few times or "evaluate all" and you're back in shape.
This is much better than having to restart the command prompt
(e.g., with a terminal interface), then paste back in all your
setup code, etc.,   Also, you can save variables as you go
easily (via the "save" command), and get back to where you
were quickly.    -- This is probably very similar to Mathematica/Maple
which have "restart kernel" commands, etc.

"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import os
import socket

# SAGE libraries
from   sage.ext.sage_object import SageObject, load
from   sage.misc.viewer     import BROWSER
from   sage.misc.misc       import alarm

# SAGE Notebook
import css          # style
import js           # javascript
import server       # web server
import workbook     # individual workbooks (which make up a notebook)

MAX_WORKBOOKS = 65535

class Notebook(SageObject):
    def __init__(self, dir='sage_notebook'):
        self.__dir = dir
        self.__workbooks = {}
        self.__load_defaults()
        self.__filename     = '%s/nb.sobj'%dir
        self.__workbook_dir = '%s/workbooks'%dir
        self.__object_dir   = '%s/objects'%dir
        self.__makedirs()
        self.__next_workbook_id = 0
        W = self.create_new_workbook('*scratch*')
        self.__default_workbook = W
        self.create_new_workbook('linear algebra')
        self.create_new_workbook('elliptic curves')
        self.create_new_workbook('modular forms')
        self.create_new_workbook('elliptic curves modular forms, and fermat')
        self.create_new_workbook('test of multi workbooks')
        self.create_new_workbook('many new ideas')
        self.create_new_workbook('l-series')
        self.create_new_workbook('modular degrees')
        self.create_new_workbook('Dirichlet Characters')
        self.create_new_workbook('SAGE tests')
        self.save()

    def default_workbook(self):
        return self.__default_workbook

    def directory(self):
        return self.__dir

    def DIR(self):
        """
        Return the absolute path to the directory that contains
        the SAGE Notebook directory.
        """
        return os.path.abspath('%s/..'%self.__dir)

    def __load_defaults(self):
        # in future this will allow override by a file, and
        # can be set by user via web interface
        self.__defaults = {'cell_input_color':'#0000000',
                           'cell_output_color':'#0000EE',
                           'word_wrap_cols':85}

    def workbook_directory(self):
        return self.__workbook_dir

    def object_directory(self):
        return self.__object_dir

    def objects(self):
        L = [x[:-5] for x in os.listdir(self.__object_dir)]
        L.sort()
        return L

    def object_list_html(self):
        m = max([len(x) for x in self.objects()] + [30])
        s = []
        a = '<a href="/%s.sobj" class="object_name">'
        for name in self.objects():
            s.append(a%name + name + '&nbsp;'*(m-len(x)) + '</a>')
        return '<br>'.join(s)

    def defaults(self):
        return self.__defaults

    def __makedirs(self):
        os.makedirs(self.__dir)
        os.makedirs(self.__workbook_dir)
        os.makedirs(self.__object_dir)

    def create_new_workbook(self, name='untitled'):
        if name in self.__workbooks.keys():
            raise ValueError, 'name (=%s) already taken.'%name
        name = str(name)
        id = self.__next_workbook_id
        if id >= MAX_WORKBOOKS:
            raise ValueError, 'there can be at most %s workbooks'%MAX_WORKBOOKS
        self.__next_workbook_id += 1
        W = workbook.Workbook(name, self, id)
        self.__workbooks[name] = W
        return W

    def workbook_names(self):
        W = self.__workbooks.keys()
        W.sort()
        return W

    def get_workbook_with_id(self, id):
        for W in self.__workbooks.itervalues():
            if W.id() == id:
                return W
        raise KeyError, "no workbook with id %s"%id

    def get_workbook_that_has_cell_with_id(self, id):
        workbook_id = id // MAX_WORKBOOKS
        return self.get_workbook_with_id(workbook_id)

    def save(self, filename=None):
        if filename is None:
            SageObject.save(self, os.path.abspath(self.__filename))
        else:
            SageObject.save(self, os.path.abspath(filename))

    def start(self, port=8000, address='localhost',
                    max_tries=10, open_viewer=False):
        tries = 0
        while True:
            try:
                ns = server.NotebookServer(self, port, address)
            except socket.error, msg:
                print msg
                port += 1
                tries += 1
                if tries > max_tries:
                    print "Not trying any more ports.  Probably your network is down."
                    break
                print "Trying next port (=%s)"%port
            else:
                break

        s = "SAGE Notebook http://%s:%s"%(address, port)
        t = len(s)
        if t%2:
            t += 1
            s += ' '
        n = max(t+4, 50)
        k = n - t  - 1
        j = k/2
        print '*'*n
        print '*'+ ' '*(n-2) + '*'
        print '*' + ' '*j + s + ' '*j + '*'
        print '*'+ ' '*(n-2) + '*'
        print '*'*n
        print "WARNING!!! Currently the SAGE Notebook *only* works with Firefox."

        if open_viewer:
            os.system('%s http://%s:%s 1>&2 >/dev/null &'%(BROWSER, address, port))
        ns.serve()
        self.save()

    def workbook_list_html(self, current_workbook):
        s = []
        names = self.workbook_names()
        m = max([len(x) for x in names] + [30])
        for n in names:
            W = self.__workbooks[n]
            if W == current_workbook:
                cls = 'workbook_current'
            else:
                cls = 'workbook_other'
            if W.computing():
                cls += '_computing' # actively computing
            name = W.name().replace(' ','&nbsp;')
            txt = '<a class="%s" onClick="switch_to_workbook(%s)" href="/%s">%s</a>'%(
                cls,W.id(),W.id(),name + '&nbsp;'*(m-len(W.name())))
            s.append(txt)
        return '<br>'.join(s)

    def _html_head(self):
        head = '<title>SAGE Notebook: %s</title>'%self.directory()
        head += '<style>' + css.css() + '</style>\n'
        head += '<script language=javascript>' + js.javascript() + '</script>\n'
        return head

    def _html_body(self, workbook_id):
        workbook = self.get_workbook_with_id(workbook_id)
        if workbook.computing():
            interrupt_class = "interrupt"
        else:
            interrupt_class = "interrupt_grey"

        vbar = '<span class="vbar"></span>'

        body = ''
        body += self._help_window()
        body += '<div class="top_control_bar">\n'
        body += '  <span class="banner">SAGE</span>\n'
        body += '  <span class="control_commands">\n'
        body += '    <a class="evaluate" onClick="evaluate_all()">Evaluate All</a>' + vbar
        body += '<a class="%s" onClick="interrupt()" id="interrupt">Interrupt</a>'%interrupt_class + vbar
        body += '<a class="hide" onClick="hide_all()">Hide Output</a>' + vbar
        body += '<a class="help" onClick="show_help_window()">Help</a>'
        body += '  </span>'
        body += '</div>'
        body += '\n<div class="workbook" id="workbook">\n' + workbook.html() + '\n</div>\n'

        body += '<span class="pane"><table bgcolor="white"><tr><td>\n'
        body += '  <div class="variables_topbar">Variables</div>\n'
        body += '  <div class="variables_list" id="variable_list">%s</div>\n'%workbook.variables_html()
        body += '  <div class="attached_topbar">Attached Files</div>\n'
        body += '  <div class="attached_list" id="attached_list">(not implemented)</div><br>\n'
        body += '  <div class="workbooks_topbar">Workbooks&nbsp;&nbsp;&nbsp;'
        body += '     <a onClick="add_new_workbook()" class="new_workbook">Add</a>\n'
        body += '     <a onClick="upload_workbook()" class="upload_workbook">Upload</a>\n'
        body += '  </div>\n'
        body += '  <div class="workbook_list">%s</div>\n'%self.workbook_list_html(workbook)
        body += '  <div class="objects_topbar">Saved Objects</div>\n'
        body += '  <div class="object_list" id="object_list">%s</div>\n'%self.object_list_html()
        body += '</td></tr></table></span>\n'
        body += '<script language=javascript>focus_on(%s)</script>\n'%(workbook[0].id())

        if workbook.computing():
            # Set the update checking back in motion.
            body += '<script language=javascript> check_for_cell_output() </script>\n'

        return body

    def _help_window(self):
        help = [('<b>Evaluate</b> an input cell', 'Press shift-enter.  You can start several calculations at once.'),
                ('Evaluate cell with <b>timing</b>', 'Press control-enter.  The CPU and wall time are printed.'),
                ('Evaluate <b>all</b> cells', 'Click the evaluate all link in the upper right.'),
                ('<b>Move</b> between cells', 'Use tab and shift-tab.'),
                ('<b>Interrupt</b> running calculations',
                 'Click the interrupt button in the upper right or press ctrl-c in any input cell. This will (attempt) to interrupt SAGE by sending control-c for several seconds; if this fails, it restarts SAGE (your workbook remains unchanged).'),
                ('<b>Completions</b> of last object in cell', 'Press the escape key.  Press escape in an empty cell for a list of all commands and objects.'),
                ('<b>Completion</b> of word anywhere in input',
                  'Type ?c immediately after the word and press escape, e.g., "Ell?c"'),
                ('<b>Help</b> about an object',
                 'Put ? right after the object and press escape.'),
                ('<b>Detailed help</b> about an object',
                 'Type "help(name of object)" and press shift-return.'),
                ('<b>Source code</b> of a command',
                 'Put ?? after the command and press escape.'),
                ('<b>Insert</b> a new cell',
                 'Put mouse between an output and input until the horizontal line appears and click.'),
                ('<b>Delete</b> a cell',
                 'Delete cell contents using backspace and the cell will be removed.'),
                ('<b>Hide/Expand</b> Output', 'Click on the left side of output to toggle between hidden, shown with word wrap, and shown without word wrap.'),
                ('<b>Hide</b> Output', 'Click on the hide output button in the upper right to hide <i>all</i> output.'),
                ('<b>Variables</b>',
                 'All variables with a new name that you create during this session are listed on the left.  If you overwrite a predefined variable, e.g., ZZ, it will not appear.'),
                ('<b>Objects</b>',
                 'All objects that you save in <i>any workbook</i> are listed on the left.  Use "save(obj, name)" and "obj = load(name)" to load and save objects.'),
                ('<b>Saving Workbooks</b>',
                 '<i>Everything</i> that has been submitted is automatically saved to disk, and is there for you next time.  You do not have to do anything special to save a workbook.'),
                ('<b>Typeset</b> version of object', 'Type "view(objname)". This requires that latex and convert are installed.  Type "latex(objname)" for latex that you can paste into your paper.'),
                ('<b>Input rules</b>', "Code is evaluated by exec'ing (after preparsing).  Only the output of the last line of the cell is implicitly printed.  If any line starts with \"sage:\" or \">>>\" the entire block is assumed to contain text and examples, so only lines that begin with a prompt are executed.   Thus you can paste in complete examples from the docs without any editing, and you can write input cells that contains non-evaluated plain text mixed with examples by starting the block with \">>>\" or including an example."),
                ('<b>Working directory</b>', 'Each block of code is run from its own directory.  The variable DIR contains the directory from which you started the SAGE notebook; to open a file in that directory, do "open(DIR+\'filename\')".'),
                ('<b>Acknowledgement</b>', 'The design of SAGE notebook was influenced by Mathematica, GMail, GNU Emacs, and IPython.  AUTHORS: William Stein, Alex Clemesha, and Tom Boothy')
                ]
        s = """
        <div class="help_window" id="help_window">
        <div class="help_window_title">&nbsp;&nbsp;&nbsp;Help</div>
        <div class="help_window_close" onClick="hide_help_window()">&#215;&nbsp;</div>
        <table class="help_window">
        """
        for x, y in help:
            s += '<tr><td class="help_window_cmd">%s</td><td class="help_window_how">%s</td></tr>'%(x,y)
        s += '</table></div>'
        return s


    def html(self, workbook_id=None):
        if workbook_id is None:
            W = self.default_workbook()
            workbook_id = W.id()
        else:
            try:
                self.get_workbook_with_id(workbook_id)
            except KeyError:
                W = self.default_workbook()
                workbook_id = W.id()
        head = self._html_head()
        body = self._html_body(workbook_id)
        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        <script language=javascript>workbook_id=%s</script>
        """%(head, body, workbook_id)


def notebook(dir       ='sage_notebook',
             port      = 8000,
             address   = 'localhost',
             open_viewer    = True,
             max_tries = 10):
    r"""
    Start a SAGE notebook web server at the given port.

    INPUT:
        dir -- (default: 'sage_notebook') name of the server directory; your
                sessions are saved in a directory with that name.  If
                you restart the server with that same name then it will
                restart in the state you left it, but with none of the
                variables defined (you have to re-eval blocks).

        port -- (default: 8000) port on computer where the server is served

        address -- (default: 'localhost') address that the server
                   will listen on
        viewer -- bool (default:True); if True, pop up a web browser at the URL
        max_tries -- (default: 10) maximum number of ports > port to try in
                     case given port can't be opened.

    NOTES:

    When you type \code{notebook(...)}  you start a web server on the
    machine you type that command on.  You don't connect to another
    machine.  So do this if you want to start a SAGE notebook
    accessible from anywhere:

    \begin{enumerate}
    \item Figure out the external IP address of your server, say
       128.208.160.191  (that's sage.math's, actually).
    \item On your server, type
       server_http1('mysession', address='128.208.160.191')
    \item Assuming you have permission to open a port on that
       machine, it will startup and display a URL, e.g.,
           \url{http://128.208.160.191:8000}
       Note this URL.
    \item Go to any computer in the world (!), or at least
       behind your firewall, and use any web browser to
       visit the above URL.  You're using \sage.
    \end{enumerate}

    \note{There are no security precautions in place yet.  If
    you open a server as above, and somebody figures this out, they
    could use their web browser to connect to the same sage session,
    and type something nasty like \code{os.system('cd; rm -rf *')}
    and your home directory would be hosed.   I'll be adding an
    authentication screen in the near future.  In the meantime
    (and even then), you should consider creating a user with
    very limited privileges (e.g., empty home directory).}
    """
    print "WARNING -- the SAGE Notebook is currently in alpha"
    print "testing.  It also only looks right on Firefox on"
    print "William Stein's Linux laptop!  Please give me feedback."
    if os.path.exists(dir):
        if not os.path.isdir(dir):
            raise RuntimeError, '"%s" is not a valid SAGE notebook directory (it is not even a directory).'%dir
        if not os.path.exists('%s/nb.sobj'%dir):
            raise RuntimeError, '"%s" is not a valid SAGE notebook directory (missing nb.sobj).'%dir
        nb = load('%s/nb.sobj'%dir)
    else:
        nb = Notebook(dir)

    nb.start(port, address, max_tries, open_viewer)
    alarm(1)
    from sage.interfaces.quit import expect_quitall
    expect_quitall(verbose=False)
    from sage.misc.misc import delete_tmpfiles
    delete_tmpfiles()






