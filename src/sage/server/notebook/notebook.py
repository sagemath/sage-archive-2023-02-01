r"""
SAGE Notebook Interface

AUTHORS:
    -- William Stein (2006-05-06): initial version
    -- Alex Clemesha
    -- Tom Boothby

\subsection{Supported Browsers}

The SAGE notebook currently is fully supported with Firefox only.
Support is planned for Opera, Konqueror and Safari, and Internet
Explorer.

\subsection{Tutorial}
Here are some things to try in the the notebook to get a feeling
for it.

Type "2+2" in the blank box and press "shift-enter".
The line below"2+2" will turn a different color for a moment while a SAGE kernel
fires up and computes the answer.

Your cursor should now be in the next box down.   Type \code{a = 2\^1000}
and press return, then "a" alone on the second line, then shift-return.
You'll see a big number.   Also, "a" will appear in the variable
browser in the left of the screen.    Next, click just to the
left of the big number in the blue-ish area.  The number will shrink
a little and change to occupy only one line.  You can see the whole
number using your browser's horizontal scroll bar.  Click again and
the number vanishes, to be replaced by a horizontal bar.  Click on the
bar and the number is back.  If you click "Hide Output" in the upper
right, all output disappears.

Next try graphics!  Type "show(plot(sin,0,10))" into an empty
box and hit shift-enter.   You'll get a graph of sin.   Try another
function, e.g.,
\begin{verbatim}
   show(plot(lambda x: sin(x)^2 - cos(2*x)^3, -5,5))
\end{verbatim}
Click on the left side of the figure (twice) to make it disappear.

One important feature of the SAGE notebook, is that you can
"queue up" a bunch of calculations in a row, *while* still editing the
notebook!  As an example, consider computing factorials, which takes a
while (but not forever).  First, enter the following in a blank box and
press"shift-return":
\begin{verbatim}
def f(n):
    return len(str(factorial(10^n)))
\end{verbatim}
This defines a function that takes a while to compute.   For example,
time the execution of "f(5)", by typing (in a new box), "time f(5)".
It should take a few seconds.   Next try
"f(6)", which takes quite a while (about 21 seconds on sage.math).
While f(6) is being computed, note that the output line for f(6) is a
different color, indicating that it is being computed
While f(6) is computing (if it finishes first, restart it by
just hitting shift-enter in the box where "f(6)" is), try typing
"f(4)" in the next box.  You're allowed to give input, but the
result doesn't get computed immediately.  You can enter several more
lines as well, etc.; when the f(6) finally finishes, SAGE goes on
to compute "f(4)".   You can queue up dozens of calculations.  For
example, if you hit the "Evaluate" link in the upper right, the
whole worksheet is queued up for computation.  Try it.   When the
computation gets stuck on "f(6)", hit the interrupt button (or press escape)
and the queued up calculations are cancelled.

Click "Hide Output" in the upper right.   You'll see just your
input and some little boxes; clicking on the boxes reveals output.

You can also embed nicely typeset math.  Try this:
\begin{verbatim}
f = maxima('sin(x^2)')
g = f.integrate('x')
view(g)
\end{verbatim}

If this silently fails, type "view(g, debug=True)" instead.
You need latex and the "convert" and "gs" commands, (use
an "apt-get install imagemagick gs").  Anyways, you get
a nicely typeset formula.  Try a matrix next:
\begin{verbatim}
A = MatrixSpace(QQ, 5).random_element()
view(A)
\end{verbatim}
Try typing this into a new box:
\begin{verbatim}
%latex
Consider the matrix $$A = \sage{A},$$
which has square $$A^2 = \sage{A^2}.$$
\end{verbatim}
If you would like to typeset a slide (suitable for presentation),
use \%slide instead.
Here is another example:
\begin{verbatim}
%latex
The first ten squares are
$$
\sage{', '.join([str(sq(i)) for i in range(1,11)])}
$$

The primes up to 100 are
$$
\sage{', '.join(str(p) for p in prime_range(100))}
$$
\end{verbatim}

\subsubsection{Using Gap, Magma, GP/PARI}
Make the first line of the input block \code{\%gap}
\code{\%magma}, or \code{\%gp}, etc.  The rest of the block
is fed directly to the corresponding interpreter.
In this way you can make a single session that has input blocks
that work with a range of different systems.

(Note -- there is currently no support for
pulling in objects and evaluating code in SAGE by typing
"sage(...)" inside the input block.  This is planned.)

\subsubsection{Typesetting}
If you have latex, gv, and the imagemagick programs (e.g., convert)
installed on your system, you can do nice latex typesetting from
within SAGE.
\begin{enumerate}
\item As usual the command \code{latex(obj)} outputs latex code
to typeset obj.
\item The command \code{view(obj)} creates an image representing
the object, which you can copy and paste into other documents.
\item If you preface a block with \code{\%latex} the rest of the
block is typeset and the corresponding image appears.
The input is also (mostly) hidden.  Use {\%latex_debug} to debug
latex problems.
\item If you preface a block with \code{\%slide} the rest of the
block is typeset as a slide (bigger san serif font)
and the corresponding image appears.  The input is again hidden.
Use {\%slide_debug} for debugging.
\end{enumerate}

Make the first line of the input block \code{\%gap}
\code{\%magma}, or \code{\%gp}, etc.  The rest of the block
is fed directly to the corresponding interpreter.
In this way you can make a single session that has input blocks
that work with a range of different systems.   You can also
pull in objects and evaluate code in SAGE by typing
"sage(...)" inside the input block.


\subsubsection{Adding and Removing Cells}
To add a new cell, click on a little black line that appears when you
hover between any two cells, or above the top one.  To delete a cell
delete all its contents, then hit backspace one more time.  The cell
vanishes forever.

You can also move back and forth between cells using the up and down
arrow.  In particular, when you are at the top of a cell and press
the up arrow the cursor jumps to the previous cell.
Press control-enter in a cell to create a new cell after the
current cell.

There is no direct support for moving and reorganizing cells, though
you can copy and paste any individual cell into another one.  However,
the "Text1" and "Text2" buttons provide the full text of the
worksheet in a very convenient format for copy and paste.


\subsubsection{History}
Click the history button near the top to pop up a history of the last
1000 (or so) input cells.  After a huge amount of design discussion about
how to design a history system, a simple popup with the text of
previous commands seems like the best choice.  It's incredibly simple,
yet provides an incredible amount of functionality, especially because
that popup window can be easily searched (at least in Firefox), pasted
from, etc., and refreshed (use F5 or Ctrl-R).


\subsubsection{Introspection}
To find all completions for an identifier you are typing press
the tab key.  This should work exactly like IPython, and even
respects the \code{trait_names()} method.

To find help for any object in a line, put ? after it
and press the tab key.  The cursor must be somewhere in the identifier
with the question mark after it.   For source code, put ?? after
the identifier and press tab.  You can also put an identifier by
itself on a line with ? (or ??) after it and press shift-enter.

To get extensive help on an object, type "help(object)" and press
return.  This works, since I set the PAGER to "cat", and I strip out
control codes that appear in the output.  And this isn't annoying,
since web browsers are very good for scrolling through long output.


\subsubsection{Objects}
When you start a notebook you give a name argument
to it, and it creates a directory.  Inside that directory there
will be many worksheets (which you can use all at once and easily
flip through -- not implemented yet), and an object store.
You can save and load objects (using save and load), and they'll
be listed in the box on the bottom let, e.g., try

a = 5
save a

and you'll see the "a" appear there.   You can load and save objects
from any worksheet in any other one.  (Currently the only way to delete
objects from the list of saved objects is to remove the object from
the objects subdirectory.)

\subsubsection{Pasting in Examples}
Code is evaluated by exec'ing (after preparsing). Only the output
of the last line of the cell is implicitly printed. If any line
starts with "sage:" or ">>>" the {\em entire block} is assumed to
contain text and examples, and only lines that begin with a
prompt are executed. Thus you can paste in *complete examples*
from the docs without any editing, and you can write input
cells that contains non-evaluated plain text mixed with
examples by starting the block with ">>>" or including an example.
(NOTE: Lines beginning with ">>>" are still preparsed.)

\subsubsection{Saving and Loading}

The SAGE notebook is very persistent.  Every time you submit
a cell for computation, the state of the notebook is saved (a
few kb's file).  If you quit the notebook and reload, it will
have everything you typed from the previous session, along
with all output.
Firefox has an excellent undo function for text input cells.
Just hit control-z to have ``infinite undo'' for the input
you've entered in that particular cell.

You can save all variables in a current session by typing
\code{save_session [optional_name]}.  You can then load
those session variables into another worksheet using
\code{load_session}, or load into the same worksheet next
time you use it.

\subsubsection{Architecture}

The SAGE Notebook is an ``AJAX application'' that can run either
entirely locally on your desktop machine, or partly on
a server and via a web browser that could be located somewhere
else.
If you run the server and allow remote access (by setting
address when starting the notebook), you should also set
the username and password, so not just anybody can access
the notebook.

Anywhere, here are the components of the SAGE Notebook:

\begin{enumerate}
\item Web Server: A Python process that uses the
      Python standard library's
     BaseHTTPServer.HTTPServer to create a web server.  This
     process also handles all requests from the web browser,
     e.g., organizing computation of cells, etc.  It
     only imports a small
     subset of the SAGE library.  In particular, if you do
     "sage -notebook" at the command line, only some of
     SAGE is imported.

 \item SAGE Server:
     A Python process with all the SAGE libraries loaded; this
     is started by (1) when a web browser first requests that
     a cell be evaluated.  There's (up to) one of these
     for each worksheet.

 \item WEB Browser: The web browser runs a 1000-line javascript (plus
     800 lines of css) program that Alex, Tom and I wrote from
     scratch, which implements much of the browser-side part of the
     SAGE notebook functionality.

\end{enumerate}

When you use the SAGE Notebook, you are mainly interacting with a
javascript program.  When you do something serious, e.g., request
computation of some input, create a new cell, etc., a request is made
from your web browser to the web server telling it what is going on.
If it's a calculation, the web server tells the SAGE server to get
started on the calculation, and tells the web browser to check several
times a second whether there is anything new with the calculation.
When something new appears it fills that in.  This continues until all
calculations are done. During this time, you can edit cells, create
new cells, submit more computations, etc.  Note that output is
updated as the computation proceeds, so you can verbosely watch
a computation progressq.  For example, try the following from the SAGE
Notebook:

\begin{verbatim}
import time
for i in range(10):
    print i
    time.sleep(0.5)
\end{verbatim}

You get to watch as the integers from 1 to 10 are "computed".
Actually, getting this output to be reported as the computation
proceeds is, I think, \emph{crucial} to making a really usable SAGE
GUI--users (i.e., me) want to run huge computations and watch the
output progress.

The architecture is also good from the point of view of being able to
interrupt running computations.  What happens when you request an
interrupt is that the web browser sends a message to the web server,
which in turn tells the SAGE server to stop computing by sending it
many interrupt signals (for several seconds) until it either stops, or
if it's really frozen (due to a bug, or calling into a C function that
isn't properly wrapped in signal handling, or maybe you run an
interactive program, e.g., via "os.system('...')"), it'll just kill that SAGE server
and start a new one.  The result is that the
user doesn't get a frozen web browser or browser interface at any point,
and even if the whole SAGE process went down and froze, at least all
your input and output from your session is still there in your
browser.  The only thing you've lost is the definition of all your
variables.  Hit "shift-enter" a few times or "evaluate all" and you're
back in shape.  This is much better than having to restart the command
prompt (e.g., with a terminal interface), then paste back in all your
setup code, etc., Also, you can save variables as you go easily (via
the "save" command), and get back to where you were quickly.

"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import os
import shutil
import socket

# SAGE libraries
from   sage.ext.sage_object import SageObject, load
from   sage.misc.viewer     import browser
from   sage.misc.misc       import alarm

# SAGE Notebook
import css          # style
import js           # javascript
import server       # web server
import worksheet    # individual worksheets (which make up a notebook)
import config       # internal configuration stuff (currently, just keycodes)
import keyboards    # keyboard layouts

MAX_WORKSHEETS = 4096  # do not change this willy nilly; that would break existing notebooks (and there is no reason to).
MAX_HISTORY_LENGTH = 500
WRAP_NCOLS = 100

# Temporarily disabled while we try fix the firefox windows hang bug.
JSMATH=False

class Notebook(SageObject):
    def __init__(self, dir='sage_notebook',
                 username=None, password=None,
                 color='default', system=None):
        self.__dir = dir
        self.set_system(system)
        self.__color = color
        if not (username is None):
            self.set_auth(username,password)
        self.__worksheets = {}
        self.__load_defaults()
        self.__filename     = '%s/nb.sobj'%dir
        self.__worksheet_dir = '%s/worksheets'%dir
        self.__object_dir   = '%s/objects'%dir
        self.__makedirs()
        self.__next_worksheet_id = 0
        self.__history = []
        W = self.create_new_worksheet('_scratch_')
        self.__default_worksheet = W
        self.save()

    def system(self):
        try:
            return self.__system
        except AttributeError:
            self.__system = None
            return None

    def set_system(self, system):
        if system == 'sage':
            self.__system = None
        elif system:  # don't change if it is None
            self.__system = system

    def color(self):
        try:
            return self.__color
        except AttributeError:
            self.__color = 'default'
            return self.__color

    def set_color(self,color):
        self.__color = color

    def set_directory(self, dir):
        if dir == self.__dir:
            return
        self.__dir = dir
        self.__filename = '%s/nb.sobj'%dir
        self.__worksheet_dir = '%s/worksheets'%dir
        self.__object_dir = '%s/objects'%dir
        for W in self.__worksheets.itervalues():
            W.set_notebook(self)

    def add_to_history(self, input_text):
        H = self.history()
        H.append(input_text)
        while len(H) > self.max_history_length():
            del H[0]

    def history(self):
        try:
            s = self.__history
        except AttributeError:
            self.__history = []
            s = self.__history
        return s

    def history_text(self):
        return '\n\n'.join([H.strip() for H in self.history()])

    def history_with_start(self, start):
        n = len(start)
        return [x for x in self.history() if x[:n] == start]

    def export_worksheet(self, worksheet_filename, filename):
        W = self.get_worksheet_with_filename(worksheet_filename)
        W.save()
        cmd = 'cd %s && tar -jcf %s.sws "%s" && mv %s.sws ..'%(
            self.__worksheet_dir,
            filename, W.filename(), filename)
        print cmd
        os.system(cmd)

    def tmpdir(self):
        d = '%s/tmp'%self.__dir
        if os.path.exists(d):
            os.system('rm -rf "%s"'%d)
        if not os.path.exists(d):
            os.makedirs(d)
        return d

    def import_worksheet(self, filename):
        if not os.path.exists(filename):
            raise ValueError, "no file %s"%filename
        if filename[-4:] != '.sws':
            raise ValueError, "file %s must have extension sws."%filename
        tmp = self.tmpdir()
        cmd = 'cd %s; tar -jxf %s'%(tmp, os.path.abspath(filename))
        print cmd
        os.system(cmd)
        D = os.listdir(tmp)[0]
        print D
        worksheet = load('%s/%s/%s.sobj'%(tmp,D,D))
        S = self.__worksheet_dir
        cmd = 'rm -rf "%s/%s"'%(S,D)
        print cmd
        os.system(cmd)
        cmd = 'mv %s/%s %s/'%(tmp, D, S)
        print cmd
        os.system(cmd)
        new_id = None
        id = worksheet.id()
        for W in self.__worksheets.itervalues():
            if W.id() == id:
                new_id = self.__next_worksheet_id
                self.__next_worksheet_id += 1
                break
        worksheet.set_notebook(self, new_id)
        name = worksheet.name()
        self.__worksheets[name] = worksheet
        return worksheet

    # unpickled, no worksheets will think they are
    # being computed, since they clearly aren't (since
    # the server just started).
    def set_not_computing(self):
        for W in self.__worksheets.values():
            W.set_not_computing()

    def default_worksheet(self):
        return self.__default_worksheet

    def directory(self):
        return self.__dir

    def DIR(self):
        """
        Return the absolute path to the directory that contains
        the SAGE Notebook directory.
        """
        return os.path.abspath('%s/..'%self.__dir)

    def max_history_length(self):
        try:
            return self.__defaults['max_history_length']
        except KeyError:
            return MAX_HISTORY_LENGTH

    def __load_defaults(self):
        # in future this will allow override by a file, and
        # can be set by user via web interface
        self.__defaults = {'cell_input_color':'#0000000',
                           'cell_output_color':'#0000EE',
                           'word_wrap_cols':WRAP_NCOLS,
                           'max_history_length':MAX_HISTORY_LENGTH}

    def worksheet_directory(self):
        return self.__worksheet_dir

    def object_directory(self):
        return self.__object_dir

    def objects(self):
        L = [x[:-5] for x in os.listdir(self.__object_dir)]
        L.sort()
        return L

    def object_list_html(self):
        m = max([len(x) for x in self.objects()] + [30])
        s = []
        a = '<a href="/%s.sobj" class="object_name">\n'
        for name in self.objects():
            s.append(a%name + name + '&nbsp;'*(m-len(name)) + '</a>\n')
        return '<br>\n'.join(s)

    def defaults(self):
        return self.__defaults

    def authorize(self, auth):
        """
        Returns True if auth is the correct authorization.
        """
        a = self.auth_string()
        if a == ':':
            return True
        return a == auth

    def auth_string(self):
        try:
            return self.__auth
        except AttributeError:
            self.__auth = ":"
        return self.__auth

    def set_auth(self, username, password):
        self.__auth = '%s:%s'%(username, password)

    def __makedirs(self):
        os.makedirs(self.__dir)
        os.makedirs(self.__worksheet_dir)
        os.makedirs(self.__object_dir)

    def worksheet_ids(self):
        return set([W.id() for W in self.__worksheets.itervalues()])

    def create_new_worksheet(self, name='untitled'):
        if name in self.__worksheets.keys():
            raise KeyError, 'name (=%s) already taken.'%name
        name = str(name)
        wids = self.worksheet_ids()
        id = 0
        while id in wids:
            id += 1

        if id >= MAX_WORKSHEETS:
            raise ValueError, 'there can be at most %s worksheets'%MAX_WORKSHEETS
        self.__next_worksheet_id += 1
        W = worksheet.Worksheet(name, self, id, system=self.system())
        self.__worksheets[name] = W
        return W

    def delete_worksheet(self, name):
        if not (name in self.__worksheets.keys()):
            raise KeyError, "Attempt to delete missing worksheet"
        W = self.__worksheets[name]
        cmd = 'rm -rf "%s"'%(W.directory())
        print cmd
        os.system(cmd)

        del self.__worksheets[name]
        if len(self.__worksheets) == 0:
            return self.create_new_worksheet('_scratch_')
        else:
            return self.__worksheets[self.__worksheets.keys()[0]]

    def worksheet_names(self):
        W = self.__worksheets.keys()
        W.sort()
        return W

    def get_worksheet_with_name(self, name):
        return self.__worksheets[name]

    def get_worksheet_with_id(self, id):
        if id != None:
            for W in self.__worksheets.itervalues():
                if W.id() == id:
                    return W
        return self.__worksheets[self.__worksheets.keys()[0]]

    def get_worksheet_with_filename(self, filename):
        if id != None:
            for W in self.__worksheets.itervalues():
                if W.filename() == filename:
                    return W
        raise KeyError, "no such worksheet %s"%filename

    def get_worksheet_that_has_cell_with_id(self, id):
        worksheet_id = id // MAX_WORKSHEETS
        return self.get_worksheet_with_id(worksheet_id)

    def save(self, filename=None):
        if filename is None:
            F = os.path.abspath(self.__filename)
            try:
                shutil.copy(F, F[:-5] + '-backup.sobj')
            except IOError:
                pass
            SageObject.save(self, os.path.abspath(self.__filename))
        else:
            SageObject.save(self, os.path.abspath(filename))

    def start(self, port=8000, address='localhost',
                    max_tries=128, open_viewer=False,
                    jsmath=True):
        global JSMATH
        JSMATH = jsmath
        tries = 0
        while True:
            try:
                notebook_server = server.NotebookServer(self,
                         port, address)
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

        s = "Open your web browser to http://%s:%s"%(address, port)
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
            cmd = '%s http://%s:%s 1>&2 >/dev/null &'%(browser(), address, port)
            os.system(cmd)
        notebook_server.serve()
        self.save()
        self.quit()

    def quit(self):
        for W in self.__worksheets.itervalues():
            W.quit()

    def worksheet_list_html(self, current_worksheet=None):
        s = []
        names = self.worksheet_names()
        m = max([len(x) for x in names] + [30])
        for n in names:
            W = self.__worksheets[n]
            if W == current_worksheet:
                cls = 'worksheet_current'
            else:
                cls = 'worksheet_other'
            if W.computing():
                cls += '_computing' # actively computing
            name = W.name()
            name += ' (%s)'%len(W)
            name += ' '*(m-len(name))
            name = name.replace(' ','&nbsp;')
            txt = '<a class="%s" onClick="switch_to_worksheet(%s)" onMouseOver="show_worksheet_menu(%s)" target="_new" href="/%s">%s</a>'%(
                cls,W.id(),W.id(),W.id(),name)
            s.append(txt)
        return '<br>'.join(s)

    def _html_head(self, worksheet_id):
        worksheet = self.get_worksheet_with_id(worksheet_id)
        head = '<title>%s (%s)</title>'%(worksheet.name(), self.directory())
        head += '<style>' + css.css(self.color()) + '</style>\n'

        if JSMATH:
            head += '<script>jsMath = {Controls: {cookie: {scale: 125}}}</script>\n'
            #head += '<script src="/jsmath/plugins/spriteImageFonts.js"></script>\n'
            head +=' <script src="/jsmath/plugins/noImageFonts.js"></script>\n'
            head += '<script src="/jsmath/jsMath.js"></script>\n'
            head += "<script>jsMath.styles['#jsMath_button'] = jsMath.styles['#jsMath_button'].replace('right','left');</script>\n"
        head += '<script language=javascript>' + js.javascript() + '</script>\n'

        return head

    def _html_body(self, worksheet_id):
        worksheet = self.get_worksheet_with_id(worksheet_id)
        #if worksheet.computing():
        interrupt_class = "interrupt"
        #else:
        #    interrupt_class = "interrupt_grey"

        add_new_worksheet_menu = """
             <div class="add_new_worksheet_menu" id="add_worksheet_menu">
             <input id="new_worksheet_box" class="add_new_worksheet_menu"
                    onKeyPress="if(is_submit(event)) process_new_worksheet_menu_submit();"></input>
             <button class="add_new_worksheet_menu"  onClick="process_new_worksheet_menu_submit();">add</button>
             &nbsp;&nbsp;&nbsp;<span class="X" onClick="hide_add_new_worksheet_menu()">X</span>
             </div>
        """

        delete_worksheet_menu = """
             <div class="delete_worksheet_menu" id="delete_worksheet_menu">
             <input id="delete_worksheet_box" class="delete_worksheet_menu"
                    onKeyPress="if(is_submit(event)) process_delete_worksheet_menu_submit();"></input>
             <button class="delete_worksheet_menu" onClick="process_delete_worksheet_menu_submit();">delete</button>
             &nbsp;&nbsp;&nbsp;<span class="X" onClick="hide_delete_worksheet_menu()">X</span>
             </div>
        """

        vbar = '<span class="vbar"></span>'

        body = ''
        body += '<div class="top_control_bar">\n'
        body += '  <span class="banner"><a class="banner" href="http://modular.math.washington.edu/sage">SAGE</a> %s</span>\n'%self.__dir
        body += '  <span class="control_commands">\n'

        body += '    <a class="help" onClick="show_help_window()">Help</a>' + vbar
        body += '    <a class="history_link" onClick="history_window()">History</a> ' + vbar
        body += '    <a class="plain_text" onClick="worksheet_text_window(\'%s\')">Text</a>'%worksheet.filename() + vbar
        body += '    <a class="doctest_text" onClick="doctest_window(\'%s\')">Text2</a>'%worksheet.filename() + vbar
        body += '    <a class="doctest_text" onClick="print_window(\'%s\')">Print</a>'%worksheet.filename() + vbar
        body += '    <a class="evaluate" onClick="evaluate_all()">Evaluate</a>' + vbar
        body += '    <a class="hide" onClick="hide_all()">Hide</a>' + vbar
        body += '    <a class="hide" onClick="show_all()">Show</a>' + vbar
        body += '     <a onClick="show_upload_worksheet_menu()" class="upload_worksheet">Open</a>' + vbar
        body += '    <a class="download_sws" href="%s.sws">Save</a>'%worksheet.filename() + vbar
        body += '    <a class="%s" onClick="interrupt()" id="interrupt">Interrupt</a>'%interrupt_class + vbar
        body += '    <a class="restart_sage" onClick="restart_sage()" id="restart_sage">Restart</a>'
        body += '  </span>'
        body += '</div>'
        body += '\n<div class="worksheet" id="worksheet">\n' + worksheet.html() + '\n</div>\n'

        body += '<span class="pane"><table bgcolor="white"><tr><td>\n'
        body += '  <div class="worksheets_topbar">'
        body += '     <a onClick="show_add_new_worksheet_menu()" class="new_worksheet">New</a> '
        body += '     <a onClick="show_delete_worksheet_menu()" class="delete_worksheet">Delete</a> '
        body += '  &nbsp;Worksheets</div>\n'
        body +=    add_new_worksheet_menu
        body +=    delete_worksheet_menu
        body += '  <div class="worksheet_list" id="worksheet_list">%s</div>\n'%self.worksheet_list_html(worksheet)
        body += '  <div class="objects_topbar">Saved Objects</div>\n'
        body += '  <div class="object_list" id="object_list">%s</div>\n'%self.object_list_html()
        body += '<br>\n'
        body += '  <div class="variables_topbar">Variables</div>\n'
        body += '  <div class="variables_list" id="variable_list">%s</div>\n'%\
                worksheet.variables_html()
        body += '  <div class="attached_topbar">Attached Files</div>\n'
        body += '  <div class="attached_list" id="attached_list">%s</div><br>\n'%\
                worksheet.attached_html()
        body += '</td></tr></table></span>\n'
        body += '<script language=javascript>focus(%s)</script>\n'%(worksheet[0].id())
        body += '<script language=javascript>jsmath_init();</script>\n'

        if worksheet.computing():
            # Set the update checking back in motion.
            body += '<script language=javascript> check_for_cell_output() </script>\n'

        return body

    def help_window(self):
        help = [
            ('HTML', 'Begin an input block with %html and it will be output as HTML.  Use the &lt;sage>...&lt;/sage> tag to do computations in an HTML block and have the typeset output inserted.  Use &lt;$>...&lt;/$> and &lt;$$>...&lt;/$$> to insert typeset math in the HTML block.  This does <i>not</i> require latex.'),
            ('shell', 'Begin a block with %sh to have the rest of the block evaluated as a shell script.  The current working directory is maintained.'),
               ('Evaluate Input', 'Press shift-enter.  You can start several calculations at once.  If you press control-enter instead, then a new cell is created after the current one.'),
                ('Timed Evaluation', 'Type "time" at the beginning of the cell.'),
                ('Evaluate all cells', 'Click <u>Evaluate All</u> in the upper right.'),
                ('Evaluate cell using <b>GAP, Singular, etc.', 'Put "%gap", "%singular", etc. as the first input line of a cell; the rest of the cell is evaluated in that system.'),
                ('Typeset a cell', 'Make the first line of the cell "%latex". The rest of the cell should be the body of a latex document.  Use \\sage{expr} to access SAGE from within the latex.  Evaluated typeset cells hide their input.  Use "%latex_debug" for a debugging version.  You must have latex for this to work.'),
               ('Typeset a slide', 'Same as typesetting a cell but use "%slide" and "%slide_debug"; will use a large san serif font.  You must have latex for this to work.'),
                ('Typesetting', 'Type "latex(objname)" for latex that you can paste into your paper.  Type "view(objname)" or "show(objname)", which will display a nicely typeset image (using javascript!).  You do <i>not</i> need latex for this to work.  Type "lprint()" to make it so output is often typeset by default.'),
                ('Move between cells', 'Use the up and down arrows on your keyboard.'),
                ('Interrupt running calculations',
                 'Click <u>Interrupt</u> in the upper right or press escape in any input cell. This will (attempt) to interrupt SAGE by sending many interrupts for several seconds; if this fails, it restarts SAGE (your worksheet is unchanged, but your session is reset).'),
                ('Tab completion', 'Press tab while the cursor is on an identifier.'),
                ('Print worksheet', 'Click the print button.'),
                ('Help About',
                 'Type ? immediately after the object or function and press tab.'),
                ('Source Code',
                 'Put ?? after the object and press tab.'),
                ('Hide Input',
                 'Put %hide at the beginning of the cell.  This can be followed by %gap, %latex, %maxima, etc.  Note that %hide must be first. Put a blank line at the beginning so the "%hide" does not appear.'),
                ('Detailed Help',
                 'Type "help(object)" and press shift-return.'),
                ('Insert New Cell',
                 'Put mouse between an output and input until the horizontal line appears and click.  Also if you press control-enter in a cell, a new cell is inserted after it.'),
                ('Delete Cell',
                 'Delete cell contents the press backspace.'),
                ('Text of Worksheet', 'Click the <u>Text</u> and <u>Doctext</u> links, which are very useful if you need to cut and paste chunks of your session into email or documentation.'),
                ('History', 'Click the <u>History</u> link for a history of commands entered in any worksheet of this notebook.  This appears in a popup window, which you can search (control-F) and copy and paste from.'),
                ('Hide/Show Output', 'Click on the left side of output to toggle between hidden, shown with word wrap, and shown without word wrap.'),
                ('Hide/Show All Output', 'Click <u>Hide</u> in the upper right to hide <i>all</i> output. Click <u>Show</u> to show all output.'),
                ('Variables',
                 'All variables and functions that you create during this session are listed on the left.  Even predefined variables that you overwrite will appear.'),
                ('Objects',
                 'All objects that you save in <i>any worksheet</i> are listed on the left.  Use "save(obj, name)" and "obj = load(name)" to save and load objects.'),
                ('Loading and Saving Sessions', 'Use "save_session name" to save all variables to an object with given name (if no name is given, defaults to name of worksheet).  Use "load_session name" to <i>merge</i> in all variables from a saved session.'),
                ('Loading and Saving Objects', 'Use "save obj1 obj2 ..." and "load obj1 obj2 ...".  This allows very easy moving of objects from one worksheet to another, and saving of objects for later use.'),
                ('Loading SAGE/Python Scripts', 'Use "load filename.sage" and "load filename.py".  Load is relative to the path you started the notebook in.  The .sage files are preparsed and .py files are not.   You may omit the .sage or .py extension.  Files may load other files.'),
                ('Attaching Scripts', 'Use "attach filename.sage" or "attach filename.py".  Attached files are automatically reloaded when the file changes.  The file $HOME/.sage/init.sage is attached on startup if it exists.'),
                ('Saving Worksheets',
                 'Click <ul>Save</ul> in the upper right to download a complete worksheet to a local .sws file.  Note that uploading is not implemented yet except locally (the <i>only</i> thing left is implementing upload of binary data to the server). Note that <i>everything</i> that has been submitted is automatically saved to disk, and is there for you next time you access the notebook.'),
                ('Restart', 'Type "restart" to restart the SAGE interpreter for a given worksheet.  (You have to interrupt first.)'),
                ('Input Rules', "Code is evaluated by exec'ing (after preparsing).  Only the output of the last line of the cell is implicitly printed.  If any line starts with \"sage:\" or \">>>\" the entire block is assumed to contain text and examples, so only lines that begin with a prompt are executed.   Thus you can paste in complete examples from the docs without any editing, and you can write input cells that contains non-evaluated plain text mixed with examples by starting the block with \">>>\" or including an example."),
                ('Working Directory', 'Each block of code is run from its own directory.  The variable DIR contains the directory from which you started the SAGE notebook.  For example, to open a file in that directory, do "open(DIR+\'filename\')".'),
                ('Customizing the look', 'Learn about cascading style sheets (CSS), then create a file notebook.css in your $HOME/.sage directory.  Use "view source" on a notebook web page to see the CSS that you can override.'),
                ('Emacs Keybindings', 'If you are using GNU/Linux, you can change (or create) a <tt>.gtkrc-2.0</tt> file.  Add the line <tt>gtk-key-theme-name = "Emacs"</tt> to it.  See <a target="_blank" href="http://kb.mozillazine.org/Emacs_Keybindings_(Firefox)">this page</a> [mozillazine.org] for more details.'),
                ('More Help', 'Type "help(sage.server.notebook.notebook)" for a detailed discussion of the architecture of the SAGE notebook and a tutorial (or see the SAGE reference manual).'),
                ]

        help.sort()
        s = """

        This is the SAGE Notebook, which is the graphical interface to
        the computer algebra system SAGE (Software for Algebra and
        Geometry Exploration).
        <br><b>Use Firefox:</b> <i>It currently only works in <b>Firefox</b>, but might
        work to some extent in other browsers.</i><br>
        AUTHORS: William Stein, Tom Boothby, Alex Clemesha (with feedback from many people,
        especially Fernando Perez and Joe Wetherell).

        <br><hr>
        <style>
        div.help_window {
            background-color:white;
            border: 3px solid #3d86d0;
            top: 10ex;
            bottom:10%;
            left:25%;
            right:15%;
            padding:2ex;
        }


        table.help_window {
            background-color:white;
            width:100%;
        }

        td.help_window_cmd {
            background-color: #f5e0aa;
            width:30%;
            padding:1ex;
            font-weight:bold;
        }

        td.help_window_how {
            padding:1ex;
            width:70%;
        }
        </style>
        <div class="help_window">

        A <i>worksheet</i> is an ordered list of SAGE calculations with output.
        A <i>session</i> is a worksheet and a set of variables in some state.
        A <i>notebook</i> is a collection of worksheets and saved objects.

        <table class="help_window">
        """
        for x, y in help:
            s += '<tr><td class="help_window_cmd">%s</td><td class="help_window_how">%s</td></tr>'%(x,y)
        s += '</table></div>'
        return s

    def upload_window(self):
        return """
          <html>
            <head>
              <title>Upload File</title>
              <style>%s</style>
              <script language=javascript>%s</script>
            </head>
            <body onLoad="if(window.focus) window.focus()">
              <div class="upload_worksheet_menu" id="upload_worksheet_menu">
              <form method="POST" action="upload_worksheet" target="_new"
                    name="upload" enctype="multipart/form-data">
              <input class="upload_worksheet_menu" type="file" name="fileField" id="upload_worksheet_filename"></input><br>
              <input type="button" class="upload_worksheet_menu" value="upload" onClick="form.submit(); window.close();">
              </form><br>
              </div>
            </body>
          </html>
         """%(css.css(self.color()),js.javascript())

    def html(self, worksheet_id=None, authorized=False):
        if worksheet_id is None:
            W = self.default_worksheet()
            worksheet_id = W.id()
        else:
            try:
                W = self.get_worksheet_with_id(worksheet_id)
            except KeyError:
                W = self.default_worksheet()
                worksheet_id = W.id()

        if authorized:
            body = self._html_body(worksheet_id)
        else:
            body = self._html_authorize()

        body += '<script language=javascript>worksheet_id=%s; worksheet_filename="%s"</script>'%(worksheet_id, W.filename())

        head = self._html_head(worksheet_id)
        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def _html_authorize(self):
        return """
        <h1>SAGE Notebook Server</h1>
        <div id="mainbody" class="login">Sign in to the SAGE Notebook<br>
        <form>
        <table>
        <tr><td>
          <span class="username">Username:</span></td>
          <td><input name="username" class="username"
                      onKeyPress="if(is_submit(event)) login(username.value, password.value)"></td>
        </tr>
        <tr><td>
           <span class="password">Password:</span></td>
           <td><input name="password" class="username" type="password"
                      onKeyPress="if(is_submit(event)) login(username.value, password.value)"></td>
        </tr>
        <td>&nbsp</td>
        <td>
           <input type='button' onClick="login(username.value,password.value);" value="Sign in">
           </td></table>
                   </form></div>

        """

    def format_completions_as_html(self, cell_id, completions):
        if len(completions) == 0:
            return ''
        lists = []

        # compute the width of each column
        column_width = []
        for i in range(len(completions[0])):
            column_width.append(max([len(x[i]) for x in completions if i < len(x)]))

        for i in range(len(completions)):
            row = completions[i]
            for j in range(len(row)):
                if len(lists) <= j:
                    lists.append([])
                cell = """
   <li id='completion%s_%s_%s' class='completion_menu_two'>
    <a onClick='do_replacement(%s, "%s")'
       onMouseOver='this.focus(); select_replacement(%s,%s);'
    >%s</a>
   </li>"""%(cell_id, i, j, cell_id, row[j], i,j,
             row[j] + '&nbsp;'*(column_width[j]-len(row[j])) )
                lists[j].append(cell)

        grid = "<ul class='completion_menu_one'>"
        for L in lists:
            s = "\n   ".join(L)
            grid += "\n <li class='completion_menu_one'>\n  <ul class='completion_menu_two'>\n%s\n  </ul>\n </li>"%s

        return grid + "\n</ul>"


def notebook(dir       ='sage_notebook',
             port      = 8000,
             address   = 'localhost',
             open_viewer    = False,
             max_tries = 10,
             username  = None,
             password  = None,
             color     = None,
             system    = None,
             jsmath    = True):
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
        open_viewer -- bool (default:False); if True, pop up a web browser at the URL
        max_tries -- (default: 10) maximum number of ports > port to try in
                     case given port can't be opened.
        username -- user name used for authenticated logins
        password -- password used for authenticated logins
        color -- string or pair of html colors, e.g.,
                    'gmail'
                    'grey'
                    ('#ff0000', '#0000ff')
        system -- default computer algebra system to use for new
                  worksheets.
        jsmath -- whether not to enable javascript typset output for math.

    NOTES:

    When you type \code{notebook(...)}  you start a web server on the
    machine you type that command on.  You don't connect to another
    machine.  So do this if you want to start a SAGE notebook
    accessible from anywhere:

    \begin{enumerate}
    \item Figure out the external address of your server, say
          'sage.math.washington.edu', for example.
    \item On your server, type
       server_http1('mysession', address='sage.math.washington.edu')
    \item Assuming you have permission to open a port on that
       machine, it will startup and display a URL, e.g.,
           \url{http://sage.math.washington.edu:8000}
       Note this URL.
    \item Go to any computer in the world (!), or at least
       behind your firewall, and use any web browser to
       visit the above URL.  You're using \sage.
    \end{enumerate}

    \note{There are no security precautions in place \emph{yet}!  If
    you open a server as above, and somebody figures this out, they
    could use their web browser to connect to the same sage session,
    and type something nasty like \code{os.system('cd; rm -rf *')}
    and your home directory would be hosed.   I'll be adding an
    authentication screen in the near future.  In the meantime
    (and even then), you should consider creating a user with
    very limited privileges (e.g., empty home directory).}

    FIREFOX ISSUE:
    If your default web browser if Firefox, then notebook will
    open a copy of Firefox at the given URL.  You should
    definitely set the "open links in new tabs" option in
    Firefox, or you might loose a web page you were looking at.
    To do this, just go to

         Edit --> Preferences --> Tabs

    and in "Open links from other apps" select the middle button
    instead of the bottom button.
    """
    if os.path.exists(dir):
        if not os.path.isdir(dir):
            raise RuntimeError, '"%s" is not a valid SAGE notebook directory (it is not even a directory).'%dir
        if not (os.path.exists('%s/nb.sobj'%dir) or os.path.exists('%s/nb-backup.sobj'%dir)):
            raise RuntimeError, '"%s" is not a valid SAGE notebook directory (missing nb.sobj).'%dir
        try:
            nb = load('%s/nb.sobj'%dir)
        except:
            print "****************************************************************"
            print "  * * * WARNING   * * * WARNING   * * * WARNING   * * * "
            print "WARNING -- failed to load notebook data. Trying the backup file."
            print "****************************************************************"
            try:
                nb = load('%s/nb-backup.sobj'%dir)
            except:
                print "Recovering from last op save failed."
                print "Trying save from last startup."
                nb = load('%s/nb-older-backup.sobj'%dir)

        nb.set_directory(dir)
        if not (username is None):
            nb.set_auth(username=username, password=password)
        if not (color is None):
            nb.set_color(color)
        if not system is None:
            nb.set_system(system)
        nb.set_not_computing()
    else:
        nb = Notebook(dir,username=username,password=password, color=color, system=system)
    nb.save()
    shutil.copy('%s/nb.sobj'%dir, '%s/nb-older-backup.sobj'%dir)
    nb.start(port, address, max_tries, open_viewer, jsmath=jsmath)
    alarm(3)
    from sage.interfaces.quit import expect_quitall
    expect_quitall(verbose=False)
    from sage.misc.misc import delete_tmpfiles
    delete_tmpfiles()
    return nb






