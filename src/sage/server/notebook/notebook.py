r"""
SAGE Notebook Interface

AUTHORS:
    -- William Stein (2006-05-06): initial version
    -- Alex Clemesha
    -- Tom Boothby: support for a wide range of web browsers; refactoring of javascript code; systematic keyboard controls

The SAGE graphical user interface is unusual in that it operates via
your web browser.  It provides you with SAGE worksheets that you can
edit and evaluate, which contain scalable typeset mathematics and
beautiful antialised images.  To try it out immediately, do this:

    sage.: notebook(open_viewer=True)
    the sage notebook starts...

\subsection{Supported Browsers}

The SAGE notebook should fully work with Firefox (and Mozilla),
Safari, and Opera. The notebook works somewhat in Internet Explorer.

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

\subsubsection{Typesetting Mathematics}
SAGE \emph{includes} jsMath, which is an implementation of the TeX
math layout engine in javascript.  If you use the show or view
commands, they display a given SAGE object typeset using jsmath.
Moreover, if you put \code{\%jsmath} at the beginning of an input
cell, the whole cell will be typeset using jsmath.  Also, you can type
\code{jsmath(obj)} to typeset a given object obj using jsmath.


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
the "Text" buttons provide the full text of the worksheet in a very
convenient format for copy and paste.


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


\subsubsection{Saving and Loading Individual Objects}
When you start a notebook you give a name argument
to it, and it creates a directory.  Inside that directory there
will be many worksheets (which you can use all at once and easily
flip through -- not implemented yet), and an object store.
You can save and load individual objects (using save and load), and they'll
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

\subsubsection{Saving and Loading Notebooks and Worksheets}

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
a computation progress.  For example, try the following from the SAGE
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

## This is commented out, since it's not recommended.  I really
## don't like crap that is both potentially insecure and will
## break on some setups.
## \subsubsection{Typesetting with Latex}
## If you have latex, gv, and the imagemagick programs (e.g., convert)
## installed on your system, you can do nice latex typesetting from
## within SAGE.
## \begin{enumerate}
## \item As usual the command \code{latex(obj)} outputs latex code
## to typeset obj.
## \item The command \code{view(obj)} creates an image representing
## the object, which you can copy and paste into other documents.
## \item If you preface a block with \code{\%latex} the rest of the
## block is typeset and the corresponding image appears.
## The input is also (mostly) hidden.  Use {\%latex_debug} to debug
## latex problems.
## \item If you preface a block with \code{\%slide} the rest of the
## block is typeset as a slide (bigger san serif font)
## and the corresponding image appears.  The input is again hidden.
## Use {\%slide_debug} for debugging.
## \end{enumerate}



###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import os
import shutil
import socket
import re           # regular expressions

# SAGE libraries
from   sage.structure.sage_object import SageObject, load
from   sage.misc.viewer     import browser
from   sage.misc.misc       import alarm, cancel_alarm

# SAGE Notebook
import css          # style
import js           # javascript
import server       # web server
import worksheet    # individual worksheets (which make up a notebook)
import config       # internal configuration stuff (currently, just keycodes)
import keyboards    # keyboard layouts

MAX_WORKSHEETS = 4096  # do not change this willy nilly; that would break existing notebooks (and there is no reason to).
MAX_HISTORY_LENGTH = 500
WRAP_NCOLS = 80

# Temporarily disabled while we try fix the firefox windows hang bug.
JSMATH=False

class Notebook(SageObject):
    def __init__(self, dir='sage_notebook', username=None,
                 password=None, color='default', system=None,
                 show_debug = False, log_server=False,
                 kill_idle=False, splashpage=False):
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
        self.__history_count = 0
        self.__log_server = log_server #log all POST's and GET's
        self.__server_log = [] #server log list
        W = self.create_new_worksheet('_scratch_')
        self.__default_worksheet = W
        self.__show_debug = show_debug
        self.__kill_idle = kill_idle
        self.__splashpage = splashpage if splashpage is not None else False
        self.save()

    def kill_idle(self):
        """
        Returns the idle timeout.  0 means don't kill
        idle processes.
        """
        try:
            return self.__kill_idle
        except AttributeError:
            self.__kill_idle = 0
            return 0

    def kill_idle_every_so_often(self):
        raise NotImplementedError


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

    def splashpage(self):
        try:
            return self.__splashpage
        except AttributeError:
            self.__splashpage = True
            return self.__splashpage

    def set_splashpage(self, splashpage):
        self.__splashpage = splashpage

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

    def history_count_inc(self):
        self.__history_count += 1

    def history_count(self):
        return self.__history_count

    def server_log(self):
        return self.__server_log

    def log_server(self):
        return self.__log_server

    def set_log_server(self, log_server):
        self.__log_server = log_server

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
        try:
            D = os.listdir(tmp)[0]
        except IndexError:
            raise ValueError, "invalid worksheet"
        worksheet = load('%s/%s/%s.sobj'%(tmp,D,D), compress=False)
        names = self.worksheet_names()
        if D in names:
            m = re.match('.*?([0-9]+)$',D)
            if m is None:
                n = 0
            else:
                n = int(m.groups()[0])
            while "%s%d"%(D,n) in names:
                n += 1
            cmd = 'mv %s/%s/%s.sobj %s/%s/%s%d.sobj'%(tmp,D,D,tmp,D,D,n)
            print cmd
            os.system(cmd)
            cmd = 'mv %s/%s %s/%s%d'%(tmp,D,tmp,D,n)
            print cmd
            os.system(cmd)
            D = "%s%d"%(D,n)
            worksheet.set_name(D)
        print D
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

    def set_debug(self,show_debug):
        self.__show_debug = show_debug

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
                           'word_wrap_cols':int(WRAP_NCOLS),
                           'max_history_length':MAX_HISTORY_LENGTH}

    def worksheet_directory(self):
        return self.__worksheet_dir

    def object_directory(self):
        O = self.__object_dir
        if not os.path.exists(O):
            os.makedirs(O)
        return O

    def objects(self):
        L = [x[:-5] for x in os.listdir(self.object_directory())]
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

    def create_new_worksheet(self, name='untitled', passcode=''):
        if name in self.__worksheets.keys():
            raise KeyError, 'name (=%s) already taken.'%name
        name = str(name)
        passcode = str(passcode)
        wids = self.worksheet_ids()
        id = 0
        while id in wids:
            id += 1

        if id >= MAX_WORKSHEETS:
            raise ValueError, 'there can be at most %s worksheets'%MAX_WORKSHEETS
        self.__next_worksheet_id += 1
        W = worksheet.Worksheet(name, self, id, system=self.system(), passcode=passcode)
        self.__worksheets[name] = W
        return W

    def delete_worksheet(self, name):
        """
        Delete the given worksheet and remove its
        name from the worksheet list.
        """
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
        """
        Get the worksheet with given id, which is either a name or the id number.
        If there is no such worksheet, a KeyError is raised.

        INPUT:
            id -- something that identifies a worksheet.
                 string -- use worksheet with that name or filename.
                 None -- use the default worksheet.
                 string int -- something that coerces to an integer; worksheet with that number

        OUTPUT:
            a worksheet.
        """
        if id is None:
            return self.default_worksheet()
        try:
            id = int(id)
            for W in self.__worksheets.itervalues():
                if W.id() == id:
                    return W
        except ValueError:
            id = str(id).lower()
            for W in self.__worksheets.itervalues():
                if W.name().lower() == id or W.filename().lower() == id:
                    return W
        raise KeyError, 'no worksheet %s'%id

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
        print "-"*70

        if filename is None:
            F = os.path.abspath(self.__filename)
            try:
                shutil.copy(F, F[:-5] + '-backup.sobj')
            except IOError:
                pass
            F = os.path.abspath(self.__filename)
        else:
            F = os.path.abspath(filename)

        print "Saving notebook to '%s'..."%F
        SageObject.save(self, F, compress=False)
        print "Press control-C twice to stop notebook server."
        print "-"*70

    def start(self, port=8000, address='localhost',
                    max_tries=128, open_viewer=False,
                    jsmath=False):
        global JSMATH
        JSMATH = jsmath
        tries = 0
        port = int(port)
        max_tries = int(max_tries)
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
        print "WARNING: The SAGE Notebook works best with Firefox/Mozilla, Safari, and Opera."

        if open_viewer:
            cmd = '%s http://%s:%s 1>&2 >/dev/null &'%(browser(), address, port)
            os.system(cmd)
        while True:
            try:
                notebook_server.serve()
            except KeyboardInterrupt, msg:
                break
            except Exception, msg:
                print msg
                print "Automatically restarting server."
            else:
                break
        self.save()
        self.quit()
        self.save()

    def quit(self):
        for W in self.__worksheets.itervalues():
            W.quit()

    def worksheet_list_html(self, current_worksheet=None):
        s = []
        names = self.worksheet_names()
        m = max([len(x) for x in names] + [30])
        for n in names:
            if n == 'doc_browser': continue
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
            txt = '<a class="%s" onClick="switch_to_worksheet(\'%s\')" onMouseOver="show_worksheet_menu(%s)" href="/%s">%s</a>'%(
                #cls,W.id(),W.id(),W.id(),name)
                cls,W.id(),W.id(), W.filename(),name)
            s.append(txt)
        return '<br>'.join(s)

    def _doc_html_head(self, worksheet_id, css_href):
        if worksheet_id is not None:
            worksheet = self.get_worksheet_with_id(worksheet_id)
            head = '\n<title>%s (%s)</title>'%(worksheet.name(), self.directory())
        else:
            head = '\n<title>SAGE Notebook | Welcome</title>'
        head += '\n<script language=javascript src="/__main__.js"></script>\n'
        head += '\n<link rel=stylesheet href="/__main__.css" type="text/css" />\n'

        if css_href:
            head += '\n<link rel=stylesheet type="text/css" href=%s />\n'%(css_href)

        if JSMATH:
            head += '<script>jsMath = {Controls: {cookie: {scale: 125}}}</script>\n'
            #head += '<script src="/jsmath/plugins/spriteImageFonts.js"></script>\n'
            head +=' <script src="/jsmath/plugins/noImageFonts.js"></script>\n'
            head += '<script src="/jsmath/jsMath.js"></script>\n'
            head += "<script>jsMath.styles['#jsMath_button'] = jsMath.styles['#jsMath_button'].replace('right','left');</script>\n"
        #head += '<script language=javascript>' + js.javascript() + '</script>\n'
        return head

    def _doc_html_body(self, worksheet_id):
        worksheet = self.get_worksheet_with_id(worksheet_id)
        if worksheet.computing():
            interrupt_class = "interrupt"
        else:
            interrupt_class = "interrupt_grey"
        main_body = worksheet.html(authorized = True)

        vbar = '<span class="vbar"></span>'

        body = ''
        body += '<div class="top_control_bar">\n'
        body += '  <span class="banner"><a class="banner" href="http://sage.math.washington.edu/sage">'
        body += '  <img src="sagelogo.png"/></a></span>\n'
        body += '  <span class="control_commands" id="cell_controls">\n'
        body += '    <a class="%s" onClick="interrupt()" id="interrupt">Interrupt</a>'%interrupt_class + vbar
        body += '    <a class="restart_sage" onClick="restart_sage()" id="restart_sage">Restart</a>' + vbar
        body += '    <a class="history_link" onClick="history_window()">History</a>' + vbar
        body += '    <a class="help" onClick="show_help_window()">Help</a>' + vbar
        # body += '    <a class="slide_mode" onClick="slide_mode()">Slideshow</a>'
        body += '  </span>\n'

        #these divs appear in backwards order because they're float:right
        body += '  <div class="hidden" id="slide_controls">\n'
        body += '    <div class="slideshow_control">'
        body += '      <a class="slide_arrow" onClick="slide_next()">&gt;</a>'
        body += '      <a class="slide_arrow" onClick="slide_last()">&gt;&gt;</a>' + vbar
        body += '      <a class="cell_mode" onClick="cell_mode()">Worksheet</a>'
        body += '    </div>'
        body += '    <div class="slideshow_progress" id="slideshow_progress" onClick="slide_next()">'
        body += '      <div class="slideshow_progress_bar" id="slideshow_progress_bar">&nbsp;</div>'
        body += '      <div class="slideshow_progress_text" id="slideshow_progress_text">&nbsp;</div>'
        body += '    </div>'
        body += '    <div class="slideshow_control">'
        body += '      <a class="slide_arrow" onClick="slide_first()">&lt;&lt;</a>'
        body += '      <a class="slide_arrow" onClick="slide_prev()">&lt;</a>'
        body += '    </div>'
        body += '  </span>\n'

        body += '</div>'
        body += '\n<div class="slideshow" id="worksheet">\n'

        body += main_body + '\n</div>\n'

        # The blank space given by '<br>'*15  is needed so the input doesn't get
        # stuck at the bottom of the screen. This could be replaced by a region
        # such that clicking on it creates a new cell at the bottom of the worksheet.
        body += '<br>'*15
        body += '\n</div>\n'

        body += '<span class="pane" id="left_pane"><table bgcolor="white"><tr><td>\n'
        body += '</td></tr></table></span>\n'

        body += '  <div class="worksheet_list" id="worksheet_list">%s</div>\n'%self.worksheet_list_html(worksheet)
        body += '<script language=javascript>focus(%s)</script>\n'%(worksheet[0].id())
        body += '<script language=javascript>jsmath_init();</script>\n'
        body += '<script language=javascript>worksheet_locked=false;</script>'

        if worksheet.computing():
            # Set the update checking back in motion.
            body += '<script language=javascript> active_cell_list = %r; \n'%worksheet.queue_id_list()
            body += 'for(var i = 0; i < active_cell_list.length; i++)'
            body += '    cell_set_running(active_cell_list[i]); \n'
            body += 'start_update_check(); </script>\n'
        return body

    def _html_head(self, worksheet_id):
        if worksheet_id is not None:
            worksheet = self.get_worksheet_with_id(worksheet_id)
            head = '\n<title>%s (%s)</title>'%(worksheet.name(), self.directory())
        else:
            head = '\n<title>SAGE Notebook | Welcome</title>'
        head += '\n<script language=javascript src="/__main__.js"></script>\n'
        head += '\n<link rel=stylesheet href="/__main__.css" type="text/css" />\n'

        if JSMATH:
            head += '<script>jsMath = {Controls: {cookie: {scale: 115}}}</script>\n'
            head +=' <script src="/jsmath/plugins/noImageFonts.js"></script>\n'
            head += '<script src="/jsmath/jsMath.js"></script>\n'
            head += "<script>jsMath.styles['#jsMath_button'] = jsMath.styles['#jsMath_button'].replace('right','left');</script>\n"
        #head += '<script language=javascript>' + js.javascript() + '</script>\n'

        return head

    def _html_body(self, worksheet_id, show_debug=False, worksheet_authorized=False):
        if worksheet_id is None or worksheet_id == '':
            main_body = '<div class="worksheet_title">Welcome to the SAGE Notebook</div>\n'
            if os.path.isfile(self.directory() + "/index.html"):
                splash_file = open(self.directory() + "/index.html")
                main_body+= splash_file.read()
                splash_file.close()
            else:
                dir = os.path.abspath('%s'%self.directory())
                main_body+= "<br>&nbsp;&nbsp;&nbsp;SAGE Notebook running from <tt>%s</tt>."%dir
                main_body+= self.help_window()
                main_body += "&nbsp;&nbsp;&nbsp;Create a file <tt>%s/index.html</tt> to replace this splash page.<br>"%(dir)
            interrupt_class = "interrupt_grey"
            worksheet = None
        else:
            worksheet = self.get_worksheet_with_id(worksheet_id)
            if worksheet.computing():
                interrupt_class = "interrupt"
            else:
                interrupt_class = "interrupt_grey"
            main_body = worksheet.html(authorized = worksheet_authorized)

        add_new_worksheet_menu = """
             <div class="add_new_worksheet_menu" id="add_worksheet_menu">
             Name: <input id="new_worksheet_box" class="add_new_worksheet_menu"
                    onKeyPress="if(is_submit(event)) process_new_worksheet_menu_submit();"></input><br>
             Password: <input id="new_worksheet_pass" class="add_new_worksheet_menu"
                    onKeyPress="if(is_submit(event)) process_new_worksheet_menu_submit();"></input>

             <button class="add_new_worksheet_menu"  onClick="process_new_worksheet_menu_submit();">add</button>
             <span class="X" onClick="hide_add_new_worksheet_menu()">X</span>
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
        body += '  <span class="banner"><a class="banner" target="_new" href="http://sage.math.washington.edu/sage">'
        body += '  <img src="sagelogo.png"/></a></span>\n'
        body += '  <span class="control_commands" id="cell_controls">\n'
        body += '    <a class="%s" onClick="interrupt()" id="interrupt">Interrupt</a>'%interrupt_class + vbar
        body += '    <a class="restart_sage" onClick="restart_sage()" id="restart_sage">Restart</a>' + vbar
        body += '    <a class="history_link" onClick="history_window()">History</a>' + vbar
        #body += '     <a onClick="toggle_left_pane()" class="worksheets_button" id="worksheets_button">Worksheets</a>' + vbar
        #body += '    <a class="doc_browser" onClick="show_doc_browser()">Documentation</a>' + vbar
        body += '    <a href="/doc_browser?/?index.html">Documentation</a>' + vbar
        body += '    <a class="help" onClick="show_help_window()">Help</a>' + vbar
        body += '    <a class="slide_mode" onClick="slide_mode()">Slideshow</a>'
        body += '  </span>\n'

        #these divs appear in backwards order because they're float:right
        body += '  <div class="hidden" id="slide_controls">\n'
        body += '    <div class="slideshow_control">'
        body += '      <a class="slide_arrow" onClick="slide_next()">&gt;</a>'
        body += '      <a class="slide_arrow" onClick="slide_last()">&gt;&gt;</a>' + vbar
        body += '      <a class="cell_mode" onClick="cell_mode()">Worksheet</a>'
        body += '    </div>'
        body += '    <div class="slideshow_progress" id="slideshow_progress" onClick="slide_next()">'
        body += '      <div class="slideshow_progress_bar" id="slideshow_progress_bar">&nbsp;</div>'
        body += '      <div class="slideshow_progress_text" id="slideshow_progress_text">&nbsp;</div>'
        body += '    </div>'
        body += '    <div class="slideshow_control">'
        body += '      <a class="slide_arrow" onClick="slide_first()">&lt;&lt;</a>'
        body += '      <a class="slide_arrow" onClick="slide_prev()">&lt;</a>'
        body += '    </div>'
        body += '  </span>\n'

        body += '</div>'
        body += '\n<div class="worksheet" id="worksheet">\n'
        if self.__show_debug or show_debug:
            body += "<div class='debug_window'>"
            body += "<div class='debug_output'><pre id='debug_output'></pre></div>"
            body += "<textarea rows=5 id='debug_input' class='debug_input' "
            body += " onKeyPress='return debug_keypress(event);' "
            body += " onFocus='debug_focus();' onBlur='debug_blur();'></textarea>"
            body += "</div>"

        body += main_body + '\n</div>\n'

        # The blank space given by '<br>'*15  is needed so the input doesn't get
        # stuck at the bottom of the screen. This could be replaced by a region
        # such that clicking on it creates a new cell at the bottom of the worksheet.
        body += '<br>'*15
        body += '\n</div>\n'

        body += '<div class="left_pane_bar" id="left_pane_bar" onClick="toggle_left_pane();"></div>\n'
        body += '<span class="pane" id="left_pane"><table bgcolor="white"><tr><td>\n'
        endpanespan = '</td></tr></table></span>\n'

        body += '  <div class="worksheets_topbar">'
        body += '     <a onClick="show_add_new_worksheet_menu()" class="new_worksheet">New</a> '
        body += '     <a onClick="show_delete_worksheet_menu()" class="delete_worksheet">Delete</a> '
        body += '  &nbsp;Worksheets</div>\n'
        body +=    add_new_worksheet_menu
        body +=    delete_worksheet_menu
        body += '  <div class="worksheet_list" id="worksheet_list">%s</div>\n'%self.worksheet_list_html(worksheet)

        if worksheet is None:
            return body + endpanespan

        body += '<div class="fivepix"></div>\n'
        body += '  <div class="objects_topbar"  onClick="toggle_menu(\'object_list\');">'
        body += '     <span class="plusminus" id="object_list_hider">[-]</span>'
        body += '     Saved Objects</div>\n'
        body += '  <div class="object_list" id="object_list">%s</div>\n'%self.object_list_html()
        body += '<div class="fivepix"></div>\n'
        body += '  <div class="variables_topbar" onClick="toggle_menu(\'variable_list\');">'
        body += '     <span class="plusminus" id="variable_list_hider">[-]</span>'
        body += '     Variables</div>\n'
        body += '  <div class="variable_list" id="variable_list">%s</div>\n'%\
                worksheet.variables_html()
        body += '<div class="fivepix"></div>\n'
        body += '  <div class="attached_topbar" onClick="toggle_menu(\'attached_list\');">'
        body += '     <span class="plusminus" id="attached_list_hider">[-]</span>'
        body += '     Attached Files</div>\n'
        body += '  <div class="attached_list" id="attached_list">%s</div><br>\n'%\
                worksheet.attached_html()
        body += endpanespan
        body += '<script language=javascript>focus(%s)</script>\n'%(worksheet[0].id())
        body += '<script language=javascript>jsmath_init();</script>\n'

        if worksheet_authorized:
            body += '<script language=javascript>worksheet_locked=false;</script>'
        else:
            body += '<script language=javascript>worksheet_locked=true;</script>'

        if worksheet.computing():
            # Set the update checking back in motion.
            body += '<script language=javascript> active_cell_list = %r; \n'%worksheet.queue_id_list()
            body += 'for(var i = 0; i < active_cell_list.length; i++)'
            body += '    cell_set_running(active_cell_list[i]); \n'
            body += 'start_update_check(); </script>\n'
        return body

    def edit_window(self, worksheet):
        """
        Return a window for editing worksheet.

        INPUT:
            worksheet -- a worksheet
        """
        t = worksheet.edit_text()
        t = t.replace('<','&lt;')
        body_html = ''
        body_html += '<h1 class="edit">SAGE Notebook: Editing Worksheet "%s"</h1>\n'%worksheet.name()
        body_html += """<b>Warnings:</b> You cannot undo after you save changes (yet).  All graphics will be deleted when you save.<br><br>"""
        body_html += '<form method="post" action="%s?edit" enctype="multipart/form-data">\n'%worksheet.filename()
        body_html += '<input type="submit" value="Save Changes" name="button_save"/>\n'
        #body_html += '<input type="submit" value="Preview" name="button_preview"/>\n'
        body_html += '<input type="submit" value="Cancel" name="button_cancel"/>\n'
        body_html += '<textarea class="edit" id="cell_intext" rows="30" name="textfield">'+t+'</textarea>'
        body_html += '</form>'
        body_html += """The format is as follows: <pre>
Arbitrary HTML
{{{
Input
///
Output
}}}
Arbitrary HTML
{{{
Input
///
Output
}}}
...
</pre>"""


        s = """
        <html><head><title>SAGE Wiki cell text </title>
        <style type="text/css">

        textarea.edit {
            font-family: monospace;
            border: 1px solid #8cacbb;
            color: black;
            background-color: white;
            padding: 3px;
            width: 100%%;
            margin-top: 0.5em;
        }
        </style>

        <script language=javascript> <!--

        %s

        function get_element(id) {
            if(document.getElementById)
                return document.getElementById(id);
            if(document.all)
                return document.all[id];
            if(document.layers)
                return document.layers[id];
        }

        function get_cell_list() {
            return window.opener.get_cell_list()
        }

        function send_doc_html() {
            var cell_id_list = get_cell_list();
            var num = cell_id_list.length;
            var lastid = cell_id_list[num-1];
            var doc_intext = get_element('cell_intext').value; /*for testing doc_html*/
            window.opener.upload_doc_html(lastid,doc_intext);
        }

        function send_cell_text() {
            var cell_id_list = get_cell_list();
            var num = cell_id_list.length;
            var lastid = cell_id_list[num-1];
            var cell_intext = get_element('cell_intext').value;
            window.opener.upload_cell_text(lastid,cell_intext);
        }

        function send_to_ws(do_eval) {
            var f = send_to_ws_callback;
            if (do_eval)
                f = send_to_ws_eval_callback;
            async_request('/delete_cell_all',f, "worksheet_id="+window.opener.get_worksheet_id());
        }

        function send_to_ws_callback(status, response_text) {
            window.opener.cell_delete_all_callback(status, response_text);
            var cell_intext = get_element('cell_intext').value;
            window.opener.insert_cells_from_wiki(cell_intext, false);
        }

        function send_to_ws_eval_callback(status, response_text) {
            window.opener.cell_delete_all_callback(status, response_text);
            var cell_intext = get_element('cell_intext').value;
            window.opener.insert_cells_from_wiki(cell_intext, true);
        }


        function clear_wiki_window() {
            get_element('cell_intext').value = ' ';
        }
        --></script></head>
        <body>%s
        </body></html>"""%(js.async_lib(), body_html)

        return s


    def help_window(self):
        help = [
            ('HTML', 'Begin an input block with %html and it will be output as HTML.  Use the &lt;sage>...&lt;/sage> tag to do computations in an HTML block and have the typeset output inserted.  Use &lt;$>...&lt;/$> and &lt;$$>...&lt;/$$> to insert typeset math in the HTML block.  This does <i>not</i> require latex.'),
            ('Shell', 'Begin a block with %sh to have the rest of the block evaluated as a shell script.  The current working directory is maintained.'),
            ('Autoevaluate Cells on Load', 'Any cells with "#auto" in the input is automatically evaluated when the worksheet is first opened.'),
            ('Create New Worksheet', "Use the menu on the left, or simply put a new worksheet name in the URL, e.g., if your notebook is at http://localhost:8000, then visiting http://localhost:8000/tests will create a new worksheet named tests."),
               ('Evaluate Input', 'Press shift-enter.  You can start several calculations at once.  If you press alt-enter instead, then a new cell is created after the current one.'),
                ('Timed Evaluation', 'Type "%time" at the beginning of the cell.'),
                ('Evaluate all Cells', 'Click <u>Eval All</u> in the upper right.'),
                ('Evaluate Cell using <b>GAP, Singular, etc.', 'Put "%gap", "%singular", etc. as the first input line of a cell; the rest of the cell is evaluated in that system.'),
                ('Typeset a Cell', 'Make the first line of the cell "%latex". The rest of the cell should be the body of a latex document.  Use \\sage{expr} to access SAGE from within the latex.  Evaluated typeset cells hide their input.  Use "%latex_debug" for a debugging version.  You must have latex for this to work.'),
               ('Typeset a slide', 'Same as typesetting a cell but use "%slide" and "%slide_debug"; will use a large san serif font.  You must have latex for this to work.'),
                ('Typesetting', 'Type "latex(objname)" for latex that you can paste into your paper.  Type "view(objname)" or "show(objname)", which will display a nicely typeset image (using javascript!).  You do <i>not</i> need latex for this to work.  Type "lprint()" to make it so output is often typeset by default.'),
                ('Move between cells', 'Use the up and down arrows on your keyboard.'),
                ('Interrupt running calculations',
                 'Click <u>Interrupt</u> in the upper right or press escape in any input cell. This will (attempt) to interrupt SAGE by sending many interrupts for several seconds; if this fails, it restarts SAGE (your worksheet is unchanged, but your session is reset).'),
                ('Tab completion', 'Press tab while the cursor is on an identifier. On some web browsers (e.g., Opera) you must use control-space instead of tab.'),
                ('Print worksheet', 'Click the print button.'),
                ('Help About',
                 'Type ? immediately after the object or function and press tab.'),
                ('Source Code',
                 'Put ?? after the object and press tab.'),
                ('Hide Cell Input',
                 'Put %hide at the beginning of the cell.  This can be followed by %gap, %latex, %maxima, etc.  Note that %hide must be first.  From the edit screen, use %hideall to hide a complete cell.'),
                ('Documentation',
                 'Click on <a href="/doc_browser?/?index.html">Documentation</a> in the upper right to browse the SAGE tutorial, reference manual, and other documentation.'),
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
                ('Downloading and Uploading Worksheets',
                 'Click <u>Download</u> in the upper right to download a complete worksheet to a local .sws file, and click <a href="__upload__.html">Upload</a> to upload a saved worksheet to the notebook.  Note that <i>everything</i> that has been submitted is automatically saved to disk when you quit the notebook server (or type "%save_server" into a cell).'),
                ('Restart', 'Type "restart" to restart the SAGE interpreter for a given worksheet.  (You have to interrupt first.)'),
                ('Input Rules', "Code is evaluated by exec'ing (after preparsing).  Only the output of the last line of the cell is implicitly printed.  If any line starts with \"sage:\" or \">>>\" the entire block is assumed to contain text and examples, so only lines that begin with a prompt are executed.   Thus you can paste in complete examples from the docs without any editing, and you can write input cells that contains non-evaluated plain text mixed with examples by starting the block with \">>>\" or including an example."),
                ('Working Directory', 'Each block of code is run from its own directory.  The variable DIR contains the directory from which you started the SAGE notebook.  For example, to open a file in that directory, do "open(DIR+\'filename\')".'),
                ('Customizing the look', 'Learn about cascading style sheets (CSS), then create a file notebook.css in your $HOME/.sage directory.  Use "view source" on a notebook web page to see the CSS that you can override.'),
                ('Emacs Keybindings', 'If you are using GNU/Linux, you can change (or create) a <tt>.gtkrc-2.0</tt> file.  Add the line <tt>gtk-key-theme-name = "Emacs"</tt> to it.  See <a target="_blank" href="http://kb.mozillazine.org/Emacs_Keybindings_(Firefox)">this page</a> [mozillazine.org] for more details.'),
                ('More Help', 'Type "help(sage.server.notebook.notebook)" for a detailed discussion of the architecture of the SAGE notebook and a tutorial (or see the SAGE reference manual).'),
                ('Javascript Debugger', 'Type ?debug at the end of a worksheet url to enable the javascript debugger.  A pair of textareas will appear at the top of the worksheet -- the upper of which is for output, the lower is a direct interface to the page\'s javascript environment.  Type any eval()-able javascript into the input box and press shift+enter to execute it.  Type debug_append(str) to print to, and debug_clear() to clear the debug output window.'),
                ]

        help.sort()
        s = """
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
        <h1 align=center><font color='darkred'>SAGE</font> Notebook Quickstart</h1>
        <div class="help_window">

        A <i>worksheet</i> is an ordered list of SAGE calculations with output.
        A <i>session</i> is a worksheet and a set of variables in some state.
        A <i>notebook</i> is a collection of worksheets and saved objects.

        <table class="help_window">
        """
        for x, y in help:
            s += '<tr><td class="help_window_cmd">%s</td><td class="help_window_how">%s</td></tr>'%(x,y)
        s += '</table></div>'

        s +="""
        <br>
        AUTHORS: William Stein, Tom Boothby, and Alex Clemesha (with feedback from many people,
        especially Fernando Perez and Joe Wetherell).<br><br>
        LICENSE: All code included with the standard SAGE install is <a href="__license__.html">licensed
        either under the GPL or a GPL-compatible license</a>.
        <br>
        """
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
              <h1><font size=+3 color="darkred">SAGE</font>&nbsp;&nbsp;&nbsp;&nbsp;<font size=+1>Upload your Worksheet</font></h1>
              <hr>
              <form method="POST" action="upload_worksheet"
                    name="upload" enctype="multipart/form-data">
              <table><tr>
              <td>
              Worksheet file:&nbsp&nbsp&nbsp </td>
              <td><input class="upload_worksheet_menu" size="40" type="file" name="fileField" id="upload_worksheet_filename"></input></td>
              </tr>
              <tr><td></td><td></td></tr>
              <tr>
              <td></td><td><input type="button" class="upload_worksheet_menu" value="Upload Worksheet" onClick="form.submit(); window.close();"></td>
              </tr>
              </form><br>
              </div>
            </body>
          </html>
         """%(css.css(self.color()),js.javascript())

    def doc_html(self,worksheet_id, css_href):
        try:
            W = self.get_worksheet_with_id(worksheet_id)
        except KeyError, msg:
            W = self.create_new_worksheet(worksheet_id)
            worksheet_id = W.id()
        head = self._doc_html_head(worksheet_id, css_href)
        body = self._doc_html_body(worksheet_id)
        if worksheet_id is not None:
           body += '<script language=javascript>worksheet_id="%s"; worksheet_filename="%s"; worksheet_name="%s"; toggle_left_pane(); </script>;'%(worksheet_id, W.filename(), W.name())

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html(self, worksheet_id=None, authorized=False, show_debug=False, worksheet_authorized=False):
        if worksheet_id is None or worksheet_id == '':
            if not self.splashpage():
                W = self.default_worksheet()
                worksheet_id = W.id()
            else:
                worksheet_id = None
                W = None
        else:
            try:
                W = self.get_worksheet_with_id(worksheet_id)
            except KeyError, msg:
                W = self.create_new_worksheet(worksheet_id)
                worksheet_id = W.id()

        if authorized:
            body = self._html_body(worksheet_id, show_debug=show_debug,
                                   worksheet_authorized=worksheet_authorized)
        else:
            body = self._html_authorize()

        if worksheet_id is not None:
            body += '<script language=javascript>worksheet_id="%s"; worksheet_filename="%s"; worksheet_name="%s";</script>;'%(worksheet_id, W.filename(), W.name())

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


import sage.interfaces.sage0
import time

## IMPORTANT!!! If you add any new input variable to notebook,
## you *must* similarly modify the restart_on_crash block
## at the beginning of the definition of notebook!!
def notebook(dir         ='sage_notebook',
             port        = 8000,
             address     = 'localhost',
             open_viewer = False,
             max_tries   = 10,
             username    = None,
             password    = None,
             color       = None,
             system      = None,
             jsmath      = True,
             show_debug  = False,
             splashpage  = None,
             warn        = True,
             ignore_lock = False,
             log_server = False,
             kill_idle   = 0,
             restart_on_crash = False,
             auto_restart = 1800):
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
        system -- (string) default computer algebra system to use for new
                  worksheets, e.g., 'maxima', 'gp', 'axiom', 'mathematica', 'macaulay2',
                  'singular', 'gap', 'octave', 'maple', etc.  (even 'latex'!)
        jsmath -- whether not to enable javascript typset output for math.
        debug -- whether or not to show a javascript debugging window
        kill_idle -- if positive, kill any idle compute processes after
                     this many auto saves.  (NOT IMPLEMENTED)

        splashpage -- whether or not to show a splash page when no worksheet is specified.
                      you can place a file named index.html into the notebook directory that
                      will be shown in place of the default.

        restart_on_crash -- if True (the default is False), the server
                      will be automatically restarted if it crashes in
                      any way.  Use this on a public servers that many
                      people might use, and which might be subjected
                      to intense usage.  NOTE: Log messages are only displayed
                      every 5 seconds in this mode.
        auto_restart -- if restart_on_crash is True, always restart
                      the server every this many seconds.

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
    if restart_on_crash:
        # Start a new subprocess
        def f(x):  # format for passing on
            if x is None:
                return 'None'
            elif isinstance(x, str):
                return "'%s'"%x
            else:
                return str(x)
        while True:
            S = sage.interfaces.sage0.Sage()
            time.sleep(1)
            S.eval("from sage.server.notebook.notebook import notebook")
            cmd = "notebook(dir=%s,port=%s, address=%s, open_viewer=%s, max_tries=%s, username=%s, password=%s, color=%s, system=%s, jsmath=%s, show_debug=%s, splashpage=%s, warn=%s, ignore_lock=%s, log_server=%s, kill_idle=%s, restart_on_crash=False)"%(
                f(dir), f(port), f(address), f(open_viewer), f(max_tries), f(username),
                f(password), f(color), f(system), f(jsmath), f(show_debug), f(splashpage),
                f(warn), f(ignore_lock), f(log_server), f(kill_idle)
                )
            print cmd
            S._send(cmd)
            tm = 0
            while True:
                s = S._get()[1].strip()
                if len(s) > 0:
                    print s
                if not S.is_running():
                    break
                time.sleep(5)
                tm += 5
                if tm > auto_restart:
                    S.quit()
                    break
            # end while
        # end while
        S.quit()
        return

    if os.path.exists(dir):
        if not os.path.isdir(dir):
            raise RuntimeError, '"%s" is not a valid SAGE notebook directory (it is not even a directory).'%dir
        if not (os.path.exists('%s/nb.sobj'%dir) or os.path.exists('%s/nb-backup.sobj'%dir)):
            raise RuntimeError, '"%s" is not a valid SAGE notebook directory (missing nb.sobj).'%dir
        if os.path.exists('%s/pid'%dir) and not ignore_lock:
            f = file('%s/pid'%dir)
            p = f.read()
            f.close()
            try:
                #This is a hack to check whether or not the process is running.
                os.kill(int(p),0)
                print "\n".join([" This notebook appears to be running with PID %s.  If it is"%p,
                                 " not responding, you will need to kill that process to continue.",
                                 " If another (non-sage) process is running with that PID, call",
                                 " notebook(..., ignore_lock = True, ...). " ])
                return
            except OSError:
                pass
        f = file('%s/pid'%dir, 'w')
        f.write("%d"%os.getpid())
        f.close()
        try:
            nb = load('%s/nb.sobj'%dir, compress=False)
        except:
            print "****************************************************************"
            print "  * * * WARNING   * * * WARNING   * * * WARNING   * * * "
            print "WARNING -- failed to load notebook data. Trying the backup file."
            print "****************************************************************"
            try:
                nb = load('%s/nb-backup.sobj'%dir, compress=False)
            except:
                print "Recovering from last op save failed."
                print "Trying save from last startup."
                nb = load('%s/nb-older-backup.sobj'%dir, compress=False)

        nb.set_directory(dir)
        if not (username is None):
            nb.set_auth(username=username, password=password)
        if not (color is None):
            nb.set_color(color)
        if not system is None:
            nb.set_system(system)
        if not splashpage is None:
            nb.set_splashpage(splashpage)
        nb.set_not_computing()
    else:
        nb = Notebook(dir,username=username,password=password, color=color,
                      system=system, kill_idle=kill_idle,splashpage=splashpage)
    nb.save()
    shutil.copy('%s/nb.sobj'%dir, '%s/nb-older-backup.sobj'%dir)
    nb.set_debug(show_debug)
    nb.set_log_server(log_server)
    if warn and address!='localhost' and username==None:
        print "WARNING -- it is *extremely* dangerous to let the server listen"
        print "on an external port without at least setting a username/password!!"
    nb.start(port, address, max_tries, open_viewer, jsmath=jsmath)
    from sage.interfaces.quit import expect_quitall
    expect_quitall(verbose=False)
    from sage.misc.misc import delete_tmpfiles
    delete_tmpfiles()
    if os.path.exists('%s/pid'%dir):
        os.remove('%s/pid'%dir)
    return nb


########################################################################
# NOTES ABOUT THE restart_on_crash option to notebook.
## I made some changes to the notebook command so it has the option of running
## the server as a completely separate process, which will restart it if it
## dies in any way.   I installed this on sage.math and started those servers
## running with this.  So if anybody has any systematic way to crash the notebook
## servers, please try it!  (The only one I have is to tell a computer lab
## full of 40 high school students to all try to crash the server at once...)
##
## Basically what happens is this:
##
##   1. You start a Python process then run the notebook command with the
##      restart_on_kill option True.
##   2. The notebook command starts another Python running, and in that
##      it runs the notebook command.  This uses the Sage0 pexpect interface.
##   3. It then monitors the process started in 2 -- if the process dies
##      or terminates, then 2 occurs again.  The server log output
##      is updated every 5 seconds.   If ctrl-c is received,
##      then everything cleans up and terminates.
##
## This means that the SAGE notebook can only currently be totally
## crashed if the Python simple http server gets into a hung state.
## From looking at the server logs after past crashes, this doesn't
## seem to ever happen.  So now instead of the server crashing
## under crazy loads, etc., it will automatically reset within seconds --
## this would of course kill all running worksheets (which is bad),
## but is much better than having the whole sage notebook go down
## until it is manually restarted!   In the long run, of course, using
## Twisted for the web server should hopefully mean that it doesn't
## crash...
###############################################################
