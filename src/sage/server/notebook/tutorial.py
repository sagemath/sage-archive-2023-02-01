"""nodoctest
"""

#############################################################################
#       Copyright (C) 2006, 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

r"""
SAGE Notebook Interface

AUTHORS:
    -- William Stein (2006-05-06): initial version
    -- Alex Clemesha
    -- Tom Boothby: support for a wide range of web browsers;
       refactoring of javascript code; systematic keyboard controls
    -- Dorian Raymer
    -- Yi Qiang
    -- Bobby Moretti

The SAGE graphical user interface is unusual in that it operates via
your web browser.  It provides you with SAGE worksheets that you can
edit and evaluate, which contain scalable typeset mathematics and
beautiful antialised images.  To try it out immediately, do this:

    sage: notebook(open_viewer=True)          # not tested
    the sage notebook starts...

\subsection{Supported Browsers}

The SAGE notebook should fully work with Firefox (and Mozilla).  It
may work to some extent with Safari and Opera.  Internet Explorer is
not supported.

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
delete all its contents, then press ctrl-backspace one more time.  The cell
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
starts with "sage:" or ">>>" the \emph{entire block} is assumed to
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

You can save all variables in a current session using the
\code{save_session} command, and you can then load those session
variables using the \code{load_session} command.

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




#####################################

notebook_help = [
            ('Full Text Search of Docs and Source', 'Search the SAGE documentation by typing <pre>search_doc("my query")</pre> in an input cell and press shift-enter.  Search the source code of SAGE by typing <pre>search_src("my query")</pre> and pressing shift-enter.  Arbitrary regular expressions are allowed as queries.'),
            ('HTML', 'Shift click between cells to create a new HTML cell.  Double click on existing HTML to edit it.  Use $...$ and $$...$$ to include typeset math in the HTML block.'),
            ('Shell', 'Begin a block with %sh to have the rest of the block evaluated as a shell script.  The current working directory is maintained.'),
            ('Interactive Dynamic Widgets', 'Put @interact on the line before a function definition.  Type interact? for more details.'),
            ('Autoevaluate Cells on Load', 'Any cells with "#auto" in the input is automatically evaluated when the worksheet is first opened.'),
               ('Evaluate Input', '<b>Press shift-enter.</b>  You can start several calculations at once.  If you press alt-enter instead, then a new cell is created after the current one.  If you press control-enter then the cell is split and both pieces are evaluated separately.'),
                ('Time', 'Type "%time" at the beginning of the cell.'),
                ('Evaluate Cell using <b>GAP, Singular, etc.', 'Put "%gap", "%singular", etc. as the first input line of a cell; the rest of the cell is evaluated in that system.'),
                ('Interrupt running calculations',
                 'Click <u>Interrupt</u> or press escape in any input cell. This will (attempt) to interrupt SAGE by sending many interrupt signals.'),
                ('Tab Completion', 'Press tab while the cursor is on an identifier. On some web browsers (e.g., Opera) you must use control-space instead of tab.'),
                ('Typesetting All Output', 'Type pretty_print_default() in an input cell and press shift-enter.  All future output will be typeset automatically.'),
                ('Help About',
                 'Type ? immediately after the object or function and press tab.'),
                ('Source Code',
                 'Put ?? after the object and press tab.'),
                ('Indenting Blocks',
                 'Highlight text and press > to indent it all and < to unindent it all (works in Safari and Firefox).  In Firefox you can also press tab and shift-tab.'),
                ('Insert New Cell',
                 'Put the mouse between an output and input until the horizontal line appears and click.  If you press Alt-Enter in a cell, the cell is evaluated and a new cell is inserted after it.'),
                ('Delete Cell',
                 'Delete all cell contents, then the press ctrl-backspace. This is a special case of joining cells.'),
                ('History', 'Click <a href="/history">log</a> commands you have entered in any worksheet of this notebook.'),
                ('Hide/Show Output', 'Click on the left side of output to toggle between hidden, shown with word wrap, and shown without word wrap.'),
                ('Loading and Saving Sessions', 'Use "save_session(\'name\')" to save all variables to an object.  Use "load_session(\'name\')" to <i>merge</i> in all variables from a saved session.'),
                ('Loading and Saving Objects', 'Use "save obj1 obj2 ..." and "load obj1 obj2 ...".  This allows for easy moving of objects from one worksheet to another, and saving of objects for later use.'),
                ('Loading SAGE/Python Scripts', 'Use "load filename.sage" and "load filename.py".  Load is relative to the path you started the notebook in.  The .sage files are preparsed and .py files are not.   You may omit the .sage or .py extension.  Files may load other files.'),
                ('Attaching Scripts', 'Use "attach filename.sage" or "attach filename.py".  Attached files are automatically reloaded when the file changes.  The file $HOME/.sage/init.sage is attached on startup if it exists.'),
                ('Restart', 'Type "restart" to restart the SAGE interpreter for a given worksheet.  (You have to interrupt first.)'),
                ('Input Rules', "Code is evaluated by exec'ing (after preparsing).  Only the output of the last line of the cell is implicitly printed.  If any line starts with \"sage:\" or \">>>\" the entire block is assumed to contain text and examples, so only lines that begin with a prompt are executed.   Thus you can paste in complete examples from the docs without any editing, and you can write input cells that contains non-evaluated plain text mixed with examples by starting the block with \">>>\" or including an example."),
                ('Working Directory', 'Each block of code is run from its own directory.  If any images are created as a side effect, they will automatically be displayed.'),
            ('DIR variable', 'The variable DIR contains the directory from which you started the SAGE notebook.  For example, to open a file in that directory, do "open(DIR+\'filename\')".'),
            ('DATA variable', 'The variable DATA contains the directory with data files that you upload into the worksheet.  For example, to open a file in that directory, do "open(DATA+\'filename\')".'),
            ('Split and join cells', 'Press ctrl-; in a cell to split it into two cells, and ctrl-backspace to join them.  Press ctrl-enter to split a cell and evaluate both pieces.'),
#                ('Emacs Keybindings', 'If you are using GNU/Linux, you can change (or create) a <tt>.gtkrc-2.0</tt> file.  Add the line <tt>gtk-key-theme-name = "Emacs"</tt> to it.  See <a target="_blank" href="http://kb.mozillazine.org/Emacs_Keybindings_(Firefox)">this page</a> [mozillazine.org] for more details.'),
 #               ('More Help', 'Type "help(sage.server.notebook.notebook)" for a detailed discussion of the architecture of the SAGE notebook and a tutorial (or see the SAGE reference manual).'),
#                ('Javascript Debugger', 'Type ?debug at the end of a worksheet url to enable the javascript debugger.  A pair of textareas will appear at the top of the worksheet -- the upper of which is for output, the lower is a direct interface to the page\'s javascript environment.  Type any eval()-able javascript into the input box and press shift+enter to execute it.  Type debug_append(str) to print to, and debug_clear() to clear the debug output window.'),
                ]

notebook_help.sort()
