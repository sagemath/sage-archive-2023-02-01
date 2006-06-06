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

class Notebook(SageObject):
    def __init__(self, dir='sage_notebook'):
        self.__dir = dir
        self.__workbooks = {}
        self.__load_defaults()
        self.__filename     = '%s/nb.sobj'%dir
        self.__workbook_dir = '%s/workbooks'%dir
        self.__object_dir   = '%s/objects'%dir
        self.__makedirs()
        self.create_new_workbook()
        self.save()

    def directory(self):
        return self.__dir

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

    def objects_html(self):
        s = ''
        div = '<div class="object_name" onClick="click_on_object(\'%s\');">'
        for name in self.objects():
            s += div%name + name + '</div>'
        return s

    def defaults(self):
        return self.__defaults

    def __makedirs(self):
        os.makedirs(self.__dir)
        os.makedirs(self.__workbook_dir)
        os.makedirs(self.__object_dir)

    def create_new_workbook(self, name='untitled'):
        name = str(name)
        self.__workbooks[name] = workbook.Workbook(name, self)

    def current_workbook(self):
        return self.__workbooks[self.current_workbook_name()]

    def current_workbook_name(self):
        try:
            return self.__current_workbook_name
        except AttributeError:
            self.__current_workbook_name = self.__workbooks.keys()[0]
            return self.__current_workbook_name

    def interrupt(self):
        """
        Interrupt currently queued up calculation in current workbook.

        OUTPUT:
            bool -- return True if no problems interrupting calculation
                    return False if the SAGE interpreter had to be restarted.
        """
        W = self.current_workbook()
        return W.interrupt()

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

        if open_viewer:
            os.system('%s http://%s:%s 1>&2 >/dev/null &'%(BROWSER, address, port))
        ns.serve()
        self.save()

    def _html_head(self):
        head = '<title>SAGE Notebook: %s</title>'%self.directory()
        head += '<style>' + css.css() + '</style>\n'
        head += '<script language=javascript>' + js.javascript() + '</script>\n'
        return head

    def _html_body(self):
        workbook = self.current_workbook()
        body = ''
        body += self._help_window()
        body += '<span class="top_control_bar">\n'
        #body += '  <input class="search_fulltext"></input>\n'
        body += '  <div class="workbook_title">%s</div>\n'%workbook.name()
        body += '  <a class="evaluate" onClick="evaluate_all()">Evaluate All</a>\n'
        body += '  <a class="interrupt_grey" onClick="interrupt()" id="interrupt">Interrupt</a>\n'
        body += '  <a class="help" onClick="show_help_window()">Help</a>\n'
        body += '  <a class="hide" onClick="hide_all()">Hide Output</a>\n'
        body += '</span>'
        body += '\n<div class="workbook" id="workbook">\n' + workbook.html() + '\n</div>\n'

        body += '\n<span class="banner">SAGE</span>\n'
        body += '<span class="pane"><table bgcolor="white"><tr><td>\n'
        body += '  <br><div class="variables_topbar">Variables</div>\n'
        body += '  <div class="variables_list" id="variable_list"></div><br>\n'
        body += '  <div class="attached_topbar">Attached Files</div>\n'
        body += '  <div class="attached_list" id="attached_list"></div><br>\n'
        body += '  <div class="workbooks_topbar">Workbooks</div>\n'
        body += '  <div class="workbooks_list">untitled</div><br>\n'
        body += '  <div class="objects_topbar">Saved Objects</div>\n'
        body += '  <div class="objects_list" id="object_list">%s</div>\n'%self.objects_html()
        body += '</td></tr></table></span>\n'
        body += '<script language=javascript>focus_on(%s)</script>\n'%(workbook[0].id())

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


    def html(self):
        head = self._html_head()
        body = self._html_body()
        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

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






