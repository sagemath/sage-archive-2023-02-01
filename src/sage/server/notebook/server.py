"""
Web Server Component of SAGE Notebook

In this file the WebServer class is defined, which
inherits from Python's BaseHTTPRequestHandler.

The main purpose of the WebServer is to handle
HTTP GET and POST requests (as with any web server).

For the GET requests, the request path is examined
and 'static' files are served based on the path.
If, for example, one wanted to extend the functionality
of the Notebook by enabling the WebServer to handle and
serve and new kind of GET request, one would write a function
that writes to the 'wfile' (the outgoing file-like object)
some data depending the incoming path.

For the POST requests, the requests comes with some
post variables.  We define functions in the WebServer class
to handle these POST requests with the post variables.
If, for example, one wanted to extend POST functionality,
one would write a function that takes a input some post variables
and serves a request depending on those post vars.
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


# Python libraries
import BaseHTTPServer
import socket
import cgi
import mimetypes
from exceptions import KeyError
import os, sys
import select
from   StringIO import StringIO
import shutil
import Cookie
import cPickle
import base64
from gzip import GzipFile
import struct
from urllib import splittag

#SAGE notebook libraries
import css, js
import keyboards
from docHTMLProcessor import DocHTMLProcessor

# SAGE libraries
import sage.interfaces.sage0

# doc browser sequence number
doc_browser_seq = 0
MAX_DOC_BROWSER_SEQ = 48  # cycle around after this many worksheets -- max number
                          # of truly simult doc requests.

from sage.misc.misc import (alarm, cancel_alarm,
                            verbose, word_wrap, SAGE_EXTCODE)

import sage.misc.copying

import sage.misc.preparser
from   sage.structure.sage_object import load, SageObject

SEP = '___S_A_G_E___'

# The unique notebook that is being served by the web server.
# This is a global variable that is set by the NotebookServer
# constructor.  It is equal to the last notebook that the
# serve method was called on, and is set back to None
# after the serve method terminates.

notebook = None

import time

SAVE_INTERVAL=5   # time in seconds between saves when notebook is in use.
last_save_time = time.time()

SAGE_ROOT = os.environ['SAGE_ROOT']

static_images = ['favicon.ico', 'corner.png', 'evaluate.png', 'evaluate_over.png', 'sagelogo.png']

def safe_path(path):
    """
    Return a safe version of the given path, i.e., a relative path with no ..'s.
    The idea is to make it so the server can't just easily return arbitrary files
    that it has access to.   Of course, right now in SAGE one can execute a command
    in a cell that uses os.system to look at anything.  But in the future the
    subcommands will be run in a sandbox themselves, so protecting the server
    is still relevant.
    """
    return path.lstrip('/').replace('..','dotdot_not_allowed')



class WebServer(BaseHTTPServer.BaseHTTPRequestHandler):
    def get_postvars(self):
        r"""
        This function gets the variables corresponding to a POST.

        Specifically, it is evaluated after \code{do_POST} is
        called, which in turn calls some other more specialized
        function that handles the POST using the POST variables.

        \code{length} is the length of the header in bytes. It is used
        to read the \code{rfile} just enough to extract the POST variables
        \code{rfile} contains the input stream, and is of type <class 'socket._fileobject'>

        \code{cgi.parse_qs} return a dictionary mapping POST variable
        names to lists of their values.

        We then loop through the keys of the cgi dictionary
        and extract the actual used values.

        """
        length = int(self.headers.getheader('content-length'))
        qs = self.rfile.read(length)
        pqs =  cgi.parse_qs(qs, keep_blank_values=1)
        postvars = {}
        #extract from the cgi dict the useable values
        for var in pqs.keys():
            try:
                if var in ['id', 'cell_id']:
                    postvars[var] = int(pqs[var][0])
                else:
                    postvars[var] = pqs[var][0]
            except KeyError:
                pass
        if notebook.log_server():
            notebook.server_log().append(["POST", self.path, postvars, self.headers.getheader('content-type')])
        return postvars

    def cell_output_set(self):
        """
        Set the output type of the cell.

        This enables the type of output cell,
        such as to allowing wrapping for output
        that is very long.
        """
        C = self.get_postvars()
        id = C['id']
        typ = C['type']
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        if self.auth_worksheet(W):
            cell = W.get_cell_with_id(id)
            cell.set_cell_output_type(typ)

    def set_cell_input(self):
        C = self.get_postvars()
        id = C['cell_id']
        input = C['input']

        W = notebook.get_worksheet_that_has_cell_with_id(id)
        if self.auth_worksheet(W):
            cell = W.get_cell_with_id(id)
            cell.set_input_text(input)

    def hide_all(self):
        """
        Sets every cell's output hidden in a given worksheet.
        """
        ws_id = self.get_postvars()['worksheet_id']
        W = notebook.get_worksheet_with_id(ws_id)
        if self.auth_worksheet(W):
            W.hide_all()

    def show_all(self):
        """
        Sets every cell's output visible in a given worksheet.
        """
        ws_id = self.get_postvars()['worksheet_id']
        W = notebook.get_worksheet_with_id(ws_id)
        if self.auth_worksheet(W):
            W.show_all()


    def restart_sage(self):
        """
        Restart a given worksheet session.

        All defined variables and functions will be
        removed from the namespace of this worksheet.

        """
        ws_id = self.get_postvars()['worksheet_id']
        W = notebook.get_worksheet_with_id(ws_id)
        if self.auth_worksheet(W):
            W.restart_sage()
            self.wfile.write('done')

    def eval_cell(self, newcell=False, introspect=False):
        """
        Evaluate a cell.

        If the request is not authorized,
        (the requester did not enter the correct password
        for the given worksheet), then the request to
        evaluate or introspect the cell is ignored.

        If the cell contains either 1 or 2 question marks,
        then this is interpreted as a request for either
        introspection to the documentation of the function,
        or the documentation of the function and the
        source code of the function respectively.
        """
        C = self.get_postvars()
        id = C['id']
        input_text = C['input']
        input_text = input_text.replace('\r\n', '\n') #TB: dos make crazy
        #input_text = input_text.replace("%2B",'+')
        verbose('%s: %s'%(id, input_text))
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        if not self.auth_worksheet(W):
            return

        cell = W.get_cell_with_id(id)
        if not introspect:
            cell.set_input_text(input_text)
        cell.evaluate(introspect=introspect)

        if cell.is_last():
            new_cell = W.append_new_cell()
            self.wfile.write(str(new_cell.id()) + SEP + 'append_new_cell' + SEP + \
                             new_cell.html(div_wrap=False))
        elif newcell:
            new_cell = W.new_cell_after(id)
            self.wfile.write(str(new_cell.id()) + SEP + 'insert_cell' + SEP + \
                             new_cell.html(div_wrap=False) + SEP + str(id))
        else:
            self.wfile.write(str(cell.next_id()) + SEP + 'no_new_cell' + SEP + str(id))

    def introspect(self):
        C = self.get_postvars()
        id = C['id']
        before_cursor = C['before_cursor']
        after_cursor = C['after_cursor']
        input_text = (before_cursor+after_cursor)
        verbose('introspect -- %s: %s|%s'%(id, before_cursor, after_cursor))

        W = notebook.get_worksheet_that_has_cell_with_id(id)
        if not self.auth_worksheet(W):
            return

        cell = W.get_cell_with_id(id)
        #TB: this tends to obliterate long cells -- if the user doesn't submit between
        #introspecting and closing the browser; there's a lot of potential to lose a
        #large amount of work without warning.  I personally would not expect hitting
        #tab to save the input.
        #cell.set_input_text(before_cursor + after_cursor)
        cell.evaluate(introspect=[before_cursor, after_cursor])

        self.wfile.write(str(cell.next_id()) + SEP +
                         'no_new_cell' + SEP + str(id))


    def new_cell(self):
        """
        Add a new cell before a given cell.
        """
        id = self.get_postvars()['id']
        verbose("Adding new cell before cell with id %s"%id)
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        if not self.auth_worksheet(W):
            return

        cell = W.new_cell_before(id)
        self.wfile.write(str(cell.id()) + SEP + cell.html(div_wrap=False) + SEP + \
                         str(id))

    def new_cell_after(self):
        """
        Add a new cell after a given cell.
        """
        id = self.get_postvars()['id']
        verbose("Adding new cell after cell with id %s"%id)
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        if not self.auth_worksheet(W):
            self.wfile.write("locked")
            return

        cell = W.new_cell_after(id)
        self.wfile.write(str(cell.id()) + SEP + cell.html(div_wrap=False) + SEP + \
                         str(id) + SEP)

    def delete_cell(self):
        """
        Deletes a notebook cell.

        If there is only one cell left in a given
        worksheet, the request to delete that cell
        is ignored because there must be a least
        one cell at all times in a worksheet.
        (This requirement exists so other functions
        that evaluate relative to existing cells will
        still work ... this requirement could be removed?)

        """
        id = self.get_postvars()['id']
        verbose("Deleting cell with id %s"%id)
        W = notebook.get_worksheet_that_has_cell_with_id(id)

        if not self.auth_worksheet(W):
            return

        if len(W) <= 1:
            self.wfile.write('ignore')
        else:
            prev_id = W.delete_cell_with_id(id)
            self.wfile.write('delete' + SEP + str(id) + SEP + str(prev_id) + SEP + str(W.cell_id_list()))

    def delete_cell_all(self):
        C = self.get_postvars()
        worksheet_id = C['worksheet_id'][0]
        W = notebook.get_worksheet_with_id(worksheet_id)
        if not self.auth_worksheet(W):
            return

        cells = W.cell_id_list()[1:]
        cells.reverse()
        for cell in cells:
            W.delete_cell_with_id(cell)
            print cell
        cell = W[0]
        cell.set_input_text("")
        cell.set_output_text("", "")
        self.wfile.write("OK")

    def save_notebook_every_so_often(self):
        """
        Writes to the nb.sobj every SAVE_INTERVAL time step.

        """
        global last_save_time
        if time.time() - last_save_time > SAVE_INTERVAL:
            notebook.save()
            last_save_time = time.time()

    def kill_idle_every_so_often(self):
        notebook.kill_idle_compute_processes()

    def cell_update(self):
        """
        Updates an evaluated cell.

        If the cell contains a calculation that takes
        a long time to complete, the async javascript will
        continously request the evaluation of this function
        to serve two purposes:
        1) To see if the long calculation is done, or
        2) To interrupt the long running calculation.
        """
        C = self.get_postvars()
        worksheet_id = C['worksheet_id']
        cell_id = C['cell_id']

        worksheet = notebook.get_worksheet_with_id(worksheet_id)
        cols = notebook.defaults()['word_wrap_cols']

        # update the computation one step.
        worksheet.check_comp()
        # now get latest status on our cell
        status, cell = worksheet.check_cell(cell_id)

        #print status, cell   # debug
        if status == 'd':
            variables = worksheet.variables_html()
            objects = notebook.object_list_html()
            attached_files = worksheet.attached_html()
        else:
            variables = '...'  # not used
            objects = "..." # not used
            attached_files = '...' # not used

        if status == 'd':
            new_input = cell.changed_input_text()
            out_html = cell.output_html()
        else:
            new_input = ''
            out_html = ''
        if cell.interrupted():
            inter = 'true'
        else:
            inter = 'false'

        raw = cell.output_text(raw=True).split("\n")
        if len(raw) == 13 and raw[4][:17] == "Unhandled SIGSEGV":
            inter = 'restart'
            print "segfault!"

        msg = '%s%s %s'%(status, cell.id(),
                          SEP.join([cell.output_text(html=True),
                                    cell.output_text(cols, html=True),
                                    out_html,
                                    new_input,
                                    inter,
                                    variables,
                                    objects,
                                    attached_files,
                                    cell.introspect_html()]))

        # more comps to go ?
        worksheet.start_next_comp()
        self.wfile.write(msg)

    def interrupt(self):
        ws_id = self.get_postvars()['worksheet_id']
        W = notebook.get_worksheet_with_id(ws_id)
        if not self.auth_worksheet(W):
            return

        t = W.interrupt()
        if t:
            self.wfile.write('ok')
        else:
            self.wfile.write('restart')

    def add_worksheet(self):
        C = self.get_postvars()
        worksheet_name = C['name']
        passcode = C['passcode']
        try:
            W = notebook.create_new_worksheet(worksheet_name, passcode)
        except ValueError, msg:
            print msg
            self.wfile.write(msg)
            return
        new_worksheet_list = notebook.worksheet_list_html(W.name())
        self.wfile.write(new_worksheet_list + SEP + str(W.name()))

    def auth_worksheet(self, W):
        n = W.filename()
        passcode = ''
        if self.cookie.has_key('ws_%s_passcode'%n):
            passcode = self.cookie['ws_%s_passcode'%n].value
        return W.auth(passcode)

    def unlock_worksheet(self):
        ws_id = self.get_postvars()['worksheet_id']
        W = notebook.get_worksheet_with_id(ws_id)
        if not self.auth_worksheet(W):
            self.wfile.write('failed')
        else:
            self.wfile.write('ok')

    def delete_worksheet(self):
        C = self.get_postvars()
        worksheet_name = C['name']
        try:
            W = notebook.get_worksheet_with_name(worksheet_name)
        except KeyError:
            # it is already deleted.
            msg = "No such worksheet '%s'"%worksheet_name
            self.wfile.write(msg)
            return
        if not self.auth_worksheet(W):
            msg = "Error deleting worksheet '%s' (you must login to it first): "%worksheet_name
            self.wfile.write(msg)
            return

        try:
            W = notebook.delete_worksheet(worksheet_name)
        except KeyError, msg:
            msg = "Error deleting worksheet: " + str(msg)
            self.wfile.write(msg)
            return
        new_worksheet_list = notebook.worksheet_list_html(W.name())
        self.wfile.write(new_worksheet_list + SEP + str(W.name()))

    def import_worksheet_local_file(self, filename):
        try:
            W = notebook.import_worksheet(filename)
        except ValueError, msg:
            msg = "Error uploading worksheet: " + str(msg)
            print msg
            self.wfile.write(msg)
            return
        #self.wfile.write(notebook.worksheet_list_html())
        #self.wfile.write(notebook.html(W.id(), authorized=self.authorize()))
        self.send_response(302)
        self.send_header("Location", '/%s'%W.name())
        self.end_headers()

    #######################################################################
    #  Source -browser functionality
    #######################################################################

    def src_browser(self, path):
        self.send_head()
        i = path.find('?')
        if i == -1:
            return self.file_not_found(path)
        filename = path[i+1:]
        file = open('%s/devel/sage/sage/%s'%(SAGE_ROOT,safe_path(filename))).read()
        file = file.replace('<','&lt;')
        s = """
<html>
<head>
"""
        s += '<title>%s | SAGE Source Browser</title>' % filename
        s += '<link rel=stylesheet href="/highlight/prettify.css" type="text/css" />\n'
        s += """
</head>
<body>
"""
        s += '<h1 align=center>SAGE Source Browser</h1>\n'
        s += '<h2 align=center>devel/sage/sage%s</h2>\n'%filename
        s += '<br><hr><br>\n'
        s += '<pre id="code">%s</pre>\n'%file
        s += '<br><hr><br>\n'
        s += '<script src="/highlight/prettify.js" type="text/javascript"></script>\n'
        s += """<script type="text/javascript">
function get_element(id) {
  if(document.getElementById)
    return document.getElementById(id);
  if(document.all)
    return document.all[id];
  if(document.layers)
    return document.layers[id];
}

var x = get_element("code");
x.innerHTML = prettyPrintOne(x.innerHTML);
</script>
"""
        s += """
</body>
</html>
"""
        return self.wfile.write(s)

    #######################################################################
    #  Doc-browser functionality
    #######################################################################

    def doc_browser(self, path):
        """
        This is the first handler function for a doc-browser request.
        If no doc-browser has been opened/created, one is started.
        The doc-browser is an instance of a special worksheet; each
        page of documentation requested is formated into a worksheet
        which replaces the doc-browsers previous worksheet.
        """
        path_split = path.split('?')
        full_path = path_split[1]
        file_name = path_split[2]
        file_name, href_tag = splittag(file_name)

        # Get the documentation path from the symlink in <SAGE_ROOT>
        # doc_path is the path from <SAGE_ROOT> to the top doc folder
        #  full_path is the path from the top doc folder to the requested
        #  file
        # file_name is the name of the file requested
        doc_path = os.path.abspath(os.environ['SAGE_ROOT'] + '/doc/')
        file = open(doc_path + '/' + safe_path(full_path + file_name),'r')
        doc_page_html = file.read()
        file.close()
        docProcessStart = time.time()
        doc_page,css_href = DocHTMLProcessor().process_doc_html(doc_path,full_path,doc_page_html)
        docProcessEnd = time.time()
        docProcessTime = docProcessEnd-docProcessStart
        verbose(file_name)
        verbose(docProcessTime)
        if css_href:
            css_href = doc_path + full_path + css_href

        wnames = notebook.worksheet_names()
        global doc_browser_seq
        worksheet_name = 'doc_browser_%s'%doc_browser_seq
        doc_browser_seq += 1
        if doc_browser_seq >= MAX_DOC_BROWSER_SEQ:
            doc_brower_seq = 0
        if worksheet_name not in wnames:
            W = notebook.create_new_worksheet(name=worksheet_name, passcode='')
        else:
            W = notebook.get_worksheet_with_name(worksheet_name)
            W.interrupt()
        W.edit_save(doc_page)
        saveTime = time.time()
        saveTime = saveTime - docProcessEnd
        verbose(saveTime)
        self.load_doc_page(W, css_href)

    def load_doc_page(self, worksheet, css_href):
        W = worksheet
        Wid = W.id()
        #s = notebook.doc_html(Wid, css_href)
        #self.wfile.write(s)
        self.show_page(Wid)
        #s = notebook.html(Wid)

    #######################################################################
    #  End doc-browser functionality
    #######################################################################

    #######################################################################
    #  SAGE plain text editing functionality
    #######################################################################

    def edit_text(self, filename, prompts=False):
        """
        Return a window that allows the user to edit the text of the
        worksheet with the given filename.
        """
        self.send_head()
        W = notebook.get_worksheet_with_filename(filename)
        s = notebook.edit_window(W)
        self.wfile.write(s)

    def edit_save(self, filename, newtext):
        W = notebook.get_worksheet_with_filename(filename)
        if self.auth_worksheet(W):
            W.edit_save(newtext)
        return self.show_page(worksheet_id=W.id())

    def edit_preview(self):
        raise NotImplementedError

    def edit_cancel(self, filename):
        W = notebook.get_worksheet_with_filename(filename)
        return self.show_page(worksheet_id=W.id())


    #######################################################################
    #  End editing functionality
    #######################################################################

    def get_queue(self):
        C = self.get_postvars()
        W = notebook.get_worksheet_with_filename(C['worksheet_id'])
        a = ",".join(["%s"%q for q in W.queue_id_list()])
        self.wfile.write(a)

    def plain_text_worksheet(self, filename, prompts=True):
        self.send_head()
        W = notebook.get_worksheet_with_filename(filename)
        t = W.plain_text(prompts = prompts)
        t = t.replace('<','&lt;')
        s = '<head>\n'
        s += '<title>SAGE Worksheet: %s</title>\n'%W.name()
        s += '</head>\n'
        s += '<body>\n'
        s += '<pre>' + t + '</pre>'
        s += '</body>\n'
        self.wfile.write(s)

    def html_worksheet(self, filename, do_print=False):
        self.send_head()
        W = notebook.get_worksheet_with_filename(filename)
        s = '<head>\n'
        s += '<title>SAGE Worksheet: %s</title>\n'%W.name()
        if do_print:
            s += '<script type="text/javascript" src="/jsmath/jsMath.js"></script>\n'
        s += '<script type="text/javascript" src="/__main__.js"></script>\n'
        s += '<link rel=stylesheet href="/__main__.css">\n'
        s += '</head>\n'
        s += '<body>\n'
        s += W.html(include_title=False, do_print=do_print)
        if do_print:
            s += '<script type="text/javascript">jsMath.Process();</script>\n'
        s += '\n</body>\n'
        self.wfile.write(s)

    def input_history_text(self):
        self.send_head()
        t = notebook.history_text()
        t = t.replace('<','&lt;')
        s = '<head>\n'
        s += '<title>SAGE Input History</title>\n'
        s += '</head>\n'
        s += '<body>\n'
        s += '<pre>' + t + '</pre>\n'
        s += '<a name="bottom"></a>\n'
        s += '<script type="text/javascript"> window.location="#bottom"</script>\n'
        s += '</body>\n'
        self.wfile.write(s)

    def help_window(self):
        self.send_head()
        self.wfile.write(notebook.help_window())

    def license_window(self):
        self.send_head()
        s = "<html><head><title>SAGE License</title></head>"
        s+= "<body><font size='-1'><pre>%s</pre></font></body></html>"%sage.misc.copying.license
        self.wfile.write(s)

    def upload_window(self):
        self.send_head()
        self.wfile.write(notebook.upload_window())

##     def insert_wiki_cells(self):
##         C = self.get_postvars()
##         W = notebook.get_worksheet_with_id(int(C['worksheet_id'][0]))
##         W.insert_wiki_cells(C['text'][0])
##         response = C['eval'][0] + SEP
##         response+= "%r"%W.cell_id_list() + SEP
##         response+= SEP.join([c.html(div_wrap=False) for c in W[:-1]])
##         self.wfile.write(response)

    def download_worksheet(self, filename):
        try:
            notebook.export_worksheet(filename, filename)
        except KeyError:
            self.file_not_found(filename)
            return
        self.send_response(200)
        self.send_header("Content-type", 'application/sage')
        self.end_headers()
        binfile = open('%s/%s.sws'%(notebook.directory(), safe_path(filename)), 'rb').read()
        f = StringIO()
        f.write(binfile)
        f.seek(0)
        shutil.copyfileobj(f, self.wfile)
        f.close()

    def cell_id_list(self):
        C = self.get_postvars()
        worksheet_id = C['worksheet_id'][0]
        worksheet = notebook.get_worksheet_with_id(worksheet_id)
        L = worksheet.cell_id_list()
        s = ' '.join(str(x) for x in L)
        self.wfile.write(s)

    def get_file(self):
        """
        Examine the request and serve a file based on
        if the path matches any of the below strings.

        If the path matches none of the strings,
        then a '200 File Not Found' error is return.
        """
        compressed = False
        path = self.path.replace('%20',' ')
        if path.endswith('.sobj'):
            path = '%s/%s'%(os.path.abspath(notebook.object_directory()), path)
        else:
            path = safe_path(path)
        if path.endswith('.html') and not '/' in path and not '/jsmath' in path and not '/highlight' in path:
            worksheet_filename = path[:-5]
            if worksheet_filename == '__history__':
                self.input_history_text()
            elif worksheet_filename == '__help__':
                self.help_window()
            elif worksheet_filename == '__license__':
                self.license_window()
            elif worksheet_filename[-7:] == '__doc__':
                self.plain_text_worksheet(worksheet_filename[:-7], prompts=True)
            elif worksheet_filename[-9:] == '__plain__':
                self.plain_text_worksheet(worksheet_filename[:-9], prompts=False)
            elif worksheet_filename[-9:] == '__print__':
                self.html_worksheet(worksheet_filename[:-9], do_print=True)
            elif worksheet_filename[-10:]== '__upload__':
                self.upload_window()
            else:
                self.html_worksheet(worksheet_filename, do_print=False)
            return

        elif path == 'robots.txt':
            self.wfile.write(self.robots())
            return

        elif path.endswith('.sws'):

            worksheet_filename = path[:-4]
            self.download_worksheet(worksheet_filename)
            return

        elif path[-12:] == '__main__.css':
            self.wfile.write(css.css())
            return

        elif path[-3:] == '.js':
            try: tempv = self._js_cache
            except: self._js_cache = {}

            binfile = None
            if self._js_cache.has_key(path):
                text, comp = self._js_cache[path]
                binfile, compressed = self.send_compressed(text, comp)
            else:
                text = None
                if self.path[-18:-7] == '__keyboard_':
                    text = keyboards.get_keyboard(self.path[-7:-5])
                elif path[-13:-3] == '__main__':
                    text = js.javascript()
                else: #if path[:7] == 'jsmath/' or path[:10] == 'highlight/':
                    try:
                        text = open(SAGE_EXTCODE + "/notebook/javascript/" + path).read()
                    except: pass

                if text is not None:
                    comp = gzip_compress(text)
                    self._js_cache[path] = (text, comp)
                    binfile, compressed = self.send_compressed(text, comp)

            if binfile is None:
                print 'file not found', path
                return self.file_not_found(path)

        else:
            try:
                if path in static_images: #this list is defined at the top of this file
                    binfile = self.image(path)
                elif path[:7] == 'jsmath/' or path[:10] == 'highlight/':
                    binfile = open(SAGE_EXTCODE + "/notebook/javascript/" + path, 'rb').read()
                else:
                    binfile = open(path, 'rb').read()
            except IOError, msg:
                print 'file not found', msg
                return self.file_not_found(path)

        self.send_response(200)

        mime_type = mimetypes.guess_type(self.path)[0]
        if mime_type is None:
            mime_type = "text/plain"
        self.send_header("Content-type", mime_type)
        if compressed:
            self.send_header("Content-Encoding", 'x-gzip')
        self.send_header("Cache-control", "no-store")

        self.end_headers()

        f = StringIO(binfile)
        f.write(binfile)
        f.flush()
        f.seek(0)


        alarm(3)
        try:
           while 1:
               buf = f.read(128)
               if not buf:
                   break
               self.wfile.write(buf)
        except KeyboardInterrupt:
           pass
        cancel_alarm()
        return f

        # the code below should work the same as above, but locks.
        alarm(3)
        try:
           shutil.copyfileobj(f, self.wfile, length=128)
        except KeyboardInterrupt:
            pass
        cancel_alarm()
        f.close()
        return f


    def show_page(self, worksheet_id,show_debug=False):
        self.send_head()
        if worksheet_id == '':
            W = None
            worksheet_id = None
            auth = True
        else:
            try:
                W = notebook.get_worksheet_with_id(worksheet_id)
            except KeyError:
                W = notebook.create_new_worksheet(worksheet_id)
                worksheet_id = W.id()
            auth = self.auth_worksheet(W)

        self.wfile.write(notebook.html(worksheet_id=worksheet_id,
                                       authorized=self.authorize(),
                                       show_debug=show_debug,
                                       worksheet_authorized = auth))


    def file_not_found(self, filename):
        self.send_response(404)
        self.send_header("Content-type", 'text/plain')
        self.end_headers()
        self.wfile.write("SAGE Server: File '%s' not found"%filename)

    def do_GET(self):
        """
        Handle a HTTP GET request.

        The basic operation of this function is to
        examine the request path and decide to serve
        either a static file based on the file extension,
        or in the case that a '?' found in the path,
        modify the path and continue.

        """
        if notebook.log_server():
            notebook.server_log().append(["GET", self.path])
        self.get_cookie()

        # Catch any doc_browser requests here
        if self.path.startswith('/doc_browser'):
            if self.path.endswith('.html'):
                return self.doc_browser(self.path)

        if self.path.startswith('/src_browser'):
            return self.src_browser(self.path)

        # The question mark trick here is so that images will be reloaded when
        # the async request requests the output text for a computation.
        # This is inspired by http://www.irt.org/script/416.htm/.
        show_debug=False
        i = self.path.rfind('?')
        if i != -1:
            if self.path[i+1:i+6] == 'debug':
                show_debug = True
            elif self.path[i+1:i+5] == 'edit':
                j = self.path.rfind('/')
                worksheet_filename = self.path[j+1:i]
                self.edit_text(worksheet_filename,prompts=False)
                return
            self.path = self.path[:i]



        if has_valid_file_extension(self.path) or \
               ('/jsmath/' in self.path and self.path.endswith('.js')):
            try:
                return self.get_file()
            except:
                print "Cancelled getting %s"%self.path
                return
        path = self.path.strip('/')
        i = path.find('/')
        if i == -1:
            i = len(path)

        try:
            worksheet_id = path[:i]
        except ValueError:
            worksheet_id = None
        path = path[i+1:]
        if path == '':
            return self.show_page(worksheet_id=worksheet_id,show_debug=show_debug)
        else:
            self.file_not_found(path)

    def get_cookie(self):
        """
        Get a cookie from the html header.
        """
        self.cookie=Cookie.SimpleCookie()
        if self.headers.has_key('cookie'):
            try:
                self.cookie=Cookie.SimpleCookie(self.headers.getheader("cookie"))
            except Cookie.CookieError, msg:
                print msg
                pass

    def authorize(self):
        """
        Set the permissions for a given worksheet.

        If the password is not given correctly, then
        the user is still allowed to view the given
        worksheet, just not evaluate any cells.

        This is not a rigorous security method,
        including the fact that the password is sent
        in clear text.

        """
        username = password = "";
        if self.cookie.has_key('username'):
            username = self.cookie['username'].value
        if self.cookie.has_key('password'):
            password = self.cookie['password'].value
        return notebook.authorize(username + ":" + password);

    def do_POST(self):
        """
        Handle a HTTP POST request.

        This differers from the HTTP GET by the fact that
        along with the POST request comes post variables
        which a used to evaluate the POST handling functions.

        """
        self.get_cookie()
        content_type, post_dict = cgi.parse_header(self.headers.getheader('content-type'))
        #verbose("POST: %s"%post_dict)
        self.save_notebook_every_so_often()

        if not self.authorize():
            self.body = {}

        elif content_type == 'multipart/form-data':
            M = cgi.parse_multipart(self.rfile, post_dict);

            if self.path[-5:] == '?edit' and self.path != '?edit':
                filename = self.path[1:-5]
                if M.has_key('button_save'):
                    self.edit_save(filename, M['textfield'][0])
                elif M.has_key('button_cancel'):
                    self.edit_cancel(filename)
                return

            if self.path == '/upload_worksheet' and M.has_key('fileField'):
                tmp = '%s/tmp.sws'%notebook.directory()
                f = file(tmp,'wb')
                f.write(M['fileField'][0])
                f.close()
                self.import_worksheet_local_file(tmp)
                os.unlink(tmp);

        elif content_type == 'application/x-www-form-urlencoded':
            self.send_response(200)
            self.send_header("Content-type", 'text/plain')
            self.end_headers()

            method = self.path[self.path.rfind('/')+1:]
            if method in ['cell_output_set', 'hide_all', 'restart_sage', 'show_all', 'introspect',
                          'new_cell', 'new_cell_after', 'delete_cell', 'cell_update', 'interrupt',
                          'cell_id_list', 'add_worksheet', 'delete_worksheet', 'unlock_worksheet',
                          'insert_wiki_cells', 'delete_cell_all', 'get_queue', 'set_cell_input']:
                eval("self.%s()"%method)
            else:
                if self.path[-8:]   == '/refresh':
                    self.show_page(worksheet_id=None, body_only=True)
                elif self.path[-6:] == '/eval0':
                    self.eval_cell(newcell=False)
                elif self.path[-6:] == '/eval1':
                    self.eval_cell(newcell=True)

        else:
            self.body = {}                   # Unknown content-type

        # some browsers send 2 more bytes...  ? which ?
        ready_to_read,x,y = select.select([self.connection],[],[],0)
        if ready_to_read:
            self.rfile.read(2)

    def do_HEAD(self):
        self.send_head()

    def send_head(self):
        self.send_response(200)
        if self.path[-4:] == '.png':
            self.send_header("Content-type", 'image/png')
        if self.path[-4:] == '.bmp':
            self.send_header("Content-type", 'image/bmp')
        elif self.path[-4:] == '.svg':
            self.send_header("Content-type", 'image/svg+xml')
        elif self.path[-4:] == '.txt':
            self.send_header("Content-type", 'text/plain')
        elif self.path[-2:] == '.c':
            self.send_header("Content-type", 'text/plain')
        elif self.path[-5:] == '.sobj':
            self.send_header("Content-type", 'application/sobj')
        else:
            self.send_header("Content-type", 'text/html')
        self.end_headers()

    def send_compressed(self, text, compressed = None):
        accept = self.headers.getheader("Accept-Encoding")
        if accept.find("gzip") >= 0:
            if compressed is None:
                return gzip_compress(text, 6), True
            else:
                return compressed, True
        else:
            return text, False

    def image(self, filename):
        try:
            return self._images[filename]
        except AttributeError:
            self._images = {}
            return self.image(filename)
        except KeyError:
	    filename = safe_path(filename)
            self._images[filename] = open(SAGE_EXTCODE +
                  '/notebook/images/' + filename, 'rb').read()
            return self._images[filename]

    def robots(self):
        """
        Return the robots.txt file contents (as a string) that should
        be used when search engines hit this site.

        The default is to not allow any robots.  This can be
        customized by creating a file robots.txt inside the
        sage_notebook directory.
        """
        try:
            return self._robots
        except AttributeError:
            pass
        filename = '%s/robots.txt'%notebook.directory()
        if os.path.exists(filename):
            robots_txt = open(filename).read()
        else:
            robots_txt = """
User-agent: *
Disallow: /
            """
            open(filename,'w').write(robots_txt)
        self._robots = robots_txt
        return robots_txt


class NotebookServer:
    def __init__(self, notebook, port, address):
        self.__notebook = notebook
        self.__httpd = BaseHTTPServer.HTTPServer((address,int(port)), WebServer)
        try:
            # this seems to be BSD specific
            self.__httpd.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEPORT, 1)
        except AttributeError:
            try:
                # on linux, this seems to do both addr and port.
                self.__httpd.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            except AttributeError:
                pass
        self.__address = address
        self.__port = port

    def auth_string(self):
        try:
            return self.__auth
        except AttributeError:
            self.__auth = ":"
        return self.__auth

    def address(self):
        return self.__address

    def port(self):
        return self.__port

    def notebook(self):
        return self.__notebook

    def url(self):
        return 'http://%s:%s"'%(self.__address, self.__port)

    def serve(self):
        global notebook
        notebook = self.notebook()
        while True:
            try:
                print "Press Control-C *TWICE* to stop the server."
                self.__httpd.serve_forever()
            except KeyboardInterrupt:
                notebook.save()
                print "Shutting down notebook server."
                notebook = None
                break
            else:
                notebook.save()

def gzip_compress(data, compresslevel=9):
    io = StringIO()
    gz = GzipFile(fileobj = io, mode='w', compresslevel=compresslevel)
    gz.write(data)
    gz.close()
    io.seek(0)
    return io.read()




VALID_EXTENSIONS = ['eps', 'pdf', 'png', 'bmp', 'svg', 'tex',
                    'dvi', 'log', 'css',
                    'txt', 'ico', 'sws',
                    'c', 'sobj', 'html',
                    'ps', 'js', 'hg', 'patch']

def has_valid_file_extension(path):
    for ext in VALID_EXTENSIONS:
        if path.endswith('.' + ext):
            return True
    return False
