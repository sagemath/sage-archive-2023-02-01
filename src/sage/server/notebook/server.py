"""
Web Server Component of SAGE Notebook
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


# Python libraries
import BaseHTTPServer
import cgi
import os, sys
import select
from   StringIO import StringIO
import shutil

import css, js

# SAGE libraries
import sage.interfaces.sage0

from   sage.misc.misc import verbose, word_wrap
import sage.misc.preparser
from   sage.misc.viewer     import BROWSER
from   sage.ext.sage_object import load, SageObject

SEP = '___S_A_G_E___'

# The unique notebook that is being served by the web server.
# This is a global variable that is set by the NotebookServer
# constructor.  It is equal to the last notebook that the
# serve method was called on, and is set back to None
# after the serve method terminates.

notebook = None


class WebServer(BaseHTTPServer.BaseHTTPRequestHandler):
    def eval_cell(self, time=False, completions=False):
        length = int(self.headers.getheader('content-length'))
        qs = self.rfile.read(length)
        C = cgi.parse_qs(qs, keep_blank_values=1)
        input_text = C['input'][0]
        id = int(C['id'][0])
        input_text = input_text.replace('__plus__','+')
        print id, input_text
        W = notebook.current_workbook()
        cell = W.get_cell_with_id(id)
        cell.set_input_text(input_text)
        notebook.save()

        cell.evaluate(time, completions)
        self.send_response(200)
        self.send_header("Content-type", 'text/plain')
        self.end_headers()

        if cell.is_last():
            new_cell = W.append_new_cell()
            self.wfile.write(str(new_cell.id()) + SEP +
                             new_cell.html(div_wrap=False) + SEP + str(id))
        else:
            self.wfile.write(str(cell.next_id()) + SEP +
                             'no_new_cell' + SEP + str(id))

    def new_cell(self):
        length = int(self.headers.getheader('content-length'))
        qs = self.rfile.read(length)
        C = cgi.parse_qs(qs, keep_blank_values=1)
        id = int(C['id'][0])
        print "Adding new cell before cell with id %s"%id
        W = notebook.current_workbook()
        C = W.new_cell_before(id)
        notebook.save()
        self.send_response(200)
        self.send_header("Content-type", 'text/plain')
        self.end_headers()
        self.wfile.write(str(C.id()) + SEP + C.html(div_wrap=False) + SEP + str(id))

    def delete_cell(self):
        length = int(self.headers.getheader('content-length'))
        qs = self.rfile.read(length)
        C = cgi.parse_qs(qs, keep_blank_values=1)
        id = int(C['id'][0])
        print "Possibly deleting cell with id %s"%id
        self.send_response(200)
        self.send_header("Content-type", 'text/plain')
        self.end_headers()

        W = notebook.current_workbook()
        if len(W) <= 1 or W.is_last_id_and_previous_is_nonempty(id):
            self.wfile.write('ignore')
        else:
            next_id = W.delete_cell_with_id(id)
            notebook.save()
            self.wfile.write('delete' + SEP + str(id) + SEP + str(next_id))


    def update_cells(self):
        workbook = notebook.current_workbook()
        cols = notebook.defaults()['word_wrap_cols']
        status, cell = workbook.check_comp()
        if status == 'd':
            notebook.save()
            variables = workbook.variables_html()
            objects = notebook.objects_html()
        else:
            variables = '...'  # not used
            objects = "..." # note used
        if cell is None:
            msg = 'empty'
        else:
            msg = '%s%s %s'%(status, cell.id(),
                              SEP.join([cell.output_text(),
                                        cell.output_text(cols),
                                        variables,
                                        objects]))
            # more comps to go.
            workbook.start_next_comp()
        self.wfile.write(msg)

    def get_file(self):
        try:
            binfile = open(self.path[1:], 'rb').read()
        except IOError:
            return self.file_not_found()
        self.send_response(200)
        if self.path[-4:] == '.png':
            self.send_header("Content-type", 'image/png')
        elif self.path[-3:] == '.ps':
            self.send_header("Content-type", 'application/postscript')
        elif self.path[-4:] == '.eps':
            print 'eps!!'
            self.send_header("Content-type", 'image/x-eps')
        elif self.path[-4:] == '.svg':
            self.send_header("Content-type", 'image/svg+xml')
        elif self.path[-4:] == '.txt':
            self.send_header("Content-type", 'text/plain')
        elif self.path[-5:] == '.sobj':
            self.send_header("Content-type", 'application/sobj')
        self.end_headers()
        f = StringIO()
        f.write(binfile)
        f.seek(0)
        shutil.copyfileobj(f, self.wfile)
        f.close()
        return f

    def show_page(self):
        self.send_head()
        self.wfile.write(notebook.html())

    def file_not_found(self):
        self.send_response(404)
        self.send_header("Content-type", 'text/plain')
        self.end_headers()
        self.wfile.write("SAGE Server: File not found")

    def do_GET(self):
        print "GET: " + self.path

        # The question mark hack here is so that images will be reloaded when
        # the async request requests the output text for a computation.
        # This is a total hack, inspired by http://www.irt.org/script/416.htm/.
        i = self.path.rfind('?')
        if i != -1:
            self.path = self.path[:i]

        print self.path

        if self.path[-13:] == '/update_cells':
            self.update_cells()

        elif self.path[-10:] == '/interrupt':
            self.send_head()
            t = notebook.interrupt()
            if t:
                self.wfile.write('ok')
            else:
                self.wfile.write('restart')
            return

        elif self.path[-13:] == '/cell_id_list':
            L = notebook.current_workbook().cell_id_list()
            self.send_head()
            self.wfile.write(' '.join(str(x) for x in L))
            return

        elif self.path[-4:] in ['.eps', '.png', '.svg', '.txt'] or \
             self.path[-5:] == '.sobj' or self.path[-3:] == '.ps':

            return self.get_file()

        else:
            self.show_page()

    def do_POST(self):
        content_type, post_dict = cgi.parse_header(self.headers.getheader('content-type'))
        print "POST: %s"%post_dict

        if content_type == 'multipart/form-data':
            self.body = cgi.parse_multipart(self.rfile, post_dict)

        elif content_type == 'application/x-www-form-urlencoded':
            if self.path[-6:] == '/eval0':
                self.eval_cell(time=False)
            elif self.path[-6:] == '/eval1':
                self.eval_cell(time=True)
            elif self.path[-6:] == '/eval2':
                self.eval_cell(completions=True)
            elif self.path[-9:]  == '/new_cell':
                self.new_cell()
            elif self.path[-12:] == '/delete_cell':
                self.delete_cell()
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
        elif self.path[-4:] == '.svg':
            self.send_header("Content-type", 'image/svg+xml')
        elif self.path[-4:] == '.txt':
            self.send_header("Content-type", 'text/plain')
        elif self.path[-5:] == '.sobj':
            self.send_header("Content-type", 'application/sobj')
        else:
            self.send_header("Content-type", 'text/html')
        self.end_headers()




class NotebookServer:
    def __init__(self, notebook, port, address):
        self.__notebook = notebook
        self.__httpd = BaseHTTPServer.HTTPServer((address,int(port)), WebServer)
        self.__address = address
        self.__port = port

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
        try:
            print "Press Control-C to stop the server."
            self.__httpd.serve_forever()
        except KeyboardInterrupt:
            print "Shutting down notebook server."
            notebook = None



