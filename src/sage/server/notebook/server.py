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
import Cookie
import cPickle
import base64
import css, js

# SAGE libraries
import sage.interfaces.sage0

from   sage.misc.misc import verbose, word_wrap
import sage.misc.preparser
from   sage.ext.sage_object import load, SageObject

SEP = '___S_A_G_E___'

# The unique notebook that is being served by the web server.
# This is a global variable that is set by the NotebookServer
# constructor.  It is equal to the last notebook that the
# serve method was called on, and is set back to None
# after the serve method terminates.

notebook = None


class WebServer(BaseHTTPServer.BaseHTTPRequestHandler):
    def get_postvars(self):
        length = int(self.headers.getheader('content-length'))
        qs = self.rfile.read(length)
        return cgi.parse_qs(qs, keep_blank_values=1)


    def cell_output_set(self):
        C = self.get_postvars()
        id = int(C['id'][0])
        typ = C['type'][0]
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        cell = W.get_cell_with_id(id)
        cell.set_cell_output_type(typ)

    def hide_all(self):
        C = self.get_postvars()
        id = int(C['worksheet_id'][0])
        W = notebook.get_worksheet_with_id(id)
        W.hide_all()

    def restart_sage(self):
        C = self.get_postvars()
        id = int(C['worksheet_id'][0])
        W = notebook.get_worksheet_with_id(id)
        W.restart_sage()
        self.wfile.write('done')

    def show_all(self):
        C = self.get_postvars()
        id = int(C['worksheet_id'][0])
        W = notebook.get_worksheet_with_id(id)
        W.show_all()

    def eval_cell(self, newcell=False, introspect=False):
        C = self.get_postvars()
        id = int(C['id'][0])
        input_text = C['input'][0]
        input_text = input_text.replace('__plus__','+')
        notebook.add_to_history(input_text)
        verbose('%s: %s'%(id, input_text))
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        cell = W.get_cell_with_id(id)
        cell.set_input_text(input_text)
        notebook.save()

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
            self.wfile.write(str(cell.next_id()) + SEP +
                             'no_new_cell' + SEP + str(id))

    def introspect(self):
        C = self.get_postvars()
        id = int(C['id'][0])
        before_cursor = C['before_cursor'][0].replace('__plus__','+')
        after_cursor = C['after_cursor'][0].replace('__plus__','+')
        input_text = (before_cursor+after_cursor)
        verbose('introspect -- %s: %s|%s'%(id, before_cursor, after_cursor))

        W = notebook.get_worksheet_that_has_cell_with_id(id)
        cell = W.get_cell_with_id(id)
        cell.set_input_text(before_cursor + after_cursor)
        cell.evaluate(introspect=[before_cursor, after_cursor])

        self.wfile.write(str(cell.next_id()) + SEP +
                         'no_new_cell' + SEP + str(id))


    def new_cell(self):
        C = self.get_postvars()
        id = int(C['id'][0])
        verbose("Adding new cell before cell with id %s"%id)
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        cell = W.new_cell_before(id)
        notebook.save()
        self.wfile.write(str(cell.id()) + SEP + cell.html(div_wrap=False) + SEP + \
                         str(id))

    def new_cell_after(self):
        C = self.get_postvars()
        id = int(C['id'][0])
        verbose("Adding new cell after cell with id %s"%id)
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        cell = W.new_cell_after(id)
        notebook.save()
        self.wfile.write(str(cell.id()) + SEP + cell.html(div_wrap=False) + SEP + \
                         str(id) + SEP)

    def delete_cell(self):
        C = self.get_postvars()
        id = int(C['id'][0])
        verbose("Deleting cell with id %s"%id)
        W = notebook.get_worksheet_that_has_cell_with_id(id)
        if len(W) <= 1 or W.is_last_id_and_previous_is_nonempty(id):
            self.wfile.write('ignore')
        else:
            prev_id = W.delete_cell_with_id(id)
            notebook.save()
            self.wfile.write('delete' + SEP + str(id) + SEP + str(prev_id) + SEP + str(W.cell_id_list()))

    def cell_update(self):
        C = self.get_postvars()
        worksheet_id = int(C['worksheet_id'][0])
        cell_id = int(C['cell_id'][0])

        worksheet = notebook.get_worksheet_with_id(worksheet_id)
        cols = notebook.defaults()['word_wrap_cols']

        # update the computation one step.
        worksheet.check_comp()
        # now get latest status on our cell
        status, cell = worksheet.check_cell(cell_id)

        #print status, cell   # debug
        if status == 'd':
            try:
                notebook.save()
            except:
                print "WARNING -- failure to pickle the notebook"
            variables = worksheet.variables_html()
            objects = notebook.object_list_html()
            attached_files = worksheet.attached_html()
        else:
            variables = '...'  # not used
            objects = "..." # note used
            attached_files = '...' # not used
        if status == 'd':
            new_input = cell.changed_input_text()
            out_html = cell.output_html()
        else:
            new_input = ''
            out_html = ''
        msg = '%s%s %s'%(status, cell.id(),
                          SEP.join([cell.output_text().replace('<','&lt;'),
                                    cell.output_text(cols).replace('<','&lt;'),
                                    out_html,
                                    new_input,
                                    variables,
                                    objects,
                                    attached_files]))

        # more comps to go ?
        worksheet.start_next_comp()
        self.wfile.write(msg)

##     def update_cells(self):
##         C = self.get_postvars()
##         worksheet_id = int(C['worksheet_id'][0])
##         worksheet = notebook.get_worksheet_with_id(worksheet_id)
##         cols = notebook.defaults()['word_wrap_cols']
##         status, cell = worksheet.check_comp()
##         print status, cell   # debug
##         if status == 'd':
##             try:
##                 notebook.save()
##             except:
##                 print "WARNING -- failure to pickle the notebook"
##             variables = worksheet.variables_html()
##             objects = notebook.object_list_html()
##             attached_files = worksheet.attached_html()
##         else:
##             variables = '...'  # not used
##             objects = "..." # note used
##             attached_files = '...' # not used
##         if cell is None or status == 'e':
##             msg = 'empty'
##         else:
##             if status == 'd':
##                 new_input = cell.changed_input_text()
##                 out_html = cell.output_html()
##             else:
##                 new_input = ''
##                 out_html = ''
##             msg = '%s%s %s'%(status, cell.id(),
##                               SEP.join([cell.output_text().replace('<','&lt;'),
##                                         cell.output_text(cols).replace('<','&lt;'),
##                                         out_html,
##                                         new_input,
##                                         variables,
##                                         objects,
##                                         attached_files]))
##             # more comps to go.
##             worksheet.start_next_comp()
##         self.wfile.write(msg)

    def interrupt(self):
        C = self.get_postvars()
        worksheet_id = int(C['worksheet_id'][0])
        worksheet = notebook.get_worksheet_with_id(worksheet_id)
        t = worksheet.interrupt()
        if t:
            self.wfile.write('ok')
        else:
            self.wfile.write('restart')

    def add_worksheet(self):
        C = self.get_postvars()
        worksheet_name = C['name'][0]
        try:
            W = notebook.create_new_worksheet(worksheet_name)
        except ValueError, msg:
            print msg
            self.wfile.write(msg)
            return
        new_worksheet_list = notebook.worksheet_list_html(W.name())
        self.wfile.write(new_worksheet_list + SEP + str(W.name()))

    def delete_worksheet(self):
        C = self.get_postvars()
        worksheet_name = C['name'][0]
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
        self.send_header("Location", '/%d'%W.id())
        self.end_headers()

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
        s += '<script language=javascript>' + js.javascript() + '</script>\n'
        s += '<style>' + css.css() + '</style>\n'
        s += '</head>\n'
        s += '<body>\n'
        s += W.html(include_title=False, do_print=do_print)
        #if do_print:
        #    s += '<script language=javascript>window.print()</script>\n'
        s += '\n</body>\n'
        self.wfile.write(s)

    def input_history_text(self):
        self.send_head()
        t = notebook.history_text()
        t = t.replace('<','&gt;')
        s = '<head>\n'
        s += '<title>SAGE Input History</title>\n'
        s += '</head>\n'
        s += '<body>\n'
        s += '<pre>' + t + '</pre>\n'
        s += '<a name="bottom"></a>\n'
        s += '<script language=javascript> window.location="#bottom"</script>\n'
        s += '</body>\n'
        self.wfile.write(s)

    def help_window(self):
        self.send_head()
        self.wfile.write(notebook.help_window())

    def upload_window(self):
        self.send_head()
        self.wfile.write(notebook.upload_window())

    def download_worksheet(self, filename):
        try:
            notebook.export_worksheet(filename, filename)
        except KeyError:
            self.file_not_found()
            return
        self.send_response(200)
        self.send_header("Content-type", 'application/sage')
        self.end_headers()
        binfile = open('%s/%s.sws'%(notebook.directory(), filename), 'rb').read()
        f = StringIO()
        f.write(binfile)
        f.seek(0)
        shutil.copyfileobj(f, self.wfile)
        f.close()

    def cell_id_list(self):
        C = self.get_postvars()
        worksheet_id = int(C['worksheet_id'][0])
        worksheet = notebook.get_worksheet_with_id(worksheet_id)
        L = worksheet.cell_id_list()
        s = ' '.join(str(x) for x in L)
        self.wfile.write(s)

    def get_file(self):
        path = self.path
        if path[-5:] == '.sobj':
            path = '%s/%s'%(os.path.abspath(notebook.object_directory()), path)
        else:
            path = path[1:]
        if path[-5:] == '.html' and not '/' in path:
            worksheet_filename = path[:-5]
            if worksheet_filename == '__history__':
                self.input_history_text()
            elif worksheet_filename == '__help__':
                self.help_window()
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

        elif path[-4:] == '.sws':

            worksheet_filename = path[:-4]
            self.download_worksheet(worksheet_filename)
            return

        try:
            if path[-11:] == 'favicon.ico':
                binfile = self.favicon()
            else:
                binfile = open(path, 'rb').read()
        except IOError, msg:
            print 'file not found', msg
            return self.file_not_found()
        self.send_response(200)
        if self.path[-4:] == '.png':
            self.send_header("Content-type", 'image/png')
        elif self.path[-3:] == '.ps':
            self.send_header("Content-type", 'application/postscript')
        elif self.path[-4:] == '.tex':
            self.send_header("Content-type", 'application/latex')
        elif self.path[-4:] == '.dvi':
            self.send_header("Content-type", 'application/x-dvi')
        elif self.path[-4:] == '.log':
            self.send_header("Content-type", 'text/plain')
        elif self.path[-4:] == '.eps':
            self.send_header("Content-type", 'image/x-eps')
        elif self.path[-4:] == '.ico':
            self.send_header("Content-type", 'image/x-icon')
        elif self.path[-4:] == '.svg':
            self.send_header("Content-type", 'image/svg+xml')
        elif self.path[-4:] == '.txt':
            self.send_header("Content-type", 'text/plain')
        elif self.path[-5:] == '.sobj':
            self.send_header("Content-type", 'application/sage')
        elif self.path[-5:] == '.html':
            self.send_header("Content-type", 'text/html')
        self.end_headers()
        f = StringIO()
        f.write(binfile)
        f.seek(0)
        shutil.copyfileobj(f, self.wfile)
        f.close()
        return f

    def show_page(self, worksheet_id=None):
        self.send_head()
        self.wfile.write(notebook.html(worksheet_id=worksheet_id,
                                       authorized=self.authorize()))


    def file_not_found(self):
        self.send_response(404)
        self.send_header("Content-type", 'text/plain')
        self.end_headers()
        self.wfile.write("SAGE Server: File not found")

    def do_GET(self):
        verbose("GET: " + self.path)

        # The question mark hack here is so that images will be reloaded when
        # the async request requests the output text for a computation.
        # This is a total hack, inspired by http://www.irt.org/script/416.htm/.
        i = self.path.rfind('?')
        if i != -1:
            self.path = self.path[:i]

        verbose(self.path)

        if self.path[-4:] in ['.eps', '.png', '.svg', '.tex', '.dvi', '.log', \
                              '.txt', '.ico', '.sws'] or \
               self.path[-5:] in ['.sobj', '.html'] or self.path[-3:] == '.ps':
            return self.get_file()

        path = self.path.strip('/')
        i = path.find('/')
        if i == -1:
            i = len(path)

        try:
            worksheet_id = int(path[:i])
        except ValueError:
            worksheet_id = None
        path = path[i+1:]
        if path == '':
            return self.show_page(worksheet_id=worksheet_id)
        else:
            self.file_not_found()

    def authorize(self):
        self.cookie=Cookie.SimpleCookie()
        if self.headers.has_key('cookie'):
            self.cookie=Cookie.SimpleCookie(self.headers.getheader("cookie"))
        username = password = "";
        if self.cookie.has_key('username'):
            username = self.cookie['username'].value
        if self.cookie.has_key('password'):
            password = self.cookie['password'].value
        return notebook.authorize(username + ":" + password);

    def do_POST(self):
        content_type, post_dict = cgi.parse_header(self.headers.getheader('content-type'))
        verbose("POST: %s"%post_dict)

        if not self.authorize():
            self.body = {}

        elif content_type == 'multipart/form-data':
            M = cgi.parse_multipart(self.rfile, post_dict);

            if self.path == '/upload_worksheet' and M.has_key('fileField'):
                f = file('temp.sws','wb')
                f.write(M['fileField'][0])
                f.close()
                self.import_worksheet_local_file('temp.sws')
                os.unlink('temp.sws');

        elif content_type == 'application/x-www-form-urlencoded':
            self.send_response(200)
            self.send_header("Content-type", 'text/plain')
            self.end_headers()
            if self.path[-8:]   == '/refresh':
                self.show_page(worksheet_id=None, body_only=True)
            elif self.path[-6:] == '/eval0':
                self.eval_cell(newcell=False)
            elif self.path[-6:] == '/eval1':
                self.eval_cell(newcell=True)
            elif self.path[-16:] == '/cell_output_set':
                self.cell_output_set()
            elif self.path[-9:]  == '/hide_all':
                self.hide_all()
            elif self.path[-13:] == '/restart_sage':
                self.restart_sage()
            elif self.path[-9:]  == '/show_all':
                self.show_all()
            elif self.path[-11:] == '/introspect':
                self.introspect()
            elif self.path[-9:]  == '/new_cell':
                self.new_cell()
            elif self.path[-15:] == '/new_cell_after':
                self.new_cell_after()
            elif self.path[-12:] == '/delete_cell':
                self.delete_cell()
            elif self.path[-12:] == '/cell_update':
                self.cell_update()
            #elif self.path[-13:] == '/update_cells':
            #    self.update_cells()
            elif self.path[-10:] == '/interrupt':
                self.interrupt()
            elif self.path[-13:] == '/cell_id_list':
                self.cell_id_list()
            elif self.path[-14:] == '/add_worksheet':
                self.add_worksheet()
            elif self.path[-17:] == '/delete_worksheet':
                self.delete_worksheet()
#            elif self.path == '/upload_worksheet':
#                self.upload_worksheet_local_file()
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



    def favicon(self):
        s = """
            AAABAAEAEBAAAAEACABoBQAAFgAAACgAAAAQAAAAIAAAAAEACAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAABDQ0MARUVFAEpKSgBLS0sATU1NAFRUVABWVlYAV1dXAFlZWQBaWloAM4IyAGRkZABlZWUA
            bm5uAG9vbwBzc3MAAcQAAHV1dQB3eHcAeXl5AHd8dwB+fn4AhYWFAIeHhwCSkpIAmJiYAJmZmQCe
            np4AoKCgAKWlpQCmpqYAp6enALa2tgC4uLgAubm5AL6+vgDAwMAAwsLCAMPDwwDGxsYA0dHRANLS
            0gDa2toA3d3dAOTk5ADm5uYA6urqAPDw8ADz8/MA9PT0APX19QD29vYA9/f3APj4+AD5+fkA+vr6
            APv7+wD8/PwA/f39AP7+/gD///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAOzsqLDs7Ozs7Ozs7Ozs7OzsQEB4jOzs7Ozs7Ozs7Ozs7Oig0CxAnJCU5Ozs7Ozs7Ozs7EBAR
            EBAKGTg6OTs7Ozs7OxAMNQkgEAUuKx0mOzs7Ozs6OBAPEBQEKRAQBxs7Ozs7OzsQFzIQABAIIRAF
            Hzs7Ozs7EA00MSIQDjEyEBo7Ozs7OzgQFTIQEAMlOTQuOzs7Ozs7LyYtEBIQBSE5Ozs7Ozs7Ozs4
            EAIkMRAWMzs7Ozs7Ozs7OBABHC82Ljs7Ozs7Ozs7OzsQEAYTJTk7Ozs7Ozs7Ozs7OBAQEBg3Ozs7
            PDw7Ozs7Ozs5NTA0Ozs7Ozs8PDs7Ozs7Ozs7Ozs7OwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA=
            """
        return base64.decodestring(s)



class NotebookServer:
    def __init__(self, notebook, port, address):
        self.__notebook = notebook
        self.__httpd = BaseHTTPServer.HTTPServer((address,int(port)), WebServer)
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
        try:
            print "Press Control-C to stop the server."
            self.__httpd.serve_forever()
        except KeyboardInterrupt:
            print "Shutting down notebook server."
            notebook = None



