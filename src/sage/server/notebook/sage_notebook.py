"""
Web Server Component of SAGE Notebook Interface
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
import shutil
import socket
from   StringIO import StringIO

import css, js

# SAGE libraries
import sage.interfaces.sage0

from   sage.misc.misc import word_wrap
import sage.misc.preparser
from   sage.misc.viewer     import BROWSER
from   sage.ext.sage_object import load, SageObject


def get_doc(query):
    query = query.replace('\n','').strip()
    if query[:9] == '?__last__':
        return get_docstring_last(int(query[9:])/15)
    if len(query) > 1:
        if query[:2] == '??':
            return get_source_code(query[2:])
        elif query[-2:] == '??':
            return get_source_code(query[:-2])
    if len(query) > 0:
        if query[0] == '?':
            return get_docstring(query[1:])
        elif query[-1] == '?':
            return get_docstring(query[:-1])
    return get_completions(query)

def get_docstring(query):
    cmd = '_support_.docstring("%s", globals())'%query
    z = sage0.eval(cmd)
    z = z.replace('\\n','\n').replace('\\t','        ')[1:-1]
    z = word_wrap(z, ncols=numcols)
    return z

def get_source_code(query):
    cmd = '"".join(_support_.source_code("%s", globals()))'%query
    z = sage0.eval(cmd)
    z = z.replace('\\n','\n').replace("\\","").replace('\\t','        ')[1:-1]
    return z

def get_completions(query):
    cmd = '"<br>".join(_support_.completions("%s", globals()))'%query
    z = sage0.eval(cmd)
    _last_ = z
    return z[1:-1]


class HTML_Interface(BaseHTTPServer.BaseHTTPRequestHandler):
    def __files(self,number):
        global directory
        return os.listdir("%s/cells/%s"%(directory,number))

    def __show_page(self, number):
        global current_workbook
        f = self.send_head()
        if f:
            f = StringIO()
            f.write('<html><head><title>SAGE: %s</title></head>\n\n'%save_name)
            f.write('<style>' + css.css() + '</style>\n\n')
            f.write('<script language=javascript>' + js.javascript(number) + '</script>\n\n')
            f.write("""
            <body onload="javascript:scroll_to_bottom()">

            <div class="search" id="search">
                <input type="text" id="search_input" class="search" onKeyUp="search_box();">

                <span id="search_doc_topbar" class="search_doc_topbar">
                    <table bgcolor="73a6ff" width="100%">
                        <tr>
                            <td align=center class="menubar">
                                Completions
                            </td>
                        </tr>
                    </table>
                </span>
                <span id="search_doc" class="search_doc">
                </span>
            </div>

            <span class="cell_collection">
                <table width="100%" height="80%"><tr><td align=center>
                <h2 class="banner">
                <a class="banner" href="http://modular.math.washington.edu/sage">SAGE:
                       Software for Algebra and Geometry Experimentation</a>
                </h2>
            """)
            #f.write('<a href="logfile.txt">Log File</a>')
            if len(current_workbook) == 0 or current_workbook[-1].cmd != '':
                I = Cell()
                current_workbook.append(I)
            current_workbook.save(workbook_sobj_file)
            f.write(current_workbook.html(numcols))
            f.write('</td></tr></table>\n')
            f.write('</span></body></html>\n')
            f.seek(0)
            self.copyfile(f, self.wfile)

            f.seek(0)
            self.copyfile(f, open('%s/index.html'%directory,'w'))

            f.close()
            return f

    def do_GET(self):
        if self.path[-4:] in ['.png', '.svg', '.txt'] or self.path[-5:] == '.sobj':
            if os.path.exists("%s%s"%(directory,self.path)):
                f = self.send_head()
                if f:
                    binfile = open("%s%s"%(directory,self.path), 'rb').read()
                    f.write(binfile)
                    f.seek(0)
                    self.copyfile(f,self.wfile)
                    f.close()
                    return f
            else:
                self.file_not_found()

        elif self.path[-8:] == 'walltime':
            self.send_response(200)
            self.send_header("Content-type", 'text/plain')
            self.end_headers()
            f = StringIO()
            f.write(sage0._eval_line('int(walltime(_start_time_))'))
            f.seek(0)
            self.copyfile(f, self.wfile)
            f.close()
            return f

        else:
            self.__show_page(0)

    def do_POST(self):
        global current_workbook
        ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))
        length = int(self.headers.getheader('content-length'))
        print "POST: ", self.path
        if ctype == 'multipart/form-data':
            self.body = cgi.parse_multipart(self.rfile, pdict)
        elif ctype == 'application/x-www-form-urlencoded':
            if self.path == '/search':
                qs = self.rfile.read(length)
                C = cgi.parse_qs(qs, keep_blank_values=1)
                d = get_doc(C['query'][0])
                d = d.replace(' ','&nbsp;')
                self.send_response(200)
                self.send_header("Content-type", 'text/plain')
                self.end_headers()
                f = StringIO()
                f.write('%s'%d)
                f.seek(0)
                self.copyfile(f, self.wfile)
                f.close()
                return f
            else:
                qs = self.rfile.read(length)
                C = cgi.parse_qs(qs, keep_blank_values=1)
                number = eval(C.keys()[0])
                print "Workbook '%s': evaluating cell number %s"%(
                    save_name, number)
                current_dir = "%s/cells/%d"%(directory,number)
                code_to_eval = C[C.keys()[0]][0]
                open(log_file,'a').write('\n#### INPUT \n' + code_to_eval + '\n')
                try:
                    if number > len(current_workbook)-1:
                        current_workbook.set_last_cmd(code_to_eval)
                        number = len(current_workbook)-1
                    else:
                        # re-evaluating a code block
                        current_workbook[number].cmd = code_to_eval
                    #code_to_eval = code_to_eval.replace('\\','')
                    s = sage.misc.preparser.preparse_file(code_to_eval, magic=False,
                                                          do_time=True, ignore_prompts=True)
                    s = [x for x in s.split('\n') if len(x.split()) > 0 and \
                          x.lstrip()[0] != '#']   # remove all blank lines and comment lines
                    if len(s) > 0:
                        t = s[-1]
                        if len(t) > 0 and not ':' in t and \
                               not t[0].isspace() and not t[:3] == '"""':
                            t = t.replace("'","\\'")
                            s[-1] = "exec compile('%s', '', 'single')"%t

                    s = '\n'.join(s) + '\n'

                    open('%s/_temp_.py'%directory, 'w').write(s)

                    if not os.path.exists(current_dir):
                        os.makedirs(current_dir)

                    for F in os.listdir(current_dir):
                        os.unlink("%s/%s"%(current_dir,F))

                    sage0._eval_line('os.chdir("%s")'%current_dir)

                    try:
                        o = sage0._eval_line('execfile("%s/_temp_.py")'%directory)
                        if not 'Traceback (most recent call last):' in o:
                            open(log_file,'a').write('#### OUTPUT \n' + o + '\n')
                        else:
                            open(log_file,'a').write('#### OUTPUT (error)\n')
                    except KeyboardInterrupt, msg:
                        print "Keyboard interrupt!"
                        o = msg

                    o = word_wrap(o, ncols=numcols)

                    #while True:
                    #    print "waiting for output"
                    #    o = sage0._get()
                    #    if o is None:
                    #        print "output isn't ready yet"
                    #    else:
                    #        print "o = ", o
                    #        break
                    #    import time
                    #    time.sleep(0.5)

                    current_workbook[number].out = o
                    current_workbook[number].file_list = self.__files(number)
                    self.__show_page(number)

                except (RuntimeError, TypeError), msg:
                    print "ERROR!!!", msg
                    self.__show_page(0)

        else:
            self.body = {}                   # Unknown content-type

        # some browsers send 2 more bytes...
        [ready_to_read,x,y] = select.select([self.connection],[],[],0)

        if ready_to_read:
            self.rfile.read(2)


    def do_HEAD(self):
        f = self.send_head()
        if f:
            f.close()

    def file_not_found(self):
        self.send_response(404)
        self.send_header("Content-type", 'text/plain')
        self.end_headers()
        f = StringIO()
        f.write("File not found")
        f.seek(0)
        self.copyfile(f, self.wfile)
        f.close()
        return f

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
        f = StringIO()
        #print "URL Path: %s\n" % self.path
        f.seek(0)
        return f

    def copyfile(self, source, outputfile):
        shutil.copyfileobj(source, outputfile)

def sage_notebook(dir       ='sage_notebook',
                  port      = 8000,
                  address   = 'localhost',
                  ncols     = 80,
                  nrows     = 8,
                  viewer    = True,
                  workbook  = None,
                  max_tries = 10):
    r"""
    Start a SAGE http server at the given port.

    Typical usage:
        \code{server_http1('mysession')}

    Use it.  To start it later, just type \code{server_http1('mysession')}
    again.

    INPUT:
        dir -- (default: 'sage_server') name of the server directory; your
                session is saved in a directory with that name.  If you restart
                the server with that same name then it will restart
                in the state you left it, but with none of the
                variables defined (you have to re-eval blocks).
        port -- (default: 8000) port on computer where the server is served
        address -- (default: 'localhost') address that the server will listen on
        ncols -- (default: 90) default number of columns for input boxes
        nrows -- (default: 8) default number of rows for input boxes
        viewer -- bool (default:True); if True, pop up a web browser at the URL
        workbook    -- (default: None) resume from a previous workbook
        max_tries -- (default: 10) maximum number of ports > port to try in
                     case given port can't be opened.

    NOTES:

    When you type

      server_http1(...)

    you start a web server on the machine you type
    that command on.  You don't connect to another
    machine.  It's kind of like starting apache
    on a server.  So do this if you want to start
    a SAGE notebook accessible from anywher.

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
    very limited privileges (e.g., empty home directory), and
    running \code{server_http1} as that user.}
    """
    global directory, current_workbook, workbook_sobj_file, \
           numcols, numrows, sage0, save_name, log_file
    save_name = dir
    dir = '_'.join(dir.split())  # no white space
    directory = os.path.abspath(dir)

    if os.path.isfile(directory):
        raise RuntimeError, 'Please delete or move the file "%s" and run the server again'%directory

    if not os.path.exists(directory):
        print 'Creating directory "%s"'%directory
        os.makedirs(directory)

    log_file = '%s/logfile.txt'%directory
    open(log_file,'a').close()  # touch the file

    numcols = int(ncols)
    numrows = int(nrows)
    workbook_sobj_file = '%s/cells.sobj'%directory

    if workbook is None and os.path.exists(workbook_sobj_file):
        try:
            workbook = load(workbook_sobj_file)
        except IOError:
            print "Unable to load log %s (creating new log)"%workbook_sobj_file

    if workbook is None:
        current_workbook = Workbook()
    else:
        current_workbook = workbook
    sage0 = sage.interfaces.sage0.Sage()
    if not os.path.exists('%s/sobj'%directory):
        os.makedirs('%s/sobj'%directory)

    # Initialize the sage0 object for use with the web notebook interface.
    sage0.eval('import sage.server.support as _support_')
    sage0.eval('_support_.init("%s")'%directory)

    #sage0._eval_line('_start_time_=walltime()')
    HTML_Interface.protocol_version = "HTTP/1.0"

    tries = 0
    while True:
        try:
            server_address = (address, int(port))
            httpd = BaseHTTPServer.HTTPServer(server_address,
                                              HTML_Interface)
            sa = httpd.socket.getsockname()
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

    print "********************************************************"
    print "*                                                      *"
    print "* This address should open in your web browser:        *"
    print "*                                                      *"
    print "*          http://%s:%s"%(address, port)
    print "*                                                      *"
    print "********************************************************"
    print "Running log at %s"%log_file

    try:

        if viewer:
            os.system('%s http://%s:%s 1>&2 >/dev/null &'%(BROWSER, address, port))

        print "Press Control-C to interrupt a running calculation."
        print "If no calculation is running, press Control-C to return to SAGE."
        httpd.serve_forever()

    except KeyboardInterrupt, msg:

        print msg
        print "Shutting down server."

    current_workbook.save(workbook_sobj_file)

    return current_workbook




