"""
SAGE Server 1 (of many)

AUTHOR:
    -- William Stein (2006-05-06): initial version
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


import BaseHTTPServer
from StringIO import StringIO
import os, sys
import shutil
import cgi
import select
import sage.interfaces.sage0
import sage.misc.preparser
import sage.misc.misc
import sage.misc.banner
from sage.misc.log import BROWSER

class IO_Line:
    def __init__(self, cmd='', out=''):
        self.cmd = cmd
        self.out = out
        self.file_list = []

    def __repr__(self):
        return 'INPUT: %s\nOUTPUT: %s\nFILES: %s'%(self.cmd, self.out, self.file_list)

    def html(self, number, nrows, ncols, out_nrows=4,
             cmd_color='#FFFFFF', out_color='#FFFFFF'):
        if self.cmd != '':
            cmd_nrows = len(self.cmd.split('\n'))
        else:
            cmd_nrows = nrows


        w = max((ncols - 15), 3)
        html_in = """
        <table width=95%% align=center border=0 cellpadding=2 bgcolor='#FFFFFF'>
        <tr><td align=center>
         <textarea name='input'  bgcolor='%s' rows='%s' cols='%s' id='input'>%s</textarea>
         </td></tr></table>
         """%(cmd_color, cmd_nrows, ncols, self.cmd)

        button = '<input align=center name="exec" type="submit" id="with4" \
                    value="%sEnter %s(%s)">'%(' '*w, ' '*w, number)

        if len(self.file_list) == 0:
            files = ''
        else:
            files = '<br>' + ('&nbsp'*3).join(['<a href="%s">%s</a>'%(F,F) for F in self.file_list])

        out = self.out

        if len(out) >  0:
            out_nrows = min(out_nrows, len(out.split('\n')))
            html_out = """
             <table width=95%%  align=center border=0 bgcolor='%s' cellpadding=2><tr>
             <td bgcolor='#FFFFFF' align=center>
             <textarea rows="%s" cols="%s">%s</textarea>
             %s</td></tr></table><br>
             """%(out_color, out_nrows, ncols, out, files)
        else:
            html_out = ''

        c = """<div align=center>
        <form name="" method=post action="">
        %s
        %s
        %s
        </form></div>"""%(button, html_in, html_out)
        return c

class Log:
    def __init__(self):
        self._log = []

    def __repr__(self):
        return str(self._log)

    def __len__(self):
        return len(self._log)

    def __getitem__(self, n):
        return self._log[n]

    def append(self, L):
        self._log.append(L)

    def html(self, ncols):
        n = len(self._log)
        s = ''
        for i in range(len(self._log)):
            if i == 0:
                color = "#FFFFFF"
            else:
                color = "#FFFFFF"
            j = n-i-1
            L = self._log[j]
            s += L.html(j, numrows, ncols, cmd_color=color) + '<hr>'
        return s

    def set_last_cmd(self, cmd):
        if len(self._log) == 0:
            self._log = [IO_Line()]
        self._log[-1].cmd = str(cmd)

    def set_last_out(self, out):
        if len(self._log) == 0:
            self._log = [IO_Line()]
        self._log[-1].out = str(out)

    def set_last_filelist(self, files):
        if len(self._log) == 0:
            self._log = [IO_Line()]
        self._log[-1].file_list = files



class HTML_Interface(BaseHTTPServer.BaseHTTPRequestHandler):
    def __new_files(self):
        global directory, files
        new = []
        for k in os.listdir(directory):
            if not k in files and k != '_temp_.py':
                new.append(k)
                files.append(k)
        new.sort()
        return new

    def __show_page(self):
        global log
        f = self.send_head()
        if f:
            f = StringIO()
            f.write("<html><head><title>SAGE Calculator</title></head>\n")
            f.write('<body><div align=center><H1><font color="darkgreen">SAGE</font></H1>')
            f.write('<h2<tr>%s</h2></div>'%\
                    sage.misc.banner.version())
            if len(log) == 0 or log[-1].cmd != '':
                I = IO_Line()
                log.append(I)
            f.write(log.html(numcols))
            f.write('<br><hr><h2 align=center>Complete Session Log</h2>')
            f.write("""<table width=90%% align=center bgcolor='#CCCCCC' cellpadding=10>
            <tr><td bgcolor='#FFFFFF'>
              <pre>%s</pre>
              </td></tr></table>
              """%fulltext_log)
            f.seek(0)
            self.copyfile(f, self.wfile)
            f.close()
            return f

    def do_GET(self):
        self.__show_page()

    def do_POST(self):
        global log, fulltext_log
        ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))
        length = int(self.headers.getheader('content-length'))
        if ctype == 'multipart/form-data':
            self.body = cgi.parse_multipart(self.rfile, pdict)
        elif ctype == 'application/x-www-form-urlencoded':
            qs = self.rfile.read(length)
            C = cgi.parse_qs(qs, keep_blank_values=1)
            code_to_eval = C['input'][0]
            fulltext_log += '\n#%s\n'%('-'*70) + '\n' + code_to_eval + '\n\n'
            number = eval(C['exec'][0].split()[-1])
            #print "INPUT:\n%s"%code_to_eval
            try:
                if number > len(log)-1:
                    log.set_last_cmd(code_to_eval)
                    number = len(log)-1
                else:
                    log[number].cmd = code_to_eval
                s = sage.misc.preparser.preparse_file(code_to_eval, magic=False,
                                                      do_time=True, ignore_prompts=True)
                s = [x for x in s.split('\n') if len(x.split()) > 0]   # remove all blank lines
                if len(s) > 0:
                    t = s[-1]
                    if len(t) > 0 and t[0] != ' ' and t[0] != '\t' and t[:5] != 'print' \
                       and t[:4] != 'time' and not ('=' in t):
                        s[-1] = 'print %s'%s[-1]
                s = '\n'.join(s)

                open('%s/_temp_.py'%directory, 'w').write(s)
                try:
                    o = sage0._eval_line('execfile("%s/_temp_.py")'%directory)
                except KeyboardInterrupt, msg:
                    print "Keyboard interrupt!"
                    o = msg

                #print 'OUTPUT:\n%s'%o
                o = sage.misc.misc.word_wrap(o, ncols=numcols)

                fulltext_log += '# ' + '\n# '.join(o.split('\n')) + '\n'

                log[number].out = o
                log[number].file_list = self.__new_files()
                self.__show_page()

            except (RuntimeError, TypeError), msg:
                print "ERROR!!!", msg
                self.__show_page()

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


    def send_head(self):
        self.send_response(200)
        self.send_header("Content-type", 'text/html')
        self.end_headers()
        f = StringIO()
        #print "URL Path: %s\n" % self.path
        f.seek(0)
        return f

    def copyfile(self, source, outputfile):
        shutil.copyfileobj(source, outputfile)

sage0=None
def server_http1(port=8000, address='localhost', ncols=90,
                nrows=8, dir=None, viewer=False):
    global directory, fulltext_log, log, files, numcols, numrows, sage0
    remove_dir = False
    if dir is None:
        remove_dir = True
        directory = sage.misc.misc.tmp_dir('server')
    else:
        directory = os.path.abspath(dir)
    logfile = '%s/logfile.txt'%directory
    open(logfile,'w').close()  # touch the file
    #os.system('tail -f %s&'%logfile)
    numcols = int(ncols)
    numrows = int(nrows)
    files = os.listdir(directory)
    fulltext_log = ''
    log = Log()
    sage0 = sage.interfaces.sage0.Sage(logfile=logfile)
    sage0.eval('os.chdir("%s")'%directory)
    server_address = (address, int(port))
    HTML_Interface.protocol_version = "HTTP/1.0"
    httpd = BaseHTTPServer.HTTPServer(server_address,
                                      HTML_Interface)
    sa = httpd.socket.getsockname()
    #print "Serving HTTP on", sa[0], "port", sa[1], "..."
    print "Web interface at http://%s:%s"%(address, port)
    print "Running log at %s"%logfile

    try:

        if viewer:
            os.system('%s file:///%s&'%(BROWSER, logfile))
            os.system('%s http://%s:%s&'%(BROWSER, address, port))
        print "Press Control-C to interrupt a running calculation."
        print "If no calculation is running, press Control-C to return to SAGE."

        httpd.serve_forever()

    except KeyboardInterrupt, msg:

        print msg
        print "Shutting down server."






