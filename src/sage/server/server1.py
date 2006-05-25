"""
SAGE Server 1 (of many)

AUTHOR:
    -- William Stein (2006-05-06): initial version

TODO:
   [] input one form shouldn't delete data from any other forms;
       e.g., you could be editing one form and submit another!
   [] The whole interface needs to be slimmed down so bunches of single
      line input (and output) will work.
   [] Ability to switch from one log (=workbook) to another via
      the web interface.
   [x] The "move to the current input box" javascript *only* works
      with firefox (not opera, not konqueror); also this should
      just keep the page position where it is rather than move it.
      Moving to a more AJAX-ish model would alternatively fix this, maybe.
   [x] A. Clemesha: shrink/expand input/output blocks
   [] Add plain text annotation that is not evaluated
      between blocks (maybe in html?)
      E.g., just make ctrl-enter on a block by HTML-it.
   [] Ability to interrupt running calculations directly
      from the web interface (no console access)
   [] Nice animation while a computation is proceeding.
   [] Some way to show output as it is computed.
   [] Option to delete blocks
   [] Make block expand if enter a lot of text into it.
   [] Evaluate the entire worksheet
   [] Theme-able / skin-able
   [] Downloading and access to exact log of IO to client SAGE process
   [] Save session objects as to log objects so don't have to re-eval?
   [] The entire page is resent/updated every time you hit shift-enter;
      using 'AJAX' this flicker/lag could be completely eliminated.
   [] When pressing shift-enter a line feed is inserted temporarily
      into the inbox, which is unnerving.
   [] Add authentication
   [] Embed the log object in the html file, so server session
      can be restared directly using the html file!  E.g., embed
      pickled Log object in a comment at end of the .html file.
   [] Ability to upload and download source files (to be run)
      via web interface; maybe ability to edit them too, via some
      'rich' code editing javascript 'widget'.
   [] rewrite tables using CSS
   [] load and attaching scripts.
   [] a way to interactively watch the output of a running computation
      (in verbose mode).
   [] undo -- have infinite undo and redo of the SAGE *log*, i.e.,
      text I/O (and possibly graphics).  does not save *state* of
      the "sage kernel".

DONE
   [x] A. Clemesha: When hit shift-enter the next text box should be made
      into focus.
   [x] Embedded graphics from plots;
       also embed png from latexing of math objects (so they look nice).
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


# Standard Python libraries
import BaseHTTPServer
import cgi
import os, sys
import select
import shutil
import socket
from   StringIO import StringIO


# SAGE libraries
import sage.interfaces.sage0

import sage.misc.banner
import sage.misc.misc
import sage.misc.preparser
from   sage.misc.viewer     import BROWSER
from   sage.ext.sage_object import load, SageObject


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

        #if not render_if_empty and len(self.cmd.strip()) == 0 \
        #   and (self.out.strip()) == 0:
        #    return ''

        w = max((ncols - 15), 3)
	#html_in is the html/javascript that starts a 'cell'
	#this probably needs some refactoring by a css/javascript pro :)
	html_in ="""
           <table  border=0 cellspacing=2 cellpadding=2 bgcolor='#FFFFFF' width='80%%'>
	    <td align=center>
	"""
	#below in the javascript that enables the +/- resizing of the textarea
	controls ="""
        <div class="controlarea">
        <table cellpadding=0 cellspacing=0>
         <tr>
          <td><span class="control"><input type='submit' class="btn" value="&gt;"></span></td>
          <td><span class="control"><a class="cs" href="javascript:changeAreaSize(-1,'in%s')"><b>-</b></a></span></td>
         </tr>
         <tr>
          <td><span class="control"><a class="cs" href="javascript:toggleVisibility(%s);"><b id='tog%s'>H</b></a></span></td>
          <td><span class="control"><a class="cs" href="javascript:changeAreaSize(1,'in%s')"><b>+</b></a></span></td>
         </tr>
        </table>
        </div>
	"""%(number, number,number,number)
	#below is where the textarea form is added
	html_in +="""
        </td>
        <tr><td align=center bgcolor='#DDDDDD' cellpadding=1 cellspacing=0>
        <table border=0 cellpadding=4 cellspacing=0 bgcolor='#FFFFFF' width='100%%'><tr><td>
                 <textarea style="border: 0px;width: 100%%;"
                   name='%s' bgcolor='%s' rows='%s'
                   cols='%s' id='in%s' onkeypress='ifShiftEnter(%s,event);'>%s</textarea>
        </td>
        </tr></table></td>
        <td align='right' valign='top'>%s</td>
        </tr></table>
        """%(number, cmd_color, cmd_nrows, ncols, number, number, self.cmd,controls)

        button = '<input align=center name="exec" type="submit" id="with4" \
                    value="%sEnter %s(%s)">'%(' '*w, ' '*w, number)

        if len(self.file_list) == 0:
            files  = ''
            images = ''
        else:
            images = []
            files  = []
            for F in self.file_list:
                if(F[-4:] == '.png'):
                    images.append('<img src="%s/%s">'%(number,F))
                else:
                    files.append('<a href="%s/%s">%s</a>'%(number,F,F))

            if len(images) == 0:
                images = ''
            else:
                images = "<br>%s"%'<br>'.join(images)
            if len(files)  == 0:
                files  = ''
            else:
                files  = ('&nbsp'*3).join(files)


        out = self.out

        if len(out) >  0:
            #out_nrows = min(out_nrows, len(out.split('\n')))
            out_nrows = len(out.split('\n'))
            #<textarea style="color:blue" readonly rows="%s" cols="%s">%s</textarea>
            html_out = """
             <table   border=0 bgcolor='%s' cellpadding=0><tr>
             <td bgcolor='#FFFFFF' align=left>
             <font color='blue'><pre>%s</pre>%s</font>
             </td></tr></table>
             """%(out_color, out.replace('<','&lt;'), files)
        else:
            html_out = ''

        c = """
        <form name="io%s" method=post action="" id="%s">
        %s
        <div id='out%s'>
        %s
        %s
        %s
        </div>
        </form>"""%(number, number, html_in, number, html_out, files, images)
        return c

class Log(SageObject):
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
            L = self._log[i]
            s += L.html(i, numrows, ncols, cmd_color=color)
        return s

    def set_last_cmd(self, cmd):
        if len(self._log) == 0:
            self._log = [IO_Line()]
        self._log[-1].cmd = str(cmd)

    def set_last_out(self, out):
        if len(self._log) == 0:
            self._log = [IO_Line()]
        self._log[-1].out = str(out)



class HTML_Interface(BaseHTTPServer.BaseHTTPRequestHandler):
    def __files(self,number):
        global directory
        return os.listdir("%s/%s"%(directory,number))

    def __show_page(self, number):
        global current_log
        f = self.send_head()
        if f:
            f = StringIO()
            f.write("""
            <html><head><title>SAGE Calculator %s</title></head>\n

	    <style>
	    div.controlarea{
                text-align:right;
            }
            span.control {
            border:1px solid white;
            font-family:Arial,Verdana, sans-serif;
            font-size:10pt;
            }

            input.btn {
              color:#999999;
              text-decoration:none;
              background: white;
              padding:0px;
              margin:0px;
              border:1px solid white;
            }
            input.btn:hover {
              color:black;
              text-decoration:none;
              background: white;
              padding:0px;
              margin:0px;
              border:1px solid #333333;
            }

            span.control a.cs {
            color:#999999;
            text-decoration:none;
            border:1px solid white;
            }
            span.control:hover a.cs, span.control a:hover.cs {
            color:black;
            border:1px solid #333333;
            }

            </style>

            <script language=javascript>
            function scroll_to_bottom() {
                document.getElementById(%s).scrollIntoView();
                document.io%s.elements[0].focus();
            }

            function ifShiftEnter(number, event) {
                var theCode = event.keyCode ? event.keyCode :
                               event.which ? event.which : event.charCode;
                if (theCode == 13 && event.shiftKey) {
                   document.forms[number].submit();
                   return false;
                }
                else
                   return true;
            }

	    function changeAreaSize(val,id) {
                var el = document.getElementById(id);
                if (val==1)
                    el.rows = el.rows + 1;
                else
                    el.rows = el.rows - 1;
            }

            function toggleVisibility(id) {
                var outBox = document.getElementById('out'+id);
                if(outBox.style.display == 'none') {
                    outBox.style.display = 'block';
                    document.getElementById('tog'+id).innerHTML='H';
                } else {
                    outBox.style.display = 'none';
                    document.getElementById('tog'+id).innerHTML='S';
                }
            }
            function showBox(id) {
                var el = document.getElementById(id);
                el.style.display='none';
            }
            </script>
            """%(save_name, number, number+1))

            f.write('<body onload="javascript:scroll_to_bottom()"><div align=center><H2><font color="darkgreen"><a href="http://modular.math.washington.edu/sage">SAGE</a>: Software for Algebra and Geometry Experimentation</font></H2>')
            f.write('<h3<tr>%s</h3>'%sage.misc.banner.version())
            f.write('(enter input and press shift-enter)')
            if len(current_log) == 0 or current_log[-1].cmd != '':
                I = IO_Line()
                current_log.append(I)
            #f.write('<br><hr><h2 align=center>Complete Session Log</h2>')
            #f.write("""<table width=90%% align=center bgcolor='#CCCCCC' cellpadding=10>
            #<tr><td bgcolor='#FFFFFF'>
            #  <pre>%s</pre>
            #  </td></tr></table>
            #  """%fulltext_log)
            f.write('<hr>')
            if save_name:
                current_log.save(save_name)
            f.write(current_log.html(numcols))
            f.write('</div><br><br><br><br></body></html>')
            f.seek(0)
            self.copyfile(f, self.wfile)
            f.close()
            return f

    def do_GET(self):
        if self.path[-4:] == '.png' or self.path[-5:] == '.sobj':
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
        else:
            self.__show_page(0)

    def do_POST(self):
        global current_log, fulltext_log
        ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))
        length = int(self.headers.getheader('content-length'))
        if ctype == 'multipart/form-data':
            self.body = cgi.parse_multipart(self.rfile, pdict)
        elif ctype == 'application/x-www-form-urlencoded':
            qs = self.rfile.read(length)
            C = cgi.parse_qs(qs, keep_blank_values=1)
            number = eval(C.keys()[0])
            current_dir = "%s/%d"%(directory,number)
            code_to_eval = C[C.keys()[0]][0]
            fulltext_log += '\n#%s\n'%('-'*70) + '\n' + code_to_eval + '\n\n'
            try:
                if number > len(current_log)-1:
                    current_log.set_last_cmd(code_to_eval)
                    number = len(current_log)-1
                else:
                    # re-evaluating a code block
                    current_log[number].cmd = code_to_eval
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
                    os.mkdir(current_dir)

                for F in os.listdir(current_dir):
                    os.unlink("%s/%s"%(current_dir,F))

                sage0._eval_line('os.chdir("%s")'%current_dir)

                try:
                    o = sage0._eval_line('execfile("%s/_temp_.py")'%directory)
                    #sage0._send('execfile("%s/_temp_.py")'%directory)
                except KeyboardInterrupt, msg:
                    print "Keyboard interrupt!"
                    o = msg

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

                o = sage.misc.misc.word_wrap(o, ncols=numcols)

                fulltext_log += '\n'.join(o.split('\n')) + '\n'

                current_log[number].out = o
                current_log[number].file_list = self.__files(number)
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

sage0=None
def server_http1(name='default_server',
                 port=8000, address='localhost', ncols=90,
                 nrows=8, dir=None, viewer=True, log=None,
                 max_tries=10):
    """
    Start a SAGE http server at the given port.

    Typical usage:
        server_http1('mysession')

    Use it.  To start it later, just type server_http1('mysession') again.

    INPUT:
        name -- name of the server; all I/O is saved in the file
                with that name in current directory.  If you restart
                the server with that same name then it will restart
                in the state you left it (though of course none of the
                blocks will have been evaluated).
        port -- port on computer where the server is served
        max_tries -- maximum number of ports > port to try in case
                port can't be opened.
    """
    global directory, fulltext_log, current_log, \
           files, numcols, numrows, sage0, save_name
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
    if isinstance(name, str):
        if name[-5:] != '.sobj':
            name += '.sobj'
    save_name = name[:-5]
    fulltext_log = ''
    if log is None and os.path.exists(name):
        try:
            log = load(name)
        except IOError:
            print "Unable to load log %s (creating new log)"%name

    if log is None:
        current_log = Log()
    else:
        current_log = log
    sage0 = sage.interfaces.sage0.Sage(logfile=logfile)
    sage0.eval('os.chdir("%s")'%directory)
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
    print "*     Open the following address in your web browser:  *"
    print "*                                                      *"
    print "*       http://%s:%s"%(address, port)
    print "*                                                      *"
    print "********************************************************"
    print "Running log at %s"%logfile

    try:

        if viewer:
            #os.system('%s file:///%s&'%(BROWSER, logfile))
            os.system('%s http://%s:%s 1>&2 >/dev/null &'%(BROWSER, address, port))
        print "Press Control-C to interrupt a running calculation."
        print "If no calculation is running, press Control-C to return to SAGE."
        httpd.serve_forever()

    except KeyboardInterrupt, msg:

        print msg
        print "Shutting down server."

    current_log.save(name)

    return current_log




