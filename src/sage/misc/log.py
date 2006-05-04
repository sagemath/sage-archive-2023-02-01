r"""
Logging of SAGE sessions.

TODO: Pressing "control-D" can mess up the I/O sequence because of
a known bug.

You can create a log of your SAGE session as a web page and/or as a
latex document.  Just type \code{log_html()} to create an HTML log, or
\code{log_dvi()} to create a dvi (latex) log.  Your complete session
so far up until when you type the above command will be logged, along
with any future input.  Thus you can view the log system as a way to
print or view your entire session so far, along with a way to see
nicely typeset incremental updates as you work.

If \code{L=log_dvi()} or \code{L=log_html()} is a logger, you can type
\code{L.stop()} and \code{L.start()} to stop and start logging.

The environment variables \code{BROWSER} and \code{DVI\_VIEWER}
determine which web browser or dvi viewer is used to display your
running log.

For both log systems you must have a tex system installed on your
computer.  For HTML logging, you must have the convert command, which
comes with the free ImageMagick tools.

\note{The HTML output is done via Latex and png images right now,
sort of like how latex2html works.  Obviously it would be interesting
to do something using MathML in the long run.}

AUTHOR:
    -- William Stein (2006-02): initial version
    -- William Stein (2006-02-27): changed html generation so log directory
                                   is relocatable (no hardcoded paths).
    -- William Stein (2006-03-04): changed environment variable to BROWSER.
"""

# Note: there is a webbrowser module standard with Python.
# But it seems so dated as to be useless.

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import time

import interpreter
import latex
import misc


offset = 0
loggers = []

def update():
    for X in loggers:
        X._update()

try:
    PDF_VIEWER = os.environ['PDF_VIEWER']
except KeyError:
    PDF_VIEWER = 'acroread'

# Set the browsers and viewers.

if os.uname()[0] == "Darwin":

    # Simple on OS X, since there is an open command that opens anything,
    # using the user's preferences.
    DVI_VIEWER = 'open'
    BROWSER = 'open'
    PDF_VIEWER = 'open'

else:

    # Try to get something from the environment on other OS's.

    try:
        DVI_VIEWER = os.environ['DVI_VIEWER']
    except KeyError:
        DVI_VIEWER = 'xdvi -s 5'
    try:
	BROWSER = os.environ['BROWSER']
    except KeyError:
	if os.system('which firefox 1>/dev/null 2>/dev/null') == 0:
	    BROWSER = 'firefox'
	else:
	    BROWSER = 'konqueror'



REFRESH = ''

#REFRESH = '<meta http-equiv="REFRESH" content="4;">'

# The following javascript is from
#   http://grizzlyweb.com/webmaster/javascripts/refresh.asp#version2
# It isn't any good because it jumps the position of the browser
# view to the top of the page.

## REFRESH = """
## <noscript>
## <!--
##     We have the "refresh" meta-tag in case the user's browser does
##     not correctly support JavaScript or has JavaScript disabled.
##     Notice that this is nested within a "noscript" block.
## -->
## <meta http-equiv="refresh" content="2">

## </noscript>

## <script language="JavaScript">
## <!--

## var sURL = unescape(window.location.pathname);

## function doLoad()
## {
##     // the timeout value should be the same as in the "refresh" meta-tag
##     setTimeout( "refresh()", 2*1000 );
## }

## function refresh()
## {
##     //  This version of the refresh function will cause a new
##     //  entry in the visitor's history.  It is provided for
##     //  those browsers that only support JavaScript 1.0.
##     //
##     window.location.href = sURL;
## }
## //-->
## </script>

## <script language="JavaScript1.1">
## <!--
## function refresh()
## {
##     //  This version does NOT cause an entry in the browser's
##     //  page view history.  Most browsers will always retrieve
##     //  the document from the web-server whether it is already
##     //  in the browsers page-cache or not.
##     //
##     window.location.replace( sURL );
## }
## //-->
## </script>

## <script language="JavaScript1.2">
## <!--
## function refresh()
## {
##     //  This version of the refresh function will be invoked
##     //  for browsers that support JavaScript version 1.2
##     //
##     //  The argument to the location.reload function determines
##     //  if the browser should retrieve the document from the
##     //  web-server.  In our example all we need to do is cause
##     //  the JavaScript block in the document body to be
##     //  re-evaluated.  If we needed to pull the document from
##     //  the web-server again (such as where the document contents
##     //  change dynamically) we would pass the argument as 'true'.

##     window.location.reload( true );
## }
## //-->
## </script>
## </head>

## <!--
##     Use the "onload" event to start the refresh process.
## -->
## <body onload="doLoad()">
## """

class Log:
    """
    This is the base logger class.  The two classes that you
    actually instantiate are derived from this one.
    """
    def __init__(self, dir=None, debug=False, viewer=None):
        if dir is None:
            dir = misc.DOT_SAGE + 'log/log'
        self._time = time.strftime('%Y-%m-%d-%H%M%S')
        dir = os.path.abspath(dir) + '-' + self._time
        if not os.path.exists(dir):
            os.makedirs(dir)
        self._debug = debug
        self._n = 0
        self._dir = dir
        self._filename = '%s/%s'%(dir, self._filename())
        self._output = __IPYTHON__.output_hist
        self._input  = __IPYTHON__.input_hist_raw
        self._text = ''
        self._after_output = False
        self._init()
        self._input_text  = ''
        self._stopped = False
        self._viewer = viewer
        loggers.append(self)
        self._update()
        self.view()

    def __repr__(self):
        return "Logger"

    def _latex_(self):
        return "\\mbox{\\rm %s}"%self

    def dir(self):
        """
        Return the directory that contains the log files.
        """
        return self._dir

    def stop(self):
        """
        Stop the logger.  To restart use the start function.
        """
        self._stopped = True

    def start(self):
        """
        Start the logger.  To stop use the stop function.
        """
        self._stopped = False

    def _update(self):
        if self._stopped:
            return
        (O, I) = (self._output, self._input)
        #print "O:", O
        #print "I:", I
        K = O.keys()
        while self._n < max(len(I), max(K + [-1])):
            n = self._n
            followed_by_output = (n in K)
            if n < len(I):
                L = I[n]
            else:
                L = "Sorry; raw input log out of sync"
            Lstrip = L.strip()
            if len(Lstrip) > 0 and Lstrip[0] != '%':
                self._write(self._get_input(n, followed_by_output))
                self._input_text += L
            #m = n - offset
            m = n
            if m in K:
                self._write(self._get_output(m))
                s = '# ' + '\n# '.join(str(O[m]).split('\n')) + '\n\n'
            self._n += 1
        A = open(self._filename,'w')
        A.write(self._header() + '\n' + self._text + '\n' + self._footer())
        A.close()
        self._update_plain()
        self._build()

    def _write(self, lines):
        self._text += lines

    def _plain_text(self):
        return self._input_text

    def _input_log_name(self):
        return '%s/input-%s.sage'%(self._dir, self._time)

    def _update_plain(self):
        open(self._input_log_name(),'w').write(self._input_text)


class log_html(Log):
    r"""
    Create a running log of your SAGE session as a web page.

    Easy usage: \code{log_html()}

    TODO: Pressing "control-D" can mess up the I/O sequence because of
    a known bug.

    Use \code{L=log_html([optional directory])} to create an HTML log.
    Your complete session so far up until when you type the above
    command will be logged, along with any future input.  Thus you can
    view the log system as a way to print or view your entire session
    so far, along with a way to see nicely typeset incremental updates
    as you work.

    If L is a logger, you can type \code{L.stop()} and \code{L.start()} to
    stop and start logging.

    The environment variable \code{WEB\_BROWSER} determines which web
    browser or dvi viewer is used to display your running log.

    You must have a tex system installed on your computer, and you
    must have the convert command, which comes with the free
    ImageMagick tools.
    """

    def _init(self):
        self._images = '%s/images/'%self._dir
        if not os.path.exists(self._images):
            os.makedirs(self._images)


    def view(self):
        if not self._viewer is None:
            viewer = self._viewer
        else:
            viewer = BROWSER
        os.system('%s %s&'%(viewer, self._filename))

    def _build(self):
        return

    def __repr__(self):
        return "HTML Logger"

    def _get_input(self, n, followed_by_output):
        if n >= len(self._input):
            return
        return """<font color=darkblue>     %s %s:</font> %s"""%(
            n, interpreter._prompt, self._input[n])

    def _get_output(self, n):
        x = self._output[n]
        try:
            L = latex.latex(x)
        except:
            L = "\\mbox{error texing object}"
        single_png = '%s/%s.png'%(self._images, n)
        try:
            x.png(single_png)
        except AttributeError:
            latex.png(x, single_png, debug=self._debug)
        oi = '%s/images/o%s.html'%(self._dir, n)
        open(oi,'w').write('<pre>OUTPUT:\n%s\n\n\nLATEX:\n%s</pre>'%(x, L))
        return """<center> <table border=0 cellpadding=20 cellspacing=2
                bgcolor=lightgrey>
               <tr><td bgcolor=white>
               <a href="%s">
               <img src="%s" alt="%s">
                </a>
             </td></tr></table> </center>\n<hr>\n"""%('images/o%s.html'%n,
                                                      'images/%s.png'%n, L)

    def _filename(self):
        return 'index.html'

    def _header(self):
        T = self._title()
        inlog = os.path.split(self._input_log_name())[1]
        return '<html>\n%s\n<title>%s</title>\n<body><h1 align=center>%s</h1>\n<h2 align=center><a href="%s">%s</a></h2><pre>'%(REFRESH,
            T,T, inlog, inlog)

    def _footer(self):
        return "</pre>\n</body>\n</html>"

    def _title(self):
        return 'SAGE Log %s'%self._time


class log_dvi(Log):
    """
    Create a running log of your SAGE session as a nicely typeset dvi
    file.

    Easy usage: \code{log_dvi()}

    TODO: Pressing "control-D" can mess up the I/O sequence because of
    a known bug.

    Use \code{L=log\_dvi([optional directory])} to create a dvi log.
    Your complete session so far up until when you type the above
    command will be logged, along with any future input.  Thus you can
    view the log system as a way to print or view your entire session
    so far, along with a way to see nicely typeset incremental updates
    as you work.

    If L is a logger, you can type \code{L.stop()} and
    \code{L.start()} to stop and start logging.

    The environment variable \code{DVI\_VIEWER} determines which web
    browser or dvi viewer is used to display your running log.

    You must have a latex system installed on your computer and a dvi
    viewer.
    """
    def _init(self):
        SAGE_ROOT = os.environ['SAGE_ROOT']
        os.system('ln -sf %s/devel/doc/commontex/macros.tex %s/macros.tex'%(SAGE_ROOT, self._dir))
        self._in_verbatim = False

    def __repr__(self):
        return "DVI Logger"

    def _build(self):
        cmd = 'cd %s; latex \\\\nonstopmode \\\\input %s '%(
            self._dir, self._filename)
        cmd += " 2>/dev/null 1>/dev/null"
        dvifile = '%s.dvi'%os.path.splitext(self._filename)[0]
        if os.path.exists(dvifile):
            cmd += ' & '
        os.system(cmd)

    def view(self):
        if not self._viewer is None:
            viewer = self._viewer
        else:
            viewer = DVI_VIEWER
        self._build()
        F = os.path.splitext(self._filename)[0] + '.dvi'
        cmd = 'cd %s; %s %s '%(
            self._dir, viewer, F)
        os.system(cmd + " 2>/dev/null 1>/dev/null &")

    def _get_input(self, n, followed_by_output):
        if n >= len(self._input):
            return
        s = ''
        if not self._in_verbatim:
            s += '\\begin{verbatim}'
            self._in_verbatim = True
        I = self._input[n]
        s += "%s %s: %s"%(n,  interpreter._prompt, I)
        if followed_by_output:
            s += '\\end{verbatim}'
            self._in_verbatim = False
        self._after_output = False
        return s

    def _get_output(self, n):
        s = ''
        self._after_output = True
        if self._in_verbatim:
            s += '\\end{verbatim}\n'
            self._in_verbatim = False
        L = latex.latex(self._output[n])
        s += '\n\\begin{center}$\\displaystyle %s $\\end{center}\n'%L
        return s

    def _filename(self):
        return 'sagelog.tex'

    def _header(self):
        return """
\\documentclass{article}
\\input{macros}
\\title{%s}\\author{}
\\begin{document}
\\maketitle
"""%self._title()

    def _footer(self):
        if self._in_verbatim:
            return r"\end{verbatim}\end{document}"
        else:
            return r"\end{document}"


    def _title(self):
        return '\\SAGE Log %s'%self._time


#
# Coarse MathML logger
class log_html_mathml(Log):
    r"""
    Create a running log of your SAGE session as a mathml web page.

    Easy usage: \code{log_html_mathml()}

    TODO: Pressing "control-D" can mess up the I/O sequence because of
    a known bug.
    """
    def _init(self):
        SAGE_ROOT = os.environ['SAGE_ROOT']
        os.system('ln -sf %s/devel/doc/commontex/macros.tex %s/macros.tex'%(SAGE_ROOT, self._dir))
        os.system('ln -sf $SAGE_ROOT/data/hermes/dlt.tex %s/dlt.tex'% self._dir)
        os.system('ln -sf $SAGE_ROOT/data/hermes/*.xsl %s/'% self._dir)
        self._in_verbatim = False

    def __repr__(self):
        return "mathml Logger"

    def _build(self):

        seed = self._filename[:-4]+'.s.tex'
        dvi = self._filename[:-4]+'.s.dvi'
        lib = self._filename[:-4]+'.lib.xml'
        pub = self._filename[:-4]+'.pub.xml'

        cmd = 'cd %s; seed %s >/dev/null; latex \\\\nonstopmode \\\\input %s '%(
            self._dir,self._filename,seed)
        cmd += " 2>/dev/null 1>/dev/null"

        # turn to xml
        cmd += '; hermes %s > %s' %(dvi, lib)
        pubfile = '$SAGE_ROOT/data/hermes/pub.xslt'
        cmd += ';xsltproc %s %s > %s '%(pubfile, lib, pub)
        os.system(cmd)

        # Bug in XSLT generation forces us to do this
        f = open(pub,'r')
        pubstr = f.read()
        import re
        repl= re.sub("<math xmlns=\"http://www.w3.org/1998/Math/MathML\">",
               "<math xmlns=\"http://www.w3.org/1998/Math/MathML\" display=\"block\">",
               pubstr)

        open(pub,'w').write(repl)

    def view(self):
        if not self._viewer is None:
            viewer = self._viewer
        else:
            viewer = BROWSER
        self._build()

        os.system('%s  %s &' %(BROWSER, self._filename[:-4]+'.pub.xml'))

    def _get_input(self, n, followed_by_output):
        if n >= len(self._input):
            return
        s = ''
        if not self._in_verbatim:
            s += '\\begin{verbatim}'
            self._in_verbatim = True
        I = self._input[n]
        s += "%s %s: %s"%(n,  interpreter._prompt, I)
        if followed_by_output:
            s += '\\end{verbatim}'
            self._in_verbatim = False
        self._after_output = False
        return s

    def _get_output(self, n):
        s = ''
        self._after_output = True
        if self._in_verbatim:
            s += '\\end{verbatim}\n'
            self._in_verbatim = False
        L = latex.latex(self._output[n])
        s += '\n\\begin{center}$\\displaystyle %s $\\end{center}\n'%L
        return s

    def _filename(self):
        return 'sagelog.pub.xml'

    def _header(self):
        return """
\\documentclass{article}
\\input{macros}
\\title{%s}\\author{}
\\begin{document}
\\maketitle
"""%self._title()

    def _footer(self):
        if self._in_verbatim:
            return r"\end{verbatim}\end{document}"
        else:
            return r"\end{document}"

    def _title(self):
        return '\\SAGE Log %s'%self._time



