"""nodoctest
A Cell.

A cell is a single input/output block.  Worksheets are built out of a
list of cells.
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

# Maximum number of characters allowed in output.  This is
# needed avoid overloading web browser.  For example, it
# should be possible to gracefully survive:
#    while True:
#       print "hello world"
# On the other hand, we don't want to loose the output of big matrices
# and numbers, so don't make this too small.

MAX_OUTPUT = 32000
MAX_OUTPUT_LINES = 120

TRACEBACK = 'Traceback (most recent call last):'

import re

# This regexp matches "cell://blah..." in a non-greedy way (the ?), so
# we don't get several of these combined in one.
re_cell = re.compile('"cell://.*?"')
re_cell_2 = re.compile("'cell://.*?'")   # same, but with single quotes

import os, shutil

from   sage.misc.misc import word_wrap
from   sage.misc.html import math_parse
from   sage.misc.preparser import strip_string_literals

import notebook

import worksheet

class Cell_generic:
    def is_interactive_cell(self):
        return False

class TextCell(Cell_generic):
    def __init__(self, id, text, worksheet):
        self.__id = int(id)
        self.__text = text
        self.__worksheet = worksheet

    def set_input_text(self, input_text):
        self.__text = input_text

    def set_worksheet(self, worksheet, id=None):
        self.__worksheet = worksheet
        if not id is None:
            self.__id = id

    def html(self, ncols, do_print=False, do_math_parse=True):
        """
        INPUT:
            do_math_parse -- bool (default: True)
                If True, call math_parse (defined in cell.py)
                on the html.
        """
        t = self.__text
        if do_math_parse:
            # Do dollar sign math parsing
            t = math_parse(t)
        s = '<div><font size=+1>%s</font></div>'%t
        return s

    def plain_text(self, prompts=False):
        return self.__text

    def edit_text(self):
        return self.__text

    def id(self):
        return self.__id

    def is_auto_cell(self):
        return False

    def __cmp__(self, right):
        return cmp(self.id(), right.id())

    def set_cell_output_type(self, typ='wrap'):
        pass # ignored


class Cell(Cell_generic):
    def __init__(self, id, input, out, worksheet):
        self.__id    = int(id)
        self.__in    = str(input).replace('\r','')
        self.__out   = str(out).replace('\r','')
        self.__worksheet = worksheet
        self.__interrupted = False
        self.__completions = False
        self.has_new_output = False
        self.__version = 0
        self.__no_output_cell = False
        self.__asap = False

    def set_asap(self, asap):
        self.__asap = bool(asap)

    def is_asap(self):
        """
        Return True if this is an asap cell, i.e., evaluation of it is
        done as soon as possible.
        """
        try:
            return self.__asap
        except AttributeError:
            self.__asap = False
            return self.__asap

    def set_no_output(self, no_output):
        self.__no_output = bool(no_output)

    def is_no_output(self):
        """
        Return True if this is an no_output cell, i.e., a cell for
        which we don't care at all about the output.
        """
        try:
            return self.__no_output
        except AttributeError:
            self.__no_output = False
            return self.__no_output

    def set_cell_output_type(self, typ='wrap'):
        self.__type = typ

    def cell_output_type(self):
        try:
            return self.__type
        except AttributeError:
            self.__type = 'wrap'
            return 'wrap'

    def set_worksheet(self, worksheet, id=None):
        self.__worksheet = worksheet
        if not id is None:
            self.set_id(id)

    def update_html_output(self, output=''):
        """
        Update the list of files with html-style links or embeddings
        for this cell.

        For interactive cells the html output section is always empty,
        mainly because there is no good way to distinguish content
        (e.g., images in the current directory) that goes into the
        interactive template and content that would go here.
        """
        if self.is_interactive_cell():
            self.__out_html = ""
        else:
            self.__out_html = self.files_html(output)

    def id(self):
        return self.__id

    def set_id(self, id):
        self.__id = int(id)

    def worksheet(self):
        return self.__worksheet

    def worksheet_filename(self):
        return self.__worksheet.filename()

    def notebook(self):
        return self.__worksheet.notebook()

    def directory(self):
        dir = self._directory_name()
        if not os.path.exists(dir):
            os.makedirs(dir)
        return dir

    def _directory_name(self):
        return '%s/cells/%s'%(self.__worksheet.directory(), self.id())


    def __cmp__(self, right):
        return cmp(self.id(), right.id())

    def __del__(self):
        dir = self._directory_name()
        if os.path.exists(dir):
            shutil.rmtree(dir, ignore_errors=True)

    def __repr__(self):
        return 'Cell %s; in=%s, out=%s'%(self.__id, self.__in, self.__out)

    def word_wrap_cols(self):
        try:
            return self.notebook().conf()['word_wrap_cols']
        except AttributeError:
            return 70

    def plain_text(self, ncols=0, prompts=True, max_out=None, wiki_out=False):
        if ncols == 0:
            ncols = self.word_wrap_cols()
        s = ''

        input_lines = self.__in
        pr = 'sage: '

        if prompts:
            input_lines = input_lines.splitlines()
            has_prompt = False
            if pr == 'sage: ':
                for v in input_lines:
                    w = v.lstrip()
                    if w[:5] == 'sage:' or w[:3] == '>>>' or w[:3] == '...':
                        has_prompt = True
                        break
            else:
                # discard first line since it sets the prompt
                input_lines = input_lines[1:]

            if has_prompt:
                s += '\n'.join(input_lines) + '\n'
            else:
                in_loop = False
                for v in input_lines:
                    if len(v) == 0:
                        pass
                    elif len(v.lstrip()) != len(v):  # starts with white space
                        in_loop = True
                        s += '...   ' + v + '\n'
                    elif v[:5] == 'else:':
                        in_loop = True
                        s += '...   ' + v + '\n'
                    else:
                        if in_loop:
                            s += '...\n'
                            in_loop = False
                        s += pr + v + '\n'
        else:
            s += self.__in

        if prompts:
            msg = TRACEBACK
            if self.__out.strip().startswith(msg):
                v = self.__out.strip().splitlines()
                w = [msg, '...']
                for i in range(1,len(v)):
                    if not (len(v[i]) > 0 and v[i][0] == ' '):
                        w = w + v[i:]
                        break
                out = '\n'.join(w)
            else:
                out = self.output_text(ncols, raw=True, html=False)
        else:
            out = self.output_text(ncols, raw=True, html=False, allow_interact=False)
            out = '///\n' + out

        if not max_out is None and len(out) > max_out:
            out = out[:max_out] + '...'

        # Get rid of spurious carriage returns
        s = s.strip('\n')
        out = out.strip('\n').strip('\r').strip('\r\n')
        s = s + '\n' + out

        if not prompts:
            s = s.rstrip('\n')
        return s

    def edit_text(self, ncols=0, prompts=False, max_out=None):
        s = self.plain_text(ncols,prompts,max_out,wiki_out=True)
        return '{{{id=%s|\n%s\n}}}'%(self.id(), s)

    def is_last(self):
        return self.__worksheet.cell_list()[-1] == self

    def next_id(self):
        L = self.__worksheet.cell_list()
        try:
            k = L.index(self)
        except ValueError:
            print "Warning -- cell %s no longer exists"%self.id()
            return L[0].id()
        try:
            return L[k+1].id()
        except IndexError:
            return L[0].id()

    def interrupt(self):
        self.__interrupted = True

    def interrupted(self):
        return self.__interrupted

    def computing(self):
        return self in self.__worksheet.queue()

    def is_interactive_cell(self):
        """
        Return True if this cell contains the use of interact either
        as a function call or a decorator.
        """
        # Do *not* cache
        s, _ = strip_string_literals(self.input_text())
        return bool(re.search('(?<!\w)interact\s*\(.*\).*', s) or re.search('\s*@\s*interact\s*\n', s))

    def is_interacting(self):
        return hasattr(self, 'interact')

    def stop_interacting(self):
        if self.is_interacting():
            del self.interact

    def set_input_text(self, input):
        # Stuff to deal with interact
        if input.startswith('%__sage_interact__'):
            self.interact = input[len('%__sage_interact__')+1:]
            self.__version = 1+self.version()
            return
        elif self.is_interacting():
            try:
                del self.interact
                del self._interact_output
            except AttributeError:
                pass

        self.__version = 1+self.version()
        self.__in = input
        if hasattr(self, '_html_cache'):
            del self._html_cache

    def input_text(self):
        return self.__in

    def is_auto_cell(self):
        return '#auto' in self.__in

    def changed_input_text(self):
        try:
            t = self.__changed_input
            del self.__changed_input
            return t
        except AttributeError:
            return ''

    def set_changed_input_text(self, new_text):
        self.__changed_input = new_text
        self.__in = new_text

    def set_output_text(self, output, html, sage=None):
        if output.count('<?__SAGE__TEXT>') > 1:
            html = '<h3><font color="red">WARNING: multiple @interacts in one cell disabled (not yet implemented).</font></h3>'
            output = ''

        # In interacting mode, we just save the computed output
        # (do not overwrite).
        if self.is_interacting():
            self._interact_output = (output, html)
            return

        if hasattr(self, '_html_cache'):
            del self._html_cache

        output = output.replace('\r','')
        if len(output) > MAX_OUTPUT or output.count('\n') > MAX_OUTPUT_LINES:
            url = ""
            if not self.computing():
                file = "%s/full_output.txt"%self.directory()
                open(file,"w").write(output)
                url = "<a target='_new' href='%s/full_output.txt' class='file_link'>full_output.txt</a>"%(
                    self.url_to_self())
                html+="<br>" + url
            lines = output.splitlines()
            start = '\n'.join(lines[:MAX_OUTPUT_LINES/2])[:MAX_OUTPUT/2]
            end = '\n'.join(lines[-MAX_OUTPUT_LINES/2:])[-MAX_OUTPUT/2:]
            warning = 'WARNING: Output truncated!  '
            if url:
                # make the link to the full output appear at the top too.
                warning += '\n<html>%s</html>\n'%url
            output = warning + '\n\n' + start + '\n\n...\n\n' + end
        self.__out = output
        if not self.is_interactive_cell():
            self.__out_html = html
        self.__sage = sage

    def sage(self):
        try:
            return self.__sage
        except AttributeError:
            return None

    def set_introspect_html(self, html, completing=False):
        if completing:
            self.__introspect_html = html
        else:
            html = html.replace('<','&lt;').strip()
            self.__introspect_html = '<pre class="introspection">'+html+'</pre>'

    def output_html(self):
        try:
            return self.__out_html
        except AttributeError:
            self.__out_html = ''
            return ''

    def introspect_html(self):
        if not self.introspect():
            return ''
        try:
            return self.__introspect_html
        except AttributeError:
            self.__introspect_html = ''
            return ''

    def process_cell_urls(self, x):
        end = '?%d"'%self.version()
        begin = self.url_to_self()
        for s in re_cell.findall(x) + re_cell_2.findall(x):
            x = x.replace(s,begin + s[7:-1] + end)
        return x

    def output_text(self, ncols=0, html=True, raw=False, allow_interact=True):
        if allow_interact and hasattr(self, '_interact_output'):
            # Get the input template
            z = self.output_text(ncols, html, raw, allow_interact=False)
            if not '<?__SAGE__TEXT>' in z or not '<?__SAGE__HTML>' in z:
                return z
            if ncols:
                # Get the output template
                try:
                    # Fill in the output template
                    output,html = self._interact_output
                    output = self.parse_html(output, ncols)
                    z = z.replace('<?__SAGE__TEXT>', output)
                    z = z.replace('<?__SAGE__HTML>', html)
                    return z
                except (ValueError, AttributeError), msg:
                    print msg
                    pass
            else:
                # Get rid of the interact div to avoid updating the wrong output location
                # during interact.
                return ''

        is_interact = self.is_interactive_cell()
        if is_interact and ncols == 0:
            if 'Traceback (most recent call last)' in self.__out:
                s = self.__out.replace('cell-interact','')
                is_interact=False
            else:
                return '<h2>Click to the left again to hide and once more to show the dynamic interactive window</h2>'
        else:
            s = self.__out

        if raw:
            return s

        if html:
            s = self.parse_html(s, ncols)

        if not is_interact and not self.is_html() and len(s.strip()) > 0:
            s = '<pre class="shrunk">' + s.strip('\n') + '</pre>'

        return s.strip('\n')

    def parse_html(self, s, ncols):
        def format(x):
            return word_wrap(x.replace('<','&lt;'), ncols=ncols)

        def format_html(x):
            return self.process_cell_urls(x)

        # if there is an error in the output,
        # specially format it.
        if not self.is_interactive_cell():
            s = format_exception(format_html(s), ncols)

        # Everything not wrapped in <html> ... </html>
        # should have the <'s replaced by &lt;'s
        # and be word wrapped.
        t = ''
        while len(s) > 0:
            i = s.find('<html>')
            if i == -1:
                t += format(s)
                break
            j = s.find('</html>')
            if j == -1:
                t += format(s[:i])
                break
            t += format(s[:i]) + format_html(s[i+6:j])
            s = s[j+7:]
        t = t.replace('</html>','')
        return t


    def has_output(self):
        return len(self.__out.strip()) > 0

    def is_html(self):
        try:
            return self.__is_html
        except AttributeError:
            return False

    def set_is_html(self, v):
        self.__is_html = v

    def introspect(self):
        try:
            return self.__introspect
        except AttributeError:
            return False

    def unset_introspect(self):
        self.__introspect = False

    def set_introspect(self, before_prompt, after_prompt):
        self.__introspect = [before_prompt, after_prompt]

    def evaluate(self, introspect=False, time=False, username=None):
        """
        INPUT:
            username -- name of user doing the evaluation
            time -- if True return time computation takes
            introspect -- either False or a pair [before_cursor, after_cursor] of strings.
        """
        self.__interrupted = False
        self.__time = time
        self.__introspect = introspect
        self.__worksheet.enqueue(self, username=username)
        self.__type = 'wrap'
        dir = self.directory()
        for D in os.listdir(dir):
            F = dir + '/' + D
            try:
                os.unlink(F)
            except OSError:
                try:
                    shutil.rmtree(F)
                except:
                    pass

    def version(self):
        try:
            return self.__version
        except AttributeError:
            self.__version = 0
            return self.__version

    def time(self):
        try:
            return self.__time
        except AttributeError:
            return "?"

    def do_time(self):
        self.__time = True

    def doc_html(self, wrap=None, div_wrap=True, do_print=False):
        """Modified version of \code{self.html} for the doc browser. This is a hack and needs to be improved.
        The problem is how to get the documentation html to display nicely between the example cells.
        The type setting (jsMath formating) needs attention too.
        """
        self.evaluate()
        if wrap is None:
            wrap = self.notebook().conf()['word_wrap_cols']
        evaluated = (self.worksheet().sage() is self.sage()) and not self.interrupted()
        if evaluated:
            cls = 'cell_evaluated'
        else:
            cls = 'cell_not_evaluated'

        html_in  = self.html_in(do_print=do_print)
        introspect = "<div id='introspect_div_%s' class='introspection'></div>"%self.id()
        #html_out = self.html_out(wrap, do_print=do_print)
        html_out = self.html()
        s = html_out
        if div_wrap:
            s = '\n\n<div id="cell_outer_%s" class="cell_visible"><div id="cell_%s" class="%s">'%(self.id(), self.id(), cls) + s + '</div></div>'
        return s

    def html(self, wrap=None, div_wrap=True, do_print=False):
        if do_print:
            wrap = 68
            div_wrap = 68
        key = (wrap,div_wrap,do_print)

        if wrap is None:
            wrap = self.notebook().conf()['word_wrap_cols']
        if self.worksheet().compute_process_has_been_started():
            evaluated = (self.worksheet().sage() is self.sage()) and not self.interrupted()
        else:
            evaluated = False
        if evaluated or do_print:
            cls = 'cell_evaluated'
        else:
            cls = 'cell_not_evaluated'

        html_in  = self.html_in(do_print=do_print)
        introspect = "<div id='introspect_div_%s' class='introspection'></div>"%self.id()
        html_out = self.html_out(wrap, do_print=do_print)

        if self.__in.lstrip()[:8] == '%hideall':
            s = html_out
        else:
            s = html_in  + introspect + html_out

        if div_wrap:
            s = '\n\n<div id="cell_outer_%s" class="cell_visible"><div id="cell_%s" class="%s">'%(self.id(), self.id(), cls) + s + '</div></div>'

        #self._html_cache[key] = s
        return s

    def html_in(self, do_print=False, ncols=80):
        s = ''
        id = self.__id
        t = self.__in.rstrip()

        if t.lstrip().startswith('%hide'):
            cls = "cell_input_hide"
        else:
            cls = "cell_input"

##         if do_print:
##             if 'hide' in cls:
##                 return ''
##             else:
##                 s = '<pre class="cell_input">%s</pre>'%(self.__in.replace('<','&lt;'))
##                 return s

        if not do_print:
            s += self.html_new_cell_before()

        #if do_print:
        #    ncols = 70

        r = max(1, number_of_rows(t.strip(), ncols))

        s += """
           <textarea class="%s" rows=%s cols=%s
              id         = 'cell_input_%s'
              onKeyPress = 'return input_keypress(%s,event);'
              onKeyDown  = 'return input_keydown(%s,event);'
              onKeyUp    = 'return cell_input_resize(this);'
              onBlur     = 'cell_blur(%s); return true;'
              onFocus    = 'cell_focused(this,%s); return true;'
              %s
           >%s</textarea>
        """%(cls, r, ncols, id, id, id, id, id, 'readonly=1' if do_print else '', t)

        if not do_print:
           s+= '<a href="javascript:evaluate_cell(%s,0)" class="eval_button" id="eval_button%s" alt="Click here or press shift-return to evaluate">evaluate</a>'%(id,id)

        t = t.replace("<","&lt;")+" "

        #s += """
        #   <pre class="%s"
        #      id         = 'cell_display_%s'
        #      onClick  = 'cell_focus(%s, false); return false;'
        #   >%s</pre>
        #"""%(cls, id, id, t)

        return s

    def html_new_cell_before(self):
        return """<div class="insert_new_cell" id="insert_new_cell_%s"
                   onmousedown="insert_new_cell_before(%s);">
                 </div>
              """%(self.id(), self.id())
    def html_new_cell_after(self):
        return """<div class="insert_new_cell" id="insert_new_cell_%s"
                   onmousedown="insert_new_cell_after(%s);">
                 </div>
              """%(self.id(), self.id())

    def url_to_self(self):
        try:
            return self.__url_to_self
        except AttributeError:
            self.__url_to_self = '/home/%s/cells/%s'%(self.worksheet_filename(), self.id())
            return self.__url_to_self

    def files(self):
        dir = self.directory()
        D = os.listdir(dir)
        return D

    def files_html(self, out):
        D = self.files()
        D.sort()
        if len(D) == 0:
            return ''
        images = []
        files  = []
        # The question mark trick here is so that images will be reloaded when
        # the async request requests the output text for a computation.
        # This is inspired by http://www.irt.org/script/416.htm/.
        for F in D:
            if 'cell://%s'%F in out:
                continue
            url = "%s/%s"%(self.url_to_self(), F)
            if F.endswith('.png') or F.endswith('.bmp') or \
                    F.endswith('.jpg') or F.endswith('.gif'):
                images.append('<img src="%s?%d">'%(url, self.version()))
            elif F.endswith('.obj'):
                images.append("""<a href="javascript:sage3d_show('%s', '%s_%s', '%s');">Click for interactive view.</a>"""%(url, self.__id, F, F[:-4]))
            elif F.endswith('.mtl') or F.endswith(".objmeta"):
                pass # obj data
            elif F.endswith('.svg'):
                images.append('<embed src="%s" type="image/svg+xml" name="emap">'%url)
            elif F.endswith('.jmol'):
                # If F ends in -size500.jmol then we make the viewer applet with size 500.
                i = F.rfind('-size')
                if i != -1:
                    size = F[i+5:-5]
                else:
                    size = 500

                #popup  = """<br><a href="javascript:jmol_popup('%s');">Enlarge</a>"""%url
                #script = '<script>jmol_applet(%s, "%s");</script>%s' % (size, url, popup)
                #script = '<script>jmol_popup("%s");</script>' % (url)

                script = '<div><script>jmol_applet(%s, "%s?%d");</script></div>' % (size, url, self.version())
                images.append(script)
            elif F.endswith('.jmol.zip'):
                pass # jmol data
            else:
                link_text = str(F)
                if len(link_text) > 40:
                    link_text = link_text[:10] + '...' + link_text[-20:]
                files.append('<a target="_new" href="%s" class="file_link">%s</a>'%(url, link_text))
        if len(images) == 0:
            images = ''
        else:
            images = "%s"%'<br>'.join(images)
        if len(files)  == 0:
            files  = ''
        else:
            files  = ('&nbsp'*3).join(files)
        return images + files

    def html_out(self, ncols=0, do_print=False):
        out_nowrap = self.output_text(0, html=True)

        out_html = self.output_html()
        if self.introspect():
            out_wrap = out_nowrap
        else:
            out_wrap = self.output_text(ncols, html=True)

        typ = self.cell_output_type()

        if self.computing():
            cls = "cell_div_output_running"
        else:
            cls = 'cell_div_output_' + typ

        top = '<div class="%s" id="cell_div_output_%s">'%(
                         cls, self.__id)

        if do_print:
            prnt = "print_"
        else:
            prnt = ""

        out_wrap   = '<div class="cell_output_%s%s" id="cell_output_%s">%s</div>'%(
            prnt, typ,self.__id, out_wrap)
        out_nowrap = '<div class="cell_output_%snowrap_%s" id="cell_output_nowrap_%s">%s</div>'%(
            prnt, typ, self.__id, out_nowrap)
        out_html   = '<div class="cell_output_html_%s" id="cell_output_html_%s">%s </div>'%(
            typ, self.__id, out_html)

        out = "%s%s%s"%(out_wrap, out_nowrap, out_html)
        s = top + out + '</div>'

        r = ''
        r += '&nbsp;'*(7-len(r))
        tbl = """
               <table class="cell_output_box"><tr>
               <td class="cell_number" id="cell_number_%s" onClick="cycle_cell_output_type(%s);">
                 %s
               </td>
               <td class="output_cell">%s</td></tr></table>"""%(
                   self.__id, self.__id, r, s)

        return tbl



########

def format_exception(s0, ncols):
    s = s0.lstrip()
    if TRACEBACK not in s:
        return s0
    if ncols > 0:
        s = s.strip()
        w = s.splitlines()
        for k in range(len(w)):
            if TRACEBACK in w[k]:
                break
        s = '\n'.join(w[:k]) + '\nTraceback (click to the left for traceback)' + '\n...\n' + w[-1]
    else:
        s = s.replace("exec compile(ur'","")
        s = s.replace("' + '\\n', '', 'single')", "")
    return s

ComputeCell=Cell


def number_of_rows(txt, ncols):
    rows = txt.splitlines()
    nrows = len(rows)
    for i in range(nrows):
        nrows += int(len(rows[i])/ncols)
    return nrows
