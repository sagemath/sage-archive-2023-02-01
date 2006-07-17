r"""
A Cell.

A cell is a single input/output block.  Worksheets are built out of a
list of cells.
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
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
MAX_OUTPUT = 65536



c = 0

import os, shutil

from   sage.misc.misc import word_wrap

import notebook

import worksheet

class Cell:
    def __init__(self, id, input, out, worksheet):
        self.__id    = int(id)
        self.__in    = str(input)
        self.__out   = str(out)
        self.__worksheet = worksheet
        self.__interrupted = False
        self.__completions = False
        self.has_new_output = False
        self.__dir   = '%s/cells/%s'%(worksheet.directory(), self.relative_id())

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
        self.__dir = '%s/cells/%s'%(worksheet.directory(), self.relative_id())
        if not id is None:
            self.set_id(id)
        self.__out_html = self.files_html()

    def __cmp__(self, right):
        return cmp(self.__id, right.__id)

    def __del__(self):
        if os.path.exists(self.__dir):
            shutil.rmtree(self.__dir, ignore_errors=True)

    def __repr__(self):
        return 'Cell %s; in=%s, out=%s'%(self.__id, self.__in, self.__out)

    def plain_text(self, ncols=0, prompts=True, max_out=None):
        if ncols == 0:
            ncols = self.notebook().defaults()['word_wrap_cols']
        s = ''

        input_lines = self.__in.strip()
        if input_lines[:1] == '%':
            pr = '%s> '%(input_lines.split()[0])[1:]
        else:
            pr = 'sage: '

        if prompts:
            input_lines = input_lines.split('\n')
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
                    #    s += '<BLANKLINE>\n'
                    elif v[0] == ' ':
                        in_loop = True
                        s += '...' + v + '\n'
                    else:
                        if in_loop:
                            s += '...\n'
                            in_loop = False
                        s += pr + v + '\n'
        else:
            s += self.__in.strip() + '\n'

        if prompts:
            msg = 'Traceback (most recent call last):'
            if self.__out[:len(msg)] == msg:
                v = self.__out.split('\n')
                w = [v[0], '...']
                for i in range(1,len(v)):
                    if not (len(v[i]) > 0 and v[i][0] == ' '):
                        w = w + v[i:]
                        break
                out = '\n'.join(w)
            else:
                out = self.output_text(ncols, html=False)
        else:
            out = self.output_text(ncols, html=False).strip().split('\n')
            out = [x for x in out if x.strip() != '']
            if len(out) > 0:
                out = '# ' + '\n# '.join(out)
            else:
                out = ''

        if not max_out is None and len(out) > max_out:
            out = out[:max_out] + '...'

        s = s.strip() + '\n' + out.strip()

        return s

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

    def directory(self):
        if not os.path.exists(self.__dir):
            os.makedirs(self.__dir)
        return self.__dir

    def id(self):
        return self.__id

    def relative_id(self):
        return self.__id - self.__worksheet.id()*notebook.MAX_WORKSHEETS

    def set_id(self, id):
        self.__id = int(id)

    def worksheet(self):
        return self.__worksheet

    def notebook(self):
        return self.__worksheet.notebook()

    def set_input_text(self, input):
        self.__in = input

    def input_text(self):
        return self.__in

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
        i = output.find(worksheet.SAGE_VARS)
        if i != -1:
            output = output[:i]
        if len(output) > MAX_OUTPUT:
            output = 'WARNING: Output truncated!\n' + output[:MAX_OUTPUT] + '\n(truncated)'
        self.__out = output.strip()
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

    def output_text(self, ncols=0, html=True):
        s = self.__out

        def format(x):
            return word_wrap(x.replace('<','&lt;'), ncols=ncols)

        if html:
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
                    t += format(s)
                    break
                t += format(s[:i]) + s[i+6:j]
                s = s[j+7:]
            s = t
            if not self.is_html() and len(s.strip()) > 0:
                s = '<pre class="shrunk">' + s + '</pre>'
        return s

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

    def evaluate(self, introspect=False, time=False):
        """
        INPUT:
            time -- if True return time computation takes
            introspect -- either False or a pair [before_curse, after_cursor] of strings.
        """
        self.__interrupted = False
        self.__time = time
        self.__introspect = introspect
        self.__worksheet.enqueue(self)
        self.__type = 'wrap'
        dir = self.directory()
        for D in os.listdir(dir):
            os.unlink(dir + '/' + D)

    def time(self):
        return self.__time

    def do_time(self):
        self.__time = True

    def html(self, wrap=None, div_wrap=True, do_print=False):
        if wrap is None:
            wrap = self.notebook().defaults()['word_wrap_cols']
        evaluated = (self.worksheet().sage() is self.sage()) and not self.interrupted()
        if evaluated:
            cls = 'cell_evaluated'
        else:
            cls = 'cell_not_evaluated'

        html_in  = self.html_in(do_print=do_print)
        introspect = "<div id='introspect_div_%s' class='introspection'></div>"%self.id()
        html_out = self.html_out(wrap, do_print=do_print)
        s = html_in  + introspect + html_out
        if div_wrap:
            s = '\n\n<div id="cell_%s" class="%s">'%(self.id(), cls) + s + '</div>'
        return s

    def html_in(self, do_print=False):
        id = self.__id
        t = self.__in.rstrip()

        if t.lstrip()[:5] == '%hide':
            cls = "cell_input_hide"
        else:
            cls = "cell_input"

        if do_print:
            if 'hide' in cls:
                return ''
            else:
                s = '<pre class="shrunk">%s</pre>'%(self.__in.replace('<','&lt;'))
                return s

        s = """<div class="insert_new_cell" id="insert_new_cell_%s"
                   onmousedown="insert_new_cell_before(%s);">
                 </div>
              """%(id, id)

        r = len(t.split('\n'))

        s += """
           <textarea class="%s" rows=%s cols=100000 columns=100000
              id         = 'cell_input_%s'
              onKeyPress = 'return input_keypress(%s,event);'
              oninput   = 'cell_input_resize(%s);'
v v v v v v v
              onFocus = 'cell_focus(%s)'
              onBlur  = 'cell_blur(%s);'
*************
              onFocus = 'return cell_focus(%s)'
              onBlur  = 'return cell_blur(%s)'
^ ^ ^ ^ ^ ^ ^
           >%s</textarea>
        """%(cls, r, id, id, id, id, id, t)
        return s

    def files_html(self):
        dir = self.directory()
        D = os.listdir(dir)
        if len(D) == 0:
            return ''
        images = []
        files  = []
        # The c and question mark hack here is so that images will be reloaded when
        # the async request requests the output text for a computation.
        # This is a total hack, inspired by http://www.irt.org/script/416.htm/.
        global c
        c += 1
        for F in D:
            if F[-4:] == '.png':
                images.append('<img src="%s/%s?%s">'%(dir,F,c))
            elif F[-4:] == '.svg':
                images.append('<embed src="%s/%s" type="image/svg+xml" name="emap">'%(dir,F))
            else:
                files.append('<a href="%s/%s">%s</a>'%(dir, F, F))
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
        out_nowrap = self.output_text(0, html=True).strip()
        if self.introspect():
            out_wrap = out_nowrap
        else:
            out_wrap = self.output_text(ncols, html=True).strip()

        typ = self.cell_output_type()

        if self.computing():
            cls = "cell_output_running"
        else:
            cls = 'cell_output_' + typ

        top = '<div class="%s" id="cell_div_output_%s">'%(
                         cls, self.__id)

        out_html = self.output_html()

        out = """<span class="cell_output_%s" id="cell_output_%s">%s </span>
                 <span class="cell_output_nowrap_%s" id="cell_output_nowrap_%s">%s </span>
                 <span class="cell_output_html_%s" id="cell_output_html_%s">%s </span>
                 """%(typ, self.__id, out_wrap,
                      typ, self.__id, out_nowrap,
                      typ, self.__id, out_html)


        s = top + out + '\n</div>'

        r = '[%s]'%self.relative_id()
        r += '&nbsp;'*(5-len(r))
        tbl = """<table class="cell_output_box"><tr>
               <td class="cell_number" id="cell_number_%s" onClick="cycle_cell_output_type(%s);">%s</td>
               <td class="output_cell">%s</td></tr></table>"""%(
                   self.__id, self.__id, r, s)

        return tbl

