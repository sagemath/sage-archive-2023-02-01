"""
A Cell.

A cell is a single input/output block.  Workbooks are built out of a
list of cells.
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

c = 0

import os, shutil

from   sage.misc.misc import word_wrap

class Cell:
    def __init__(self, id, input, out, workbook):
        self.__id    = int(id)
        self.__in    = str(input)
        self.__out   = str(out)
        self.__workbook = workbook
        self.__dir   = '%s/cells/%s'%(workbook.directory(), self.__id)
        self.__interrupted = False
        self.__completions = False
        self.has_new_output = False

    def __cmp__(self, right):
        return cmp(self.__id, right.__id)

    def __del__(self):
        if os.path.exists(self.__dir):
            shutil.rmtree(self.__dir, ignore_errors=True)

    def __repr__(self):
        return 'Cell %s'%self.__id

    def is_last(self):
        return self.__workbook.cell_list()[-1] == self

    def next_id(self):
        L = self.__workbook.cell_list()
        k = L.index(self)
        if k == len(L):
            return L[0]
        return L[k+1].id()

    def interrupt(self):
        self.__interrupted = True

    def interrupted(self):
        return self.__interrupted

    def computing(self):
        return self in self.__workbook.queue()

    def directory(self):
        if not os.path.exists(self.__dir):
            os.makedirs(self.__dir)
        return self.__dir

    def id(self):
        return self.__id

    def workbook(self):
        return self.__workbook

    def notebook(self):
        return self.__workbook.notebook()

    def set_input_text(self, input):
        self.__in = input

    def input_text(self):
        return self.__in

    def set_output_text(self, output):
        self.__out = output

    def output_text(self, ncols=0):
        if ncols:
            return word_wrap(self.__out, ncols=ncols)
        return self.__out

    def completions(self):
        return self.__completions

    def unset_completions(self):
        self.__completions = False

    def evaluate(self, time=False, completions=False):
        self.__interrupted = False
        self.__time = time
        self.__completions = completions
        self.__workbook.enqueue(self)
        dir = self.directory()
        for D in os.listdir(dir):
            os.unlink(dir + '/' + D)

    def time(self):
        return self.__time

    def html(self, wrap=None, div_wrap=True):
        if wrap is None:
            wrap = self.notebook().defaults()['word_wrap_cols']
        html_in  = self.html_in()
        html_out = self.html_out(wrap)
        s = html_in + html_out
        if div_wrap:
            s = '<div id="cell_%s">'%self.id() + s + '</div>'
        return s

    def html_in(self):
        id = self.__id
        t = self.__in
        r = len(t.split('\n'))
        if r <= 1:
            style = 'style = "height:1.5em"'
        else:
            style = ''
            #     onClick="if (event.shiftKey) {
            #           insert_new_cell_before(%s)}">
            #
        new = """<div class="insert_new_cell" id="insert_new_cell_%s"
                      onClick="insert_new_cell_before(%s)">
                 </div>
              """%(id, id)
        return """%s
           <textarea class="cell_input" rows=%s
              id         = 'cell_input_%s'
              onKeyPress = 'return cell_input_key_event(%s,event);'
              oninput   = 'cell_input_resize(%s);'
              onFocus   = 'cell_input_resize(%s);'
              onBlur    = 'cell_input_minimize_size(%s);'
              %s
           >%s</textarea>
        """%(new, r, id, id, id, id, id, style, t)

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
            images = "<br>%s"%'<br>'.join(images)
        if len(files)  == 0:
            files  = ''
        else:
            files  = ('&nbsp'*3).join(files)
        return images + files

    def html_out(self, ncols=0):
        out_wrap = self.output_text(ncols)
        out_no_wrap = self.output_text(0)
        s = """<div class="cell_output" id="cell_div_output_%s"
                 onClick="cell_output_click(%s, event);">
                 <table class="cell_output"><tr><td>
                 <pre class="cell_output" id="cell_output_%s">%s</pre>
                 <pre class="cell_output_nowrap" id="cell_output_nowrap_%s">%s</pre>
                 </tr></td></table>
               </div>"""%(self.__id, self.__id,
                          self.__id, out_wrap,
                          self.__id, out_no_wrap)
        return s

