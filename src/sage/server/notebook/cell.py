"""
A Cell.

A cell is a single input/output block. Worksheets are built out of
a list of cells.
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
from   sage.misc.package   import is_package_installed

from cgi import escape

if is_package_installed("tinyMCE"):
    JEDITABLE_TINYMCE = True
else:
    JEDITABLE_TINYMCE = False


class Cell_generic:
    def is_interactive_cell(self):
        """
        Returns True if this cell contains the use of interact either as a
        function call or a decorator.

        EXAMPLES::

            sage: from sage.server.notebook.cell import Cell_generic
            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: Cell_generic.is_interactive_cell(C)
            False
        """
        return False

    def delete_output(self):
        """
        Delete all output in this cell. This is not executed - it is an
        abstract function that must be overwritten in a derived class.

        EXAMPLES: This function just raises a NotImplementedError, since it
        most be defined in derived class.

        ::

            sage: C = sage.server.notebook.cell.Cell_generic()
            sage: C.delete_output()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class TextCell(Cell_generic):
    def __init__(self, id, text, worksheet):
        """
        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C == loads(dumps(C))
            True
        """
        self.__id = int(id)
        self.__text = text
        self.__worksheet = worksheet

    def __repr__(self):
        """
        String representation of this text cell.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C.__repr__()
            'TextCell 0: 2+3'
        """
        return "TextCell %s: %s"%(self.__id, self.__text)

    def delete_output(self):
        """
        Delete all output in this cell. This does nothing since text cells
        have no output.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C
            TextCell 0: 2+3
            sage: C.delete_output()
            sage: C
            TextCell 0: 2+3
        """
        pass # nothing to do -- text cells have no output

    def set_input_text(self, input_text):
        """
        Sets the input text of self to be input_text.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C
            TextCell 0: 2+3
            sage: C.set_input_text("3+2")
            sage: C
            TextCell 0: 3+2
        """
        self.__text = input_text

    def set_worksheet(self, worksheet, id=None):
        """
        Sets the worksheet object of self to be worksheet and optionally
        changes the id of self.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: W = "worksheet object"
            sage: C.set_worksheet(W)
            sage: C.worksheet()
            'worksheet object'
            sage: C.set_worksheet(None, id=2)
            sage: C.id()
            2
        """
        self.__worksheet = worksheet
        if id is not None:
            self.__id = id

    def worksheet(self):
        """
        Returns the worksheet object associated to self.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', 'worksheet object')
            sage: C.worksheet()
            'worksheet object'
        """
        return self.__worksheet

    def html(self, ncols=0, do_print=False, do_math_parse=True, editing=False):
        """
        Returns an HTML version of self as a string.

        INPUT:

        - ``do_math_parse`` - bool (default: True)
          If True, call math_parse (defined in cell.py)
          on the html.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C.html()
            '<div class="text_cell" id="cell_text_0">2+3...'
            sage: C.set_input_text("$2+3$")
            sage: C.html(do_math_parse=True)
            '<div class="text_cell" id="cell_text_0"><span class="math">2+3</span>...'
        """

        s = """<div class="text_cell" id="cell_text_%s">%s</div>"""%(self.__id,self.html_inner(ncols=ncols, do_print=do_print, do_math_parse=do_math_parse, editing=editing))

        if JEDITABLE_TINYMCE and hasattr(self.worksheet(),'is_published') and not self.worksheet().is_published():
            s += """<script>$("#cell_text_%s").unbind('dblclick').editable(function(value,settings) {
evaluate_text_cell_input(%s,value,settings);
return(value);
}, {
      tooltip   : "",
      placeholder : "",
//      type   : 'textarea',
      type   : 'mce',
      onblur : 'ignore',
      select : false,
      submit : 'Save changes',
      cancel : 'Cancel changes',
      event  : "dblclick",
      style  : "inherit",
      data   : %r
  });
</script>"""%(self.__id,self.__id,self.__text)


        if editing:
            s += """<script>$("#cell_text_%s").trigger('dblclick');</script>"""%self.__id

        return s

    def html_inner(self,ncols=0, do_print=False, do_math_parse=True, editing=False):
        """
        Returns an HTML version of the content of self as a string.

        INPUT:

        - ``do_math_parse`` - bool (default: True)
          If True, call math_parse (defined in cell.py)
          on the html.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C.html_inner()
            '2+3...'
            sage: C.set_input_text("$2+3$")
            sage: C.html_inner(do_math_parse=True)
            '<span class="math">2+3</span>...'
        """
        t = self.__text
        if do_math_parse:
            # Do dollar sign math parsing
            t = math_parse(t)
        s = """%s"""%t

        return s


    def plain_text(self, prompts=False):
        """
        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C.plain_text()
            '2+3'
        """
        return self.__text

    def edit_text(self):
        """
        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C.edit_text()
            '2+3'
        """
        return self.__text

    def id(self):
        """
        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C.id()
            0
        """
        return self.__id

    def is_auto_cell(self):
        """
        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C.is_auto_cell()
            False
        """
        return False

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: C1 = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C2 = sage.server.notebook.cell.TextCell(0, '3+2', None)
            sage: C3 = sage.server.notebook.cell.TextCell(1, '2+3', None)
            sage: C1 == C1
            True
            sage: C1 == C2
            True
            sage: C1 == C3
            False
        """
        return cmp(self.id(), right.id())

    def set_cell_output_type(self, typ='wrap'):
        """
        This does nothing for TextCells.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.TextCell(0, '2+3', None)
            sage: C.set_cell_output_type("wrap")
        """
        pass # ignored


class Cell(Cell_generic):
    def __init__(self, id, input, out, worksheet):
        """
        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C == loads(dumps(C))
            True
        """
        self.__id    = int(id)
        self.__out   = str(out).replace('\r','')
        self.__worksheet = worksheet
        self.__interrupted = False
        self.__completions = False
        self.has_new_output = False
        self.__no_output_cell = False
        self.__asap = False
        self.__version = -1
        self.set_input_text(str(input).replace('\r',''))

    def set_asap(self, asap):
        """
        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.is_asap()
            False
            sage: C.set_asap(True)
            sage: C.is_asap()
            True
        """
        self.__asap = bool(asap)

    def is_asap(self):
        """
        Return True if this is an asap cell, i.e., evaluation of it is done
        as soon as possible.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.is_asap()
            False
            sage: C.set_asap(True)
            sage: C.is_asap()
            True
        """
        try:
            return self.__asap
        except AttributeError:
            self.__asap = False
            return self.__asap

    def delete_output(self):
        """
        Delete all output in this cell.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None); C
            Cell 0; in=2+3, out=5
            sage: C.delete_output()
            sage: C
            Cell 0; in=2+3, out=
        """
        self.__out = ''
        self.__out_html = ''
        self.__evaluated = False

    def evaluated(self):
        r"""
        Return True if this cell has been successfully evaluated in a
        currently running session.

        This is not about whether the output of the cell is valid given the
        input.

        OUTPUT:


        -  ``bool`` - whether or not this cell has been
           evaluated in this session


        EXAMPLES: We create a worksheet with a cell that has wrong output::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: W.edit_save('Sage\n{{{\n2+3\n///\n20\n}}}')
            sage: C = W.cell_list()[0]
            sage: C
            Cell 0; in=2+3, out=20

        We re-evaluate that input cell::

            sage: C.evaluate()
            sage: W.check_comp(wait=9999)
            ('d', Cell 0; in=2+3, out=
            5
            )

        Now the output is right::

            sage: C
            Cell 0; in=2+3, out=
            5

        And the cell is considered to have been evaluated.

        ::

            sage: C.evaluated()
            True

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        # Cells are never considered evaluated in a new session.
        if not self.worksheet().compute_process_has_been_started():
            self.__evaluated = False
            return False

        # Figure out if the worksheet is using the same sage
        # session as this cell.  (I'm not sure when this would
        # be False.)
        same_session = self.worksheet().sage() is self.sage()
        try:
            # Always not evaluated if sessions are different.
            if not same_session:
                self.__evaluated = False
                return False
            return self.__evaluated
        except AttributeError:
            # Default assumption is that cell has not been evaluated.
            self.__evaluated = False
            return False

    def set_no_output(self, no_output):
        """
        Sets whether or not this is an no_output cell, i.e., a cell for
        which we don't care at all about the output.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.is_no_output()
            False
            sage: C.set_no_output(True)
            sage: C.is_no_output()
            True
        """
        self.__no_output = bool(no_output)

    def is_no_output(self):
        """
        Return True if this is an no_output cell, i.e., a cell for which
        we don't care at all about the output.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.is_no_output()
            False
            sage: C.set_no_output(True)
            sage: C.is_no_output()
            True
        """
        try:
            return self.__no_output
        except AttributeError:
            self.__no_output = False
            return self.__no_output

    def set_cell_output_type(self, typ='wrap'):
        """
        Sets the cell output type.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.cell_output_type()
            'wrap'
            sage: C.set_cell_output_type('nowrap')
            sage: C.cell_output_type()
            'nowrap'
        """
        self.__type = typ

    def cell_output_type(self):
        """
        Returns the cell output type.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.cell_output_type()
            'wrap'
            sage: C.set_cell_output_type('nowrap')
            sage: C.cell_output_type()
            'nowrap'
        """
        try:
            return self.__type
        except AttributeError:
            self.__type = 'wrap'
            return self.__type

    def set_worksheet(self, worksheet, id=None):
        """
        Sets the worksheet object of self to be worksheet and optionally
        changes the id of self.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: W = "worksheet object"
            sage: C.set_worksheet(W)
            sage: C.worksheet()
            'worksheet object'
            sage: C.set_worksheet(None, id=2)
            sage: C.id()
            2
        """
        self.__worksheet = worksheet
        if id is not None:
            self.set_id(id)

    def worksheet(self):
        """
        Returns the worksheet object associated to self.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', 'worksheet object')
            sage: C.worksheet()
            'worksheet object'
        """
        return self.__worksheet

    def update_html_output(self, output=''):
        """
        Update the list of files with html-style links or embeddings for
        this cell.

        For interactive cells the html output section is always empty,
        mainly because there is no good way to distinguish content (e.g.,
        images in the current directory) that goes into the interactive
        template and content that would go here.
        """
        if self.is_interactive_cell():
            self.__out_html = ""
        else:
            self.__out_html = self.files_html(output)

    def id(self):
        """
        Returns the id of self.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.id()
            0
        """
        return self.__id

    def set_id(self, id):
        """
        Sets the id of self to id.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.set_id(2)
            sage: C.id()
            2
        """
        self.__id = int(id)

    def worksheet(self):
        """
        Returns the workseet associated to self.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', W)
            sage: C.worksheet() is W
            True

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        return self.__worksheet

    def worksheet_filename(self):
        """
        Returns the filename of the worksheet associated to self.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', W)
            sage: C.worksheet_filename()
            'sage/0'

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        return self.__worksheet.filename()


    def notebook(self):
        """
        Returns the notebook object associated to self.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', W)
            sage: C.notebook() is nb
            True

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        return self.__worksheet.notebook()

    def directory(self):
        """
        Returns the directory associated to self. If the directory doesn't
        already exist, then this method creates it.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', W)
            sage: C.directory()
            '.../worksheets/sage/0/cells/0'

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        dir = self._directory_name()
        if not os.path.exists(dir):
            os.makedirs(dir)
        return dir

    def _directory_name(self):
        """
        Returns a string of the directory associated to self.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', W)
            sage: C._directory_name()
            '.../worksheets/sage/0/cells/0'

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        return '%s/cells/%s'%(self.__worksheet.directory(), self.id())


    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: C1 = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C2 = sage.server.notebook.cell.Cell(0, '3+2', '5', None)
            sage: C3 = sage.server.notebook.cell.Cell(1, '2+3', '5', None)
            sage: C1 == C1
            True
            sage: C1 == C2
            True
            sage: C1 == C3
            False
        """
        return cmp(self.id(), right.id())

    def __repr__(self):
        """
        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None); C
            Cell 0; in=2+3, out=5
        """
        return 'Cell %s; in=%s, out=%s'%(self.__id, self.__in, self.__out)

    def word_wrap_cols(self):
        """
        Returns the number of columns for word wrapping. This defaults to
        70, but the default setting for a notebook is 72.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.word_wrap_cols()
            70

        ::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', W)
            sage: C.word_wrap_cols()
            72

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        try:
            return self.notebook().conf()['word_wrap_cols']
        except AttributeError:
            return 70

    def plain_text(self, ncols=0, prompts=True, max_out=None):
        """
        Returns the plain text version of self.

        TODO: Add more comprehensive doctests.
        """
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
        r"""
        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.edit_text()
            '{{{id=0|\n2+3\n///\n5\n}}}'
        """
        s = self.plain_text(ncols, prompts, max_out)
        return '{{{id=%s|\n%s\n}}}'%(self.id(), s)

    def is_last(self):
        """
        Returns True if self is the last cell in the worksheet.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = W.new_cell_after(0, "2^2"); C
            Cell 1; in=2^2, out=
            sage: C.is_last()
            True
            sage: C = W.get_cell_with_id(0)
            sage: C.is_last()
            False

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        return self.__worksheet.cell_list()[-1] == self

    def next_id(self):
        """
        Returns the id of the next cell in the worksheet associated to
        self. If self is not in the worksheet or self is the last cell in
        the cell_list, then the id of the first cell is returned.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = W.new_cell_after(0, "2^2")
            sage: C = W.get_cell_with_id(0)
            sage: C.next_id()
            1
            sage: C = W.get_cell_with_id(1)
            sage: C.next_id()
            0

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
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
        """
        Record that the calculation running in this cell was interrupted.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = W.new_cell_after(0, "2^2")
            sage: C.interrupt()
            sage: C.interrupted()
            True
            sage: C.evaluated()
            False

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        self.__interrupted = True
        self.__evaluated = False

    def interrupted(self):
        """
        Returns True if the evaluation of this cell has been interrupted.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = W.new_cell_after(0, "2^2")
            sage: C.interrupt()
            sage: C.interrupted()
            True

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        return self.__interrupted

    def computing(self):
        """
        Returns True if self is in its worksheet's queue.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = W.new_cell_after(0, "2^2")
            sage: C.computing()
            False

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        return self in self.__worksheet.queue()

    def is_interactive_cell(self):
        r"""
        Return True if this cell contains the use of interact either as a
        function call or a decorator.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = W.new_cell_after(0, "@interact\ndef f(a=slider(0,10,1,5):\n    print a^2")
            sage: C.is_interactive_cell()
            True
            sage: C = W.new_cell_after(C.id(), "2+2")
            sage: C.is_interactive_cell()
            False

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        # Do *not* cache
        s = strip_string_literals(self.input_text())[0]
        return bool(re.search('(?<!\w)interact\s*\(.*\).*', s) or re.search('\s*@\s*interact\s*\n', s))

    def is_interacting(self):
        """
        Returns True
        """
        return hasattr(self, 'interact')

    def stop_interacting(self):
        """

        """
        if self.is_interacting():
            del self.interact

    def set_input_text(self, input):
        """
        Sets the input text of self to be the string input.

        TODO: Add doctests for the code dealing with interact.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = W.new_cell_after(0, "2^2")
            sage: C.evaluate()
            sage: W.check_comp(wait=9999)
            ('d', Cell 1; in=2^2, out=
            4
            )
            sage: C.version()
            0

        ::

            sage: C.set_input_text('3+3')
            sage: C.input_text()
            '3+3'
            sage: C.evaluated()
            False
            sage: C.version()
            1

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        # Stuff to deal with interact
        if input.startswith('%__sage_interact__'):
            self.interact = input[len('%__sage_interact__')+1:]
            self.__version = self.version() + 1
            return
        elif self.is_interacting():
            try:
                del self.interact
                del self._interact_output
            except AttributeError:
                pass

        # We have updated the input text so the cell can't have
        # been evaluated.
        self.__evaluated = False
        self.__version = self.version() + 1
        self.__in = input
        if hasattr(self, '_html_cache'):
            del self._html_cache

        #Run get the input text with all of the percent
        #directives parsed
        self._cleaned_input = self.parse_percent_directives()

    def input_text(self):
        """
        Returns self's input text.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.input_text()
            '2+3'
        """
        return self.__in

    def cleaned_input_text(self):
        r"""
        Returns the input text with all of the percent directives
        removed.  If the cell is interacting, then the interacting
        text is returned.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '%hide\n%maxima\n2+3', '5', None)
            sage: C.cleaned_input_text()
            '2+3'

        """
        if self.is_interacting():
            return self.interact
        else:
            return self._cleaned_input

    def parse_percent_directives(self):
        r"""
        Returns a string which consists of the input text of this cell
        with the percent directives at the top removed.  As it's doing
        this, it computes a list of all the directives and which
        system (if any) the cell should be run under.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '%hide\n%maxima\n2+3', '5', None)
            sage: C.parse_percent_directives()
            '2+3'
            sage: C.percent_directives()
            ['hide', 'maxima']

        """
        self._system = None
        text = self.input_text().split('\n')
        directives = []
        for i, line in enumerate(text):
            if not line.startswith('%'):
                #Handle the #auto case here for now
                if line == "#auto":
                    pass
                else:
                    break
            elif line in ['%auto', '%hide', '%hideall', '%save_server', "%time", "%timeit"]:
                #We do not consider any of the above percent
                #directives as specifying a system.
                pass
            else:
                self._system = line[1:]

            directives.append(line[1:])

        self._percent_directives = directives
        return "\n".join(text[i:]).strip()

    def percent_directives(self):
        r"""
        Returns a list of all the percent directives that appear
        in this cell.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '%hide\n%maxima\n2+3', '5', None)
            sage: C.percent_directives()
            ['hide', 'maxima']

        """
        return self._percent_directives

    def system(self):
        r"""
        Returns the system used to evaluate this cell. The system
        is specified by a percent directive like '%maxima' at
        the top of a cell.

        If no system is explicitly specified, then None is returned
        which tells the notebook to evaluate the cell using the
        worksheet's default system.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '%maxima\n2+3', '5', None)
            sage: C.system()
            'maxima'
            sage: prefixes = ['%hide', '%time', '']
            sage: cells = [sage.server.notebook.cell.Cell(0, '%s\n2+3'%prefix, '5', None) for prefix in prefixes]
            sage: [(C, C.system()) for C in cells if C.system() is not None]
            []
        """
        return self._system


    def is_auto_cell(self):
        r"""
        Returns True if self is an auto cell.

        An auto cell is a cell that is automatically evaluated when the
        worksheet starts up.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.is_auto_cell()
            False
            sage: C = sage.server.notebook.cell.Cell(0, '#auto\n2+3', '5', None)
            sage: C.is_auto_cell()
            True
        """
        return 'auto' in self.percent_directives()

    def changed_input_text(self):
        """
        Returns the changed input text for the cell. If there was any
        changed input text, then it is reset to " before this method
        returns.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.changed_input_text()
            ''
            sage: C.set_changed_input_text('3+3')
            sage: C.input_text()
            '3+3'
            sage: C.changed_input_text()
            '3+3'
            sage: C.changed_input_text()
            ''
            sage: C.version()
            0
        """
        try:
            t = self.__changed_input
            del self.__changed_input
            return t
        except AttributeError:
            return ''

    def set_changed_input_text(self, new_text):
        """
        Note that this does not update the version of the cell. This is
        typically used for things like tab completion.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.set_changed_input_text('3+3')
            sage: C.input_text()
            '3+3'
            sage: C.changed_input_text()
            '3+3'
        """
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
        # We do not truncate if "notruncate" or "Output truncated!" already
        # appears in the output.  This notruncate tag is used right now
        # in sage.server.support.help.
        if 'notruncate' not in output and 'Output truncated!' not in output and \
               (len(output) > MAX_OUTPUT or output.count('\n') > MAX_OUTPUT_LINES):
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
        """
        TODO: Figure out what exactly this does.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.sage() is None
            True
        """
        try:
            return self.__sage
        except AttributeError:
            return None

    def output_html(self):
        try:
            return self.__out_html
        except AttributeError:
            self.__out_html = ''
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
            return word_wrap(escape(x), ncols=ncols)

        def format_html(x):
            return self.process_cell_urls(x)

        # if there is an error in the output,
        # specially format it.
        if not self.is_interactive_cell():
            s = format_exception(format_html(s), ncols)

        # Everything not wrapped in <html> ... </html> should be
        # escaped and word wrapped.
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

        # Get rid of the <script> tags, since we do not want them to
        # be evaluated twice.  They are only evaluated in the wrapped
        # version of the output.
        if ncols == 0:
            while True:
                i = t.lower().find('<script>')
                if i == -1: break
                j = t[i:].lower().find('</script>')
                if j == -1: break
                t = t[:i] + t[i+j+len('</script>'):]

        return t


    def has_output(self):
        """
        Returns True if there is output for this cell.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.has_output()
            True
            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '', None)
            sage: C.has_output()
            False
        """
        return len(self.__out.strip()) > 0

    def is_html(self):
        r"""
        Returns True if this is an HTML cell. An HTML cell whose system is
        'html' and is typically specified by ``%html``.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, "%html\nTest HTML", None, None)
            sage: C.system()
            'html'
            sage: C.is_html()
            True
            sage: C = sage.server.notebook.cell.Cell(0, "Test HTML", None, None)
            sage: C.is_html()
            False

        """
        try:
            return self.__is_html
        except AttributeError:
            return self.system() == 'html'

    def set_is_html(self, v):
        """
        Sets whether or not this cell is an HTML cell.

        This is called by check_for_system_switching in worksheet.py.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.is_html()
            False
            sage: C.set_is_html(True)
            sage: C.is_html()
            True
        """
        self.__is_html = v

    #################
    # Introspection #
    #################
    def set_introspect_html(self, html, completing=False):
        if completing:
            self.__introspect_html = html
        else:
            html = escape(html).strip()
            self.__introspect_html = '<pre class="introspection">'+html+'</pre>'

    def introspect_html(self):
        if not self.introspect():
            return ''
        try:
            return self.__introspect_html
        except AttributeError:
            self.__introspect_html = ''
            return ''

    def introspect(self):
        """
        TODO: Figure out what the __introspect method is for and write a
        better doctest.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.introspect()
            False
            sage: C.set_introspect("a", "b")
            sage: C.introspect()
            ['a', 'b']
        """
        try:
            return self.__introspect
        except AttributeError:
            return False

    def unset_introspect(self):
        """
        TODO: Figure out what the __introspect method is for and write a
        better doctest.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.set_introspect("a", "b")
            sage: C.introspect()
            ['a', 'b']
            sage: C.unset_introspect()
            sage: C.introspect()
            False
        """
        self.__introspect = False

    def set_introspect(self, before_prompt, after_prompt):
        """
        TODO: Figure out what the __introspect method is for and write a
        better doctest.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.set_introspect("a", "b")
            sage: C.introspect()
            ['a', 'b']
        """
        self.__introspect = [before_prompt, after_prompt]

    def evaluate(self, introspect=False, time=None, username=None):
        r"""
        INPUT:


        -  ``username`` - name of user doing the evaluation

        -  ``time`` - if True return time computation takes

        -  ``introspect`` - either False or a pair
           [before_cursor, after_cursor] of strings.


        EXAMPLES: We create a notebook, worksheet, and cell and evaluate it
        in order to compute `3^5`::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: W.edit_save('Sage\n{{{\n3^5\n}}}')
            sage: C = W.cell_list()[0]; C
            Cell 0; in=3^5, out=
            sage: C.evaluate(username='sage')
            sage: W.check_comp(wait=9999)
            ('d', Cell 0; in=3^5, out=
            243
            )
            sage: C
            Cell 0; in=3^5, out=
            243

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        self.__interrupted = False
        self.__evaluated = True
        if time is not None:
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
        """
        Returns the version number of this cell.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.version()
            0
            sage: C.set_input_text('2+3')
            sage: C.version()
            1
        """
        try:
            return self.__version
        except AttributeError:
            self.__version = 0
            return self.__version

    def time(self):
        r"""
        Returns True if the time it takes to evaluate this cell should be
        printed.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: C.time()
            False
            sage: C = sage.server.notebook.cell.Cell(0, '%time\n2+3', '5', None)
            sage: C.time()
            True
        """
        return ('time' in self.percent_directives() or
                'timeit' in self.percent_directives() or
                getattr(self, '__time', False))

    def doc_html(self, wrap=None, div_wrap=True, do_print=False):
        """
        Modified version of ``self.html`` for the doc browser.
        This is a hack and needs to be improved. The problem is how to get
        the documentation html to display nicely between the example cells.
        The type setting (jsMath formating) needs attention too.
        """
        self.evaluate()
        if wrap is None:
            wrap = self.notebook().conf()['word_wrap_cols']
        evaluated = self.evaluated()
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
        evaluated = self.evaluated()
        if evaluated or do_print:
            cls = 'cell_evaluated'
        else:
            cls = 'cell_not_evaluated'

        html_in  = self.html_in(do_print=do_print)
        introspect = "<div id='introspect_div_%s' class='introspection'></div>"%self.id()
        html_out = self.html_out(wrap, do_print=do_print)

        if 'hideall' in self.percent_directives():
            s = html_out
        else:
            s = html_in + introspect + html_out

        if div_wrap:
            s = '\n\n<div id="cell_outer_%s" class="cell_visible"><div id="cell_%s" class="%s">'%(self.id(), self.id(), cls) + s + '</div></div>'

        #self._html_cache[key] = s
        return s

    def html_in(self, do_print=False, ncols=80):
        """
        Returns the HTML code for the input of this cell.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: print C.html_in()
            <div class="insert_new_cell" id="insert_new_cell_0"...</a>
        """
        s = ''
        id = self.__id
        t = self.__in.rstrip()

        cls = "cell_input_hide" if 'hide' in self.percent_directives() else "cell_input"

        if not do_print:
            s += self.html_new_cell_before()

        r = max(1, number_of_rows(t.strip(), ncols))

        if do_print:
            tt = escape(t).replace('\n','<br>').replace('  ',' &nbsp;') + '&nbsp;'
            s += '<div class="cell_input_print">%s</div>'%tt
        else:
            s += """
               <textarea class="%s" rows=%s cols=%s
                  id         = 'cell_input_%s'
                  onKeyPress = 'return input_keypress(%s,event);'
                  onKeyDown  = 'return input_keydown(%s,event);'
                  onKeyUp    = 'return input_keyup(%s, event);'
                  onBlur     = 'cell_blur(%s); return true;'
                  onFocus    = 'cell_focused(this,%s); return true;'
               >%s</textarea>
            """%(cls, r, ncols, id, id, id, id, id, id, t)

        if not do_print:
           s+= '<a href="javascript:evaluate_cell(%s,0)" class="eval_button" id="eval_button%s" alt="Click here or press shift-return to evaluate">evaluate</a>'%(id,id)

        t = escape(t)+" "

        return s

    def html_new_cell_before(self):
        """
        Returns the HTML code for inserting a new cell before self.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: print C.html_new_cell_before()
            <div class="insert_new_cell" id="insert_new_cell_0">...
        """
        return """<div class="insert_new_cell" id="insert_new_cell_%(id)s">
                 </div>
<script type="text/javascript">
$("#insert_new_cell_%(id)s").plainclick(function(e) {insert_new_cell_before(%(id)s);});
$("#insert_new_cell_%(id)s").shiftclick(function(e) {insert_new_text_cell_before(%(id)s);});
</script>"""%{'id': self.id()}
    def html_new_cell_after(self):
        """
        Returns the HTML code for inserting a new cell after self.

        EXAMPLES::

            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', None)
            sage: print C.html_new_cell_after()
            <div class="insert_new_cell" id="insert_new_cell_0">...
        """
        return """<div class="insert_new_cell" id="insert_new_cell_%(id)s">
                 </div>
<script type="text/javascript">
$("#insert_new_cell_%(id)s").plainclick(function(e) {insert_new_cell_after(%(id)s);});
$("#insert_new_cell_%(id)s").shiftclick(function(e) {insert_new_text_cell_after(%(id)s);});
</script>"""%{'id': self.id()}


    def url_to_self(self):
        """
        Returns a notebook URL for this cell.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, '2+3', '5', W)
            sage: C.url_to_self()
            '/home/sage/0/cells/0'

        """
        try:
            return self.__url_to_self
        except AttributeError:
            self.__url_to_self = '/home/%s/cells/%s'%(self.worksheet_filename(), self.id())
            return self.__url_to_self

    def files(self):
        """
        Returns a list of all the files in self's directory.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, 'plot(sin(x),0,5)', ", W)
            sage: C.evaluate()
            sage: W.check_comp(wait=9999)
            ('d', Cell 0; in=plot(sin(x),0,5), out=
            <BLANKLINE>
            )
            sage: C.files()
            ['sage0.png']

        ::

            sage: import shutil; shutil.rmtree(nb.directory())
        """
        dir = self.directory()
        D = os.listdir(dir)
        return D

    def delete_files(self):
        """
        Deletes all of the files associated with this cell.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: C = sage.server.notebook.cell.Cell(0, 'plot(sin(x),0,5)', ", W)
            sage: C.evaluate()
            sage: W.check_comp(wait=9999)
            ('d', Cell 0; in=plot(sin(x),0,5), out=
            <BLANKLINE>
            )
            sage: C.files()
            ['sage0.png']
            sage: C.delete_files()
            sage: C.files()
            []
        """
        try:
            dir = self._directory_name()
        except AttributeError:
            return
        if os.path.exists(dir):
            shutil.rmtree(dir, ignore_errors=True)



    def files_html(self, out):
        import time
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
                images.append('<img src="%s?%d">'%(url, time.time()))
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

                script = '<div><script>jmol_applet(%s, "%s?%d");</script></div>' % (size, url, time.time())
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
        if do_print and self.cell_output_type() == 'hidden':
            return '<pre>\n</pre>'

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
    r"""
    Make it so excpetions don't appear expanded by default.

    INPUT:


    -  ``s0`` - string

    -  ``ncols`` - integer


    OUTPUT: string

    If s0 contains "notracebacks" then this function always returns s0

    EXAMPLES::

        sage: sage.server.notebook.cell.format_exception(sage.server.notebook.cell.TRACEBACK,80)
        '\nTraceback (click to the left for traceback)\n...\nTraceback (most recent call last):'
        sage: sage.server.notebook.cell.format_exception(sage.server.notebook.cell.TRACEBACK + "notracebacks",80)
        'Traceback (most recent call last):notracebacks'
    """
    s = s0.lstrip()
    # Add a notracebacks option -- if it is in the string then tracebacks aren't shrunk.
    # This is currently used by the sage.server.support.help command.
    if TRACEBACK not in s or 'notracebacks' in s:
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
    r"""
    Returns the number of rows needed to display the string in txt if
    there are a maximum of ncols columns per row.

    EXAMPLES::

        sage: from sage.server.notebook.cell import number_of_rows
        sage: s = "asdfasdf\nasdfasdf\n"
        sage: number_of_rows(s, 8)
        2
        sage: number_of_rows(s, 5)
        4
        sage: number_of_rows(s, 4)
        4
    """
    rows = txt.splitlines()
    nrows = len(rows)
    for i in range(nrows):
        nrows += int((len(rows[i])-1)/ncols)
    return nrows
