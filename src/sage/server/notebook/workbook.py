"""
A Workbook.

A workbook is embedded in a webpage that is served by the SAGE server.
It is a linearly-ordered collections of numbered cells, where a
cell is a single input/output block.
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import os
import re
import string

import pexpect

from sage.ext.sage_object  import SageObject
from sage.interfaces.sage0 import Sage
from sage.misc.preparser   import preparse_file
from sage.misc.misc        import verbose

from cell import Cell

INTERRUPT_TRIES = 40

class Workbook:
    def __init__(self, name, notebook):
        name = ' '.join(name.split())
        self.__next_id = 0
        self.__name = name
        self.__notebook = notebook
        self.__dir = '%s/%s'%(notebook.workbook_directory(), '_'.join(name.split()))
        self.__comp_is_running = False
        if not os.path.exists(self.__dir):
            os.makedirs(self.__dir)
        self.__queue = []
        self.__cells = [ ]
        self.append_new_cell()
        #for i in range(20):
        #    self.__cells.append(self.new_cell())

    # The following setstate method is here
    # so that when this object is pickled and
    # unpickled, the self.__sage attribute
    # will not be set, so it will properly initialized.
    def __setstate__(self, state):
        self.__dict__ = state
        try:
            del self.__sage
            self.__queue = []
        except AttributeError:
            pass

    def cell_id_list(self):
        return [C.id() for C in self.__cells]

    def cell_list(self):
        return self.__cells

    def append_new_cell(self):
        """
        Create an append a new cell to the list of cells.
        """
        C = self._new_cell()
        self.__cells.append(C)
        return C

    def new_cell_before(self, id):
        """
        Insert a new cell into the cell list before the cell
        with the given integer id.  If the id is not the
        id of any cell, inserts a new cell at the end of the
        cell list.

        INPUT:
            id -- integer

        OUTPUT:
            an empty cell
        """
        cells = self.__cells
        for i in range(len(cells)):
            if cells[i].id() == id:
                C = self._new_cell()
                cells.insert(i, C)
                return C
        C = self._new_cell()
        cells.append(C)
        return C

    def delete_cell_with_id(self, id):
        """
        Remove the cell with given id and return the next cell after it.
        """
        cells = self.__cells
        for i in range(len(cells)):
            if cells[i].id() == id:
                del cells[i]
                if i < len(cells):
                    return cells[i].id()
        return cells[-1].id()

    def directory(self):
        return self.__dir

    def sage(self):
        try:
            return self.__sage
        except AttributeError:
            verbose("Initializing SAGE.")
            os.environ['PAGER'] = 'cat'
            self.__sage = Sage(logfile='%s/sage.log'%self.directory())
            S = self.__sage
            S.eval('import sage.server.support as _support_')
            S.eval('__SAGENB__globals = set(globals().keys())')
            object_directory = os.path.abspath(self.__notebook.object_directory())
            verbose(object_directory)
            S.eval('_support_.init("%s")'%object_directory)
            return S


    def _new_cell(self):
        D = self.__notebook.defaults()
        id = self.__next_id
        self.__next_id += 1
        return Cell(id, '', '', self)

    def __repr__(self):
        return str(self.__cells)

    def __len__(self):
        return len(self.__cells)

    def __getitem__(self, n):
        return self.__cells[n]

    def get_cell_with_id(self, id):
        for c in self.__cells:
            if c.id() == id:
                return c
        raise IndexError

    def queue(self):
        return self.__queue

    def enqueue(self, C):
        if not isinstance(C, Cell):
            raise TypeError
        if C.workbook() != self:
            raise ValueError, "C must be have self as workbook."
        if not (C in self.__queue):
            self.__queue.append(C)
        self.start_next_comp()

    def start_next_comp(self):
        if len(self.__queue) == 0:
            return
        if self.__comp_is_running:
            return

        self.__comp_is_running = True

        C = self.__queue[0]
        if C.interrupted():
            # don't actually compute
            return

        D = C.directory()
        S = self.sage()
        S._eval_line('os.chdir("%s")'%os.path.abspath(D))
        if C.time():
            S._eval_line('__t__=cputime(); __w__=walltime()')


        tmp = '%s/tmp.py'%self.directory()
        i = self.preparse_input(C.input_text(), C.completions())
        open(tmp,'w').write(i + '\nprint ""\n')
        S._send('execfile("%s")'%os.path.abspath(tmp))


    def check_comp(self):
        if len(self.__queue) == 0:
            return 'e', None
        S = self.sage()
        C = self.__queue[0]
        if C.interrupted():
            self.__comp_is_running = False
            del self.__queue[0]
            return 'd', C

        try:
            done, out = S._so_far()
        except RuntimeError:
            verbose("Computation was interrrupted or failed. Restarting.")
            self.__comp_is_running = False
            self.start_next_comp()
            return 'w', C

        out = self.postprocess_output(out, C.completions())
        if not done:
            # Still computing
            C.set_output_text(out)
            return 'w', C

        # Finished a computation.
        self.__comp_is_running = False
        if C.time():
            tm = S._eval_line('print "CPU time: %.2f s,  Wall time: %.2f s"%(cputime(__t__), walltime(__w__))')
            out = tm + '\n' + out
        out = self._process_output(out)
        C.set_output_text(out + C.files_html())
        del self.__queue[0]
        return 'd', C

    def _process_output(self, s):
        s = s.replace('<','&lt;')
        s = re.sub('\x08.','',s)
        return s

    def is_last_id_and_previous_is_nonempty(self, id):
        if self.__cells[-1].id() != id:
            return False
        if len(self.__cells) == 1:
            return False
        if len(self.__cells[-2].output_text(ncols=0)) == 0:
            return False
        return True

    def interrupt(self):
        """
        Interrupt all currently queued up calculations.

        OUTPUT:
            bool -- return True if no problems interrupting calculation
                    return False if the SAGE interpreter had to be restarted.
        """
        if len(self.__queue) == 0:
            # nothing to do
            return True

        # stop the current computation in the running SAGE
        S = self.sage()
        E = S._expect
        t = E.timeout
        E.timeout = 0.2
        success = False
        for i in range(INTERRUPT_TRIES):
            E.sendline('q')
            E.sendline(chr(3))
            try:
                E.expect(S._prompt)
                E.expect(S._prompt)
                success = True
                E.timeout = t
                break
            except (pexpect.TIMEOUT, pexpect.EOF), msg:
                verbose("Trying again to interrupt SAGE (try %s)..."%i)

        if not success:
            del self.__sage

        # empty the queue
        for C in self.__queue:
            C.interrupt()

        return success

    def postprocess_output(self, out, completions=False):
        i = out.find('\r\n')
        out = out[i+2:]
        out = out.rstrip()
        return out

    def _get_last_identifier(self, s):
        X = string.ascii_letters + string.digits + '._'
        i = len(s)-1
        while i >= 0 and s[i] in X:
            i -= 1
        return s[i+1:]

    def _input_contains_question_mark_query_not_in_quotes(self, x):
        """
        Return either False and None or True and the substring up to
        and including the question marks.
        """
        quote1 = False
        quote2 = False
        quote3 = False
        i = 0
        while i < len(x):
            if x[i] == "'":
                if not quote2 and not quote3:
                    quote1 = not quote1
            elif x[i] == '"':
                if not quote1 and not quote3:
                    quote2 = not quote2
            elif x[i] == '"""':
                if not quote1 and not quote2:
                    quote3 = not quote3
            elif x[i] == '?' and not quote1 and not quote2 and not quote3:
                if i + 1 < len(x) and (x[i+1] == '?' or x[i+1] == 'c'):
                    return True, x[:i+2]
                else:
                    return True, x[:i+1]
            i += 1
        return False, None


    def preparse_input(self, input, completions=False):
        input = input.strip()

        contains, new_input = self._input_contains_question_mark_query_not_in_quotes(input)
        if contains:
            input = new_input

        if len(input) > 1 and (input[-2:] == '??' or input[:2] == '??'):
            # source code request
            input = input.strip('?')
            input = self._get_last_identifier(input)
            input = 'print _support_.source_code("%s", globals())'%input

        elif len(input) > 0 and (input[-1] == '?' or input[0] == '?'):
            # docstring request
            input = input.strip('?')
            input = self._get_last_identifier(input)
            input = 'print _support_.docstring("%s", globals())'%input

        elif input[-2:] == '?c':
            input = self._get_last_identifier(input[:-2])
            input = 'print _support_.completions("%s", globals(), format=True)'%input

        elif completions:
            input = self._get_last_identifier(input)
            input = 'print _support_.completions("%s", globals(), format=True)'%input

        input = ignore_prompts_and_output(input)

        s = preparse_file(input, magic=False, do_time=True, ignore_prompts=True)
        s = [x for x in s.split('\n') if len(x.split()) > 0 and \
               x.lstrip()[0] != '#']   # remove all blank lines and comment lines
        if len(s) > 0:
            t = s[-1]
            if len(t) > 0 and not ':' in t and \
                   not t[0].isspace() and not t[:3] == '"""': # \
                   #and not t[:5] == 'print' and not '=' in t:
                #s[-1] = 'print %s'%t
                t = t.replace("'","\\'")
                s[-1] = "exec compile('%s', '', 'single')"%t
        s = '\n'.join(s) + '\n'
        return s

    def notebook(self):
        return self.__notebook

    def name(self):
        return self.__name

    def append(self, L):
        self.__cells.append(L)

    def variables(self):
        cmd = '",".join(["%s-%s"%(__x__,type(globals()[__x__])) for __x__ in globals().keys() if not __x__ in __SAGENB__globals and __x__[0] != "_"])'
        S = self.sage()
        v = S.eval(cmd)[1:-1]
        v = v.replace("<type '","").replace("<class '","").replace("'>","")
        w = v.split(',')
        w.sort()
        return w

    def variables_html(self):
        s = ''
        div = '<div class="variable_name" onClick="inspect_variable(\'%s\');">'
        for v in self.variables():
            try:
                name, typ = v.split('-')
            except ValueError:
                name = v; typ = ''
            if name:
                s += div%name + '<span class="varname">%s</span>&nbsp;<span class="vartype">(%s)</span></div>'%(name, typ)
        return s

    def html(self):
        n = len(self.__cells)
        s = ''
        D = self.__notebook.defaults()
        ncols = D['word_wrap_cols']
        for i in range(n):
            cell = self.__cells[i]
            s += cell.html(ncols) + '\n'
        return s


def ignore_prompts_and_output(s):
    """
    Given a string s that defines an input block of code,
    if any line begins in "sage:" (or ">>>"), strip out all lines
    that don't begin in either "sage:" (or ">>>") or"...", and
    remove all "sage:" (or ">>>") and "..." from the beginning
    of the remaining lines.
    """
    t = s.split('\n')
    do_strip = False
    for I in t:
        I2 = I.lstrip()
        if I2[:5] == 'sage:' or I2[:3] == '>>>':
            do_strip = True
            break
    if not do_strip:
        return s
    s = ''
    for I in t:
        I2 = I.lstrip()
        if I2[:5] == 'sage:':
            s += I2[5:].lstrip() + '\n'
        elif I2[:3] == '>>>':
            s += I2[3:].lstrip() + '\n'
        elif I2[:3] == '...':
            s += I2[3:] + '\n'
    return s



