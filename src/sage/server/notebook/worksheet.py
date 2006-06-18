"""
A Worksheet.

A worksheet is embedded in a webpage that is served by the SAGE server.
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

from sage.ext.sage_object  import load, save
from sage.interfaces.sage0 import Sage
from sage.misc.preparser   import preparse_file
from sage.misc.misc        import alarm, cancel_alarm, verbose, DOT_SAGE
import sage.server.support as support
from cell import Cell

INTERRUPT_TRIES = 20
INITIAL_NUM_CELLS = 1
import notebook as _notebook

SAGE_BEGIN='__SAGE_BEGIN__'
SAGE_END='__SAGE_END__'
SAGE_ERROR='error' + SAGE_END
SAGE_VARS='__SAGE_VARS__'

class Worksheet:
    def __init__(self, name, notebook, id):
        name = ' '.join(name.split())
        self.__id = id
        self.__next_id = (_notebook.MAX_WORKSHEETS) * id
        self.__name = name
        self.__notebook = notebook
        dir = list(name)
        for i in range(len(dir)):
            if not dir[i].isalnum() and dir[i] != '_':
                dir[i]='_'
        dir = ''.join(dir)
        self.__filename = dir
        self.__dir = '%s/%s'%(notebook.worksheet_directory(), dir)
        while os.path.exists(self.__dir):
            self.__dir += "_"
            self.__filename += '_'
        self.__comp_is_running = False
        if not os.path.exists(self.__dir):
            os.makedirs(self.__dir)
        self.__queue = []
        self.__cells = [ ]
        for i in range(INITIAL_NUM_CELLS):
            self.append_new_cell()

    def set_notebook(self, notebook, new_id=None):
        self.__notebook = notebook
        self.__dir = '%s/%s'%(notebook.worksheet_directory(), self.__filename)
        if not new_id is None:
            for C in self.__cells:
                i = C.relative_id()
                C.set_worksheet(self, new_id * _notebook.MAX_WORKSHEETS + i)
            self.__id = new_id
        else:
            for C in self.__cells:
                C.set_worksheet(self)

    def filename(self):
        return self.__filename

    def save(self, filename=None):
        if filename is None:
            save(self, os.path.abspath(self.__dir + '/' + self.__filename))
        else:
            save(self, filename)

    def __cmp__(self, other):
        try:
            return cmp((self.__id, self.__name), (other.__id, other.__name))
        except AttributeError:
            return -1

    def computing(self):
        return self.__comp_is_running

    def set_not_computing(self):
        self.__comp_is_running = False
        self.__queue = []

    def plain_text(self, prompts=False):
        """
        Return a plain-text version of the worksheet.

        prompts -- if True format for inclusion in docstrings.
        """
        s  = "#"*80 + '\n'
        s += "# Worksheet: %s"%self.name() + '\n'
        s += "#"*80+ '\n\n'
        for C in self.__cells:
            t = C.plain_text(prompts=prompts).strip()
            if t != '':
                s += '\n\n' + t
        return s

    # The following setstate method is here
    # so that when this object is pickled and
    # unpickled, the self.__sage attribute
    # will not be set, so it will properly initialized.
    def __setstate__(self, state):
        self.__dict__ = state
        try:
            del self.__sage
            del self.__variables
            self.__queue = []
        except AttributeError:
            pass

    def id(self):
        return self.__id

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

    def new_cell_after(self, id):
        """
        Insert a new cell into the cell list after the cell
        with the given integer id.
        """
        cells = self.__cells
        for i in range(len(cells)):
            if cells[i].id() == id:
                C = self._new_cell()
                cells.insert(i+1, C)
                return C
        C = self._new_cell()
        cells.append(C)
        return C

    def delete_cell_with_id(self, id):
        """
        Remove the cell with given id and return the cell before it.
        """
        cells = self.__cells
        for i in range(len(cells)):
            if cells[i].id() == id:
                del cells[i]
                if i > 0:
                    return cells[i-1].id()
                else:
                    return cells[0].id()
        return cells[0].id()

    def directory(self):
        return self.__dir

    def DIR(self):
        return self.__notebook.DIR()

    def sage(self):
        try:
            return self.__sage
        except AttributeError:
            verbose("Initializing SAGE.")
            os.environ['PAGER'] = 'cat'
            self.__sage = Sage()
            try:
                del self.__variables
            except AttributeError:
                pass
            S = self.__sage
            print "Starting SAGE server for worksheet %s..."%self.name()
            S.eval('__DIR__="%s/"; DIR=__DIR__'%self.DIR())
            S.eval('from sage.all_notebook import *')
            S.eval('import sage.server.support as _support_')
            S.eval('__SAGENB__globals = set(globals().keys())')
            object_directory = os.path.abspath(self.__notebook.object_directory())
            verbose(object_directory)
            S.eval('_support_.init("%s", globals())'%object_directory)
            print "(done)"

            A = self.attached_files()
            for F in A.iterkeys():
                A[F] = 0  # expire all

            return S


    def _new_cell(self, id=None):
        D = self.__notebook.defaults()
        if id is None:
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
        return self._new_cell(id)

    def queue(self):
        return self.__queue

    def enqueue(self, C):
        if not isinstance(C, Cell):
            raise TypeError
        if C.worksheet() != self:
            raise ValueError, "C must be have self as worksheet."
        if not (C in self.__queue):
            self.__queue.append(C)
        self.start_next_comp()

    def synchronize(self, s):
        try:
            i = (self.__synchro + 1)%65536
        except AttributeError:
            i = 0
        self.__synchro = i
        return 'print "%s%s"\n'%(SAGE_BEGIN,i) + s + '\nprint "%s%s"\n'%(SAGE_END,i)

    def synchro(self):
        try:
            return self.__synchro
        except AttributeError:
            return 0

    def start_next_comp(self):
        if len(self.__queue) == 0:
            return

        if self.__comp_is_running:
            return

        C = self.__queue[0]
        if C.interrupted():
            # don't actually compute
            return

        D = C.directory()
        V = self.known_variables()
        I = C.input_text().strip()
        if I in ['restart', 'quit', 'exit'] and not I in V:
            self.restart_sage()
            C.set_output_text('Restarted SAGE','')
            return
        elif I[:5] in ['time ', 'time\n', 'time\t'] and not 'time' in V:
            C.do_time()
            C.set_input_text(I[5:].lstrip())

        S = self.sage()


        tmp = '%s/tmp.py'%self.directory()
        input = 'os.chdir("%s")\n'%os.path.abspath(D)
        if C.time():
            input += '__SAGE_t__=cputime()\n__SAGE_w__=walltime()\n'
        if C.input_text()[-1:] == '?':
            C.set_introspect(C.input_text(), '')
        input += self.preparse_input(C)

        if C.time():
            input += 'print "CPU time: %.2f s,  Wall time: %.2f s"%(cputime(__SAGE_t__), walltime(__SAGE_w__))\n'

        if not C.introspect():
            input += 'print "%s'%SAGE_VARS + '=%s"%_support_.variables(True)'

        input = self.synchronize(input)

        open(tmp,'w').write(input)
        e = 'execfile("%s")\n'%os.path.abspath(tmp)
        # just in case, put an extra end...
        cmd = e + 'print "%s"+"%s"'%(SAGE_ERROR,self.synchro())
        self.__comp_is_running = True
        S._send(cmd)

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
            done, out, new = S._so_far(alternate_prompt=SAGE_END+str(self.synchro()))
        except RuntimeError:
            verbose("Computation was interrupted or failed. Restarting.")
            self.__comp_is_running = False
            self.start_next_comp()
            return 'w', C

        out = self.postprocess_output(out, C.introspect())
        if not done:
            # Still computing
            out = self._strip_synchro_from_start_of_output(out)
            C.set_output_text(out, C.files_html())
            return 'w', C

        # Finished a computation.
        self.__comp_is_running = False
        out = self._process_output(out)
        C.set_output_text(out, C.files_html())
        if C.introspect():
            before_prompt, after_prompt = C.introspect()
            if before_prompt[-1] != '?':
                # completions
                c = self.best_completion(out, C._word_being_completed)
                C.set_changed_input_text(before_prompt + c + after_prompt)

        del self.__queue[0]
        return 'd', C

    def best_completion(self, s, word):
        if 'no completions of' in s:
            return ''
        completions = s.split()
        if len(completions) == 0:
            return ''
        n = len(word)
        i = n
        m = min([len(x) for x in completions])
        while i <= m:
            word = completions[0][:i]
            for w in completions[1:]:
                if w[:i] != word:
                    return w[n:i-1]
            i += 1
        return completions[0][n:m]

    def _strip_synchro_from_start_of_output(self, s):
        z = SAGE_BEGIN+str(self.synchro())
        i = s.find(z)
        if i == -1:
            return s
        return s[i+len(z):]

    def _process_output(self, s):
        s = re.sub('\x08.','',s)
        s = self._strip_synchro_from_start_of_output(s)
        if SAGE_ERROR in s:
            i = s.rfind('>>>')
            if i >= 0:
                return s[:i-1]
        else:
            i = s.rfind(SAGE_VARS)
            if i != -1:
                t = s[i+len(SAGE_VARS)+1:]
                t = t.replace("<type '","").replace("<class '","").replace("'>","")
                self.__variables = eval(t)
                s = s[:i-1]
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
        E.timeout = 0.1
        success = False
        alarm(INTERRUPT_TRIES * E.timeout)
        try:
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
        except:
            print "Interrupted (escape via alarm)!"
            success = False
        else:
            # Turn off the alarm.
            cancel_alarm()

        if not success:
            del self.__sage

        # empty the queue
        for C in self.__queue:
            C.interrupt()

        return success

    def restart_sage(self):
        """
        Restart SAGE kernel.
        """
        # stop the current computation in the running SAGE
        try:
            S = self.__sage
        except AttributeError:
            # no sage running anyways!
            pass
        else:
            del self.__sage

        # empty the queue
        for C in self.__queue:
            C.interrupt()
        self.__comp_is_running = False

    def postprocess_output(self, out, introspect=False):
        i = out.find('\r\n')
        out = out[i+2:]
        out = out.rstrip()
        return out

    def _get_last_identifier(self, s):
        return support.get_rightmost_identifier(s)

    def preparse(self, s):
        return preparse_file(s, magic=False, do_time=False, ignore_prompts=False)

    def load_any_changed_attached_files(self, s):
        """
        Modify s by prepending any necessary load commands
        corresponding to attached files that have changed.
        """
        A = self.attached_files()
        for F, tm in A.iteritems():
            try:
                new_tm = os.path.getmtime(F)
            except OSError:
                pass
            else:
                if new_tm > tm:
                    A[F] = new_tm
                    s = 'load %s\n'%F + s
        return s

    def attached_files(self):
        try:
            A = self.__attached
        except AttributeError:
            A = {}

            init_sage = DOT_SAGE + 'init.sage'
            if os.path.exists(init_sage):
                A[init_sage] = 0

            self.__attached = A
        return A

    def attach(self, filename):
        A = self.attached_files()
        try:
            A[filename] = os.path.getmtime(filename)
        except OSError:
            print "WARNING: File %s vanished"%filename
            pass

    def detach(self, filename):
        A = self.attached_files()
        try:
            A.pop(filename)
        except KeyError:
            pass

    def _normalized_filenames(self, L):
        a = []
        OBJECTS = os.path.abspath(self.notebook().object_directory())
        for filename in L.split():
            filename = filename.strip('"').strip("'")
            if os.path.exists(OBJECTS + '/' + filename):
                filename = OBJECTS + '/' + filename
            elif os.path.exists(OBJECTS + '/' + filename + '.sobj'):
                filename = OBJECTS + '/' + filename + '.sobj'
            else:
                if len(filename) > 0 and filename[0] != '/':
                    filename = '%s/%s'%(self.DIR(), filename)
                if filename[-3:] != '.py' and filename[-5:] != '.sage' and \
                   filename[-5:] != '.sobj' and not os.path.exists(filename):
                    if os.path.exists(filename + '.sage'):
                        filename = filename + '.sage'
                    elif os.path.exists(filename + '.py'):
                        filename = filename + '.py'
                    elif os.path.exists(filename + '.sobj'):
                        filename = filename + '.sobj'
            a.append(filename)
        return a

    def _load_file(self, filename, files_seen_so_far, this_file):
        if filename[-5:] == '.sobj':
            i = filename.rfind('/')
            if i != -1:
                name = filename[i+1:-5]
            else:
                name = filename[:-5]
            return '%s = load("%s");'%(name, filename)

        if filename in files_seen_so_far:
            return "print 'WARNING: Not loading %s -- would create recursive load'"%filename
        try:
            F = open(filename).read()
        except IOError:
            t = "print 'Error loading %s -- file not found'"%filename
        else:
            if filename[-3:] == '.py':
                t = F
            elif filename[-5:] == '.sage':
                t = self.preparse(F)
        t = self.do_sage_extensions_preparsing(t,
                          files_seen_so_far + [this_file], filename)
        return t

    def _save_objects(self, s):
        s = s.replace(',',' ').replace('(',' ').replace(')',' ')
        v = s.split()
        return ';'.join(['save(%s,"%s")'%(x,x) for x in v])


    def do_sage_extensions_preparsing(self, s, files_seen_so_far=[], this_file=''):
        u = []
        for t in s.split('\n'):
            if t[:5] == 'load ':
                z = ''
                for filename in self._normalized_filenames(t[5:]):
                    z += self._load_file(filename, files_seen_so_far, this_file) + '\n'
                t = z

            elif t[:7] == 'attach ':
                z = ''
                for filename in self._normalized_filenames(t[7:]):
                    if not os.path.exists(filename):
                        z += "print 'Error attaching %s -- file not found'\n"%filename
                    else:
                        self.attach(filename)
                        z += self._load_file(filename, files_seen_so_far, this_file) + '\n'
                t = z

            elif t[:7]  == 'detach ':
                for filename in self._normalized_filenames(t[7:]):
                    self.detach(filename)
                t = ''

            elif t[:12] in ['save_session', 'load_session']:
                F = t[12:].strip().strip('(').strip(')').strip("'").strip('"').split(',')[0]
                if len(F) == 0:
                    filename = self.__filename
                else:
                    filename = F
                if t[:4] == 'save':
                    t = '_support_.save_session("%s")'%filename
                else:
                    t = 'load_session(locals(), "%s")'%filename

            elif t[:5] == 'save ':
                t = self._save_objects(t[5:])

            u.append(t)

        return '\n'.join(u)


    def check_for_system_switching(self, s):
        # s = '%gap\n' + s
        if len(s) == 0:
            return False, s
        if s[0] == '%':
            i = s.find('\n')
            if i == -1:
                # nothing to evaluate
                return True, ''
            j = s.find(' ')
            if j == -1:
                j = i
            else:
                j = min(i,j)
            sys = s[1:j]
            s = s[i+1:]
            if sys in ['latex', 'latex_debug', 'slide', 'slide_debug']:
                t = 'print %s.eval(r"""%s""", vars=globals())'%(sys,s)
            else:
                t = 'print %s.eval(r"""%s""")'%(sys,s)
            return True, t
        return False, s

    def preparse_input(self, C):
        input = C.input_text().strip()
        introspect = C.introspect()
        if introspect:
            before_prompt, after_prompt = introspect
            i = 0
            while i < len(after_prompt):
                if after_prompt[i] == '?':
                    if i < len(after_prompt)-1 and after_prompt[i+1] == '?':
                        i += 1
                    before_prompt += after_prompt[:i+1]
                    after_prompt = after_prompt[i+1:]
                    C.set_introspect(before_prompt, after_prompt)
                    break
                elif after_prompt[i] in ['"', "'", ' ', '\t', '\n']:
                    break
                i += 1
            if before_prompt[-2:] == '??':
                input = self._get_last_identifier(before_prompt[:-2])
                input = 'print _support_.source_code("%s", globals())'%input
            elif before_prompt[-1:] == '?':
                input = self._get_last_identifier(before_prompt[:-1])
                input = 'print _support_.docstring("%s", globals())'%input
            else:
                input = self._get_last_identifier(before_prompt)
                C._word_being_completed = input
                input = 'print _support_.completions("%s", globals(), format=True)\n'%input

        else:
            switched, input = self.check_for_system_switching(input)

            if not switched:
                input = ignore_prompts_and_output(input)
                input = self.preparse(input)
                input = self.load_any_changed_attached_files(input)
                input = self.do_sage_extensions_preparsing(input)

                input = [x for x in input.split('\n') if len(x.split()) > 0 and \
                         x.lstrip()[0] != '#']   # remove all blank lines and comment lines

                if len(input) > 0:
                    t = input[-1]
                    if not switched and len(t) > 0 and not ':' in t and \
                       not t[0].isspace() and not t[:3] == '"""' and not t[:3] == "'''":
                        t = t.replace("'", "\\u0027")
                        input[-1] = "exec compile(ur'%s', '', 'single')"%t
                input = '\n'.join(input)

            input += '\n'

        return input

    def notebook(self):
        return self.__notebook

    def name(self):
        return self.__name

    def append(self, L):
        self.__cells.append(L)

    def known_variables(self):
        try:
            return self.__variables
        except AttributeError:
            return []

    def variables(self, with_types=True):
        try:
            self.__sage
        except AttributeError:
            try:
                del self.__variables
            except AttributeError:
                pass
            return []
        try:
            v = self.__variables
        except AttributeError:
            return []
        if with_types:
            return v
        else:
            return [x.split('-')[0] for x in v]

    def variables_html(self):
        s = ''
        div = '<div class="variable_name">'
        for v in self.variables():
            try:
                name, typ = v.split('-')
            except ValueError:
                name = v; typ = ''
            if name:
                s += div + '<span class="varname">%s</span>&nbsp;<span class="vartype">(%s)</span></div>'%(name, typ)
        return s

    def attached_html(self):
        s = ''
        div = '<div class="attached_filename" onClick="inspect_attached_file(\'%s\')">'
        A = self.attached_files()
        D = self.DIR()
        for F, tm in A.iteritems():
            # uncomment this to remove some absolute path info...
            # if F[:len(D)] == D: F = F[len(D)+1:]
            s += div%F + '%s</div>'%F
        return s

    def html(self):
        n = len(self.__cells)
        s = ''

        s += '<div class="worksheet_title">Worksheet: %s</div>\n'%self.name()
        D = self.__notebook.defaults()
        ncols = D['word_wrap_cols']
        s += '<div class="worksheet_cell_list" id="worksheet_cell_list">\n'
        for i in range(n):
            cell = self.__cells[i]
            s += cell.html(ncols) + '\n'
        s += '\n</div>\n'
        s += '<div class="worksheet_bottom_padding"></div>\n'
        s += '<script language=javascript>cell_id_list=%s</script>\n'%self.cell_id_list()
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



