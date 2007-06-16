"""
A Worksheet.

A worksheet is embedded in a webpage that is served by the SAGE server.
It is a linearly-ordered collections of numbered cells, where a
cell is a single input/output block.
"""

from __future__ import with_statement

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import os
import re
import string
import traceback
import time
import crypt
import pexpect

from sage.structure.sage_object  import load, save
from sage.interfaces.sage0 import Sage
from sage.misc.preparser   import preparse_file
from sage.misc.misc        import alarm, cancel_alarm, verbose, DOT_SAGE
import sage.server.support as support

from cell import Cell, TextCell

INTERRUPT_TRIES = 3

INITIAL_NUM_CELLS = 1
HISTORY_MAX_OUTPUT = 92*5
HISTORY_NCOLS = 90


import notebook as _notebook

#If you make any changes to this, be sure to change the
# error line below that looks like this:
#         cmd += 'print "\\x01r\\x01e%s"'%self.synchro()
SC='\x01'
SAGE_BEGIN=SC+'b'
SAGE_END=SC+'e'
SAGE_ERROR=SC+'r'


# The default is for there to be one sage session for
# each worksheet.  If this is False, then there is just
# one global SAGE session, like with Mathematica.
# This variable gets sets when the notebook function
# in notebook.py is called.
multisession = True
def initialized_sage():
    S = Sage(maxread = 1)
    S._start(block_during_init=False)
    E = S.expect()
    E.sendline('\n')
    E.expect('>>>')
    cmd = 'from sage.all_notebook import *;'
    cmd += 'import sage.server.support as _support_; '
    E.sendline(cmd)
    return S



_a_sage = None
def init_sage_prestart():
    global _a_sage
    _a_sage = initialized_sage()

def one_prestarted_sage():
    global _a_sage
    X = _a_sage
    if multisession:
        init_sage_prestart()
    return X

class Worksheet:
    def __init__(self, name, notebook, id, system=None, passcode = ''):
        name = ' '.join(name.split())
        self.__id = id
        self.__system = system
        self.__next_id = (_notebook.MAX_WORKSHEETS) * id
        self.set_name(name)
        self.__notebook = notebook
        self.__passcode = crypt.crypt(passcode, self.salt())
        self.__passcrypt= True
        self.__dir = '%s/%s'%(notebook.worksheet_directory(), dir)
        self.clear()

    def clear(self):
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
            nid = self.__next_id  - _notebook.MAX_WORKSHEETS*self.__id
            self.__id = new_id
            self.__next_id = _notebook.MAX_WORKSHEETS*new_id + nid
        else:
            for C in self.__cells:
                C.set_worksheet(self)

    def salt(self):
        try:
            return self.__salt
        except AttributeError:
            self.__salt = "%f"%time.time()
            return self.__salt

    def passcode(self):
        try:
            c = self.__passcrypt
        except AttributeError:
            self.__passcrypt = False
        try:
            if not self.__passcrypt:
                self.__passcode = crypt.crypt(self.__passcode, self.salt())
                self.__passcrypt = True
            return self.__passcode
        except AttributeError:
            self.__passcode = crypt.crypt('', self.salt())
            self.__passcrypt = True
            return self.__passcode

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
        try:
            return self.__comp_is_running
        except AttributeError:
            return False

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
            t = C.plain_text(prompts=prompts).strip('\n')
            if t != '':
                s += '\n' + t
        return s

    def edit_text(self, prompts=False):
        """
        Returns a plain-text version of the worksheet with \{\{\{\}\}\} wiki-formatting,
        suitable for hand editing.
        """
        s = ''
        for C in self.__cells:
            t = C.edit_text(prompts=prompts).strip()
            if t != '':
                s += '\n\n' + t
        return s

    def edit_save(self, text):
        text.replace('\r\n','\n')
        # This is where we would save the last version in a history.
        cells = []
        while True:
            plain_text = extract_text_before_first_compute_cell(text).strip()
            if len(plain_text) > 0:
                T = self._new_text_cell(plain_text)
                cells.append(T)
            try:
                input, output, graphics, i = extract_first_compute_cell(text)
            except EOFError:
                cells.append(self._new_cell())  # make sure last cell is a compute cell
                break
            text = text[i:]
            C = self._new_cell()
            C.set_input_text(input)
            C.set_output_text(output, '')
            cells.append(C)
        if len(cells) == 0:   # there must be at least one cell.
            cells = [self._new_cell()]
        self.__cells = cells

##         lines = text.split('\n')
##         input = ""
##         output = ""
##         in_cell = False
##         in_output = False
##         old_first = self.__cells[0].id()
##         for line in lines:
##             if not in_cell:
##                 if line[:3] == '{{{':
##                     in_cell = True
##             elif line != '}}}':
##                 if not in_output:
##                     if line == '///':
##                         in_output = True
##                     else:
##                         input += line+'\n'
##                 else:
##                     output += line+'\n'
##             else:
##                 C = self.new_cell_before(old_first)
##                 C.set_input_text(input)
##                 C.set_output_text(output,output)
##                 C.set_cell_output_type()
##                 input = ""
##                 output = ""
##                 in_cell = False
##                 in_output = False

    def input_text(self):
        """
        Return text version of the input to the worksheet.
        """
        return '\n\n---\n\n'.join([C.input_text() for C in self.__cells])

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

    def compute_cell_id_list(self):
        return [C.id() for C in self.__cells if isinstance(C, Cell)]

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
                    break
        return cells[0].id()

    def directory(self):
        if not os.path.exists(self.__dir):
            # prevent "rm -rf" accidents.
            os.makedirs(self.__dir)
        return self.__dir

    def DIR(self):
        return self.__notebook.DIR()

    def quit(self):
        self.restart_sage()

    def next_block_id(self):
        try:
            i = self.__next_block_id
        except AttributeError:
            i = 0
        i += 1
        self.__next_block_id = i
        return i

    def compute_process_has_been_started(self):
        """
        Return True precisely if the compute process has been started,
        irregardless of whether or not it is currently churning away
        on a computation.
        """
        try:
            S = self.__sage
            if S._expect is None:
                return False
        except AttributeError:
            return False
        return True

    def initialize_sage(self):
        print "Starting SAGE server for worksheet %s..."%self.name()
        self.delete_cell_input_files()
        object_directory = os.path.abspath(self.__notebook.object_directory())
        S = self.__sage
        try:
            cmd = '__DIR__="%s/"; DIR=__DIR__;'%self.DIR()
            cmd += '_support_.init("%s", globals()); '%object_directory
            S._send(cmd)   # non blocking
        except Exception, msg:
            print "ERROR initializing compute process:\n"
            print msg
            del self.__sage
            raise RuntimeError
        print "Done starting"
        A = self.attached_files()
        for F in A.iterkeys():
            A[F] = 0  # expire all
        return S

    def sage(self):
        try:
            S = self.__sage
            if not S._expect is None:
                return S
        except AttributeError:
            pass
        self.__sage = one_prestarted_sage()
        verbose("Initializing SAGE.")
        os.environ['PAGER'] = 'cat'
        try:
            del self.__variables
        except AttributeError:
            pass
        self.__next_block_id = 0
        self.initialize_sage()
        return self.__sage

    def _enqueue_auto_cells(self):
        for c in self.__cells:
            if c.is_auto_cell():
                self.enqueue(c)

    def _new_text_cell(self, plain_text, id=None):
        if id is None:
            id = self.__next_id
            self.__next_id += 1
        return TextCell(id, plain_text, self)

    def _new_cell(self, id=None):
        if id is None:
            id = self.__next_id
            self.__next_id += 1
        return Cell(id, '', '', self)

    def __repr__(self):
        return str(self.__cells)

    def __len__(self):
        return len(self.__cells)

    def __getitem__(self, n):
        try:
            return self.__cells[n]
        except IndexError:
            if n >= 0:  # this should never happen -- but for robustness we cover this case.
                for k in range(len(self.__cells),n+1):
                    self.__cells.append(self._new_cell())
                return self.__cells[n]
            raise IndexError


    def get_cell_with_id(self, id):
        for c in self.__cells:
            if c.id() == id:
                return c
        return self._new_cell(id)

    def queue(self):
        return list(self.__queue)

    def queue_id_list(self):
        return [c.id() for c in self.__queue]

    def _enqueue_auto(self):
        for c in self.__cells:
            if c.is_auto_cell():
                self.__queue.append(c)


    def enqueue(self, C):
        if not isinstance(C, Cell):
            raise TypeError
        if C.worksheet() != self:
            raise ValueError, "C must be have self as worksheet."
        # If the SAGE server hasn't started and the queue is empty,
        # first enqueue the auto cells:
        if len(self.__queue) == 0:
            try:
                self.__sage
            except AttributeError:
                self._enqueue_auto()

        # Now enqueue the requested cell.
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

    def delete_cell_input_files(self):
        """
        Delete all the files code_%s.py and code_%s.spyx that are created
        when evaluating cells.  We do this when we first start the notebook
        to get rid of clutter.
        """
        D = self.directory() + '/code/'
        if os.path.exists(D):
            for X in os.listdir(D):
                os.unlink('%s/%s'%(D,X))
        else:
            os.makedirs(D)

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
        if not C.introspect():
            I = C.input_text().strip()
            if I in ['restart', 'quit', 'exit'] and not I in V:
                self.restart_sage()
                S = self.system()
                if S is None: S = 'SAGE'
                C.set_output_text('Exited %s process'%S,'')
                return
            if I[:5] == '%time':
                C.do_time()
                I = I[5:].lstrip()
            elif I[:5] in ['time ', 'time\n', 'time\t'] and not 'time' in V:
                C.do_time()
                I = I[5:].lstrip()
        else:
            I = C.introspect()[0]

        S = self.sage()

        id = self.next_block_id()

        # prevent directory disappear problems
        if not os.path.exists('%s/code'%self.directory()):
            os.makedirs('%s/code'%self.directory())
        tmp = '%s/code/%s.py'%(self.directory(), id)
        input = 'os.chdir("%s")\n'%os.path.abspath(D)

        if C.time():
            input += '__SAGE_t__=cputime()\n__SAGE_w__=walltime()\n'
        if I[-1:] == '?':
            C.set_introspect(I, '')
        I = I.replace('\\\n','')
        C._before_preparse = input + I
        input += self.preparse_input(I, C)

        try:
            compile(input, '', 'exec')
        except SyntaxError, msg:
            t = traceback.format_exc()
            s = 'File "<string>",'
            i = t.find(s)
            if i != -1:
                t = t[i+len(s):]
            i = t.find('\n')
            try:
                n = int(t[t[:i].rfind(' '):i])  # line number of the exception
                try:
                    t = 'Syntax Error:\n    %s'%C._before_preparse.split('\n')[n-1]
                except IndexError:
                    pass
                if False:
                    if i != -1:
                        t = t[i:]
                    v = [w for w in t.split('\n') if w]
                    t = '\n'.join(['Syntax Error:'] + v[0:-1])
                C.set_output_text(t, '')
                del self.__queue[0]
                return
            except ValueError:
                pass

        if C.time() and not C.introspect():
            input += 'print "CPU time: %.2f s,  Wall time: %.2f s"%(cputime(__SAGE_t__), walltime(__SAGE_w__))\n'


        input = self.synchronize(input)
        # Unfortunately, this has to go here at the beginning of the file until Python 2.6,
        # in order to support use of the with statement in the notebook.  Very annoying.
        input = 'from __future__ import with_statement\n' + input


        open(tmp,'w').write(input)

        cmd = 'execfile("%s")\n'%os.path.abspath(tmp)
        # Signal an end (which would only be seen if there is an error.)
        cmd += 'print "\\x01r\\x01e%s"'%self.synchro()
        self.__comp_is_running = True
        try:
            S._send(cmd)
        except OSError, msg:
            self.restart_sage()
            C.set_output_text('The SAGE compute process quit (possibly SAGE crashed?).\nPlease retry your calculation.','')



    def check_cell(self, id):
        """
        Check the status on computation of the cell with given id.

        INPUT:
            id -- an integer

        OUTPUT:
            status -- a string, either 'd' (done) or 'w' (working)
            cell -- the cell with given id
        """
        cell = self.get_cell_with_id(id)

        if cell in self.__queue:
            status = 'w'
        else:
            status = 'd'
        return status, cell

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
            done, out, new = S._so_far(wait=0.2, alternate_prompt=SAGE_END+str(self.synchro()))
        except RuntimeError, msg:
            verbose("Computation was interrupted or failed. Restarting.\n%s"%msg)
            self.__comp_is_running = False
            self.start_next_comp()
            return 'w', C

        out = self.postprocess_output(out, C)
        if not done:
            # Still computing
            out = self._process_output(out)
            if not C.introspect():
                C.set_output_text(out, '')
            return 'w', C

        # Finished a computation.
        self.__comp_is_running = False
        del self.__queue[0]
        out = self._process_output(out)
        if C.introspect():
            before_prompt, after_prompt = C.introspect()
            if len(before_prompt) == 0:
                return
            if before_prompt[-1] != '?':
                # completions
                c = self.best_completion(out, C._word_being_completed)
                C.set_changed_input_text(before_prompt + c + after_prompt)
                out = self.completions_html(C.id(), out)
                C.set_introspect_html(out, completing=True)
            else:
                C.set_introspect_html(out, completing=False)
        else:
            C.set_output_text(out, C.files_html(out), sage=self.sage())
            C.set_introspect_html('')
            history = "Worksheet '%s' (%s)\n"%(self.name(), time.strftime("%Y-%m-%d at %H:%M",time.localtime(time.time())))
            history += C.edit_text(ncols=HISTORY_NCOLS, prompts=False,
                                    max_out=HISTORY_MAX_OUTPUT)
            self.notebook().add_to_history(history)

        return 'd', C

    def best_completion(self, s, word):
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

    def completions_html(self, id, s, cols=3):
        if 'no completions of' in s:
            return ''

        completions = s.split()

        n = len(completions)
        l = n/cols + n%cols

        if n == 1:
            return '' # don't show a window, just replace it

        rows = []
        for r in range(0,l):
            row = []
            for c in range(cols):
                try:
                    cell = completions[r + l*c]
                    row.append(cell)
                except:
                    pass
            rows.append(row)
        return self.__notebook.format_completions_as_html(id, rows)

    def auth(self, passcode):
        return self.passcode() == crypt.crypt(passcode, self.salt())

    def _strip_synchro_from_start_of_output(self, s):
        z = SAGE_BEGIN+str(self.synchro())
        i = s.find(z)
        if i == -1:
            # did not find any synchronization info in the output stream
            j = s.find('Traceback')
            if j != -1:
                # Probably there was an error; better not hide it.
                return s[j:]
            else:
                # Maybe we just read too early -- supress displaying anything yet.
                return ''
        else:
            return s[i+len(z):]

    def _process_output(self, s):
        s = re.sub('\x08.','',s)
        s = self._strip_synchro_from_start_of_output(s)
        if SAGE_ERROR in s:
            i = s.rfind('>>>')
            if i >= 0:
                return s[:i-1]
        # Remove any control codes that might have not got stripped out.
        return s.replace(SAGE_BEGIN,'').replace(SAGE_END,'').replace(SC,'')

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

        success = False
        # stop the current computation in the running SAGE
        try:
            S = self.__sage
        except AttributeError:
            pass
        else:
            success = S.interrupt(INTERRUPT_TRIES, timeout=0.3, quit_on_fail=False)

        if success:
            self.clear_queue()

        return success

    def clear_queue(self):
        # empty the queue
        for C in self.__queue:
            C.interrupt()
        self.__queue = []
        self.__comp_is_running = False

    def restart_sage(self):
        """
        Restart SAGE kernel.
        """
        print "restarting"

        try:
            S = self.__sage
        except AttributeError:
            # no sage running anyways!
            return

        try:
            pid = S._expect.pid
            os.killpg(pid, 9)
            os.kill(pid, 9)
            S._expect = None
            del self.__sage
        except AttributeError, msg:
            print "WARNING: %s"%msg
        except Exception, msg:
            print msg
            print "WARNING: Error deleting SAGE object!"

        try:
            del self.__variables
        except AttributeError:
            pass

        # We do this to avoid getting a stale SAGE that uses old code.
        self.clear_queue()
        self.__sage = initialized_sage()
        self.initialize_sage()
        self._enqueue_auto_cells()
        self.start_next_comp()

    def postprocess_output(self, out, C):
        i = out.find('\r\n')
        out = out[i+2:]
        #out = out.rstrip()
        if C.introspect():
            return out

        # this isn't needed anymore !
        # the python error message for list indices is not good enough.
        # out = out.replace('indices must be integers', 'indices must be of type Python int.\n(Hint: Use int(n) to make n into a Python int.)')

        out = out.replace("NameError: name 'os' is not defined", "NameError: name 'os' is not defined\nTHERE WAS AN ERROR LOADING THE SAGE LIBRARIES.  Try starting SAGE from the command line to see what the error is.")

        try:
            tb = 'Traceback (most recent call last):'
            i = out.find(tb)
            if i != -1:
                t = '.py", line'
                j = out.find(t)
                z = out[j+5:].find(',')
                n = int(out[j+len(t):z + j+5])
                k = out[j:].find('\n')
                if k != -1:
                    k += j
                    l = out[k+1:].find('\n')
                    if l != -1:
                        l += k+1
                        I = C._before_preparse.split('\n')
                        out = out[:i + len(tb)+1] + '    ' + I[n-2] + out[l:]
        except (ValueError, IndexError), msg:
            pass
        return out

    def _get_last_identifier(self, s):
        t = support.get_rightmost_identifier(s)
        S = self.system()
        # If we are tab completing in a worksheet for another
        # system, e.g., mathematica, this is like typing mathematica.[tab].
        # However, for SAGE and Python, we want regular tab completion.
        if S and not (S in ['sage', 'python']):
            t = S + '.' + t
        return t

    def preparse(self, s):
        s = preparse_file(s, magic=False, do_time=True,
                          ignore_prompts=False)
        return s

    def load_any_changed_attached_files(self, s):
        """
        Modify s by prepending any necessary load commands
        corresponding to attached files that have changed.
        """
        A = self.attached_files()
        init_sage = DOT_SAGE + 'init.sage'
        if not init_sage in A.keys() and os.path.exists(init_sage):
            A[init_sage] = 0

        for F, tm in A.iteritems():
            try:
                new_tm = os.path.getmtime(F)
            except OSError:
                del A[F]
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
        i = L.find('#')
        if i != -1:
            L = L[:i]
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
            t = "print 'WARNING: Not loading %s -- would create recursive load'"%filename
        try:
            F = open(filename).read()
        except IOError:
            return "print 'Error loading %s -- file not found'"%filename
        else:
            if filename[-3:] == '.py':
                t = F
            elif filename[-5:] == '.sage':
                t = self.preparse(F)
            else:
                t = "print 'Loading of file \"%s\" has type not implemented.'"%filename

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

    def system(self):
        try:
            return self.__system
        except AttributeError:
            self.__system = None
            return None

    def set_system(self, system=None):
        self.__system = system

    def _eval_cmd(self, system, cmd):
        cmd = cmd.replace("'", "\\u0027")
        return "print _support_.syseval(%s, ur'''%s''')"%(system, cmd)

    def sagex_import(self, cmd, C):
        # Choice: Can use either C.relative_id() or self.next_block_id().
        # C.relative_id() has the advantage that block evals are cached, i.e.,
        # no need to recompile.  On the other hand tracebacks don't work if
        # you change a cell and create a new function in it.  Caching is
        # also annoying since the linked .c file disappears.
        id = self.next_block_id()
        # id = C.relative_id()
        spyx = os.path.abspath('%s/code/sage%s.spyx'%(self.directory(), id))
        if not (os.path.exists(spyx) and open(spyx).read() == cmd):
            open(spyx,'w').write(cmd)
        s  = '_support_.sagex_import_all("%s", globals())'%spyx
        return s

    def check_for_system_switching(self, s, C):
        r"""
        Check for input cells that start with \code{\%foo},
        where foo is an object with an eval method.
        """
        z = s
        s = s.lstrip()
        S = self.system()
        if not (S is None):
            if s[:5] != '%sage':
                return True, self._eval_cmd(self.__system, s)
            else:
                s = s[5:].lstrip()
                z = s

        if len(s) == 0 or s[0] != '%':
            return False, z
        if s[:5] == '%hide':
            t = s[5:].lstrip()
            if len(t) == 0 or t[0] != '%':
                return False, t
            s = t
        if s[:12] == '%save_server':
            self.notebook().save()
            t = s[12:].lstrip()
            if len(t) == 0 or t[0] != '%':
                return False, t
            s = t
        if s[:6] == "%pyrex" or s[:6] == "%sagex":  # a block of Sagex code.
            return True, self.sagex_import(s[6:].lstrip(), C)

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
        cmd = self._eval_cmd(sys, s)
        if sys == 'html':
            C.set_is_html(True)
        return True, cmd

    def preparse_input(self, input, C):
        C.set_is_html(False)
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
                input = 'print "\\n".join(_support_.completions("%s", globals()))'%(input)

        else:
            switched, input = self.check_for_system_switching(input, C)

            if not switched:
                input = ignore_prompts_and_output(input).rstrip()
                input = self.preparse(input)
                input = self.load_any_changed_attached_files(input)
                input = self.do_sage_extensions_preparsing(input)
                input = input.split('\n')

                # The following is all so the last line (or single lines)
                # will implicitly print as they should, unless they are
                # an assignment.   "display hook"  It's very complicated,
                # but it has to be...
                i = len(input)-1
                if i >= 0:
                    while len(input[i]) > 0 and input[i][0] in ' \t':
                        i -= 1
                    t = '\n'.join(input[i:])
                    if t[:4] != 'def ':
                        try:
                            compile(t+'\n', '', 'single')
                            t = t.replace("'", "\\u0027").replace('\n','\\u000a')
                            # IMPORTANT: If you change this line, also change
                            # the function format_exception in cell.py
                            input[i] = "exec compile(ur'%s' + '\\n', '', 'single')"%t
                            input = input[:i+1]
                        except SyntaxError, msg:
                            pass
                input = '\n'.join(input)

            input += '\n'

        #print input
        return input

    def notebook(self):
        return self.__notebook

    def name(self):
        return self.__name

    def set_name(self, name):
        self.__name = name
        dir = list(name)
        for i in range(len(dir)):
            if not dir[i].isalnum() and dir[i] != '_':
                dir[i] = '_'
        dir = ''.join(dir)
        self.__filename = dir

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

    def html(self, include_title=True, do_print=False,
             authorized=False, confirm_before_leave=False):
        n = len(self.__cells)
        s = ''
        if include_title:
            if self.computing():
                interrupt_class = "interrupt"
            else:
                interrupt_class = "interrupt_grey"
            S = self.system()
            if not (S is None):
                system = ' (%s mode)'%S
            else:
                system =''
            if not authorized:
                lock_text = '&nbsp;&nbsp;<span id="worksheet_lock" class="locked" onClick="unlock_worksheet()">[locked]</span>'
            else:
                lock_text = ''

            vbar = '<span class="vbar"></span>'

            name = self.filename()
            menu  = '  <span class="worksheet_control_commands">'
            menu += '    <a class="%s" onClick="interrupt()" id="interrupt">Interrupt</a>'%interrupt_class + vbar
            menu += '    <a class="restart_sage" onClick="restart_sage()" id="restart_sage">Restart</a>' +vbar
            menu += '    <a class="plain_text" href="edit">Edit</a>' + vbar
            menu += '    <a class="doctest_text" onClick="doctest_window(\'%s\')">Text</a>'%name + vbar
            menu += '    <a class="doctest_text" onClick="print_window(\'%s\')">Print</a>'%name + vbar
            menu += '    <a class="evaluate" onClick="evaluate_all()">Eval All</a>' + vbar
            menu += '    <a class="hide" onClick="hide_all()">Hide</a>/<a class="hide" onClick="show_all()">Show</a>' + vbar
            menu += '    <a class="slide_mode" onClick="slide_mode()">Focus</a>' + vbar
            menu += '    <a class="download_sws" href="download/%s.sws">Download</a>'%self.filename() + vbar
            menu += '    <a class="delete" href="delete">Delete</a>'
            menu += '  </span>'

            s += '<div class="worksheet_title">'
            s += '%s%s%s%s</div>\n'%(self.name(),system,lock_text,menu)

        D = self.__notebook.defaults()
        ncols = D['word_wrap_cols']
        s += '<div class="worksheet_cell_list" id="worksheet_cell_list">\n'
        for i in range(n):
            cell = self.__cells[i]
            s += cell.html(ncols,do_print=do_print) + '\n'

        s += '\n</div>\n'
        s += '\n<div class="insert_new_cell" id="insert_last_cell" onmousedown="insert_new_cell_after(cell_id_list[cell_id_list.length-1]);"> </div>\n'
        s += '<div class="worksheet_bottom_padding"></div>\n'

        if not do_print:
            s += '<script language=javascript>cell_id_list=%s;\n'%self.compute_cell_id_list()
            s += 'for(i=0;i<cell_id_list.length;i++) prettify_cell(cell_id_list[i]);</script>\n'
        else:
            s += '<script language=javascript>jsMath.ProcessBeforeShowing();</script>\n'

        if not do_print and confirm_before_leave:
            s += """<script type="text/javascript">
            window.onbeforeunload = confirmBrowseAway;
            function confirmBrowseAway()
            {
            return "Unsubmitted cells will be lost.";
            }
            </script>
            """
        return s

    def show_all(self):
        for C in self.__cells:
            try:
                C.set_cell_output_type('wrap')
            except AttributeError:   # for backwards compatibility
                pass

    def hide_all(self):
        for C in self.__cells:
            try:
                C.set_cell_output_type('hidden')
            except AttributeError:
                pass


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


def extract_text_before_first_compute_cell(text):
    """
    OUTPUT:
        Everything in text up to the first \{\{\{.
    """
    i = text.find('{{{')
    if i == -1:
        return text
    return text[:i]

def extract_first_compute_cell(text):
    """
    INPUT:
        a block of wiki-like marked up text
    OUTPUT:
        input -- string, the input text
        output -- string, the output text
        graphics -- string, text that describes any embedded graphics
        end -- integer, first position after }}} in text.
    """
    # Find the input block
    i = text.find('{{{')
    if i == -1:
        raise EOFError
    j = text.find('\n}}}')
    if j <= i:
        j = len(text)
    k = text.find('\n///')
    if k == -1 or k > j:
        input = text[i+3:j]
        output = ''
        graphics = ''
    else:
        input = text[i+3:k].strip()
        # Find the graphics block, if there is one.
        l = text[k+4:].find('\n///')
        if l != -1 and l+k+4 < j:
            graphics = text[l+k+4+3:j]
        else:
            graphics = ''
            l = j
        output = text[k+4:l].strip()
    return input.strip(), output, graphics, j+4

