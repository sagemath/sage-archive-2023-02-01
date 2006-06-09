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

from sage.ext.sage_object  import SageObject
from sage.interfaces.sage0 import Sage
from sage.misc.preparser   import preparse_file
from sage.misc.misc        import verbose, DOT_SAGE

from cell import Cell

INTERRUPT_TRIES = 20
import notebook as _notebook

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
        self.append_new_cell()

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
            self.__sage = Sage(logfile='%s/sage.log'%self.directory())
            S = self.__sage
            S.eval('__DIR__="%s/"; DIR=__DIR__'%self.DIR())
            S.eval('import sage.server.support as _support_')
            S.eval('__SAGENB__globals = set(globals().keys())')
            object_directory = os.path.abspath(self.__notebook.object_directory())
            verbose(object_directory)
            S.eval('_support_.init("%s", globals())'%object_directory)
            S.eval('print ""')
            S.eval('print ""')
            S.eval('print ""')

            A = self.attached_files()
            for F in A.iterkeys():
                A[F] = 0  # expire all

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
        if C.worksheet() != self:
            raise ValueError, "C must be have self as worksheet."
        if not (C in self.__queue):
            self.__queue.append(C)
        self.start_next_comp()

    def start_next_comp(self):
        if len(self.__queue) == 0:
            return
        if self.__comp_is_running:
            return

        self.__comp_is_running = True
        try:
            del self.__variables   # cache could change...
        except AttributeError:
            pass

        C = self.__queue[0]
        if C.interrupted():
            # don't actually compute
            return

        D = C.directory()
        S = self.sage()


        S._eval_line('os.chdir("%s")'%os.path.abspath(D))
        if C.time():
            S._eval_line('__t__=cputime(); __w__=walltime()')

        tmp = '%s/_temp_.py'%self.directory()
        input = self.preparse_input(C.input_text(), C.completions())
        open(tmp,'w').write(input + '\nprint ""\n')   # the print "" is a hack.
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
        E.timeout = 0.1
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

    def preparse(self, s):
        return preparse_file(s, magic=False, do_time=False, ignore_prompts=False)

    def load_any_changed_attached_files(self, s):
        """
        Modify s by prepending any necessary load commands
        corresponding to attached files that have changed.
        """
        A = self.attached_files()
        for F, tm in A.iteritems():
            new_tm = os.path.getmtime(F)
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
                    d = '{' + ','.join(["'%s':%s"%(v,v) for v in self.variables(with_types=False)]) + '}'
                    t = '_support_.save_session(%s, "%s")'%(d, filename)
                else:
                    t = 'load_session(locals(), "%s")'%filename

            elif t[:5] == 'save ':
                t = self._save_objects(t[5:])

            u.append(t)

        return '\n'.join(u)

    def check_for_system_switching(self, s):
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
            s = 'print %s.eval(r"""%s""")'%(sys,s[i+1:].replace('\n',' '))
            print s
            return True, s
        return False, s

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
        input = '\n'.join(input) + '\n'
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
            return []
        if self.__comp_is_running:
            return []  # can't get the info.
        if with_types:
            try:
                return self.__variables
            except AttributeError:
                pass
            cmd = '_support_.variables(True)'
            S = self.sage()
            v = S.eval(cmd)[1:-1]
            v = v.replace("<type '","").replace("<class '","").replace("'>","").replace('"','')
        else:
            S = self.sage()
            cmd = '_support_.variables(False)'
            v = S.eval(cmd)[1:-1]
        w = v.split(',')
        w.sort()
        if w[0] == '':
            del w[0]
        if with_types:
            self.__variables = w
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

        s += '<span class="worksheet_title">%s</span>\n'%self.name()
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



