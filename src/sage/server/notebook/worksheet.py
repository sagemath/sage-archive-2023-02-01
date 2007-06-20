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

import re
whitespace = re.compile('\s')  # Match any whitespace character
non_whitespace = re.compile('\S')

import pexpect

from sage.structure.sage_object  import load, save
from sage.interfaces.sage0 import Sage
from sage.misc.preparser   import preparse_file
from sage.misc.misc        import alarm, cancel_alarm, verbose, DOT_SAGE, walltime
import sage.server.support as support

import worksheet_conf

import twist

from cell import Cell, TextCell

INTERRUPT_TRIES = 3
INITIAL_NUM_CELLS = 1


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
def initialized_sage(server, ulimit):
    S = Sage(server=server, ulimit=ulimit, maxread = 1, python=True, verbose_start=True)
    S._start(block_during_init=False)
    E = S.expect()
    E.sendline('\n')
    E.expect('>>>')
    cmd = 'from sage.all_notebook import *;'
    cmd += 'import sage.server.support as _support_; '
    E.sendline(cmd)
    return S



_a_sage = None
def init_sage_prestart(server, ulimit):
    global _a_sage
    _a_sage = initialized_sage(server, ulimit)

def one_prestarted_sage(server, ulimit):
    global _a_sage
    X = _a_sage
    if multisession:
        init_sage_prestart(server, ulimit)
    return X

import notebook as _notebook
def worksheet_filename(name, owner):
    return owner + '/' + _notebook.clean_name(name)

class Worksheet:
    def __init__(self, name, dirname, notebook, system, owner, docbrowser=False):

        # Record the basic properties of the worksheet
        self.__system   = system
        self.__owner         = owner
        self.__viewers       = []
        self.__collaborators = [owner]
        self.__docbrowser = docbrowser

        # Initialize the cell id counter.
        self.__next_id = 0

        self.set_name(name)

        # set the directory in which the worksheet files will be stored.
        # We also add the hash of the name, since the cleaned name loses info, e.g.,
        # it could be all _'s if all characters are funny.
        filename = '%s/%s'%(owner, dirname)
        self.__filename = filename
        self.__dir = '%s/%s'%(notebook.worksheet_directory(), filename)

        self.clear()

    def __cmp__(self, other):
        try:
            return cmp(self.name(), other.name())
        except AttributeError:
            return cmp(type(self), type(other))

    def __repr__(self):
        return str(self.__cells)

    def __len__(self):
        return len(self.__cells)

    def docbrowser(self):
        try:
            return self.__docbrowser
        except AttributeError:
            return False

    ##########################################################
    # Configuration
    ##########################################################
    def conf(self):
        try:
            return self.__conf
        except AttributeError:
            file = '%s/conf.sobj'%self.directory()
            if os.path.exists(file):
                C = load(file)
            else:
                C = worksheet_conf.WorksheetConfiguration()
            self.__conf = C
            return C

    ##########################################################
    # Basic properties
    # TODO: Need to create a worksheet configuration object
    ##########################################################
    def collaborators(self):
        return self.__collaborators

    def viewers(self):
        return self.__viewers

    def delete_notebook_specific_data(self):
        self.__attached = {}
        self.__collaborators = [self.owner()]
        self.__viewer = []

    def name(self):
        return self.__name

    def set_name(self, name):
        if len(name.strip()) == 0:
            name = 'Untitled'
        self.__name = name

    def set_filename_without_owner(self, nm):
        filename = '%s/%s'%(self.owner(), nm)
        self.set_filename(filename)

    def set_filename(self, filename):
        old_filename = self.__filename
        self.__filename = filename
        self.__dir = '%s/%s'%(self.notebook().worksheet_directory(), filename)
        self.notebook().change_worksheet_key(old_filename, filename)

    def filename(self):
        return self.__filename

    def filename_without_owner(self):
        """
        Return the part of the worksheet filename after the last /, i.e.,
        without any information about the owner of this worksheet.
        """
        return os.path.split(self.__filename)[-1]

    def directory(self):
        return self.__dir

    def data_directory(self):
        return self.directory() + '/data/'

    def notebook(self):
        return twist.notebook

    def DIR(self):
        return self.notebook().DIR()

    def system(self):
        try:
            return self.__system
        except AttributeError:
            self.__system = None
            return None

    def set_system(self, system=None):
        self.__system = system

    ##########################################################
    # Owner/viewer/user management
    ##########################################################

    def owner(self):
        try:
            return self.__owner
        except AttributeError:
            self.__owner = 'pub'
            return 'pub'

    def set_owner(self, owner):
        self.__owner = owner
        if not owner in self.__collaborators:
            self.__collaborators.append(owner)

    def user_is_only_viewer(self, user):
        try:
            return user in self.__viewers
        except AttributeError:
            return False

    def user_is_viewer(self, user):
        try:
            return user in self.__viewers or user in self.__collaborators
        except AttributeError:
            return True

    def user_is_collaborator(self, user):
        try:
            return user in self.__collaborators
        except AttributeError:
            return True


    def add_viewer(self, user):
        try:
            if not user in self.__viewers:
                self.__viewers.append(user)
        except AttributeError:
            self.__viewers = [user]

    def add_collaborator(self, user):
        try:
            if not user in self.__collaborators:
                self.__collaborators.append(user)
        except AttributeError:
            self.__collaborators = [user]

    ##########################################################
    # Saving
    ##########################################################
    def save(self):
        path = os.path.abspath(self.__dir)
        open(path + '/worksheet.txt','w').write(self.edit_text())
        save(self.conf(), path + '/conf.sobj')

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

    ##########################################################
    # Exporting the worksheet in plain text command-line format
    ##########################################################
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

    def input_text(self):
        """
        Return text version of the input to the worksheet.
        """
        return '\n\n---\n\n'.join([C.input_text() for C in self.__cells])

    ##########################################################
    # Editing the worksheet in plain text format (export and import)
    ##########################################################
    def edit_text(self):
        """
        Returns a plain-text version of the worksheet with \{\{\{\}\}\} wiki-formatting,
        suitable for hand editing.
        """
        s = self.name() + '\n'

        for C in self.__cells:
            t = C.edit_text().strip()
            if t != '':
                s += '\n\n' + t
        return s

    def edit_save(self, text):
        text.replace('\r\n','\n')
        name, i = extract_name(text)
        text = text[i:]
        self.set_name(name)

        data = []
        while True:
            plain_text = extract_text_before_first_compute_cell(text).strip()
            if len(plain_text) > 0:
                T = plain_text
                data.append(('plain', T))
            try:
                meta, input, output, i = extract_first_compute_cell(text)
                data.append(('compute', (meta,input,output)))
            except EOFError, msg:
                print msg
                break
            text = text[i:]

        ids = set([x[0]['id'] for typ, x in data if typ == 'compute' and  x[0].has_key('id')])
        used_ids = set([])

        cells = []
        for typ, T in data:
            if typ == 'plain':
                if len(T) > 0:
                    id = next_available_id(ids)
                    ids.add(id)
                    cells.append(self._new_text_cell(T, id=id))
                    used_ids.add(id)
            elif typ == 'compute':
                meta, input, output = T
                if meta.has_key('id'):
                    id = meta['id']
                    if id in used_ids:
                        # In this case don't reuse, since ids for must be unique.
                        id = next_available_id(ids)
                    html = True
                else:
                    id = next_available_id(ids)
                    ids.add(id)
                    html = False
                used_ids.add(id)
                C = self.get_cell_with_id(id = id)
                if isinstance(C, TextCell):
                    C = self._new_cell(id)
                C.set_input_text(input)
                C.set_output_text(output, '')
                if html:
                    print C.directory()
                    C.update_html_output()
                cells.append(C)

        if len(cells) == 0:   # there must be at least one cell.
            cells = [self._new_cell()]

        self.set_cell_counter()

        self.__cells = cells

    ##########################################################
    # HTML rendering of the wholea worksheet
    ##########################################################
    def html(self, include_title=True, do_print=False,
             confirm_before_leave=False, read_only=False):
        s = ''
        if include_title:
            s += self.html_title()
            s += self.html_menu()

        s += self.html_worksheet_body(do_print=do_print)

        if do_print:
            s += self.javascript_for_jsmath_rendering()
        else:
            s += self.javascript_for_being_active_worksheet()

        if not do_print and confirm_before_leave:
            s += self.javascript_confirm_before_leave()

        return s

    def html_title(self):
        name = self.name()

        s = ''
        s += '<div class="worksheet_title">'
        s += '<a id="worksheet_title" class="worksheet_title" onClick="rename_worksheet(); return false;" title="Click to rename this worksheet">%s</a>'%(name)
        s += '&nbsp;'*5 + self.html_time_last_edited()
        s += '</div>'
        return s

    def html_menu(self):
        name = self.filename()
        if self.computing():
            interrupt_class = "interrupt"
        else:
            interrupt_class = "interrupt_grey"
        vbar = '<span class="vbar"></span>'
        menu = ''

        menu  += '  <span class="worksheet_control_commands">'
        menu += '    <a class="%s" onClick="interrupt()" id="interrupt">Interrupt</a>'%interrupt_class + vbar
        menu += '    <a class="restart_sage" onClick="restart_sage()" id="restart_sage">Restart</a>' +vbar
        menu += '    <a class="plain_text" href="%s">Edit</a>'%self.worksheet_command('edit') + vbar
        menu += '    <a class="doctest_text" href="%s">Text</a>'%self.worksheet_command('plain') + vbar
        menu += '    <a class="doctest_text" href="%s">Preview</a>'%self.worksheet_command('print') + vbar
        menu += '    <a class="evaluate" onClick="evaluate_all()">Eval All</a>' + vbar
        menu += '    <a class="hide" onClick="hide_all()">Hide</a>/<a class="hide" onClick="show_all()">Show</a>' + vbar
        menu += '    <a class="slide_mode" onClick="slide_mode()">Focus</a>' + vbar

        filename = os.path.split(self.filename())[-1]
        download_name = _notebook.clean_name(self.name())

        menu += '    <a class="download_sws" href="download/%s.sws">Download</a>'%download_name + vbar
        menu += '    <a class="delete" href="delete">Delete</a>'
        menu += '  </span>'

        return menu

    def html_worksheet_body(self, do_print):
        n = len(self.__cells)

        s = ''
        D = self.notebook().conf()
        ncols = D['word_wrap_cols']
        s += '<div class="worksheet_cell_list" id="worksheet_cell_list">\n'
        for i in range(n):
            cell = self.__cells[i]
            s += cell.html(ncols, do_print=do_print) + '\n'

        s += '\n</div>\n'
        s += '\n<div class="insert_new_cell" id="insert_last_cell" onmousedown="insert_new_cell_after(cell_id_list[cell_id_list.length-1]);"> </div>\n'
        s += '<div class="worksheet_bottom_padding"></div>\n'
        return s

    def javascript_for_being_active_worksheet(self):
        s =  '<script language=javascript>cell_id_list=%s;\n'%self.compute_cell_id_list()
        s += 'for(i=0;i<cell_id_list.length;i++) prettify_cell(cell_id_list[i]);</script>\n'
        return s

    def javascript_for_jsmath_rendering(self):
        return '<script language=javascript>jsMath.ProcessBeforeShowing();</script>\n'

    def javascript_confirm_before_leave(self):
        return """<script type="text/javascript">
            window.onbeforeunload = confirmBrowseAway;
            function confirmBrowseAway()
            {
            return "Unsubmitted cells will be lost.";
            }
            </script>
            """



    ##########################################################
    # Last edited
    ##########################################################
    def last_edited(self):
        try:
            return self.__last_edited[0]
        except AttributeError:
            t = time.time()
            self.__last_edited = (t, self.owner())
            return t

    def last_to_edit(self):
        try:
            return self.__last_edited[1]
        except AttributeError:
            return self.owner()

    def record_edit(self, user):
        self.__last_edited = (time.time(), user)

    def time_since_last_edited(self):
        return time.time() - self.last_edited()

    def html_time_since_last_edited(self):
        t = self.time_since_last_edited()
        return convert_seconds_to_meaningful_time_span(t)

    def html_time_last_edited(self):
        lt = time.localtime(self.last_edited())
        tm = time.strftime('%B %d, %Y %I:%M %p', lt)
        who = self.last_to_edit()
        t = '<span class="lastedit">last edited on %s by %s</span>'%(tm, who)
        return t

    ##########################################################
    # Managing cells and groups of cells in this worksheet
    ##########################################################

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

    ##########################################################
    # Managing whether computing is happening: stop, start, clear, etc.
    ##########################################################
    def clear(self):
        self.__comp_is_running = False
        self.__queue = []
        self.__cells = [ ]
        for i in range(INITIAL_NUM_CELLS):
            self.append_new_cell()

    def computing(self):
        try:
            return self.__comp_is_running
        except AttributeError:
            return False

    def set_not_computing(self):
        self.__comp_is_running = False
        self.__queue = []

    def quit(self):
        try:
            S = self.__sage
        except AttributeError:
            # no sage running anyways!
            return

        try:
            pid = S._expect.pid
            print "PID = ", pid
            os.killpg(pid, 9)
            os.kill(pid, 9)
            S._expect = None
        except AttributeError, msg:
            print "WARNING: %s"%msg
        except Exception, msg:
            print msg
            print "WARNING: Error deleting SAGE object!"

        try:
            os.kill(pid, 9)
        except:
            pass

        del self.__sage

        # We do this to avoid getting a stale SAGE that uses old code.
        self.clear_queue()


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
        object_directory = os.path.abspath(self.notebook().object_directory())
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
        self.__sage = one_prestarted_sage(server = self.notebook().get_server(),
                                          ulimit = self.notebook().get_ulimit())
        verbose("Initializing SAGE.")
        os.environ['PAGER'] = 'cat'
        self.__next_block_id = 0
        self.initialize_sage()
        return self.__sage

    def start_next_comp(self):
        if len(self.__queue) == 0:
            return

        if self.__comp_is_running:
            self._record_that_we_are_computing()
            return

        C = self.__queue[0]
        if C.interrupted():
            # don't actually compute
            return

        D = C.directory()
        if not C.introspect():
            I = C.input_text().strip()
            if I in ['restart', 'quit', 'exit']:
                self.restart_sage()
                S = self.system()
                if S is None: S = 'SAGE'
                C.set_output_text('Exited %s process'%S,'')
                return
            if I.startswith('%time'):
                C.do_time()
                I = after_first_word(I).lstrip()
            elif first_word(I) == 'time':
                C.do_time()
                I = after_first_word(I).lstrip()
        else:
            I = C.introspect()[0]

        S = self.sage()

        id = self.next_block_id()

        # prevent directory disappear problems
        dir = self.directory()
        if not os.path.exists('%s/code'%dir):
            os.makedirs('%s/code'%dir)
        if not os.path.exists('%s/cells'%dir):
            os.makedirs('%s/cells'%dir)
        tmp = '%s/code/%s.py'%(dir, id)

        absD = os.path.abspath(D)
        input = 'os.chdir("%s")\n'%absD

        # TODOss
        os.system('chmod -R a+rw "%s"'%absD)

        if C.time():
            input += '__SAGE_t__=cputime()\n__SAGE_w__=walltime()\n'
        if I.endswith('?'):
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
            self._record_that_we_are_computing()
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
            #history = "Worksheet '%s' (%s)\n"%(self.name(), time.strftime("%Y-%m-%d at %H:%M",time.localtime(time.time())))
            #history += C.edit_text(ncols=HISTORY_NCOLS, prompts=False,
            #                        max_out=HISTORY_MAX_OUTPUT)
            #self.notebook().add_to_history(history)

        return 'd', C

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
        self.quit()

        self.__sage = initialized_sage(server = self.notebook().get_server(),
                                       ulimit = self.notebook().get_ulimit())
        self.initialize_sage()
        self._enqueue_auto_cells()
        self.start_next_comp()


    def worksheet_command(self, cmd):
        return '/home/%s/%s'%(self.filename(), cmd)

    ##########################################################
    # Idle timeout
    ##########################################################
    def quit_if_idle(self, timeout):
        """
        Quit the worksheet process if it has been "idle" for more than timeout seconds,
        where idle is by definition that the worksheet has not reported back that it
        is actually computing.  I.e., an ignored worksheet process (since the user closed
        their browser) is also considered idle, even if code is running.
        """
        if self.time_idle() > timeout:
            print "Quitting idle or ignored worksheet process for '%s'."%self.name()
            self.quit()

    def time_idle(self):
        return walltime() - self.last_compute_walltime()

    def last_compute_walltime(self):
        try:
            return self.__last_compute_walltime
        except AttributeError:
            t = walltime()
            self.__last_compute_walltime = t
            return t

    def _record_that_we_are_computing(self, username=None):
        self.__last_compute_walltime = walltime()
        if username:
            self.record_edit(username)

    ##########################################################
    # Enqueuing cells
    ##########################################################
    def queue(self):
        return list(self.__queue)

    def queue_id_list(self):
        return [c.id() for c in self.__queue]

    def _enqueue_auto(self):
        for c in self.__cells:
            if c.is_auto_cell():
                self.__queue.append(c)


    def enqueue(self, C, username=None):
        self._record_that_we_are_computing(username)
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

    def _enqueue_auto_cells(self):
        for c in self.__cells:
            if c.is_auto_cell():
                self.enqueue(c)

    def set_cell_counter(self):
        self.__next_id = 1 + max([C.id() for C in self.__cells])

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

    def append(self, L):
        self.__cells.append(L)

    ##########################################################
    # Accessing existing cells
    ##########################################################
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

    def is_last_id_and_previous_is_nonempty(self, id):
        if self.__cells[-1].id() != id:
            return False
        if len(self.__cells) == 1:
            return False
        if len(self.__cells[-2].output_text(ncols=0)) == 0:
            return False
        return True


    ##########################################################
    # (Tab) Completions
    ##########################################################
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
        return format_completions_as_html(id, rows)

    ##########################################################
    # Processing of input and output to worksheet process.
    ##########################################################
    def preparse_input(self, input, C):
        C.set_is_html(False)
        introspect = C.introspect()
        if introspect:
            input = self.preparse_introspection_input(input, C, introspect)
        else:
            switched, input = self.check_for_system_switching(input, C)
            if not switched:
                input = self.preparse_nonswitched_input(input)
            input += '\n'
        return input

    def preparse_introspection_input(self, input, C, introspect):
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
        if before_prompt.endswith('??'):
            input = self._get_last_identifier(before_prompt[:-2])
            input = 'print _support_.source_code("%s", globals())'%input
        elif before_prompt.endswith('?'):
            input = self._get_last_identifier(before_prompt[:-1])
            input = 'print _support_.docstring("%s", globals())'%input
        else:
            input = self._get_last_identifier(before_prompt)
            C._word_being_completed = input
            input = 'print "\\n".join(_support_.completions("%s", globals()))'%(input)
        return input

    def preparse_nonswitched_input(self, input):
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
            if not t.startswith('def '):
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
        return input


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

    ##########################################################
    # Loading and attaching files
    ##########################################################
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
                if not filename.endswith('.py') and not filename.endswith('.sage') and \
                       not filename.endswith('.sobj') and not os.path.exists(filename):
                    if os.path.exists(filename + '.sage'):
                        filename = filename + '.sage'
                    elif os.path.exists(filename + '.py'):
                        filename = filename + '.py'
                    elif os.path.exists(filename + '.sobj'):
                        filename = filename + '.sobj'
            a.append(filename)
        return a

    def _load_file(self, filename, files_seen_so_far, this_file):
        if filename.endswith('.sobj'):
            name = os.path.splitext(filename)[0]
            name = os.path.split(name)[-1]
            return '%s = load("%s");'%(name, filename)

        if filename in files_seen_so_far:
            t = "print 'WARNING: Not loading %s -- would create recursive load'"%filename
        try:
            F = open(filename).read()
        except IOError:
            return "print 'Error loading %s -- file not found'"%filename
        else:
            if filename.endswith('.py'):
                t = F
            elif filename.endswith('.sage'):
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
            if t.startswith('load '):
                z = ''
                for filename in self._normalized_filenames(after_first_word(t)):
                    z += self._load_file(filename, files_seen_so_far, this_file) + '\n'
                t = z

            elif t.startswith('attach '):
                z = ''
                for filename in self._normalized_filenames(after_first_word(t)):
                    if not os.path.exists(filename):
                        z += "print 'Error attaching %s -- file not found'\n"%filename
                    else:
                        self.attach(filename)
                        z += self._load_file(filename, files_seen_so_far, this_file) + '\n'
                t = z

            elif t.startswith('detach '):
                for filename in self._normalized_filenames(after_first_word(t)):
                    self.detach(filename)
                t = ''

            elif t.startswith('save_session') or t.startswith('load_session'):
                F = after_first_word(t).strip().strip('(').strip(')').strip("'").strip('"').split(',')[0]
                if len(F) == 0:
                    filename = self.__filename
                else:
                    filename = F
                if t.startswith('save'):
                    t = '_support_.save_session("%s")'%filename
                else:
                    t = 'load_session(locals(), "%s")'%filename

            elif t.startswith('save '):
                t = self._save_objects(after_first_word(t))

            u.append(t)

        return '\n'.join(u)

    def _eval_cmd(self, system, cmd):
        cmd = cmd.replace("'", "\\u0027")
        return "print _support_.syseval(%s, ur'''%s''')"%(system, cmd)

    ##########################################################
    # Parsing the %sagex, %jsmath, %python, etc., extension.
    ##########################################################

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
            if s.startswith('%sage'):
                s = after_first_word(s).lstrip()
                z = s
            else:
                return True, self._eval_cmd(self.__system, s)

        if len(s) == 0 or s[0] != '%':
            return False, z
        if s.startswith('%hide'):
            t = after_first_word(s).lstrip()
            if len(t) == 0 or t[0] != '%':
                return False, t
            s = t
        if s.startswith('%save_server'):
            self.notebook().save()
            t = after_first_word(s).lstrip()
            if len(t) == 0 or t[0] != '%':
                return False, t
            s = t
        if s.startswith("%pyrex") or s.startswith("%sagex"):  # a block of Sagex code.
            return True, self.sagex_import(after_first_word(s).lstrip(), C)

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

    ##########################################################
    # List of attached files.
    ##########################################################
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

    ##########################################################
    # Showing and hiding all cells
    ##########################################################
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
        if I2.startswith('sage:') or I2.startswith('>>>'):
            do_strip = True
            break
    if not do_strip:
        return s
    s = ''
    for I in t:
        I2 = I.lstrip()
        if I2.startswith('sage:'):
            s += after_first_word(I2).lstrip() + '\n'
        elif I2.startswith('>>>'):
            s += after_first_word(I2).lstrip() + '\n'
        elif I2.startswith('...'):
            s += after_first_word(I2) + '\n'
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
        meta -- meta information about the cell (as a dictionary)
        input -- string, the input text
        output -- string, the output text
        end -- integer, first position after }}} in text.
    """
    # Find the input block
    i = text.find('{{{')
    if i == -1:
        raise EOFError
    j = text[i:].find('\n')
    if j == -1:
        raise EOFError
    k = text[i:].find('|')
    if k != -1 and k < j:
        try:
            meta = dictify(text[i+3:i+k])
        except TypeError:
            meta = {}
        i += k + 1
    else:
        meta = {}
        i += 3

    j = text[i:].find('\n}}}')
    if j == -1:
        j = len(text)
    else:
        j += i
    k = text[i:].find('\n///')
    if k == -1 or k+i > j:
        input = text[i:j]
        output = ''
    else:
        input = text[i:i+k].strip()
        output = text[i+k+4:j].strip()

    return meta, input.strip(), output, j+4

def after_first_word(s):
    """
    Return everything after the first whitespace in the string s.
    Returns the empty string if there is nothing after the
    first whitespace.

    INPUT:
        s -- string
    OUTPUT:
        a string
    """
    i = whitespace.search(s)
    if i is None:
        return ''
    return s[i.start()+1:]

def first_word(s):
    i = whitespace.search(s)
    if i is None:
        return s
    return s[:i.start()]



def format_completions_as_html(cell_id, completions):
    if len(completions) == 0:
        return ''
    lists = []

    # compute the width of each column
    column_width = []
    for i in range(len(completions[0])):
        column_width.append(max([len(x[i]) for x in completions if i < len(x)]))

    for i in range(len(completions)):
        row = completions[i]
        for j in range(len(row)):
            if len(lists) <= j:
                lists.append([])
            cell = """
<li id='completion%s_%s_%s' class='completion_menu_two'>
<a onClick='do_replacement(%s, "%s"); return false;'
   onMouseOver='this.focus(); select_replacement(%s,%s);'
>%s</a>
</li>"""%(cell_id, i, j, cell_id, row[j], i,j,
         row[j])
         #row[j] + '&nbsp;'*(column_width[j]-len(row[j])) )

            lists[j].append(cell)

    grid = "<ul class='completion_menu_one'>"
    for L in lists:
        s = "\n   ".join(L)
        grid += "\n <li class='completion_menu_one'>\n  <ul class='completion_menu_two'>\n%s\n  </ul>\n </li>"%s

    return grid + "\n</ul>"


def extract_name(text):
    # The first line is the title
    i = non_whitespace.search(text)
    if i is None:
        name = 'Untitled'
        n = 0
    else:
        i = i.start()
        j = text[i:].find('\n')
        if j != -1:
            name = text[i:i+j]
            n = j+1
        else:
            name = text[i:]
            n = len(text)-1
    return name.strip(), n


def dictify(s):
    """
    INPUT:
        s -- a string like 'in=5, out=7'
    OUTPUT:
        dict -- such as {'in':5, 'out':7}
    """
    w = []
    try:
        for v in s.split(','):
            a, b = v.strip().split('=')
            try:
                b = eval(b)
            except:
                pass
            w.append([a, b])
    except ValueError:
        return {}
    return dict(w)


def next_available_id(v):
    """
    Return smallest positive integer not in v.
    """
    i = 0
    while i in v:
        i += 1
    return i


def convert_seconds_to_meaningful_time_span(t):
    if t < 60:
        return "%d seconds"%t
    if t < 3600:
        return "%d minutes"%(t/60)
    if t < 3600*24:
        return "%d hours"%(t/3600)
    return "%d days"%(t/(3600*24))
