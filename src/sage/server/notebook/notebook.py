"""
The SAGE Notebook object
"""

#############################################################################
#       Copyright (C) 2006, 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

import os
import shutil
import socket
import re           # regular expressions

# SAGE libraries
from   sage.structure.sage_object import SageObject, load
from   sage.misc.viewer     import browser
from   sage.misc.misc       import alarm, cancel_alarm
from   sage.server.misc import print_open_msg

# SAGE Notebook
import css          # style
import js           # javascript
import server       # web server
import worksheet    # individual worksheets (which make up a notebook)
import config       # internal configuration stuff (currently, just keycodes)
import keyboards    # keyboard layouts

MAX_HISTORY_LENGTH = 500
WRAP_NCOLS = 80

PUBLIC_USER = 'pub'

JSMATH = True

def open_page(address, port):
    cmd = '%s http://%s:%s 1>&2 >/dev/null &'%(browser(), address, port)
    os.system(cmd)

class Notebook(SageObject):
    def __init__(self,
                 dir='sage_notebook',
                 system=None,
                 show_debug = False,
                 log_server=False,
                 address='localhost',
                 port=8000,
                 secure=True,
                 server_pool = []):
        self.__dir = dir
        self.__server_pool = server_pool
        self.set_system(system)
        self.__worksheets = {}
        self.__load_defaults()
        self.__filename      = '%s/nb.sobj'%dir
        self.__worksheet_dir = '%s/worksheets'%dir
        self.__object_dir    = '%s/objects'%dir
        self.__makedirs()
        self.__history = []
        self.__history_count = 0
        self.__log_server = log_server #log all POST's and GET's
        self.__server_log = [] #server log list
        self.__show_debug = show_debug
        self.save()
        self.__admins = []

    def _migrate_old(self):
        """
        Migrate all old worksheets, i.e., ones with no owner
        to /pub.
        """
        for w in self.__worksheets.itervalues():
            if not '/' in w.filename():
                print "Moving worksheet ", w.name()
                w.set_owner(PUBLIC_USER)
                self.rename_worksheet_filename(w, w.filename())

    def user_is_admin(self, user):
        # todo -- make this use the password file !!!
        return user in ['a', 'admin']
        try:
            return user in self.__admins
        except AttributeError:
            self.__admins = []
            return False

    def add_admin(self, user):
        try:
            if not user in self.__admins:
                self.__admins.append(user)
        except AttributeError:
            self.__admins = [user]


    def server_pool(self):
        try:
            return self.__server_pool
        except AttributeError:
            self.__server_pool = []
            return []

    def set_server_pool(self, servers):
        self.__server_pool = servers

    def get_ulimit(self):
        try:
            return self.__ulimit
        except AttributeError:
            self.__ulimit = ''
            return ''

    def set_ulimit(self, ulimit):
        self.__ulimit = ulimit

    def get_server(self):
        P = self.server_pool()
        if len(P) == 0:
            return None
        try:
            i = (self.__server_number + 1)%len(P)
        except AttributeError:
            self.__server_number = 0
            i = 0
        return P[i]

    def system(self, username=None):
        try:
            return self.__system
        except AttributeError:
            self.__system = None
            return None

    def set_system(self, system):
        if system == 'sage':
            self.__system = None
        elif system:  # don't change if it is None
            self.__system = system

    def color(self):
        try:
            return self.__color
        except AttributeError:
            self.__color = 'default'
            return self.__color

    def set_color(self,color):
        self.__color = color

    def set_directory(self, dir):
        if dir == self.__dir:
            return
        self.__dir = dir
        self.__filename = '%s/nb.sobj'%dir
        self.__worksheet_dir = '%s/worksheets'%dir
        self.__object_dir = '%s/objects'%dir
        for W in self.__worksheets.itervalues():
            W.set_notebook(self)

    def add_to_history(self, input_text):
        H = self.history()
        H.append(input_text)
        while len(H) > self.max_history_length():
            del H[0]

    def history_count_inc(self):
        self.__history_count += 1

    def history_count(self):
        return self.__history_count

    def server_log(self):
        return self.__server_log

    def log_server(self):
        return self.__log_server

    def set_log_server(self, log_server):
        self.__log_server = log_server

    def history(self):
        try:
            s = self.__history
        except AttributeError:
            self.__history = []
            s = self.__history
        return s

    def history_text(self):
        return '\n\n'.join([H.strip() for H in self.history()])

    def history_html(self):
        t = self.history_text()
        t = t.replace('<','&lt;')
        s = '<head>\n'
        s += '<title>SAGE Input History</title>\n'
        s += '</head>\n'
        s += '<body>\n'
        s += '<pre>' + t + '</pre>\n'
        s += '<a name="bottom"></a>\n'
        s += '<script type="text/javascript"> window.location="#bottom"</script>\n'
        s += '</body>\n'
        return s


    def history_with_start(self, start):
        n = len(start)
        return [x for x in self.history() if x[:n] == start]

    def export_worksheet(self, worksheet_filename, filename):
        W = self.get_worksheet_with_filename(worksheet_filename)
        W.save()
        cmd = 'cd %s && tar -jcf %s.sws "%s" && mv %s.sws ..'%(
            self.__worksheet_dir,
            filename, W.filename(), filename)
        print cmd
        os.system(cmd)

    def plain_text_worksheet_html(self, name, prompts=True):
        W = self.get_worksheet_with_filename(name)
        t = W.plain_text(prompts = prompts)
        t = t.replace('<','&lt;')
        s = '<head>\n'
        s += '<title>SAGE Worksheet: %s</title>\n'%W.name()
        s += '</head>\n'
        s += '<body>\n'
        s += '<h1><a href=".">SAGE Worksheet: %s</a></h1>\n'%W.name()
        s += '<pre>' + t + '</pre>'
        s += '</body>\n'
        return s

    def tmpdir(self):
        d = '%s/tmp'%self.__dir
        if os.path.exists(d):
            os.system('rm -rf "%s"'%d)
        if not os.path.exists(d):
            os.makedirs(d)
        return d

    def import_worksheet(self, filename):
        # TODO -- this is broken -- it does *not* work if you
        # upload a worksheet with the same name as an existing worksheet,
        # though somebody thinks they wrote it to do that.
        # The problem is that changing the worksheet name is not
        # enough -- one must also change the worksheet id.
        # One should get rid of the id stuff, maybe.  For now,
        # we raise an error if the worksheet name is already used.
        if not os.path.exists(filename):
            raise ValueError, "no file %s"%filename
        if filename[-4:] != '.sws':
            raise ValueError, "file %s must have extension sws."%filename
        tmp = self.tmpdir()
        cmd = 'cd %s; tar -jxf %s'%(tmp, os.path.abspath(filename))
        print cmd
        os.system(cmd)
        try:
            D = os.listdir(tmp)[0]
        except IndexError:
            raise ValueError, "invalid worksheet"
        worksheet = load('%s/%s/%s.sobj'%(tmp,D,D), compress=False)
        names = self.worksheet_names()
        if D in names:
            raise ValueError, "Worksheet with given name already defined."
            m = re.match('.*?([0-9]+)$',D)
            if m is None:
                n = 0
            else:
                n = int(m.groups()[0])
            while "%s%d"%(D,n) in names:
                n += 1
            cmd = 'mv %s/%s/%s.sobj %s/%s/%s%d.sobj'%(tmp,D,D,tmp,D,D,n)
            print cmd
            os.system(cmd)
            cmd = 'mv %s/%s %s/%s%d'%(tmp,D,tmp,D,n)
            print cmd
            os.system(cmd)
            D = "%s%d"%(D,n)
            worksheet.set_name(D)
        print D
        S = self.__worksheet_dir
        cmd = 'rm -rf "%s/%s"'%(S,D)
        print cmd
        os.system(cmd)
        cmd = 'mv %s/%s %s/'%(tmp, D, S)
        print cmd
        os.system(cmd)
        new_id = None
        worksheet.set_notebook(self)
        filename = worksheet.filename()
        self.__worksheets[filename] = worksheet
        return worksheet

    # unpickled, no worksheets will think they are
    # being computed, since they clearly aren't (since
    # the server just started).
    def set_not_computing(self):
        for W in self.__worksheets.values():
            W.set_not_computing()

    def set_debug(self,show_debug):
        self.__show_debug = show_debug

    def directory(self):
        if not os.path.exists(self.__dir):
            # prevent "rm -rf" accidents.
            os.makedirs(self.__dir)
        return self.__dir

    def DIR(self):
        """
        Return the absolute path to the directory that contains
        the SAGE Notebook directory.
        """
        P = os.path.abspath('%s/..'%self.__dir)
        if not os.path.exists(P):
            # prevent "rm -rf" accidents.
            os.makedirs(P)
        return P

    def max_history_length(self):
        try:
            return self.__defaults['max_history_length']
        except KeyError:
            return MAX_HISTORY_LENGTH

    def __load_defaults(self):
        # in future this will allow override by a file, and
        # can be set by user via web interface
        self.__defaults = {'cell_input_color':'#0000000',
                           'cell_output_color':'#0000EE',
                           'word_wrap_cols':int(WRAP_NCOLS),
                           'max_history_length':MAX_HISTORY_LENGTH}

    def worksheet_directory(self):
        return self.__worksheet_dir

    def object_directory(self):
        O = self.__object_dir
        if not os.path.exists(O):
            os.makedirs(O)
        return O

    def objects(self):
        L = [x[:-5] for x in os.listdir(self.object_directory())]
        L.sort()
        return L

    def object_list_html(self):
        m = max([len(x) for x in self.objects()] + [30])
        s = []
        a = '<a href="/%s.sobj" class="object_name">\n'
        for name in self.objects():
            s.append(a%name + name + '</a>\n')  # '&nbsp;'*(m-len(name)) +
        return '<br>\n'.join(s)

    def defaults(self):
        return self.__defaults

    def authorize(self, auth):
        """
        Returns True if auth is the correct authorization.
        """
        a = self.auth_string()
        if a == ':':
            return True
        return a == auth

    def auth_string(self):
        try:
            return self.__auth
        except AttributeError:
            self.__auth = ":"
        return self.__auth

    def set_auth(self, username, password):
        self.__auth = '%s:%s'%(username, password)

    def __makedirs(self):
        if not os.path.exists(self.__dir):
            os.makedirs(self.__dir)
        if not os.path.exists(self.__worksheet_dir):
            os.makedirs(self.__worksheet_dir)
        if not os.path.exists(self.__object_dir):
            os.makedirs(self.__object_dir)

    def create_new_worksheet(self, worksheet_name, username):
        filename = worksheet.worksheet_filename(worksheet_name, username)
        if self.__worksheets.has_key(filename):
            return self.__worksheets[filename]
        W = worksheet.Worksheet(worksheet_name, self,
                        system = self.system(username), owner=username)
        self.__worksheets[W.filename()] = W
        return W

    def delete_worksheet(self, filename):
        """
        Delete the given worksheet and remove its name from the
        worksheet list.
        """
        if not (filename in self.__worksheets.keys()):
            print self.__worksheets.keys()
            raise KeyError, "Attempt to delete missing worksheet '%s'"%filename
        W = self.__worksheets[filename]
        W.quit()
        cmd = 'rm -rf "%s"'%(W.directory())
        print cmd
        os.system(cmd)

        self.deleted_worksheets()[filename] = W
        del self.__worksheets[filename]

    def deleted_worksheets(self):
        try:
            return self.__deleted_worksheets
        except AttributeError:
            self.__deleted_worksheets = {}
            return self.__deleted_worksheets

    def worksheet_names(self):
        W = self.__worksheets.keys()
        W.sort()
        return W

    def worksheet_html(self, name, do_print=False):
        W = self.get_worksheet_with_filename(name)
        s = '<head>\n'
        s += '<title>SAGE Worksheet: %s</title>\n'%W.name()
        if do_print:
            s += '<script type="text/javascript" src="/javascript/jsmath/jsMath.js"></script>\n'
        s += '<script type="text/javascript" src="/javascript/main.js"></script>\n'
        s += '<link rel=stylesheet href="/css/main.css">\n'
        s += '</head>\n'
        s += '<body>\n'
        if do_print:
            s += '<h1><a href=".">SAGE Worksheet: %s</a></h1>'%W.name()
        s += W.html(include_title=False, do_print=do_print)
        if do_print:
            s += '<script type="text/javascript">jsMath.Process();</script>\n'
        s += '\n</body>\n'
        return s

    def get_worksheets_with_collaborator(self, user):
        return [w for w in self.__worksheets.itervalues() if w.user_is_collaborator(user)]

    def get_worksheet_names_with_collaborator(self, user):
        return [W.name() for W in self.get_worksheets_with_collaborator(user)]

    def get_worksheets_with_viewer(self, user):
        return [w for w in self.__worksheets.itervalues() if w.user_is_viewer(user)]

    def get_worksheets_with_owner(self, owner):
        return [w for w in self.__worksheets.itervalues() if w.owner() == owner]

    def get_worksheets_with_owner_that_are_viewable_by_user(self, owner, user):
        return [w for w in self.get_worksheets_with_owner(owner) if w.user_is_viewer(user)]

    def get_worksheet_names_with_viewer(self, user):
        return [W.name() for W in self.get_worksheets_with_viewer(user)]

    def get_worksheet_with_name(self, name):
        return self.__worksheets[name]

    def get_worksheet_with_id(self, id):
        try:
            id = int(id)
            for W in self.__worksheets.itervalues():
                if W.id() == id:
                    return W
        except ValueError:
            id = str(id).lower()
            for W in self.__worksheets.itervalues():
                if W.name().lower() == id or W.filename().lower() == id:
                    return W
        raise KeyError, 'no worksheet %s'%id

    def get_worksheet_with_filename(self, filename):
        """
        Get the worksheet with given filename.  If there is no such
        worksheet, raise a KeyError.

        INPUT:
            string
        OUTPUT:
            a worksheet or KeyError
        """
        if self.__worksheets.has_key(filename):
            return self.__worksheets[filename]
        if not filename is None:
            for W in self.__worksheets.itervalues():
                if W.filename() == filename:
                    return W
        raise KeyError, "no such worksheet %s"%filename

    ###########################################################
    def html_worksheet_list_for_user(self, user):
        add_new_worksheet_menu = """
             <div class="add_new_worksheet_menu" id="add_worksheet_menu">
             <input id="new_worksheet_box" class="add_new_worksheet_menu"
                    onKeyPress="if(is_submit(event)) process_new_worksheet_menu_submit();"></input><br>
             <button class="add_new_worksheet_menu"  onClick="process_new_worksheet_menu_submit();">New</button>
             </div>
        """
        W = self.get_worksheets_with_viewer(user)
        s = '<html><body> <ol>\n'

        # This is stupid -- just used for add_new_worksheet_menu -- get rid of this.
        s += '<script type="text/javascript" src="/javascript/main.js"></script>\n'
        s += '<script type="text/javascript">user_name="%s"; </script>'%user

        s += '<br>'*2
        s += add_new_worksheet_menu
        s += '<br>'*2
        s += '<h2>Logged in as: %s</h2>'%user
        s += '<h2>Active Worksheets</h2>'
        s += '<br>'*2
        for w in W:
            s += '<li> <a href="/home/%s">%s</a>\n'%(w.filename(), w.name())
        s += '</body></html>'
        return s



    ###########################################################

    def save(self, filename=None):
        #print "-"*70

        if filename is None:
            F = os.path.abspath(self.__filename)
            try:
                shutil.copy(F, F[:-5] + '-backup.sobj')
            except IOError:
                pass
            F = os.path.abspath(self.__filename)
        else:
            F = os.path.abspath(filename)

        print "Saving notebook to '%s'..."%F
        D, _ = os.path.split(F)
        if not os.path.exists(D):
            os.makedirs(D)
        SageObject.save(self, F, compress=False)
        #print "Press control-C to stop the notebook server."
        #print "-"*70

    def quit(self):
        for W in self.__worksheets.itervalues():
            W.quit()

    def delete_doc_browser_worksheets(self):
        names = self.worksheet_names()
        for n in self.__worksheets.keys():
            if n.startswith('doc_browser'):
                self.delete_worksheet(n)

    def worksheet_list_html(self, current_worksheet, username):
        print current_worksheet, username
        s = []
        names = self.get_worksheet_names_with_viewer(username)
        m = max([len(x) for x in names] + [30])
        for n in names:
            if n.startswith('doc_browser'): continue
            W = self.__worksheets[n]
            if W == current_worksheet:
                cls = 'worksheet_current'
            else:
                cls = 'worksheet_other'
            if W.computing():
                cls += '_computing' # actively computing
            name = W.name()
            name += ' (%s)'%len(W)
            name = name.replace(' ','&nbsp;')
        return '<br>'.join(s)

    def _html_head(self, worksheet_filename, username):
        if worksheet_filename is not None:
            worksheet = self.get_worksheet_with_filename(worksheet_filename)
            head = '\n<title>%s (%s)</title>'%(worksheet.name(), self.directory())
        else:
            head = '\n<title>SAGE Notebook | Welcome</title>'
        head += '\n<script type="text/javascript" src="/javascript/main.js"></script>\n'
        head += '\n<link rel=stylesheet href="/css/main.css" type="text/css">\n'

        if JSMATH:
            head += '<script type="text/javascript">jsMath = {Controls: {cookie: {scale: 115}}}</script>\n'
            head +=' <script type="text/javascript" src="/javascript/jsmath/plugins/noImageFonts.js"></script>\n'
            head += '<script type="text/javascript" src="/javascript/jsmath/jsMath.js"></script>\n'
            head += "<script type='text/javascript'>jsMath.styles['#jsMath_button'] = jsMath.styles['#jsMath_button'].replace('right','left');</script>\n"

        head +=' <script type="text/javascript" src="/javascript/highlight/prettify.js"></script>\n'
        head += '<link rel=stylesheet href="/css/highlight/prettify.css" type="text/css">\n'

        return head

    def _html_body(self, worksheet_filename, show_debug=False, username=''):
        if worksheet_filename is None or worksheet_filename == '':
            main_body = '<div class="worksheet_title">Welcome %s to the SAGE Notebook</div>\n'%username
            if os.path.isfile(self.directory() + "/index.html"):
                splash_file = open(self.directory() + "/index.html")
                main_body+= splash_file.read()
                splash_file.close()
            else:
                dir = os.path.abspath('%s'%self.directory())
                main_body+= "<br>&nbsp;&nbsp;&nbsp;SAGE Notebook running from <tt>%s</tt>."%dir
                main_body+= self.help_window()
                main_body += "&nbsp;&nbsp;&nbsp;Create a file <tt>%s/index.html</tt> to replace this splash page.<br>"%(dir)
            interrupt_class = "interrupt_grey"
            worksheet = None
        else:

            worksheet = self.get_worksheet_with_filename(worksheet_filename)
            if worksheet.computing():
                interrupt_class = "interrupt"
            else:
                interrupt_class = "interrupt_grey"
            main_body = worksheet.html()

##         add_new_worksheet_menu = """
##              <div class="add_new_worksheet_menu" id="add_worksheet_menu">
##              <input id="new_worksheet_box" class="add_new_worksheet_menu"
##                     onKeyPress="if(is_submit(event)) process_new_worksheet_menu_submit();"></input><br>
##              <button class="add_new_worksheet_menu"  onClick="process_new_worksheet_menu_submit();">New</button>
##              </div>
##         """

##         delete_worksheet_menu = """
##              <div class="delete_worksheet_menu" id="delete_worksheet_menu">
##              <input id="delete_worksheet_box" class="delete_worksheet_menu"
##                     onKeyPress="if(is_submit(event)) process_delete_worksheet_menu_submit();"></input>
##              <button class="delete_worksheet_menu" onClick="process_delete_worksheet_menu_submit();">delete</button>
##              &nbsp;&nbsp;&nbsp;<span class="X" onClick="hide_delete_worksheet_menu()">X</span>
##              </div>
##         """

        vbar = '<span class="vbar"></span>'

        body = ''

        body += '<div class="top_control_bar">\n'
        body += '  <span class="banner"><a class="banner" target="_new" href="http://www.sagemath.org">'
        body += '  <img src="/images/sagelogo.png" alt="SAGE"></a></span>\n'
        body += '  <span class="control_commands" id="cell_controls">\n'
        body += '    <a class="help" href="/home/%s">Worksheets</a>'%username + vbar
        body += '    <a class="history_link" onClick="history_window()">History</a>' + vbar
        body += '    <a class="help" onClick="show_help_window()">Help</a>' + vbar
        body += '    <a href="/doc">Documentation</a>' + vbar
        body += '     <a href="/upload" class="upload_worksheet">Upload</a>' + vbar
        body += '     <a href="/logout" class="help">Sign Out</a>'
        body += '  </span>\n'

        #these divs appear in backwards order because they're float:right
        body += '  <div class="hidden" id="slide_controls">\n'
        body += '    <div class="slideshow_control">\n'
        body += '      <a class="slide_arrow" onClick="slide_next()">&gt;</a>\n'
        body += '      <a class="slide_arrow" onClick="slide_last()">&gt;&gt;</a>\n' + vbar
        body += '      <a class="cell_mode" onClick="cell_mode()">Show Full Worksheet</a>\n'
        body += '    </div>\n'
        body += '    <div class="slideshow_progress" id="slideshow_progress" onClick="slide_next()">\n'
        body += '      <div class="slideshow_progress_bar" id="slideshow_progress_bar">&nbsp;</div>\n'
        body += '      <div class="slideshow_progress_text" id="slideshow_progress_text">&nbsp;</div>\n'
        body += '    </div>\n'
        body += '    <div class="slideshow_control">\n'
        body += '      <a class="slide_arrow" onClick="slide_first()">&lt;&lt;</a>\n'
        body += '      <a class="slide_arrow" onClick="slide_prev()">&lt;</a>\n'
        body += '    </div>\n'
        body += '  </div>\n'

        body += '</div>\n'
        body += '\n<div class="worksheet" id="worksheet">\n'
        if self.__show_debug or show_debug:
            body += "<div class='debug_window'>"
            body += "<div class='debug_output'><pre id='debug_output'></pre></div>"
            body += "<textarea rows=5 id='debug_input' class='debug_input' "
            body += " onKeyPress='return debug_keypress(event);' "
            body += " onFocus='debug_focus();' onBlur='debug_blur();'></textarea>"
            body += "</div>"

        body += main_body + '\n</div>\n'

        # The blank space given by '<br>'*15  is needed so the input doesn't get
        # stuck at the bottom of the screen. This could be replaced by a region
        # such that clicking on it creates a new cell at the bottom of the worksheet.
        body += '<br>'*15
        body += '\n</div>\n'

#        body += '<div class="left_pane_bar" id="left_pane_bar" onClick="toggle_left_pane();"></div>\n'
#        body += '<span class="pane" id="left_pane"><table bgcolor="white"><tr><td>\n'
        endpanespan = '</td></tr></table></span>\n'

##         body += '  <div class="worksheets_topbar">'
##         body += '     <b>Worksheets</b> '
##         body += '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a class="left_panel_hide" onClick="toggle_left_pane()" class="worksheets_button" id="worksheets_button">Hide</a>'
##         body += '<br></div>'
##         body +=    add_new_worksheet_menu
##         body += '  <div class="worksheet_list" id="worksheet_list">%s</div>\n'%self.worksheet_list_html(worksheet, username)

        if worksheet is None:
             return body + endpanespan

##         body += '<script type="text/javascript">focus(%s)</script>\n'%(worksheet[0].id())
##         body += '<script type="text/javascript">jsmath_init();</script>\n'

        if worksheet.user_is_only_viewer(username):
            body += '<script type="text/javascript">worksheet_locked=true;</script>'
        else:
            body += '<script type="text/javascript">worksheet_locked=false;</script>'

        if worksheet.computing():
            # Set the update checking back in motion.
            body += '<script type="text/javascript"> active_cell_list = %r; \n'%worksheet.queue_id_list()
            body += 'for(var i = 0; i < active_cell_list.length; i++)'
            body += '    cell_set_running(active_cell_list[i]); \n'
            body += 'start_update_check(); </script>\n'

        #body += '<script type="text/javascript">toggle_left_pane()</script>'
        return body

    def edit_window(self, worksheet):
        """
        Return a window for editing worksheet.

        INPUT:
            worksheet -- a worksheet
        """
        t = worksheet.edit_text()
        t = t.replace('<','&lt;')
        body_html = ''
        body_html += '<h1 class="edit">SAGE Notebook: Editing Worksheet "%s"</h1>\n'%worksheet.name()
        body_html += """<b>Warnings:</b> You cannot undo after you save changes (yet).  All graphics will be deleted when you save.<br><br>"""
        body_html += '<form method="post" action="save" enctype="multipart/form-data">\n'
        body_html += '<input type="submit" value="Save Changes" name="button_save"/>\n'
        #body_html += '<input type="submit" value="Preview" name="button_preview"/>\n'
        body_html += '<input type="submit" value="Cancel" name="button_cancel"/>\n'
        body_html += '<textarea class="edit" id="cell_intext" rows="30" name="textfield">'+t+'</textarea>'
        body_html += '</form>'
        body_html += """The format is as follows: <pre>
Arbitrary HTML
{{{
Input
///
Output
}}}
Arbitrary HTML
{{{
Input
///
Output
}}}
...
</pre>"""

        s = """
        <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
        <html><head><title>SAGE Wiki cell text </title>
        <style type="text/css">

        textarea.edit {
            font-family: courier, monospace;
            font-size:12pt;
            border: 1px solid #8cacbb;
            color: black;
            background-color: white;
            padding: 3px;
            width: 100%%;
            margin-top: 0.5em;
        }
        </style>
        <body>%s
        </body></html>"""%body_html

        return s

    def help_window(self):
        try:
            return self._help_window
        except AttributeError:
            pass
        from tutorial import notebook_help
        s = """
        <br><hr>
        <style>
        div.help_window {
            background-color:white;
            border: 3px solid #3d86d0;
            top: 10ex;
            bottom:10%;
            left:25%;
            right:15%;
            padding:2ex;
        }


        table.help_window {
            background-color:white;
            width:100%;
        }

        td.help_window_cmd {
            background-color: #f5e0aa;
            width:30%;
            padding:1ex;
            font-weight:bold;
        }

        td.help_window_how {
            padding:1ex;
            width:70%;
        }
        </style>
        <h1 align=center><font color='darkred'>SAGE</font> Notebook Quickstart</h1>
        <div class="help_window">

        A <i>worksheet</i> is an ordered list of SAGE calculations with output.
        A <i>session</i> is a worksheet and a set of variables in some state.
        A <i>notebook</i> is a collection of worksheets and saved objects.
        <br>
        <br>
        To get started with SAGE, <a href="doc_browser?/tut/?tut.html">view the tutorial</a>.
        <br><br>

        <table class="help_window">
        """
        for x, y in notebook_help:
            s += '<tr><td class="help_window_cmd">%s</td><td class="help_window_how">%s</td></tr>\n'%(x,y)
        s += '</table></div>'

        s +="""
        <br>
        AUTHORS: Tom Boothby, Alex Clemesha, Bobby Moretti, Yi Qiang, Dorian Ramier, and William Stein<br><br>
        LICENSE: SAGE is <a href="/license.html">GPL-compatible</a>.
        <br>
        """
        self._help_window = s
        return s

    def upload_window(self):
        return """
          <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
          <html>
            <head>
              <title>Upload File</title>
              <style>%s</style>
              <script type='text/javascript'>%s</script>
            </head>
            <body onLoad="if(window.focus) window.focus()">
              <div class="upload_worksheet_menu" id="upload_worksheet_menu">
              <h1><font size=+3 color="darkred">SAGE</font>&nbsp;&nbsp;&nbsp;&nbsp;<font size=+1>Upload your Worksheet</font></h1>
              <hr>
              <form method="POST" action="upload_worksheet"
                    name="upload" enctype="multipart/form-data">
              <table><tr>
              <td>
              Worksheet file:&nbsp&nbsp&nbsp </td>
              <td><input class="upload_worksheet_menu" size="40" type="file" name="fileField" id="upload_worksheet_filename"></input></td>
              </tr>
              <tr><td></td><td></td></tr>
              <tr>
              <td></td><td><input type="button" class="upload_worksheet_menu" value="Upload Worksheet" onClick="form.submit(); window.close();"></td>
              </tr>
              </form><br>
              </div>
            </body>
          </html>
         """%(css.css(self.color()),js.javascript())

    def html(self, worksheet_filename=None, username=None, show_debug=False, admin=False):
        if worksheet_filename is None or worksheet_filename == '':
            worksheet_filename = None
            W = None
        else:
            try:
                W = self.get_worksheet_with_filename(worksheet_filename)
            except KeyError:
                W = None

        head = self._html_head(worksheet_filename=worksheet_filename, username=username)
        body = self._html_body(worksheet_filename=worksheet_filename, username=username, show_debug=show_debug)

        head += '<script type="text/javascript">user_name="%s"; </script>'%username

        if worksheet_filename is not None:
            head += '<script  type="text/javascript">worksheet_filename="%s"; worksheet_name="%s"; </script>'%(worksheet_filename, W.name())

        return """
        <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def _html_authorize(self):
        return """
        <h1>SAGE Notebook Server</h1>
        <div id="mainbody" class="login">Sign in to the SAGE Notebook<br>
        <form>
        <table>
        <tr><td>
          <span class="username">Username:</span></td>
          <td><input name="username" class="username"
                      onKeyPress="if(is_submit(event)) login(username.value, password.value)"></td>
        </tr>
        <tr><td>
           <span class="password">Password:</span></td>
           <td><input name="password" class="username" type="password"
                      onKeyPress="if(is_submit(event)) login(username.value, password.value)"></td>
        </tr>
        <td>&nbsp</td>
        <td>
           <input type='button' onClick="login(username.value,password.value);" value="Sign in">
           </td></table>
                   </form></div>

        """

    def format_completions_as_html(self, cell_id, completions):
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


import sage.interfaces.sage0
import time

def load_notebook(dir, address=None, port=None, secure=None):
    sobj = '%s/nb.sobj'%dir
    if os.path.exists(sobj):
        try:
            nb = load(sobj, compress=False)
        except:
            print "****************************************************************"
            print "  * * * WARNING   * * * WARNING   * * * WARNING   * * * "
            print "WARNING -- failed to load notebook data. Trying the backup file."
            print "****************************************************************"
            try:
                nb = load('%s/nb-backup.sobj'%dir, compress=False)
            except:
                print "Recovering from last op save failed."
                print "Trying save from last startup."
                nb = load('%s/nb-older-backup.sobj'%dir, compress=False)

        nb.delete_doc_browser_worksheets()
        nb.set_directory(dir)
        nb.set_not_computing()
    else:
        nb = Notebook(dir)

    nb.address = address
    nb.port = port
    nb.secure = secure
    return nb

## IMPORTANT!!! If you add any new input variable to notebook,
## you *must* similarly modify the restart_on_crash block
## at the beginning of the definition of notebook!!
def notebook(dir         ='sage_notebook',
             port        = 8000,
             address     = 'localhost',
             open_viewer = True,
             max_tries   = 10,
             username    = None,
             password    = None,
             color       = None,
             system      = None,
             jsmath      = True,
             show_debug  = False,
             splashpage  = True,
             warn        = True,
             ignore_lock = False,
             log_server = False,
             restart_on_crash = False,
             auto_restart = 1800,
             multisession = True):
    r"""
    Start a SAGE notebook web server at the given port.

    INPUT:
        dir -- (default: 'sage_notebook') name of the server directory; your
                sessions are saved in a directory with that name.  If
                you restart the server with that same name then it will
                restart in the state you left it, but with none of the
                variables defined (you have to re-eval blocks).
        port -- (default: 8000) port on computer where the server is served
        address -- (default: 'localhost') address that the server
                   will listen on
        open_viewer -- bool (default: True); if True, pop up a web browser at the URL
        max_tries -- (default: 10) maximum number of ports > port to try in
                     case given port can't be opened.
        username -- user name used for authenticated logins
        password -- password used for authenticated logins
        color -- string or pair of html colors, e.g.,
                    'gmail'
                    'grey'
                    ('#ff0000', '#0000ff')
        system -- (string) default computer algebra system to use for new
                  worksheets, e.g., 'maxima', 'gp', 'axiom', 'mathematica', 'macaulay2',
                  'singular', 'gap', 'octave', 'maple', etc.  (even 'latex'!)
        jsmath -- whether not to enable javascript typset output for math.
        debug -- whether or not to show a javascript debugging window
        splashpage -- whether or not to show a splash page when no worksheet is specified.
                      you can place a file named index.html into the notebook directory that
                      will be shown in place of the default.

        restart_on_crash -- if True (the default is False), the server
                      will be automatically restarted if it crashes in
                      any way.  Use this on a public servers that many
                      people might use, and which might be subjected
                      to intense usage.  NOTE: Log messages are only displayed
                      every 5 seconds in this mode.
        auto_restart -- if restart_on_crash is True, always restart
                      the server every this many seconds.
        multisession -- (default: True) The default is for there to be
                       one sage session for each worksheet.  If this
                       is False, then there is just one global SAGE
                       session, like with Mathematica.

    NOTES:

    When you type \code{notebook(...)}  you start a web server on the
    machine you type that command on.  You don't connect to another
    machine.  So do this if you want to start a SAGE notebook
    accessible from anywhere:

    \begin{enumerate}
    \item Figure out the external address of your server, say
          'www.sagemath.org', for example.
    \item On your server, type
        notebook('mysession', address='www.sagemath.org')
    \item Assuming you have permission to open a port on that
       machine, it will startup and display a URL, e.g.,
           \url{http://www.sagemath.org:8000}
       Note this URL.
    \item Go to any computer in the world (!), or at least
       behind your firewall, and use any web browser to
       visit the above URL.  You're using \sage.
    \end{enumerate}

    \note{There are no security precautions in place \emph{yet}!  If
    you open a server as above, and somebody figures this out, they
    could use their web browser to connect to the same sage session,
    and type something nasty like \code{os.system('cd; rm -rf *')}
    and your home directory would be hosed.   I'll be adding an
    authentication screen in the near future.  In the meantime
    (and even then), you should consider creating a user with
    very limited privileges (e.g., empty home directory).}

    FIREFOX ISSUE:
    If your default web browser if Firefox, then notebook will
    open a copy of Firefox at the given URL.  You should
    definitely set the "open links in new tabs" option in
    Firefox, or you might loose a web page you were looking at.
    To do this, just go to

         Edit --> Preferences --> Tabs

    and in "Open links from other apps" select the middle button
    instead of the bottom button.
    """
    assert 0, "deprecated"
    import worksheet
    worksheet.init_sage_prestart()
    worksheet.multisession = multisession

    if '/' in dir:
	# change current working directory and make the notebook
	# directory a subdirectory of the working directory.
        base, dir = os.path.split(dir)
        os.chdir(base)

    if restart_on_crash:
        # Start a new subprocess
        def f(x):  # format for passing on
            if x is None:
                return 'None'
            elif isinstance(x, str):
                return "'%s'"%x
            else:
                return str(x)
        while True:
            S = sage.interfaces.sage0.Sage()
            time.sleep(1)
            S.eval("from sage.server.notebook.notebook import notebook")
            cmd = "notebook(dir=%s,port=%s, address=%s, open_viewer=%s, max_tries=%s, username=%s, password=%s, color=%s, system=%s, jsmath=%s, show_debug=%s, splashpage=%s, warn=%s, ignore_lock=%s, log_server=%s, restart_on_crash=False, multisession=%s)"%(
                f(dir), f(port), f(address), f(open_viewer), f(max_tries), f(username),
                f(password), f(color), f(system), f(jsmath), f(show_debug), f(splashpage),
                f(warn), f(ignore_lock), f(log_server), f(multisession)
                )
            print cmd
            S._send(cmd)
            tm = 0
            while True:
                s = S._get()[1].strip()
                if len(s) > 0:
                    print s
                if not S.is_running():
                    break
                time.sleep(5)
                tm += 5
                if tm > auto_restart:
                    S.quit()
                    break
            # end while
        # end while
        S.quit()
        return

    if os.path.exists(dir):
        if not os.path.isdir(dir):
            raise RuntimeError, '"%s" is not a valid SAGE notebook directory (it is not even a directory).'%dir
        if not (os.path.exists('%s/nb.sobj'%dir) or os.path.exists('%s/nb-backup.sobj'%dir)):
            raise RuntimeError, '"%s" is not a valid SAGE notebook directory (missing nb.sobj).'%dir
        pidfile = '%s/pid'%dir

        if os.path.exists(pidfile) and not ignore_lock:
            f = file(pidfile)
            try:
                p, oldport = f.readlines()
            except ValueError:
                p = file(pidfile).read()
                oldport = port
            f.close()
            try:
                #This is a hack to check whether or not the process is running.
                os.kill(int(p),0)
                print "\n".join([" This notebook appears to be running with PID %s.  If it is"%p,
                                 " not responding, you will need to kill that process to continue.",
                                 " If another (non-sage) process is running with that PID, call",
                                 " notebook(..., ignore_lock = True, ...). " ])
                if open_viewer:
                   open_page(address, int(oldport))
                return
            except OSError:
                pass

    nb = load_notebook(dir, username, password, color,system, splashpage)

    nb.save()
    shutil.copy('%s/nb.sobj'%dir, '%s/nb-older-backup.sobj'%dir)
    nb.set_debug(show_debug)
    nb.set_log_server(log_server)
    if warn and address!='localhost' and username==None:
        print "WARNING -- it is *extremely* dangerous to let the server listen"
        print "on an external port without at least setting a username/password!!"
    nb.start(port, address, max_tries, open_viewer, jsmath=jsmath)
    from sage.interfaces.quit import expect_quitall
    expect_quitall(verbose=False)
    from sage.misc.misc import delete_tmpfiles
    delete_tmpfiles()
    if os.path.exists('%s/pid'%dir):
        os.remove('%s/pid'%dir)
    return nb


#######
# Misc

def clean_name(name):
    return ''.join([x if (x.isalnum() or x == '_') else '_' for x in name])

