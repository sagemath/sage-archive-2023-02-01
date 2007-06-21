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
import random
import re
import shutil
import socket
import time

# SAGE libraries
from   sage.structure.sage_object import SageObject, load
from   sage.misc.viewer     import browser
from   sage.misc.misc       import alarm, cancel_alarm, tmp_dir
from   sage.server.misc import print_open_msg

# SAGE Notebook
import css          # style
import js           # javascript
import server       # web server
import worksheet    # individual worksheets (which make up a notebook)
import config       # internal configuration stuff (currently, just keycodes)
import keyboards    # keyboard layouts
import server_conf  # server configuration
import user_conf    # user configuration
import user         # users


SYSTEMS = ['sage', 'axiom', 'gap', 'gp', 'jsmath', 'kash', 'latex', 'lisp', 'macaulay2', 'magma', 'maple', 'mathematica', 'matlab', 'maxima', 'mupad', 'mwrank', 'octave', 'python', 'sage', 'sh', 'singular']

JSMATH = True

vbar = '<span class="vbar"></span>'

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
        self.__conf = server_conf.ServerConfiguration()
        self.default_users()

    ##########################################################
    # Users
    ##########################################################
    def default_users(self):
        self.add_user('pub', '', '', account_type='guest')
        self.add_user('guest', '', '', account_type='guest')
        pw = hex(random.randint(1,2**24))[2:]
        self.add_user('admin', pw, '', account_type='admin')
        print "*"*60
        print "*"*60
        print ""
        print "Setting Admin password to '%s'"%pw
        print ""
        print "*"*60
        print "*"*60

    def users(self):
        try:
            return self.__users
        except AttributeError:
            self.__users = {}
            return self.__users

    def user(self, username):
        if '/' in username:
            raise ValueError
        try:
            return self.__users[username]
        except KeyError:
            raise KeyError, "no user '%s'"%username

    def user_is_admin(self, user):
        return self.user(user).is_admin()

    def user_is_guest(self, username):
        try:
            return self.user(username).is_guest()
        except KeyError:
            return False

    def user_list(self):
        try:
            return list(self.__users.itervalues())
        except AttributeError:
            self.__users = {}
            return []

    def usernames(self):
        U = self.users()
        return U.keys()

    def add_user(self, username, password, email, account_type="user"):
        us = self.users()
        if us.has_key(username):
            raise ValueError, "User '%s' already exists"%username
        U = user.User(username, password, email, account_type)
        us[username] = U

    def change_password(self, username, password):
        self.user(username).set_password(password)

    def del_user(self, username):
        us = self.users()
        if us.has_key(username):
            del us[username]

    def passwords(self):
        """
        Return the username:password dictionary.
        """
        return dict([(user.username(), user.password()) for user in self.user_list()])

    def user_conf(self, username):
        return self.users()[username].conf()

    ##########################################################
    # Publishing worksheets
    ##########################################################
    def _initialize_worksheet(self, src, W):
        """
        src and W are worksheets and W is brand new.
        """
        # Copy over images and other files
        data = src.data_directory()
        if os.path.exists(data):
            shutil.copytree(data, W.directory() + '/data')
        cells = src.cells_directory()
        if os.path.exists(cells):
            shutil.copytree(cells, W.directory() + '/cells')
        W.edit_save(src.edit_text())

    def publish_worksheet(self, worksheet, username):
        """
        Publish the given worksheet.

        This creates a new worksheet in the pub directory with the
        same contents as worksheet.
        """

        for X in self.__worksheets.itervalues():
            if X.owner() == 'pub' and X.worksheet_that_was_published() == worksheet:
                # Update X based on worksheet instead of creating something new
                # 1. delete cells and data directories
                # 2. copy them over
                # 3. update worksheet text
                if os.path.exists(X.data_directory()):
                    shutil.rmtree(X.data_directory())
                if os.path.exists(X.cells_directory()):
                    shutil.rmtree(X.cells_directory())
                self._initialize_worksheet(worksheet, X)
                X.set_worksheet_that_was_published(worksheet)
                X.move_to_archive(username)
                worksheet.set_published_version(X.filename())
                return X

        # Have to create a new worksheet
        W = self.create_new_worksheet(worksheet.name(), 'pub')
        self._initialize_worksheet(worksheet, W)
        W.set_worksheet_that_was_published(worksheet)
        W.move_to_archive(username)
        worksheet.set_published_version(W.filename())
        return W

    ##########################################################
    # Moving, copying, creating, renaming, and listing worksheets
    ##########################################################

    def scratch_worksheet(self):
        try:
            return self.__scratch_worksheet
        except AttributeError:
            W = self.create_new_worksheet('scratch', '_sage_', add_to_list=False)
            self.__scratch_worksheet = W
            return W

    def create_new_worksheet(self, worksheet_name, username, docbrowser=False, add_to_list=True):
        if username!='pub' and self.user_is_guest(username):
            raise ValueError, "guests cannot create new worksheets"

        filename = worksheet.worksheet_filename(worksheet_name, username)
        if self.__worksheets.has_key(filename):
            return self.__worksheets[filename]
        i = 0
        dir = self.worksheet_directory() + '/' + username
        if os.path.exists(dir):
            D = os.listdir(dir)
            D.sort()
            dirname = str(i)
            while dirname in D:
                i += 1
                dirname = str(i)
        else:
            dirname = '0'

        W = worksheet.Worksheet(worksheet_name, dirname, self,
                                system = self.system(username),
                                owner=username,
                                docbrowser = docbrowser)

        if add_to_list:
            self.__worksheets[W.filename()] = W
        return W

    def copy_worksheet(self, ws, owner):
        W = self.create_new_worksheet('default', owner)
        print W.name(), W.filename()
        self._initialize_worksheet(ws, W)
        name = "Copy of %s"%ws.name()
        W.set_name(name)
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

    def migrate_old(self):
        """
        Migrate all old worksheets, i.e., ones with no owner
        to /pub.
        """
        raise NotImplementedError
        for w in self.__worksheets.itervalues():
            if not '/' in w.filename():
                print "Moving worksheet ", w.name()
                w.set_owner('old')
                self.rename_worksheet_filename(w, w.filename())

    ##########################################################
    # Information about the pool of worksheet compute servers
    ##########################################################

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
        if P is None or len(P) == 0:
            return None
        try:
            i = (self.__server_number + 1)%len(P)
        except AttributeError:
            self.__server_number = 0
            i = 0
        return P[i]

    ##########################################################
    # The default math software system for new worksheets for
    # a given user or the whole notebook (if username is None).
    ##########################################################

    # TODO -- only implemented for the notebook right now
    def system(self, username=None):
        return self.user(username).conf()['default_system']

    def set_system(self, system):
        self.__system = system

    ##########################################################
    # The default color scheme for the notebook.
    ##########################################################
    def color(self):
        try:
            return self.__color
        except AttributeError:
            self.__color = 'default'
            return self.__color

    def set_color(self,color):
        self.__color = color

    ##########################################################
    # The directory the notebook runs in.
    ##########################################################
    def set_directory(self, dir):
        if dir == self.__dir:
            return
        self.__dir = dir
        self.__filename = '%s/nb.sobj'%dir
        self.__worksheet_dir = '%s/worksheets'%dir
        self.__object_dir = '%s/objects'%dir
        for W in self.__worksheets.itervalues():
            W.set_notebook(self)

    ##########################################################
    # The notebook history.
    ##########################################################
    def user_history(self, username):
        U = self.user(username)
        try:
            return U.history
        except AttributeError:
            U.history = []
            return U.history

    def create_new_worksheet_from_history(self, name, username, maxlen=None):
        W = self.create_new_worksheet(name, username)
        W.edit_save('Log Worksheet\n' + self.user_history_text(username, maxlen=None))
        return W

    def user_history_text(self, username, maxlen=None):
        H = self.user_history(username)
        if maxlen:
            H = H[-maxlen:]
        return '\n\n'.join([L.strip() for L in H])

    def user_history_html(self, username):
        t = self.user_history_text(username)
        t = t.replace('<','&lt;')
        s = """
        <html>
        <head>
           <link rel=stylesheet href="/css/main.css">
           <title>SAGE: History for %s</title>
        </head>
        <body>
        %s
        <pre>
        %s
        </pre>
        <hr class="usercontrol">
        <a title="Click here to turn the above into a SAGE worksheet" href="/live_history">Create a new SAGE worksheet version of the last 100 commands in the above log.</a>
        <a name="bottom"></a>
        <script type="text/javascript"> window.location="#bottom"</script>
        </body>
        </html>
        """%(username, self.html_worksheet_list_top(username, actions=False), t)
        return s


    def add_to_user_history(self, entry, username):
        H = self.user_history(username)
        H.append(entry)
        maxlen = self.user_conf(username)['max_history_length']
        while len(H) > maxlen:
            del H[0]

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

    def max_history_length(self):
        try:
            return self.conf()['max_history_length']
        except KeyError:
            return MAX_HISTORY_LENGTH

    def history_html(self):
        t = self.history_text()
        t = t.replace('<','&lt;')
        s = '<head>\n'
        s += '<title>Command History</title>\n'
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

    ##########################################################
    # Importing and exporting worksheets to files
    ##########################################################
    def export_worksheet(self, worksheet_filename, output_filename):
        W = self.get_worksheet_with_filename(worksheet_filename)
        W.save()
        path = W.filename_without_owner()
        cmd = 'cd "%s/%s/" && tar -jcf "%s" "%s"'%(
            self.__worksheet_dir, W.owner(),
            os.path.abspath(output_filename), path)
        e = os.system(cmd)
        if e:
            print "Failed to execute command to export worksheet:\n'%s'"%cmd

    def new_worksheet_with_title_from_text(self, text, owner):
        name, _ = worksheet.extract_name(text)
        W = self.create_new_worksheet(name, owner)
        return W

    def change_worksheet_key(self, old_key, new_key):
        ws = self.__worksheets
        W = ws[old_key]
        ws[new_key] = W
        del ws[old_key]

    def import_worksheet(self, filename, owner):
        """
        Upload the worksheet with name filename and make it have the
        given owner.
        """
        if not os.path.exists(filename):
            raise ValueError, "no file %s"%filename

        # Decompress the worksheet to a temporary directory.
        tmp = tmp_dir()
        cmd = 'cd "%s"; tar -jxf "%s"'%(tmp, os.path.abspath(filename))
        print cmd
        e = os.system(cmd)
        if e:
            raise ValueError, "Error decompressing saved worksheet."

        # Find the worksheet text representation and load it into memory.
        try:
            D = os.listdir(tmp)[0]
        except IndexError:
            raise ValueError, "invalid worksheet"
        text_filename = '%s/%s/worksheet.txt'%(tmp,D)
        worksheet_txt = open(text_filename).read()
        worksheet = self.new_worksheet_with_title_from_text(worksheet_txt, owner)
        worksheet.set_owner(owner)
        name = worksheet.filename_without_owner()

        # Change the filename of the worksheet, if necessary
        names = [w.filename_without_owner() for w in self.get_worksheets_with_owner(owner)]
        if name in names:
            name = 0
            while str(name) in names:
                name += 1
            name = str(name)
            worksheet.set_filename_without_owner(name)

        # Change the display name of the worksheet if necessary
        name = worksheet.name()
        display_names = [w.name() for w in self.get_worksheets_with_owner(owner)]
        if name in display_names:
            j = name.rfind('(')
            if j != -1:
                name = name[:j].rstrip()
            i = 2
            while name + " (%s)"%i in display_names:
                i += 1
            name = name + " (%s)"%i
            worksheet.set_name(name)


        # Put the worksheet files in the target directory.
        S = self.__worksheet_dir
        target = '%s/%s'%(os.path.abspath(S), worksheet.filename())
        if not os.path.exists(target):
            os.makedirs(target)
        cmd = 'rm -rf "%s"/*; mv "%s/%s/"* "%s/"'%(target, tmp, D, target)
        print cmd
        if os.system(cmd):
            raise ValueError, "Error moving over files when loading worksheet."

        worksheet.edit_save(worksheet_txt)

        shutil.rmtree(tmp)

        return worksheet


    ##########################################################
    # Importing and exporting worksheets to a plain text format
    ##########################################################

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

    ##########################################################
    # Directories for worksheets, etc.
    ##########################################################
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

    def worksheet_directory(self):
        return self.__worksheet_dir

    def __makedirs(self):
        if not os.path.exists(self.__dir):
            os.makedirs(self.__dir)
        if not os.path.exists(self.__worksheet_dir):
            os.makedirs(self.__worksheet_dir)
        if not os.path.exists(self.__object_dir):
            os.makedirs(self.__object_dir)

    ##########################################################
    # Server configuration
    ##########################################################
    def conf(self):
        try:
            return self.__conf
        except AttributeError:
            C = server_conf.ServerConfiguration()
            self.__conf = C
            return C


    def set_debug(self,show_debug):
        self.__show_debug = show_debug

    def number_of_backups(self):
        return self.conf()['number_of_backups']

    def backup_directory(self):
        try:
            D = self.__backup_dir
        except AttributeError:
            D = self.__dir + "/backups/"
            self.__backup_dir = D
        if not os.path.exists(D):
            os.makedirs(D)
        return D


    ##########################################################
    # The object store for the notebook.
    ##########################################################
    # Todo: like with worksheets, objects should belong to
    # users, some should be published, rateable, etc.
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

    ##########################################################
    # Computing control
    ##########################################################
    def set_not_computing(self):
        # unpickled, no worksheets will think they are
        # being computed, since they clearly aren't (since
        # the server just started).
        for W in self.__worksheets.values():
            W.set_not_computing()

    def quit(self):
        for W in self.__worksheets.itervalues():
            W.quit()

    def quit_idle_worksheet_processes(self):
        timeout = self.conf()['idle_timeout']
        for W in self.__worksheets.itervalues():
            if W.compute_process_has_been_started():
                W.quit_if_idle(timeout)


    ##########################################################
    # Worksheet HTML generation
    ##########################################################
    def list_window_javascript(self, worksheet_filenames):
        s = """
           <script type="text/javascript" src="/javascript/main.js"></script>
           <script type="text/javascript">
           var worksheet_filenames = %s;
           </script>
        """%(worksheet_filenames)

        return s

    def worksheet_html(self, filename, do_print=False):
        W = self.get_worksheet_with_filename(filename)
        s = '<head>\n'
        s += '<title>SAGE Worksheet: %s</title>\n'%W.name()
        if do_print:
            s += '<script type="text/javascript" src="/javascript/jsmath/jsMath.js"></script>\n'
        s += '<script type="text/javascript" src="/javascript/main.js"></script>\n'
        s += '<link rel=stylesheet href="/css/main.css">\n'
        s += '</head>\n'
        s += '<body>\n'
        #if do_print:
        #    s += '<h1><a href=".">SAGE Worksheet: %s</a></h1>'%W.name()
        s += W.html(include_title=False, do_print=do_print)
        if do_print:
            s += '<script type="text/javascript">jsMath.Process();</script>\n'
        s += '\n</body>\n'
        return s


    def html_worksheet_list_public(self, username,
                                  sort='last_edited', reverse=False, search=None):
        W = [x for x in self.__worksheets.itervalues() if x.is_published() and not x.is_trashed(user)]

        if search:
            W = [x for x in W if x.satisfies_search(search)]

        sort_worksheet_list(W, sort, reverse)  # changed W in place
        worksheet_filenames = [x.filename() for x in W]

        top  = self.html_worksheet_list_top(username, False, pub=True, search=search)
        list = self.html_worksheet_list(W, username, False, sort=sort, reverse=reverse, typ="all", pub=True)

        s = """
        <html>
           <link rel=stylesheet href="/css/main.css">
           <title>SAGE: Published Worksheets</title>
           %s
        <body>
        %s
        %s
        </body>
        </html>
        """%(self.list_window_javascript(worksheet_filenames), top, list)

        return s

    def html_worksheet_list_for_user(self, user,
                                     typ="active",
                                     sort='last_edited',
                                     reverse=False,
                                     search=None):

        X = [x for x in self.get_worksheets_with_viewer(user) if not x.docbrowser()]
        if typ == "trash":
            worksheet_heading = "Trash"
            W = [x for x in X if x.is_trashed(user)]
        elif typ == "active":
            worksheet_heading = "Active Worksheets"
            W = [x for x in X if x.is_active(user)]
        else: # typ must be archived or "all"
            worksheet_heading = "Archived and Active"
            W = [x for x in X if not x.is_trashed(user)]
        if search:
            W = [x for x in W if x.satisfies_search(search)]

        sort_worksheet_list(W, sort, reverse)  # changed W in place

        top = self.html_worksheet_list_top(user, typ=typ, search=search)
        list = self.html_worksheet_list(W, user, worksheet_heading, sort=sort, reverse=reverse, typ=typ)
        worksheet_filenames = [x.filename() for x in W]

        s = """
        <html>
           <link rel=stylesheet href="/css/main.css">
           <title>SAGE: Worksheet List</title>
           %s
        <body>
        %s
        %s
        </body>
        </html>
        """%(self.list_window_javascript(worksheet_filenames), top, list)

        return s

    def html_topbar(self, user, pub=False):
        s = ''
        entries = []
        if self.user_is_guest(user):
            entries.append(('/', 'Log in', 'Please log in to the SAGE notebook'))
        else:
            entries.append(('/home/%s'%user, 'Home', 'Back to your personal worksheet list'))
            #entries.append(('/settings', 'Settings', 'Change user settings'))   # TODO -- settings
            entries.append(('/help', 'Help', 'Documentation'))

        ## TODO -- settings
        #if self.user(user).is_admin():
        #    entries.insert(1, ('/notebook_settings', 'Server', 'Change general SAGE notebook server configuration'))
        if not pub:
            entries.insert(1, ('history_window()', 'Log', 'View a log of recent computations'))
        if not self.user_is_guest(user):
            entries.append(('/logout', 'Sign out', 'Logout of the SAGE notebook'))

        s += self.html_user_control(user, entries)
        s += self.html_banner()
        s += '<hr class="usercontrol">'
        return s

    def html_worksheet_list_top(self, user, actions=True, typ='active', pub=False, search=None):
        s = self.html_topbar(user, pub)
        s += self.html_new_or_upload()
        s += self.html_search(search, typ)
        s += '<br>'
        s += '<hr class="usercontrol">'
        if actions:
            s += self.html_worksheet_actions(user, typ=typ)
        return s


    def html_user_control(self, user, entries):
        s = ''
        s += '<div class="flush-right">'
        s += '<span class="username">%s</span>'%user
        for href, name, title in entries:
            if '(' in href:
                action = 'onClick="%s"'%href
            else:
                action = 'href="%s"'%href
            x = '<a title="%s" class="usercontrol" %s>%s</a>\n'%(title, action, name)
            s += vbar + x
        s += '</div>'
        return s

    def html_banner(self):
        s = """
        <span class="banner">
        <a class="banner" href="http://www.sagemath.org"><img align="top" src="/images/sagelogo.png" alt="SAGE"> Mathematics Software</a> <span class="ping" id="ping">Unable to contact SAGE server!</span>
        </span>
        """
        return s

    def html_search(self, search, typ):
        s = """
        <span class="flush-right">
        <input id="search_worksheets" size=20 onkeypress="return entsub(event);" value="%s"></input>
        <button class="add_new_worksheet_menu" onClick="search_worksheets('%s');">Search Worksheets</button>
        &nbsp;&nbsp;&nbsp;
        </span>
        """%('' if search is None else search, typ)
        return s

    def html_new_or_upload(self):
        s = """
        <a class="boldusercontrol" href="/new_worksheet">New Worksheet</a>\n
        <a class="boldusercontrol" href="/upload">Upload</a>\n
        """
        return s

    def html_worksheet_actions(self, user, typ):
##
##          <option onClick="save_worksheets('sws');" title="Save the selected worksheets to disk">Save ...</option>
##          <option onClick="save_worksheets('html');" title="Save the selected worksheets as a single HTML web page">Save as HTML (zipped) ... </option>
##          <option onClick="save_worksheets('latex');" title="Save the selected worksheets as a single LaTeX document">Save as LaTeX (zipped) ... </option>
##          <option onClick="save_worksheets('pdf');" title="Save the selected worksheets as a single PDF document">Save as PDF...</option>
##          <option onClick="save_worksheets('txt');" title="Save the selected worksheets to a single text file">Save as Text...</option>
##         s = """
##          <select class="worksheet_list">
##          <option onClick="archive_button();" title="Archive selected worksheets so they do not appear in the default worksheet list">Archive</option>
##          <option onClick="make_active_button();" title="Unarchive this worksheet so it appears in the default worksheet list">Unarchive</option>
##          <option onClick="uncollaborate_me();" title="Remove myself from collaboration or viewing of this worksheet">Un-collaborate me</option>
##         </select>
##         """
        s = ''

        if not self.user_is_guest(user):
            if typ == 'archive':
                s += '<button onClick="make_active_button();" title="Unarchive selected worksheets so it appears in the default worksheet list">Unarchive</button>'
            else:
                s += '<button onClick="archive_button();" title="Archive selected worksheets so they do not appear in the default worksheet list">Archive</button>'

            if typ != 'trash':
                s += '&nbsp;&nbsp;<button onClick="delete_button();" title="Move the selected worksheets to the trash">Delete</button>'
            else:
                s += '&nbsp;&nbsp;<button onClick="make_active_button();" title="Move the selected worksheets out of the trash">Undelete</button>'

            s += '<span class="flush-right">'
            s += '<a class="control" href="/pub" title="Browse everyone\'s published worksheets">Published</a>'
            s += '&nbsp;'*10
            s += '&nbsp;<a class="usercontrol" href=".">Active</a>'
            s += '&nbsp;<a class="usercontrol" href=".?typ=archive">Archive</a>'
            s += '&nbsp;<a class="usercontrol" href=".?typ=trash">Trash</a>&nbsp;&nbsp;'
            s += '</span>'
        return s


    def html_worksheet_list(self, worksheets, user, worksheet_heading, sort, reverse, typ, pub=False):
        s = ''

        s = '<br><br>'
        s += '<table width=100% border=0 cellspacing=0 cellpadding=0>'
        s += '<tr class="greybox"><td colspan=4><div class="thinspace"></div></td></tr>'
        s += '<tr  class="greybox">'

        if not pub:
            s += '<td>&nbsp;<input id="controlbox" onClick="set_worksheet_list_checks();" class="entry" type=checkbox></td>'
        else:
            s += '<td>&nbsp;&nbsp;<a class="listcontrol" href=".?sort=rating">Rating</a></td>'
            worksheet_heading = "Published Worksheets"

        s += '<td><a class="listcontrol" href=".?typ=%s&sort=name%s">%s</a> </td>'%(typ,
            '' if sort != 'name' or reverse else '&reverse=True', worksheet_heading)
        s += '<td><a class="listcontrol" href=".?typ=%s&sort=owner%s">Owner / Collaborators</a> </td>'%(typ,
            '' if sort != 'owner' or reverse else '&reverse=True')
        s += '<td><a class="listcontrol" href=".?typ=%s&%s">Last Edited</a> </td>'%(typ,
            '' if sort != 'last_edited' or reverse else 'reverse=True')
        s += '</tr>'
        s += '<tr class="greybox"><td colspan=4><div class="thinspace"></div></td></tr>'

        v = []
        for w in worksheets:
            k = '<tr>'
            k += '<td class="entry">%s</td>'%self.html_check_col(w, user, pub)
            if w.is_active(user):
                k += '<td class="worksheet_link">%s</td>'%self.html_worksheet_link(w, pub)
            else:
                k += '<td class="archived_worksheet_link">%s</td>'%self.html_worksheet_link(w, pub)
            k += '<td class="owner_collab">%s</td>'%self.html_owner_collab_view(w, user, typ)
            k += '<td class="last_edited">%s</td>'%self.html_last_edited(w, user)
            k += '</tr>'
            k += '<tr class="thingreybox"><td colspan=4><div class="ultrathinspace"></div></td></tr>'
            v.append(k)

        s += ''.join(v)
        s += '</table>'

        return s

    def html_check_col(self, worksheet, user, pub):
        def doc_options(name):
            if pub:
                rating = worksheet.rating()
                if rating == 0:
                    r = "----"
                else:
                    r = "%.1f"%rating
                if not self.user_is_guest(user) and not worksheet.is_rater(user):
                    r = '<i>%s</i>'%r
                return r

            return """
            <select class="worksheet_edit">
            <option onClick="list_edit_worksheet('%s');" title="Open this worksheet and edit it">Edit</option>
            <option onClick="list_copy_worksheet('%s');" title="Copy this worksheet">Copy Worksheet</option>
            <option onClick="list_share_worksheet('%s');" title="Share this worksheet with others">Collaborate</option>
            <option onClick="list_publish_worksheet('%s');" title="Publish this worksheet on the internet">Publish</option>
            <option onClick="list_revisions_of_worksheet('%s');" title="See all revisions of this worksheet">Revisions</option>
            <option onClick="list_preview_worksheet('%s');" title="Preview this worksheet">Preview</option>
            </select>
            """%(name, name,name,name,name,name)

        k = ''
        if not pub:
            k += '<input type=checkbox unchecked id="%s">'%worksheet.filename()
        k += '&nbsp;'*4
        k += doc_options(worksheet.filename())
        k += '&nbsp;'*4
        return k

    def html_worksheet_link(self, worksheet, pub):
        name = worksheet.truncated_name(35)
        if not pub and worksheet.is_published():
            name += ' (published)'
        return '<a title="%s" id="name/%s" class="worksheetname" href="/home/%s">%s</a>\n'%(
            worksheet.name(), worksheet.filename(), worksheet.filename(), name)

    def html_owner_collab_view(self, worksheet, user, typ):
        v = []

        owner = worksheet.owner()
        pub = False
        if owner == 'pub':
            pub = True
            owner = worksheet.worksheet_that_was_published().owner()
        elif owner == user:
            owner = "Me"

        v.append(owner)

        collab = worksheet.collaborators()

        if not pub and typ != 'trash' and (owner == "Me" or self.user(user).is_admin()):
            if len(collab) <= 1:
                share = '<a class="share" href="/home/%s/share">Share now</a>'%(worksheet.filename())
            else:
                collaborators = ', '.join([x for x in collab if x != user])
                if len(collaborators) > 21:
                    collaborators = collaborators[:21] + '...'
                v.append(collaborators)
                share = '<a class="share" href="/home/%s/share">Add or Delete</a>'%(worksheet.filename())
        else:
            share = ''
        if worksheet.has_published_version():
            pub_ver = worksheet.published_version().filename()
            share += ' <a href="/home/%s">(published)'%pub_ver

        viewers = worksheet.viewers()
        if len(viewers) > 0:
            viewers = '<i>' + ', '.join(viewers) + '</i>'
            v.append(viewers)

        s = ' / '.join(v) + ' ' + share

        return s

    def html_last_edited(self, worksheet, user):
        s = worksheet.html_time_since_last_edited()
        who = worksheet.last_to_edit()
        if who == user:
            who = 'Me'
        return s + ' ago by ' + who

    ##########################################################
    # Revision history for a worksheet
    ##########################################################
    def html_worksheet_revision_list(self, username, worksheet):
        head, body = self.html_worksheet_page_template(worksheet, username, "Revision History")
        data = worksheet.snapshot_data()  # pairs ('how long ago', key)
        rows = []
        i = 0
        for desc, key in data:
            rows.append('<tr><td></td><td><a href="revisions?rev=%s">Revision %s</a></td><td><span class="revs">%s</span></td></tr>'%
                        (key, i, desc))
            i += 1

        rows = list(reversed(rows))
        rows = '\n'.join(rows)
        body += """
        <hr class="usercontrol">
<table width=100%%>
<tr><td width=1%%></td><td width=20%%><b>Revision</b></td> <td width=20%%><b>Last Edited</b></td><td width=30%%></td>
%s
</table>
"""%rows

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)


    def html_specific_revision(self, username, ws, rev):
        t = time.time() - float(rev[:-4])
        when = worksheet.convert_seconds_to_meaningful_time_span(t)
        head, body = self.html_worksheet_page_template(ws, username, "Revision from %s ago"%when)

        filename = ws.get_snapshot_text_filename(rev)
        txt = open(filename).read()
        W = self.scratch_worksheet()
        W.delete_cells_directory()
        W.edit_save(txt)
        html = W.html_worksheet_body(do_print=True, publish=True)

        actions = """
        <a class="listcontrol" href="revisions?rev=%s&action=publish">Publish this one</a>&nbsp;&nbsp;
        <a class="listcontrol" href="revisions?rev=%s&action=revert">Revert to this one</a>  (note that images are note recorded)
        """%(rev, rev)

        s = """
        %s
        <hr class="usercontrol">
<table width=100%%>
%s
        <hr class="usercontrol">
%s
</table>
"""%(actions, html, actions)
        body += s

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_worksheet_page_template(self, worksheet, username, title):
        head = self._html_head(worksheet_filename=worksheet.filename(), username=username)
        head += '<script  type="text/javascript">worksheet_filename="%s"; worksheet_name="%s"; server_ping_while_alive(); </script>'%(worksheet.filename(), worksheet.name())
        body = self._html_body(worksheet.filename(), top_only=True, username=username)
        body += self.html_worksheet_topbar(worksheet)
        body += '<hr class="usercontrol">'
        body += '<span class="sharebar">%s</span>'%title
        body += '<br>'*3
        return head, body


    def html_share(self, worksheet, username):
        head, body = self.html_worksheet_page_template(worksheet, username, "Share this document")

        body += 'This SAGE Worksheet is currently shared with the people listed in the box below.<br>'
        body += 'You may add or remove collaborators (separate user names by commas).<br><br>'

        collabs = ', '.join(worksheet.collaborators())
        body += '<form width=70% method="post" action="invite_collab">\n'
        body += '<textarea name="collaborators" rows=5 cols=70 class="edit" id="collaborators">%s</textarea><br><br>'%collabs
        body += '<input type="submit" title="Give access to your worksheet to the above collaborators" value="Invite Collaborators">'
        body += '</form>'

        body += '<br>'*2
        body += '<hr class="usercontrol">'
        body += '<span class="username">SAGE Users:</span>'
        U = self.users()
        K = [x for x, u in U.iteritems() if not u.is_guest() and u.username() != username]
        def mycmp(x,y):
            return cmp(x.lower(), y.lower())
        K.sort(mycmp)
        body += '<span class="users">%s</span>'%(', '.join(K))


        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)



    ##########################################################
    # Accessing all worksheets with certain properties.
    ##########################################################
    def get_all_worksheets(self):
        return list(self.__worksheets.itervalues())

    def get_worksheets_with_collaborator(self, user):
        if user == 'admin': return self.get_all_worksheets()
        return [w for w in self.__worksheets.itervalues() if w.user_is_collaborator(user)]

    def get_worksheet_names_with_collaborator(self, user):
        if user == 'admin': return [W.name() for W in self.get_all_worksheets()]
        return [W.name() for W in self.get_worksheets_with_collaborator(user)]

    def get_worksheets_with_viewer(self, user):
        if user == 'admin': return self.get_all_worksheets()
        return [w for w in self.__worksheets.itervalues() if w.user_is_viewer(user)]

    def get_worksheets_with_owner(self, owner):
        return [w for w in self.__worksheets.itervalues() if w.owner() == owner]

    def get_worksheets_with_owner_that_are_viewable_by_user(self, owner, user):
        return [w for w in self.get_worksheets_with_owner(owner) if w.user_is_viewer(user)]

    def get_worksheet_names_with_viewer(self, user):
        if user == 'admin': return [W.name() for W in self.get_all_worksheets()]
        return [W.name() for W in self.get_worksheets_with_viewer(user)]

    def get_worksheet_with_name(self, name):
        for W in self.__worksheets.itervalues():
            if W.name() == name:
                return W
        raise KeyError, "No worksheet with name '%s'"%name

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
        raise KeyError, "No worksheet with filename '%s'"%filename

    ###########################################################
    # Saving the whole notebook
    ###########################################################

    def save(self, filename=None):
        #print "-"*70

        if filename is None:
            F = os.path.abspath(self.__filename)
            backup_dir = self.backup_directory()
            backup = backup_dir + '/nb-backup-'
            for i in range(self.number_of_backups()-1,0,-1):
                a = padzeros(i-1); b = padzeros(i)
                try:
                    shutil.move(backup + '%s.sobj'%a, backup + '%s.sobj'%b)
                except IOError, msg:
                    pass
            a = '%s.sobj'%padzeros(0)
            try:
                shutil.copy(F, backup + a)
            except Exception, msg:
                pass
            F = os.path.abspath(self.__filename)
        else:
            F = os.path.abspath(filename)

        #print "Saving notebook to '%s'..."%F
        D, _ = os.path.split(F)
        if not os.path.exists(D):
            os.makedirs(D)
        SageObject.save(self, F, compress=False)
        #print "Press control-C to stop the notebook server."
        #print "-"*70

    def delete_doc_browser_worksheets(self):
        names = self.worksheet_names()
        for n in self.__worksheets.keys():
            if n.startswith('doc_browser'):
                self.delete_worksheet(n)

    ###########################################################
    # HTML -- generate most html related to the whole notebook page
    ###########################################################
    def html_slide_controls(self):
        return """
          <div class="hidden" id="slide_controls">
            <div class="slideshow_control">
             <a class="slide_arrow" onClick="slide_next()">&gt;</a>
              <a class="slide_arrow" onClick="slide_last()">&gt;&gt;</a> %s
              <a class="cell_mode" onClick="cell_mode()">Exit</a>
            </div>
            <div class="slideshow_progress" id="slideshow_progress" onClick="slide_next()">
              <div class="slideshow_progress_bar" id="slideshow_progress_bar">&nbsp;</div>
              <div class="slideshow_progress_text" id="slideshow_progress_text">&nbsp;</div>
            </div>
            <div class="slideshow_control">
              <a class="slide_arrow" onClick="slide_first()">&lt;&lt;</a>
              <a class="slide_arrow" onClick="slide_prev()">&lt;</a>
            </div>
          </div>
          """%vbar

    def html_debug_window(self):
        return """
        <div class='debug_window'>
        <div class='debug_output'><pre id='debug_output'></pre></div>
        <textarea rows=5 id='debug_input' class='debug_input'
         onKeyPress='return debug_keypress(event);'
         onFocus='debug_focus();' onBlur='debug_blur();'></textarea>
        </div>"""


    def _html_head(self, worksheet_filename, username):
        if worksheet_filename is not None:
            worksheet = self.get_worksheet_with_filename(worksheet_filename)
            head = '\n<title>SAGE: %s</title>'%(worksheet.name())
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

    def html_worksheet_topbar(self, worksheet):
        body = ''
        body += worksheet.html_title()
        body += worksheet.html_save_discard_buttons() + '<br>'

        body += '<hr class="greybar">'
        menu = worksheet.html_menu()
        share = '<span class=thin-right>%s</span>'%worksheet.html_share_publish_buttons()

        body += '<table><tr><td>%s</td><td>%s</td></tr></table>'%(menu, share)

        body += self.html_slide_controls()
        return body

    def _html_body(self, worksheet_filename, show_debug=False, username='', top_only=False):
        worksheet = self.get_worksheet_with_filename(worksheet_filename)
        worksheet_html = worksheet.html()

        body = ''

        if worksheet.is_published() or self.user_is_guest(username):
            original_worksheet = worksheet.worksheet_that_was_published()
            body += '<h1 align=center>%s</h1>'%original_worksheet.name()
            body += '<h2 align=center>%s</h2>'%worksheet.html_time_last_edited()
            body += worksheet_html
            body += '<hr class="usercontrol">'
            if original_worksheet.user_is_collaborator(username) or original_worksheet.is_owner(username):
                s = "Edit this worksheet."
                url = 'edit_published_page'
            elif self.user_is_guest(username):
                s = 'Log in if you would like to edit a copy of this worksheet.'
                url = '/'
            else:
                s = 'Edit a copy of this worksheet.'
                url = 'edit_published_page'
            body += '<a class="usercontrol" href="%s">%s</a>'%(url, s)
            body += '&nbsp;'*10

            r = worksheet.rating()
            if r == 0:
                rating = 'This page has not been rated.'
            else:
                rating = '<a class="usercontrol" href="rating_info">This page is rated %.1f.</a>'%r
            if not self.user_is_guest(username) and not worksheet.is_rater(username) \
                   and not worksheet.is_publisher(username):
                rating += '&nbsp;&nbsp; Rate it: '
                rating += '  '.join(['<a class="usercontrol" href="rate%s">&nbsp;%s&nbsp;</a>'%(i,i) for
                                   i in range(1,5)])
            body += '<span class="ratingmsg">%s</span>'%rating

            body += '<hr class="usercontrol">'
            body += '<br>'
            body += '<span class="pubmsg">This document was published using <a href="/">SAGE</a>.'
            body += '  Browser <a href="/pub/">other published documents.</a></span>'

        else:

            entries = [('/', 'Home', 'Back to your personal worksheet list'),
                       ('history_window()', 'Log', 'View a log of recent computations'),
                       #('settings', 'Settings', 'Worksheet settings'),  # TODO -- settings
                       ('/help', 'Help', 'Documentation')]

            if not self.user_is_guest(username):
                entries.append(('/logout', 'Sign out', 'Logout of the SAGE notebook'))

            body += self.html_user_control(username, entries)
            body += self.html_banner()+ '<br><br>'


            if top_only:
                return body

            if worksheet_filename:
                body += self.html_worksheet_topbar(worksheet)

            if self.__show_debug or show_debug:
                body += self.html_debug_window()


            body += '<div class="worksheet" id="worksheet">%s</div>'%worksheet_html

        # The blank space given by '<br>'*15  is needed so the input doesn't get
        # stuck at the bottom of the screen. This could be replaced by a region
        # such that clicking on it creates a new cell at the bottom of the worksheet.
        body += '<br>'*15
        body += '\n</div>\n'
        endpanespan = '</td></tr></table></span>\n'


        if worksheet is None:
             return body + endpanespan

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

        return body

    def edit_window(self, worksheet, username):
        """
        Return a window for editing worksheet.

        INPUT:
            worksheet -- a worksheet
        """
        head, body = self.html_worksheet_page_template(worksheet, username, 'Edit Text &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<input type="submit" value="Save Changes" name="button_save"> <input type="submit" value="Cancel" name="button_cancel">')

        t = worksheet.edit_text()
        t = t.replace('<','&lt;')
        body = '        <form method="post" action="save" enctype="multipart/form-data">' + body
        body += """
        <textarea class="plaintextedit" id="cell_intext" name="textfield">%s</textarea>
        </form>
        """%t

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_notebook_help_window(self, username):
        top = self._html_head(None, username) + self.html_topbar(username)

        from tutorial import notebook_help
        s = """
        <html>
        <title>SAGE Documentation</title>

        <body>
        """ + top + \
        """
        <style>

        div.help_window {
            font-family: sans-serif;
            background-color:white;
            border: 3px solid #3d86d0;
            top: 10ex;
            bottom:10%;
            left:25%;
            right:15%;
            padding:2ex;
            width:80%;
        }


        table.help_window {
            background-color:white;
            width:95%;
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

        <center>
        <br>
        <a class="control" href="/doc/live/tut/index.html">Tutorial</a>
        &nbsp;&nbsp;
        <a class="control" href="/doc/live/ref/index.html">Reference Manual</a>
        &nbsp;&nbsp;
        <a class="control" href="/doc/live/prog/index.html">Programming Guide</a>
        &nbsp;&nbsp;
        <a class="control" href="/doc/live/const/const.html">Constructions</a>
        <br><br>
        <hr class="usercontrol">
        <br>

        <div class="help_window">
        <h2>How to use the SAGE Notebook</h2>

        <br><br>

        A <i>worksheet</i> is an ordered list of SAGE calculations with output. <br>
        A <i>session</i> is a worksheet and a set of variables in some state.<br>
        The <i>SAGE notebook</i> is a collection of worksheets, saved objects, and user information.<br>
        <br>
        <br>
        To get started with SAGE, <a href="/doc/live/tut/tut.html">work through the tutorial</a> (if
        you have trouble with it, view the <a href="/doc/static/tut/tut.html">static version</a>).
        <br><br>

        <table class="help_window">
        """

        for x, y in notebook_help:
            s += '<tr><td class="help_window_cmd">%s</td><td class="help_window_how">%s</td></tr>\n'%(x,y)
        s += '</table></div>'

        s +="""
        <br>        <br>
        The SAGE Notebook was written by Tom Boothby, Timothy Clemens, Alex Clemesha, Bobby Moretti, Yi Qiang, Dorian Ramier, and William Stein.
        </center>
        </body>
        </html>
        """
        return s

    def upload_window(self):
        return """
          <html>
            <head>
              <title>Upload File</title>
              <style>%s</style>
            </head>
            <body>
              <div class="upload_worksheet_menu" id="upload_worksheet_menu">
              <h1><font size=+3 color="darkred">SAGE</font>&nbsp;&nbsp;&nbsp;&nbsp;<font size=+1>Upload your Worksheet</font></h1>
              <hr>
              <form method="POST" action="upload_worksheet"
                    name="upload" enctype="multipart/form-data">
              <table><tr>
              <td>
              <b>Browse your computer to select a worksheet file to upload:</b><br>
              <input class="upload_worksheet_menu" size="50" type="file" name="fileField" id="upload_worksheet_filename"></input><br><br>
              <b>Or enter the url of a worksheet file on the web:</b><br>

              <input class="upload_worksheet_menu" size="50" type="text" name="urlField" id="upload_worksheet_url"></input></br>
              <br><br>
              <b>What do you want to call it? (if different than the original name)</b><br>
              <input class="upload_worksheet_menu" size="50" type="text" name="nameField" id="upload_worksheet_name"></input></br>
              </td>
              </tr>
              <tr>
              <td><br><input type="button" class="upload_worksheet_menu" value="Upload Worksheet" onClick="form.submit();"></td>
              </tr>
              </form><br>
              </div>
            </body>
          </html>
         """%css.css(self.color())

    def html(self, worksheet_filename=None, username='guest', show_debug=False, admin=False):
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
            head += '<script  type="text/javascript">worksheet_filename="%s"; worksheet_name="%s"; server_ping_while_alive(); </script>'%(worksheet_filename, W.name())

        return """
        <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
        <html>
        <head>%s</head>
        <body>%s</body>
        <script type="text/javascript">jsmath_init()</script>
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


    ####################################################################
    # Configuration html.
    # In each case the settings html is a form that when submitted
    # pulls up another web page and sets the corresponding options.
    ####################################################################
    def html_system_select_form_element(self, ws):
        system = ws.system()
        options = ''
        for S in SYSTEMS:
            if S == system:
                selected = "selected=1"
            else:
                selected = ''
            options += '<option %s onClick="system_select(\'%s\');">%s</option>\n'%(selected, S,S)
        s = """<select class="worksheet">
            %s
            </select>"""%options
        return s

    def html_worksheet_settings(self, ws, username):
        head, body = self.html_worksheet_page_template(ws, username, 'Worksheet Settings &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<button name="button_save">Save Settings</button>  <input type="submit" value="Cancel" name="button_cancel"/>')

        body = '<form width=70%% method="post" action="input_settings"  enctype="multipart/form-data">' + body
        body += '</form>'

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_settings(self):
        s = """
        <h1>Settings</h1>
        """
        return s

    def html_user_settings(self, username):
        s = self.html_settings()
        return s

    def html_notebook_settings(self):
        s = self.html_settings()
        return s

    def html_doc(self, username):
        top = self._html_head(None, username) + self.html_topbar(username)
        body = """
        <br>
        <div class="docidx">
        <h1>SAGE Documentation</h1>
        <br>
        <hr class="usercontrol">
        <br><br>
        <font size=+2>
        <a href="/doc/live/">Documentation</a><br><br>
        <a href="/help/">SAGE Notebook Howto</a><br><br>
        <br><br>
        <br>
        <hr class="usercontrol">
        </font>
        </div>
        """
        #(<a href="/doc/static/">static</a>)

        s = """
        <html>
        %s
        <body>
        %s
        </body>
        </html>
        """%(top, body)

        return s

####################################################################

def load_notebook(dir, address=None, port=None, secure=None):
    """
    Load the notebook from the given directory, or create one in that directory
    if one isn't already there.

    INPUT:
        dir -- a string that defines a directory name
        address -- the address that the notebook server will listen on
        port -- the port the server listens on
        secure -- whether or not the notebook is secure
    """
    sobj = '%s/nb.sobj'%dir
    nb = None
    if os.path.exists(dir):
        try:
            nb = load(sobj, compress=False)
        except:
            backup = '%s/backups/'%dir
            if os.path.exists(backup):
                print "****************************************************************"
                print "  * * * WARNING   * * * WARNING   * * * WARNING   * * * "
                print "WARNING -- failed to load notebook object. Trying backup files."
                print "****************************************************************"
                for F in os.listdir(backup):
                    file = backup + '/' + F
                    try:
                        nb = load(file, compress=False)
                    except Exception, msg:
                        print "Failed to load backup '%s'"%file
                    else:
                        print "Successfully loaded backup '%s'"%file
                        nb.save()
                        break
                if nb is None:
                    print "Unable to restore notebook from *any* auto-saved backups."
                    print "This is a serious problem."
                nb.delete_doc_browser_worksheets()
                nb.set_directory(dir)
                nb.set_not_computing()
    if nb is None:
        nb = Notebook(dir)

    nb.address = address
    nb.port = port
    nb.secure = secure
    return nb


##########################################################
# Misc
##########################################################

def clean_name(name):
    return ''.join([x if (x.isalnum() or x == '_') else '_' for x in name])

def padzeros(s):
    return "0"*(3-len(str(s))) + str(s)

def sort_worksheet_list(v, sort, reverse):
    """
    INPUT:
        sort -- 'last_edited', 'owner', or 'name'
        reverse -- if True, reverse the order of the sort.
    """
    f = None
    if sort == 'last_edited':
        def c(a, b):
            return -cmp(a.last_edited(), b.last_edited())
        f = c
    elif sort == 'name':
        def c(a,b):
            return cmp((a.name().lower(), -a.last_edited()), (b.name().lower(), -b.last_edited()))
        f = c
    elif sort == 'owner':
        def c(a,b):
            return cmp((a.owner().lower(), -a.last_edited()), (b.owner().lower(), -b.last_edited()))
        f = c
    elif sort == "rating":
        def c(a,b):
            return -cmp((a.rating(), -a.last_edited()), (b.rating(), -b.last_edited()))
        f = c
    else:
        raise ValueError, "invalid sort key '%s'"%sort
    v.sort(cmp = f, reverse=reverse)
