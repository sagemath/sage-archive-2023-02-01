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

PUBLIC_USER = 'pub'

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


    ##########################################################
    # Users
    ##########################################################
    def users(self):
        try:
            return self.__users
        except AttributeError:
            self.__users = {}
            return self.__users

    def user(self, username):
        try:
            return self.__users[username]
        except KeyError:
            U = user.User(username)
            self.__users[username] = U
            return U
        except AttributeError:
            self.__users = {}
            raise KeyError, "no user '%s'"%username

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

    def passwords(self):
        """
        Return the username:password dictionary.
        """
        return dict([(user.username(), user.password()) for user in self.user_list()])

    def user_conf(self, username):
        return self.users()[username].conf()

    ##########################################################
    # Moving, copying, creating, renaming, and listing worksheets
    ##########################################################

    def create_new_worksheet(self, worksheet_name, username):
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

    ##########################################################
    # Information about users
    ##########################################################
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
        if len(P) == 0:
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

    def user_history_text(self, username):
        H = self.user_history(username)
        return '\n\n'.join([L.strip() for L in H])

    def user_history_html(self, username):
        t = self.user_history_text(username)
        t = t.replace('<','&lt;')
        s = '<head>\n'
        s += '<title>Command History for %s</title>\n'%username
        s += '</head>\n'
        s += '<body>\n'
        s += '<pre>' + t + '</pre>\n'
        s += '<a name="bottom"></a>\n'
        s += '<script type="text/javascript"> window.location="#bottom"</script>\n'
        s += '</body>\n'
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
        if do_print:
            s += '<h1><a href=".">SAGE Worksheet: %s</a></h1>'%W.name()
        s += W.html(include_title=False, do_print=do_print)
        if do_print:
            s += '<script type="text/javascript">jsMath.Process();</script>\n'
        s += '\n</body>\n'
        return s


    def html_worksheet_list_for_user(self, user, active_only=True, sort='last_edited', reverse=False):
        W = self.get_worksheets_with_viewer(user)
        sort_worksheet_list(W, sort, reverse)  # changed W in place

        top = self.html_worksheet_list_top(user, active_only=active_only)
        list = self.html_worksheet_list(W, user, active_only=active_only, sort=sort, reverse=reverse)

        s = """
        <html>
           <link rel=stylesheet href="/css/main.css">
        <body>
        %s
        %s
        </body>
        </html>
        """%(top, list)

        return s


    def html_worksheet_list_top(self, user, active_only):
        s = ''

        s += self.html_user_control(user)
        s += self.html_banner()
        s += '<hr class="usercontrol">'
        s += self.html_new_or_upload()
        s += self.html_search()
        s += '<br>'
        s += '<hr class="usercontrol">'
        s += self.html_worksheet_acts()

        return s

    def html_user_control(self, user):
        s = ''
        s += '<div class="flush-right">'
        s += '<span class="username">%s</span>'%user
        s += vbar + '<a class="usercontrol" href="/settings">Settings</a>\n'
        s += vbar + '<a class="usercontrol" href="/doc">Help</a>\n'
        s += vbar + '<a class="usercontrol" href="/logout" class="help">Sign Out</a>\n'
        s += '</div>'
        return s

    def html_banner(self):
        s = """
        <span class="banner">
        <a class="banner" href="http://www.sagemath.org"><img align="top" src="/images/sagelogo.png" alt="SAGE"> Mathematics Software</a>
        </span>
        """
        return s

    def html_search(self):
        s = """
        <span class="flush-right">
        <input id="search_worksheets"></input>
        <button class="add_new_worksheet_menu" onClick="process_new_worksheet_menu_submit();">Search Worksheets</button>
        </span>
        """
        return s

    def html_new_or_upload(self):
        s = """
        <a class="boldusercontrol" href="/new_worksheet">New Worksheet</a>\n
        <a class="boldusercontrol" href="/upload">Upload</a>\n
        """
        return s

    def html_worksheet_acts(self):
        s = """
        <select>
         <option>Save</option>
         <option>Save as HTML (zipped) ... </option>
         <option>Save as LaTeX (zipped) ... </option>
         <option>Save as PDF...</option>
         <option>Save as Text...</option>
         <option>Copy Worksheet</option>
         <option>Archive</option>
         <option>Unarchive</option>
         <option>Un-collaborate me</option>
        </select>
        """
        s += '<button>Archive</button>'
        s += '&nbsp;&nbsp;<button>Delete</button>'
        return s


    def html_worksheet_list(self, worksheets, user, active_only, sort, reverse):
        s = ''

        s = '<br><br>'
        s += '<table width=100% border=0 cellspacing=0 cellpadding=0>'
        s += '<tr class="greybox"><td colspan=4><div class="thinspace"></div></td></tr>'
        s += '<tr  class="greybox">'
        s += '<td>&nbsp;<input class="entry" type=checkbox></td>'
        s += '<td><a class="listcontrol" href=".?sort=name%s">Active Worksheets</a> </td>'%(
            '' if sort != 'name' or reverse else '&reverse=True')
        s += '<td><a class="listcontrol" href=".?sort=owner%s">Owner / Collaborators / <i>Viewers</i></a> </td>'%(
            '' if sort != 'owner' or reverse else '&reverse=True')
        s += '<td><a class="listcontrol" href=".%s">Last Edited</a> </td>'%(
            '' if sort != 'last_edited' or reverse else '?reverse=True')
        s += '</tr>'
        s += '<tr class="greybox"><td colspan=4><div class="thinspace"></div></td></tr>'

        v = []
        for w in worksheets:
            k = '<tr>'
            k += '<td class="entry">%s</td>'%self.html_check_col(w)
            k += '<td>%s</td>'%self.html_worksheet_link(w)
            k += '<td>%s</td>'%self.html_owner_collab_view(w, user)
            k += '<td>%s</td>'%self.html_last_edited(w, user)
            k += '</tr>'
            k += '<tr class="thingreybox"><td colspan=4><div class="ultrathinspace"></div></td></tr>'
            v.append(k)

        s += ''.join(v)
        s += '</table>'

        return s

    def html_check_col(self, worksheet):
        def doc_options(name):
            return """
            <select>
            <option>Edit</option>
            <option>Collaborate</option>
            <option>Publish</option>
            <option>Revisions</option>
            <option>Preview</option>
            </select>
        """

        k = ''
        k += '<input type=checkbox>'
        k += '&nbsp;'*4
        k += doc_options(worksheet.filename())
        k += '&nbsp;'*4
        return k

    def html_worksheet_link(self, worksheet):
        return '<a class="worksheetname" href="/home/%s" target="_new">%s</a>\n'%(
              worksheet.filename(), worksheet.name())

    def html_owner_collab_view(self, worksheet, user):
        v = []

        owner = worksheet.owner()
        if owner == user:
            owner = "Me"

        v.append(owner)

        collab = worksheet.collaborators()

        if owner == "Me" or self.user(user).account_type() == 'admin':
            if len(collab) <= 1:
                share = '<a class="share" href="%s/share">Share now</a>'%(worksheet.filename_without_owner())
            else:
                collaborators = ', '.join([x for x in collab if x != user])
                v.append(collaborators)
                share = '<a class="share" href="%s/collaborate">Add</a>'%(worksheet.filename_without_owner())
        else:
            share = ''

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

        print "Saving notebook to '%s'..."%F
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


        body = ''

        body += '<div class="top_control_bar">\n'
        body += self.html_banner()
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
        body_html += """The format is as follows: <pre>
... Arbitrary HTML with latex formulas (in $ and $$)...
{{{meta info about cell|
Input
///
Output
}}}
</pre>"""
        body_html += '<form method="post" action="save" enctype="multipart/form-data">\n'
        body_html += '<input type="submit" value="Save Changes" name="button_save"/>\n'
        #body_html += '<input type="submit" value="Preview" name="button_preview"/>\n'
        body_html += '<input type="submit" value="Cancel" name="button_cancel"/>\n'
        body_html += '<textarea class="edit" id="cell_intext" rows="22" name="textfield">'+t+'</textarea>'
        body_html += '</form>'

        s = """
        <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
        <html><head><title>SAGE Wiki cell text </title>
        <style type="text/css">

        textarea.edit {
            font-family: courier, monospace;
            font-size:10pt;
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
    else:
        raise ValueError, "invalid sort key '%s'"%sort
    v.sort(cmp = f, reverse=reverse)
