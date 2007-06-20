#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

"""
SAGE Notebook (Twisted Version)
"""
import os, time

from twisted.web2 import server, http, resource, channel
from twisted.web2 import static, http_headers, responsecode

import css, js, keyboards

import notebook as _notebook

HISTORY_MAX_OUTPUT = 92*5
HISTORY_NCOLS = 90

from sage.misc.misc import SAGE_EXTCODE, walltime, tmp_filename
from sage.misc.remote_file import get_remote_file

p = os.path.join
css_path        = p(SAGE_EXTCODE, "notebook/css")
image_path      = p(SAGE_EXTCODE, "notebook/images")
javascript_path = p(SAGE_EXTCODE, "notebook/javascript")

# the list of users waiting to register
waiting = {}

# the user database
from user_db import UserDatabase
users = UserDatabase()

_cols = None
def word_wrap_cols():
    global _cols
    if _cols is None:
        _cols = notebook.conf()['word_wrap_cols']
    return _cols

############################
# Encoding data to go between the server and client
############################
SEP = '___S_A_G_E___'

def encode_list(v):
    return SEP.join([str(x) for x in v])



############################
# Notebook autosave.
############################
# save if make a change to notebook and at least some seconds have elapsed since last save.
def init_updates():
    global save_interval, idle_interval, last_save_time, last_idle_time

    save_interval = notebook.conf()['save_interval']
    idle_interval = notebook.conf()['idle_check_interval']
    last_save_time = walltime()
    last_idle_time = walltime()

def notebook_save_check():
    global last_save_time
    t = walltime()
    if t > last_save_time + save_interval:
        notebook.save()
        last_save_time = t

def notebook_idle_check():
    notebook.quit_idle_worksheet_processes()
    global last_idle_time
    t = walltime()
    if t > last_idle_time + idle_interval:
        notebook.save()
        last_idle_time = t

def notebook_updates():
    notebook_save_check()
    notebook_idle_check()

######################################################################################
# RESOURCES
######################################################################################

############################
# Create a SAGE worksheet from a latex2html'd file
############################
from docHTMLProcessor import DocHTMLProcessor

def doc_worksheet():
    wnames = notebook.worksheet_names()
    name = 'doc_browser_0'
    if name in wnames:
        W = notebook.get_worksheet_with_name(name)
        W.restart_sage()
        W.clear()
        return W
    else:
        return notebook.create_new_worksheet(name, username, docbrowser=True)


class WorksheetFile(resource.Resource):
    addSlash = False

    def __init__(self, path):
        self.docpath = path

    def render(self, ctx=None):
        # Create a live SAGE worksheet out of self.path and render it.
        doc_page_html = open(self.docpath).read()
        directory = os.path.split(self.docpath)[0]
        doc_page, css_href = DocHTMLProcessor().process_doc_html(DOC,
                               directory, doc_page_html)
        if css_href:
            css_href = DOC + directory + css_href
        W = doc_worksheet()
        W.edit_save(doc_page)
        s = notebook.html(worksheet_filename = W.filename(),  username=username)
        return http.Response(stream=s)

    def childFactory(self, request, name):
        path = self.docpath + '/' + name
        if name.endswith('.html'):
            return WorksheetFile(path)
        else:
            return static.File(path)

############################
# The documentation browsers
############################

DOC = os.path.abspath(os.environ['SAGE_ROOT'] + '/doc/')
class DocStatic(resource.Resource):
    addSlash = True
    def render(self, ctx):
        return static.File('%s/index.html'%DOC)

    def childFactory(self, request, name):
        return static.File('%s/%s'%(DOC,name))

class DocLive(resource.Resource):
    addSlash = True

    def render(self, ctx):
        return static.File('%s/index.html'%DOC)

    def childFactory(self, request, name):
        return WorksheetFile('%s/%s'%(DOC,name))

class Doc(resource.Resource):
    addSlash = True
    child_static = DocStatic()
    child_live = DocLive()

    def render(self, ctx):
        s = """
        <html>
           <link rel=stylesheet href="/css/main.css">
           <title>SAGE Help</title>
        <body>
        <a href="/">Home</a>
        <br>
        <center>
        <h1><font color="darkred">SAGE Documentation</font></h1>
        <br>
        <hr class="usercontrol">
        <br><br>
        <font size=+2>
        <a href="/help/">SAGE Notebook Quickstart</a><br><br>
        <a href="/doc/static/">Static Documentation</a><br><br>
        <a href="/doc/live/">Interactive Live Documentation</a><br>
        <br>
        <hr class="usercontrol">
        </font>
        </center>
        </body>
        </html>
        """
        return http.Response(stream=s)

############################
# A New Worksheet
############################
class NewWorksheet(resource.Resource):
    def render(self, ctx):
        W = notebook.create_new_worksheet("Untitled", username)
        return http.RedirectResponse('/home/'+W.filename())


############################
# Uploading a saved worksheet file
############################

def redirect(url):
    return '<html><head><meta http-equiv="REFRESH" content="0; URL=%s"></head></html>'%url

class Upload(resource.Resource):
    def render(self, ctx):
        return http.Response(stream = notebook.upload_window())

class UploadWorksheet(resource.PostableResource):
    def render(self, ctx):
        url = ctx.args['urlField'][0].strip()
        if url != '':
            tmp = get_remote_file(url, verbose=True)
        else:
            tmp = '%s/tmp.sws'%notebook.directory()
            f = file(tmp,'wb')
            f.write(ctx.files['fileField'][0][2].read())
            f.close()

        try:
            W = notebook.import_worksheet(tmp, username)
            os.unlink(tmp)
        except ValueError, msg:
            s = "<html>Error uploading worksheet '%s'.  <a href='/'>continue</a></html>"%msg
            return http.Response(stream = s)

        name = ctx.args['nameField'][0].strip()
        if len(name) > 0:
            W.set_name(name)

        return http.RedirectResponse('/home/'+W.filename())



############################
# A resource attached to a given worksheet.
#
# This has the name of the worksheet and the
# worksheet object itself set as attributes.
# It's much better to do it once-for-all here
# instead of doing it in the derived classes
# over and over again.
############################
class WorksheetResource:
    def __init__(self, name):
        self.name = name
        self.worksheet = notebook.get_worksheet_with_filename(name)
        if not username in self.worksheet.collaborators() and user_type(username) != 'admin':
            raise RuntimeError, "illegal worksheet access"

    def id(self, ctx):
        return int(ctx.args['id'][0])


###############################################
# Worksheet data -- a file that
# is associated with a cell in some worksheet.
# The file is stored on the filesystem.
#      /home/worksheet_name/data/cell_number/filename
##############################################
class CellData(resource.Resource):
    def __init__(self, worksheet, number):
        self.worksheet = worksheet
        self.number = number

    def childFactory(self, request, name):
        dir = self.worksheet.directory()
        path = '%s/cells/%s/%s'%(dir, self.number, name)
        return static.File(path)

class Worksheet_data(WorksheetResource, resource.Resource):
    def childFactory(self, request, number):
        return CellData(self.worksheet, number)

########################################################
# Use this to wrap a worksheet operation in a confirmation
# request.  See WorksheetDelete and WorksheetAdd for
# examples.
########################################################

class FastRedirect(resource.Resource):
    def __init__(self, dest):
        self.dest = dest
    def render(self, ctx):
        return http.RedirectResponse(self.dest)

class YesNo(resource.Resource):
    addSlash = True

    def __init__(self, mesg, yes, no):
        self.mesg = mesg
        self.yes = yes
        self.no  = no

    def render(self, ctx):
        from sage.server.notebook.template import yes_no_template
        lt = yes_no_template(mesg=self.mesg)
        return http.Response(stream = lt)

        s = '<html><body>%s<br>'%self.mesg
        s += '<a href="yes">Yes</a> or <a href="no">No</a></body></html>'
        return http.Response(stream = s)

    def childFactory(self, request, op):
        if op == 'yes':
            return FastRedirect(self.yes())
        elif op == 'no':
            return FastRedirect(self.no())


########################################################
# Completely delete the worksheet from the notebook
# server.  It is assumed that the javascript has already
# done all relevant confirmation, and of course in the
# future we'll check that the users is authenticated to
# touch the worksheet.  Delete should also, in the
# future, really just put the result in a trash can.
########################################################

def Worksheet_delete(name):
    def yes():
        notebook.delete_worksheet(name)
        return '/'
    def no():
        return '..'

    return YesNo('Do you want to delete the worksheet "%s"?'%name,
                 yes=yes, no=no)

########################################################
# Create a new worksheet.
########################################################

## def Worksheet_create(name):
##     def yes():
##         W = notebook.create_new_worksheet(name, username)
##         return '/home/' + W.filename()
##     def no():
##         return '/'

##     wsname = _notebook.clean_name(name)
##     return YesNo('Do you want to create the worksheet "%s"?'%name,
##                  yes = yes, no = no)

########################################################
# Worksheet configuration.
########################################################
class Worksheet_conf(WorksheetResource, resource.Resource):
    def render(self, ctx):
        conf = self.worksheet.conf()
        s = str(conf)
        # TODO: This should be a form that allows for configuring all options
        # of a given worksheet, saves the result,
        return http.Response(stream = s)

########################################################
# Cell introspection
########################################################
class Worksheet_introspect(WorksheetResource, resource.PostableResource):
    """
    Cell introspection.  This is called when the user presses the tab
    key in the browser in order to introspect.
    """
    def render(self, ctx):
        try:
            id = int(ctx.args['id'][0])
        except (KeyError,TypeError):
            return http.Response(stream = 'Error in introspection -- invalid cell id.')
        try:
            before_cursor = ctx.args['before_cursor'][0]
        except KeyError:
            before_cursor = ''
        try:
            after_cursor = ctx.args['after_cursor'][0]
        except KeyError:
            after_cursor = ''
        C = self.worksheet.get_cell_with_id(id)
        C.evaluate(introspect=[before_cursor, after_cursor])
        return http.Response(stream = encode_list([C.next_id(),'no_new_cell',id]))

########################################################
# Edit the entire worksheet
########################################################
class Worksheet_edit(WorksheetResource, resource.Resource):
    """
    Return a window that allows the user to edit the text of the
    worksheet with the given filename.
    """
    def render(self, ctx):
        return http.Response(stream = notebook.edit_window(self.worksheet))


########################################################
# Preview what the worksheet will look like when published.
########################################################
class Worksheet_preview(WorksheetResource, resource.PostableResource):
    def render(self, ctx):
        s = notebook.html(worksheet_filename = self.name,  username='pub')
        return http.Response(stream=s)


########################################################
# Copy a worksheet
########################################################
class Worksheet_copy(WorksheetResource, resource.PostableResource):
    def render(self, ctx):
        W = notebook.copy_worksheet(self.worksheet, username)
        return http.RedirectResponse('/home/' + W.filename())

########################################################
# Get a copy of a published worksheet and start editing it
########################################################
class Worksheet_edit_published_page(WorksheetResource, resource.Resource):
    def render(self, ctx):
        W = notebook.copy_worksheet(self.worksheet, username)
        W.set_name(self.worksheet.name())
        return http.RedirectResponse('/home/' + W.filename())

########################################################
# Save a worksheet
########################################################
class Worksheet_save(WorksheetResource, resource.PostableResource):
    """
    Save the contents of a worksheet after editing it in plain-text edit mode.
    """
    def render(self, ctx):
        if ctx.args.has_key('button_save'):
            self.worksheet.edit_save(ctx.args['textfield'][0])
        return http.RedirectResponse('/home/'+self.worksheet.filename())


class Worksheet_save_snapshot(WorksheetResource, resource.PostableResource):
    """
    Save a snapshot of a worksheet.
    """
    def render(self, ctx):
        self.worksheet.save_snapshot()
        return http.Response(stream="saved")

class Worksheet_revert_to_last_saved_state(WorksheetResource, resource.PostableResource):
    def render(self, ctx):
        self.worksheet.revert_to_last_saved_state()
        return http.Response(stream="reverted")

class Worksheet_save_and_close(WorksheetResource, resource.PostableResource):
    """
    Save a snapshot of a worksheet then quit it.
    """
    def render(self, ctx):
        self.worksheet.save_snapshot()
        self.worksheet.quit()
        return http.Response(stream="saved")

########################################################
# Collaborate with others
########################################################
class Worksheet_share(WorksheetResource, resource.Resource):
    def render(self, ctx):
        s = notebook.html_share(self.worksheet, username)
        return http.Response(stream = s)

class Worksheet_invite_collab(WorksheetResource, resource.PostableResource):
    def render(self, ctx):
        collab = ctx.args['collaborators'][0]
        v = [x.strip() for x in collab.split(',')]
        self.worksheet.add_collaborators(v)
        return http.RedirectResponse('share')


#################################
# Revisions
#################################

class PublishWorksheetRevision(resource.Resource):
    def __init__(self, worksheet, rev):
        self.worksheet = worksheet
        self.rev = rev

    def render(self, ctx):
        W = notebook.publish_worksheet(self.worksheet)
        txt = open(self.worksheet.get_snapshot_text_filename(self.rev)).read()
        W.delete_cells_directory()
        W.edit_save(txt)
        return http.RedirectResponse('/home/'+W.filename())

class RevertToWorksheetRevision(resource.Resource):
    def __init__(self, worksheet, rev):
        self.worksheet = worksheet
        self.rev = rev

    def render(self, ctx):
        self.worksheet.save_snapshot()
        txt = open(self.worksheet.get_snapshot_text_filename(self.rev)).read()
        self.worksheet.delete_cells_directory()
        self.worksheet.edit_save(txt)
        return http.RedirectResponse('/home/'+self.worksheet.filename())

def worksheet_revision_publish(worksheet, rev):
    W = notebook.publish_worksheet(worksheet)
    txt = open(worksheet.get_snapshot_text_filename(rev)).read()
    W.delete_cells_directory()
    W.edit_save(txt)
    return http.RedirectResponse('/home/'+W.filename())

def worksheet_revision_revert(worksheet, rev):
    worksheet.save_snapshot()
    txt = open(worksheet.get_snapshot_text_filename(rev)).read()
    worksheet.delete_cells_directory()
    worksheet.edit_save(txt)
    return http.RedirectResponse('/home/'+worksheet.filename())


class Worksheet_revisions(WorksheetResource, resource.PostableResource):
    """
    Show a list of revisions of this worksheet.
    """
    def render(self, ctx):
        if not ctx.args.has_key('action'):
            if ctx.args.has_key('rev'):
                rev = ctx.args['rev'][0]
                s = notebook.html_specific_revision(username, self.worksheet, rev)
            else:
                s = notebook.html_worksheet_revision_list(username, self.worksheet)
        else:
            rev = ctx.args['rev'][0]
            action = ctx.args['action'][0]
            print action
            if action == 'revert':
                return worksheet_revision_revert(self.worksheet, rev)
            elif action == 'publish':
                return worksheet_revision_publish(self.worksheet, rev)
            else:
                s = 'Error'
        return http.Response(stream = s)


########################################################
# Worksheet/User/Notebooks settings and configuration
########################################################

class Worksheet_input_settings(WorksheetResource, resource.PostableResource):
    def render(self, ctx):
        if self.worksheet.owner() != username:
            s = 'You must be the owner of this worksheet to configure it.'
            return http.Response(stream = s)
        else:
            system = ctx.args['system'][0].strip()
            self.worksheet.set_system(system)
            return http.RedirectResponse('/home/'+ self.worksheet.filename())

class Worksheet_settings(WorksheetResource, resource.Resource):
    def render(self, ctx):
        if self.worksheet.owner() != username:
            s = 'You must be the owner of this worksheet to configure it.'
        else:
            s = notebook.html_worksheet_settings(self.worksheet, username)
        return http.Response(stream = s)

class ProcessUserSettings(resource.PostableResource):
    def render(self, ctx):
        pass

class UserSettings(resource.Resource):
    child_process = ProcessUserSettings()

    def render(self, ctx):
        s = notebook.html_user_settings(username)
        return http.Response(stream = s)

class ProcessNotebookSettings(resource.PostableResource):
    def render(self, ctx):
        pass

class NotebookSettings(resource.Resource):
    child_process = ProcessNotebookSettings()
    def render(self, ctx):
        if user_type(username) != 'admin':
            s = 'You must an admin to configure the notebook.'
        else:
            s = notebook.html_notebook_settings()
        return http.Response(stream = s)


########################################################
# Set output type of a cell
########################################################
class Worksheet_set_cell_output_type(WorksheetResource, resource.PostableResource):
    """
    Set the output type of the cell.

    This enables the type of output cell, such as to allowing wrapping
    for output that is very long.
    """
    def render(self, ctx):
        id = self.id(ctx)
        typ = ctx.args['type'][0]
        W = self.worksheet
        W.get_cell_with_id(id).set_cell_output_type(typ)
        return http.Response(stream = '')

########################################################
# The new cell command: /home/worksheet/new_cell?id=number
########################################################
class Worksheet_new_cell(WorksheetResource, resource.PostableResource):
    """
    Adds a new cell before a given cell.
    """
    def render(self, ctx):
        id = self.id(ctx)
        cell = self.worksheet.new_cell_before(id)
        s = encode_list([cell.id(), cell.html(div_wrap=False), id])
        return http.Response(stream = s)


########################################################
# The delete cell command: /home/worksheet/delete_cell?id=number
########################################################
class Worksheet_delete_cell(WorksheetResource, resource.PostableResource):
    """
    Deletes a notebook cell.

    If there is only one cell left in a given worksheet, the request
    to delete that cell is ignored because there must be a least one
    cell at all times in a worksheet.  (This requirement exists so
    other functions that evaluate relative to existing cells will
    still work, and so one can add new cells.)
    """
    def render(self, ctx):
        id = self.id(ctx)
        W = self.worksheet
        if len(W) <= 1:
            s = 'ignore'
        else:
            prev_id = W.delete_cell_with_id(id)
            s = encode_list(['delete', id, prev_id, W.cell_id_list()])
        return http.Response(stream = s)


############################
# Get the latest update on output appearing
# in a given output cell.
############################
class Worksheet_cell_update(WorksheetResource, resource.PostableResource):
    def render(self, ctx):
        id = self.id(ctx)

        worksheet = self.worksheet

        # update the computation one "step".
        worksheet.check_comp()

        # now get latest status on our cell
        status, cell = worksheet.check_cell(id)

        if status == 'd':
            new_input = cell.changed_input_text()
            out_html = cell.output_html()
            H = "Worksheet '%s' (%s)\n"%(worksheet.name(), time.strftime("%Y-%m-%d at %H:%M",time.localtime(time.time())))
            H += cell.edit_text(ncols=HISTORY_NCOLS, prompts=False,
                                max_out=HISTORY_MAX_OUTPUT)
            notebook.add_to_user_history(H, username)
        else:
            new_input = ''
            out_html = ''

        if cell.interrupted():
            inter = 'true'
        else:
            inter = 'false'

        raw = cell.output_text(raw=True).split("\n")
        if "Unhandled SIGSEGV" in raw:
            inter = 'restart'
            print "Segmentation fault detected in output!"

        msg = '%s%s %s'%(status, cell.id(),
                       encode_list([cell.output_text(html=True),
                                    cell.output_text(word_wrap_cols(), html=True),
                                    out_html,
                                    new_input,
                                    inter,
                                    cell.introspect_html()]))

        # There may be more computations left to do, so start one if there is one.
        worksheet.start_next_comp()

        return http.Response(stream=msg)


class Worksheet_eval(WorksheetResource, resource.PostableResource):
    """
    Evaluate a worksheet cell.

    If the request is not authorized, (the requester did not enter the
    correct password for the given worksheet), then the request to
    evaluate or introspect the cell is ignored.

    If the cell contains either 1 or 2 question marks at the end (not
    on a comment line), then this is interpreted as a request for
    either introspection to the documentation of the function, or the
    documentation of the function and the source code of the function
    respectively.
    """
    def render(self, ctx):
        newcell = int(ctx.args['newcell'][0])  # whether to insert a new cell or not
        id = self.id(ctx)
        if not ctx.args.has_key('input'):
            input_text = ''
        else:
            input_text = ctx.args['input'][0]
            input_text = input_text.replace('\r\n', '\n')   # DOS

        W = self.worksheet
        cell = W.get_cell_with_id(id)
        cell.set_input_text(input_text)
        cell.evaluate(username=username)

        if cell.is_last():
            new_cell = W.append_new_cell()
            s = encode_list([new_cell.id(), 'append_new_cell', new_cell.html(div_wrap=False)])
        elif newcell:
            new_cell = W.new_cell_after(id)
            s = encode_list([new_cell.id(), 'insert_cell', new_cell.html(div_wrap=False), str(id)])
        else:
            s = encode_list([cell.next_id(), 'no_new_cell', str(id)])

        notebook_updates()
        return http.Response(stream=s)


########################################################
# Publication and rating of a worksheet
########################################################

class Worksheet_publish(WorksheetResource, resource.Resource):
    def render(self, ctx):
        W = notebook.publish_worksheet(self.worksheet)
        addr = '/home/' + W.filename()
        return http.RedirectResponse(addr)


class Worksheet_rating_info(WorksheetResource, resource.Resource):
    def render(self, ctx):
        ratings = self.worksheet.ratings()
        r = '\n'.join(['<tr><td>%s</td><td align=center>%s</td></tr>'%(x,y) for x,y in sorted(ratings)])
        return http.Response(stream="""
        <html>
        <h2 align=center>Ratings for %s</h2>
        <h3 align=center><a href='/home/%s'>Return to the worksheet.</a>
        <br><br>
        <table width=30%% align=center border=1 cellpadding=10 cellspacing=0>
        <tr bgcolor="lightgray"><td width=85%%>User</td><td width=10em align=center>Rating</td></tr>
        %s
        </table>
        <br><br>
        </html>
        """%(self.worksheet.name(), self.worksheet.filename(), r))


class WorksheetRating(WorksheetResource, resource.Resource):
    def render(self, ctx):
        if self.worksheet.is_rater(username):
            return http.Response(stream="You have already rated the worksheet '%s'."%self.worksheet.name())
        self.do_rating()
        return http.Response(stream="""
        <html>
        <body>
        <br><br><br>
        Thank you for rating the worksheet '%s'!
        <br><br>
        <a href='/home/%s'>Click here to return to the worksheet</a> or <a href="/pub">browse other published worksheets</a>.
        </body>
        </html>
        """%(self.worksheet.name(), self.worksheet.filename()))

class Worksheet_rate1(WorksheetRating):
    def do_rating(self):
        self.worksheet.rate(1, username)

class Worksheet_rate2(WorksheetRating):
    def do_rating(self):
        self.worksheet.rate(2, username)

class Worksheet_rate3(WorksheetRating):
    def do_rating(self):
        self.worksheet.rate(3, username)

class Worksheet_rate4(WorksheetRating):
    def do_rating(self):
        self.worksheet.rate(4, username)

########################################################
# Downloading, moving around, renaming, etc.
########################################################


class Worksheet_download(WorksheetResource, resource.Resource):
    def childFactory(self, request, name):
        worksheet_name = self.name
        filename = tmp_filename() + '.sws'
        try:
            notebook.export_worksheet(worksheet_name, filename)
        except KeyError:
            return http.Response(stream='No such worksheet.')
        r = open(filename, 'rb').read()
        os.unlink(filename)
        return static.Data(r, 'application/sage')
        #return static.File(filename)

class Worksheet_rename(WorksheetResource, resource.PostableResource):
    def render(self, ctx):
        self.worksheet.set_name(ctx.args['name'][0])
        return http.Response(stream='done')

class Worksheet_restart_sage(WorksheetResource, resource.Resource):
    def render(self, ctx):
        # TODO -- this must not block long (!)
        self.worksheet.restart_sage()
        return http.Response(stream='done')

class Worksheet_interrupt(WorksheetResource, resource.Resource):
    def render(self, ctx):
        # TODO -- this must not block long (!)
        s = self.worksheet.interrupt()
        return http.Response(stream='ok' if s else 'failed')

class Worksheet_plain(WorksheetResource, resource.Resource):
    def render(self, ctx):
        s = notebook.plain_text_worksheet_html(self.name)
        return http.Response(stream=s)

class Worksheet_print(WorksheetResource, resource.Resource):
    def render(self, ctx):
        s = notebook.worksheet_html(self.name, do_print=True)
        return http.Response(stream=s)


class NotImplementedWorksheetOp(resource.Resource):
    def __init__(self, op):
        self.op = op

    def render(self, ctx):
        return http.Response(stream = 'The worksheet operation "%s" is not implemented.'%self.op)


class Worksheet(WorksheetResource, resource.Resource):
    addSlash = True

    def render(self, ctx):
        s = notebook.html(worksheet_filename = self.name,  username=username)
        self.worksheet.sage()
        return http.Response(stream=s)

    def childFactory(self, request, op):
        notebook_updates()
        try:
            # Rather than a bunch of if-else statements, we wrap
            # any Worksheet_... class as a subresource of a worksheet
            # using the following  statement:
            R = globals()['Worksheet_%s'%op]
            return R(self.name)
        except KeyError:
            file = self.worksheet.data_directory() + '/' + op
            if os.path.exists(file):
                return static.File(file)

            return NotImplementedWorksheetOp(op)

def render_worksheet_list(args, pub=False):
    if args.has_key('typ'):
        typ = args['typ'][0]
    else:
        typ = 'active'
    if args.has_key('search'):
        search = args['search'][0]
    else:
        search = None
    if not args.has_key('sort'):
        sort = 'last_edited'
    else:
        sort = args['sort'][0]
    if args.has_key('reverse'):
        reverse = (args['reverse'][0] == 'True')
        print args['reverse'], reverse
    else:
        reverse = False

    if pub:
        if username is None or username == tuple([]):
            user = 'pub'
        else:
            user = username
        s = notebook.html_worksheet_list_public(
            user, sort=sort, reverse=reverse, search=search)
    else:
        s = notebook.html_worksheet_list_for_user(
            username, typ=typ, sort=sort, reverse=reverse, search=search)
    return s


class WorksheetsByUser(resource.Resource):
    addSlash = True

    def __init__(self, user):
        self.user = user

    def render_list(self, ctx):
        s = render_worksheet_list(ctx.args)
        return http.Response(stream = s)

    def render(self, ctx):
        if self.user == username:
            return self.render_list(ctx)
        else:
            return http.Response(stream = "<html><br><br><br><h2>You are logged in as '%s' so you do not have permission to view the home page of '%s'.</h2></html>."%(
                username, self.user))

    def childFactory(self, request, name):
        if name == "trash":
            return TrashCan(self.user)

        filename = self.user + '/' + name
        try:
            return Worksheet(filename)
        except KeyError:
            return http.Response(stream = "The user '%s' has no worksheet '%s'."%(self.user, name))


############################
# Trash can, archive and active
############################
class SendWorksheetToFolder(resource.PostableResource):
    def action(self, W):
        raise NotImplementedError
    def render(self, ctx):
        filename = ctx.args['filename'][0]
        W = notebook.get_worksheet_with_filename(filename)
        X = notebook.user(username)
        if X.is_admin() or W.owner() == username or W.publisher() == username:
            self.action(W)
            s = ''
        else:
            s = "You are not authorized to move '%s'"%W.name()
        return http.Response(stream = s)

class SendWorksheetToTrash(SendWorksheetToFolder):
    def action(self, W):
        W.move_to_trash()

class SendWorksheetToArchive(SendWorksheetToFolder):
    def action(self, W):
        W.move_to_archive()

class SendWorksheetToActive(SendWorksheetToFolder):
    def action(self, W):
        W.set_active()

############################
# Publically Available Worksheets
############################
class PublicWorksheets(resource.Resource):
    addSlash = True

    def render(self, ctx):
        s = render_worksheet_list(ctx.args, pub=True)
        return http.Response(stream = s)

    def childFactory(self, request, name):
        return http.Response(stream = "Not implemented")

############################
# Resource that gives access to worksheets
############################

class Worksheets(resource.Resource):
    def render(self, ctx):
        return http.Response(stream = "Please request a specific worksheet")

    def childFactory(self, request, name):
        return WorksheetsByUser(name)


class WorksheetsByUserAdmin(WorksheetsByUser):
    def render(self, ctx):
        return self.render_list(ctx)

class WorksheetsAdmin(Worksheets):
    def childFactory(self, request, name):
        return WorksheetsByUserAdmin(name)

############################
# Adding a new worksheet
############################

class NotebookConf(Worksheets):
    def render(self, ctx):
        s = '<html>' + notebook.conf().html_conf_form('submit') + '</html>'
        return http.Response(stream = s)



############################
# Adding a new worksheet
############################

class AddWorksheet(resource.Resource):
    def render(self, ctx):
        name = ctx.args['name'][0]
        W = notebook.create_new_worksheet(name)
        v = notebook.worksheet_list_html(W.name())
        return http.Response(stream = encode_list([v, W.name()]))

class Notebook(resource.Resource):
    child_add_worksheet = AddWorksheet()


############################

class Help(resource.Resource):
    addSlash = True
    def render(self, ctx):
        try:
            s = self._cache
        except AttributeError:
            s = notebook.help_window()
            self._cache = s
        return http.Response(stream=s)


############################

############################

class History(resource.Resource):
    def render(self, ctx):
        s = notebook.user_history_html(username)
        return http.Response(stream=s)

class LiveHistory(resource.Resource):
    def render(self, ctx):
        W = notebook.create_new_worksheet_from_history('Log', username, 100)
        return http.RedirectResponse('/home/'+W.filename())


############################

class Main_css(resource.Resource):
    def render(self, ctx):
        s = css.css()
        return http.Response(stream=s)

class CSS(resource.Resource):
    addSlash = True

    def render(self, ctx):
        return static.File(css_path)

    def childFactory(self, request, name):
        return static.File(css_path + "/" + name)

setattr(CSS, 'child_main.css', Main_css())

############################


############################
# Javascript resources
############################

class Main_js(resource.Resource):
    def render(self, ctx):
        s = js.javascript()
        return http.Response(stream=s)

class Keyboard_js_specific(resource.Resource):
    def __init__(self, browser_os):
        self.s = keyboards.get_keyboard(browser_os)

    def render(self, ctx):
        return http.Response(stream = self.s)


class Keyboard_js(resource.Resource):
    def childFactory(self, request, browser_os):
        return Keyboard_js_specific(browser_os)

class Javascript(resource.Resource):
    addSlash = True
    child_keyboard = Keyboard_js()

    def render(self, ctx):
        return static.File(javascript_path)

    def childFactory(self, request, name):
        return static.File(javascript_path + "/" + name)

setattr(Javascript, 'child_main.js', Main_js())

############################
# Logout
############################
class Logout(resource.Resource):
    def render(self, ctx):
        # TODO -- actually log out.
        notebook.save()
        return http.Response(stream="thank you for using sage")

############################
# Image resource
############################

class Images(resource.Resource):
    addSlash = True

    def render(self, ctx):
        return static.File(image_path)

    def childFactory(self, request, name):
        return static.File(image_path + "/" + name)

#####################################
# Confirmation of registration
####################################
class RegConfirmation(resource.Resource):
    def render(self, request):
        key = request.args['key'][0]
        global notebook
        url_prefix = "https" if notebook.secure else "http"
        invalid_confirm_key = """\
<html>
<h1>Invalid confirmation key</h1>
<p>You are reporting a confirmation key that has not been assigned by this
server. Please <a href="%s://%s:%s/register">register</a> with the server.</p>
</html>""" % (url_prefix, notebook.address, notebook.port)
        key = int(key)
        global waiting
        try:
            username = waiting[key]
        except KeyError:
            return http.Response(stream=invalid_confirm_key)
        success = """\
<html>
<h1>Hello, %s. Thank you for registering!</h1>
</html>""" % username
        del waiting[key]
        return http.Response(stream=success)

############################
# Registration page
############################

class RegistrationPage(resource.PostableResource):
    def __init__(self, userdb):
        self.userdb = userdb

    def render(self, request):
        if request.args.has_key('email'):
            if request.args['email'][0] is not None:

                s = ''
                try:
                    username = request.args['username'][0]
                except KeyError:
                    s += "You must specify a username."
                try:
                    passwd  = request.args['password'][0]
                except KeyError:
                    s += "  You must specify a password."
                else:
                    if len(passwd) == 0:
                        s = "  Password must be nonempty."
                if s:
                    return http.Response(stream=s)


                destaddr = """%s""" % request.args['email'][0]
                from sage.server.notebook.smtpsend import send_mail
                from sage.server.notebook.register import make_key, build_msg
                # TODO: make this come from the server settings
                key = make_key()
                listenaddr = notebook.address
                port = notebook.port
                fromaddr = 'no-reply@%s' % listenaddr
                body = build_msg(key, username, listenaddr, port,
                                 notebook.secure)

                # Send a confirmation message to the user.
                try:
                    send_mail(self, fromaddr, destaddr, "SAGE Notebook Registration",body)
                except ValueError:
                    # the email address is invalid
                    return http.Response(stream="Registration failed -- the email address '%s' is invalid."%destaddr)

                # Store in memory that we are waiting for the user to respond
                # to their invitation to join the SAGE notebook.
                waiting[key] = username

            # Add the user to passwords.txt
            try:
                self.userdb.add_user(username, passwd, destaddr)
                # now say that the user has been registered.
                s = """
                <html>
                <h1>Registration information received</h1>
                <p>Thank you for registering with the SAGE notebook. A
                confirmation message will be sent to %s.</p>
                <br>
                <p><a href="/">Click here to login with your new account.</a></p>
                </html>
                """%destaddr
            except ValueError:
                s = """
                <html>
                <h1>Username is already taken, please choose another one.</h1>
                </html>
                """
        else:
            url_prefix = "https" if notebook.secure else "http"
            s = """<html><h1>This is the registration page.</h1>
            <form method="POST" action="%s://%s:%s/register">
            Username: <input type="text" name="username" size="15" />  Password:
                <input type="password" name="password" size="15" /><br /> Email
                Address: <input type="text" name="email" size="15" /><br /> <div align="center">  <p><input type="submit" value="Register" /></p>  </div> </form><br /><br />
            </html>""" % (url_prefix, notebook.address, notebook.port)
        return http.Response(stream=s)

class Toplevel(resource.PostableResource):
    def __init__(self, cookie, _username):
        self.cookie = cookie
        global username, admin
        username = _username

setattr(Toplevel, 'child_favicon.ico', static.File(image_path + '/favicon.ico'))



from sage.server.notebook.template import login_template
from sage.server.notebook.template import failed_login_template

class LoginResourceClass(resource.Resource):
    def render(self, ctx):
        return http.Response(stream =  login_template())

    def childFactory(self, request, name):
        return LoginResource

LoginResource = LoginResourceClass()

class AnonymousToplevel(Toplevel):
    from sage.server.notebook.avatars import PasswordChecker
    addSlash = True
    child_register = RegistrationPage(PasswordChecker())
    child_confirm = RegConfirmation()
    child_images = Images()
    child_css = CSS()
    child_pub = PublicWorksheets()

    def render(self, ctx):
        return http.Response(stream =  login_template())

    def childFactory(self, request, name):
        return LoginResource

class FailedToplevel(Toplevel):
    def __init__(self, info, problem):
        self.info = info
        self.problem= problem

    def render(self, ctx):
        return http.Response(stream=failed_login_template(problem=self.problem))

class UserToplevel(Toplevel):
    addSlash = True

    child_images = Images()
    child_css = CSS()
    child_javascript = Javascript()
    child_home = Worksheets()
    child_notebook = Notebook()
    child_doc = Doc()
    child_upload = Upload()
    child_upload_worksheet = UploadWorksheet()
    child_new_worksheet = NewWorksheet()
    child_logout = Logout()
    child_pub = PublicWorksheets()
    child_live_history = LiveHistory()
    child_history = History()
    child_help = Help()
    child_send_to_trash = SendWorksheetToTrash()
    child_send_to_archive = SendWorksheetToArchive()
    child_send_to_active = SendWorksheetToActive()
    child_notebook_settings = NotebookSettings()
    child_settings = UserSettings()

    def render(self, ctx):
        s = render_worksheet_list(ctx.args)
        return http.Response(responsecode.OK,
                             {'content-type': http_headers.MimeType('text',
                                                                    'html'),
                             'set-cookie':[http_headers.Cookie("sid",
                                                            self.cookie)]},
                             stream=s)

class AdminToplevel(UserToplevel):
    child_home = WorksheetsAdmin()
    child_conf = NotebookConf()

    def render(self, ctx):
        s = render_worksheet_list(ctx.args)
        return http.Response(responsecode.OK,
                             {'content-type': http_headers.MimeType('text',
                                                                    'html'),
                             'set-cookie':[http_headers.Cookie("sid",
                                                            self.cookie)]},
                             stream=s)



notebook = None  # this gets set on startup.
username = None  # This is set when a request comes in.
OPEN_MODE = None # this gets set on startup.


def user_type(user):
    return notebook.user(user).account_type()
