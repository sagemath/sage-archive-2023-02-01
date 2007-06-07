"""
SAGE Notebook (Twisted Version)
"""
from twisted.web2 import server, http, resource, channel
from twisted.web2 import static, http_headers, responsecode

import css, js

from sage.misc.misc import SAGE_EXTCODE
javascript_path = SAGE_EXTCODE + "/notebook/javascript/"
css_path = SAGE_EXTCODE + "/notebook/css/"

class Main_css(resource.Resource):
    def render(self, ctx):
        s = css.css()
        return http.Response(stream=s)

class CSS(resource.Resource):
    def childFactory(self, request, name):
        return static.File(css_path + "/" + name)

setattr(CSS, 'child_main.css', Main_css())


class Main_js(resource.Resource):
    def render(self, ctx):
        s = js.javascript()
        return http.Response(stream=s)

class Javascript(resource.Resource):
    def childFactory(self, request, name):
        return static.File(javascript_path + "/" + name)

setattr(Javascript, 'child_main.js', Main_js())

class Toplevel(resource.Resource):
    addSlash = True

    child_javascript = Javascript()
    child_css = CSS()

    def render(self, ctx):
        s = notebook.html(authorized=True, worksheet_authorized=True)
        return http.Response(stream=s)



site = server.Site(Toplevel())
notebook = None  # this gets set on startup.


##########################################################
# This actually serves up the notebook.
  ##########################################################

from   sage.server.misc import print_open_msg
import os, socket

def notebook_twisted(directory='sage_notebook',
                     port=8000,
                     address='localhost',
                     port_tries=10):
    r"""
    Experimental twisted version of the SAGE Notebook.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    port = int(port)
    conf = '%s/twistedconf.py'%directory

    def run(port):
        ## Create the config file
        config = open(conf, 'w')
        config.write("""
import sage.server.notebook.notebook as notebook
import sage.server.notebook.twist as twist
twist.notebook = notebook.load_notebook('%s')

from twisted.web2 import channel
from twisted.application import service, strports
application = service.Application("SAGE Notebook")
s = strports.service('tcp:%s', channel.HTTPFactory(twist.site))
s.setServiceParent(application)
"""%(directory, port))
        config.close()

        ## Start up twisted
        print_open_msg(address, port)
        e = os.system('cd "%s" && sage -twistd -ny twistedconf.py'%directory)
        if e == 256:
            raise socket.error


    for i in range(port_tries):
        try:
            run(port + i)
        except socket.error:
            print "Port %s is already in use.  Trying next port..."%port
        else:
            break

    return True
