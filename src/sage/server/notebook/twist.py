"""
SAGE Notebook (Twisted Version)
"""

from twisted.web2 import server, http, resource, channel
from twisted.web2 import static, http_headers, responsecode

import css

class MainCSS(resource.Resource):
    def render(self, ctx):
        s = css.css()
        return http.Response(stream=s)

class Toplevel(resource.Resource):
    addSlash = True
    def render(self, ctx):
        s = notebook.html(authorized=True)
        return http.Response(stream=s)

setattr(Toplevel, 'child___main__.css', MainCSS())


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
