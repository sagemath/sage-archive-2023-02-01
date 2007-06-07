"""
SAGE Notebook (Twisted Version)
"""

import notebook

from twisted.web2 import server, http, resource, channel
from twisted.web2 import static, http_headers, responsecode

class Foo(resource.Resource):
    def render(self, ctx):
        print ctx
        return http.Response(stream="Child/Foo")

class Sage(resource.Resource):
    def render(self, ctx):
        print ctx.args
        try:
          print 1
          #from sage.misc.all import sage_eval
          ans = str(eval(ctx.args['expr'][0]))
          print 2
        except Exception, msg:
          ans = msg
        return http.Response(stream=ans)

class Bar(resource.Resource):
    def render(self, ctx):
        print ctx
        return http.Response(stream="Bar")

class Child(resource.Resource):
    child_foo = Foo()

    def locateChild(self, request, segments):
        return self, ()

    def render(self, ctx):
        Toplevel.child_bar = Bar()
        return http.Response(stream="Child" + str(ctx.prepath))


class Toplevel(resource.Resource):
    addSlash = True
    child_child = Child()
    child_dir = static.File('/home/was/')
    child_sage = Sage()
    nb = None
    def render(self, ctx):
        print ctx
        print ctx.method
        print ctx.headers
        print ctx.args
        return http.Response(stream="Toplevel")


site = server.Site(Toplevel())



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

    Toplevel.nb = notebook.load_notebook(directory)

    def run(port):
        ## Create the config file
        config = open(conf, 'w')
        config.write("""
from sage.server.notebook.twist import site
from twisted.web2 import channel
from twisted.application import service, strports
application = service.Application("SAGE Notebook")
s = strports.service('tcp:%s', channel.HTTPFactory(site))
s.setServiceParent(application)
"""%(port))
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
