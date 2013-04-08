#####################################################################
# Copyright (C)  2007 Alex Clemesha <clemesha@gmail.com>
#                 and Dorian Raymer <deldotdr@gmail.com>
#####################################################################

from knoboo import config
# init of config can take options
kbconfig = config.Config()
CONFIG = kbconfig.read_config()

#indexer for searching notebooks
from knoboo import indexer
import sys

#knoboo.external.twisted.web2
from knoboo.external.twisted.web2 import log, server, channel
from twisted.application import service, strports
from twisted.cred import portal, checkers, credentials
#from twisted.python import log

#for periodic indexing
from twisted.application import internet

#twisted database acccess
from twisted.enterprise import adbapi

# access control / security
from knoboo.authority import guard, check, avatars

def knoboo(port, secure, proxy, appname="Knoboo"):
    print "port = ", port
    dbargs = {"database":CONFIG['db_main_path']}
    dbConnection = adbapi.ConnectionPool(CONFIG['db_driver'], **dbargs)

    realm = avatars.LoginSystem(dbConnection)
    #create portal and attach checkers
    p = portal.Portal(realm)
    p.registerChecker(check.HashedPasswordDataBaseChecker(dbConnection))
    p.registerChecker(checkers.AllowAnonymousAccess(), credentials.IAnonymous)

    # wrap the resources to make them guarded
    rsrc = guard.SessionWrapper(p)

    #rsrc = log.LogWrapperResource(rsrc) #XXX This enable verbose http logging.
    #log.DefaultCommonAccessLoggingObserver().start()

    # build website resources
    site = server.Site(rsrc)
    factory = channel.HTTPFactory(site)

    #create the application
    application = service.Application(appname)

    #create the service
    # If proxy, serve only tcp on local interface
    if int(proxy):
        # backend server behind apache or nginx
        s = 'tcp:%s:interface=127.0.0.1' % str(port)
        srv = strports.service(s, factory)
        srv.setServiceParent(application)
    else: # Not proxy; standalone server open to the outside
        if int(secure):
            from twisted.internet import ssl
            # serve ssl, secure webserver
            s = 'ssl:%s:privateKey=ssl/privkey.pem:certKey=ssl/cacert.pem' % str(port)
            srv = strports.service(s, factory)
            srv.setServiceParent(application)
        else:
            # serve tcp, webserver
            s = 'tcp:%s' % str(port)
            srv = strports.service(s, factory)
            srv.setServiceParent(application)
    return application


application = knoboo(CONFIG['port'], CONFIG['secure'], CONFIG['proxy'])

if sys.version_info > (2, 4):
    class LoopingIndexer(object):
        def __init__(self):
            indexer.indexer()
            self.lastest = indexer.latest()

        def update_index(self):
            indexer.update(self.lastest)
            self.lastest = indexer.latest()

    li = LoopingIndexer()
    loop = internet.TimerService(10, li.update_index)
    loop.setServiceParent(application)

