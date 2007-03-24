# Twisted imports
from twisted.application import internet, service

from twisted.web import script, static, server, vhost, resource, util

from twisted.internet import threads, reactor

from twisted.internet import ssl

# Enable threads
from twisted.python import threadable
threadable.init(1)

class WikiResource(resource.Resource):
    """ Wiki resource """
    isLeaf = 1

    def render(self, request):
        return server.NOT_DONE_YET


class WikiRoot(resource.Resource):
    """ Wiki root resource """

    def getChild(self, name, request):
        # Serve images and css from '/wiki'
        if request.prepath == [] and name == 'wiki':
            return resource.Resource.getChild(self, name, request)

        # Serve special 'root' files from '/wiki'
        elif name in ['favicon.ico', 'robots.txt'] and request.postpath == []:
            return self.children['wiki'].getChild(name, request)

        # All other through moin

        # TODO: fix profile code to include the request init and ignore
        # first request. I'm not doing this now since its better done
        # with the new twisted code waiting in fix branch. --Nir
        else:
            if config.memoryProfile:
                config.memoryProfile.addRequest()
            req = RequestTwisted(request, name, reactor,
                                 properties=config.properties)
            if config.hotshotProfile:
                threads.deferToThread(config.hotshotProfile.runcall, req.run)
            else:
                threads.deferToThread(req.run)
            return WikiResource()


class MoinRequest(server.Request):
    """ MoinMoin request

    Enable passing of file-upload filenames
    """

    def requestReceived(self, command, path, version):
        """ Called by channel when all data has been received.

        Override server.Request method for POST requests, to fix file
        uploads issue.
        """
        if command == 'POST':
            self.requestReceivedPOST(path, version)
        else:
            server.Request.requestReceived(self, command, path, version)

    def requestReceivedPOST(self, path, version):
        """ Handle POST requests

        This is a modified copy of server.Request.requestRecived,
        modified to use cgi.FieldStorage to handle file uploads
        correctly.

        Creates an extra member extended_args which also has
        filenames of file uploads ( FIELDNAME__filename__ ).
        """
        import cgi

        self.content.seek(0,0)
        self.args = {}
        self.extended_args = {}
        self.stack = []

        self.method = 'POST'
        self.uri = path
        self.clientproto = version
        x = self.uri.split('?')

        argstring = ""
        if len(x) == 1:
            self.path = self.uri
        else:
            if len(x) != 2:
                from twisted.python import log
                log.msg("May ignore parts of this invalid URI: %s"
                        % repr(self.uri))
            self.path, argstring = x[0], x[1]

        # cache the client and server information, we'll need this later to be
        # serialized and sent with the request so CGIs will work remotely
        self.client = self.channel.transport.getPeer()
        self.host = self.channel.transport.getHost()

        # create dummy env for cgi.FieldStorage
        env = {
            'REQUEST_METHOD': self.method,
            'QUERY_STRING': argstring,
            }
        form = cgi.FieldStorage(fp=self.content,
                                environ=env,
                                headers=self.received_headers)

        # Argument processing

        args = self.args
        try:
            keys = form.keys()
        except TypeError:
            pass
        else:
            for key in keys:
                values = form[key]
                if not isinstance(values, list):
                    values = [values]
                fixedResult = []
                for i in values:
                    if isinstance(i, cgi.MiniFieldStorage):
                        fixedResult.append(i.value)
                    elif isinstance(i, cgi.FieldStorage):
                        fixedResult.append(i.value)
                        # multiple uploads to same form field are stupid!
                        if i.filename:
                            args[key + '__filename__'] = i.filename
                args[key] = fixedResult

        self.process()


class MoinSite(server.Site):
    """ Moin site """
    requestFactory = MoinRequest

    def startFactory(self):
        """ Setup before starting """
        # Memory profile
        if config.memoryProfile:
            config.memoryProfile.sample()

        # hotshot profile
        if config.hotshotProfile:
            import hotshot
            config.hotshotProfile = hotshot.Profile(config.hotshotProfile)
        server.Site.startFactory(self)

    def stopFactory(self):
        """ Cleaup before stoping """
        server.Site.stopFactory(self)
        if config.hotshotProfile:
            config.hotshotProfile.close()


class TwistedConfig(Config):
    """ Twisted server default config """

    name = 'mointwisted'
    properties = {}
    docs = '/usr/share/moin/htdocs'
    user = 'www-data'
    group = 'www-data'
    port = 8080
    interfaces = ['']
    threads = 10
    timeout = 15*60 # 15 minutes
    logPath = None
    virtualHosts = None
    memoryProfile = None
    hotshotProfile = None

    # sslcert = ('/whereever/cert/sitekey.pem', '/whereever/cert/sitecert.pem')
    sslcert = None

    def __init__(self):
        Config.__init__(self)

        # Check for '' in interfaces, then ignore other
        if '' in self.interfaces:
            self.interfaces = ['']


