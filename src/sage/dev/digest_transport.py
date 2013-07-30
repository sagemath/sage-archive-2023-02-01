r"""
HTTP transport to the trac server

AUTHORS:

- TODO: add authors from github's history and trac's history

"""
#*****************************************************************************
#       Copyright (C) 2013 TODO
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

class DigestTransport(object, SafeTransport):
    """
    Handles an HTTP transaction to an XML-RPC server.

    EXAMPLES::

        sage: sage.dev.digest_transport.DigestTransport()
        <sage.dev.digest_transport.DigestTransport object at ...>

    """
    def __init__(self, **kwds):
        """
        Initialization.

        EXAMPLES::

            sage: type(sage.dev.digest_transport.DigestTransport())
            <class 'sage.dev.digest_transport.DigestTransport'>
            sage: type(sage.dev.digest_transport.DigestTransport(realm='realm', url='url', username='username', password='password'))
            <class 'sage.dev.digest_transport.DigestTransport'>

        """
        auth = (kwds[k] if k in kwds else None for k in ('realm', 'url', 'username', 'password'))

        SafeTransport.__init__(self, **kwds)

        authhandler = urllib2.HTTPDigestAuthHandler()
        if all(x is not None for x in auth):
            authhandler.add_password(*auth)

        self.opener = urllib2.build_opener(authhandler)

    def single_request(self, host, handler, request_body, verbose):
        """
        Issue an XML-RPC request.

        EXAMPLES::

            sage: from sage.env import TRAC_SERVER_URI
            sage: import urlparse
            sage: url = urlparse.urlparse(TRAC_SERVER_URI).netloc
            sage: d = sage.dev.digest_transport.DigestTransport()
            sage: d.single_request(url, 'xmlrpc',         # optional: internet
            ....: '''<?xml version='1.0'?>
            ....: <methodCall>
            ....: <methodName>ticket.get</methodName>
            ....: <params>
            ....: <param>
            ....: <value><int>1000</int></value>
            ....: </param>
            ....: </params>
            ....: </methodCall>
            ....: ''', 0)
            ([1000,
              <DateTime '20071025T16:48:05' at ...>,
              <DateTime '20080110T08:28:40' at ...>,
              {'status': 'closed',
               'changetime': <DateTime '20080110T08:28:40' at ...>,
               'description': '',
               'reporter': 'was',
               'cc': '',
               'type': 'defect',
               'milestone': 'sage-2.10',
               '_ts': '1199953720000000',
               'component': 'distribution',
               'summary': 'Sage does not have 10000 users yet.',
               'priority': 'major',
               'owner': 'was',
               'time': <DateTime '20071025T16:48:05' at ...>,
               'keywords': '',
               'resolution': 'fixed'}],)

        """
        req = urllib2.Request(
                urlparse.urlunparse(('http', host, handler, '', '', '')),
                request_body, {'Content-Type': 'text/xml',
                    'User-Agent': self.user_agent})

        response = self.opener.open(req)

        self.verbose = verbose
        return self.parse_response(response)

