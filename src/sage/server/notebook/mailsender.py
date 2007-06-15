from twisted.application import service
from twisted.application import internet
from twisted.internet import protocol
from twisted.mail import smtp, relaymanager
from StringIO import StringIO
from email.Generator import Generator

# Make an instance of this class if you need to send an email
class SMTPMessage:
    """
    A message to be sent off to an SMTP server.

    INPUT:
        mail_from   -- the SMTP MAIL FROM: ...
        rcpt        -- the SMTP RCPT TO: ...
        data        -- the entire message. This should be a StringIO object.
        identity    -- the SMTP HELO identity. For some SMTP servers, the
                       SAGE email sender's ip must resolve to this address.
        secret      -- Dunno?
    """
    def __init__(self, mail_from, rcpt, data, identity, secret=None):

        self._mail_from = mail_from
        self._rcpt = rcpt
        self._id = identity
        self._data = data
        self._secret = secret
        self._factory = SMTPClientFactory(self)
        # the server that we will deliver the mail to
        try:
            self._rcpt_domain = [(r.split('@'))[1] for r in rcpt]
        except ValueError:
            raise ValueError, "mal-formed recipient email address"

        print self._rcpt_domain

    def exchange_mail(self, _):
        smtp_client = internet.TCPClient(self._exchange, 25, self._factory)
        smtp_client.setServiceParent(self._app)

    def get_mx(self, host):
        def on_found_record(record):
            self._exchange = str(record.name)
        return relaymanager.MXCalculator().getMX(host).addCallback(on_found_record)

    def run_from(self, app):
        self._app = app
        self.get_mx(self._rcpt_domain[0]).addCallback(self.exchange_mail)

class MailClient(smtp.ESMTPClient):
    def __init__(self, mesg, **kwds):
        smtp.ESMTPClient.__init__(self,**kwds)
        self._mesg = mesg

    # these are all overridden from the super class, and called by Twisted to do
    # the real work
    getMailFrom = lambda self: self._mesg._mail_from
    getMailTo = lambda self: self._mesg._rcpt

    def getMailData(self):
        return StringIO(self._mesg._data)

    def sentMail(self, code, resp, numOk, addresses, log):
        print '<< dest SMTP server >> %s: %s' % (code, resp)
        print 'Sent %s messages.' % numOk
        from twisted.internet import reactor
        reactor.stop()

class SMTPClientFactory(protocol.ClientFactory):
    def __init__(self, mesg):
        self._mesg = mesg
        self._protocol = MailClient

    def buildProtocol(self, addr):
        mesg = self._mesg
        return self._protocol(mesg, secret=mesg._secret, identity=mesg._id)

application = service.Application("SAGE SMTP Client")
_from = "moretti@sage.math.washington.edu"
_to =  ["bobmoretti@gmail.com"]
_id = "sage.math.washington.edu"
data = "BLAH 2"

m = SMTPMessage(_from, _to, data, _id)
m.run_from(application)
