# This file is part of the OLD Sage notebook and is NOT actively developed,
# maintained, or supported.  As of Sage v4.1.2, all notebook development has
# moved to the separate Sage Notebook project:
#
# http://nb.sagemath.org/
#
# The new notebook is installed in Sage as an spkg (e.g., sagenb-0.3.spkg).
#
# Please visit the project's home page for more information, including directions on
# obtaining the latest source code.  For notebook-related development and support,
# please consult the sage-notebook discussion group:
#
# http://groups.google.com/group/sage-notebook

"""nodoctest
"""

#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

from twisted.application import service
from twisted.application import internet
from twisted.internet import protocol
from twisted.mail import smtp, relaymanager
from StringIO import StringIO
from email.MIMEBase import MIMEBase
from email.MIMEMultipart import MIMEMultipart


class MailMessage(MIMEMultipart):
    """
    Represents an email's data.
    """
    def __init__(self, fromaddr, toaddr, subject, body):
        MIMEMultipart.__init__(self)
        self['From'] = fromaddr
        self['To'] = toaddr
        self['Subject'] = subject
        text_part = MIMEBase('text', 'plain')
        text_part.set_payload(body)
        self.attach(text_part)

# Make an instance of this class if you need to send an email
class SMTPInput:
    """
    A message to be sent off to an SMTP server.

    INPUT:
        mail_from   -- the SMTP MAIL FROM: ...
        rcpt        -- the SMTP RCPT TO: ...
        data        -- the entire message. This should be a StringIO object.
        identity    -- the SMTP HELO identity. For some SMTP servers, the
                       Sage email sender's IP must resolve to this address.
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

    def exchange_mail(self, exchange):
        smtp_client = internet.TCPClient(exchange, 25, self._factory)
        smtp_client.setServiceParent(self._app)

    def get_mx(self, host):
        on_found_record = lambda rec: str(rec.name)
        # return a deffered that will call on_found_record when it's done.
        # on_found_record's return value gets passed to exchange_mail
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
_from = "moretti@u.math.washington.edu"
_to =  ["wstein@gmail.com"]
_id = "sage.math.washington.edu"
subject = "Progress"
body = \
"""
William,

I'm sending this to you from Twisted. I think I'm ready to try to plug this in
from SAGE, as soon as Yi has the login stuff working (which he almost does.)

~Bobby
"""
data = MailMessage(_from, _to[0], subject, body)
m = SMTPInput(_from, _to, data.as_string(unixfrom=False), _id)
m.run_from(application)
