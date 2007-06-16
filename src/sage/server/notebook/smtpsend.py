from twisted.protocols import smtp
from twisted.internet import reactor, defer
from twisted.mail import relaymanager
from email.MIMEBase import MIMEBase
from email.MIMEMultipart import MIMEMultipart
from email import Encoders
import sys, mimetypes, os

def buildMessage(fromaddr, toaddr, subject, body):
    message = MIMEMultipart()
    message['From'] = fromaddr
    message['To'] = toaddr
    message['Subject'] = subject
    textPart = MIMEBase('text', 'plain')
    textPart.set_payload(body)
    message.attach(textPart)
    return message

def sendComplete(result):
    print "Message sent."

def handleError(error):
    print >> sys.stderr, "Error", error.getErrorMessage()

def send_mail(self, fromaddr, toaddr, subject, body):
    try:
        recpt_domain = toaddr.split('@')[1]
    except ValueError:
        raise ValueError, "mal-formed destination address"
    message = buildMessage(fromaddr, toaddr, subject, body)
    messageData = message.as_string(unixfrom=False)

    def on_found_record(mx_rec):
        smtp_server = str(mx_rec.name)
        sending = smtp.sendmail(smtp_server, fromaddr, [toaddr], messageData)
        sending.addCallback(sendComplete).addErrback(handleError)

    relaymanager.MXCalculator().getMX(recpt_domain).addCallback(on_found_record)
