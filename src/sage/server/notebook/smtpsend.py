from twisted.protocols import smtp
from twisted.internet import reactor, defer
from email.MIMEBase import MIMEBase
from email.MIMEMultipart import MIMEMultipart
from email import Encoders
import sys, mimetypes, os

def buildMessage(fromaddr, toaddr, subject, body, filenames):
    message = MIMEMultipart()
    message['From'] = fromaddr
    message['To'] = toaddr
    message['Subject'] = subject
    textPart = MIMEBase('text', 'plain')
    textPart.set_payload(body)
    message.attach(textPart)
    for filename in filenames:
        # guess the mimetype
        mimetype = mimetypes.guess_type(filename)[0]
        if not mimetype: mimetype = 'application/octet-stream'
        maintype, subtype = mimetype.split('/')
        attachment = MIMEBase(maintype, subtype)
        attachment.set_payload(file(filename).read())
        # base64 encode for safety
        Encoders.encode_base64(attachment)
        # include filename info
        attachment.add_header('Content-Disposition', 'attachment',
                              filename=os.path.split(filename)[1])
        message.attach(attachment)
    return message

def sendComplete(result):
    print "Message sent."
    # reactor.stop()

def handleError(error):
    print >> sys.stderr, "Error", error.getErrorMessage()
    reactor.stop()

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: %s smtphost fromaddr toaddr file1 [file2, ...]" % (
            sys.argv[0])
        sys.exit(1)

    smtphost = sys.argv[1]
    fromaddr = sys.argv[2]
    toaddr = sys.argv[3]
    filenames = sys.argv[4:]
    subject = raw_input("Subject: ")
    body = raw_input("Message (one line): ")
    message = buildMessage(fromaddr, toaddr, subject, body, filenames)
    messageData = message.as_string(unixfrom=False)
    sending = smtp.sendmail(smtphost, fromaddr, [toaddr], messageData)
    sending.addCallback(sendComplete).addErrback(handleError)
    reactor.run()

