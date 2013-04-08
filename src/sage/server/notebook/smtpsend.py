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

"""
Sending mail using Twisted

AUTHOR:
    -- Bobby Moretti
"""

from twisted.mail import smtp, relaymanager
from email.MIMEBase import MIMEBase
from email.MIMEMultipart import MIMEMultipart
import sys

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

def send_mail(fromaddr, toaddr, subject, body, on_success=sendComplete, on_failure=handleError):
    try:
        recpt_domain = toaddr.split('@')[1]
    except (ValueError, IndexError):
        raise ValueError, "mal-formed destination address"
    message = buildMessage(fromaddr, toaddr, subject, body)
    messageData = message.as_string(unixfrom=False)

    def on_found_record(mx_rec):
        smtp_server = str(mx_rec.name)
        sending = smtp.sendmail(smtp_server, fromaddr, [toaddr], messageData)
        sending.addCallback(on_success).addErrback(on_failure)

    relaymanager.MXCalculator().getMX(recpt_domain).addCallback(on_found_record)

