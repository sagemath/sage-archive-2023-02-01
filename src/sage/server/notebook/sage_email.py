"""
Send an Email

Sage supports very easily sending an email from Sage to notify
yourself when some event occurs.  This does not require configuring an
email server or anything else, since Sage already includes by default
a sophisticated email server (which is part of Twisted).

EXAMPLES::

    sage: email('xxxsageuser@gmail.com', 'The calculation finished!')  # not tested
    Child process ... is sending email to xxxsageuser@gmail.com

AUTHOR:
    -- William Stein (2008-12-13)
"""

#############################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


import os

def default_email_address():
    """
    Get the hostname and username of the user running this
    program. This provides a default email from address.

    OUTPUT:
        string

    EXAMPLES::

        sage: sage.server.notebook.sage_email.default_email_address()
        '...@...'
    """
    import socket
    hostname = socket.gethostname()
    username = os.popen('whoami').read().strip()
    return '%s@%s'%(username, hostname)

def email(to, subject, body = '', from_address = None, verbose = True, block = False, kill_on_exit = False):
    """
    Send an email message.

    INPUT:
        to           -- string; address of recipient
        subject      -- string; subject of the email
        body         -- string (default: ''); body of the email
        from_address -- string (default: username@hostname); address
                        email will appear to be from
        verbose      -- whether to print status information when the email is sent
        block        -- bool (default: False); if True this function doesn't
                        return until the email is actually sent.  if
                        False, the email gets sent in a background
                        thread.
        kill_on_exit -- bool (default: False): if True, guarantee that
                        the sending mail subprocess is killed when you
                        exit sage, even if it failed to send the
                        message.  If False, then the subprocess might
                        keep running for a while.  This should never
                        be a problem, but might be useful for certain
                        users.

    EXAMPLES::

        sage: email('xxxsageuser@gmail.com', 'The calculation finished!')  # not tested
        Child process ... is sending email to xxxsageuser@gmail.com

    NOTE: This function does not require configuring an email server
          or anything else at all, since Sage already includes by
          default a sophisticated email server (which is part of
          Twisted).
    """

    # We use Fork to make this work, because we have to start the
    # Twisted reactor in order to use it's powerful sendmail
    # capabilities.  Unfortunately, Twisted is purposely designed so
    # that its reactors cannot be restarted.  Thus if we don't fork,
    # one could send at most one email.  Of course, forking means this
    # won't work on native Windows.  It might be possible to get this
    # to work using threads instead, but I did not do so, since Python
    # threading with Twisted is not fun, and would likely have many
    # of the same problems.  Plus the below works extremely well.

    try:
        pid = os.fork()
    except:
        print "Fork not possible -- the email command is not supported on this platform."
        return

    if from_address is None:
        # Use a default email address as the from: line.
        from_address = default_email_address()

    if pid: # We're the parent process
        if kill_on_exit:
            # Tell the Sage cleaner about this subprocess, just in case somehow it fails
            # to properly quit (e.g., smtp is taking a long time), so it will get killed
            # no matter what when sage exits.  Zombies are bad bad bad, no matter what!
            import sage.interfaces.cleaner
            sage.interfaces.cleaner.cleaner(pid)  # register pid of forked process with cleaner
        if verbose:
            print "Child process %s is sending email to %s..."%(pid,to)
        # Now wait for the fake subprocess to finish.
        os.waitpid(pid,0)
        return

    if not block:
        # Do a non-block sendmail, which is typically what a user wants, since it can take
        # a while to send an email.

        # Use the old "double fork" trick -- otherwise there would *definitely* be a zombie
        # every time.  Here's a description from the web of this trick:
        # "If you can't stand zombies, you can get rid of them with a double fork().
        #  The forked child immediately forks again while its parent calls waitpid().
        #  The first forked process exits, and the parent's waitpid() returns, right
        #  away.  That leaves an orphaned process whose parent reverts to 1 ("init")."
        pid = os.fork()
        if pid:
            # OK, we're in the subprocess of the subprocess -- we
            # again register the subprocess we just spawned with the
            # zombie cleaner just in case, then we kill ourself, as
            # explained above.
            if kill_on_exit:
                import sage.interfaces.cleaner
                sage.interfaces.cleaner.cleaner(pid)   # register with cleaner
            os.kill(os.getpid(),9)                 # suicide

    # Now we're the child process.  Let's do stuff with Twisetd!
    from smtpsend import send_mail, reactor

    # First define two callback functions.  Each one optionally prints
    # some information, then kills the subprocess dead.
    def on_success(result):
        """
        Callback in case of a successfully sent email.
        """
        if verbose:
            print "Successfully sent an email to %s."%to
        reactor.stop()
        os.kill(os.getpid(),9)                     # suicide

    def on_failure(error):
        """
        Callback in case of a failure sending an email.
        """
        if verbose:
            print "Failed to send email to %s."%to
            print "-"*70
            print error.getErrorMessage()
            print "-"*70
        reactor.stop()
        os.kill(os.getpid(),9)                    # suicide

    # Finally, call the send_mail function.  This is code that sets up
    # a twisted deferred, which actually happens when we run the
    # reactor.
    send_mail(from_address, to, subject, body, on_success, on_failure)

    # Start the twisted reactor.
    reactor.run()
