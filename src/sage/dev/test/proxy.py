class Proxy(object):
    """
    A fake trac proxy for doctesting the functionality in this file which would require authentication by trac.

    EXAMPLES::

        sage: sage.dev.trac_interface.DoctestServerProxy(dev.trac)
        <sage.dev.trac_interface.DoctestServerProxy object at ...>
    """
    _ticket_data = dict()
    def __init__(self, trac):
        """
        Initialization.

        EXAMPLES::

            sage: type(sage.dev.trac_interface.DoctestServerProxy(dev.trac))
            <class 'sage.dev.trac_interface.DoctestServerProxy'>
        """
        self._trac = trac
        self._sshkeys = {}

    @property
    def sshkeys(self):
        """
        fake sshkeys trac plugin for doctest

        TESTS::

            sage: from sage.dev.trac_interface import DoctestServerProxy
            sage: proxy = DoctestServerProxy(dev.trac)
            sage: sshkeys = proxy.sshkeys
            sage: type(sshkeys)
            <class 'sage.dev.trac_interface.SshKeys'>
            sage: sshkeys.listusers()
            []
            sage: sshkeys.getkeys()
            []
            sage: sshkeys.setkeys(["foo", "bar"])
            0
            sage: sshkeys.getkeys()
            ['foo', 'bar']
        """
        try:
            return self._sshkeys_impl
        except AttributeError:
            pass
        class SshKeys(object):
            def _user(this):
                return self._trac._username
            def _setdefault(this):
                return self._sshkeys.setdefault(this._user(), set())
            def __setitem__(this, key, value):
                self._sshkeys[key] = value
            def setkeys(this, keys):
                this[this._user()] = this._setdefault().union(set(keys))
                return 0
            def getkeys(this):
                return list(this._setdefault())
            def listusers(this):
                return self._sshkeys.keys()

        self._sshkeys_impl = SshKeys()
        return self._sshkeys_impl

    @property
    def ticket(self):
        """
        fake ticket methods for trac xmlrpc plugin

        TESTS::

            sage: from sage.dev.trac_interface import DoctestServerProxy
            sage: proxy = DoctestServerProxy(dev.trac)
            sage: ticket = proxy.ticket
            sage: type(ticket)
            <class 'sage.dev.trac_interface.Ticket'>
            sage: ticket.create('a comment', 'a description', {}, False)
            14366
            sage: ticket.update(5614, 'a comment', {})
            Traceback (most recent call last):
            ...
            AssertionError
            sage: ticket.update(int(5614), 'a comment', {})
            (5614,)
        """
        class Ticket(object):
            def create(self, summary, description, attributes, notify):
                ticketnum = 14366
                DoctestServerProxy._ticket_data[ticketnum] = [ticketnum, summary, description, attributes, notify]
                return ticketnum
            def update(self, ticketnum, comment, attributes):
                assert isinstance(ticketnum, int)
                try:
                    DoctestServerProxy._ticket_data[ticketnum][3] = attributes
                except KeyError:
                    DoctestServerProxy._ticket_data[ticketnum] = [ticketnum, None, None, attributes, None]
                return (ticketnum,)
            def get(self, ticketnum):
                assert isinstance(ticketnum, int)
                print("returning attributes %s" %  DoctestServerProxy._ticket_data[ticketnum])
                return DoctestServerProxy._ticket_data[ticketnum]
                

        return Ticket()

