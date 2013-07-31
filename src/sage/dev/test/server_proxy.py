class DoctestServerProxy(object):
    def __init__(self, server):
        self._server = server

        self.ticket = DoctestTicketProxy(self)

    def _check_authentication(self, privilege):
        import xmlrpclib
        raise xmlrpclib.Fault(403, "%s privileges are required to perform this operation. You don't have the required permissions."%privilege)

class AuthenticatedDoctestServerProxy(DoctestServerProxy):
    def __init__(self, server, username, password):
        DoctestServerProxy.__init__(self, server)

        self._username = username
        self._password = password

    def _check_authentication(self, privilege):
        pass

class DoctestTicketProxy(object):
    def __init__(self, server_proxy):
        self._server_proxy = server_proxy

    def create(self, summary, description, attributes):
        self._server_proxy._check_authentication("TICKET_CREATE")

        from trac_server import Ticket
        ticket = len(self._server_proxy._server.tickets)+1
        self._server_proxy._server.tickets[ticket] = Ticket(ticket, summary, description, attributes)
        return ticket

    def update(self, ticket, comment, attributes):
        self._server_proxy._check_authentication("TICKET_MODIFY")

        ticket = self._server_proxy._server.tickets[ticket]
        ticket.comments.append(comment)
        ticket.attributes = attributes

    def get(self, ticket):
        ticket = self._server_proxy._server.tickets[ticket]
        return [ticket.id, ticket.time_created, ticket.time_changed, ticket.attributes]
