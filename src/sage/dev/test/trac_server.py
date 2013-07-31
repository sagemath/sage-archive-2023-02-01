class DoctestTracServer(object):
    def __init__(self):
        self.tickets = {}

class Ticket(object):
    def __init__(self, ticket, summary, description, attributes):
        attributes['summary'] = summary
        attributes['description'] = description
        self.id = ticket
        self.attributes = attributes
        self.time_created = 'not implemented'
        self.time_changed = 'not implemented'
        self.comments = []
