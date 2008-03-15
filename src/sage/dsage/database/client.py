class Client(object):
    """
    A client of the dsage server.

    """

    def __init__(self, username, public_key):
        """
        :type username: string
        :param username: username

        :type public_key: string
        :param public_key: public key of user

        """

        self.username = username
        self.public_key = public_key
        self.creation_time = None
        self.access_time = None
        self.last_login = None
        self.connected = None
        self.enabled = None

    def get_username(self):
        return self.username

    def get_public_key(self):
        return self.public_key

    def is_connected(self):
        return self.connected

    def is_enabled(self):
        return self.enabled