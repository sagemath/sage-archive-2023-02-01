class UserDatabase(UserDict):
    def add_user(self, user, passwd, email):
        self[user] = UserRecord(user, passwd, email)

    def remove_user(self, user):
        del self[user]

class UserRecord:
    def __init__(self, username, passwd, email):
        self.username = username
        self.passwd = passwd
        self.email = email
