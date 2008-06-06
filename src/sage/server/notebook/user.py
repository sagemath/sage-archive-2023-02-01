"""nodoctest
"""

import crypt

SALT = 'aa'

import user_conf

class User:
    def __init__(self, username, password='', email='', account_type='admin'):
        self.__username = username
        self.__password = crypt.crypt(password, SALT)
        self.__email = email
        if not account_type in ['admin', 'user', 'guest']:
            raise ValueError, "account type must be one of admin, user, or guest"
        self.__account_type = account_type
        self.__conf = user_conf.UserConfiguration()

    def username(self):
        return self.__username

    def __repr__(self):
        return self.__username

    def conf(self):
        return self.__conf

    def __getitem__(self, *args):
        return self.__conf.__getitem__(*args)

    def __setitem__(self, *args):
        self.__conf.__setitem__(*args)

    def password(self):
        return self.__password

    def set_password(self, password):
        if password == '':
            self.__password = 'x'   # won't get as a password -- i.e., this account is closed.
        else:
            self.__password = crypt.crypt(password, SALT)

    def set_hashed_password(self, password):
        self.__password = password

    def set_email(self, email):
        self.__email = email

    def password_is(self, password):
        if self.__username == "pub":
            return False
        return self.__password == crypt.crypt(password, SALT)

    def account_type(self):
        return self.__account_type

    def is_admin(self):
        return self.__account_type == 'admin'

    def is_guest(self):
        return self.__account_type == 'guest'
