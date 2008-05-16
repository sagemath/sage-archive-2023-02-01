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
        """
        sage: from sage.server.notebook.user import User
        sage: User('andrew', 'tEir&tiwk!', 'andrew@matrixstuff.com', 'user').username()
        'andrew'
        sage: User('sarah', 'Miaasc!', 'sarah@ellipticcurvesrock.org', 'user').username()
        'sarah'
        sage: User('bob', 'Aisfa!!', 'bob@pizzaisyummy.net', 'admin').username()
        'bob'
        """
        return self.__username

    def __repr__(self):
        return self.__username

    def conf(self):
        """
        sage: from sage.server.notebook.user import User
        sage: config = User('bob', 'Aisfa!!', 'bob@pizzaisyummy.net', 'admin').conf(); config
        Configuration: {}
        sage: config['max_history_length']
        500
        sage: config['default_system']
        'sage'
        sage: config['autosave_interval']
        180
        sage: config['default_pretty_print']
        False
        """
        return self.__conf

    def __getitem__(self, *args):
        return self.__conf.__getitem__(*args)

    def __setitem__(self, *args):
        self.__conf.__setitem__(*args)

    def password(self):
        """
        sage: from sage.server.notebook.user import User
        sage: user = User('bob', 'Aisfa!!', 'bob@pizzaisyummy.net', 'admin')
        sage: user.password()
        'aamxw5LCYcWY.'
        """
        return self.__password

    def set_password(self, password):
        """
        sage: from sage.server.notebook.user import User
        sage: user = User('bob', 'Aisfa!!', 'bob@pizzaisyummy.net', 'admin')
        sage: old = user.password()
        sage: user.set_password(self, 'Crrc!')
        sage: old != user.password()
        True
        """
        if password == '':
            self.__password = 'x'   # won't get as a password -- i.e., this account is closed.
        else:
            self.__password = crypt.crypt(password, SALT)

    def set_hashed_password(self, password):
        """
        sage: from sage.server.notebook.user import User
        sage: user = User('bob', 'Aisfa!!', 'bob@pizzaisyummy.net', 'admin')
        sage: user.set_hashed_password(self, 'Crrc!')
        sage: user.password()
        'Crrc!'
        """
        self.__password = password

    def get_email(self):
        """
        sage: from sage.server.notebook.user import User
        sage: user = User('bob', 'Aisfa!!', 'bob@pizzaisyummy.net', 'admin')
        sage: user.get_email()
        'bob@pizzaisyummy.net'
        """
        return self.__email

    def set_email(self, email):
        """
        sage: from sage.server.notebook.user import User
        sage: user = User('bob', 'Aisfa!!', 'bob@pizzaisyummy.net', 'admin').conf()
        sage: user.get_email()
        'bob@pizzaisyummy.net'
        sage: user.set_email('bob@ilovepizza.gov')
        sage: user.get_email()
        'bob@ilovepizza.gov'
        """
        self.__email = email

    def password_is(self, password):
        """
        sage: from sage.server.notebook.user import User
        sage: user = User('bob', 'Aisfa!!', 'bob@pizzaisyummy.net', 'admin')
        sage: user.password_is('ecc')
        False
        sage: user.password_is('Aisfa!!')
        True
        """
        if self.__username == "pub":
            return False
        return self.__password == crypt.crypt(password, SALT)

    def account_type(self):
        """
        sage: from sage.server.notebook.user import User
        sage: User('A', account_type='admin').account_type()
        'admin'
        sage: User('B', account_type='user').account_type()
        'user'
        sage: User('C', account_type='guest').account_type()
        'guest'
        """
        return self.__account_type

    def is_admin(self):
        """
        sage: from sage.server.notebook.user import User
        sage: User('A', account_type='admin').is_admin()
        True
        sage: User('B', account_type='user').is_admin()
        False
        """
        return self.__account_type == 'admin'

    def is_guest(self):
        """
        sage: from sage.server.notebook.user import User
        sage: User('A', account_type='guest').is_guest()
        True
        sage: User('B', account_type='user').is_guest()
        False
        """
        return self.__account_type == 'guest'
