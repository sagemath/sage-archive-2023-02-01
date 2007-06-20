import user_conf

class User:
    def __init__(self, username, password='', email='', account_type='admin'):
        self.__username = username
        self.__password = password
        self.__email = email
        self.__account_type = account_type
        self.__conf = user_conf.UserConfiguration()

    def username(self):
        return self.__username

    def __repr__(self):
        return self.__username

    def conf(self):
        return self.__conf

    def password(self):
        return self.__password

    def account_type(self):
        return self.__account_type

    def is_admin(self):
        return self.__account_type == 'admin'
