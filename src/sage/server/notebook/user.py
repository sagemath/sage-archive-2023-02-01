"""nodoctest
"""

import copy
import crypt
import cPickle
import os

SALT = 'aa'

import user_conf
import twist

class User:
    def __init__(self, username, password='', email='', account_type='admin'):
        self.__username = username
        self.__password = crypt.crypt(password, SALT)
        self.__email = email
        if not account_type in ['admin', 'user', 'guest']:
            raise ValueError, "account type must be one of admin, user, or guest"
        self.__account_type = account_type
        self.__conf = user_conf.UserConfiguration()

    def __getstate__(self):
        d = copy.copy(self.__dict__)

        # Some old worksheets have this attribute, which we do *not* want to save.
        if d.has_key('history'):
            try:
                self.save_history()
                del d['history']
            except Exception, msg:
                print msg
                print "Unable to dump history of user %s to disk yet."%self.__username
        return d

    def history_list(self):
        try:
            return self.history
        except AttributeError:
            if twist.notebook is None: return []
            history_file = "%s/worksheets/%s/history.sobj"%(twist.notebook.directory(), self.__username)
            if os.path.exists(history_file):
                try:
                    self.history = cPickle.load(open(history_file))
                except:
                    print "Error loading history for user %s"%self.__username
            else:
                self.history = []
            return self.history

    def save_history(self):
        if not hasattr(self, 'history'):
            return
        if twist.notebook is None: return
        history_file = "%s/worksheets/%s/history.sobj"%(twist.notebook.directory(), self.__username)
        try:
            print "Dumping %s history to '%s'"%(self.__username, history_file)
            his = cPickle.dumps(self.history)
        except AttributeError:
            his = cPickle.dumps([])
        open(history_file,'w').write(his)

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
