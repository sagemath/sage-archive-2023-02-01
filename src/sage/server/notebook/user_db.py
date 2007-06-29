"""nodoctest
"""

#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

class UserDatabase():
    def add_user(self, user, passwd, email):
        self[user] = UserRecord(user, passwd, email)

    def remove_user(self, user):
        del self[user]

class UserRecord:
    def __init__(self, username, passwd, email):
        self.username = username
        self.passwd = passwd
        self.email = email
