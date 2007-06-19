import conf

defaults = {'max_history_length':500,
            }

class UserConfiguration(conf.Configuration):
    def defaults(self):
        return defaults

