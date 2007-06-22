import conf

defaults = {'max_history_length':500,
            'default_system':'sage',
            }

class UserConfiguration(conf.Configuration):
    def defaults(self):
        return defaults

