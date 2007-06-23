import conf

defaults = {'max_history_length':500,
            'default_system':'sage',
            'autosave_interval':3*60,   # (in seconds)
            }

class UserConfiguration(conf.Configuration):
    def defaults(self):
        return defaults

