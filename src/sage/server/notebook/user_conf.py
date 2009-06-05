"""nodoctest
"""

import conf

defaults = {'max_history_length':1000,
            'default_system':'sage',
            'autosave_interval':60*60,   # 1 hour in seconds
            'default_pretty_print': False
            }

class UserConfiguration(conf.Configuration):
    def defaults(self):
        return defaults

