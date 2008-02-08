"""nodoctest
"""

import conf

defaults = {'max_history_length':500,
            'default_system':'sage',
            'autosave_interval':3*60,   # (in seconds)
            'default_pretty_print': False
            }

class UserConfiguration(conf.Configuration):
    def defaults(self):
        return defaults

