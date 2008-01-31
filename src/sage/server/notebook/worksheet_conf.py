"""nodoctest
"""

import conf

defaults = {'system':'sage',
            'pretty_print':False
           }

class WorksheetConfiguration(conf.Configuration):
    def defaults(self):
        return defaults
