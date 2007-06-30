"""nodoctest
"""

import conf

defaults = {'system':'sage',
           }

class WorksheetConfiguration(conf.Configuration):
    def defaults(self):
        return defaults
