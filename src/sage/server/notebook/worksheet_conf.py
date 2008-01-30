"""nodoctest
"""

import conf

defaults = {'system':'sage',
            'prettyprint':False
           }

class WorksheetConfiguration(conf.Configuration):
    def defaults(self):
        return defaults
