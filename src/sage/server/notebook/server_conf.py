# This file is part of the OLD Sage notebook and is NOT actively developed,
# maintained, or supported.  As of Sage v4.1.2, all notebook development has
# moved to the separate Sage Notebook project:
#
# http://nb.sagemath.org/
#
# The new notebook is installed in Sage as an spkg (e.g., sagenb-0.3.spkg).
#
# Please visit the project's home page for more information, including directions on
# obtaining the latest source code.  For notebook-related development and support,
# please consult the sage-notebook discussion group:
#
# http://groups.google.com/group/sage-notebook

"""nodoctest
"""

import conf

defaults = {'cell_input_color':'#0000000',
            'cell_output_color':'#0000EE',
            'word_wrap_cols':72,
            'max_history_length':250,
            'number_of_backups':3,

            'idle_timeout':120,        # 2 minutes
            'idle_check_interval':360,

            'save_interval':360,        # seconds

            'doc_pool_size':128,
            'email':False
           }

class ServerConfiguration(conf.Configuration):
    def defaults(self):
        return defaults
