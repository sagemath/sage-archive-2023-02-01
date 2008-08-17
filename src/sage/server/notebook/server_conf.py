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
