import conf

defaults = {'cell_input_color':'#0000000',
            'cell_output_color':'#0000EE',
            'word_wrap_cols':80,
            'max_history_length':500,
            'number_of_backups':3,

            'idle_timeout':1800,
            'idle_check_interval':30,

            'save_interval':15,
           }

class ServerConfiguration(conf.Configuration):
    def defaults(self):
        return defaults
