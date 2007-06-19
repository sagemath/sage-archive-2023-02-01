import conf

defaults = {'cell_input_color':'#0000000',
            'cell_output_color':'#0000EE',
            'word_wrap_cols':80,
            'max_history_length':500,
            'number_of_backups':3
           }

class ServerConfiguration(conf.Configuration):
    def __init__(self):
        conf.Configuration.__init__(self, defaults)
