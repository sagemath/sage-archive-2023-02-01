import conf

defaults = {}

class UserConfiguration(conf.Configuration):
    def __init__(self):
        conf.Configuration.__init__(self, defaults)


