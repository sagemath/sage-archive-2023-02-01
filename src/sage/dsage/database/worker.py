class Worker(object):
    def __init__(self, host_info):
        for k, v in host_info.iteritems():
            setattr(self, k, v)
