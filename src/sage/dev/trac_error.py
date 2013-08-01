class TracError(RuntimeError):
    pass

class TracConnectionError(TracError):
    def __init__(self):
        TracError.__init__(self, "Connection to trac server failed.")

class TracInternalError(TracError):
    def __init__(self, fault):
        self._fault = fault
        self.faultCode = fault.faultCode

    def __str__(self):
        return str(self._fault)

