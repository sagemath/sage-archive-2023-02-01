class TracError(RuntimeError):
    pass

class TracConnectionError(TracError):
    pass

class TracInternalError(TracError):
    def __init__(self, fault):
        self._fault = fault
        self.faultCode = fault.faultCode
