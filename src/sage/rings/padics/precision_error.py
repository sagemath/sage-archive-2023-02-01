class PrecisionError(Exception):
    pass

class HaltingError(PrecisionError):
    pass

class PrecisionLimitError(PrecisionError):
    pass
