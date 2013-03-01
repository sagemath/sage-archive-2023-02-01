The ``NameError`` raised on the second line should be displayed, even
if we crash immediately afterwards::

    sage: import time
    sage: import signal
    sage: this_gives_a_NameError
    sage: os.kill(os.getpid(), signal.SIGKILL)
