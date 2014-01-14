Interrupt the doctester (the parent process)::

    sage: import signal
    sage: import time
    sage: os.kill(os.getppid(), signal.SIGINT)
    sage: time.sleep(10)
    sage: os._exit(0)
