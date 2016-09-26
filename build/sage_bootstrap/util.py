import time


def is_url(url):
    """
    Test whether argument is url
    """
    url = url.rstrip()
    if len(url.splitlines()) > 1:
        return False
    if url.find(' ') >= 0:
        return False
    return (
        url.startswith('http://') or
        url.startswith('https://') or
        url.startswith('ftp://')
    )


def retry(func, exc=Exception, retries=2, delay=1):
    """
    Call ``func()``, and if an exception is raised retry it up to ``retries``
    times, with a ``delay`` in seconds between each retry.

    To be clear, the function is called up to ``retries + 1`` times.

    If given, ``exc`` can be either an `Exception` or a tuple of `Exception`s
    in which only those exceptions result in a retry, and all other exceptions
    are raised.
    """

    while retries >= 0:
        try:
            func()
        except exc:
            retries -= 1
            time.sleep(delay)
        else:
            break
