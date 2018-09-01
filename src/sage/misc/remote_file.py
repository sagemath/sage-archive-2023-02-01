"get_remote_file"

from __future__ import absolute_import

import os
import sys


def get_remote_file(filename, verbose=True):
    """
    INPUT:

    - ``filename`` -- the URL of a file on the web, e.g.,
      ``"http://modular.math.washington.edu/myfile.txt"``

    - ``verbose`` -- whether to display download status

    OUTPUT:

    creates a file in the temp directory and returns the absolute path
    to that file.

    EXAMPLES::

        sage: g = get_remote_file("http://sagemath.org/ack.html", verbose=False)   # optional - internet
        sage: len(open(g).read())   # optional - internet; random
        10198
    """
    if verbose:
        print("Attempting to load remote file: " + filename)

    from sage.misc.temporary_file import tmp_filename
    temp_name = tmp_filename() + '.' + os.path.splitext(filename)[1][1:]
    # IMPORTANT -- urllib takes a long time to load,
    # so do not import it in the module scope.

    # import compatible with py2 and py3
    from six.moves.urllib.request import urlretrieve, Request, urlopen
    req = Request(filename, headers={"User-Agent":"sage-doctest"})

    global cur
    cur = 0
    if verbose:
        sys.stdout.write("Loading: [")
        sys.stdout.flush()
        #urlretrieve(filename, temp_name, report_hook)
        content = urlopen(req, timeout=1)
        with open(temp_name, 'w') as f:
            f.write(content.read())
        print("]")
    else:
        content = urlopen(req, timeout=1)
        with open(temp_name, 'w') as f:
            f.write(content.read())
    return temp_name

cur = 0
def report_hook(block, size, total):
     global cur
     n = block*size*50/total
     if n > cur:
          cur = n
          sys.stdout.write('.')
          sys.stdout.flush()
