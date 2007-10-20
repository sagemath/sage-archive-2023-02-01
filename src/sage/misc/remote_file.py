import os, sys

def get_remote_file(filename, verbose=True):
    """
    INPUT:
        filename -- the URL of a file on the web, e.g.,
             "http://modular.math.washington.edu/myfile.txt"
        verbose -- whether to display download status

    OUTPUT:
        creates a file in the temp directory and returns the
        absolute path to that file.

    EXAMPLES:
        sage: g = get_remote_file("http://sage.math.washington.edu/home/novoselt/lp/K3vertices", verbose=False)   # optional -- requires an internet connection
        sage: L = lattice_polytope.read_all_polytopes(g, "K3 %4d")                    # optional
        sage: L[0]                                                                    # optional
        K3    0: 3-dimensional, 4 vertices.
    """
    if verbose:
        print "Attempting to load remote file: " + filename
    import misc

    temp_name = misc.tmp_filename() + '.' + os.path.splitext(filename)[1][1:]
    # IMPORTANT -- urllib takes a long time to load,
    # so do not import it in the module scope.
    import urllib
    global cur
    cur = 0
    if verbose:
        sys.stdout.write("Loading: [")
        sys.stdout.flush()
        urllib.urlretrieve(filename, temp_name, report_hook)
        print "]"
    else:
        urllib.urlretrieve(filename, temp_name)
    return temp_name

cur = 0
def report_hook(block, size, total):
     global cur
     n = block*size*50/total
     if n > cur:
          cur = n
          sys.stdout.write('.')
          sys.stdout.flush()

