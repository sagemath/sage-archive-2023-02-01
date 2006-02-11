import bz2
import cPickle


def dump_to_string(obj, protocol=0):
    """
    Return a compressed binary string representation for obj from
    which obj can be reconstructed using the load_from_string
    function.
    """
    return bz2.compress(cPickle.dumps(obj, protocol))

def load_from_string(string):
    """
    Given a string as output by dump_to_string, returns an object.
    """
    return cPickle.loads(bz2.decompress(string))


def dump_to_file(obj, file, protocol=0):
    """
    Writes to the file with given name, a compressed binary string
    representation for obj from which obj can be reconstructed using
    the load_from_file function.  Only one object can be written to a
    file.
    INPUT:
        obj -- Python object
        file -- string (a filename)
    OUTPUT:
        creates a file
    """
    return open(file,'w').write(bz2.compress(cPickle.dumps(obj, protocol)))

def load_from_file(file):
    """
    Reads an object saved using dump_to_file.
    INPUT:
        file -- string (a filename)
    """
    return cPickle.loads(bz2.decompress(open(file).read()))

dumps = dump_to_string
loads = load_from_string
dump = dump_to_file
load = load_from_file

format_version=cPickle.format_version




