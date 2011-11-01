"""
Interface to Sloane On-Line Encyclopedia of Integer Sequences

To look up sequence A060843, type one of the following::

    sage: sloane_sequence(60843)       # optional - internet
    Searching Sloane's online database...
    [60843, 'Busy Beaver problem: maximal number of steps that an n-state Turing machine can make on an initially blank tape before eventually halting.', [1, 6, 21, 107]]

::

    sage: sloane_sequence("60843")     # optional - internet
    Searching Sloane's online database...
    [60843, 'Busy Beaver problem: maximal number of steps that an n-state Turing machine can make on an initially blank tape before eventually halting.', [1, 6, 21, 107]]

::

    sage: sloane_sequence("060843")    # optional - internet
    Searching Sloane's online database...
    [60843, 'Busy Beaver problem: maximal number of steps that an n-state Turing machine can make on an initially blank tape before eventually halting.', [1, 6, 21, 107]]

Do not prefix an integer with a 0 or it will be interpreted in
octal. Results are of the form [number, description, list], and
invalid numbers will cause ``sloane_sequence`` to
raise an ValueError exception::

    sage: sloane_sequence('sage')     # optional - internet
    Traceback (most recent call last):
    ...
    ValueError: sequence 'sage' not found

To look up the sequence "2, 3, 5, 7", simply put the numbers in a list.
The second argument specifies that at most 2 results will be returned.

::

    sage: sloane_find([2,3,5,7], 2)       # optional - internet
    Searching Sloane's online database...
    [[40, 'The prime numbers.', [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271]], [41, 'a(n) = number of partitions of n (the partition numbers).', [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627, 792, 1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, 12310, 14883, 17977, 21637, 26015, 31185, 37338, 44583, 53174, 63261, 75175, 89134]]]

To return no more than 3 results (default is 30), type

::

    sage: len(sloane_find([1,2,3,4,5], 3))      # optional - internet
    Searching Sloane's online database...
    3

Sequence A137443 includes Sage code, and a "b file" that was computed
with Sage::

    sage: sloane_find([7, 71, 281, 4523, 74713])  # optional - internet
    Searching Sloane's online database...
    [[137443, 'First n-digit prime in consecutive digits of e.', [7, 71, 281, 4523, 74713, 904523, 6028747, 72407663, 360287471, 7427466391, 75724709369, 749669676277, 8284590452353, 99959574966967, 724709369995957, 2470936999595749, 28459045235360287, 571382178525166427]]]

Note that the OEIS (http://oeis.org/classic/) claims to limit the number of
results to 100. Results are lists of the form [[number, description, list]],
and invalid input will cause sloane_find to return [].

In some cases, these functions may return [] even though the inputs
are legal. These cases correspond to errors from the OEIS server,
and calling the functions again may fix the problem.

Alternatively, the SloaneEncyclopedia object provides access to a
local copy of the database containing only the sequences and their
names. To use this you must download and install the database using
``SloaneEncyclopedia.install()``, or
``SloaneEncyclopedia.install_from_gz()`` if you have already
downloaded the database manually.

To look up a sequence, type

::

    sage: SloaneEncyclopedia[60843]               # optional - sloane_database
    [1, 6, 21, 107]

To get the name of a sequence, type

::

    sage: SloaneEncyclopedia.sequence_name(1)     # optional - sloane_database
    'Number of groups of order n.'

To search locally for a particular subsequence, type

::

    sage: SloaneEncyclopedia.find([1,2,3,4,5], 1)    # optional - sloane_database
    [(15, [1, 2, 3, 4, 5, 7, 7, 8, 9, 11, 11, 13, 13, 16, 16, 16, 17, 19, 19, 23, 23, 23, 23, 25, 25, 27, 27, 29, 29, 31, 31, 32, 37, 37, 37, 37, 37, 41, 41, 41, 41, 43, 43, 47, 47, 47, 47, 49, 49, 53, 53, 53, 53, 59, 59, 59, 59, 59, 59, 61, 61, 64, 64, 64, 67, 67, 67, 71, 71, 71, 71, 73])]

The default maximum number of results is 30, but to return up to
100, type

::

    sage: SloaneEncyclopedia.find([1,2,3,4,5], 100)    # optional - sloane_database
    [(15, [1, 2, 3, 4, 5, 7, 7, 8, 9, 11, 11, ...

Results in either case are of the form [ (number, list) ].

TODO:

- When this program gets a sloane sequence from the database it
  actually downloads a huge amount of information about it, then
  throws most of it away. Also, it returns the data to the user as a
  very simple tuple. It would be much better to return an instance of
  a class::

      class SloaneSequence: ...

  and the class should have methods for each of the things that
  Sloane records about a sequence. Also, when possible, it should be
  able to compute more terms.

AUTHORS:

- Steven Sivek (2005-12-22): first version

- Steven Sivek (2006-02-07): updated to correctly handle the new
  search form on the Sloane website, and it's now also smarter about
  loading the local database in that it doesn't convert a sequence
  from string form to a list of integers until absolutely necessary.
  This seems to cut the loading time roughly in half.

- Steven Sivek (2009-12-22): added the SloaneEncyclopedia functions
  install() and install_from_gz() so users can get the latest versions
  of the OEIS without having to get an updated spkg; added
  sequence_name() to return the description of a sequence; and changed
  the data type for elements of each sequence from int to Integer.

"""

#*****************************************************************************
#
#      Sage: Copyright (C) 2005-2006 William Stein <wstein@gmail.com>
#                               and  Steven Sivek  <ssivek@mit.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import bz2, os, re, urllib

from sage.misc.all import verbose
import sage.rings.integer_ring
ZZ = sage.rings.integer_ring.IntegerRing()

class SloaneEncyclopediaClass:
    """
    A local copy of the Sloane Online Encyclopedia of Integer Sequences
    that contains only the sequence numbers and the sequences
    themselves.
    """
    def __init__(self):
        """
        Initialize the database but do not load any of the data.
        """
        self.__path__ = "%s/data/sloane/"%os.environ["SAGE_ROOT"]
        self.__file__ = "%ssloane-oeis.bz2"%self.__path__
        self.__file_names__ = "%ssloane-names.bz2"%self.__path__
        self.__loaded__ = False
        self.__loaded_names__ = False

    def __repr__(self):
        """
        String representation of this database. OUTPUT: str
        """
        return "Sloane Online Encyclopedia of Integer Sequences"

    def __iter__(self):
        """
        Returns an iterator through the encyclopedia. Elements are of the
        form [number, sequence].
        """
        for i in self.__data__:
            yield [i, self[i]]

    def __getitem__(self, N):
        """
        Return sequence N in the encyclopedia. If sequence N does not
        exist, return [].

        INPUT:


        -  ``N`` - int


        OUTPUT: list
        """
        self.load()
        if not N in self.__data__: # sequence N does not exist
            return []
        if self.__data__[N][1] is None: # list N has not been created yet
            list = self.__data__[N][2].strip(',').split(',')
            self.__data__[N][1] = [ZZ(n) for n in list]
        return self.__data__[N][1]

    def __len__(self):
        """
        Return the number of sequences in the encyclopedia.
        """
        self.load()
        return len(self.__data__)


    def find(self, seq, maxresults=30):
        """
        Return a list of all sequences which have seq as a subsequence, up
        to maxresults results. Sequences are returned in the form (number,
        list).

        INPUT:


        -  ``seq`` - list

        -  ``maxresults`` - int


        OUTPUT: list of 2-tuples (i, v), where v is a sequence with seq as
        a subsequence.
        """
        self.load()

        answer, nanswer = [], 0
        pattern = re.sub(r'[\[\]]', ',', str(seq).replace(' ',''))
        for i in self.__data__:
            if self.__data__[i][2].find(pattern) != -1:
                answer.append((i, self[i]))
                nanswer = nanswer + 1
                if nanswer == maxresults:
                    return answer

        return answer

    def install(self, oeis_url="http://oeis.org/classic/stripped.gz", names_url="http://oeis.org/classic/names.gz", overwrite=False):
        """
        Download and install the online encyclopedia, raising an IOError if
        either step fails.

        INPUT:

        - ``oeis_url`` - string (default: "http://www.research.att.com...")
          The URL of the stripped.gz encyclopedia file.

        - ``names_url`` - string (default: "http://www.research.att.com...")
          The URL of the names.gz encyclopedia file.  If you do not want to
          download this file, set names_url=None.

        - ``overwrite`` - boolean (default: False) If the encyclopedia is
          already installed and overwrite=True, download and install the latest
          version over the installed one.
        """
        ### See if the encyclopedia already exists
        if not overwrite and os.path.exists(self.__file__):
            raise IOError, "Sloane encyclopedia is already installed"

        tm = verbose("Downloading stripped version of Sloane encyclopedia")
        try:
            fname, _ = urllib.urlretrieve(oeis_url);
        except IOError, msg:
            raise IOError, "%s\nError fetching the following website:\n    %s\nTry checking your internet connection."%(msg, oeis_url)

        if not names_url is None:
            try:
                nname, _ = urllib.urlretrieve(names_url);
            except IOError, msg:
                raise IOError, "%s\nError fetching the following website:\n    %s\nTry checking your internet connection."%(msg, names_url)
        else:
            nname = None
        verbose("Finished downloading", tm)

        self.install_from_gz(fname, nname, overwrite)
        # Delete the temporary downloaded files
        os.remove(fname)
        if not nname is None:
            os.remove(nname)

    def install_from_gz(self, stripped_file, names_file, overwrite=False):
        """
        Install the online encyclopedia from a local stripped.gz file.

        INPUT:

        - ``stripped_file`` - string. The name of the stripped.gz OEIS file.

        - ``names_file`` - string.  The name of the names.gz OEIS file, or
          None if the user does not want it installed.

        - ``overwrite`` - boolean (default: False) If the encyclopedia is
          already installed and overwrite=True, install 'filename' over the
          old encyclopedia.
        """
        if not overwrite and os.path.exists(self.__file__):
            raise IOError, "Sloane encyclopedia is already installed"

        copy_gz_file(stripped_file, self.__file__)

        if not names_file is None:
            copy_gz_file(names_file, self.__file_names__)
        else:
            # Delete old copies of names.gz since their sequence numbers
            # probably won't match the newly installed stripped.gz
            if os.path.exists(self.__file_names__):
                os.remove(self.__file_names__)

        # Remove the old database from memory so the new one will be
        # automatically loaded next time the user tries to access it
        self.unload()


    def load(self):
        """
        Load the entire encyclopedia into memory from a file. This is done
        automatically if the user tries to perform a lookup or a search.
        """
        if self.__loaded__ == True:
            return
        try:
            file_seq = bz2.BZ2File(self.__file__, 'r')
        except IOError:
            raise IOError, "The Sloane Encyclopedia database must be installed.  Use e.g. 'SloaneEncyclopedia.install()' to download and install it."

        self.__data__ = {}

        tm = verbose("Loading Sloane encyclopedia from disk")
        entry = re.compile(r'A(?P<num>\d{6}) ,(?P<body>.*),$');
        for L in file_seq:
            if len(L) == 0:
                continue
            m = entry.search(L)
            if m:
                seqnum = int(m.group('num'))
                msg = m.group('body').strip()
                self.__data__[seqnum] = [seqnum, None, ','+msg+',', None]
        file_seq.close()

        try:
            file_names = bz2.BZ2File(self.__file_names__, 'r')
            entry = re.compile(r'A(?P<num>\d{6}) (?P<body>.*)$');
            for L in file_names:
                if len(L) == 0: continue
                m = entry.search(L)
                if m:
                    seqnum = int(m.group('num'))
                    self.__data__[seqnum][3] = m.group('body').strip()
            file_names.close()
            self.__loaded_names__ = True
        except KeyError:
            ### Some sequence in the names file isn't in the database
            raise KeyError, "Sloane OEIS sequence and name files do not match.  Try reinstalling, e.g. SloaneEncyclopedia.install(overwrite=True)."
        except IOError, msg:
            ### The names database is not installed
            self.__loaded_names__ = False

        verbose("Finished loading", tm)
        self.__loaded__ = True

    def sequence_name(self, N):
        """
        Return the name of sequence N in the encyclopedia. If sequence N
        does not exist, return ''.  If the names database is not installed,
        raise an IOError.

        INPUT:

        -  ``N`` - int

        OUTPUT: string

        EXAMPLES:

        sage: SloaneEncyclopedia.sequence_name(1) # optional - sloane_database
        'Number of groups of order n.'
        """
        self.load()
        if not self.__loaded_names__:
            raise IOError, "The Sloane OEIS names file is not installed.  Try reinstalling, e.g. SloaneEncyclopedia.install(overwrite=True)."

        if not N in self.__data__: # sequence N does not exist
            return ''
        return self.__data__[N][3]

    def unload(self):
        """
        Remove the database from memory.
        """
        if self.__loaded__ == False:
            return
        del self.__data__
        self.__loaded__ = False
        self.__loaded_names__ = False

SloaneEncyclopedia = SloaneEncyclopediaClass()

def copy_gz_file(gz_source, bz_destination):
    """
    Decompress a gzipped file and install the bzipped verson.  This is
    used by SloaneEncyclopedia.install_from_gz to install several
    gzipped OEIS database files.

    INPUT:

    - ``gz_source`` - string. The name of the gzipped file.

    - ``bz_destination`` - string.  The name of the newly compressed file.
    """
    import gzip
    from sage.misc.misc import sage_makedirs

    # Read the gzipped input
    try:
        gz_input = gzip.open(gz_source, 'r')
        db_text = gz_input.read()
        gz_input.close()
    except IOError, msg:
        raise IOError, "Error reading gzipped input file:\n%s"%msg

    # Write the bzipped output
    try:
        sage_makedirs(os.path.dirname(bz_destination))
        bz2_output = bz2.BZ2File(bz_destination, 'w')
        bz2_output.write(db_text)
        bz2_output.close()
    except IOError, msg:
        raise IOError, "Error writing bzipped output file:\n%s"%msg


def parse_sequence(text):
    entry = re.compile(r'%(?P<letter>[A-Za-z]) A(?P<num>\d{6}) (?P<body>.*)$')
    unsigned, signed, list = [], [], []
    description = ''
    seqnum = -1

    # Fix broken lines: the next line is indented.
    text = text.replace('\n               ', '');

    for line in re.split(r'[\s,]*\n', text):
        m = entry.search(line)
        if m:
            seqnum = ZZ(m.group('num').lstrip('0'))
            type = m.group('letter')
            msg = m.group('body').lstrip().rstrip()
            if type == 'S' or type == 'T' or type == 'U':
                unsigned += msg.split(',')
            elif type == 'V' or type == 'W' or type == 'X':
                signed += msg.split(',')
            elif type == 'N':
                description = msg

    if signed == []:
        list = unsigned
    else:
        list = signed
    return [seqnum, description, [ZZ(n) for n in list]]

def sloane_sequence(number, verbose = True):
    """
    Returns a list with the number, name, and values for the sequence
    ``number`` in Sloane's online database of integer
    sequences.

    INPUT:

    - ``number`` (integer) corresponding to the sequence.

    -  ``verbose`` - (boolean) print a string to let the
       user know that it is working and not hanging.

       Set to ``True`` by default.


    EXAMPLES::

        sage: sloane_sequence(22) # optional - internet
        Searching Sloane's online database...
        [22,
         'Number of centered hydrocarbons with n atoms.',
         [0,
          1,
          0,
          1,
          ...
          36201693122]]

    The input must not be a sequence itself::

        sage: sloane_sequence(prime_range(100))
        Traceback (most recent call last):
        ...
        TypeError: input must be an integer or string that specifies the id of the Sloane sequence to download
    """
    if not isinstance(number, str):
        try:
            number = str(ZZ(number))
        except TypeError:
            raise TypeError, "input must be an integer or string that specifies the id of the Sloane sequence to download"
    results = sloane_find('id:A%s'%number, verbose = verbose)
    if len(results) == 0:
        raise ValueError, "sequence '%s' not found"%number
    return results[0]

def sloane_find(list, nresults = 30, verbose=True):
    """
    Searches Sloane's Online Encyclopedia of Integer Sequences for a
    sequence containing the number provided in ``list``.

    INPUT:


    -  ``list`` - (list) a list of integers to search
       Sloane's for

    -  ``nresults`` - (integer) the maximum number of
       results to return default: 30

    -  ``verbose`` - (boolean) print a string to let the
       user know that it is working and not hanging. default: True


    OUTPUT: A list of matches in Sloane's database. Each match consists
    of a list of the sequence number, the name of the sequence, and
    some initial terms of the sequence.

    EXAMPLES::

        sage: sloane_find([1,1,2,3,5,8,13,21], nresults=1) #optional - internet
        Searching Sloane's online database...
        [[45,
          'Fibonacci numbers: F(n) = F(n-1) + F(n-2), F(0) = 0, F(1) = 1, F(2) = 1, ...',
          [0,
           1,
           1,
           2,
           3,
           5,
           8,
           13,
           21,
           ...
           39088169]]]
    """
    liststr = re.sub(r'[\[\] ]', '', str(list))
    urlparams = urllib.urlencode({'q': liststr,
                                  'p': 1,
                                  'n': nresults,
                                  'fmt': 2,
                                  'sort': 0});


    try:
        if verbose:
            print "Searching Sloane's online database..."
        url = "http://oeis.org/classic/"
        f = urllib.urlopen(url+'?'+urlparams);
        s = f.read()
        f.close()
    except IOError, msg:
        raise IOError, "%s\nError fetching the following website:\n    %s\nTry checking your internet connection."%(msg, url)

    t = s.lower()
    i = t.find("<pre>")
    j = t.find("</pre>")
    if i == -1 or j == -1:
        #raise IOError, "Error parsing data (missing pre tags)."
        return []
    text = s[i+5:j].strip()

    results = []
    for line in text.split('\n\n'):
        if line.find('%') != -1: # ignore text at top of results
            results.append(parse_sequence(line))

    return results
