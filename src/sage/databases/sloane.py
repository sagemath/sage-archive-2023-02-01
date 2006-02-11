"""nodoctest  -- it takes too long!
Interface to Sloane On-Line Encyclopedia of Integer Sequences

AUTHOR: Steven Sivek (ssivek@mit.edu), 2005-12-22

To look up sequence A060843, type one of the following:

    sage: sloane_sequence(60843)
    [60843, 'Busy Beaver problem: maximal number of steps that an n-state Turing machine can make on an initially blank tape before eventually halting.', [1, 6, 21, 107]]
    sage: sloane_sequence("60843")
    [60843, 'Busy Beaver problem: maximal number of steps that an n-state Turing machine can make on an initially blank tape before eventually halting.', [1, 6, 21, 107]]
    sage: sloane_sequence("060843")
    [60843, 'Busy Beaver problem: maximal number of steps that an n-state Turing machine can make on an initially blank tape before eventually halting.', [1, 6, 21, 107]]

Do not prefix an integer with a 0 or it will be interpreted in octal.
Results are of the form [number, description, list], and invalid
numbers will cause \code{sloane_sequence} to raise an ValueError exception:

    sage: sloane_sequence('sage')
    []

To look up the sequence
    sage: sloane_find([2,3,5,7])
    [[53874, 'Triangle T(n,k) = number of Boolean functions mapping {0,1}^n to {0,1}^n with image of size k (k = 0...2^n) under action of GL(n,2).', [1, 2, 1, 1, 2, 2, 2, 1, 1, 2, 2, 3, 4, 3, 2, 2, 1, 1, 2, 2, 3, 5, 7, 9, 11, 12, 11, 9, 7, 5, 3, 2, 2, 1, 1, 2, 2, 3, 5, 8, 14, 23, 35, 55, 84, 117, 158, 204, 242, 274, 290, 274, 242, 204, 158, 117, 84, 55, 35, 23, 14, 8, 5, 3, 2, 2, 1, 1, 2, 2, 3, 5, 8, 15, 29, 54, 107, 227, 495, 1131]], [58398, 'Partition triangle A008284 read from right to left.', ...

To return no more than 2 results (default is 30), type
    sage: sloane_find([1,2,3,4,5], 2)

Note that the OEIS (???: todo) claims to limit the number of results
to 100.  Results are lists of the form [ [number, description, list]
], and invalid input will cause sloane_find to return [].

In some cases, these functions may return [] even though the inputs are legal.
These cases correspond to errors from the OEIS server, and calling the
functions again may fix the problem.


Alternatively, the SloaneEncyclopedia object provides access to a local copy
of the database containing only the sequences.   To use this you must install
the optional database_sloane_oeis-2005-12 package using
\code{sage -i database_sloane_oeis-2005-12}.

To look up a sequence, type
    sage: SloaneEncyclopedia[60843]
    [1, 6, 21, 107]

To search locally for a particular subsequence, type
    sage: SloaneEncyclopedia.find([1,2,3,4,5])

The default maximum number of results is 30, but to return up to 200, type
    sage: SloaneEncyclopedia.find([1,2,3,4,5], 200)

Results in either case are of the form [ (number, list) ].
"""

#*****************************************************************************
#
#      SAGE: Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#                          and  Steven Sivek  <ssivek@mit.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import bz2, os, re, urllib

from sage.misc.all import verbose

class SloaneEncyclopediaClass:
    """
    A local copy of the Sloane Online Encyclopedia of Integer Sequences
    that contains only the sequence numbers and the sequences themselves.
    """
    def __init__(self):
        """
        Initialize the database but do not load any of the data.
        """
        self.__file__ = "%s/data/sloane/sloane-oeis.bz2"%os.environ["SAGE_ROOT"]
        self.__arraysize__ = 114751 # maximum sequence number + 1
        self.__loaded__ = False

    def __repr__(self):
        """
        String representation of this database.
        OUTPUT:
           str
        """
        return "Sloane Online Encyclopedia of Integer Sequences"

    def __iter__(self):
        """
        Returns an iterator through the encyclopedia.  Elements are of the
        form [number, sequence].
        """
        for i in range(1, self.__arraysize__):
            if self[i] != []:
                yield [i, self[i]]

    def __getitem__(self, N):
        """
        Return sequence N in the encyclopedia.  If sequence N does not
        exist, return [].

        INPUT:
            N -- int
        OUTPUT:
            list
        """
        self.load()
        if self.__data__[N] == None:
            return []
        return self.__data__[N][1]

    def __len__(self):
        """
        Return the number of sequences in the encyclopedia.
        """
        self.load()
        return len(self.__data__)


    def find(self, seq, maxresults=30):
        """
        Return a list of all sequences which have seq as a subsequence,
        up to maxresults results.  Sequences are returned in the form
        [number, list].

        INPUT:
            seq -- list
            maxresults -- int
        OUTPUT:
            list of 2-tuples (i, v), where v is a sequence with seq as a subsequence.
        """
        self.load()

        answer, nanswer = [], 0
        pattern = re.sub(r'[\[\]]', ',', str(seq).replace(' ',''))
        for i in range(1, self.__arraysize__):
            if self.__data__[i] != None:
                if self.__data__[i][2].find(pattern) != -1:
                    answer.append((i, self[i]))
                    nanswer = nanswer + 1
                    if nanswer == maxresults:
                        return answer

        return answer

    def load(self):
        """
        Load the entire encyclopedia into memory from a file.  This is done
        automatically if the user tries to perform a lookup or a search.
        """
        if self.__loaded__ == True:
            return
        try:
            file = bz2.BZ2File(self.__file__, 'r')
        except IOError:
            raise IOError, "The Sloane Encyclopedia optional SAGE package must be installed.  Usage e.g., 'sage -i database_sloane_oeis-2005-12' to download and install it."

        entry = re.compile(r'A(?P<num>\d{6}) ,(?P<body>.*),$');
        self.__data__ = [None] * self.__arraysize__
        tm = verbose("Loading Sloane encyclopedia from disk")
        for L in file:
            if len(L) == 0:
                continue
            m = entry.search(L)
            if m:
                seqnum = int(m.group('num'));
                msg = m.group('body').strip();
                list = [int(n) for n in msg.split(',')]
                self.__data__[seqnum] = [seqnum, list, ','+msg+',']
        verbose("Finished loading", tm)
        self.__loaded__ = True

    def unload(self):
        """
        Remove the database from memory.
        """
        if self.__loaded__ == False:
            return
        del self.__data__
        self.__loaded__ = False

SloaneEncyclopedia = SloaneEncyclopediaClass()

def parse_sequence(text):
    entry = re.compile(r'%(?P<letter>[A-Za-z]) A(?P<num>\d{6}) (?P<body>.*)$')
    unsigned, signed, list = [], [], []
    description = ''
    seqnum = -1

    # Fix broken lines: they end with '\' and the next line is indented.
    text = text.replace('\\\n  ', '');

    for line in re.split(r'[\s,]*\n', text):
        m = entry.search(line)
        if m:
            seqnum = int(m.group('num'))
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
    return [seqnum, description, [int(n) for n in list]]

def sloane_sequence(number):
    try:
        print "Looking up in Sloane's online database (this requires a net connection and is slow)..."
        f = urllib.urlopen("http://www.research.att.com/cgi-bin/access.cgi/as/njas/sequences/eisA2.cgi?Anum=A%s"%number)
        s = f.read()
        f.close()
    except IOError:
        s = ''

    i = s.find("<PRE>")
    j = s.find("</PRE>")
    if i == -1 or j == -1:
        raise IOError, "Error parsing data (missing pre tags)."
    text = s[i+5:j].strip()

    return parse_sequence(text)

def sloane_find(list, nresults = 30):
    liststr = re.sub(r'[\[\] ]', '', str(list))
    print liststr
    urlparams = urllib.urlencode({'sequence': liststr,
                                  'choice': 1,
                                  'maxhit': nresults,
                                  'noold': 2});

    try:
        print "Searching Sloane's online database (this requires a net connection and is slow)..."
        f = urllib.urlopen("http://www.research.att.com/cgi-bin/access.cgi/as/njas/sequences/eismum2.cgi", urlparams);
        s = f.read()
        f.close()
    except IOError:
        s = ''

    i = s.find("<PRE>")
    j = s.find("</PRE>")
    if i == -1 or j == -1:
        raise IOError, "Error parsing data (missing pre tags)."
    text = s[i+5:j].strip()

    results = []
    for line in text.split('\n\n'):
        if line.find('%') != -1: # ignore text at top of results
            results.append(parse_sequence(line))

    return results
