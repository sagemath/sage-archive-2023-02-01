r"""
Huffman Encoding
"""
from string import join

def frequency_table(string):
    r"""
    Return the frequency table corresponding to the given
    string.

    INPUT:

    - ``string`` -- a string

    EXAMPLE::

        sage: from sage.coding.source_coding.huffman import frequency_table
        sage: str = "Sage is my most favorite general purpose computer algebra system"
        sage: frequency_table(str)
        {'a': 5, ' ': 9, 'c': 1, 'b': 1, 'e': 8, 'g': 3, 'f': 1, 'i': 2, 'm': 4, 's': 5, 'o': 4, 'n': 1, 'p': 3, 'S': 1, 'r': 5, 'u': 2, 't': 4, 'v': 1, 'y': 2, 'l': 2}

    """
    d = {}
    for l in string:
        d[l] = d.get(l,0) + 1

    return d


class Huffman():
    r"""
    Huffman Encoding

    This class implements the basic functionalities
    of Huffman's encoding.

    It can build a Huffman code from a given string, or
    from the information of a dictionary associating
    to each key (the elements of the alphabet) a weight
    (most of the time, a probability value or a number
    of occurrences). For example ::

        sage: from sage.coding.source_coding.huffman import Huffman, frequency_table
        sage: h1 = Huffman("There once was a french fry")
        sage: for letter, code in h1.encoding_table().iteritems():
        ...       print "'"+ letter + "' : " + code
        'a' : 0111
        ' ' : 00
        'c' : 1010
        'e' : 100
        'f' : 1011
        'h' : 1100
        'o' : 11100
        'n' : 1101
        's' : 11101
        'r' : 010
        'T' : 11110
        'w' : 11111
        'y' : 0110

    We could have obtained the same result by "training" the Huffman
    code on the following table of frequency ::

        sage: ft = frequency_table("There once was a french fry"); ft
        {'a': 2, ' ': 5, 'c': 2, 'e': 4, 'f': 2, 'h': 2, 'o': 1, 'n': 2, 's': 1, 'r': 3, 'T': 1, 'w': 1, 'y': 1}
        sage: h2 = Huffman(frequencies = ft)

    Once ``h1`` has been trained, and hence possesses an encoding code,
    it is possible to obtain the Huffman encoding of any string
    (possibly the same) using this code::

        sage: encoded = h1.encode("There once was a french fry"); encoded
        '11110110010001010000111001101101010000111110111111010001110010110101001101101011000010110100110'

    Which can be decoded the following way::

        sage: h1.decode(encoded)
        'There once was a french fry'

    Obviously, if we try to decode a string using a Huffman instance which
    has been trained on a different sample (and hence has a different encoding
    table), we are likely to get some random-looking string ::

        sage: h3 = Huffman("There once were two french fries")
        sage: h3.decode(encoded)
        ' wehnefetrhft ne ewrowrirTc'

    ... precisely what we deserved :-)

    INPUT:

    One among the following:

    - ``string`` -- a string from which the Huffman encoding should
      be created

    - ``frequencies`` -- a dictionary associating its frequency or
     its number of occurrences to each letter of the alphabet.

    """

    def __init__(self, string = None, frequencies = None):
        r"""
        Constructor for Huffman

        INPUT:

        One among the following:

        - ``string`` -- a string from which the Huffman encoding should
          be created

        - ``frequencies`` -- a dictionary associating its frequency or
         its number of occurrences to each letter of the alphabet.

        EXAMPLE::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)

        If both arguments are supplied, an exception is raised ::

            sage: Huffman(string=str, frequencies={'a':8})
            Traceback (most recent call last):
            ...
            ValueError: Exactly one of `string` or `frequencies` parameters must be defined

        """

        self._character_to_code = []

        if sum([string is not None, frequencies is not None]) != 1:
            raise ValueError("Exactly one of `string` or `frequencies` parameters must be defined")

        if string is not None:
            self._build_code(frequency_table(string))
        elif frequencies is not None:
            self._build_code(frequencies)

    def _build_code_from_tree(self, tree, d, prefix=''):
        r"""
        Builds the code corresponding to a given tree and prefix

        INPUT:

        - ``tree`` -- integer, or list of size `2`

        - ``d`` -- the dictionary to fill

        - ``prefix`` (string) -- binary string which is the prefix
          of any element of the tree

        EXAMPLE::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: d = {}
            sage: h._build_code_from_tree(h._tree, d)

        """
        try:
            self._build_code_from_tree(tree[0], d, prefix=prefix+'0')
            self._build_code_from_tree(tree[1], d, prefix=prefix+'1')
        except TypeError:
            d[tree] = prefix

    def _build_code(self, dic):
        r"""
        Returns a Huffman code for each one of the given elements.

        INPUT:

        - ``dic`` (dictionary) -- associates to each letter of the alphabet
          a frequency or a number of occurrences.

        EXAMPLE::

            sage: from sage.coding.source_coding.huffman import Huffman, frequency_table
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: d = {}
            sage: h._build_code(frequency_table(str))
        """

        from heapq import heappush, heappop

        index = dic.items()
        heap = []

        for i,(e,w) in enumerate(index):
            heappush(heap, (w, i) )

        while len(heap)>=2:
            (w1, i1) = heappop(heap)
            (w2, i2) = heappop(heap)
            heappush(heap, (w1+w2,[i1,i2]))


        d = {}
        self._tree = heap[0][1]
        self._build_code_from_tree(self._tree, d)
        self._index = dict([(i,e) for i,(e,w) in enumerate(index)])
        self._character_to_code = dict([(e,d[i]) for i,(e,w) in enumerate(index)])


    def encode(self, string):
        r"""
        Returns an encoding of the given string based
        on the current encoding table

        INPUT:

        - ``string`` (string)

        EXAMPLE:

        This is how a string is encoded then decoded ::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: encoded = h.encode(str); encoded
            '00000110100010101011000011101010011100101010011011011100111101110010110100001011011111000001110101010001010110011010111111011001110100101000111110010011011100101011100000110001100101000101110101111101110110011000101011000111111101101111010010111001110100011'
            sage: h.decode(encoded)
            'Sage is my most favorite general purpose computer algebra system'

        """
        if self._character_to_code:
            return join(map(lambda x:self._character_to_code[x],string), '')


    def decode(self, string):
        r"""
        Returns a decoded version of the given string
        corresponding to the current encoding table.

        INPUT:

        - ``string`` (string)


        EXAMPLE:

        This is how a string is encoded then decoded ::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: encoded = h.encode(str); encoded
            '00000110100010101011000011101010011100101010011011011100111101110010110100001011011111000001110101010001010110011010111111011001110100101000111110010011011100101011100000110001100101000101110101111101110110011000101011000111111101101111010010111001110100011'
            sage: h.decode(encoded)
            'Sage is my most favorite general purpose computer algebra system'

        Of course, the string one tries to decode has to be a binary one. If
        not, an exception is raised ::

            sage: h.decode('I clearly am not a binary string')
            Traceback (most recent call last):
            ...
            ValueError: The given string does not only contain 0 and 1
        """
        chars = []
        tree = self._tree
        index = self._index
        for i in string:

            if i == '0':
                tree = tree[0]
            elif i == '1':
                tree = tree[1]
            else:
                raise ValueError('The given string does not only contain 0 and 1')

            if not isinstance(tree,list):
                chars.append(index[tree])
                tree = self._tree

        return join(chars, '')

    def encoding_table(self):
        r"""
        Returns the current encoding table

        OUTPUT:

        A dictionary associating its code to each trained letter of
        the alphabet

        EXAMPLE::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: h.encoding_table()
            {'S': '00000', 'a': '1101', ' ': '101', 'c': '110000', 'b': '110001', 'e': '010', 'g': '0001', 'f': '110010', 'i': '10000', 'm': '0011', 'l': '10011', 'o': '0110', 'n': '110011', 'p': '0010', 's': '1110', 'r': '1111', 'u': '10001', 't': '0111', 'v': '00001', 'y': '10010'}
        """
        return self._character_to_code.copy()

    def tree(self):
        r"""
        Returns the Huffman tree corresponding to the current encoding

        OUTPUT:

        A tree

        EXAMPLE::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: T = h.tree(); T
            Digraph on 39 vertices
            sage: T.show(figsize=[20,20])
        """

        from sage.graphs.digraph import DiGraph
        g = DiGraph()
        g.add_edges(self._generate_edges(self._tree))
        return g

    def _generate_edges(self, tree, father='', id=''):
        if father=='':
            u = 'root'
        else:
            u = father
        try:
            return self._generate_edges(tree[0], father=father+id, id='0') + \
                self._generate_edges(tree[1], father=father+id, id='1') + \
                ([(u, father+id)] if (father+id) != '' else [])

        except TypeError:
            return [(u, self.decode(father+id)+' : '+(father+id))]
