"""
TODO:

* check which import statements are needed
* include tests (a la OEIS)
* include full text search
* truncate input by size - how could we do that for distribution search?
* get code for a composition of maps
* for statistics: do we want to "call" them on an object, or rather have a getitem method?
  currently, it's call:

            sage: f=findstat(2); f
            St000002: The number of occurrences of the pattern $[1,2,3]$ inside a permutation of length at least 3.

            sage: f(Permutation([1,2,4,3,5]))        # optional -- internet
            7

"""
from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.databases.oeis import FancyTuple

from string import join
from ast import literal_eval
from collections import defaultdict
from urllib import urlopen, urlencode
import re
import webbrowser

# Combinatoral collections
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix
from sage.combinat.binary_tree import BinaryTree
from sage.combinat.core import Core
from sage.combinat.dyck_word import DyckWord
from sage.combinat.root_system.cartan_type import CartanType_abstract
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern
from sage.graphs.graph import Graph
from sage.combinat.composition import Composition
from sage.combinat.partition import Partition
from sage.combinat.ordered_tree import OrderedTree
from sage.combinat.parking_functions import ParkingFunction, ParkingFunction_class
from sage.combinat.perfect_matching import PerfectMatching
from sage.combinat.permutation import Permutation
from sage.combinat.posets.posets import Poset, FinitePoset
from sage.combinat.tableau import SemistandardTableau, StandardTableau
from sage.combinat.set_partition import SetPartition

findstat_url = 'http://www.findstat.org/'
findstat_url_result = findstat_url + "StatisticFinder/Result/"
findstat_url_submit = findstat_url + 'StatisticsDatabase/NewStatistic/'
findstat_url_edit = findstat_url + 'StatisticsDatabase/EditStatistic/'
findstat_url_downloads = 'http://downloads.findstat.org/'
# the number of values FindStat allows to search for at most
findstat_max_values = 200

def _fetch(url):
    r"""
    Fetch the given ``url``.

    INPUT:

    - ``url`` - a string corresponding to the URL to be fetched.

    OUTPUT:

    - a string representing the fetched web page.

    TESTS::

        sage: from sage.databases.findstat import _fetch, findstat_url
        sage: _fetch(findstat_url + 'hints.html')[-8:-1]            # optional -- internet
        '</html>'
    """
    try:
        _ = verbose("Fetching URL %s ..." %url, caller_name='FindStat')
        f = urlopen(url)
        result = f.read()
        f.close()
        return result
    except IOError as msg:
        raise IOError("%s\nError fetching %s." % (msg, url))

class FindStat:
    r"""
    The Combinatorial Statistic Finder.

    ``FindStat`` is a class representing the FindStat database.  Its
    behaviour is modeled after the ``OEIS`` class.  You can query it
    using its methods, but ``FindStat`` can also be called directly
    with three arguments:

    - ``query`` - it can be:

      - a string representing a FindStat ID (e.g. 'St000045').

      - an integer representing a FindStat ID (e.g. 45).

      - a dictionary from sage objects to integer values 
        representing the statistic on a combinatorial
        collection, and optionally a string naming the collection.

      - a list of pairs of the form (object, value) or (list of
        objects, list of values), where in the latter case, there
        should be as many values as objects, and optionally a string
        naming the collection

      - a single pair of the form (list of objects, list of values),
        with as many objects as values, and optionally a string
        naming the collection

    TODO:
      - a string, representing the ID of a collection.
      - a string, representing a text search.

    TODO:
    - ``max_results`` - (integer, default: 30) the maximum number of
      results to return, they are sorted according to their relevance. In
      any cases, the OEIS website will never provide more than 100 results.

    - ``first_result`` - (integer, default: 0) allow to skip the
      ``first_result`` first results in the search, to go further.
      This is useful if you are looking for a sequence that may appear
      after the 100 first found sequences.

    OUTPUT:

    - if ``query`` is an integer or a FindStat ID (e.g. 'St000045'), returns
      the associated FindStat statistic.

    TODO:
    - if ``query`` is a string, returns a tuple of OEIS sequences whose
      description corresponds to the query. Those sequences can be used
      without the need to fetch the database again.

    - if ``query`` is a list of integers, returns a tuple of OEIS sequences
      containing it as a subsequence. Those sequences can be used without
      the need to fetch the database again.

    EXAMPLES::

        sage: FindStat
        The Combinatorial Statistic Finder (http://findstat.org/)

    A particular statistic can be called by its St-number or number::

        sage: findstat('St000045')                   # optional -- internet
        St000045: The number of linear extensions of the tree. 

        sage: findstat(3)                            # optional -- internet
        St000003: The number of [[/StandardTableaux|standard Young tableaux]] of the partition.

    The database can be searched by providing a dictionary::

        sage: stat = {pi: pi.length() for pi in Permutations(3)}
        sage: search = findstat(stat) ; search    # optional -- internet
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 6)
        1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., ['inversion-number to major-index bijection'], 6)
        2: (St000067: The inversion number of the alternating sign matrix., ['to alternating sign matrix'], 6)
        3: (St000224: The sorting index of a permutation., ['first fundamental transformation'], 6)
        4: (St000008: The major index of the composition., ['inversion-number to major-index bijection', 'descent composition'], 6)
        5: (St000059: The inversion number of a standard Young tableau as defined by Haglund and Stevens. [1]  , ['inversion-number to major-index bijection', 'Robinson-Schensted recording tableau'], 6)
        6: (St000153: The number of adjacent cycles of a permutation., ['reverse', 'first fundamental transformation'], 6)
        7: (St000156: The Denert index of a permutation., ['inverse first fundamental transformation', 'foata_bijection'], 6)
        8: (St000174: The flush statistic on semistandard tableaux., ['to alternating sign matrix', 'to semistandard tableau'], 6)

        Each result is a triple consisting of the statistic, a list
        of strings naming certain maps, and a number which says how
        many of the values submitted agree with the values in the
        database, when applying the maps in the given order to the
        object (here: the permutation) and then computing the
        statistic on the result.

        We can access the description of a statistic as follows:

        sage: s = search[5][0]; print s.description()                        # optional -- internet
        The inversion number of a standard Young tableau as defined by Haglund and Stevens. [1]  

        Their inversion number is the total number of inversion pairs for the tableau.  An inversion pair is defined as a pair of cells (a,b), (x,y) such that the content of (x,y) is greater than the content of (a,b) and (x,y) is north of the inversion path of (a,b), where the inversion path is defined in detail in [1].

        We can access references provided as follows:
        sage: s.references()
        0: [1] J. Haglund and L. Stevens, An Extension of the Foata Map to Standard Young Tableaux, October, 2006

        TODO:
        A full-text search should also be available

    """
    def __call__(self, query, collection=None, depth=2):
        r"""
        See the documentation of :class:`FindStat`.

        TESTS::

            sage: findstat()
            Traceback (most recent call last):
            ...
            TypeError: __call__() takes at least 2 arguments (1 given)
        """
        if isinstance(query, str):
            if collection is None:
                if re.match('^St[0-9]{6}$',query):
                    return self.find_by_id(query)
                else:
                    raise ValueError("The value %s is not a valid statistic identifier." %query)
                    #                 return self.find_by_description(query, max_results, first_result)
            else:
                raise ValueError("When providing an identifier, do not provide a collection.")
        elif isinstance(query, (int, Integer)):
            if collection is None:
                return self.find_by_id(query)
            else:
                raise ValueError("When providing an identifier, do not provide a collection.")

        elif isinstance(query, dict):
            # values must be lists because otherwise we get a trailing comma
            stat = [(key, list(value)) if isinstance(value, (list, tuple)) else ([key], [value]) for (key, value) in query.iteritems()]
            return self.find_by_values(stat, collection, depth)

        elif isinstance(query, (list, tuple)):
            if (len(query) == 2 and isinstance(query[1], (list, tuple)) and
                query[1] != "" and isinstance(query[1][0], (int, Integer))):
                stat = [(query[0], list(query[1]))]
            else:
                stat = [(key, list(value)) if isinstance(value, (list, tuple)) else ([key], [value]) for (key, value) in query]
            return self.find_by_values(stat, collection, depth)

    def __repr__(self):
        r"""
        Return the representation of ``self``.

        TESTS::

            sage: findstat
            The Combinatorial Statistics Finder (http://findstat.org/)
        """
        return "The Combinatorial Statistics Finder (%s)" % findstat_url

    def find_by_id(self, ident):
        r"""

        INPUT:

        - ``ident`` - a string representing the St-number of the statistic
          or an integer representing its number.

        OUTPUT:

        - The FindStat statistic whose St-number or number corresponds to
          ``ident``.

        """
        if not isinstance(ident, str):
            ident = str(ident)
            ident = 'St000000'[:-len(ident)] + ident
        url = findstat_url_downloads + "statistics/" + ident + ".txt"
        statistic = _fetch(url) # we do not do any parsing here!
        return FindStatStatistic(statistic)

    def _parse(self, text):
        r"""
        parses the result of a findstat query
        """
        stat_re = r"Identifier (St[0-9]*) in the following ([0-9]*) positions"
        maps_re = r"the map \*(.*)\*"
        result = []
        start = 0
        stats = [(m.start(), m.groups()) for m in re.finditer(stat_re, text)]
        for (end, (stat, quality)) in stats:
            maps = [m.groups()[0] for m in re.finditer(maps_re, text[start:end])]
            result += [(findstat(stat), maps, Integer(quality))]
            start = end
        return result

    def find_by_values(self, statistic, collection, depth):
        r"""
        Search for FindStat statistics with the given values applying at most ``depth`` maps.

        INPUT:

        - ``statistic`` - a list of pairs of the form (list of objects, list of values),
                          each containing as many values as objects.
        - ``collection`` - a collection (possibly as a string) or None

        OUTPUT:

        - a tuple (with fancy formatting) of FindStat statistics.
          Those statistics can be used without the need to fetch the
          database again.

        """
        import urllib, urllib2

        if sum(len(values) for (keys, values) in statistic) > findstat_max_values:
            raise ValueError, "FindStat allows to search for at most " + str(findstat_max_values) + " values."

        if collection is None:
            collection = FindStatCollection(statistic[0][0][0])
        else:
            collection = FindStatCollection(collection)

        url = findstat_url_result + collection.url_name() + "/"

        to_str = collection.to_string()
        stat = [(map(to_str, keys), str(values)[1:-1]) for (keys, values) in statistic]

        stat_str = join([join(keys, "\n") + "\n====> " + values for (keys, values) in stat], "\n")

# the following was for single values
#
# the spaces in " => " are important!
#        stat_str = join([a + " => " + b for (keys, values) in stat], " \n")

        _ = verbose("Sending data to Findstat %s" %stat_str, caller_name='FindStat')
        values = urllib.urlencode({"freedata": stat_str, "depth": str(depth), "caller": "Sage"})

        _ = verbose("Fetching URL %s with encoded data %s" %(url, values), caller_name='FindStat')

        request = urllib2.Request(url, data=values)

        _ = verbose("Requesting %s" %request, caller_name='FindStat')

        response = urllib2.urlopen(request)

        _ = verbose("Response was %s" %response.info(), caller_name='FindStat')

        result_html = response.read()

        _ = verbose("Result was %s" %result_html, level = 2, caller_name='FindStat')

        download_idx = result_html.find("Download your result")

        _ = verbose("Download link at position %s" %download_idx, caller_name='FindStat')

        if download_idx == -1:
            return []
        else:
            link = findstat_url + result_html[download_idx-76:download_idx-18]
            
            _ = verbose("Download link is %s..." %link, caller_name='FindStat')
            result_text = urllib2.urlopen(link).read()
            return FancyTuple(self._parse(result_text))

    def browse(self):
        r"""
        Open the FindStat web page in a browser.

        EXAMPLES::

            sage: findstat.browse()                         # optional -- webbrowser
        """
        webbrowser.open(findstat_url)

    def submit(code=None):
        r"""
        Open the FindStat web page for submissions in a browser.
        """
        webbrowser.open(findstat_url_submit)
        

#     def _imaginary_entry(self, keywords=''):
#         r"""
#         This is an imaginary entry of an OEIS sequence for offline tests.
# 
#         INPUT:
# 
#         - ``keywords`` - a string corresponding to the keyword field of the
#           sequence.
# 
#         OUTPUT:
# 
#         - a string representing the entry of the sequence.
# 
#         TESTS::
# 
#             sage: oeis._imaginary_entry().split('\n')[0]
#             '%I A999999 M9999 N9999'
# 
#             sage: from sage.databases.oeis import OEISSequence
#             sage: keywords = 'simon,cussonet'
#             sage: s = OEISSequence(oeis._imaginary_entry(keywords))
#             sage: ','.join(s.keywords()) == keywords
#             True
# 
#         """
#         return ('%I A999999 M9999 N9999\n'
#                 '%S A999999 1,1,1,1,1,1,1,1,\n'
#                 '%T A999999 1,1,1,1,1,1,1,1,1,\n'
#                 '%U A999999 1,1,1,1,1,1,1,1,1\n'
#                 '%V A999999 1,1,1,1,-1,1,1,1,\n'
#                 '%W A999999 1,1,1,1,1,1,1,1,1,\n'
#                 '%X A999999 1,1,1,1,1,1,1,1,1\n'
#                 '%N A999999 The opposite of twice the characteristic sequence of 42 plus one, starting from 38.\n'
#                 '%D A999999 Lewis Carroll, Alice\'s Adventures in Wonderland.\n'
#                 '%D A999999 Lewis Carroll, The Hunting of the Snark.\n'
#                 '%D A999999 Deep Thought, The Answer to the Ultimate Question of Life, The Universe, and Everything.\n'
#                 '%H A999999 Wikipedia, <a href="http://en.wikipedia.org/wiki/42_(number)">42 (number)</a>\n'
#                 '%H A999999 See. also <a href="http://trac.sagemath.org/sage_trac/ticket/42">trac ticket #42</a>\n'
#                 '%H A999999 Do not confuse with the sequence <a href="/A000042">A000042</a> or the sequence <a href="/A000024">A000024</a>\n'
#                 '%H A999999 The string http://42.com is not a link.\n'
#                 '%F A999999 For n big enough, s(n+1) - s(n) = 0.\n'
#                 '%Y A999999 Related sequences are A000042 and its friend A000024.\n'
#                 '%A A999999 Anonymous.\n'
#                 '%O A999999 38,4\n'
#                 '%E A999999 This sequence does not contain errors.\n'
#                 '%e A999999 s(42) + s(43) = 0.\n'
#                 '%p A999999 Do not even try, Maple is not able to produce such a sequence.\n'
#                 '%t A999999 Mathematica neither.\n'
#                 '%o A999999 (Python)\n'
#                 '%o A999999 def A999999(n):\n'
#                 '%o A999999     assert(isinstance(n, (int, Integer))), "n must be an integer."\n'
#                 '%o A999999     if n < 38:\n'
#                 '%o A999999         raise ValueError("The value %s is not accepted." %str(n)))\n'
#                 '%o A999999     elif n == 42:\n'
#                 '%o A999999         return -1\n'
#                 '%o A999999     else:\n'
#                 '%o A999999         return 1\n'
#                 '%K A999999 ' + keywords + '\n'
#                 '%C A999999 42 is the product of the first 4 prime numbers, except 5 and perhaps 1.\n'
#                 '%C A999999 Apart from that, i have no comment.')
# 
#     def _imaginary_sequence(self, keywords='sign,easy'):
#         r"""
#         This is the OEIS sequence corresponding to the imaginary entry.
#         Its main purpose is to allow offline doctesting.
# 
#         INPUT:
# 
#         - ``keywords`` - string (default: 'sign,easy'), a list of words
#           separated by commas.
# 
#         OUTPUT:
# 
#         - OEIS sequence.
# 
#         TESTS::
# 
#             sage: s = oeis._imaginary_sequence()
#             sage: s
#             A999999: The opposite of twice the characteristic sequence of 42 plus one, starting from 38.
#             sage: s[4]
#             -1
#             sage: s(42)
#             -1
#         """
#         return OEISSequence(self._imaginary_entry(keywords))

class FindStatStatistic(SageObject):
    r"""
    The class of FindStat statistics.

    This class implements FindStat statistics. Such statistics are
    produced from a string in the FindStat format. They are usually
    produced by calls to the FindStat database, represented by the
    class :class:`FindStat`.

    .. automethod:: __call__
    """

    # separates different fields
    _field_separator1 = "-----------------------------------------------------------------------------\n"
    # separates field name from field content
    _field_separator2 = ":"
    # separates name from description
    _field_separator3 = "\r\n"
    # separates references
    _field_separator4 = "\r\n\r\n"

    # the fields of the FindStat database we expect
    _field_code = 'Code'
    _field_desc = 'Description'
    _field_copy = 'Copyright'
    _field_id   = 'Statistic identifier'
    _field_coll = 'Collection'
    _field_refs = 'References'
    _field_vals = 'Statistic values'

    def __init__(self, entry):
        r"""
        Initializes a FindStat statistic.

        TESTS::

            sage: stat = findstat('St000045')               # optional -- internet

            sage: s = findstat._imaginary_sequence()
        """
        # we delete leading and trailing whitespace for every field
        self._raw = entry
        self._fields = dict()
        split_entry = entry.split(FindStatStatistic._field_separator1)
        # first the copyright
        self._fields["Copyright"] = split_entry[0].strip()
        # the last element of the split entry should be empty
        for line in split_entry[1:-1:]:
            pair = line.partition(FindStatStatistic._field_separator2)
            if pair[1] == "":
                raise ValueError, line

            self._fields[pair[0]] = pair[2].strip()

    def id(self, format='St'):
        r"""
        The ID of the statistic ``self`` is the St-number that identifies
        ``self``.

        INPUT:

        - ``format`` - (string, default: 'St').

        OUTPUT:

        - if ``format`` is set to 'St', returns a string of the form 'St000123'.
        - if ``format`` is set to 'int' returns an integer of the form 123.

        """
        if format == 'St':
            return self._fields[FindStatStatistic._field_id]
        elif format == 'int':
            return Integer(self._fields[FindStatStatistic._field_id][2:].lstrip("0"))

    def raw_entry(self):
        r"""
        Return the raw entry of the statistic ``self``, in the FindStat format.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: print findstat(45).raw_entry()                   # optional -- internet
            *****************************************************************************
            *       www.FindStat.org - The Combinatorial Statistic Finder               *
            *                                                                           *
            *       Copyright (C) 2013 Chris Berg <chrisjamesberg@gmail.com>            *
            *                          Christian Stump <christian.stump@gmail.com>      *
            *                                                                           *
            *    This information is distributed in the hope that it will be useful,    *
            *    but WITHOUT ANY WARRANTY; without even the implied warranty of         *
            *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
            *****************************************************************************
            
            -----------------------------------------------------------------------------
            Statistic identifier: St000045
            
            -----------------------------------------------------------------------------
            Collection: Binary trees
            
            -----------------------------------------------------------------------------
            Description: The number of linear extensions of the tree. 
            
            Also, the number of increasing / decreasing binary trees labelled by $1, \dots, n$ of this shape.
            
            Also, the size of the sylvester class corresponding to this tree when the Tamari order is seen as a quotient poset of the right weak order on permutations. 
            
            Also, the number of permutations which give this tree shape when inserted in a binary search tree.
            
            Also, the number of permutations which increasing / decreasing tree is of this shape.
            
            -----------------------------------------------------------------------------
            References: [1] D. E. Knuth. The art of computer programming. Volume 3. Addison-Wesley Publishing
            Co., Reading, Mass.-London-Don Mills, Ont., 1973. Sorting and searching, Addison-Wesley
            Series in Computer Science and Information Processing
            
            [2] A. Bjorner and M. Wachs, q-Hook length formulas for forests, J. Combin. Theory Ser. A 52 (1989) 165–187.
            
            [3] F. Hivert, J.-C. Novelli, and J.-Y. Thibon. The algebra of binary search trees. Theoret. Comput. Sci.,339(1):129–165, 2005. [[arXiv:abs/math/0401089]]
            
            
            -----------------------------------------------------------------------------
            Code:
            def hook_product(tree):
                if(not tree):
                    return 1
                hl = hook_product(tree[0])
                hr = hook_product(tree[1])
                return hl*hr*tree.node_number()
                
            def linear_extensions(tree):
                return factorial(tree.node_number()-1)/hook_product(tree[0])/hook_product(tree[1])
            
            -----------------------------------------------------------------------------
            Statistic values:
            
            [.,[.,[.,.]]]             => 1
            [.,[[.,.],.]]             => 1



        """
        return self._raw

    def description(self):
        r"""
        Return the description of the statistic ``self``.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: print findstat(45).description()                              # optional -- internet
            The number of linear extensions of the tree. 
            
            Also, the number of increasing / decreasing binary trees labelled by $1, \dots, n$ of this shape.
            
            Also, the size of the sylvester class corresponding to this tree when the Tamari order is seen as a quotient poset of the right weak order on permutations. 
            
            Also, the number of permutations which give this tree shape when inserted in a binary search tree.
            
            Also, the number of permutations which increasing / decreasing tree is of this shape.

        """
        return self._fields[FindStatStatistic._field_desc]

    def name(self):
        r"""
        Return the name of the statistic ``self``, i.e., the first line of the description.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: findstat(45).name()                              # optional -- internet
            'The number of linear extensions of the tree. '

        """
        return self._fields[FindStatStatistic._field_desc].partition(FindStatStatistic._field_separator3)[0]


    def author(self):
        r"""
        Returns the author of the statistic in the encyclopedia.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: findstat(45).author()                            # optional -- internet

        """
        raise NotImplementedError
        return self._fields["XXX"]

    @cached_method
    def collection(self):
        r"""
        Returns the FindStat collection corresponding to the combinatorial collection, or `None`.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: findstat(45).collection()                      # optional -- internet
            Binary trees
        """
        return FindStatCollection(self._fields[FindStatStatistic._field_coll])

    @cached_method
    def first_terms(self):
        r"""

        OUTPUT:

        - dictionary.

        EXAMPLES::

            sage: f = findstat(2); f                  # optional -- internet
            St000002: The number of occurrences of the pattern $[1,2,3]$ inside a permutation of length at least 3.
            sage: f.first_terms()[Permutation([1, 2, 4, 3])]
            2

        """
        vals = self._fields[FindStatStatistic._field_vals].split()
        # elements with index
        # 0 mod 3 are the objects
        vals0 = vals[0::3]
        # 1 mod 3 should always be "=>"
        # 2 mod 3 are the values
        vals2 = vals[2::3]
        from_str = self.collection().from_string()
        if from_str is None:
            return {obj: Integer(val) for (obj, val) in zip(vals0, vals2)}
        else:
            return {from_str(obj): Integer(val) for (obj, val) in zip(vals0, vals2)}

    def code(self, security_warning=True):
        r"""
        Returns the contents of the field for the code for the statistic as a string.

        Warning: although it may be desirable to execute this code,
        be aware of the fact that it is user contributed code.  It
        could kill your computer.  Therefore, this method prepends
        the following two lines to the contents of the field, when
        ``security_warning`` is ``True``:

        print "Warning: this is user contributed code, potentially killing your computer"
        raise StandardError

        OUTPUT:

        - string.

        EXAMPLES::

            sage: print findstat(45).code()                            # optional -- internet
            print "Warning: this is user contributed code, potentially killing your computer"
            raise StandardError

            def hook_product(tree):
                if(not tree):
                    return 1
                hl = hook_product(tree[0])
                hr = hook_product(tree[1])
                return hl*hr*tree.node_number()

            def linear_extensions(tree):
                return factorial(tree.node_number()-1)/hook_product(tree[0])/hook_product(tree[1])

        """
        if security_warning:
            s = 'raise StandardError("This is user contributed code, potentially killing your computer!")'
            return s + '\r\n' + self._fields[FindStatStatistic._field_code]
        else:
            return self._fields[FindStatStatistic._field_code]

    def __repr__(self):
        r"""
        Prints the statistic number and a short summary of this statistic.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: findstat(45)                      # optional -- internet
            St000045: The number of linear extensions of the tree. 
        """
        return "%s: %s" % (self.id(), self.name())

    def __call__(self, obj):
        r"""
        Returns the value of the statistic ``self`` for ``obj``.

        INPUT:

        - ``obj`` - an element of the combinatorial collection.

        OUTPUT:

        - integer.

        EXAMPLES::

            sage: f=findstat(2); f
            St000002: The number of occurrences of the pattern $[1,2,3]$ inside a permutation of length at least 3.

            sage: f(Permutation([1,2,4,3,5]))        # optional -- internet
            7

        ::

        TESTS::

            sage: s = findstat._imaginary_sequence()
            sage: s(38)
            1
            sage: s(42)
            -1
            sage: s(2)
            Traceback (most recent call last):
            ...
            ValueError: Sequence A999999 is not defined (or known) for index 2
        """
        vals = self.first_terms()
        if isinstance(obj, str):
            from_str = self.collection().from_string()
            if from_str is None:
                return vals[obj]
            else:
                return vals[from_str(obj)]
        else:
            return vals[obj]

    def __eq__(self,other):
        r"""
        Return ``True`` if ``self`` is equal to ``other`` and
        ``False`` otherwise.  Two statistics are considered equal if
        they have the same FindStat id.

        INPUT:

        - ``other`` -- a FindStat statistic

        OUTPUT:

        - boolean.

        EXAMPLES::

            sage: findstat([1,2,3,5,8,13])[0] == findstat(45)   # optional -- internet
            True

        TESTS::

            sage: s = findstat._imaginary_sequence()
            sage: s == findstat._imaginary_sequence()
            True
        """
        return self.id() == other.id()

    def __ne__(self, other):
        r"""
        Return ``True`` if ``self`` has a different FindStat id than ``other``
        and ``False`` otherwise.

        INPUT:

        - ``other`` -- a findstat statistic

        OUTPUT:

        - boolean.

        EXAMPLES::

            sage: findstat([1,2,3,5,8,13])[0] != findstat(40)   # optional -- internet
            True

        TESTS::

            sage: s = findstat._imaginary_sequence()
            sage: s != findstat._imaginary_sequence()
            False
        """
        return not self.__eq__(other)

    def references(self):
        r"""
        Return a tuple of references associated to the statistic ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting)

        EXAMPLES::

            sage: w = findstat(21) ; w                # optional -- internet
            St000021: The number of descents of a permutation.

            sage: w.references()                        # optional -- internet
            0: [[/Permutations/Descents]]
            1: [[oeis:A008292]]

            sage: _[0]                                  # optional -- internet
            'A. H. Beiler, Recreations in the Theory of Numbers, Dover, NY, 1964, p. 52.'

        """
        return FancyTuple(self._fields[FindStatStatistic._field_refs].split(FindStatStatistic._field_separator4))

    def url(self):
        r"""
        Return the URL of the page associated to the statistic ``self``.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: findstat(45).url()                               # optional -- internet
            'http://www.findstat.org/St000045'

        """
        return findstat_url + self.id()

    def browse(self):
        r"""
        Open the FindStat web page associated to the statistic ``self`` in a browser.

        EXAMPLES::

            sage: findstat(45).browse()                            # optional -- internet, webbrowser
        """
        webbrowser.open(self.url())

    def edit(self):
        r"""
        Open the FindStat web page for editing the statistic ``self`` in a browser.

        EXAMPLES::

            sage: findstat(45).edit()                            # optional -- internet, webbrowser

        """
        webbrowser.open(findstat_url_edit + self.id())

class FindStatCollection():
    r"""
    The class of FindStat collections.

    .. automethod:: __call__
    """
    # we set up a dictionary from FindStat collections to 
    # * url's as needed by FindStat, 
    # * sage objects
    # * a method to get the size of the sage object
    # * the (FindStat) string representations of the sage object
    # * sage constructors of the FindStat string representation
    _findstat_collections = {
        'Alternating sign matrices': ("AlternatingSignMatrices", AlternatingSignMatrix, 
                                      lambda x: x.to_matrix().nrows(),
                                      lambda x: str(map(list,list(x._matrix))), 
                                      lambda x: AlternatingSignMatrix(literal_eval(x))),
        'Binary trees':              ("BinaryTrees",             BinaryTree,            
                                      lambda x: x.node_number(),
                                      str, 
                                      lambda x: BinaryTree(x)),
        'Cores':                     ("Cores",                   Core,
                                      lambda x: x.k(),
                                      lambda X: "( "+X._repr_()+", "+str(X.k())+" )",
                                      lambda x: Core(literal_eval(x)[0],literal_eval(x)[1])),
        'Dyck paths':                ("DyckPaths",               DyckWord,              
                                      lambda x: x.length(),
                                      lambda x: str(list(DyckWord(x))), 
                                      lambda x: DyckWord(literal_eval(x))),
        'Finite Cartan types':       ("FiniteCartanType",        CartanType_abstract,            
                                      lambda x: x.rank(),
                                      str,
                                      lambda x: CartanType(*literal_eval(str(x)))),
        'Gelfand-Tsetlin patterns':  ("GelfandTsetlinPatterns",  GelfandTsetlinPattern, 
                                      len,
                                      str,                              
                                      lambda x: GelfandTsetlinPattern(literal_eval(x))),
        'Graphs':                    ("Graphs",                  Graph,                 
                                      lambda x: x.num_verts(),
                                      lambda X: str((sorted(X.canonical_label().edges(False)), X.num_verts())),
                                      lambda x: Graph([range(literal_eval(x)[1]), lambda i,j: (i,j) in literal_eval(x)[0] or (j,i) in literal_eval(x)[0]], immutable=True)),
        'Integer compositions':      ("IntegerCompositions",     Composition,           
                                      lambda x: x.size(),
                                      str,                              
                                      lambda x: Composition(literal_eval(x))),
        'Integer partitions':        ("IntegerPartitions",       Partition,             
                                      lambda x: x.size(),
                                      str,                              
                                      lambda x: Partition(literal_eval(x))),
        'Ordered trees':             ("OrderedTrees",            OrderedTree,           
                                      lambda x: x.node_number(),
                                      str,                              
                                      lambda x: OrderedTree(literal_eval(x))),
        'Parking functions':         ("ParkingFunctions",        ParkingFunction_class, 
                                      len,
                                      str,                              
                                      lambda x: ParkingFunction(literal_eval(x))),
        'Perfect matchings':         ("PerfectMatchings",        PerfectMatching,       
                                      lambda x: x.size(),
                                      str,                              
                                      lambda x: PerfectMatching(literal_eval(x))),
        'Permutations':              ("Permutations",            Permutation,           
                                      lambda x: x.size(),
                                      str,                              
                                      lambda x: Permutation(literal_eval(x))),
        'Posets':                    ("Posets",                  FinitePoset,           
                                      lambda x: x.cardinality(),
                                      lambda X: str((sorted(X._hasse_diagram.canonical_label().cover_relations()), len(X._hasse_diagram.vertices()))),                                                                       
                                      lambda x: Poset((range(literal_eval(x)[1]), literal_eval(x)[0]))),
        'Semistandard tableaux':     ("SemistandardTableaux",    SemistandardTableau,   
                                      lambda x: x.size(),
                                      str,                              
                                      lambda x: SemistandardTableau(literal_eval(x))),
        'Set partitions':            ("SetPartitions",           SetPartition,          
                                      lambda x: x.size(),
                                      str,                              
                                      lambda x: SetPartition(literal_eval(x.replace('{','[').replace('}',']')))),
        'Standard tableaux':         ("StandardTableaux",        StandardTableau,       
                                      lambda x: x.size(),
                                      str,                              
                                      lambda x: StandardTableau(literal_eval(x)))}

    def __init__(self, entry):
        r"""
        Initializes a FindStat collection from
        
        * a string, eg. 'Dyck paths'
        * a sage object

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection('Dyck paths'); c               # optional -- internet
            Dyck paths

        """
        if isinstance(entry, str):
            # we should make this case insensitive, and generally more flexible 
            value = FindStatCollection._findstat_collections.get(entry, None)
            if value == None:
                raise ValueError, "unknown collection '" + entry + "'"
            else:
                self._name = entry
                (self._url_name, self._sagename, self._to_size, self._to_str, self._from_str) = value
        else:
            bad = True
            for (key, value) in FindStatCollection._findstat_collections.iteritems():
                if isinstance(entry, value[1]):
                    self._name = key
                    (self._url_name, self._sagename, self._to_size, self._to_str, self._from_str) = value
                    bad = False
                    break
            if bad:
                raise ValueError, "Could not find FindStat collection for " + str(entry)

    def url_name(self):
        return self._url_name

    def to_size(self):
        r"""

        returns a function that returns the FindStat size of a sage object
        """
        return self._to_size


    def to_string(self):
        r"""

        returns a function that converts a sage object into a FindStat string
        """
        return self._to_str

    def from_string(self):
        r"""
        Returns the sage constructor corresponding to the combinatorial collection, or `None`.

        OUTPUT:

        - string or `None`.

        """
        return self._from_str

    def __repr__(self):
        r"""
        Prints the name of the collection

        OUTPUT:

        - string.

        EXAMPLES::

            sage: FindStatCollection("Binary trees")                      # optional -- internet
            Binary trees
        """
        return self._name

    def name(self):
        r"""
        returns the name of the collection

        OUTPUT:

        - string.

        EXAMPLES::

            sage: FindStatCollection("Binary trees")                      # optional -- internet
            Binary trees
        """
        return self._name


    def __eq__(self,other):
        r"""
        Return ``True`` if ``self`` is equal to ``other`` and
        ``False`` otherwise.  Two collections are considered equal if
        they have the same FindStat name.

        INPUT:

        - ``other`` -- a FindStat collection or a string

        OUTPUT:

        - boolean.

        EXAMPLES::

            sage: findstat(1).collection() == "Permutations"   # optional -- internet
            True

        TESTS::

            sage: s = findstat._imaginary_sequence()
            sage: s == findstat._imaginary_sequence()
            True
        """
        return self._name == str(other)

    def __ne__(self, other):
        r"""
        Return ``True`` if ``self`` has a different FindStat id than ``other``
        and ``False`` otherwise.

        INPUT:

        - ``other`` -- a findstat statistic

        OUTPUT:

        - boolean.

        EXAMPLES::

            sage: findstat([1,2,3,5,8,13])[0] != findstat(40)   # optional -- internet
            True

        TESTS::

            sage: s = findstat._imaginary_sequence()
            sage: s != findstat._imaginary_sequence()
            False
        """
        return not self.__eq__(other)



findstat = FindStat()
