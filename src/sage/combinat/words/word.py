# coding=utf-8
r"""
Word classes

AUTHORS:

- Arnaud Bergeron
- Amy Glen
- Sébastien Labbé
- Franco Saliola

"""
#*****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#                          Amy Glen <amy.glen@gmail.com>,
#                          Sébastien Labbé <slabqc@gmail.com>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.combinat.words.abstract_word import Word_class
from sage.combinat.words.finite_word import FiniteWord_class
from sage.combinat.words.infinite_word import InfiniteWord_class
from word_datatypes import (WordDatatype_str,
                            WordDatatype_list,
                            WordDatatype_tuple)
from word_infinite_datatypes import (
                            WordDatatype_iter_with_caching,
                            WordDatatype_iter,
                            WordDatatype_callable_with_caching,
                            WordDatatype_callable)

# TODO. Word needs to be replaced by Word. Consider renameing
# Word_class to Word and imbedding Word as its __call__ method.

def Word(data=None, alphabet=None, length=None, datatype=None, caching=True):
    r"""
    Construct a word.

    INPUT:

    -  ``data`` - (default: None) list, string, tuple, iterator, None
       (shorthand for []), or a callable defined on [0,1,...,length].

    -  ``alphabet`` - any argument accepted by Words

    -  ``length`` - (default: None) This is dependent on the type of data.
       It is ignored for words defined by lists, strings, tuples,
       etc., because they have a naturally defined length.
       For callables, this defines the domain of definition,
       which is assumed to be [0, 1, 2, ..., length-1].
       For iterators: Infinity if you know the iterator will not
       terminate (default); "unknown" if you do not know whether the
       iterator terminates; "finite" if you know that the iterator
       terminates, but do know know the length.

    -  ``datatype`` - (default: None) None, "list", "str", "tuple", "iter",
       "callable". If None, then the function
       tries to guess this from the data.
    -  ``caching`` - (default: True) True or False. Whether to keep a cache
       of the letters computed by an iterator or callable.

    .. note::

       Be careful when defining words using callables and iterators. It
       appears that islice does not pickle correctly causing various errors
       when reloading. Also, most iterators do not support copying and
       should not support pickling by extension.

    EXAMPLES:

    Empty word::

        sage: Word()
        word:

    Word with string::

        sage: Word("abbabaab")
        word: abbabaab

    Word with string constructed from other types::

        sage: Word([0,1,1,0,1,0,0,1], datatype="str")
        word: 01101001
        sage: Word((0,1,1,0,1,0,0,1), datatype="str")
        word: 01101001

    Word with list::

        sage: Word([0,1,1,0,1,0,0,1])
        word: 01101001

    Word with list constructed from other types::

        sage: Word("01101001", datatype="list")
        word: 01101001
        sage: Word((0,1,1,0,1,0,0,1), datatype="list")
        word: 01101001

    Word with tuple::

        sage: Word((0,1,1,0,1,0,0,1))
        word: 01101001

    Word with tuple constructed from other types::

        sage: Word([0,1,1,0,1,0,0,1], datatype="tuple")
        word: 01101001
        sage: Word("01101001", datatype="str")
        word: 01101001

    Word with iterator::

        sage: from itertools import count
        sage: Word(count())
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: Word(iter("abbabaab")) # iterators default to infinite words
        word: abbabaab
        sage: Word(iter("abbabaab"), length="unknown")
        word: abbabaab
        sage: Word(iter("abbabaab"), length="finite")
        word: abbabaab

    Word with function (a 'callable')::

        sage: f = lambda n : add(Integer(n).digits(2)) % 2
        sage: Word(f)
        word: 0110100110010110100101100110100110010110...
        sage: Word(f, length=8)
        word: 01101001

    Word over a string with a parent::

        sage: w = Word("abbabaab", alphabet="abc"); w
        word: abbabaab
        sage: w.parent()
        Words over Ordered Alphabet ['a', 'b', 'c']

    The default parent is the combinatorial class of all words::

        sage: w = Word("abbabaab"); w
        word: abbabaab
        sage: w.parent()
        Words

    TESTS::

        sage: Word(5)
        Traceback (most recent call last):
        ...
        ValueError: Cannot guess a datatype from data (=5); please specify one

    ::

        sage: w  = Word('abc')
        sage: w is Word(w)
        True
        sage: w is Word(w, alphabet='abc')
        False
    """
    # Create the parent object
    from words import Words
    parent = Words() if alphabet is None else Words(alphabet)

    return parent(data=data, length=length, datatype=datatype, caching=caching)

###########################################################################
##### DEPRECATION WARNINGS ################################################
##### Added July 2009 #####################################################
###########################################################################

def is_Word(obj):
    r"""
    Returns True if obj is a word, and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.word import is_Word
        sage: is_Word(33)
        doctest:1: DeprecationWarning: is_Word is deprecated, use isinstance(your_object, Word_all) instead!
        False
        sage: is_Word(Word('abba'))
        True
    """
    from sage.misc.misc import deprecation
    deprecation("is_Word is deprecated, use isinstance(your_object, Word_all) instead!")
    return isinstance(obj, Word_class)

def is_FiniteWord(obj):
    r"""
    Returns True if obj is a finite word, and False otherwise.

    EXAMPLES::

        sage: from sage.combinat.words.word import is_FiniteWord
        sage: is_FiniteWord(33)
        doctest:1: DeprecationWarning: is_Word is deprecated, use isinstance(your_object, Word_all) instead!
        False
        sage: is_FiniteWord(Word('baab'))
        True
    """
    from sage.misc.misc import deprecation
    deprecation("is_Word is deprecated, use isinstance(your_object, Word_all) instead!")
    if isinstance(obj, Word_class):
        try:
            len(obj)
        except:
            return False
        return True
    else:
        return False

#######################################################################
#                                                                     #
#                    Concrete word classes                            #
#                                                                     #
#######################################################################

##### Finite Words #####

class FiniteWord_list(WordDatatype_list, FiniteWord_class):
    r"""
    Finite word represented by a Python list.

    For any word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word(range(10))
        sage: w.iterated_right_palindromic_closure()
        word: 0102010301020104010201030102010501020103...

    TESTS::

        sage: w = Word([0,1,1,0])
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_str(WordDatatype_str, FiniteWord_class):
    r"""
    Finite word represented by a Python str.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word('abcdef')
        sage: w.is_square()
        False

    TESTS::

        sage: w = Word('abba')
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_tuple(WordDatatype_tuple, FiniteWord_class):
    r"""
    Finite word represented by a Python tuple.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word(())
        sage: w.is_empty()
        True

    TESTS::

        sage: w = Word((0,1,1,0))
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_iter_with_caching(WordDatatype_iter_with_caching, FiniteWord_class):
    r"""
    Finite word represented by an iterator (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word(iter('abcdef'))
        sage: w.conjugate(2)
        word: cdefab

    TESTS::

        sage: w = Word(iter(range(10)))
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_list'>
    """
    pass

class FiniteWord_iter(WordDatatype_iter, FiniteWord_class):
    r"""
    Finite word represented by an iterator.

    For such word `w`, type  ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word(iter(range(10)), caching=False)
        sage: w
        word: 0123456789
        sage: w.finite_differences()
        word: 111111111

    TESTS::

        sage: w = Word(iter(range(10)), caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_iter'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_list'>
    """
    pass

class FiniteWord_callable_with_caching(WordDatatype_callable_with_caching, FiniteWord_class):
    r"""
    Finite word represented by a callable (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: f = lambda n : n % 3
        sage: w = Word(f, length=32)
        sage: w
        word: 01201201201201201201201201201201
        sage: w.border()
        word: 01201201201201201201201201201

    TESTS::

        sage: w = Word(lambda n:n, length=10)
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>

    Pickle also works for concatenation of words::

        sage: w = Word(range(10)) * Word('abcdef')
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_list'>

    Pickle also works for power of words::

        sage: w = Word(range(10)) ^ 2
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_list'>
    """
    pass

class FiniteWord_callable(WordDatatype_callable, FiniteWord_class):
    r"""
    Finite word represented by a callable.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: f = lambda n : 3 if n > 8 else 6
        sage: w = Word(f, length=30, caching=False)
        sage: w
        word: 666666666333333333333333333333
        sage: w.is_symmetric()
        True

    TESTS::

        sage: w = Word(lambda n:n, length=10, caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_callable'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_callable'>
    """
    pass

##### Infinite Words #####

class InfiniteWord_iter_with_caching(WordDatatype_iter_with_caching, InfiniteWord_class):
    r"""
    Infinite word represented by an iterable (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: from itertools import cycle
        sage: w = Word(cycle([9,8,4]))
        sage: w
        word: 9849849849849849849849849849849849849849...
        sage: prefix = w[:23]
        sage: prefix
        word: 98498498498498498498498
        sage: prefix.minimal_period()
        3

    TESTS::

        sage: from itertools import count
        sage: w = Word(count())
        sage: type(w)
        <class 'sage.combinat.words.word.InfiniteWord_iter_with_caching'>

    Pickle is not supported for infinite word defined by an iterator::

        sage: dumps(w)
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle <type 'generator'>: attribute lookup __builtin__.generator failed
    """
    pass

class InfiniteWord_iter(WordDatatype_iter, InfiniteWord_class):
    r"""
    Infinite word represented by an iterable.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: from itertools import chain, cycle
        sage: w = Word(chain('letsgo', cycle('forever')), caching=False)
        sage: w
        word: letsgoforeverforeverforeverforeverforeve...
        sage: prefix = w[:100]
        sage: prefix
        word: letsgoforeverforeverforeverforeverforeve...
        sage: prefix.is_lyndon()
        False

    TESTS::

        sage: from itertools import count
        sage: w = Word(count(), caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.InfiniteWord_iter'>

    Pickle is not supported for infinite word defined by an iterator::

        sage: dumps(w)
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle <type 'generator'>: attribute lookup __builtin__.generator failed
    """
    pass

class InfiniteWord_callable_with_caching(WordDatatype_callable_with_caching, InfiniteWord_class):
    r"""
    Infinite word represented by a callable (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: w = Word(lambda n:n)
        sage: factor = w[4:13]
        sage: factor
        word: 4,5,6,7,8,9,10,11,12

    TESTS::

        sage: w = Word(lambda n:n)
        sage: type(w)
        <class 'sage.combinat.words.word.InfiniteWord_callable_with_caching'>
        sage: z = loads(dumps(w))
        sage: z
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: type(z)
        <class 'sage.combinat.words.word.InfiniteWord_callable_with_caching'>
    """
    pass

class InfiniteWord_callable(WordDatatype_callable, InfiniteWord_class):
    r"""
    Infinite word represented by a callable.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: w = Word(lambda n:n, caching=False)
        sage: w
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: w.iterated_right_palindromic_closure()
        word: 0102010301020104010201030102010501020103...

    TESTS::

        sage: w = Word(lambda n:n, caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.InfiniteWord_callable'>
        sage: z = loads(dumps(w))
        sage: z
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: type(z)
        <class 'sage.combinat.words.word.InfiniteWord_callable'>
    """
    pass

##### Words of unknown length #####

class Word_iter_with_caching(WordDatatype_iter_with_caching, Word_class):
    r"""
    Word of unknown length (finite or infinite) represented by an
    iterable (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: w = Word(iter([1,2,3]*1000), length='unknown')
        sage: w
        word: 1231231231231231231231231231231231231231...
        sage: w.finite_differences(mod=2)
        word: 1101101101101101101101101101101101101101...

    TESTS::

        sage: w = Word(iter('abcd'*100), length='unknown')
        sage: type(w)
        <class 'sage.combinat.words.word.Word_iter_with_caching'>
        sage: w
        word: abcdabcdabcdabcdabcdabcdabcdabcdabcdabcd...

    Pickle is not supported for word of unknown length defined by an iterator::

        sage: dumps(w)
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle <type 'generator'>: attribute lookup __builtin__.generator failed
    """
    pass

class Word_iter(WordDatatype_iter, Word_class):
    r"""
    Word of unknown length (finite or infinite) represented by an
    iterable.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: w = Word(iter([1,1,4,9]*1000), length='unknown', caching=False)
        sage: w
        word: 1149114911491149114911491149114911491149...
        sage: w.delta()
        word: 2112112112112112112112112112112112112112...

    TESTS::

        sage: w = Word(iter('abcd'*100), length='unknown', caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.Word_iter'>
        sage: w
        word: abcdabcdabcdabcdabcdabcdabcdabcdabcdabcd...

    Pickle is not supported for word of unknown length defined by an iterator::

        sage: dumps(w)
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle <type 'generator'>: attribute lookup __builtin__.generator failed
    """
    pass

#######################################################################
#                                                                     #
#                    Old word classes                                 #
#                                                                     #
#######################################################################

def _word_from_word_content(data, parent):
    r"""
    Create a word with the given parent from a :class:`WordContent` object.

    .. note::

       :class:`WordContent` is deprecated and will be deleted in a future
       version of Sage. This function exists to maintain backwards
       compatibility for unpickling objects saved with older versions of Sage.

    TESTS::

        sage: from sage.combinat.words.word import _word_from_word_content
        sage: import sage.combinat.words.word_content as word_content
        sage: cl = word_content.WordContentFromList([0, 1, 1, 2, 1])
        doctest:...: DeprecationWarning: WordContentFromList is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        sage: _word_from_word_content(cl,Words([0,1,2]))
        word: 01121
        sage: _word_from_word_content(cl,Words([2,1,0]))
        word: 21101

    ::

        sage: cf = word_content.WordContentFromFunction(lambda x : x)[:10]
        doctest:...: DeprecationWarning: WordContentFromFunction is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        doctest:...: DeprecationWarning: WordContentFromFunction is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        sage: _word_from_word_content(cf, Words())
        word: 0123456789

    ::

        sage: from itertools import count
        sage: ci = word_content.WordContentFromIterator(count(3))[:5]
        doctest:...: DeprecationWarning: WordContentFromIterator is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        doctest:...: DeprecationWarning: WordContentFromIterator is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        sage: _word_from_word_content(ci, Words())
        word: 34567

    ::

        sage: cc = word_content.ConcatenateContent((cl, cf, ci))
        doctest:...: DeprecationWarning: ConcatenateContent is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        doctest:...: DeprecationWarning: ConcatenateContent is deprecated and will be deleted in a future version of Sage; see sage.combinat.words.word_datatypes for alternatives
        sage: _word_from_word_content(cc, Words())
        word: 01121012345678934567

    """
    from sage.combinat.words import word_content
    # extract data from WordContent object
    length = None
    if isinstance(data, word_content.WordContentFromList):
        unranker = parent.alphabet().unrank
        data = map(unranker, data._list)
        datatype = 'list'
    elif isinstance(data, word_content.WordContentFromIterator):
        data = data._get_it()
        datatype = 'iter'
    elif isinstance(data, word_content.WordContentFromFunction):
        length = data._len
        data = data._func
        datatype = 'callable'
    elif isinstance(data, word_content.ConcatenateContent):
        l = map(_word_from_word_content, data)
        data = CallableFromListOfWords(l)
        datatype = 'callable'
    else:
        raise TypeError, 'data is not an instance of WordContent'
    return Word(data=data, alphabet=parent, datatype=datatype, length=length)

class DeprecatedWordClass(SageObject):
    def __setstate__(self, data):
        from sage.misc.misc import deprecation
        cls = self.__class__
        deprecation('Your word object is saved in an old file format '
            'since %s is deprecated and will be deleted in a future version '
            'of Sage (you can use %s instead). You can re-save your word by '
            'typing "word.save(filename)" to ensure that it will load in '
            'future versions of Sage.''' %
            (cls.__name__, cls.__bases__[1].__name__))
        parent = data['_parent']
        del data['_parent']
        wordcontent = data.get('_word_content', None)
        del data['_word_content']
        w = _word_from_word_content(wordcontent, parent)
        self.__class__ = type(w)
        self.__init__(parent, list(w))
        for key, item in data.iteritems():
            setattr(self, key, item)

class AbstractWord(DeprecatedWordClass, FiniteWord_list):
    pass

class AbstractFiniteWord(DeprecatedWordClass, FiniteWord_list):
    pass

class AbstractInfiniteWord(DeprecatedWordClass, InfiniteWord_class):
    pass

class Word_over_Alphabet(DeprecatedWordClass, FiniteWord_list):
    pass

class Word_over_OrderedAlphabet(DeprecatedWordClass, FiniteWord_list):
    pass

class FiniteWord_over_Alphabet(DeprecatedWordClass, FiniteWord_list):
    pass

class FiniteWord_over_OrderedAlphabet(DeprecatedWordClass, FiniteWord_list):
    pass

class InfiniteWord_over_Alphabet(DeprecatedWordClass, InfiniteWord_class):
    pass

class InfiniteWord_over_OrderedAlphabet(DeprecatedWordClass, InfiniteWord_class):
    pass

